#!/usr/bin/env python3
"""
nist-openfoam-wallboiling.py

Generates OpenFOAM-style isobaric property tables from the NIST WebBook
for a species, a pressure range, and a temperature range of choice.

For LIQUID table, for T > Tsat(P):
    - rho  : ramps LINEARLY from saturated-liquid density at Tsat down to a
             floor value, then holds constant.
    - hs, Cp, Cv, mu, kappa: held CONSTANT at saturated-liquid value at Tsat(P).

For VAPOUR table, for T < Tsat(P):
    - rho  : extrapolated with the ideal-gas isobaric scaling law:
                 rho(T) = rho_sat_vapour(P) * Tsat(P) / T
    - hs, Cp, Cv, mu, kappa: held CONSTANT at saturated-vapour value at Tsat(P).

Units for input:
    Temperature  : K
    Pressure     : MPa

Output is converted to SI system of units for writing OpenFOAM files.

Usage:
    python nist-openfoam-wallboiling.py --species R134a \
        --p-low 2.2 --p-high 3.2 --p-inc 0.04 \
        --t-low 325 --t-high 425 --t-inc 4
"""

import argparse
import io
import os
import sys
import time

import pandas as pd
import requests

BASE_FLUID_URL = "https://webbook.nist.gov/cgi/fluid.cgi"

HEADERS = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}

SESSION = requests.Session()
SESSION.headers.update(HEADERS)

REQUEST_DELAY_S = 0.35
DEBUG = False

CAS_LOOKUP = {
    "water": "7732-18-5",
    "nitrogen": "7727-37-9",
    "oxygen": "7782-44-7",
    "carbon dioxide": "124-38-9",
    "co2": "124-38-9",
    "methane": "74-82-8",
    "ethane": "74-84-0",
    "propane": "74-98-6",
    "n-butane": "106-97-8",
    "butane": "106-97-8",
    "isobutane": "75-28-5",
    "ammonia": "7664-41-7",
    "carbon monoxide": "630-08-0",
    "hydrogen": "1333-74-0",
    "helium": "7440-59-7",
    "argon": "7440-37-1",
    "ethanol": "64-17-5",
    "methanol": "67-56-1",
    "acetone": "67-64-1",
    "benzene": "71-43-2",
    "toluene": "108-88-3",
    "r134a": "811-97-2",
    "sulfur dioxide": "7446-09-5",
    "hydrogen sulfide": "7783-06-4",
    "nitrous oxide": "10024-97-2",
    "xenon": "7440-63-3",
    "krypton": "7439-90-9",
    "neon": "7440-01-9",
}

COMMON_PARAMS = dict(
    TUnit="K",
    PUnit="MPa",
    DUnit="kg/m3",
    HUnit="kJ/kg",
    WUnit="m/s",
    VisUnit="Pa*s",
    STUnit="N/m",
)
PROPS = ["rho", "hs", "Cp", "Cv", "mu", "kappa"]


def frange(low, high, inc):
    if inc <= 0:
        return [low]
    n = int(round((high - low) / inc))
    return [
        round(low + i * inc, 10)
        for i in range(n + 1)
        if low + i * inc <= high + 1e-9
    ]


def _get_table(params, digits=6):
    p = dict(Action="Data", Wide="on", Digits=digits)
    p.update(COMMON_PARAMS)
    p.update(params)

    try:
        resp = SESSION.get(BASE_FLUID_URL, params=p, timeout=30)
        time.sleep(REQUEST_DELAY_S)
    except requests.RequestException as exc:
        if DEBUG:
            print(f"    [debug] Request failed: {exc}")
        return pd.DataFrame()

    if DEBUG:
        print(f"    [debug] GET {resp.url}")
        print(
            f"    [debug] status={resp.status_code} bytes={len(resp.content)}"
        )

    try:
        resp.raise_for_status()
    except requests.HTTPError as exc:
        if DEBUG:
            print(f"    [debug] HTTP error: {exc}")
        return pd.DataFrame()

    try:
        tables = pd.read_html(io.BytesIO(resp.content))
        if tables:
            df = max(tables, key=lambda t: t.shape[0]).copy()
            if isinstance(df.columns, pd.MultiIndex):
                df.columns = [
                    " ".join(
                        str(c) for c in col if "Unnamed" not in str(c)
                    ).strip()
                    for col in df.columns.values
                ]
            if DEBUG:
                print(f"    [debug] parsed HTML table shape={df.shape}")
            return df
    except Exception as exc:
        if DEBUG:
            print(f"    [debug] pd.read_html failed: {exc}")

    try:
        df = pd.read_csv(io.BytesIO(resp.content), sep="\t")
        df.columns = [str(c).strip() for c in df.columns]
        if not df.empty and df.shape[1] > 1:
            if DEBUG:
                print(f"    [debug] parsed TSV table shape={df.shape}")
                print(f"    [debug] columns={list(df.columns)}")
            return df
    except Exception as exc:
        if DEBUG:
            print(f"    [debug] pd.read_csv TSV failed: {exc}")

    if DEBUG:
        snippet = resp.text[:800].replace("\n", " ")
        print(f"    [debug] response snippet: {snippet!r}")

    return pd.DataFrame()


def _find_col(df, *substrings):
    for c in df.columns:
        cl = str(c).lower()
        if all(s.lower() in cl for s in substrings):
            return c
    return None


def _find_phase_col(df, phase, *prop_substrings):
    phase_keys = {
        "liquid": ["liquid", "(l,", "l,", "(l)", " l "],
        "vapor": ["vapor", "vapour", "(v,", "v,", "(v)", " v "],
    }[phase]

    for c in df.columns:
        cl = str(c).lower()
        if all(ps.lower() in cl for ps in prop_substrings):
            if any(pk.lower() in cl for pk in phase_keys):
                return c
    return None


def fetch_isobaric(species_id, p_mpa, t_low, t_high, t_inc):
    if t_high < t_low:
        return pd.DataFrame()
    
    safe_t_inc = max(t_inc, 0.1)

    df = _get_table(
        dict(
            ID=species_id,
            Type="IsoBar",
            P=p_mpa,
            TLow=t_low,
            THigh=t_high,
            TInc=safe_t_inc,
            RefState="DEF",
        )
    )
    if df.empty:
        return df

    def get_num(col_name):
        if col_name is None or col_name not in df.columns:
            return pd.Series([None] * len(df))
        return pd.to_numeric(df[col_name], errors="coerce")

    out = pd.DataFrame(
        {
            "T": get_num(_find_col(df, "temp")),
            "rho": get_num(_find_col(df, "density")),
            "hs": get_num(_find_col(df, "enthalpy")) * 1000.0,
            "Cp": get_num(_find_col(df, "cp")) * 1000.0,
            "Cv": get_num(_find_col(df, "cv")) * 1000.0,
            "mu": get_num(_find_col(df, "visc")),
            "kappa": get_num(_find_col(df, "cond")),
        }
    )
    return out.dropna(subset=["T"]).reset_index(drop=True)


def fetch_saturation_table(species_id, p_low, p_high, p_inc):
    df = _get_table(
        dict(
            ID=species_id,
            Type="SatP",
            RefState="DEF",
            PLow=p_low,
            PHigh=p_high,
            PInc=p_inc,
        )
    )
    if df.empty:
        raise RuntimeError(
            "NIST WebBook returned no saturation data for this range."
            "Run with --debug to see HTTP status / response snippet."
        )

    def get_num(col_name):
        if col_name is None or col_name not in df.columns:
            return pd.Series([None] * len(df))
        return pd.to_numeric(df[col_name], errors="coerce")

    out = pd.DataFrame(
        {
            "T": get_num(_find_col(df, "temp")),
            "P": get_num(_find_col(df, "press")),
            "rho_l": get_num(_find_phase_col(df, "liquid", "density")),
            "rho_v": get_num(_find_phase_col(df, "vapor", "density")),
            "hs_l": get_num(_find_phase_col(df, "liquid", "enthalpy")) * 1000.0,
            "hs_v": get_num(_find_phase_col(df, "vapor", "enthalpy")) * 1000.0,
            "Cp_l": get_num(_find_phase_col(df, "liquid", "cp")) * 1000.0,
            "Cp_v": get_num(_find_phase_col(df, "vapor", "cp")) * 1000.0,
            "Cv_l": get_num(_find_phase_col(df, "liquid", "cv")) * 1000.0,
            "Cv_v": get_num(_find_phase_col(df, "vapor", "cv")) * 1000.0,
            "mu_l": get_num(_find_phase_col(df, "liquid", "visc")),
            "mu_v": get_num(_find_phase_col(df, "vapor", "visc")),
            "kappa_l": get_num(_find_phase_col(df, "liquid", "cond")),
            "kappa_v": get_num(_find_phase_col(df, "vapor", "cond")),
            "sigma": get_num(_find_col(df, "surf")),
        }
    )
    return out.dropna(subset=["T", "P"]).reset_index(drop=True)


def build_liquid_row(species_id, P, T_grid, sat_row):
    Tsat = sat_row["T"]
    real_T = [t for t in T_grid if t <= Tsat]
    values = {p: {} for p in PROPS}

    if real_T:
        t_inc = (real_T[1] - real_T[0]) if len(real_T) > 1 else 1.0
        real = fetch_isobaric(species_id, P, real_T[0], real_T[-1], t_inc)
        for _, r in real.iterrows():
            for p in PROPS:
                values[p][round(r["T"], 6)] = r[p]

    anchor = {
        "rho": sat_row["rho_l"],
        "hs": sat_row["hs_l"],
        "Cp": sat_row["Cp_l"],
        "Cv": sat_row["Cv_l"],
        "mu": sat_row["mu_l"],
        "kappa": sat_row["kappa_l"],
    }

    floor_T = T_grid[-1]
    floor_val = sat_row["rho_v"]

    rows = {p: [] for p in PROPS}
    for t in T_grid:
        val_found = False
        if t <= Tsat:
            for k in values["rho"]:
                if abs(k - t) < 1e-3:
                    for p in PROPS:
                        rows[p].append(values[p][k])
                    val_found = True
                    break

        if not val_found:
            for p in PROPS:
                if p == "rho":
                    if t >= floor_T:
                        rows[p].append(floor_val)
                    else:
                        frac = (
                            (t - Tsat) / (floor_T - Tsat)
                            if floor_T != Tsat
                            else 1.0
                        )
                        rows[p].append(
                            anchor["rho"] + frac * (floor_val - anchor["rho"])
                        )
                else:
                    rows[p].append(anchor[p])
    return rows


def build_vapour_row(species_id, P, T_grid, sat_row):
    Tsat = sat_row["T"]
    real_T = [t for t in T_grid if t >= Tsat]
    values = {p: {} for p in PROPS}

    if real_T:
        t_inc = (real_T[1] - real_T[0]) if len(real_T) > 1 else 1.0
        real = fetch_isobaric(species_id, P, real_T[0], real_T[-1], t_inc)
        for _, r in real.iterrows():
            for p in PROPS:
                values[p][round(r["T"], 6)] = r[p]

    anchor = {
        "rho": sat_row["rho_v"],
        "hs": sat_row["hs_v"],
        "Cp": sat_row["Cp_v"],
        "Cv": sat_row["Cv_v"],
        "mu": sat_row["mu_v"],
        "kappa": sat_row["kappa_v"],
    }

    rows = {p: [] for p in PROPS}
    for t in T_grid:
        val_found = False
        if t >= Tsat:
            for k in values["rho"]:
                if abs(k - t) < 1e-3:
                    for p in PROPS:
                        rows[p].append(values[p][k])
                    val_found = True
                    break

        if not val_found:
            for p in PROPS:
                if p == "rho":
                    rows[p].append(anchor["rho"] * Tsat / t)
                else:
                    rows[p].append(anchor[p])
    return rows


def format_row(vals):
    return "            ( " + " ".join(f"{v:.10g}" for v in vals) + " )"


def write_foam_table(path, p_grid, t_grid, table):
    p_low_pa = (
        int(round(p_grid[0] * 1e6))
        if (p_grid[0] * 1e6).is_integer()
        else p_grid[0] * 1e6
    )
    p_high_pa = (
        int(round(p_grid[-1] * 1e6))
        if (p_grid[-1] * 1e6).is_integer()
        else p_grid[-1] * 1e6
    )
    t_low = (
        int(round(t_grid[0])) if float(t_grid[0]).is_integer() else t_grid[0]
    )
    t_high = (
        int(round(t_grid[-1])) if float(t_grid[-1]).is_integer() else t_grid[-1]
    )
    nr, nc = len(p_grid), len(t_grid)

    def block(name, rows):
        lines = [
            f"    {name}",
            "    {",
            f"        low             ( {p_low_pa} {t_low} );",
            f"        high            ( {p_high_pa} {t_high} );",
            f"        values          {nr} {nc}",
            "        (",
        ]
        for r in rows:
            lines.append(format_row(r))
        lines += ["        );", "    }"]
        return "\n".join(lines)

    with open(path, "w") as f:
        f.write("/*--------------------------------*- C++ -*----------------------------------*\\\n")
        f.write("  =========                 |\n")
        f.write("  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n")
        f.write("   \\\\    /   O peration     | Website:  https://openfoam.org\n")
        f.write("    \\\\  /    A nd           | Version:  dev\n")
        f.write("     \\\\/     M anipulation  |\n")
        f.write("\\*---------------------------------------------------------------------------*/\n\n")

        f.write("equationOfState\n{\n")
        f.write(block("rho", table["rho"]) + "\n}\n\n")

        f.write("thermodynamics\n{\n")
        f.write("    hf              0;\n    sf              0;\n")
        f.write(block("hs", table["hs"]) + "\n")
        f.write(block("Cp", table["Cp"]) + "\n")
        f.write(block("Cv", table["Cv"]) + "\n}\n\n")

        f.write("transport\n{\n")
        f.write(block("mu", table["mu"]) + "\n")
        f.write(block("kappa", table["kappa"]) + "\n}\n\n\n")
        f.write("// ************************************************************************* //\n")


def main():
    global DEBUG

    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--species")
    ap.add_argument("--p-low", type=float)
    ap.add_argument("--p-high", type=float)
    ap.add_argument("--p-inc", type=float)
    ap.add_argument("--t-low", type=float)
    ap.add_argument("--t-high", type=float)
    ap.add_argument("--t-inc", type=float)
    ap.add_argument("--debug", action="store_true", help=("Print HTTP status codes and response diagnostics for each NIST WebBook request."))

    args = ap.parse_args()
    DEBUG = args.debug

    species = (
        args.species or input("Species of interest (e.g. R134a): ").strip()
    )
    p_low = (
        args.p_low
        if args.p_low is not None
        else float(input("Lowest pressure [MPa]: "))
    )
    p_high = (
        args.p_high
        if args.p_high is not None
        else float(input("Highest pressure [MPa]: "))
    )
    p_inc = (
        args.p_inc
        if args.p_inc is not None
        else float(input("Pressure increment [MPa]: "))
    )
    t_low = (
        args.t_low
        if args.t_low is not None
        else float(input("Lowest temperature [K]: "))
    )
    t_high = (
        args.t_high
        if args.t_high is not None
        else float(input("Highest temperature [K]: "))
    )
    t_inc = (
        args.t_inc
        if args.t_inc is not None
        else float(input("Temperature increment [K]: "))
    )

    species_key = species.strip().lower()
    if species_key not in CAS_LOOKUP:
        sys.exit(f"Unknown species '{species}'.")
    species_id = "C" + CAS_LOOKUP[species_key].replace("-", "")
    print(f"Using NIST WebBook ID: {species_id}")

    p_grid = frange(p_low, p_high, p_inc)
    t_grid = frange(t_low, t_high, t_inc)
    print(
        f"Pressure grid: {len(p_grid)} points, Temperature grid:"
        f" {len(t_grid)} points"
    )

    print("Fetching saturation table...")
    sat = fetch_saturation_table(species_id, p_low, p_high, p_inc)
    sat = sat.sort_values("P").reset_index(drop=True)

    sat_out = sat.rename(
        columns={
            "T": "Temperature [K]",
            "P": "Pressure [MPa]",
            "rho_l": "Liquid Density [kg/m3]",
            "mu_l": "Liquid Viscosity [Pa s]",
            "kappa_l": "Liquid Thermal Conductivity [W/m K]",
            "sigma": "Surface Tension [N/m]",
            "rho_v": "Vapour Density [kg/m3]",
            "mu_v": "Vapour Viscosity [Pa s]",
            "kappa_v": "Vapour Thermal Conductivity [W/m K]",
        }
    )
    keep = [
        "Temperature [K]",
        "Pressure [MPa]",
        "Liquid Density [kg/m3]",
        "Liquid Viscosity [Pa s]",
        "Liquid Thermal Conductivity [W/m K]",
        "Surface Tension [N/m]",
        "Vapour Density [kg/m3]",
        "Vapour Viscosity [Pa s]",
        "Vapour Thermal Conductivity [W/m K]",
    ]
    sat_out[[c for c in keep if c in sat_out.columns]].to_csv("saturation.csv", index=False)

    liquid_table = {p: [] for p in PROPS}
    vapour_table = {p: [] for p in PROPS}

    for P in p_grid:
        srow = sat.iloc[(sat["P"] - P).abs().idxmin()]
        print(
            f"P = {P:g} MPa -> Tsat = {srow['T']:.3f} K ... fetching liquid/vapour branches"
        )

        lrows = build_liquid_row(
            species_id,
            P,
            t_grid,
            srow,
        )
        vrows = build_vapour_row(species_id, P, t_grid, srow)
        for p in PROPS:
            liquid_table[p].append(lrows[p])
            vapour_table[p].append(vrows[p])

    write_foam_table("liquid", p_grid, t_grid, liquid_table)
    write_foam_table("vapour", p_grid, t_grid, vapour_table)
    print("Files written successfully.")


if __name__ == "__main__":
    main()
