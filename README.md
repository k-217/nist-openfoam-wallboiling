# nist-openfoam-wallboiling

`nist-openfoam-wallboiling.py` is a Python utility that fetches thermophysical and transport property data directly from the [NIST Chemistry WebBook](https://webbook.nist.gov/chemistry/fluid/) and formats it into OpenFOAM-style property tables. It can be directly used for wall boiling or other multiphase simulations (`multiphaseEulerFoam`), where phase properties be evaluated across subcooled, saturated, and superheated regimes.

Reference files for the OpenFOAM expected input are obtained from [OpenFOAM 12 Resources](https://github.com/OpenFOAM/OpenFOAM-12/tree/master/tutorials/resources/thermoData). 

---

## Usage

Clone the repository: 

```bash
git clone https://github.com/k-217/nist-openfoam-wallboiling.git
cd nist-openfoam-wallboiling
```

Install the required dependencies listed in `requirements.txt`:

```bash
pip install -r requirements.txt
```

### Method 1: Python interactive mode

```bash
python nist-openfoam-wallboiling.py
```

Example output:

Species of interest (e.g. R134a): R134a
Lowest pressure [MPa]: 2.2
Highest pressure [MPa]: 3.2
Pressure increment [MPa]: 0.04
Lowest temperature [K]: 325
Highest temperature [K]: 425
Temperature increment [K]: 4
Using NIST WebBook ID: C811972
Pressure grid: 26 points, Temperature grid: 26 points
Fetching saturation table...
P = 2.2 MPa -> Tsat = 344.918 K ... fetching liquid/vapour branches
...
Files written successfully.

### Method 2: Command-line interface

Example:

```bash
python nist-openfoam-wallboiling.py \
    --species R134a \
    --p-low 2.2 \
    --p-high 3.2 \
    --p-inc 0.04 \
    --t-low 325 \
    --t-high 425 \
    --t-inc 4 \
```

Run with ```--debug``` option to verify the HTTP status codes and responses.

## Output

Three files are created ```saturation.csv```, ```liquid```, and ```vapor```. The last two may be used directly for the OpenFOAM simulation. 
