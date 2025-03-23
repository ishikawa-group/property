# Property
* A Python package for calculating band gap and hydration properties of crystalline materials.

## Features
* Band gap calculation
  - Calculate using either VASP or M3GNet
  - M3GNet requires no additional setup
  - VASP requires environment variable setup

* Hydration energy calculation
  - Calculate using VASP
  - Oxygen vacancy formation energy
  - Proton formation energy
  - Hydration energy

## Requirements
* Python ~= 3.12
* ase ~= 3.24
* pymatgen == 2025.3.10
* matgl == 1.2.1

## Installation

```bash
pip install .
```

## Usage

### Band gap calculation

```python
from ase.io import read
from property.bandgap import get_bandgap

# Load structure
atoms = read("structure.cif")

# Calculate using M3GNet
bandgap = get_bandgap(atoms=atoms, calculator="m3gnet")
print(f"bandgap = {bandgap:5.3f} eV")

# Calculate using VASP (requires environment setup)
bandgap = get_bandgap(atoms=atoms, calculator="vasp")
print(f"bandgap = {bandgap:5.3f} eV")
```

### Hydration energy calculation

```python
from ase.io import read
from property.hydration import get_hydration_energy

# Load structure
atoms = read("structure.cif")

# Calculate hydration energy, oxygen vacancy formation energy, and proton formation energy
hydration, vacancy, oh = get_hydration_energy(atoms)
print(f"hydration energy = {hydration:5.3f} eV")
print(f"oxygen vacancy formation energy = {vacancy:5.3f} eV")
print(f"proton formation energy = {oh:5.3f} eV")
```

## VASP environment setup

When using VASP, the following environment variables must be set:

1. VASP execution command (one of the following)
   - `ASE_VASP_COMMAND`
   - `VASP_COMMAND`
   - `VASP_SCRIPT`

2. Pseudopotential path
   - `VASP_PP_PATH`

## License

Not specified
