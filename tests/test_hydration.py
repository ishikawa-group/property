import sys
sys.path.append("../")
from ase.io import read
from property.hydration import get_hydration_energy
import matgl
from matgl.ext.ase import PESCalculator

cif = "BaZrO3.cif"  # should be ~3.3 eV
atoms = read(cif)
atoms *= [2, 2, 2]

pot = matgl.load_model("M3GNet-MP-2021.2.8-PES")
calc = PESCalculator(pot)
atoms.calc = calc

hydration, vacancy, oh = get_hydration_energy(atoms)
print(f"hydration energy = {hydration:5.3f} eV")
print(f"oxygen vacation formation energy = {vacancy:5.3f} eV")
print(f"proton formation energy = {oh:5.3f} eV")
