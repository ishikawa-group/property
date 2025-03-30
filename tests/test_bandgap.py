import sys
sys.path.append("../")
from ase.io import read
from property.bandgap import get_bandgap
import matgl
from matgl.ext.ase import PESCalculator

cif = "./BaZrO3.cif"   # should be ~3.3 eV

atoms = read(cif)
atoms *= [2, 2, 2]

pot = matgl.load_model("M3GNet-MP-2021.2.8-PES")
calc = PESCalculator(pot)
atoms.calc = calc

bandgap = get_bandgap(atoms=atoms)
print(f"bandgap = {bandgap:5.3f} eV")
