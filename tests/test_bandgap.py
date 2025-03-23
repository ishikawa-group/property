import sys
sys.path.append("../")
from ase.io import read
from property.bandgap import get_bandgap

calculator = "m3gnet"  # ( "vasp" | "m3gnet" )
cif = "./BaZrO3.cif"   # should be ~3.3 eV

atoms = read(cif)
atoms *= [2, 2, 2]

bandgap = get_bandgap(atoms=atoms, calculator=calculator)
print(f"bandgap = {bandgap:5.3f} eV")

