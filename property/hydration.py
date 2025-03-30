import os
import copy
from ase.io import read
from ase import Atom
from ase.io.vasp import read_vasp_out
from ase.build import molecule
from ase.calculators.vasp import Vasp
import matgl
from matgl.ext.ase import PESCalculator
import warnings
warnings.simplefilter("ignore")  # To suppress warnings for clearer output


def set_vasp_calculator(atoms=None, directory=None):
    """
    Set VASP calculator for Atoms.

    Args:
        atoms: ASE Atoms object.
    """

    calc = Vasp(prec="normal", xc="pbe", ispin=2, lorbit=10,
                encut=520, ediff=1e-6, algo="Normal", nelm=40, nelmin=5,
                ibrion=2, nsw=50, isif=2, potim=0.1, ediffg=-0.05,
                ismear=0, sigma=0.05,
                kpts=[4, 4, 4], kgamma=True,
                npar=4, nsim=4,
                directory=directory,
                lwave=False, lcharg=False,
                lreal=False,
                )
    atoms.calc = calc
    return None


def get_molecule_energy(atoms=None) -> float:
    """
    """
    mol = copy.deepcopy(atoms)
    mol.cell = [20]*3
    mol.pbc = True
    mol.center()
    formula = mol.get_chemical_formula()

    if "vasp" in mol.calc.name:
        ispin = 2 if formula == "O2" else 1
        directory = "tmpdir_" + formula
        mol.calc = Vasp(prec="normal", xc="pbe", encut=400, ispin=ispin, ibrion=2,
                        nsw=50, ismear=0, kpts=[1, 1, 1], directory=directory)
    else:
        pass

    energy = mol.get_potential_energy()

    return energy


def get_hydration_energy(atoms=None) -> float:
    """
    Calculates and prints defect formation and hydration energies with ZPE corrections.

    Args:
        atoms: ASE Atoms object.
    Returns:
        delta_E_hydr (float): Hydration energy
        delta_E_VOx (float): Oxygen vacancy formation energy
        delta_E_OHx (float): Proton formation energy
    """
    # constant energy values for H2O, H2, O2 (eV)
    use_reference = False
    use_vasp = False

    if use_reference:
        E_H2O = -14.891
        E_H2 = -6.715
        E_O2
    else:
        h2  = molecule("H2")
        h2o = molecule("H2O")

        if use_vasp:
            h2.calc = Vasp()
            h2o.calc = Vasp()
        else:
            pot = matgl.load_model("M3GNet-MP-2021.2.8-PES")
            calc = PESCalculator(pot)
            h2.calc = calc
            h2o.calc = calc

        E_H2  = get_molecule_energy(atoms=h2)
        E_H2O = get_molecule_energy(atoms=h2o)

    # constant ZPE values (eV)
    ZPE_H2O = 0.56
    ZPE_H2 = 0.27
    ZPE_O2 = 0.10

    # apply ZPE corrections
    E_H2_corr = E_H2 + ZPE_H2
    E_H2O_corr = E_H2O + ZPE_H2O

    # H2O formation enthalpy (experimental value, 241.8 (gas) or 285.8 (liquid) kJ/mol -> eV)
    E_exp_HOF = -2.507  # -2.962

    E_O2_corr = 2 * ((E_H2O + ZPE_H2O) - (E_H2 + ZPE_H2) - E_exp_HOF) - ZPE_O2

    # make the hydrated structure
    atoms_with_Vox = copy.deepcopy(atoms)

    # Step1: Remove an oxygen atom to simulate the vacancy(Vox)

    # Find the oxygen atom 
    oxygen_indices = [atom.index for atom in atoms_with_Vox if atom.symbol == "O"]
    # Record the oxygen atom position
    oxygen_position = atoms[oxygen_indices[0]].position

    # Step2: Add hydrogen atom to original lattice, to simulate proton defect (OHx)
    atoms_with_OHx = copy.deepcopy(atoms)
    atoms_with_OHx.append(Atom("H", position=oxygen_position + [0.6, 0.0, 0.6]))

    if "vasp" in atoms.calc.name:
        tmpdir = "tmpdir_hydration"
        tmpdir_vox = "tmpdir_hydration_vox"
        tmpdir_ohx = "tmpdir_hydration_ohx"

        set_vasp_calculator(atoms, directory=tmpdir)
        atoms.calc.set(isif=8)
        set_vasp_calculator(atoms_with_Vox, directory=tmpdir_vox)
        set_vasp_calculator(atoms_with_OHx, directory=tmpdir_ohx)
    else:
        pass

    E_pristine = atoms.get_potential_energy()
    E_Vox = atoms_with_Vox.get_potential_energy()
    E_OHx = atoms_with_OHx.get_potential_energy()

    # calculate chemical potential
    mu_O = 0.5*E_O2_corr
    mu_H2O = E_H2O_corr
    mu_H = 0.5*(mu_H2O - mu_O)

    # calculate hydration energy
    delta_E_Vox = E_Vox - (E_pristine - mu_O)
    delta_E_OHx = E_OHx - (E_pristine + mu_H)
    delta_E_hydr = 2*E_OHx - (E_pristine + E_Vox + mu_H2O)

    # return results
    return delta_E_hydr, delta_E_Vox, delta_E_OHx
