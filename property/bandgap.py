import os
from ase.calculators.vasp import Vasp
from ase.dft.bandgap import bandgap
from ase.io import read, write
import warnings
import torch
from pymatgen.core import Lattice, Structure
from pymatgen.io.ase import AseAtomsAdaptor
import matgl

warnings.simplefilter("ignore")  # To suppress warnings for clearer output


def get_bandgap(atoms=None, verbose=False) -> float:
    """
    Calculate band gap.

    Args:
        atoms: ASE Atoms object.
        calculator: VASP or M3GNet
        verbose: Verbose printing or not.
    Returns:
        badngap: Calculated band gap value.
    """
    if "vasp" in atoms.calc.name:
        # Check VASP command environment variable is set
        if os.getenv("ASE_VASP_COMMAND"):
            if verbose:
                print("ASE_VASP_COMMAND is set to:", os.getenv("ASE_VASP_COMMAND"))
        elif os.getenv("VASP_COMMAND"):
            if verbose:
                print("VASP_COMMAND is set to:", os.getenv("VASP_COMMAND"))
        elif os.getenv("VASP_SCRIPT"):
            if verbose:
                print("VASP_SCRIPT is set to:", os.getenv("VASP_SCRIPT"))
        else:
            raise EnvironmentError("Error: One of ASE_VASP_COMMAND, VASP_COMMAND, VASP_SCRIPT should be set.")

        # Check if VASP_PP_PATH is set
        if not os.getenv("VASP_PP_PATH"):
            raise EnvironmentError("Error: VASP_PP_PATH is not set.")

        if verbose:
            print("VASP_PP_PATH is set to:", os.getenv("VASP_PP_PATH"))

        # Set up VASP calculator with standard settings
        tmpdir = "tmpdir_bandgap"
        atoms.calc = Vasp(prec="normal", xc="pbe", ispin=2, lorbit=10,
                          ibrion=2, nsw=10, isif=8,
                          encut=520, ediff=1e-6, algo="Normal", nelm=50, nelmin=5,
                          kpts=[4, 4, 4], kgamma=True,
                          ismear=0, sigma=0.05,
                          lwave=False, lcharg=False,
                          npar=4, nsim=4,
                          directory=tmpdir,
                          lreal=False,
                          )

        # Calculate total energy (needed for band structure calculation)
        atoms.get_potential_energy()

        # Calculate band gap
        gap, _, _ = bandgap(atoms.calc, direct=False)

    elif "pescalculator" in atoms.calc.name:  # matgl
        model = matgl.load_model("MEGNet-MP-2019.4.1-BandGap-mfi")
        imethod = 0   # 0: PBE, 1: GLIB-SC, 2: HSE, 3: SCAN
        state_attr = torch.tensor([imethod])
        structure = AseAtomsAdaptor.get_structure(atoms)
        gap = model.predict_structure(structure=structure, state_attr=state_attr)
        gap = float(gap)

    else:
        raise ValueError("Unknown calculator.")

    return gap
