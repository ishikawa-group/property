import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator
from matgl.models import MEGNetCalculator as MatGLMLPotential


class CustomMatGLCalculator(Calculator):
    """
    A custom ASE calculator for MatGL machine learning potential

    Implements the standard ASE calculator interface with MatGL
    """

    implemented_properties = ['energy', 'forces', 'stress']

    def __init__(self, model=None, **kwargs):
        """
        Initialize the MatGL calculator

        Parameters:
        -----------
        model : MEGNetCalculator, optional
            Pre-trained MatGL model. If None, uses default pretrained model.
        **kwargs : dict
            Additional parameters to pass to the base Calculator
        """
        super().__init__(**kwargs)

        # Use default pretrained model if not provided
        if model is None:
            model = MatGLMLPotential.load()

        self.model = model

    def calculate(self, atoms=None, properties=['energy'], system_changes=None):
        """
        Calculate energy, forces, and stress for the given atomic configuration

        Parameters:
        -----------
        atoms : Atoms object
            The atomic configuration to calculate
        properties : list
            Properties to calculate (energy, forces, stress)
        system_changes : list
            List of system changes
        """
        # Always call the parent class's calculate method
        super().calculate(atoms, properties, system_changes)

        # Convert ASE Atoms to MatGL compatible format
        pred = self.model.predict(atoms)

        # Store results
        self.results['energy'] = pred['energy'].numpy()[0]
        self.results['forces'] = pred['forces'].numpy()[0]

        # Calculate stress (optional, might require additional processing)
        # Note: MatGL doesn't directly provide stress, so this is a placeholder
        self.results['stress'] = np.zeros((3, 3))

    def get_potential_energy(self, atoms=None, force_consistent=False):
        """
        Convenience method to get potential energy

        Parameters:
        -----------
        atoms : Atoms object, optional
            Atomic configuration to calculate energy for
        force_consistent : bool
            Flag for force-consistent energy (not typically used with ML potentials)

        Returns:
        --------
        float
            Potential energy of the system
        """
        self.calculate(atoms)
        return self.results['energy']


# Example usage
def example_usage():
    # Create a sample atomic structure
    atoms = Atoms('Cu4',
                  positions=[[0, 0, 0],
                             [0, 2.5, 2.5],
                             [2.5, 0, 2.5],
                             [2.5, 2.5, 0]],
                  cell=[5, 5, 5])

    # Initialize the custom MatGL calculator
    calculator = CustomMatGLCalculator()

    # Attach calculator to atoms
    atoms.calc = calculator

    # Calculate and print energy and forces
    print(f"Potential Energy: {atoms.get_potential_energy()} eV")
    print(f"Forces:\n{atoms.get_forces()}")


# Run the example
if __name__ == "__main__":
    example_usage()
