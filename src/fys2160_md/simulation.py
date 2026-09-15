"""Run dynamics on an explicitly constructed System."""
from .system import System


class Simulation:
    """Course solver entry point; geometry and evolving state belong to System."""

    @staticmethod
    def run(system, timestep=None, steps=1000, **options):
        """Advance system in place and return a saved Run.

        Omitted physical settings retain their values on this system. Use
        storage_name to choose the output folder; reusing a name overwrites it.
        """
        if not isinstance(system,System):
            raise TypeError('Pass a System, such as FCC(rho=0.1, N=500).')
        return system._run(timestep=timestep,steps=steps,**options)
