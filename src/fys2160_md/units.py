"""Unit systems: reduced Lennard-Jones units (default) and LAMMPS ``real`` units."""
from dataclasses import dataclass


@dataclass(frozen=True)
class Units:
    name: str
    boltzmann: float  # energy per temperature unit
    mvv2e: float      # mass * velocity**2 -> energy
    qqr2e: float      # Coulomb constant: energy * length / charge**2
    pressure: float   # pressure unit per energy/length**3
    description: str
    time: str
    temperature: str
    pressure_unit: str
    energy: str


REDUCED = Units('reduced', 1., 1., 1., 1., 'fixed reduced reference units: length=energy=mass=kB=1',
                't*', 'T*', 'P*', 'reduced energy')
# Constants as in LAMMPS units real: Å, fs, g/mol, kcal/mol, K, atm, e.
REAL = Units('real', 0.0019872067, 48.88821291**2, 332.06371, 68568.415,
             'LAMMPS real units: length Å, time fs, mass g/mol, energy kcal/mol, temperature K, pressure atm, charge e',
             'fs', 'K', 'atm', 'kcal/mol')
UNITS = {'reduced': REDUCED, 'real': REAL}
AVOGADRO_PER_CM3_A3 = 0.602214076  # (g/mol)/Å^3 -> g/cm^3 divisor


def units(name):
    if name not in UNITS:
        raise ValueError('units must be "reduced" or "real".')
    return UNITS[name]
