"""Unit systems: reduced Lennard-Jones units (default) and LAMMPS ``real`` and ``metal`` units."""
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
# Constants as in LAMMPS units metal: Å, ps, g/mol, eV, K, bar, e.
METAL = Units('metal', 8.617343e-5, 1.0364269e-4, 14.399645, 1.6021765e6,
              'LAMMPS metal units: length Å, time ps, mass g/mol, energy eV, temperature K, pressure bar, charge e',
              'ps', 'K', 'bar', 'eV')
UNITS = {'reduced': REDUCED, 'real': REAL, 'metal': METAL}
AVOGADRO_PER_CM3_A3 = 0.602214076  # (g/mol)/Å^3 -> g/cm^3 divisor


def units(name):
    if name not in UNITS:
        raise ValueError('units must be "reduced", "real" or "metal".')
    return UNITS[name]
