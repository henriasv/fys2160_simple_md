"""Pair potentials given for each pair of species: the Buckingham potential."""
from collections.abc import Mapping
from dataclasses import dataclass
import math
import numpy as np


@dataclass(frozen=True)
class Buckingham:
    """The Buckingham potential between two species, U(r) = A exp(-r/rho) - C/r^6.

    A is an energy, rho a length and C an energy times length^6, in the units of the
    simulation (eV and Å in ``units="metal"``). The repulsion A exp(-r/rho) stays finite
    as r goes to 0, so the attraction -C/r^6 wins at short distances: the potential
    has a barrier at ``barrier`` and falls without bound inside it. Two atoms that get
    that close (at a high temperature, or with a long timestep) stop the run.
    """
    A: float
    rho: float
    C: float

    def __post_init__(self):
        for name in ('A', 'rho', 'C'):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)) or not math.isfinite(value):
                raise ValueError(f'Buckingham {name} must be a finite number.')
            object.__setattr__(self, name, float(value))
        if self.A < 0 or self.C < 0 or self.rho <= 0:
            raise ValueError('Buckingham needs A >= 0, rho > 0 and C >= 0.')
        if self.C > 0 and self.barrier is None:
            raise ValueError('This Buckingham potential is attractive at every distance: its repulsion A exp(-r/rho) '
                             'never outweighs C/r^6. Check A, rho and C.')

    def energy(self, r):
        """U(r) at the distances r (a number or an array)."""
        r = np.asarray(r, dtype=float)
        return self.A*np.exp(-r/self.rho) - self.C/r**6

    @property
    def barrier(self):
        """The distance of the potential's inner maximum (0 without the C/r^6 term).

        Inside it the force pulls the atoms together. Charges are left out: opposite
        charges attract as well, so the true barrier lies a little further out.
        """
        if self.C == 0:
            return 0.
        # dU/dr = 0 where log(A/rho) - r/rho = log(6C) - 7 log r. The difference h(r)
        # of the two sides rises up to r = 7 rho, so it has at most one root below that.
        h = lambda r: math.log(self.A/self.rho) - r/self.rho - math.log(6*self.C) + 7*math.log(r) if self.A > 0 else -math.inf
        low, high = 1e-6*self.rho, 7*self.rho
        if h(high) <= 0:
            return None
        for _ in range(200):
            middle = .5*(low+high)
            low, high = (middle, high) if h(middle) < 0 else (low, middle)
        return high


def pair_table(pairs):
    """Check {(species, species): Buckingham} and return it with sorted name pairs."""
    if not isinstance(pairs, Mapping) or not pairs:
        raise ValueError('pair_potential must be "lj", "lj96", "gravity" or a mapping from pairs of species to '
                         'Buckingham(A=..., rho=..., C=...), such as {("Si", "O"): Buckingham(...)}.')
    table = {}
    for key, term in pairs.items():
        if (not isinstance(key, tuple) or len(key) != 2 or not all(isinstance(name, str) and name.strip() for name in key)):
            raise ValueError(f'Each key of pair_potential is a pair of species names, such as ("Si", "O"); got {key!r}.')
        if not isinstance(term, Buckingham):
            raise TypeError(f'pair_potential[{key!r}] must be a Buckingham(A=..., rho=..., C=...).')
        name = tuple(sorted(key))
        if name in table:
            raise ValueError(f'The pair {key!r} is given twice (the order of the two species does not matter).')
        table[name] = term
    return table


def buckingham_array(table, names):
    """(types, types, 4) array of A, rho, C and barrier; pairs left out do not interact."""
    index = {name: k for k, name in enumerate(names)}
    unknown = sorted({name for pair in table for name in pair} - set(names))
    if unknown:
        raise ValueError(f'pair_potential names species {unknown} that this state does not have; it has {list(names)}.')
    array = np.zeros((len(names), len(names), 4))
    array[:, :, 1] = 1.  # rho for pairs without a term; A = C = 0 skips them
    for (a, b), term in table.items():
        row = (term.A, term.rho, term.C, term.barrier)
        array[index[a], index[b]] = row
        array[index[b], index[a]] = row
    return np.ascontiguousarray(array)
