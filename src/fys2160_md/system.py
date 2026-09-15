"""Physical particle state, independent of integration and experiment settings."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .bonds import resolve_model, bonds_for_species, coefficients

MODELS = {'atomic': 0, 'diatomic-flexible': 1, 'diatomic-rigid': 2}


def _positive(value, name, *, zero=False):
    value = float(value)
    if not np.isfinite(value) or (value < 0 if zero else value <= 0):
        raise ValueError(f'{name} must be finite and {"nonnegative" if zero else "positive"}.')
    return value


def _integer(value, name, minimum=1):
    if isinstance(value, (bool, np.bool_)) or int(value) != value or value < minimum:
        raise ValueError(f'{name} must be an integer >= {minimum}.')
    return int(value)


def _readonly(array):
    if array is None:
        return None
    result = array.copy()
    result.flags.writeable = False
    return result


@dataclass(frozen=True)
class Box:
    lengths: np.ndarray

    @property
    def volume(self):
        return float(np.prod(self.lengths))


@dataclass(frozen=True)
class Atoms:
    positions: np.ndarray
    velocities: np.ndarray | None
    masses: np.ndarray
    molecule_ids: np.ndarray
    species: np.ndarray

    def __len__(self):
        return len(self.positions)


class System:
    """Periodic box and particle state. Generated systems have no velocities yet.

    MDSimulation copies this state, initializes missing velocities, and evolves
    its own copy. All public array accessors return read-only snapshots.
    """

    def __init__(self, positions, box, *, velocities=None, masses=None, molecule=None, bond=None, model=None, species=None):
        model, bond = resolve_model(model, molecule, bond)
        self._bond = bond
        x=np.array(positions,dtype=float,order='C',copy=True)
        lengths=np.broadcast_to(np.asarray(box,dtype=float),(3,)).copy()
        if x.ndim!=2 or x.shape[1]!=3 or not np.isfinite(x).all():
            raise ValueError('positions must be a finite (N,3) array.')
        if not np.isfinite(lengths).all() or np.any(lengths<=0):
            raise ValueError('box lengths must be finite and positive.')
        if model not in MODELS:
            raise ValueError(f'model must be one of {tuple(MODELS)}')
        if model!='atomic' and len(x)%2:
            raise ValueError('Diatomic models require an even number of atoms.')
        labels = np.asarray('A' if species is None else species)
        if labels.ndim == 0: labels = np.repeat(labels, len(x))
        if labels.shape != (len(x),) or labels.dtype.kind not in 'US':
            raise ValueError('species must be a name or one string label per atom.')
        if any(not str(name).strip() or str(name) != str(name).strip() for name in labels):
            raise ValueError('Species names must be nonempty strings without surrounding spaces.')
        self._species = labels.astype(str, copy=True)
        if isinstance(bond,dict) and np.any(self._species[::2] != self._species[1::2]):
            raise ValueError('Per-species bonds require homonuclear pairs; use a shared bond for heteronuclear molecules.')
        self._bond_specs = () if model == 'atomic' else bonds_for_species(bond, self._species[::2])
        self._bond_parameters = np.array([coefficients(b) for b in self._bond_specs],dtype=float).reshape(-1,4)
        if len(self._bond_parameters) and np.any(self._bond_parameters[:,0]>=min(lengths)/2):
            raise ValueError('Bond length must be less than half the shortest box side.')
        self._model=model
        self._N=_integer(len(x) if model=='atomic' else len(x)//2,'N',2)
        self._box=lengths
        self._x=np.ascontiguousarray(x%lengths)
        self._mass=np.ones(len(x)) if masses is None else np.array(masses,dtype=float,order='C',copy=True)
        if self._mass.shape!=(len(x),) or not np.isfinite(self._mass).all() or np.any(self._mass<=0):
            raise ValueError('masses must be a positive finite (N,) array.')

        self._v=None if velocities is None else np.array(velocities,dtype=float,order='C',copy=True)
        if self._v is not None and (self._v.shape!=x.shape or not np.isfinite(self._v).all()):
            raise ValueError('velocities must be a finite (N,3) array.')
        if model=='diatomic-rigid':
            bond=self._x[::2]-self._x[1::2]
            bond-=lengths*np.rint(bond/lengths)
            if not np.allclose(np.linalg.norm(bond,axis=1),self._bond_parameters[:,0],atol=1e-10,rtol=0):
                raise ValueError('Rigid bonds must have their specified lengths.')
            if self._v is not None and np.max(np.abs((bond*(self._v[::2]-self._v[1::2])).sum(axis=1)))>1e-10:
                raise ValueError('Rigid bond relative velocities must be perpendicular to the bonds.')

    @property
    def species(self):
        return tuple(dict.fromkeys(self._species.tolist()))

    @property
    def bond(self):
        return self._bond.copy() if isinstance(self._bond,dict) else self._bond

    @property
    def molecule(self):
        return None if self.model == 'atomic' else 'diatomic'

    @property
    def bonds(self):
        """Connected atom-index pairs (zero-based); empty for an atomic system."""
        pairs = np.empty((0, 2), dtype=int) if self.bond is None else np.arange(2*self.N).reshape(-1, 2)
        return _readonly(pairs)

    @property
    def model(self):
        return self._model

    @property
    def N(self):
        return self._N

    @property
    def rho(self):
        return self.N/self.box.volume

    @property
    def box(self):
        return Box(_readonly(self._box))

    @property
    def atoms(self):
        ids=np.arange(len(self._x))//(1 if self.model=='atomic' else 2)
        return Atoms(_readonly(self._x),_readonly(self._v),_readonly(self._mass),_readonly(ids),_readonly(self._species))

    @property
    def dof(self):
        return 3*len(self._x)-3-(self.N if self.model=='diatomic-rigid' else 0)

    def copy(self):
        """Independent copy of geometry and any supplied velocities."""
        return System(self._x,self._box,velocities=self._v,masses=self._mass,model=self.model,bond=self.bond,species=self._species)

    @classmethod
    def from_arrays(cls, positions, velocities, box, *, masses=None, molecule=None, bond=None, model=None, species=None):
        return cls(positions,box,velocities=velocities,masses=masses,model=model,molecule=molecule,bond=bond,species=species)

    def __repr__(self):
        return f'System(model={self.model!r}, N={self.N}, rho={self.rho:g}, velocities={self._v is not None})'
