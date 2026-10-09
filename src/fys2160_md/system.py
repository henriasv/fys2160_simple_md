"""Physical particle state, independent of integration and experiment settings."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .bonds import resolve_model, bonds_for_species, coefficients, RigidWater
from .units import units as _units

MODELS = {'atomic': 0, 'diatomic-flexible': 1, 'diatomic-rigid': 2, 'water-rigid': 3}
_SIZE = {'atomic': 1, 'diatomic-flexible': 2, 'diatomic-rigid': 2, 'water-rigid': 3}


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
    """The state of a system: box, atoms, optional velocities and a clock.

    ``MDSimulation.run(system, ...)`` evolves a copy and never changes this one;
    the run's final state is ``run.state``. A state from a run also carries the
    simulation time, step count and thermostat/barostat state, so continuing from
    it is identical to one longer run. All public arrays are read-only snapshots.
    """

    def __init__(self, positions, box, *, velocities=None, masses=None, molecule=None, bond=None, model=None, species=None, boundary="periodic", units="reduced",
                 unwrapped=None, time=0., step=0):
        if boundary not in ("periodic", "open"):
            raise ValueError("boundary must be periodic or open.")
        self._boundary = boundary
        self._units = _units(units).name
        model, bond = resolve_model(model, molecule, bond)
        if boundary == "open" and model != "atomic":
            raise ValueError("Open boundaries currently support unbonded particles only.")
        self._bond = bond
        x=np.array(positions,dtype=float,order='C',copy=True)
        lengths=np.broadcast_to(np.asarray(box,dtype=float),(3,)).copy()
        if x.ndim!=2 or x.shape[1]!=3 or not np.isfinite(x).all():
            raise ValueError('positions must be a finite (N,3) array.')
        if not np.isfinite(lengths).all() or np.any(lengths<=0):
            raise ValueError('box lengths must be finite and positive.')
        if model not in MODELS:
            raise ValueError(f'model must be one of {tuple(MODELS)}')
        if model.startswith('diatomic') and len(x)%2:
            raise ValueError('Diatomic models require an even number of atoms.')
        if model=='water-rigid' and len(x)%3:
            raise ValueError('Water requires O, H, H atom triplets.')
        labels = np.asarray('A' if species is None else species)
        if labels.ndim == 0: labels = np.repeat(labels, len(x))
        if labels.shape != (len(x),) or labels.dtype.kind not in 'US':
            raise ValueError('species must be a name or one string label per atom.')
        if any(not str(name).strip() or str(name) != str(name).strip() for name in labels):
            raise ValueError('Species names must be nonempty strings without surrounding spaces.')
        self._species = labels.astype(str, copy=True)
        if isinstance(bond,dict) and np.any(self._species[::2] != self._species[1::2]):
            raise ValueError('Per-species bonds require homonuclear pairs; use a shared bond for heteronuclear molecules.')
        if model=='water-rigid' and (np.any(self._species[1::3] != self._species[2::3]) or np.any(self._species[0::3] == self._species[1::3])):
            raise ValueError('Each water molecule must list O first, then two H atoms of one other species.')
        self._bond_specs = () if model in ('atomic', 'water-rigid') else bonds_for_species(bond, self._species[::2])
        self._bond_parameters = np.array([coefficients(b) for b in self._bond_specs],dtype=float).reshape(-1,4)
        if len(self._bond_parameters) and np.any(self._bond_parameters[:,0]>=min(lengths)/2):
            raise ValueError('Bond length must be less than half the shortest box side.')
        self._model=model
        self._N=_integer(len(x)//_SIZE[model],'N',2)
        self._box=lengths
        self._x=np.ascontiguousarray(x%lengths if boundary == "periodic" else x)
        self._mass=np.ones(len(x)) if masses is None else np.array(masses,dtype=float,order='C',copy=True)
        if self._mass.shape!=(len(x),) or not np.isfinite(self._mass).all() or np.any(self._mass<=0):
            raise ValueError('masses must be a positive finite (N,) array.')

        self._v=None if velocities is None else np.array(velocities,dtype=float,order='C',copy=True)
        if self._v is not None and (self._v.shape!=x.shape or not np.isfinite(self._v).all()):
            raise ValueError('velocities must be a finite (N,3) array.')
        self._u=None if unwrapped is None else np.array(unwrapped,dtype=float,order='C',copy=True)
        if self._u is not None and (self._u.shape!=x.shape or not np.isfinite(self._u).all()):
            raise ValueError('unwrapped positions must be a finite (N,3) array.')
        self._time=float(time)
        if not np.isfinite(self._time):
            raise ValueError('time must be finite.')
        self._step=_integer(step,'step',0)
        # Thermostat random-number state, barostat rate and PPPM mesh choice
        # (state array, kspace_accuracy, coulomb_cutoff), carried between runs.
        self._rng=None
        self._rate=0.
        self._pppm=None
        self._heat=0.  # total energy added by heat_rate along the runs leading here
        if model=='diatomic-rigid':
            bond=self._x[::2]-self._x[1::2]
            bond-=lengths*np.rint(bond/lengths)
            if not np.allclose(np.linalg.norm(bond,axis=1),self._bond_parameters[:,0],atol=1e-10,rtol=0):
                raise ValueError('Rigid bonds must have their specified lengths.')
            if self._v is not None and np.max(np.abs((bond*(self._v[::2]-self._v[1::2])).sum(axis=1)))>1e-10:
                raise ValueError('Rigid bond relative velocities must be perpendicular to the bonds.')
        if model=='water-rigid':
            for (a, b), length in zip(((0, 1), (0, 2), (1, 2)), self._bond.constraints):
                if length >= min(lengths)/2:
                    raise ValueError('Water must be smaller than half the shortest box side.')
                dr=self._x[a::3]-self._x[b::3]
                dr-=lengths*np.rint(dr/lengths)
                if not np.allclose(np.linalg.norm(dr,axis=1),length,atol=1e-9,rtol=0):
                    raise ValueError('Water molecules must have the RigidWater O-H lengths and H-O-H angle.')
                if self._v is not None and np.max(np.abs((dr*(self._v[a::3]-self._v[b::3])).sum(axis=1)))>1e-9:
                    raise ValueError('Rigid-water relative velocities must be perpendicular to the constrained distances.')

    @property
    def boundary(self):
        """Open systems use box only as a viewing frame, without walls or images."""
        return self._boundary

    @property
    def time(self):
        """Simulation time of this state, in the system's time unit."""
        return self._time

    @property
    def step(self):
        """Number of MD steps taken to reach this state."""
        return self._step

    @property
    def units(self):
        """Unit system: "reduced" Lennard-Jones units or LAMMPS "real" units."""
        return self._units

    @property
    def species(self):
        return tuple(dict.fromkeys(self._species.tolist()))

    @property
    def bond(self):
        return self._bond.copy() if isinstance(self._bond,dict) else self._bond

    @property
    def molecule(self):
        return {'atomic': None, 'water-rigid': 'water'}.get(self.model, 'diatomic')

    @property
    def bonds(self):
        """Connected atom-index pairs (zero-based); empty for an atomic system."""
        if self.model == 'water-rigid':
            o = 3*np.arange(self.N)
            pairs = np.stack((np.stack((o, o+1), axis=1), np.stack((o, o+2), axis=1)), axis=1).reshape(-1, 2)
        else:
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
        return self.N/self.box.volume if self.boundary == "periodic" else float("nan")

    @property
    def box(self):
        return Box(_readonly(self._box))

    @property
    def atoms(self):
        ids=np.arange(len(self._x))//_SIZE[self.model]
        return Atoms(_readonly(self._x),_readonly(self._v),_readonly(self._mass),_readonly(ids),_readonly(self._species))

    @property
    def dof(self):
        return 3*len(self._x)-3-{'diatomic-rigid': self.N, 'water-rigid': 3*self.N}.get(self.model, 0)

    def copy(self):
        """Independent copy of the complete state."""
        other=System(self._x,self._box,velocities=self._v,masses=self._mass,model=self.model,bond=self.bond,species=self._species,boundary=self.boundary,units=self.units,
                     unwrapped=self._u,time=self._time,step=self._step)
        other._rng,other._rate,other._heat=self._rng,self._rate,self._heat
        other._pppm=None if self._pppm is None else (self._pppm[0].copy(),)+tuple(self._pppm[1:])
        return other

    @classmethod
    def from_arrays(cls, positions, velocities, box, *, masses=None, molecule=None, bond=None, model=None, species=None, boundary="periodic", units="reduced"):
        return cls(positions,box,velocities=velocities,masses=masses,model=model,molecule=molecule,bond=bond,species=species,boundary=boundary,units=units)

    def __repr__(self):
        return f'System(model={self.model!r}, N={self.N}, rho={self.rho:g}, velocities={self._v is not None}, time={self._time:g}, step={self._step})'
