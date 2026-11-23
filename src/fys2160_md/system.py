"""Physical particle state, independent of integration and experiment settings."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .bonds import resolve_model, bonds_for_species, coefficients, RigidWater
from .units import units as _units

MODELS = {'atomic': 0, 'diatomic-flexible': 1, 'diatomic-rigid': 2, 'water-rigid': 3}
_SIZE = {'atomic': 1, 'diatomic-flexible': 2, 'diatomic-rigid': 2, 'water-rigid': 3}
# The six faces of the box, two per axis: at x = 0 and x = Lx, y = 0 and y = Ly,
# z = 0 and z = Lz. x runs to the right, y upwards and z towards the viewer.
WALLS = ('left', 'right', 'bottom', 'top', 'back', 'front')
_BEYOND = 1.  # how far beyond a wall an atom of a given state may be (reduced length)


def _walled(boundary, dimension):
    """Which axes have walls, from one word for all axes or one word per axis."""
    words = (boundary,)*dimension if isinstance(boundary, str) else tuple(boundary)
    if boundary != 'open' and (len(words) != dimension or any(word not in ('periodic', 'walls') for word in words)):
        raise ValueError('boundary must be "periodic", "walls" or "open", or one of "periodic" and "walls" '
                         'for each axis, such as ("walls", "periodic").')
    return tuple(word == 'walls' for word in words)


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


class Box(np.ndarray):
    """Box lengths, one per axis: an array, so ``state.box[0] *= 2`` doubles x."""

    @property
    def lengths(self):
        return self.view(np.ndarray)

    @property
    def volume(self):
        """Product of the lengths: the area of a flat (two-dimensional) box."""
        return float(np.prod(self))


def _check_tilt(lengths, tilt):
    """LAMMPS's limits for a tilted box: each tilt at most half the length it leans along."""
    limits = np.array([lengths[0], lengths[0], lengths[1]])/2
    if np.any(np.abs(tilt) > limits*(1+1e-12)):
        raise ValueError('Each tilt may be at most half the box length it leans along: |xy| and |xz| at most '
                         f'{limits[0]:g} (half of x), |yz| at most {limits[2]:g} (half of y).')


def wrap(x, lengths, tilt):
    """Positions (N, 3) moved by whole box vectors into the box: the parallelepiped with
    edges a = (lx, 0, 0), b = (xy, ly, 0), c = (xz, yz, lz) when tilted."""
    lx, ly, lz = lengths
    xy, xz, yz = tilt
    edges = np.array([[lx, 0, 0], [xy, ly, 0], [xz, yz, lz]])  # rows a, b, c
    fractional = np.linalg.solve(edges.T, np.asarray(x, dtype=float).T).T
    return np.asarray(x, dtype=float)-np.floor(fractional) @ edges


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
    it is identical to one longer run. Its parts can be changed in place, to see
    what happens: positions, velocities, masses, species labels, box, tilt,
    boundary, time and step (``state.box[0] *= 2``, ``state.velocities *= 1.1``).
    The next run checks the result; a run never changes the state it is given.

    Positions with two columns and a box with two lengths make a flat,
    two-dimensional system: ``dimension`` is 2, velocities have two columns and
    ``box.volume`` is the area. ``boundary="walls"`` puts a wall on every face
    of the box instead of periodic images; one word per axis, such as
    ``("walls", "periodic")``, puts walls on the two faces across that axis
    only. The walls are named in ``walls``; every atom must start between them.
    """

    def __init__(self, positions, box, *, velocities=None, masses=None, molecule=None, bond=None, model=None, species=None, boundary="periodic", units="reduced",
                 unwrapped=None, time=0., step=0, tilt=None):
        x=np.array(positions,dtype=float,order='C',copy=True)
        if x.ndim!=2 or x.shape[1] not in (2,3) or not np.isfinite(x).all():
            raise ValueError('positions must be a finite (N,3) array, or (N,2) for a flat system.')
        D=x.shape[1]
        self._dim=D
        self._walled=_walled(boundary,D)
        walled=np.array(self._walled)
        self._boundary = boundary if isinstance(boundary,str) else "walls" if walled.all() else "periodic" if not walled.any() else tuple(boundary)
        self._units = _units(units).name
        model, bond = resolve_model(model, molecule, bond)
        if boundary == "open" and (model != "atomic" or D != 3):
            raise ValueError("Open boundaries currently support unbonded particles in three dimensions only.")
        if (D==2 or walled.any()) and model != "atomic":
            raise ValueError("A flat system and walls currently take unbonded atoms only.")
        self._bond = bond
        try:
            lengths=np.broadcast_to(np.asarray(box,dtype=float),(D,)).copy()
        except ValueError:
            raise ValueError(f'box must be one length or {D} lengths, to match positions with {D} columns.') from None
        if not np.isfinite(lengths).all() or np.any(lengths<=0):
            raise ValueError('box lengths must be finite and positive.')
        try:
            tilts=np.zeros(3) if tilt is None else np.broadcast_to(np.asarray(tilt,dtype=float),(3,)).copy()
        except ValueError:
            raise ValueError('tilt must be three numbers: xy, xz and yz.') from None
        if not np.isfinite(tilts).all():
            raise ValueError('tilt must be finite.')
        if tilts.any():
            if D!=3 or boundary!='periodic' or resolve_model(model, molecule, bond)[0]!='atomic':
                raise ValueError('A tilted box is periodic in three dimensions and holds unbonded atoms.')
            _check_tilt(lengths, tilts)
        # A state taken from a run may hold an atom in the middle of a bounce, a
        # little beyond a wall; further out the wall would throw it back hard.
        if np.any(x[:,walled]<-_BEYOND) or np.any(x[:,walled]>lengths[walled]+_BEYOND):
            raise ValueError('Every atom must start between the walls: coordinates from 0 to the box length. '
                             'To squeeze the atoms between walls, move a wall with a Piston.')
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
        if any(len(str(name)) > 32 for name in labels):
            raise ValueError('Species names may be at most 32 characters long.')
        # Room for any name up to 32 characters, so a label can be changed in place.
        self._species = labels.astype('<U32', copy=True)
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
        if tilts.any():
            x=wrap(x,lengths,tilts)
        elif boundary != "open":
            x[:,~walled]%=lengths[~walled]
        # The solver works in three dimensions: a flat system has z = 0 throughout
        # and a box of thickness 0. The public arrays have D columns.
        self._box=np.zeros(3)
        self._box[:D]=lengths
        # Tilt factors xy, xz, yz as in LAMMPS: the box edges are a = (lx, 0, 0),
        # b = (xy, ly, 0) and c = (xz, yz, lz); all zero for a rectangular box.
        self._tilt=tilts
        self._x=self._columns(x)
        self._mass=np.ones(len(x)) if masses is None else np.array(masses,dtype=float,order='C',copy=True)
        if self._mass.shape!=(len(x),) or not np.isfinite(self._mass).all() or np.any(self._mass<=0):
            raise ValueError('masses must be a positive finite (N,) array.')

        self._v=None if velocities is None else np.array(velocities,dtype=float,order='C',copy=True)
        if self._v is not None and (self._v.shape!=x.shape or not np.isfinite(self._v).all()):
            raise ValueError(f'velocities must be a finite (N,{D}) array.')
        self._u=None if unwrapped is None else np.array(unwrapped,dtype=float,order='C',copy=True)
        if self._u is not None and (self._u.shape!=x.shape or not np.isfinite(self._u).all()):
            raise ValueError(f'unwrapped positions must be a finite (N,{D}) array.')
        self._v=self._columns(self._v)
        self._u=self._columns(self._u)
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
        # For each wall, the impulse it has received and the work it has done on
        # the atoms by moving, along the runs leading here (None: all zero).
        self._wall_sums=None
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

    def _columns(self, array):
        """The solver's three-column array for a public one with D columns."""
        if array is None or self._dim==3:
            return None if array is None else np.ascontiguousarray(array)
        return np.ascontiguousarray(np.column_stack((array,np.zeros(len(array)))))

    @property
    def boundary(self):
        """"periodic", "walls" (a wall on every face), one of the two per axis, or "open".

        Open systems use box only as a viewing frame, without walls or images.
        It can be changed, ``state.boundary = "walls"``: every atom must then be
        inside the box, between the new walls."""
        return self._boundary

    @boundary.setter
    def boundary(self, boundary):
        self._rebuild(boundary=boundary)
        self._wall_sums = None  # new walls have received no impulse yet

    def _rebuild(self, **changes):
        """Check a change by building the state anew, then take it over in place."""
        values = dict(positions=self._x[:, :self._dim], box=self._box[:self._dim], velocities=self.velocities,
                      masses=self._mass, model=self.model, bond=self.bond, species=self._species,
                      boundary=self._boundary, units=self._units, unwrapped=None if self._u is None else self._u[:, :self._dim],
                      time=self._time, step=self._step, tilt=self._tilt)
        values.update(changes)
        other = System(values.pop('positions'), values.pop('box'), **values)
        other._rng, other._rate, other._heat, other._pppm, other._wall_sums = (
            self._rng, self._rate, self._heat, self._pppm, self._wall_sums)
        self.__dict__.update(other.__dict__)

    @property
    def dimension(self):
        """3, or 2 for a flat system: the number of columns of positions and velocities."""
        return self._dim

    @property
    def walls(self):
        """Names of the walls: "left" and "right" (at x = 0 and x = Lx), "bottom" and
        "top" (y), "back" and "front" (z). Empty for a periodic box."""
        return tuple(WALLS[2*axis+side] for axis,walled in enumerate(self._walled) if walled for side in (0,1))

    @property
    def time(self):
        """Simulation time of this state, in the system's time unit; it can be set."""
        return self._time

    @time.setter
    def time(self, value):
        value = float(value)
        if not np.isfinite(value):
            raise ValueError('time must be finite.')
        self._time = value

    @property
    def step(self):
        """Number of MD steps taken to reach this state; it can be set."""
        return self._step

    @step.setter
    def step(self, value):
        self._step = _integer(value, 'step', 0)

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
        """Number density: particles per volume, or per area for a flat system."""
        return self.N/self.box.volume if self.boundary != "open" else float("nan")

    @property
    def box(self):
        """Box lengths, one per axis (two for a flat system), which can be changed.

        ``state.box[0] *= 2`` doubles the x length; ``state.box[:] = 20, 20, 20``
        sets all three. Atoms keep their positions. The next run starts from the
        new box: across a periodic axis atoms outside it are wrapped back in, so a
        smaller box can bring atoms close together; between walls, the far wall
        moves and every atom must stay inside.
        """
        return self._box[:self._dim].view(Box)

    @property
    def tilt(self):
        """Tilt factors (xy, xz, yz) of a tilted (triclinic) box, which can be changed.

        As in LAMMPS, the box edges are a = (lx, 0, 0), b = (xy, ly, 0) and
        c = (xz, yz, lz): xy leans the y edge along x, and so on. Each tilt may be
        at most half the length it leans along. All zero for a rectangular box.
        ``state.tilt[0] = 1.5`` keeps the atoms where they are, as editing the box
        does; ``state.strain(xy=0.01)`` shears the atoms with the box.
        """
        return self._tilt

    def strain(self, *, xx=0., yy=0., zz=0., yz=0., xz=0., xy=0.):
        """Stretch or shear the box and the atoms with it, in place.

        ``xx``, ``yy``, ``zz`` change a length by that fraction; ``yz``, ``xz``, ``xy``
        are shear strains: ``xy=0.01`` leans the y edge by 0.01 times its length along
        x. Every position moves with the box (as LAMMPS ``change_box ... remap``), so a
        crystal is stretched or sheared rather than cut. Velocities are kept. For an
        elastic constant, strain by a small amount and compare the stress.
        """
        F=np.array([[1+xx,xy,xz],[0,1+yy,yz],[0,0,1+zz]],dtype=float)
        if not np.isfinite(F).all() or min(1+xx,1+yy,1+zz)<=0:
            raise ValueError('Strains must be finite, and each length must stay positive.')
        if self._dim!=3 or self.boundary!='periodic' or self.model!='atomic':
            raise ValueError('strain works on a periodic three-dimensional box of unbonded atoms.')
        lx,ly,lz=self._box
        xyt,xzt,yzt=self._tilt
        H=F@np.array([[lx,xyt,xzt],[0,ly,yzt],[0,0,lz]])
        lengths,tilts=np.diag(H).copy(),np.array([H[0,1],H[0,2],H[1,2]])
        _check_tilt(lengths,tilts)
        self._x=np.ascontiguousarray(wrap(self._x@F.T,lengths,tilts))
        if self._u is not None:
            self._u=np.ascontiguousarray(self._u@F.T)
        self._box[:]=lengths
        self._tilt[:]=tilts

    def _public(self, array):
        return None if array is None else _readonly(array[:,:self._dim])

    @property
    def positions(self):
        """Positions, one row per atom: (x, y, z), or (x, y) for a flat system.

        They can be changed: ``state.positions[0, 0] += 0.5`` moves atom 0 along x,
        and ``state.positions = new`` replaces them all. The next run checks them:
        atoms outside a periodic box are wrapped back in, and between walls every
        atom must be inside.
        """
        return self._x[:, :self._dim]

    @positions.setter
    def positions(self, value):
        value = np.asarray(value, dtype=float)
        if value.shape != (len(self._x), self._dim):
            raise ValueError(f'positions must have one row per atom: shape ({len(self._x)}, {self._dim}).')
        self._x[:, :self._dim] = value

    @property
    def velocities(self):
        """Velocities, one row per atom, or None: then a run starts with random ones at its temperature.

        They can be changed: ``state.velocities *= 1.1`` heats the system up, and
        ``state.velocities[0] += (1, 0, 0)`` kicks atom 0. A run along a periodic
        direction needs zero total momentum: after kicking one atom, subtract the
        mean, ``state.velocities -= np.average(state.velocities, axis=0, weights=state.masses)``.
        """
        return None if self._v is None else self._v[:, :self._dim]

    @velocities.setter
    def velocities(self, value):
        if value is None:
            self._v = None
            return
        value = np.asarray(value, dtype=float)
        if value.shape != (len(self._x), self._dim):
            raise ValueError(f'velocities must have one row per atom: shape ({len(self._x)}, {self._dim}).')
        if self._v is None:
            self._v = np.zeros_like(self._x)
        self._v[:, :self._dim] = value

    @property
    def masses(self):
        """Mass of each atom; they can be changed, ``state.masses[:10] = 4``."""
        return self._mass

    @masses.setter
    def masses(self, value):
        value = np.asarray(value, dtype=float)
        if value.shape != self._mass.shape:
            raise ValueError(f'masses must have one value per atom: shape ({len(self._mass)},).')
        self._mass[:] = value

    @property
    def atoms(self):
        """Positions, velocities, masses, molecule numbers and species of the atoms.

        positions, velocities, masses and species are the state's own arrays and can
        be changed in place (``state.atoms.species[:50] = "B"``); molecule numbers
        follow from the model.
        """
        ids=np.arange(len(self._x))//_SIZE[self.model]
        return Atoms(self.positions,self.velocities,self._mass,_readonly(ids),self._species)

    @property
    def _conserved(self):
        """For each of the solver's three axes: is the total momentum along it conserved?"""
        return np.array([not walled for walled in self._walled]+[False]*(3-self._dim))

    @property
    def dof(self):
        """Degrees of freedom in the kinetic temperature: one per atom and direction,
        less one for each direction that conserves the total momentum (a direction
        with walls does not), less the rigid constraints."""
        return self._dim*len(self._x)-int(self._conserved.sum())-{'diatomic-rigid': self.N, 'water-rigid': 3*self.N}.get(self.model, 0)

    def copy(self):
        """Independent copy of the complete state."""
        return self._copy(self._public(self._v))

    def _copy(self, velocities):
        other=System(self._public(self._x),self._box[:self._dim],velocities=velocities,masses=self._mass,model=self.model,bond=self.bond,species=self._species,boundary=self.boundary,units=self.units,
                     unwrapped=self._public(self._u),time=self._time,step=self._step,tilt=self._tilt)
        other._rng,other._rate,other._heat=self._rng,self._rate,self._heat
        other._pppm=None if self._pppm is None else (self._pppm[0].copy(),)+tuple(self._pppm[1:])
        other._wall_sums=None if self._wall_sums is None else self._wall_sums.copy()
        return other

    def with_velocities(self, velocities):
        """Copy of this state with other velocities, one row per atom.

        ``state.with_velocities(-state.atoms.velocities)`` runs the motion
        backwards; changing one atom's row perturbs it. Positions, box, masses
        and the clock are kept. The total momentum must stay zero for a run, and
        velocities along a rigid bond are rejected.
        """
        if velocities is None:
            raise ValueError(f'velocities must be a finite (N,{self._dim}) array, one row per atom.')
        return self._copy(velocities)

    def with_box(self, box):
        """Copy of this state in another box, with nothing stretched.

        Along a periodic direction the box can only grow, and the contents are
        centred in it. Molecules keep their positions relative to each other
        and stay whole. Enlarging one side turns a bulk liquid into a slab
        with two free surfaces and vacuum on both sides, for evaporation:
        ``slab = state.with_box([Lx, Ly, 4*Lz])``.

        Along a direction with walls every atom keeps its position: the wall at
        0 stays and the other wall is put at the new length, which must leave
        every atom inside. ``gas.with_box((200, 100))`` gives the gas of a
        100 x 100 box twice the room, as when a partition is removed.

        Velocities, masses and the clock are kept; the barostat velocity and
        PPPM mesh choice are reset.
        """
        if self.boundary=='open':
            raise ValueError('with_box needs a box: periodic or with walls.')
        if self._tilt.any():
            raise ValueError('with_box works on a rectangular box; this one is tilted.')
        D=self._dim
        old=self._box[:D]
        walled=np.array(self._walled)
        try:
            new=np.broadcast_to(np.asarray(box,dtype=float),(D,)).copy()
        except ValueError:
            raise ValueError(f'box must be one length or {D} lengths.') from None
        if not np.isfinite(new).all() or np.any(new<=0) or np.any(new[~walled]<old[~walled]*(1-1e-12)):
            raise ValueError(f'Along a periodic direction each new box length must be at least the current one, {old.tolist()}.')
        if np.any(self._x[:,:D][:,walled]>new[walled]+_BEYOND):
            raise ValueError('The new walls must leave every atom inside the box.')
        size=_SIZE[self.model]
        x=self._x[:,:D].reshape(-1,size,D)
        m=self._mass.reshape(-1,size)
        # Whole molecules: offsets from the first atom by minimum image, centre of
        # mass wrapped into the old box, then the old box centred in the new one.
        free=~walled
        offset=x-x[:,:1]
        offset[...,free]-=old[free]*np.rint(offset[...,free]/old[free])
        centre=x[:,0]+(m[:,:,None]*offset).sum(axis=1)/m.sum(axis=1)[:,None]
        centre[:,free]-=old[free]*np.floor(centre[:,free]/old[free])
        local=offset-(m[:,:,None]*offset).sum(axis=1)[:,None,:]/m.sum(axis=1)[:,None,None]
        positions=(centre[:,None,:]+local+np.where(walled,0.,.5*(new-old))).reshape(-1,D)
        unwrapped=None if self._u is None else self._u[:,:D]+(positions-self._x[:,:D])
        other=System(positions,new,velocities=self._public(self._v),masses=self._mass,model=self.model,bond=self.bond,species=self._species,
                     boundary=self.boundary,units=self.units,unwrapped=unwrapped,time=self._time,step=self._step)
        other._rng,other._heat=self._rng,self._heat
        other._wall_sums=None if self._wall_sums is None else self._wall_sums.copy()
        return other

    def draw(self, ax=None, *, view=None, projection='perspective', radius=None, arrows=True):
        """Draw the atoms as balls in a matplotlib figure and return its Axes.

        ``view`` is where you look from, with z up: "angled" (the default),
        "front" (the xz plane), "side" (the yz plane), "top" (the xy plane), or
        (elevation, azimuth) in degrees as in matplotlib's 3D plots. Without a
        view, ``projection="orthographic"`` gives "front", and a box much longer in
        z (a slab) lies across the picture. Small x, y, z arrows show the directions
        (``arrows=False`` leaves them out). A flat system is drawn as it is, with
        axes. ``radius`` is a display radius (a number, or one per species); by
        default 0.5 (sigma = 1), and the player's sizes for water. Give ``ax`` to
        draw into a figure of your own; ``ax.figure.savefig(...)`` saves it.
        """
        from .drawing import draw_atoms, display_radii
        if radius is None and self.boundary == 'open':
            radius = .035
        return draw_atoms(self._x[:, :self._dim], self._box[:self._dim], self._species,
                          display_radii(self._species, self.model, radius, self.units), bonds=self.bonds,
                          walled=tuple(self._walled), framed=self.boundary != 'open', ax=ax, view=view,
                          projection=projection, arrows=arrows, tilt=self._tilt)

    @classmethod
    def from_arrays(cls, positions, velocities, box, *, masses=None, molecule=None, bond=None, model=None, species=None, boundary="periodic", units="reduced"):
        return cls(positions,box,velocities=velocities,masses=masses,model=model,molecule=molecule,bond=bond,species=species,boundary=boundary,units=units)

    def __repr__(self):
        shape = '' if self._dim==3 and self._boundary=='periodic' else f'dimension={self._dim}, boundary={self._boundary!r}, '
        if self._tilt.any():
            shape += f'tilt={self._tilt.tolist()}, '
        return f'System({shape}model={self.model!r}, N={self.N}, rho={self.rho:g}, velocities={self._v is not None}, time={self._time:g}, step={self._step})'
