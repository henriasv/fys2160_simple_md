"""Integration, experiment controls and recording for an owned physical system."""
from __future__ import annotations
from datetime import datetime, timezone
from pathlib import Path
import csv
from functools import wraps
import json
import shutil
import time
import numpy as np
from . import _core
from .system import System, MODELS, _positive, _integer, _readonly
from .results import Run, write_json
from .species import parameter, per_atom, charges as _charges
from .bonds import bond_config, read_bond
from .units import units as _units, AVOGADRO_PER_CM3_A3

def _storage_name(value):
    if not isinstance(value, str) or not value or value != value.strip() or len(value) > 100:
        raise ValueError('storage_name must be a nonempty folder name (at most 100 characters).')
    if value.startswith('.') or any(c in value for c in '/\\\0<>:"|?*'):
        raise ValueError('storage_name must be a plain folder name, such as atomic-gas.')
    return value


def _replace_run_directory(path):
    """Overwrite only this package's own run directories, never arbitrary folders."""
    if path.is_symlink():
        raise ValueError('The storage folder must not be a symbolic link.')
    if path.exists():
        try:
            metadata = json.loads((path/'metadata.json').read_text())
        except (OSError, ValueError) as exc:
            raise ValueError(f'{path} already exists and is not a saved MD run; choose another storage_name.') from exc
        if not isinstance(metadata, dict) or metadata.get('schema_version') != 1 or 'package_version' not in metadata or 'configuration' not in metadata:
            raise ValueError(f'{path} is not a saved MD run; choose another storage_name.')
        shutil.rmtree(path)
    (path/'frames').mkdir(parents=True)


def _legacy_output_names(method):
    @wraps(method)
    def wrapped(self, *args, **kwargs):
        for old, new in [('sample_every', 'thermo_every'), ('save_every', 'trajectory_every')]:
            if old in kwargs:
                if new in kwargs:
                    raise TypeError(f'Use {new}, not both {old} and {new}.')
                kwargs[new] = kwargs.pop(old)
        return method(self, *args, **kwargs)
    return wrapped


class MDSimulation:
    """An experiment evolving its own copy of a physical System.

    Missing velocities are initialized at temperature. Explicit velocities are
    preserved; temperature then specifies only the thermostat target. Repeated
    run calls continue this simulation's clock, state and physical settings.
    """

    def __init__(self, system, *, pair_potential='lj', epsilon=1., sigma=1., mixing_rule='lorentz-berthelot', exclude_bonded_pairs=True, temperature=None, ensemble='nvt', timestep=None,
                 cutoff=None, G=1., softening=.05, skin=None, friction=None, pressure=None, pressure_time=None,
                 heat_rate=0., seed=87287, output_dir='runs', storage_name='simulation', units=None,
                 charge=None, coulomb=None, coulomb_cutoff=None, kspace_accuracy=1e-5, lj_shift=True):
        if not isinstance(system,System):
            raise TypeError('Pass a System, such as FCC(N=500, rho=0.1).')
        if len(system.atoms)>20000:
            raise ValueError('This course solver supports at most 20000 atoms.')
        if pair_potential not in ('lj', 'lj96', 'gravity'):
            raise ValueError('pair_potential must be lj, lj96 or gravity.')
        self._units = _units(system.units if units is None else units)
        if self._units.name != system.units:
            raise ValueError(f'units="{self._units.name}" does not match the system units ("{system.units}").')
        if (charge is None) != (coulomb is None):
            raise ValueError('Charges need coulomb="pppm" (long-range electrostatics), and coulomb needs charge.')
        if coulomb not in (None, 'pppm'):
            raise ValueError('coulomb must be "pppm".')
        real = self._units.name == 'real'
        self._engine = 'molecular' if (system.model == 'water-rigid' or real or charge is not None) else 'pair'
        if self._engine == 'molecular':
            if pair_potential != 'lj':
                raise ValueError('Water, charges and real units use the 12-6 Lennard-Jones potential (pair_potential="lj").')
            if system.model.startswith('diatomic'):
                raise ValueError('Diatomic molecules use reduced units without charges.')
        if not isinstance(lj_shift, bool):
            raise TypeError('lj_shift must be True or False.')
        if self._engine == 'pair' and not lj_shift:
            raise ValueError('lj_shift=False is available for water, charges and real units.')
        if self._engine == 'pair' and coulomb_cutoff is not None:
            raise ValueError('coulomb_cutoff applies only to charged systems.')
        if real:
            missing = [name for name, value in (('temperature', temperature), ('timestep', timestep), ('cutoff', cutoff)) if value is None]
            if missing:
                raise ValueError(f'Real units need explicit {", ".join(missing)} (K, fs, Å).')
        defaults = dict(skin=2., friction=.005, pressure=None, pressure_time=1000.) if real else dict(temperature=2., skin=.4, friction=1., pressure=.01, pressure_time=5.)
        temperature = defaults['temperature'] if temperature is None else temperature
        skin = defaults['skin'] if skin is None else skin
        friction = defaults['friction'] if friction is None else friction
        pressure = defaults['pressure'] if pressure is None else pressure
        pressure_time = defaults['pressure_time'] if pressure_time is None else pressure_time
        if pair_potential == 'gravity':
            if system.boundary != 'open' or system.model != 'atomic':
                raise ValueError('Gravity requires an unbonded System with boundary="open".')
            if cutoff is not None:
                raise ValueError('Isolated gravity has no cutoff; use cutoff=None.')
            if len(system.atoms) > 4096:
                raise ValueError('Direct all-pairs gravity supports at most 4096 particles.')
            if epsilon != 1. or sigma != 1.:
                raise ValueError('Gravity uses G and softening, not LJ epsilon/sigma.')
        elif system.boundary != 'periodic':
            raise ValueError('Lennard–Jones currently requires periodic boundaries.')
        self._G = _positive(G, 'G')
        self._softening = _positive(softening, 'softening', zero=True)
        if not isinstance(exclude_bonded_pairs, bool):
            raise TypeError('exclude_bonded_pairs must be True or False.')
        if mixing_rule != 'lorentz-berthelot':
            raise ValueError('mixing_rule must be "lorentz-berthelot".')
        zero = self._engine == 'molecular'
        self._epsilon = parameter(epsilon, system.species, 'epsilon', zero=zero)
        self._sigma = parameter(sigma, system.species, 'sigma', zero=zero)
        self._charge = None if charge is None else _charges(charge, system.species)
        self._coulomb = coulomb
        self._lj_shift = lj_shift
        self._pair_parameters = np.ascontiguousarray(np.column_stack((per_atom(self._epsilon, system._species), per_atom(self._sigma, system._species))))
        self._pair_potential = pair_potential
        self._exclude_bonded_pairs = exclude_bonded_pairs
        self._system=system.copy()

        self._x=self._system._x
        self._box=self._system._box
        self._mass=self._system._mass
        self.cutoff=None if pair_potential == 'gravity' else _positive(2.5 if cutoff is None else cutoff,'cutoff')
        self.skin=_positive(skin,'skin')
        self.friction=_positive(friction,'friction')
        self.pressure_time=_positive(pressure_time,'pressure_time')
        self.seed=_integer(seed,'seed')
        if self.seed>=2**64:
            raise ValueError('seed must be below 2**64.')
        self.output_dir=Path(output_dir).expanduser()
        self.storage_name=_storage_name(storage_name)
        self._settings=dict(ensemble=ensemble,
                            timestep=(.005 if self.model=='atomic' else .002) if timestep is None else timestep,
                            temperature=temperature,pressure=pressure,heat_rate=heat_rate)
        if self._engine == 'molecular':
            self._settings.update(coulomb_cutoff=self.cutoff if coulomb_cutoff is None else coulomb_cutoff,
                                  kspace_accuracy=kspace_accuracy)
            self._prepare_molecular()
        self._deform = None
        self._validate_settings(self._settings)
        if self._system._v is None:
            rng=np.random.default_rng(self.seed)
            self._v=np.ascontiguousarray(rng.normal(size=self._x.shape)/np.sqrt(self._mass[:,None]))
            self._v-=(self._v*self._mass[:,None]).sum(axis=0)/self._mass.sum()
            if self.model=='diatomic-rigid':self._project_velocities()
            if self.model=='water-rigid':_core.project_molecular_velocities(self._x,self._v,self._box,self._mass,self._constraints)
            self._v*=np.sqrt(self.dof*self._units.boltzmann*self._settings['temperature']/(self._units.mvv2e*np.sum(self._mass[:,None]*self._v**2)))
            self._system._v=self._v
        else:
            self._v=self._system._v
            momentum=self._v*self._mass[:,None]
            if np.linalg.norm(momentum.sum(axis=0))>1e-10*max(1.,np.linalg.norm(momentum)):
                raise ValueError('Remove centre-of-mass velocity before starting this solver.')
        self._u=self._x.copy()
        self._f=np.zeros_like(self._x)
        self.step=0
        self.time=0.
        self.heat_added=0.
        self._rng=self.seed
        self._rate=0.
        self._baromass=(self.dof+3)*max(self._settings['temperature'],1e-6)*self.pressure_time**2
        if self._engine == 'molecular':
            # Barostat inertia from centre-of-mass degrees of freedom (3N-3)+3.
            self._baromass=3*self.particles*self._units.boltzmann*max(self._settings['temperature'],1e-6)*self.pressure_time**2
        self._failed=False
        self._advance(0)

    def _prepare_molecular(self):
        """Type tables, per-atom charges and rigid-water geometry for the native solver."""
        system = self._system
        names = system.species
        index = {name: k for k, name in enumerate(names)}
        self._types = np.ascontiguousarray([index[s] for s in system._species], dtype=np.intc)
        value = lambda p, name: p[name] if isinstance(p, dict) else p
        eps = np.array([value(self._epsilon, n) for n in names])
        sig = np.array([value(self._sigma, n) for n in names])
        if np.any((eps > 0) & (sig <= 0)):
            raise ValueError('sigma must be positive for every species with epsilon > 0.')
        # Lorentz-Berthelot mixing; species with epsilon=0 have no LJ interactions.
        self._epsilon_table = np.ascontiguousarray(np.sqrt(np.outer(eps, eps)))
        self._sigma_table = np.ascontiguousarray(.5*(sig[:, None]+sig[None, :]))
        q = per_atom(self._charge, system._species) if self._charge is not None else np.zeros(len(system._species))
        if self._charge is not None and abs(q.sum()) > 1e-8*max(1., np.abs(q).sum()):
            raise ValueError(f'The total charge must be zero; these charges sum to {q.sum():.6g} e.')
        self._charges = np.ascontiguousarray(q)
        self._pppm_state = np.zeros(8)
        water = system.bond if system.model == 'water-rigid' else None
        self._constraints = water.constraints if water else (1., 1., 1.)
        self._alpha = water.alpha if water else 0.

    def _validate_settings(self, config):
        if self.pair_potential == 'gravity' and config['ensemble'] != 'nve':
            raise ValueError('Isolated gravity uses ensemble="nve" (no thermostat, walls or barostat).')
        if config['ensemble'] not in ('nve', 'nvt', 'nph', 'npt'):
            raise ValueError('ensemble must be nve, nvt, nph or npt.')
        if config['ensemble'] in ('nph', 'npt') and config['pressure'] is None:
            raise ValueError(f'Set pressure ({self._units.pressure_unit}) for {config["ensemble"]}.')
        for key in ('timestep', 'temperature', 'pressure'):
            if config[key] is not None:
                config[key] = _positive(config[key], key, zero=(key != 'timestep'))
        if self._engine == 'molecular':
            config['coulomb_cutoff'] = _positive(config['coulomb_cutoff'], 'coulomb_cutoff')
            config['kspace_accuracy'] = _positive(config['kspace_accuracy'], 'kspace_accuracy')
            if config['kspace_accuracy'] >= 1:
                raise ValueError('kspace_accuracy is a relative force accuracy, such as 1e-5.')
        config['heat_rate'] = float(config['heat_rate'])
        if not np.isfinite(config['heat_rate']):
            raise ValueError('heat_rate must be finite (total reduced energy per time).')
        if config['heat_rate'] and config['ensemble'] in ('nvt', 'npt'):
            raise ValueError('Turn off the thermostat (nve or nph) before adding/removing heat.')
        if self._engine == 'pair' and self.model != 'atomic' and config['ensemble'] in ('nph', 'npt'):
            raise ValueError('Pressure control is supported for atoms and rigid water only.')

    @property
    def G(self):
        return self._G

    @property
    def softening(self):
        return self._softening

    @property
    def units(self):
        return self._units.name

    @property
    def charge(self):
        return self._charge.copy() if isinstance(self._charge, dict) else self._charge

    @property
    def coulomb(self):
        return self._coulomb

    @property
    def lj_shift(self):
        return self._lj_shift

    @property
    def epsilon(self):
        return self._epsilon.copy() if isinstance(self._epsilon, dict) else self._epsilon

    @property
    def sigma(self):
        return self._sigma.copy() if isinstance(self._sigma, dict) else self._sigma

    @property
    def mixing_rule(self):
        return 'lorentz-berthelot'

    @property
    def pair_potential(self):
        return self._pair_potential

    @property
    def exclude_bonded_pairs(self):
        return self._exclude_bonded_pairs

    @property
    def system(self):
        return self._system

    @property
    def model(self):
        return self.system.model

    @property
    def particles(self):
        return self.system.N

    @property
    def N(self):
        return self.system.N

    @property
    def rho(self):
        return self.system.rho

    def _bond_vectors(self):
        dr = self._x[::2]-self._x[1::2]
        return dr-self._box*np.rint(dr/self._box)

    def _project_velocities(self):
        r = self._bond_vectors()
        dot = (r*(self._v[::2]-self._v[1::2])).sum(axis=1)
        inv_mass = 1/self._mass[::2] + 1/self._mass[1::2]
        correction = dot[:,None]*r/((r*r).sum(axis=1)*inv_mass)[:,None]
        self._v[::2] -= correction/self._mass[::2,None]
        self._v[1::2] += correction/self._mass[1::2,None]

    @property
    def dof(self):
        return self.system.dof

    @property
    def box(self):
        return self.system.box

    @property
    def atoms(self):
        return self.system.atoms

    @property
    def forces(self):
        return _readonly(self._f)

    @property
    def settings(self):
        return self._settings.copy()

    def _advance(self, steps):
        c = self._settings
        if self.pair_potential == 'gravity':
            self._potential, self._virial, done, interrupted = _core.advance_gravity(
                self._x, self._v, self._f, self._mass, steps, c['timestep'],
                self.G, self.softening, c['heat_rate'])
            self._u[:] = self._x
            self.step += done
            self.time += done*c['timestep']
            self.heat_added += done*c['timestep']*c['heat_rate']
            self._pressure_numerator = float('nan')
            return bool(interrupted)
        if self._engine == 'molecular':
            return self._advance_molecular(steps)
        result = _core.advance(self._x, self._v, self._u, self._f, self._box,
            self._mass, steps, c['timestep'], self.cutoff, self.skin,
            c['temperature'], self.friction if c['ensemble'] in ('nvt','npt') else 0.,
            self._rng, MODELS[self.model], c['heat_rate'],
            c['pressure'] if c['ensemble'] in ('nph','npt') else -1., self._baromass, self._rate,
            int(self.pair_potential == 'lj96'), int(self.exclude_bonded_pairs), .7, 0., 0., 0., 1., 1., self._pair_parameters, self.system._bond_parameters)
        self._potential, self._virial, self._rng, done, interrupted, self._rate, self._pressure_numerator = result
        self.step += done
        self.time += done*c['timestep']
        self.heat_added += done*c['timestep']*c['heat_rate']
        return bool(interrupted)

    def _advance_molecular(self, steps):
        c = self._settings
        U = self._units
        controlled = c['ensemble'] in ('nph', 'npt')
        deform = None
        if self._deform is not None:
            lengths, scale, start, total = self._deform
            deform = np.array([*lengths, scale, self.step-start, total], dtype=float)
        result = _core.advance_molecular(x=self._x, v=self._v, u=self._u, f=self._f, box=self._box,
            mass=self._mass, types=self._types, charges=self._charges, epsilon=self._epsilon_table,
            sigma=self._sigma_table, pppm_state=self._pppm_state, deform=deform, steps=steps,
            timestep=c['timestep'], cutoff=self.cutoff, coulomb_cutoff=c['coulomb_cutoff'], skin=self.skin,
            temperature=U.boltzmann*c['temperature'], friction=self.friction if c['ensemble'] in ('nvt','npt') else 0.,
            seed=self._rng, heat_rate=c['heat_rate'], pressure=c['pressure']/U.pressure if controlled else -1.,
            barostat_mass=self._baromass, rate=self._rate, mvv2e=U.mvv2e, qqr2e=U.qqr2e,
            molecule_size=3 if self.model == 'water-rigid' else 1, constraints=self._constraints,
            alpha=self._alpha, lj_shift=self._lj_shift, coulomb=int(self._charge is not None),
            accuracy=c['kspace_accuracy'], order=5)
        self._terms = result
        self._potential = result['lj'] + result['coulomb'] + result['kspace']
        self._virial = result['virial']
        self._rng, self._rate = result['rng'], result['rate']
        self._pressure_numerator = result['pressure_numerator']
        done = result['done']
        self.step += done
        self.time += done*c['timestep']
        self.heat_added += done*c['timestep']*c['heat_rate']
        return bool(result['interrupted'])

    @property
    def observables(self):
        U = self._units
        K = float(.5*np.sum(self._mass[:,None]*self._v**2))*U.mvv2e
        V = float(np.prod(self._box))
        T = 2*K/(self.dof*U.boltzmann)
        P_internal = self._pressure_numerator/(3*V)
        P = P_internal*U.pressure
        E = K+self._potential
        if self.pair_potential == 'gravity':
            center = np.average(self._x, axis=0, weights=self._mass)
            radius = np.linalg.norm(self._x-center, axis=1)
            order = np.argsort(radius)
            half = np.searchsorted(np.cumsum(self._mass[order]), .5*self._mass.sum())
            return dict(step=self.step, time=self.time, temperature=T,
                        kinetic_energy=K, potential_energy=self._potential,
                        total_energy=E, virial=self._virial,
                        virial_ratio=2*K/(-self._virial) if self._virial else float('nan'),
                        half_mass_radius=float(radius[order[half]]), heat_added=self.heat_added,
                        pressure=float('nan'))
        controlled = self._settings['ensemble'] in ('nph','npt')
        H = E+self._settings['pressure']/U.pressure*V if controlled else E+P_internal*V
        B = .5*self._baromass*self._rate**2 if controlled else 0.
        values = dict(step=self.step, time=self.time, temperature=T, pressure=P,
                      kinetic_energy=K, potential_energy=self._potential,
                      total_energy=E, volume=V, density=self.particles/V,
                      atom_density=len(self._x)/V, compressibility_factor=P_internal*V/(self.particles*T*U.boltzmann) if T else float('nan'),
                      enthalpy=H, barostat_energy=B, extended_enthalpy=H+B,
                      heat_added=self.heat_added)
        if U.name == 'real':
            values['mass_density'] = float(self._mass.sum())/(AVOGADRO_PER_CM3_A3*V)
        if self._engine == 'molecular' and self._charge is not None:
            values['lj_energy'] = self._terms['lj']
            values['coulomb_energy'] = self._terms['coulomb']+self._terms['kspace']
        return values

    def _config(self):
        return dict(model=self.model, bond=bond_config(self.system.bond), units=self.units,
                    charge=self.charge, coulomb=self.coulomb, lj_shift=self.lj_shift,
                    boundary=self.system.boundary, G=self.G, softening=self.softening,
                    pair_potential=self.pair_potential, epsilon=self.epsilon, sigma=self.sigma, mixing_rule=self.mixing_rule, exclude_bonded_pairs=self.exclude_bonded_pairs, particles=self.particles,
                    density=self.rho if self.system.boundary == "periodic" else None, seed=self.seed,
                    cutoff=self.cutoff, skin=self.skin, friction=self.friction,
                    pressure_time=self.pressure_time, storage_name=self.storage_name, **self.settings)

    def _save_checkpoint(self, path):
        with (path/'checkpoint.tmp').open('wb') as stream:
            np.savez_compressed(stream, positions=self._x, velocities=self._v,
                unwrapped=self._u, masses=self._mass, box=self._box, species=self.system._species,
                rng=np.uint64(self._rng), step=self.step, time=self.time,
                heat_added=self.heat_added, barostat_rate=self._rate, barostat_mass=self._baromass,
                **({'pppm_state': self._pppm_state} if self._engine == 'molecular' else {}))
        (path/'checkpoint.tmp').replace(path/'checkpoint.npz')

    @classmethod
    def from_run(cls, run, *, output_dir=None):
        """Explicitly resume the final checkpoint of a completed/interrupted run."""
        run = run if isinstance(run, Run) else Run.load(run)
        if run.metadata['status'] not in ('completed', 'interrupted'):
            raise ValueError('Only completed or interrupted runs can be resumed.')
        c = run.metadata['configuration'].copy()
        c['output_dir'] = output_dir or run.path.parent
        c.pop('particles',None)
        c.pop('density',None)
        c.setdefault('pair_potential', 'lj' if c['model'] == 'atomic' else 'lj96')
        bond_args = {'bond': read_bond(c.pop('bond'))} if 'bond' in c else {}
        with np.load(run.path/'checkpoint.npz', allow_pickle=False) as data:
            system=System(data['positions'],data['box'],velocities=data['velocities'],masses=data['masses'],model=c.pop('model'),boundary=c.pop('boundary', 'periodic'),species=data['species'] if 'species' in data else None,units=c.get('units','reduced'),**bond_args)
            obj=cls(system,**c)
            for attr, key in (('_x','positions'),('_v','velocities'),('_u','unwrapped'),('_mass','masses'),('_box','box')):
                getattr(obj,attr)[:] = data[key]
            obj._f = np.zeros_like(obj._x)
            obj._rng = int(data['rng'])
            obj.step = int(data['step'])
            obj.time = float(data['time'])
            obj.heat_added = float(data['heat_added'])
            obj._rate = float(data['barostat_rate'])
            obj._baromass = float(data['barostat_mass'])
            if 'pppm_state' in data:
                obj._pppm_state[:] = data['pppm_state']
        obj._advance(0)
        return obj

    @_legacy_output_names
    def run(self, timestep=None, steps=1000, *, ensemble=None, temperature=None,
            heat_rate=None, pressure=None, thermo_every=100, trajectory_every=500,
            show=False, max_fps=None, frame_every=100, storage_name=None,
            coulomb_cutoff=None, kspace_accuracy=None, deform=None):
        """Advance and save a run; omitted physical settings retain their values.

        Results go to output_dir/storage_name. Reusing a name replaces the old
        result, including its trajectory. The constructor's name is the default.
        Sampling and trajectory intervals are in steps. ``trajectory_every=None``
        disables trajectory frames, but preserves thermo and a final checkpoint.
        ``show=True`` gives an interactive live 3D view in a Jupyter notebook.
        ``max_fps`` optionally paces execution at at most that many displayed
        frames/second, with ``frame_every`` MD steps/frame. None runs at full
        speed with display-only throttling. Pacing does not change physical dt.

        For charged systems, coulomb_cutoff and kspace_accuracy may change per run;
        the PPPM mesh is chosen again at the start of every run, as in LAMMPS.
        ``deform=s`` (nve/nvt only) scales every box length linearly in time to
        s times its starting value, moving whole molecules with their centres of
        mass. An interrupted ramp stops where it is; a new run starts a new ramp.
        """
        if self._failed:
            raise RuntimeError('This simulation failed; create a new one or resume a saved checkpoint.')
        storage_name = _storage_name(self.storage_name if storage_name is None else storage_name)
        steps = _integer(steps, 'steps', 0)
        thermo_every = _integer(thermo_every, 'thermo_every')
        if trajectory_every is not None:
            trajectory_every = _integer(trajectory_every, 'trajectory_every')
        frame_every = _integer(frame_every, 'frame_every')
        if max_fps is not None:
            max_fps = _positive(max_fps, 'max_fps')
            if not show:
                raise ValueError('max_fps is a live-view option; set show=True.')
        config = self.settings
        overrides = dict(timestep=timestep, ensemble=ensemble, temperature=temperature,
                         heat_rate=heat_rate, pressure=pressure)
        if self._engine == 'molecular':
            overrides.update(coulomb_cutoff=coulomb_cutoff, kspace_accuracy=kspace_accuracy)
            if self._charge is None and (coulomb_cutoff is not None or kspace_accuracy is not None):
                raise ValueError('coulomb_cutoff and kspace_accuracy apply only to charged systems.')
        elif coulomb_cutoff is not None or kspace_accuracy is not None:
            raise ValueError('coulomb_cutoff and kspace_accuracy apply only to charged systems.')
        for key, value in overrides.items():
            if value is not None:
                config[key] = value
        self._validate_settings(config)
        if deform is not None:
            deform = _positive(deform, 'deform')
            if self._engine != 'molecular':
                raise ValueError('deform is available for water, charges and real units.')
            if config['ensemble'] not in ('nve', 'nvt'):
                raise ValueError('deform needs ensemble="nve" or "nvt" (no barostat).')
        live = None
        if show:
            from .visualization import LiveView
            live = LiveView(max_fps=max_fps, frame_every=frame_every)
        was_controlled = self._settings['ensemble'] in ('nph','npt')
        self._settings = config
        if not was_controlled and config['ensemble'] in ('nph','npt'):
            self._rate = 0.
        if self._engine == 'molecular':
            self._pppm_state[0] = 0.  # choose the PPPM mesh again, as LAMMPS does at each run
        self._advance(0)
        start = self.step
        if deform is not None and steps:
            self._deform = (self._box.copy(), deform, start, steps)
        path = self.output_dir.resolve()/storage_name
        _replace_run_directory(path)
        meta = dict(schema_version=1, package_version='0.1.0', label=storage_name,
                    storage_name=storage_name,
                    created_at=datetime.now(timezone.utc).isoformat(), status='running',
                    configuration={**self._config(), 'storage_name': storage_name}, atoms=len(self._x),
                    molecules=self.particles if self.model != 'atomic' else 0,
                    atom_species=self.system._species.tolist(),
                    degrees_of_freedom=self.dof, units=self._units.description, deform=deform,
                    potential_shifted=self.pair_potential != 'gravity' and self.lj_shift, pressure_estimator=('not defined for an unconfined cluster' if self.pair_potential == 'gravity' else 'atomic virial' if self.model == 'atomic' else 'molecular COM virial'),
                    start_step=start, start_time=self.time, requested_steps=steps,
                    thermo_every=thermo_every, trajectory_every=trajectory_every,
                    live_view=bool(show), max_fps=max_fps, frame_every=frame_every if show else None)
        write_json(path/'metadata.json', meta)
        started = time.monotonic()
        frame_index = 0
        last_sample = last_frame = -1
        status = 'completed'
        def frame():
            nonlocal frame_index, last_frame
            np.savez_compressed(path/'frames'/f'{frame_index:07d}.npz',
                positions=self._x, unwrapped=self._u, velocities=self._v,
                box=self._box, step=self.step, time=self.time)
            frame_index += 1
            last_frame = self.step
        try:
            if live:
                live.start(self)
            with (path/'thermo.csv').open('w', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=list(self.observables))
                writer.writeheader()
                writer.writerow(self.observables)
                last_sample = self.step
                if trajectory_every is not None:
                    frame()
                while self.step-start < steps:
                    elapsed = self.step-start
                    chunk = min(steps-elapsed, thermo_every-elapsed%thermo_every, 256)
                    if trajectory_every is not None:
                        chunk = min(chunk, trajectory_every-elapsed%trajectory_every)
                    if live:
                        chunk = min(chunk, frame_every-elapsed%frame_every)
                    interrupted = self._advance(chunk)
                    elapsed = self.step-start
                    final = elapsed == steps or interrupted
                    if elapsed%thermo_every == 0 or final:
                        writer.writerow(self.observables)
                        stream.flush()
                        last_sample = self.step
                    if trajectory_every is not None and (elapsed%trajectory_every == 0 or final):
                        frame()
                    if live and (elapsed%frame_every == 0 or final):
                        live.update(self, force=final)
                    if interrupted:
                        status = 'interrupted'
                        break
        except KeyboardInterrupt:
            status = 'interrupted'
        except Exception as exc:
            status = 'failed'
            self._failed = True
            meta['error'] = f'{type(exc).__name__}: {exc}'
            raise
        finally:
            self._deform = None
            if live:
                live.finish(self, status)
            if status != 'failed':
                if last_sample != self.step:
                    with (path/'thermo.csv').open('a', newline='') as stream:
                        csv.DictWriter(stream, fieldnames=list(self.observables)).writerow(self.observables)
                if trajectory_every is not None and last_frame != self.step:
                    frame()
                self._save_checkpoint(path)
            meta.update(status=status, completed_steps=self.step-start, end_step=self.step,
                        end_time=self.time, wall_seconds=time.monotonic()-started)
            if self._engine == 'molecular' and self._charge is not None:
                state = self._pppm_state
                meta['pppm'] = dict(mesh=[int(n) for n in state[:3]], g_ewald=float(state[3]), order=5,
                                    estimated_relative_accuracy=float(state[7]), kspace_accuracy=self._settings['kspace_accuracy'],
                                    coulomb_cutoff=self._settings['coulomb_cutoff'])
            write_json(path/'metadata.json', meta)
        return Run(path)

    def __repr__(self):
        return f'MDSimulation(model={self.model!r}, particles={self.particles}, step={self.step}, settings={self.settings})'


# Compatibility name for existing notebooks; new examples use MDSimulation.
Simulation = MDSimulation
