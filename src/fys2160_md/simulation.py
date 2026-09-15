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
from .species import parameter, per_atom
from .bonds import bond_config, read_bond

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

    def __init__(self, system, *, pair_potential='lj', epsilon=1., sigma=1., mixing_rule='lorentz-berthelot', exclude_bonded_pairs=True, temperature=2., ensemble='nvt', timestep=None,
                 cutoff=2.5, skin=.4, friction=1., pressure=.01, pressure_time=5.,
                 heat_rate=0., seed=87287, output_dir='runs', storage_name='simulation'):
        if not isinstance(system,System):
            raise TypeError('Pass a System, such as FCC(N=500, rho=0.1).')
        if len(system.atoms)>20000:
            raise ValueError('This course solver supports at most 20000 atoms.')
        if pair_potential not in ('lj', 'lj96'):
            raise ValueError('pair_potential must be "lj" (12–6) or "lj96" (9–6).')
        if not isinstance(exclude_bonded_pairs, bool):
            raise TypeError('exclude_bonded_pairs must be True or False.')
        if mixing_rule != 'lorentz-berthelot':
            raise ValueError('mixing_rule must be "lorentz-berthelot".')
        self._epsilon = parameter(epsilon, system.species, 'epsilon')
        self._sigma = parameter(sigma, system.species, 'sigma')
        self._pair_parameters = np.ascontiguousarray(np.column_stack((per_atom(self._epsilon, system._species), per_atom(self._sigma, system._species))))
        self._pair_potential = pair_potential
        self._exclude_bonded_pairs = exclude_bonded_pairs
        self._system=system.copy()

        self._x=self._system._x
        self._box=self._system._box
        self._mass=self._system._mass
        self.cutoff=_positive(cutoff,'cutoff')
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
        self._validate_settings(self._settings)
        if self._system._v is None:
            rng=np.random.default_rng(self.seed)
            self._v=np.ascontiguousarray(rng.normal(size=self._x.shape)/np.sqrt(self._mass[:,None]))
            self._v-=(self._v*self._mass[:,None]).sum(axis=0)/self._mass.sum()
            if self.model=='diatomic-rigid':self._project_velocities()
            self._v*=np.sqrt(self.dof*self._settings['temperature']/np.sum(self._mass[:,None]*self._v**2))
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
        self._failed=False
        self._advance(0)

    def _validate_settings(self, config):
        if config['ensemble'] not in ('nve', 'nvt', 'nph', 'npt'):
            raise ValueError('ensemble must be nve, nvt, nph or npt.')
        for key in ('timestep', 'temperature', 'pressure'):
            config[key] = _positive(config[key], key, zero=(key != 'timestep'))
        config['heat_rate'] = float(config['heat_rate'])
        if not np.isfinite(config['heat_rate']):
            raise ValueError('heat_rate must be finite (total reduced energy per time).')
        if config['heat_rate'] and config['ensemble'] in ('nvt', 'npt'):
            raise ValueError('Turn off the thermostat (nve or nph) before adding/removing heat.')
        if self.model != 'atomic' and config['ensemble'] in ('nph', 'npt'):
            raise ValueError('Pressure control is supported for the atomic model only.')

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

    @property
    def observables(self):
        K = float(.5*np.sum(self._mass[:,None]*self._v**2))
        V = float(np.prod(self._box))
        T = 2*K/self.dof
        P = self._pressure_numerator/(3*V)
        E = K+self._potential
        controlled = self._settings['ensemble'] in ('nph','npt')
        H = E+self._settings['pressure']*V if controlled else E+P*V
        B = .5*self._baromass*self._rate**2 if controlled else 0.
        return dict(step=self.step, time=self.time, temperature=T, pressure=P,
                    kinetic_energy=K, potential_energy=self._potential,
                    total_energy=E, volume=V, density=self.particles/V,
                    atom_density=len(self._x)/V, compressibility_factor=P*V/(self.particles*T) if T else float('nan'),
                    enthalpy=H, barostat_energy=B, extended_enthalpy=H+B,
                    heat_added=self.heat_added)

    def _config(self):
        return dict(model=self.model, bond=bond_config(self.system.bond),
                    pair_potential=self.pair_potential, epsilon=self.epsilon, sigma=self.sigma, mixing_rule=self.mixing_rule, exclude_bonded_pairs=self.exclude_bonded_pairs, particles=self.particles,
                    density=self.particles/np.prod(self._box), seed=self.seed,
                    cutoff=self.cutoff, skin=self.skin, friction=self.friction,
                    pressure_time=self.pressure_time, storage_name=self.storage_name, **self.settings)

    def _save_checkpoint(self, path):
        with (path/'checkpoint.tmp').open('wb') as stream:
            np.savez_compressed(stream, positions=self._x, velocities=self._v,
                unwrapped=self._u, masses=self._mass, box=self._box, species=self.system._species,
                rng=np.uint64(self._rng), step=self.step, time=self.time,
                heat_added=self.heat_added, barostat_rate=self._rate, barostat_mass=self._baromass)
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
            system=System(data['positions'],data['box'],velocities=data['velocities'],masses=data['masses'],model=c.pop('model'),species=data['species'] if 'species' in data else None,**bond_args)
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
        obj._advance(0)
        return obj

    @_legacy_output_names
    def run(self, timestep=None, steps=1000, *, ensemble=None, temperature=None,
            heat_rate=None, pressure=None, thermo_every=100, trajectory_every=500,
            show=False, max_fps=None, frame_every=100, storage_name=None):
        """Advance and save a run; omitted physical settings retain their values.

        Results go to output_dir/storage_name. Reusing a name replaces the old
        result, including its trajectory. The constructor's name is the default.
        Sampling and trajectory intervals are in steps. ``trajectory_every=None``
        disables trajectory frames, but preserves thermo and a final checkpoint.
        ``show=True`` gives an interactive live 3D view in a Jupyter notebook.
        ``max_fps`` optionally paces execution at at most that many displayed
        frames/second, with ``frame_every`` MD steps/frame. None runs at full
        speed with display-only throttling. Pacing does not change physical dt.
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
        for key, value in dict(timestep=timestep, ensemble=ensemble, temperature=temperature,
                               heat_rate=heat_rate, pressure=pressure).items():
            if value is not None:
                config[key] = value
        self._validate_settings(config)
        live = None
        if show:
            from .visualization import LiveView
            live = LiveView(max_fps=max_fps, frame_every=frame_every)
        was_controlled = self._settings['ensemble'] in ('nph','npt')
        self._settings = config
        if not was_controlled and config['ensemble'] in ('nph','npt'):
            self._rate = 0.
        self._advance(0)
        start = self.step
        path = self.output_dir.resolve()/storage_name
        _replace_run_directory(path)
        meta = dict(schema_version=1, package_version='0.1.0', label=storage_name,
                    storage_name=storage_name,
                    created_at=datetime.now(timezone.utc).isoformat(), status='running',
                    configuration={**self._config(), 'storage_name': storage_name}, atoms=len(self._x),
                    molecules=self.particles if self.model != 'atomic' else 0,
                    atom_species=self.system._species.tolist(),
                    degrees_of_freedom=self.dof, units='fixed reduced reference units: length=energy=mass=kB=1',
                    potential_shifted=True, pressure_estimator='atomic virial' if self.model == 'atomic' else 'molecular COM virial',
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
            write_json(path/'metadata.json', meta)
        return Run(path)

    def __repr__(self):
        return f'MDSimulation(model={self.model!r}, particles={self.particles}, step={self.step}, settings={self.settings})'


# Compatibility name for existing notebooks; new examples use MDSimulation.
Simulation = MDSimulation
