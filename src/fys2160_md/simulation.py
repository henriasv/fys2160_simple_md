"""Small, explicit Python interface to the course's native C integrator."""
from __future__ import annotations
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import csv
import json
import shutil
import time
import numpy as np
from . import _core
from .results import Run, write_json

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


def _readonly(array):
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
    velocities: np.ndarray
    masses: np.ndarray
    molecule_ids: np.ndarray

    def __len__(self):
        return len(self.positions)


class Simulation:
    """Three-dimensional periodic LJ dynamics, in reduced units.

    ``particles`` counts atoms for the atomic model and molecules for diatomics.
    ``density`` uses that same count. Omitted run settings retain their values.
    Each run replaces the saved result with the same storage_name. See README for model conventions.
    """

    def __init__(self, *, model='atomic', particles=500, density=.001,
                 temperature=2., ensemble='nvt', timestep=None, cutoff=2.5,
                 skin=.4, friction=1., pressure=.01, pressure_time=5.,
                 heat_rate=0., seed=87287, output_dir='runs', storage_name='simulation'):
        if model not in MODELS:
            raise ValueError(f'model must be one of {tuple(MODELS)}')
        self.model = model
        self.particles = _integer(particles, 'particles', 2)
        if self.particles * (1 if model == 'atomic' else 2) > 20000:
            raise ValueError('This course solver supports at most 20000 atoms.')
        self.cutoff = _positive(cutoff, 'cutoff')
        self.skin = _positive(skin, 'skin')
        self.friction = _positive(friction, 'friction')
        self.pressure_time = _positive(pressure_time, 'pressure_time')
        self.seed = _integer(seed, 'seed', 1)
        if self.seed >= 2**64:
            raise ValueError('seed must be below 2**64.')
        self.output_dir = Path(output_dir).expanduser()
        self.storage_name = _storage_name(storage_name)
        self._settings = dict(ensemble=ensemble,
                              timestep=(.005 if model == 'atomic' else .002) if timestep is None else timestep,
                              temperature=temperature, pressure=pressure, heat_rate=heat_rate)
        self._validate_settings(self._settings)
        density = _positive(density, 'density')
        length = (self.particles / density)**(1/3)
        if length <= 2*self.cutoff:
            raise ValueError('Box must be larger than twice the cutoff; increase particles or lower density.')
        self._box = np.full(3, length, dtype=float)
        rng = np.random.default_rng(self.seed)
        cells = int(np.ceil((self.particles/4)**(1/3)))
        while 4*cells**3 < self.particles:
            cells += 1
        bases = np.array([[0,0,0], [0,.5,.5], [.5,0,.5], [.5,.5,0]])
        sites = (np.indices((cells,)*3).reshape(3,-1).T[:,None,:] + bases).reshape(-1,3)
        if len(sites) > self.particles:
            sites = sites[rng.choice(len(sites), self.particles, replace=False)]
        centers = (sites + .25)*length/cells
        if model == 'atomic':
            self._x = np.ascontiguousarray(centers)
        else:
            directions = rng.normal(size=(self.particles,3))
            directions /= np.linalg.norm(directions, axis=1)[:,None]
            self._x = np.ascontiguousarray(np.stack((centers-.35*directions, centers+.35*directions), axis=1).reshape(-1,3) % length)
        self._mass = np.ones(len(self._x))
        self._v = np.ascontiguousarray(rng.normal(size=self._x.shape))
        self._v -= self._v.mean(axis=0)
        if model == 'diatomic-rigid':
            self._project_velocities()
        self._v *= np.sqrt(self.dof*self._settings['temperature']/(self._v*self._v).sum())
        self._u = self._x.copy()
        self._f = np.zeros_like(self._x)
        self.step = 0
        self.time = 0.
        self.heat_added = 0.
        self._rng = self.seed
        self._rate = 0.
        self._baromass = (self.dof+3)*max(self._settings['temperature'], 1.e-6)*self.pressure_time**2
        self._failed = False
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

    @classmethod
    def from_arrays(cls, positions, velocities, box, *, masses=None, model='atomic', **kwargs):
        """Start from explicit atom arrays. Diatomic partners must be consecutive.

        Total momentum must be zero: temperature uses the COM-removed DOF.
        Rigid bonds must already have length 0.7 and tangent velocities.
        """
        x = np.array(positions, dtype=float, order='C', copy=True)
        v = np.array(velocities, dtype=float, order='C', copy=True)
        lengths = np.broadcast_to(np.asarray(box, dtype=float), (3,)).copy()
        if x.ndim != 2 or x.shape[1] != 3 or v.shape != x.shape or not np.isfinite(x).all() or not np.isfinite(v).all():
            raise ValueError('positions and velocities must be finite (N,3) arrays.')
        count = len(x) if model == 'atomic' else len(x)//2
        if model != 'atomic' and len(x)%2:
            raise ValueError('Diatomic models require an even number of atoms.')
        obj = cls(model=model, particles=count, density=count/np.prod(lengths), **kwargs)
        obj._box = lengths
        obj._x = x % lengths
        obj._u = x.copy()
        obj._v = v
        obj._f = np.zeros_like(x)
        obj._mass = np.ones(len(x)) if masses is None else np.array(masses, dtype=float, order='C', copy=True)
        if obj._mass.shape != (len(x),) or not np.isfinite(obj._mass).all() or np.any(obj._mass <= 0):
            raise ValueError('masses must be a positive finite (N,) array.')
        if np.linalg.norm((v*obj._mass[:,None]).sum(axis=0)) > 1.e-10*max(1., np.linalg.norm(v*obj._mass[:,None])):
            raise ValueError('Remove centre-of-mass velocity before using from_arrays.')
        if model == 'diatomic-rigid':
            bond = obj._bond_vectors()
            if not np.allclose(np.linalg.norm(bond, axis=1), .7, atol=1e-10, rtol=0):
                raise ValueError('Rigid bonds must initially have length 0.7.')
            if np.max(np.abs((bond*(v[::2]-v[1::2])).sum(axis=1))) > 1.e-10:
                raise ValueError('Rigid bond relative velocities must be perpendicular to the bonds.')
        obj._advance(0)
        return obj

    def _bond_vectors(self):
        dr = self._x[::2]-self._x[1::2]
        return dr-self._box*np.rint(dr/self._box)

    def _project_velocities(self):
        r = self._bond_vectors()
        dot = (r*(self._v[::2]-self._v[1::2])).sum(axis=1)
        correction = dot[:,None]*r/(r*r).sum(axis=1)[:,None]/2
        self._v[::2] -= correction
        self._v[1::2] += correction

    @property
    def dof(self):
        return 3*len(self._x)-3-(self.particles if self.model == 'diatomic-rigid' else 0)

    @property
    def box(self):
        return Box(_readonly(self._box))

    @property
    def atoms(self):
        ids = np.arange(len(self._x)) // (1 if self.model == 'atomic' else 2)
        return Atoms(_readonly(self._x), _readonly(self._v), _readonly(self._mass), _readonly(ids))

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
            c['pressure'] if c['ensemble'] in ('nph','npt') else -1., self._baromass, self._rate)
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
        return dict(model=self.model, particles=self.particles,
                    density=self.particles/np.prod(self._box), seed=self.seed,
                    cutoff=self.cutoff, skin=self.skin, friction=self.friction,
                    pressure_time=self.pressure_time, storage_name=self.storage_name, **self.settings)

    def _save_checkpoint(self, path):
        with (path/'checkpoint.tmp').open('wb') as stream:
            np.savez_compressed(stream, positions=self._x, velocities=self._v,
                unwrapped=self._u, masses=self._mass, box=self._box,
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
        obj = cls(**c)
        with np.load(run.path/'checkpoint.npz', allow_pickle=False) as data:
            for attr, key in (('_x','positions'),('_v','velocities'),('_u','unwrapped'),('_mass','masses'),('_box','box')):
                setattr(obj, attr, data[key].copy())
            obj._f = np.zeros_like(obj._x)
            obj._rng = int(data['rng'])
            obj.step = int(data['step'])
            obj.time = float(data['time'])
            obj.heat_added = float(data['heat_added'])
            obj._rate = float(data['barostat_rate'])
            obj._baromass = float(data['barostat_mass'])
        obj._advance(0)
        return obj

    def run(self, timestep=None, steps=1000, *, ensemble=None, temperature=None,
            heat_rate=None, pressure=None, sample_every=100, save_every=500,
            show=False, max_fps=None, frame_every=100, storage_name=None):
        """Advance and save a run; omitted physical settings retain their values.

        Results go to output_dir/storage_name. Reusing a name replaces the old
        result, including its trajectory. The constructor's name is the default.
        Sampling and trajectory intervals are in steps. ``save_every=None``
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
        sample_every = _integer(sample_every, 'sample_every')
        if save_every is not None:
            save_every = _integer(save_every, 'save_every')
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
                    degrees_of_freedom=self.dof, units='LJ reduced: sigma=epsilon=reference_mass=kB=1',
                    potential_shifted=True, pressure_estimator='atomic virial' if self.model == 'atomic' else 'molecular COM virial',
                    start_step=start, start_time=self.time, requested_steps=steps,
                    sample_every=sample_every, save_every=save_every,
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
                if save_every is not None:
                    frame()
                while self.step-start < steps:
                    elapsed = self.step-start
                    chunk = min(steps-elapsed, sample_every-elapsed%sample_every, 256)
                    if save_every is not None:
                        chunk = min(chunk, save_every-elapsed%save_every)
                    if live:
                        chunk = min(chunk, frame_every-elapsed%frame_every)
                    interrupted = self._advance(chunk)
                    elapsed = self.step-start
                    final = elapsed == steps or interrupted
                    if elapsed%sample_every == 0 or final:
                        writer.writerow(self.observables)
                        stream.flush()
                        last_sample = self.step
                    if save_every is not None and (elapsed%save_every == 0 or final):
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
                if save_every is not None and last_frame != self.step:
                    frame()
                self._save_checkpoint(path)
            meta.update(status=status, completed_steps=self.step-start, end_step=self.step,
                        end_time=self.time, wall_seconds=time.monotonic()-started)
            write_json(path/'metadata.json', meta)
        return Run(path)

    def __repr__(self):
        return f'Simulation(model={self.model!r}, particles={self.particles}, step={self.step}, settings={self.settings})'
