"""Interactions, run controls and recording. Each run evolves a copy of a given state."""
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
from .system import System, MODELS, WALLS, _positive, _integer, _readonly
from .walls import Piston, external_terms, piston_config
from .results import Run, write_json
from .species import parameter, per_atom, charges as _charges
from .bonds import bond_config
from .units import units as _units, AVOGADRO_PER_CM3_A3
from ._version import __version__


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


def _seed(value):
    value = _integer(value, 'seed')
    if value >= 2**64:
        raise ValueError('seed must be below 2**64.')
    return value


def _species_values(value, name, zero):
    """Check scalar or per-species parameters before the species are known."""
    if isinstance(value, dict):
        if not value:
            raise ValueError(f'{name} mapping must not be empty.')
        return {k: _positive(v, name, zero=zero) for k, v in value.items()}
    return _positive(value, name, zero=zero)


class MDSimulation:
    """Interactions and default run controls; it holds no particles.

    ``run = sim.run(state, steps=...)`` evolves a copy of ``state`` (a System)
    and saves the result. Continue from ``run.state``. Every run uses the
    controls given here plus its own arguments; nothing carries over between
    runs except through the state you pass in. ``live = sim.start(state)``
    gives the same run in memory, to advance a few steps at a time.

    A state with walls (``boundary="walls"``) is held in by them: an atom that
    has passed a wall by the depth d is pushed back with the force
    ``wall_stiffness*d**2`` (potential energy ``wall_stiffness*d**3/3``), so a
    bounce conserves energy; at the default an atom turns within about a tenth
    of a diameter. The force on each wall is among the observables.
    """

    def __init__(self, system=None, *, pair_potential='lj', epsilon=1., sigma=1., mixing_rule='lorentz-berthelot', exclude_bonded_pairs=True, temperature=None, ensemble='nvt', timestep=None,
                 cutoff=None, G=1., softening=.05, skin=None, friction=None, pressure=None, pressure_time=None,
                 heat_rate=0., seed=87287, output_dir='runs', storage_name='simulation', units='reduced',
                 charge=None, coulomb=None, coulomb_cutoff=None, kspace_accuracy=1e-5, lj_shift=True,
                 wall_stiffness=1e4):
        if system is not None:
            raise TypeError('MDSimulation holds interactions and controls only. Pass the state to each run: '
                            'run = sim.run(system, steps=...), then continue with sim.run(run.state, ...).')
        if pair_potential not in ('lj', 'lj96', 'gravity'):
            raise ValueError('pair_potential must be lj, lj96 or gravity.')
        self._units = _units(units)
        if (charge is None) != (coulomb is None):
            raise ValueError('Charges need coulomb="pppm" (long-range electrostatics), and coulomb needs charge.')
        if coulomb not in (None, 'pppm'):
            raise ValueError('coulomb must be "pppm".')
        real = self._units.name == 'real'
        if (real or charge is not None) and pair_potential != 'lj':
            raise ValueError('Water, charges and real units use the 12-6 Lennard-Jones potential (pair_potential="lj").')
        if not isinstance(lj_shift, bool):
            raise TypeError('lj_shift must be True or False.')
        if charge is None and coulomb_cutoff is not None:
            raise ValueError('coulomb_cutoff applies only to charged systems.')
        if real:
            missing = [name for name, value in (('temperature', temperature), ('timestep', timestep), ('cutoff', cutoff)) if value is None]
            if missing:
                raise ValueError(f'Real units need explicit {", ".join(missing)} (K, fs, Å).')
        defaults = dict(skin=2., friction=.005, pressure=None, pressure_time=1000.) if real else dict(temperature=2., skin=.4, friction=1., pressure=.01, pressure_time=5.)
        temperature = defaults['temperature'] if temperature is None else temperature
        if pair_potential == 'gravity':
            if cutoff is not None:
                raise ValueError('Isolated gravity has no cutoff; use cutoff=None.')
            if epsilon != 1. or sigma != 1.:
                raise ValueError('Gravity uses G and softening, not LJ epsilon/sigma.')
        self._G = _positive(G, 'G')
        self._softening = _positive(softening, 'softening', zero=True)
        if not isinstance(exclude_bonded_pairs, bool):
            raise TypeError('exclude_bonded_pairs must be True or False.')
        if mixing_rule != 'lorentz-berthelot':
            raise ValueError('mixing_rule must be "lorentz-berthelot".')
        zero = real or charge is not None
        self._epsilon = _species_values(epsilon, 'epsilon', zero)
        self._sigma = _species_values(sigma, 'sigma', zero)
        self._charge = None if charge is None else _charges(charge, tuple(charge) if isinstance(charge, dict) else ())
        self._coulomb = coulomb
        self._lj_shift = lj_shift
        self._pair_potential = pair_potential
        self._exclude_bonded_pairs = exclude_bonded_pairs
        self.wall_stiffness = _positive(wall_stiffness, 'wall_stiffness')
        self.cutoff = None if pair_potential == 'gravity' else _positive(2.5 if cutoff is None else cutoff, 'cutoff')
        self.skin = _positive(defaults['skin'] if skin is None else skin, 'skin')
        self.friction = _positive(defaults['friction'] if friction is None else friction, 'friction')
        self.pressure_time = _positive(defaults['pressure_time'] if pressure_time is None else pressure_time, 'pressure_time')
        self.seed = _seed(seed)
        self.output_dir = Path(output_dir).expanduser()
        self.storage_name = _storage_name(storage_name)
        self._defaults = dict(ensemble=ensemble, timestep=timestep, temperature=temperature,
                              pressure=defaults['pressure'] if pressure is None else pressure, heat_rate=heat_rate)
        if charge is not None:
            self._defaults.update(coulomb_cutoff=self.cutoff if coulomb_cutoff is None else coulomb_cutoff,
                                  kspace_accuracy=kspace_accuracy)
        _check_controls(dict(self._defaults), self)

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
    def settings(self):
        """Default run controls; each run may override them for that run only."""
        return self._defaults.copy()

    def _engine(self, state, *, timestep=None, ensemble=None, temperature=None, pressure=None,
                heat_rate=None, coulomb_cutoff=None, kspace_accuracy=None, seed=None, piston=None):
        """Working copy of state for one run, with this run's controls resolved."""
        if not isinstance(state, System):
            raise TypeError('Pass a System as the state, such as FCC(N=500, rho=0.1) or run.state.')
        if self._charge is None and (coulomb_cutoff is not None or kspace_accuracy is not None):
            raise ValueError('coulomb_cutoff and kspace_accuracy apply only to charged systems.')
        settings = self.settings
        for key, value in dict(timestep=timestep, ensemble=ensemble, temperature=temperature, pressure=pressure,
                               heat_rate=heat_rate, coulomb_cutoff=coulomb_cutoff, kspace_accuracy=kspace_accuracy).items():
            if value is not None:
                settings[key] = value
        if state.walls or piston is not None:
            settings['piston'] = piston
        return _Engine(self, state, settings, None if seed is None else _seed(seed))

    def start(self, state, *, timestep=None, ensemble=None, temperature=None, pressure=None,
              heat_rate=None, coulomb_cutoff=None, kspace_accuracy=None, seed=None, piston=None):
        """Start a run of a copy of ``state`` in memory, and return it as a LiveRun.

        Nothing is written to disk. ``live.advance(steps)`` takes steps;
        ``live.observables``, ``live.positions``, ``live.velocities`` and
        ``live.box`` read the result; ``live.set(temperature=...)`` changes a
        control for the steps that follow; ``live.state`` is a System to continue
        from, here or with ``sim.run``. The controls, the clock and the random
        sequence are treated as in ``run``, which is this loop with files:
        advancing in pieces gives the same trajectory as one call.
        ``piston=Piston("right", to=50, speed=0.05)`` moves a wall, as in ``run``.
        """
        return LiveRun(self._engine(state, timestep=timestep, ensemble=ensemble, temperature=temperature,
                                    pressure=pressure, heat_rate=heat_rate, coulomb_cutoff=coulomb_cutoff,
                                    kspace_accuracy=kspace_accuracy, seed=seed, piston=piston))

    @_legacy_output_names
    def run(self, state, steps=1000, *, timestep=None, ensemble=None, temperature=None, pressure=None,
            heat_rate=None, coulomb_cutoff=None, kspace_accuracy=None, deform=None, seed=None, piston=None,
            thermo_every=100, trajectory_every=500, show=False, max_fps=None, frame_every=100, storage_name=None):
        """Evolve a copy of ``state`` and save the run; the final state is ``run.state``.

        Controls given here apply to this run only; omitted ones take the values
        given to MDSimulation. A ``temperature`` or ``pressure`` is the target from
        the first step: a step change, not a ramp. Velocities are not reset; the
        thermostat (and barostat) relax the system towards the new target. A state
        without velocities gets random ones at ``temperature``.

        The clock continues from the state's time and step. The thermostat continues
        the state's random sequence; give ``seed`` for a different one. Results go to
        output_dir/storage_name; reusing a name replaces the old result. Sampling
        and trajectory intervals are in steps; ``trajectory_every=None`` saves no
        frames but keeps thermo data and the final state. ``show=True`` gives a live
        3D view in a notebook; ``max_fps`` optionally paces execution at
        ``frame_every`` MD steps per displayed frame without changing physical dt.

        For charged systems, coulomb_cutoff and kspace_accuracy may be set per run.
        The PPPM mesh is chosen as in LAMMPS and kept in the state; it is chosen
        again when these settings change or the box has changed by more than 5%.
        ``deform=s`` (nve/nvt only) scales every box length linearly in time to s
        times its starting value, moving whole molecules with their centres of mass.

        In a box with walls, the thermo table has for each wall the force on it at
        that step (``wall_force_right``, positive outwards) and the momentum it has
        received so far (``wall_impulse_right``); the mean force between two rows is
        the difference in impulse over the difference in time.
        ``piston=Piston("right", to=50, speed=0.05)`` moves that wall at constant
        speed to x = 50, where it stops; the box follows it, and ``work`` is the work
        the moving wall has done on the atoms.
        """
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
        advancing = self.start(state, timestep=timestep, ensemble=ensemble, temperature=temperature, pressure=pressure,
                               heat_rate=heat_rate, coulomb_cutoff=coulomb_cutoff, kspace_accuracy=kspace_accuracy, seed=seed,
                               piston=piston)
        engine = advancing._engine
        if deform is not None:
            deform = _positive(deform, 'deform')
            if engine.kind != 'molecular':
                raise ValueError('deform is available for water, charges and real units.')
            if engine.settings['ensemble'] not in ('nve', 'nvt'):
                raise ValueError('deform needs ensemble="nve" or "nvt" (no barostat).')
        live = None
        if show:
            from .visualization import LiveView
            live = LiveView(max_fps=max_fps, frame_every=frame_every)
        start = engine.step
        if deform is not None and steps:
            engine._deform = (engine._box.copy(), deform, start, steps)
        path = self.output_dir.resolve()/storage_name
        _replace_run_directory(path)
        meta = dict(schema_version=1, package_version=__version__, label=storage_name,
                    storage_name=storage_name,
                    created_at=datetime.now(timezone.utc).isoformat(), status='running',
                    configuration={**engine._config(), 'storage_name': storage_name}, atoms=len(engine._x),
                    molecules=engine.particles if engine.model != 'atomic' else 0,
                    atom_species=engine.system._species.tolist(),
                    degrees_of_freedom=engine.dof, units=self._units.description, deform=deform,
                    potential_shifted=self.pair_potential != 'gravity' and self.lj_shift, pressure_estimator=('not defined for an unconfined cluster' if self.pair_potential == 'gravity' else 'atomic virial' if engine.model == 'atomic' else 'molecular COM virial'),
                    start_step=start, start_time=engine.time, requested_steps=steps,
                    thermo_every=thermo_every, trajectory_every=trajectory_every,
                    live_view=bool(show), max_fps=max_fps, frame_every=frame_every if show else None)
        write_json(path/'metadata.json', meta)
        started = time.monotonic()
        frame_index = 0
        last_sample = last_frame = -1
        status = 'completed'
        def frame():
            nonlocal frame_index, last_frame
            D = engine.dimension
            np.savez_compressed(path/'frames'/f'{frame_index:07d}.npz',
                positions=engine._x[:, :D], unwrapped=engine._u[:, :D], velocities=engine._v[:, :D],
                box=engine._box[:D], step=engine.step, time=engine.time)
            frame_index += 1
            last_frame = engine.step
        try:
            if live:
                live.start(engine, steps)
            with (path/'thermo.csv').open('w', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=list(engine.observables))
                writer.writeheader()
                writer.writerow(engine.observables)
                last_sample = engine.step
                if trajectory_every is not None:
                    frame()
                while engine.step-start < steps:
                    elapsed = engine.step-start
                    chunk = min(steps-elapsed, thermo_every-elapsed%thermo_every, 256)
                    if trajectory_every is not None:
                        chunk = min(chunk, trajectory_every-elapsed%trajectory_every)
                    if live:
                        chunk = min(chunk, frame_every-elapsed%frame_every)
                    interrupted = advancing._advance(chunk)
                    elapsed = engine.step-start
                    final = elapsed == steps or interrupted
                    if elapsed%thermo_every == 0 or final:
                        writer.writerow(engine.observables)
                        stream.flush()
                        last_sample = engine.step
                    if trajectory_every is not None and (elapsed%trajectory_every == 0 or final):
                        frame()
                    if live and (elapsed%frame_every == 0 or final):
                        live.update(engine, force=final)
                    if interrupted:
                        status = 'interrupted'
                        break
        except KeyboardInterrupt:
            status = 'interrupted'
        except Exception as exc:
            status = 'failed'
            meta['error'] = f'{type(exc).__name__}: {exc}'
            raise
        finally:
            engine._deform = None
            if live:
                live.finish(engine, status)
            if status != 'failed':
                if last_sample != engine.step:
                    with (path/'thermo.csv').open('a', newline='') as stream:
                        csv.DictWriter(stream, fieldnames=list(engine.observables)).writerow(engine.observables)
                if trajectory_every is not None and last_frame != engine.step:
                    frame()
                engine._save_checkpoint(path)
            meta.update(status=status, completed_steps=engine.step-start, end_step=engine.step,
                        end_time=engine.time, wall_seconds=time.monotonic()-started)
            if engine.kind == 'molecular' and self._charge is not None:
                pppm = engine._pppm_state
                meta['pppm'] = dict(mesh=[int(n) for n in pppm[:3]], g_ewald=float(pppm[3]), order=5,
                                    estimated_relative_accuracy=float(pppm[7]), kspace_accuracy=engine.settings['kspace_accuracy'],
                                    coulomb_cutoff=engine.settings['coulomb_cutoff'])
            write_json(path/'metadata.json', meta)
        return Run(path)

    def __repr__(self):
        return f'MDSimulation(pair_potential={self.pair_potential!r}, units={self.units!r}, settings={self.settings})'


class LiveRun:
    """A run held in memory, which its caller advances a few steps at a time.

    ``live = sim.start(state, ...)`` copies the state; nothing is written to
    disk. Between calls to ``advance`` the caller may read the run, change a
    control with ``set``, draw a frame or stop. ``live.state`` is a System that
    ``sim.run`` or another ``sim.start`` continues from, with the clock and the
    thermostat's random sequence carried over. How the steps are divided
    between calls does not change the trajectory.
    """

    def __init__(self, engine):
        self._engine = engine
        self._failed = False

    def _advance(self, steps):
        """Take steps; True if a keyboard interrupt stopped them early."""
        self._check()
        try:
            return self._engine._advance(steps)
        except BaseException:
            self._failed = True
            raise

    def _check(self):
        if self._failed:
            raise RuntimeError('This run stopped on an error and cannot continue; start again from an earlier state.')

    def advance(self, steps):
        """Take ``steps`` MD steps (a nonnegative integer).

        A keyboard interrupt stops after a whole step and is raised again, so a
        loop around ``advance`` ends; the run stays usable.
        """
        if self._advance(_integer(steps, 'steps', 0)):
            raise KeyboardInterrupt

    def set(self, *, timestep=None, ensemble=None, temperature=None, pressure=None, heat_rate=None, piston=None):
        """Change controls for the steps that follow; omitted ones keep their values.

        The effect is that of a new run started from ``live.state`` with these
        controls: a new ``temperature`` or ``pressure`` is the target from the
        next step (a step change, not a ramp), and velocities are not reset.
        ``piston=Piston("right", to=50, speed=0.05)`` sets that wall moving from
        where it is; it stops by itself at ``to``.
        """
        self._check()
        self._engine._set(timestep=timestep, ensemble=ensemble, temperature=temperature,
                          pressure=pressure, heat_rate=heat_rate, piston=piston)

    @property
    def settings(self):
        """The controls now in force (timestep, ensemble, temperature, pressure,
        heat_rate, and piston for a system with walls)."""
        return self._engine.settings.copy()

    @property
    def step(self):
        """Number of MD steps taken to reach the current state, counted from the state's clock."""
        return self._engine.step

    @property
    def time(self):
        """Simulation time of the current state."""
        return self._engine.time

    @property
    def observables(self):
        """Current values of the quantities that ``run`` records in thermo.csv, as a dict.

        With walls it also has, for each wall, ``wall_force_<name>`` (the force
        of the atoms on that wall at this step, positive outwards) and
        ``wall_impulse_<name>`` (the momentum they have delivered to it so far),
        and ``work``, the work moving walls have done on the atoms; with a
        piston, ``piston_position``. The mean force on a wall between two
        readings is the difference in impulse over the difference in time.
        """
        return self._engine.observables

    @property
    def positions(self):
        """Read-only snapshot of the positions, one row per atom: (x, y, z), or
        (x, y) for a flat system. Wrapped into the box along periodic directions."""
        return _readonly(self._engine._x[:, :self._engine.dimension])

    @property
    def velocities(self):
        """Read-only snapshot of the velocities, one row per atom, with the columns of ``positions``."""
        return _readonly(self._engine._v[:, :self._engine.dimension])

    @property
    def box(self):
        """Read-only snapshot of the box lengths (three, or two for a flat system),
        as in a saved frame. A piston changes the length along its axis."""
        return _readonly(self._engine._box[:self._engine.dimension])

    @property
    def state(self):
        """The current state, an independent System: continue from it with ``sim.run`` or ``sim.start``."""
        self._check()
        return self._engine.state()

    def __repr__(self):
        return f'LiveRun(step={self.step}, time={self.time:g}, settings={self.settings})'


def _check_controls(config, sim, kind=None, model=None):
    """Validate run controls; kind and model are known once a state is given."""
    if sim.pair_potential == 'gravity' and config['ensemble'] != 'nve':
        raise ValueError('Isolated gravity uses ensemble="nve" (no thermostat, walls or barostat).')
    if config['ensemble'] not in ('nve', 'nvt', 'nph', 'npt'):
        raise ValueError('ensemble must be nve, nvt, nph or npt.')
    if config['ensemble'] in ('nph', 'npt') and config['pressure'] is None:
        raise ValueError(f'Set pressure ({sim._units.pressure_unit}) for {config["ensemble"]}.')
    for key in ('timestep', 'temperature', 'pressure'):
        if config[key] is not None:
            config[key] = _positive(config[key], key, zero=(key != 'timestep'))
    if 'coulomb_cutoff' in config:
        config['coulomb_cutoff'] = _positive(config['coulomb_cutoff'], 'coulomb_cutoff')
        config['kspace_accuracy'] = _positive(config['kspace_accuracy'], 'kspace_accuracy')
        if config['kspace_accuracy'] >= 1:
            raise ValueError('kspace_accuracy is a relative force accuracy, such as 1e-5.')
    config['heat_rate'] = float(config['heat_rate'])
    if not np.isfinite(config['heat_rate']):
        raise ValueError('heat_rate must be finite (total energy per time).')
    if config['heat_rate'] and config['ensemble'] in ('nvt', 'npt'):
        raise ValueError('Turn off the thermostat (nve or nph) before adding/removing heat.')
    if kind == 'pair' and model != 'atomic' and config['ensemble'] in ('nph', 'npt'):
        raise ValueError('Pressure control is supported for atoms and rigid water only.')
    return config


class _Engine:
    """One run's working copy of a state, with the native solver's arrays."""

    def __init__(self, sim, state, settings, seed):
        self.sim = sim
        self._units = sim._units
        system = state.copy()
        self._system = system
        if len(system.atoms) > 20000:
            raise ValueError('This course solver supports at most 20000 atoms.')
        if system.units != self._units.name:
            raise ValueError(f'This simulation uses units="{self._units.name}", but the state is in "{system.units}" units.')
        real = self._units.name == 'real'
        self.kind = 'molecular' if (system.model == 'water-rigid' or real or sim._charge is not None) else 'pair'
        if self.kind == 'molecular':
            if sim.pair_potential != 'lj':
                raise ValueError('Water, charges and real units use the 12-6 Lennard-Jones potential (pair_potential="lj").')
            if system.model.startswith('diatomic'):
                raise ValueError('Diatomic molecules use reduced units without charges.')
        elif not sim.lj_shift:
            raise ValueError('lj_shift=False is available for water, charges and real units.')
        if sim.pair_potential == 'gravity':
            if system.boundary != 'open' or system.model != 'atomic':
                raise ValueError('Gravity requires an unbonded System with boundary="open".')
            if len(system.atoms) > 4096:
                raise ValueError('Direct all-pairs gravity supports at most 4096 particles.')
        elif system.boundary == 'open':
            raise ValueError('Lennard–Jones needs a box: periodic or with walls. boundary="open" is for gravity.')
        elif (system.dimension == 2 or system.walls) and self.kind != 'pair':
            raise ValueError('A flat system and walls use reduced units and atoms without charges.')
        self.dimension = system.dimension
        zero = self.kind == 'molecular'
        self._epsilon = parameter(sim._epsilon, system.species, 'epsilon', zero=zero)
        self._sigma = parameter(sim._sigma, system.species, 'sigma', zero=zero)
        self._charge = None if sim._charge is None else _charges(sim._charge, system.species)
        self._pair_parameters = np.ascontiguousarray(np.column_stack((per_atom(self._epsilon, system._species), per_atom(self._sigma, system._species))))
        self._x = system._x
        self._box = system._box
        self._mass = system._mass
        self.cutoff = sim.cutoff
        self.skin = sim.skin
        self.friction = sim.friction
        self.pressure_time = sim.pressure_time
        if settings['timestep'] is None:
            settings['timestep'] = .005 if system.model == 'atomic' else .002
        if self.kind == 'pair':
            settings.pop('coulomb_cutoff', None)
            settings.pop('kspace_accuracy', None)
        self.settings = self._controls(settings)
        if self.kind == 'molecular':
            self._prepare_molecular()
        self._deform = None
        # The total momentum is conserved, and kept at zero, along periodic directions only.
        conserved = system._conserved
        if system._v is None:
            rng = np.random.default_rng(sim.seed if seed is None else seed)
            self._v = np.zeros_like(self._x)
            self._v[:, :self.dimension] = rng.normal(size=(len(self._x), self.dimension))/np.sqrt(self._mass[:, None])
            self._v -= np.where(conserved, (self._v*self._mass[:, None]).sum(axis=0)/self._mass.sum(), 0.)
            if system.model == 'diatomic-rigid':
                self._project_velocities()
            if system.model == 'water-rigid':
                _core.project_molecular_velocities(self._x, self._v, self._box, self._mass, self._constraints)
            self._v *= np.sqrt(self.dof*self._units.boltzmann*self.settings['temperature']/(self._units.mvv2e*np.sum(self._mass[:, None]*self._v**2)))
            system._v = self._v
        else:
            self._v = system._v
            momentum = self._v*self._mass[:, None]
            if np.linalg.norm(momentum.sum(axis=0)[conserved]) > 1e-10*max(1., np.linalg.norm(momentum)):
                raise ValueError('Remove centre-of-mass velocity before starting this solver.')
        self._u = self._x.copy() if system._u is None else system._u
        system._u = self._u
        self._f = np.zeros_like(self._x)
        self.step = system.step
        self.time = system.time
        self.heat_added = system._heat
        self._rng = seed if seed is not None else (system._rng if system._rng is not None else sim.seed)
        # One row per wall for the solver, and for each wall the force on it now,
        # the impulse it has received and the work it has done on the atoms.
        self._walls, self._external = external_terms(system, sim.wall_stiffness, self.settings.get('piston'))
        self._wall_sums = np.zeros((len(self._walls), 3)) if self._walls else None
        if self._walls and system._wall_sums is not None:
            self._wall_sums[:, 1:] = system._wall_sums
        self._barostat(system._rate)
        # Kept by the native solver between calls: the neighbour list and the forces.
        self._workspace = _core.molecular_workspace() if self.kind == 'molecular' else _core.workspace()
        self._advance(0)

    def _barostat(self, rate):
        """Barostat rate and inertia for the controls now in force."""
        self._rate = rate if self.settings['ensemble'] in ('nph', 'npt') else 0.
        self._baromass = (self.dof+3)*max(self.settings['temperature'], 1e-6)*self.pressure_time**2
        if self.kind == 'molecular':
            # Barostat inertia from centre-of-mass degrees of freedom (3N-3)+3.
            self._baromass = 3*self.particles*self._units.boltzmann*max(self.settings['temperature'], 1e-6)*self.pressure_time**2

    def _controls(self, settings):
        """Validate run controls for this state."""
        settings = _check_controls(settings, self.sim, self.kind, self.model)
        if (self.dimension == 2 or self._system.walls) and settings['ensemble'] in ('nph', 'npt'):
            raise ValueError('Pressure control needs a three-dimensional periodic box; use nve or nvt.')
        external_terms(self._system, self.sim.wall_stiffness, settings.get('piston'))
        return settings

    def _set(self, **controls):
        """Change run controls as a new run from the current state would have them."""
        settings = self.settings.copy()
        settings.update({key: value for key, value in controls.items() if value is not None})
        self.settings = self._controls(settings)
        self._walls, self._external = external_terms(self._system, self.sim.wall_stiffness, self.settings.get('piston'))
        self._barostat(self._rate)

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
        stored = system._pppm
        if stored is not None and self._charge is not None and stored[1:] == (self.settings['kspace_accuracy'], self.settings['coulomb_cutoff']):
            # Same electrostatics settings: keep the mesh, so continuing matches one longer run.
            self._pppm_state[:] = stored[0]
        water = system.bond if system.model == 'water-rigid' else None
        self._constraints = water.constraints if water else (1., 1., 1.)
        self._alpha = water.alpha if water else 0.

    @property
    def system(self):
        return self._system

    @property
    def model(self):
        return self._system.model

    @property
    def particles(self):
        return self._system.N

    @property
    def dof(self):
        return self._system.dof

    @property
    def atoms(self):
        return self._system.atoms

    @property
    def pair_potential(self):
        return self.sim.pair_potential

    @property
    def sigma(self):
        return self.sim.sigma

    @property
    def units(self):
        return self._units.name

    @property
    def forces(self):
        return _readonly(self._f)

    def _bond_vectors(self):
        dr = self._x[::2]-self._x[1::2]
        return dr-self._box*np.rint(dr/self._box)

    def _project_velocities(self):
        r = self._bond_vectors()
        dot = (r*(self._v[::2]-self._v[1::2])).sum(axis=1)
        inv_mass = 1/self._mass[::2] + 1/self._mass[1::2]
        correction = dot[:, None]*r/((r*r).sum(axis=1)*inv_mass)[:, None]
        self._v[::2] -= correction/self._mass[::2, None]
        self._v[1::2] += correction/self._mass[1::2, None]

    def _advance(self, steps, measure=False):
        c = self.settings
        sim = self.sim
        if sim.pair_potential == 'gravity':
            self._potential, self._virial, done, interrupted = _core.advance_gravity(
                self._x, self._v, self._f, self._mass, steps, c['timestep'],
                sim.G, sim.softening, c['heat_rate'], self._workspace)
            self._u[:] = self._x
            self.step += done
            self.time += done*c['timestep']
            self.heat_added += done*c['timestep']*c['heat_rate']
            self._pressure_numerator = float('nan')
            return bool(interrupted)
        if self.kind == 'molecular':
            return self._advance_molecular(steps)
        # The pressure of molecules costs a pass over the pairs, so it waits until it is read.
        measure = measure or self.model == 'atomic'
        result = _core.advance(self._x, self._v, self._u, self._f, self._box,
            self._mass, steps, c['timestep'], self.cutoff, self.skin,
            c['temperature'], self.friction if c['ensemble'] in ('nvt','npt') else 0.,
            self._rng, MODELS[self.model], c['heat_rate'],
            c['pressure'] if c['ensemble'] in ('nph','npt') else -1., self._baromass, self._rate,
            int(sim.pair_potential == 'lj96'), int(sim.exclude_bonded_pairs), .7, 0., 0., 0., 1., 1., self._pair_parameters, self._system._bond_parameters, self._workspace, int(measure),
            self.dimension, self._external, self._wall_sums)
        self._potential, self._virial, self._rng, done, interrupted, self._rate, numerator = result
        self._pressure_numerator = numerator if measure else None
        self.step += done
        self.time += done*c['timestep']
        self.heat_added += done*c['timestep']*c['heat_rate']
        return bool(interrupted)

    def _advance_molecular(self, steps):
        c = self.settings
        U = self._units
        controlled = c['ensemble'] in ('nph', 'npt')
        deform = None
        if self._deform is not None:
            lengths, scale, start, total = self._deform
            deform = np.array([*lengths, scale, self.step-start, total], dtype=float)
        result = _core.advance_molecular(x=self._x, v=self._v, u=self._u, f=self._f, box=self._box,
            mass=self._mass, types=self._types, charges=self._charges, epsilon=self._epsilon_table,
            sigma=self._sigma_table, pppm_state=self._pppm_state, deform=deform, steps=steps,
            timestep=c['timestep'], cutoff=self.cutoff, coulomb_cutoff=c.get('coulomb_cutoff', self.cutoff), skin=self.skin,
            temperature=U.boltzmann*c['temperature'], friction=self.friction if c['ensemble'] in ('nvt','npt') else 0.,
            seed=self._rng, heat_rate=c['heat_rate'], pressure=c['pressure']/U.pressure if controlled else -1.,
            barostat_mass=self._baromass, rate=self._rate, mvv2e=U.mvv2e, qqr2e=U.qqr2e,
            molecule_size=3 if self.model == 'water-rigid' else 1, constraints=self._constraints,
            alpha=self._alpha, lj_shift=self.sim.lj_shift, coulomb=int(self._charge is not None),
            accuracy=c.get('kspace_accuracy', 1e-5), order=5, workspace=self._workspace)
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
        if self._pressure_numerator is None:
            self._advance(0, measure=True)
        U = self._units
        K = float(.5*np.sum(self._mass[:,None]*self._v**2))*U.mvv2e
        D = self.dimension
        V = float(np.prod(self._box[:D]))  # the area of a flat system
        T = 2*K/(self.dof*U.boltzmann)
        P_internal = self._pressure_numerator/(D*V)
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
        controlled = self.settings['ensemble'] in ('nph','npt')
        H = E+self.settings['pressure']/U.pressure*V if controlled else E+P_internal*V
        B = .5*self._baromass*self._rate**2 if controlled else 0.
        values = dict(step=self.step, time=self.time, temperature=T, pressure=P,
                      kinetic_energy=K, potential_energy=self._potential,
                      total_energy=E, volume=V, density=self.particles/V,
                      atom_density=len(self._x)/V, compressibility_factor=P_internal*V/(self.particles*T*U.boltzmann) if T else float('nan'),
                      enthalpy=H, barostat_energy=B, extended_enthalpy=H+B,
                      heat_added=self.heat_added)
        if U.name == 'real':
            values['mass_density'] = float(self._mass.sum())/(AVOGADRO_PER_CM3_A3*V)
        if self.kind == 'molecular' and self._charge is not None:
            values['lj_energy'] = self._terms['lj']
            values['coulomb_energy'] = self._terms['coulomb']+self._terms['kspace']
        if self._walls:
            for name, (force, impulse, _) in zip(self._walls, self._wall_sums.tolist()):
                values[f'wall_force_{name}'] = force
                values[f'wall_impulse_{name}'] = impulse
            values['work'] = float(self._wall_sums[:, 2].sum())
            piston = self.settings.get('piston')
            if piston is not None:
                values['piston_position'] = float(self._box[WALLS.index(piston.wall)//2])
        return values

    def state(self):
        """Independent System holding this run's current state and clock."""
        state = self._system.copy()
        state._time, state._step = self.time, self.step
        state._rng = self._rng
        state._heat = self.heat_added
        state._rate = self._rate if self.settings['ensemble'] in ('nph', 'npt') else 0.
        if self._walls:
            state._wall_sums = self._wall_sums[:, 1:].copy()
        if self.kind == 'molecular' and self._charge is not None:
            state._pppm = (self._pppm_state.copy(), self.settings['kspace_accuracy'], self.settings['coulomb_cutoff'])
        return state

    def _config(self):
        sim, system = self.sim, self._system
        settings = self.settings.copy()
        if 'piston' in settings:
            settings['piston'] = piston_config(settings['piston'])
        return dict(model=self.model, bond=bond_config(system.bond), units=self.units,
                    charge=sim.charge, coulomb=sim.coulomb, lj_shift=sim.lj_shift,
                    dimension=self.dimension, walls=list(system.walls), wall_stiffness=sim.wall_stiffness,
                    boundary=system.boundary, G=sim.G, softening=sim.softening,
                    pair_potential=sim.pair_potential, epsilon=sim.epsilon, sigma=sim.sigma, mixing_rule=sim.mixing_rule,
                    exclude_bonded_pairs=sim.exclude_bonded_pairs, particles=self.particles,
                    density=system.rho if system.boundary != "open" else None, seed=sim.seed,
                    cutoff=self.cutoff, skin=self.skin, friction=self.friction,
                    pressure_time=self.pressure_time, storage_name=sim.storage_name, **settings)

    def _save_checkpoint(self, path):
        state = self.state()
        D = self.dimension
        with (path/'checkpoint.tmp').open('wb') as stream:
            np.savez_compressed(stream, positions=self._x[:, :D], velocities=self._v[:, :D],
                unwrapped=self._u[:, :D], masses=self._mass, box=self._box[:D], species=self._system._species,
                **({'wall_sums': self._wall_sums[:, 1:]} if self._walls else {}),
                rng=np.uint64(self._rng), step=self.step, time=self.time,
                heat_added=self.heat_added, barostat_rate=state._rate, barostat_mass=self._baromass,
                **({'pppm_state': self._pppm_state, 'pppm_settings': np.array([self.settings['kspace_accuracy'], self.settings['coulomb_cutoff']])}
                   if self.kind == 'molecular' and self._charge is not None else {}))
        (path/'checkpoint.tmp').replace(path/'checkpoint.npz')
