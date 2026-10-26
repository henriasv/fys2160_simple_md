"""Interactive 3D views with bounded data transfer and optional MD pacing."""
from pathlib import Path
import json
import time
import uuid
import numpy as np

# Ball-and-stick display diameters (Å) for water; they do not affect the physics.
WATER_DISPLAY = {'O': .8, 'H': .5}
VIEW_UNITS = {'real': dict(time='ps', time_scale=.001, temperature='K', pressure='atm')}


def _visible_atoms(limit, molecule_size):
    """Largest atom count <= limit that keeps whole molecules."""
    return max(molecule_size, int(limit)//molecule_size*molecule_size)


def _display(configuration, species, particle_radius):
    """Display diameters, explicit bonds and unit labels for the viewer."""
    model = configuration.get('model')
    if particle_radius is not None:
        if isinstance(particle_radius, dict):
            sigma = {k: 2*float(v) for k, v in particle_radius.items()}
            if any(not np.isfinite(v) or v <= 0 for v in sigma.values()):
                raise ValueError('particle_radius values must be positive and finite.')
        else:
            sigma = 2*float(particle_radius)
            if not np.isfinite(sigma) or sigma <= 0:
                raise ValueError('particle_radius must be positive and finite.')
    elif configuration.get('pair_potential') == 'gravity':
        sigma = .07
    elif model == 'water-rigid':
        sigma = {name: WATER_DISPLAY.get(name, 1.) for name in dict.fromkeys(species)}
    else:
        sigma = configuration.get('sigma', 1.)
    extra = {}
    if model == 'water-rigid':
        extra['bonds'] = [[i, i+k] for i in range(0, len(species), 3) for k in (1, 2)]
    if configuration.get('units') in VIEW_UNITS:
        extra['units'] = VIEW_UNITS[configuration['units']]
    return sigma, extra


class LiveView:
    """One persistent widget; only atom snapshots change during a run."""

    def __init__(self, *, max_fps=None, frame_every=100):
        from IPython import get_ipython
        if get_ipython() is None:
            raise ValueError('show=True requires IPython/Jupyter; use Run.view() for saved trajectories.')
        try:
            from ._live_widget import MDWidget
        except ImportError as exc:
            raise ImportError('Live 3D needs the notebook extras: pip install "fys2160-simple-md[notebook]" (or the local package path).') from exc
        self.widget_type = MDWidget
        self.max_fps = max_fps
        self.frame_every = frame_every
        self.last = 0.
        self.last_step = None
        self.widget = None

    @staticmethod
    def snapshot(sim):
        obs = sim.observables
        return dict(x=np.round(sim._x[:LiveView.visible(sim)], 5).tolist(), box=sim._box.tolist(),
                    step=sim.step, time=sim.time, temperature=obs['temperature'],
                    pressure=obs['pressure'] if sim.system.boundary == 'periodic' else None)

    @staticmethod
    def visible(sim):
        return min(len(sim._x), _visible_atoms(3000, 3 if sim.model == 'water-rigid' else 1))

    def start(self, sim):
        from IPython.display import display
        shown = self.visible(sim)
        note = (f'Live 3D · showing {shown} of {len(sim._x)} atoms. '
                'Drag to rotate; two-finger click-drag or Shift-drag to pan; scroll to zoom. ')
        note += (f'{self.frame_every} MD steps/frame; capped at {self.max_fps:g} fps.'
                 if self.max_fps is not None else 'Full simulation speed; display updates capped at 20 fps.')
        species = sim.system._species[:shown].tolist()
        sigma, extra = _display(dict(model=sim.model, pair_potential=sim.pair_potential, sigma=sim.sigma, units=sim.units), species, None)
        self.widget = self.widget_type(
            configuration=dict(projection='perspective', zoom=1., boundary=sim.system.boundary, sigma=sigma, species=species, molecular=sim.model!='atomic', note=note, **extra),
            frame=self.snapshot(sim))
        display(self.widget)
        self.last = time.monotonic()
        self.last_step = sim.step

    def update(self, sim, *, force=False, pace=True):
        if self.widget is None or self.last_step == sim.step:
            return
        now = time.monotonic()
        interval = 1/(self.max_fps if self.max_fps is not None else 20.)
        if self.max_fps is not None and pace:
            # Wait between displayed MD batches; never change the physical dt.
            delay = self.last + interval - now
            if delay > 0:
                time.sleep(delay)
        elif not force and now-self.last < interval:
            return
        self.widget.frame = self.snapshot(sim)
        self.last = time.monotonic()
        self.last_step = sim.step

    def finish(self, sim, status):
        if self.widget is not None:
            if status != 'failed':
                self.update(sim, force=True, pace=False)
            self.widget.status = status


def trajectory_player(run, *, max_frames=150, max_atoms=3000, projection='perspective', zoom=1.0, particle_radius=None, renderer='omovi'):
    from IPython.display import HTML
    zoom = float(zoom)
    if not np.isfinite(zoom) or not .25 <= zoom <= 3:
        raise ValueError('zoom must be between 0.25 and 3.')
    if projection not in ('orthographic', 'perspective'):
        raise ValueError('projection must be orthographic or perspective.')
    if renderer not in ('omovi', 'basic'):
        raise ValueError('renderer must be "omovi" or "basic".')
    if int(max_frames) != max_frames or max_frames < 2 or int(max_atoms) != max_atoms or max_atoms < 2:
        raise ValueError('max_frames and max_atoms must be integers >= 2.')
    paths = sorted((run.path/'frames').glob('*.npz'))
    if not paths:
        raise ValueError('No trajectory frames were saved.')
    config = run.metadata['configuration']
    max_atoms = _visible_atoms(max_atoms, 3 if config.get('model') == 'water-rigid' else 1)
    selected = np.unique(np.linspace(0, len(paths)-1, min(int(max_frames), len(paths))).astype(int))
    frames = []
    for index in selected:
        with np.load(paths[index], allow_pickle=False) as f:
            frames.append(dict(x=np.round(f['positions'][:int(max_atoms)], 5).tolist(),
                               box=f['box'].tolist(), time=float(f['time']), step=int(f['step'])))
    species = run.metadata.get('atom_species', [])[:int(max_atoms)]
    display_sigma, extra = _display(config, species, particle_radius)
    data = dict(boundary=config.get('boundary','periodic'), frames=frames, projection=projection, zoom=zoom, sigma=display_sigma, species=species, molecular=config['model'] != 'atomic', renderer=renderer, **extra)
    identifier = 'md-'+uuid.uuid4().hex
    note = f'{len(frames)} of {len(paths)} saved frames; {len(frames[0]["x"])} of {run.metadata["atoms"]} atoms. Drag to rotate; two-finger click-drag or Shift-drag to pan; scroll over the view to zoom. Arrow keys also rotate the focused view, space plays and pauses, and < and > change the speed. Display sampling does not change saved data.'
    if config.get('boundary') == 'open':
        note += ' Isolated cluster: no walls or periodic images. The view may omit distant particles; the solver retains them.'
    data.update(note=note, live=False)
    script = Path(__file__).with_name('viewer.js').read_text()
    markup = (f'<div id="{identifier}"></div><script>(()=>{{' + script +
              f'createMDView(document.getElementById("{identifier}"),{json.dumps(data)});' +
              '})();</script>')
    return HTML(markup)
