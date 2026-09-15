"""Interactive 3D views with bounded data transfer and optional MD pacing."""
from pathlib import Path
import json
import time
import uuid
import numpy as np


class LiveView:
    """One persistent widget; only atom snapshots change during a run."""

    def __init__(self, *, max_fps=None, frame_every=100):
        from IPython import get_ipython
        if get_ipython() is None:
            raise ValueError('show=True requires IPython/Jupyter; use Run.view() for saved trajectories.')
        try:
            from ._live_widget import MDWidget
        except ImportError as exc:
            raise ImportError('Live 3D needs the notebook extras: pip install "fys2160-md[notebook]" (or the local package path).') from exc
        self.widget_type = MDWidget
        self.max_fps = max_fps
        self.frame_every = frame_every
        self.last = 0.
        self.last_step = None
        self.widget = None

    @staticmethod
    def snapshot(sim):
        obs = sim.observables
        return dict(x=np.round(sim._x[:1500], 5).tolist(), box=sim._box.tolist(),
                    step=sim.step, time=sim.time, temperature=obs['temperature'],
                    pressure=obs['pressure'])

    def start(self, sim):
        from IPython.display import display
        note = (f'Live 3D · showing {min(len(sim._x),1500)} of {len(sim._x)} atoms. '
                'Drag to rotate; scroll to zoom. ')
        note += (f'{self.frame_every} MD steps/frame; capped at {self.max_fps:g} fps.'
                 if self.max_fps is not None else 'Full simulation speed; display updates capped at 20 fps.')
        self.widget = self.widget_type(
            configuration=dict(projection='perspective', zoom=1., molecular=sim.model!='atomic', note=note),
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


def trajectory_player(run, *, max_frames=150, max_atoms=1500, projection='orthographic', zoom=1.0):
    from IPython.display import HTML
    zoom = float(zoom)
    if not np.isfinite(zoom) or not .25 <= zoom <= 3:
        raise ValueError('zoom must be between 0.25 and 3.')
    if projection not in ('orthographic', 'perspective'):
        raise ValueError('projection must be orthographic or perspective.')
    if int(max_frames) != max_frames or max_frames < 2 or int(max_atoms) != max_atoms or max_atoms < 2:
        raise ValueError('max_frames and max_atoms must be integers >= 2.')
    paths = sorted((run.path/'frames').glob('*.npz'))
    if not paths:
        raise ValueError('No trajectory frames were saved.')
    selected = np.unique(np.linspace(0, len(paths)-1, min(int(max_frames), len(paths))).astype(int))
    frames = []
    for index in selected:
        with np.load(paths[index], allow_pickle=False) as f:
            frames.append(dict(x=np.round(f['positions'][:int(max_atoms)], 5).tolist(),
                               box=f['box'].tolist(), time=float(f['time']), step=int(f['step'])))
    data = dict(frames=frames, projection=projection, zoom=zoom, molecular=run.metadata['configuration']['model'] != 'atomic')
    identifier = 'md-'+uuid.uuid4().hex
    note = f'{len(frames)} of {len(paths)} saved frames; {len(frames[0]["x"])} of {run.metadata["atoms"]} atoms. Drag to rotate; scroll over the view to zoom. Arrow keys also rotate the focused view. Display sampling does not change saved data.'
    data.update(note=note, live=False)
    renderer = Path(__file__).with_name('viewer.js').read_text()
    markup = (f'<div id="{identifier}"></div><script>(()=>{{' + renderer +
              f'createMDView(document.getElementById("{identifier}"),{json.dumps(data)});' +
              '})();</script>')
    return HTML(markup)
