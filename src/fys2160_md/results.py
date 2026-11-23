"""Portable, engine-independent saved results. Reading a Run never starts a solver."""
from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:  # pandas is loaded only when Run.thermo is read.
    import pandas as pd


def _readonly(array):
    array.flags.writeable = False
    return array


def write_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


class Run:
    """One saved call to MDSimulation.run, including its configuration and samples."""

    def __init__(self, path):
        self.path = Path(path).expanduser().resolve()
        self.metadata = json.loads((self.path / 'metadata.json').read_text())
        if self.metadata.get('schema_version') != 1:
            raise ValueError('Unsupported run format.')

    @classmethod
    def load(cls, path) -> 'Run':
        return cls(path)

    @classmethod
    def find(cls, directory='runs') -> list['Run']:
        return [cls(p.parent) for p in sorted(Path(directory).glob('*/metadata.json'))]

    @property
    def state(self):
        """Final state of this run, a System: pass it to the next sim.run.

        It holds positions, velocities, box, masses and species, plus the clock
        and thermostat/barostat state, so continuing from it matches one longer run.
        """
        # One state per Run, so that a change such as run.state.box[0] *= 2 is
        # kept for the next sim.run(run.state).
        if getattr(self, '_state', None) is not None:
            return self._state
        from .system import System
        from .bonds import read_bond
        if self.metadata.get('status') not in ('completed', 'interrupted') or not (self.path/'checkpoint.npz').exists():
            raise ValueError('This run did not finish, so it has no final state.')
        c = self.metadata['configuration']
        bond = {'bond': read_bond(c['bond'])} if 'bond' in c else {}
        with np.load(self.path/'checkpoint.npz', allow_pickle=False) as data:
            boundary = c.get('boundary', 'periodic')
            state = System(data['positions'], data['box'], velocities=data['velocities'], masses=data['masses'],
                           model=c['model'], boundary=boundary if isinstance(boundary, str) else tuple(boundary),
                           species=data['species'] if 'species' in data else None, units=c.get('units', 'reduced'),
                           unwrapped=data['unwrapped'], time=float(data['time']), step=int(data['step']),
                           tilt=data['tilt'] if 'tilt' in data else None, **bond)
            state._rng = int(data['rng'])
            state._heat = float(data['heat_added'])
            rate = data['barostat_rate']
            state._rate = float(rate) if rate.ndim == 0 else rate.copy()  # one per axis if anisotropic
            if 'wall_sums' in data:
                state._wall_sums = data['wall_sums'].copy()
            if 'pppm_settings' in data:
                accuracy, cutoff = (float(v) for v in data['pppm_settings'])
                state._pppm = (data['pppm_state'].copy(), accuracy, cutoff)
        self._state = state
        return state

    @property
    def thermo(self) -> pd.DataFrame:
        """Recorded observables for the whole system, in the run's units (see metadata['units'])."""
        try:
            import pandas as pd
        except ImportError as exc:
            raise ImportError(
                'Run.thermo returns a pandas DataFrame, and pandas is not installed. Install pandas '
                '(in a browser page: await micropip.install("pandas")), or read the same table with '
                "numpy.genfromtxt(run.path/'thermo.csv', delimiter=',', names=True).") from exc
        return pd.read_csv(self.path / 'thermo.csv')

    def iter_frames(self):
        """Yield saved frames in time order, loading only one at a time."""
        for path in sorted((self.path / 'frames').glob('*.npz')):
            with np.load(path, allow_pickle=False) as frame:
                yield {k: frame[k].copy() for k in frame.files}

    def load_frames(self) -> dict[str, np.ndarray]:
        """Load and stack all saved frames. Use iter_frames for large runs."""
        frames = list(self.iter_frames())
        if not frames:
            raise ValueError('This run contains no saved trajectory frames.')
        return {key: np.stack([frame[key] for frame in frames]) for key in frames[0]}

    def _sampled(self, key):
        """One stacked frame array, loaded with the others on first use and kept."""
        if getattr(self, '_frames', None) is None:
            paths = sorted((self.path / 'frames').glob('*.npz'))
            if not paths:
                raise ValueError('This run saved no frames: run it with trajectory_every set.')
            with np.load(paths[0], allow_pickle=False) as first:
                size = len(paths) * sum(first[k].nbytes for k in first.files)
            if size > 2e9:
                raise MemoryError(f'The saved frames would take {size/1e9:.1f} GB in memory. Read them one at '
                                  'a time with iter_frames(), or save fewer (a larger trajectory_every).')
            self._frames = {key: _readonly(value) for key, value in self.load_frames().items()}
        return self._frames[key]

    @property
    def positions(self) -> np.ndarray:
        """Saved positions, one row per saved step: (saved steps, atoms, dimensions)."""
        return self._sampled('positions')

    @property
    def velocities(self) -> np.ndarray:
        """Saved velocities: (saved steps, atoms, dimensions)."""
        return self._sampled('velocities')

    @property
    def unwrapped(self) -> np.ndarray:
        """Saved positions not wrapped into the box, for displacements: (saved steps, atoms, dimensions)."""
        return self._sampled('unwrapped')

    @property
    def boxes(self) -> np.ndarray:
        """Box lengths at each saved step: (saved steps, dimensions)."""
        return self._sampled('box')

    @property
    def tilts(self) -> np.ndarray:
        """Tilt factors (xy, xz, yz) of a tilted box at each saved step; zeros for a rectangular box."""
        try:
            return self._sampled('tilt')
        except KeyError:
            return np.zeros((len(self.steps), 3))

    @property
    def steps(self) -> np.ndarray:
        """Step number of each saved frame."""
        return self._sampled('step')

    @property
    def times(self) -> np.ndarray:
        """Simulation time of each saved frame."""
        return self._sampled('time')

    def draw(self, frame=-1, ax=None, *, view=None, projection='perspective', radius=None, arrows=True):
        """Draw one saved frame as balls in a matplotlib figure and return its Axes.

        ``frame`` indexes the saved frames (``run.steps``), the last by default.
        Views, sizes and colours as in ``System.draw``.
        """
        from .drawing import draw_atoms, display_radii
        from .system import _walled
        from .visualization import _display
        config = self.metadata['configuration']
        species = np.array(self.metadata['atom_species'])
        positions, box = self.positions[frame], self.boxes[frame]
        D = positions.shape[1]
        model = config['model']
        if radius is None:
            sigma, _ = _display(config, species.tolist(), None)
            radii = (np.array([sigma[name] for name in species]) if isinstance(sigma, dict)
                     else np.full(len(species), float(sigma)))/2
        else:
            radii = display_radii(species, model, radius, config.get('units', 'reduced'))
        if model == 'water-rigid':
            bonds = [(i, i+k) for i in range(0, len(species), 3) for k in (1, 2)]
        elif model.startswith('diatomic'):
            bonds = [(i, i+1) for i in range(0, len(species), 2)]
        else:
            bonds = []
        boundary = config.get('boundary', 'periodic')
        open_space = boundary == 'open'
        walled = (False,)*D if open_space else _walled(boundary if isinstance(boundary, str) else tuple(boundary), D)
        return draw_atoms(positions, box, species, radii, bonds=bonds, walled=walled, framed=not open_space,
                          ax=ax, view=view, projection=projection, arrows=arrows,
                          tilt=self.tilts[frame] if D == 3 else None)

    @property
    def measured(self) -> dict[str, np.ndarray]:
        """The run's measurements (``measure=``), one value per step from the start; empty without them."""
        if not (self.path / 'measure.npz').exists():
            return {}
        with np.load(self.path / 'measure.npz', allow_pickle=False) as data:
            return {key: data[key].copy() for key in data.files if key != 'steps'}

    @property
    def measured_steps(self) -> np.ndarray:
        """Step numbers of the values in ``measured``."""
        if not (self.path / 'measure.npz').exists():
            return np.empty(0, dtype=int)
        with np.load(self.path / 'measure.npz', allow_pickle=False) as data:
            return data['steps'].copy()

    def plot(self, *quantities, x='time'):
        """Plot one quantity per panel, with no twin axes."""
        import matplotlib.pyplot as plt
        quantities = quantities or ('temperature', 'pressure')
        data = self.thermo
        fig, axes = plt.subplots(len(quantities), 1, sharex=True,
                                 squeeze=False, figsize=(8, 2.5 * len(quantities)))
        for ax, quantity in zip(axes[:, 0], quantities):
            ax.plot(data[x], data[quantity])
            ax.set_ylabel(quantity.replace('_', ' '))
            ax.grid(alpha=.2)
        axes[-1, 0].set_xlabel(x.replace('_', ' '))
        fig.tight_layout()
        return fig

    def view(self, *, max_frames=150, max_atoms=3000, projection='perspective', zoom=1.0, particle_radius=None, renderer='omovi'):
        """Return a notebook trajectory player with drag-to-rotate.

        renderer="omovi" draws with OMOVI (loaded from jsDelivr, with ambient
        occlusion in perspective view); without a network connection, or with
        renderer="basic", the built-in renderer is used. particle_radius sets
        display sizes only: one value, or one per species."""
        from .visualization import trajectory_player
        return trajectory_player(self, max_frames=max_frames, max_atoms=max_atoms, projection=projection, zoom=zoom,
                                 particle_radius=particle_radius, renderer=renderer)

    def __repr__(self):
        return (f"Run(storage_name={self.metadata.get('storage_name', self.metadata.get('label'))!r}, "
                f"status={self.metadata.get('status')!r}, path={str(self.path)!r})")
