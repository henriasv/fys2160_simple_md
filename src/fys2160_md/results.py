"""Portable, engine-independent saved results. Reading a Run never starts a solver."""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd


def write_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


class Run:
    """One saved call to Simulation.run, including its configuration and samples."""

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
    def thermo(self) -> pd.DataFrame:
        """Recorded observables. Energies are totals, time is reduced MD time."""
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

    def view(self, *, max_frames=150, max_atoms=1500, projection='orthographic', zoom=1.0):
        """Return a self-contained notebook trajectory player with drag-to-rotate."""
        from .visualization import trajectory_player
        return trajectory_player(self, max_frames=max_frames, max_atoms=max_atoms, projection=projection, zoom=zoom)

    def __repr__(self):
        return (f"Run(storage_name={self.metadata.get('storage_name', self.metadata.get('label'))!r}, "
                f"status={self.metadata.get('status')!r}, path={str(self.path)!r})")
