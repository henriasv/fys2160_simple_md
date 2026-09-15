"""Species labels, exact composition and positive per-species parameters."""
from collections.abc import Mapping
import numpy as np
from .system import _positive, _integer


def species_counts(N, species):
    if species is None: species = 'A'
    if isinstance(species, str):
        counts = {species: 500 if N is None else N}
    elif isinstance(species, Mapping) and species:
        counts = dict(species)
    else:
        raise ValueError('species must be a name or a nonempty mapping of names to counts.')
    if any(not isinstance(name, str) or not name.strip() or name != name.strip() for name in counts):
        raise ValueError('Species names must be nonempty strings without surrounding spaces.')
    counts = {name: _integer(count, 'species count') for name, count in counts.items()}
    total = sum(counts.values())
    if total > 20000:
        raise ValueError('This course solver supports at most 20000 atoms.')
    if N is not None and _integer(N, 'N', 2) != total:
        raise ValueError('N must equal the sum of species counts; omit N to infer it.')
    return total, np.array([name for name, count in counts.items() for _ in range(count)])


def parameter(value, names, name):
    if isinstance(value, Mapping):
        if set(value) != set(names):
            raise ValueError(f'{name} must specify exactly these species: {list(names)}.')
        return {s: _positive(value[s], name) for s in names}
    return _positive(value, name)


def per_atom(value, labels):
    return np.array([value[label] for label in labels], dtype=float) if isinstance(value, dict) else np.full(len(labels), value, dtype=float)
