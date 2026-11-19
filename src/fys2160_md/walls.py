"""Walls and pistons: the one-body potentials the solver knows, as it takes them."""
from __future__ import annotations
from dataclasses import dataclass, asdict
import numpy as np
from .system import WALLS, _positive


@dataclass(frozen=True)
class Piston:
    """A wall that moves: ``wall`` goes to the position ``to`` at ``speed``, then stops.

    ``wall`` is "right", "top" or "front": the wall at the far end of the x, y
    or z axis, whose position is that box length. It starts where it is and
    moves at constant speed, inwards or outwards, until it is at ``to``.
    Give it to a run, ``sim.run(state, piston=Piston("right", to=50, speed=0.05))``,
    to ``sim.start`` or to ``live.set``. The box length follows the wall, and the
    work it does on the atoms is the observable ``work``.
    """
    wall: str
    to: float
    speed: float

    def __post_init__(self):
        if self.wall not in WALLS[1::2]:
            raise ValueError('A piston is the wall "right", "top" or "front" (the far end of x, y or z); '
                             'the walls at 0 stay where they are.')
        object.__setattr__(self, 'to', _positive(self.to, 'to'))
        object.__setattr__(self, 'speed', _positive(self.speed, 'speed'))


def piston_config(piston):
    """A piston as it is written to metadata.json."""
    return None if piston is None else asdict(piston)


def external_terms(system, stiffness, piston):
    """The solver's table of one-body terms for a state, and their names.

    One row per wall: kind (0 at coordinate 0, 1 at the box length), axis,
    stiffness, and the speed and target of a piston (0 for a wall at rest).
    A later kind of term (gravity, a bias, a step) is another row here and
    another case in external_forces in _core.c.
    """
    names = system.walls
    if piston is not None and not isinstance(piston, Piston):
        raise TypeError('piston must be a Piston, such as Piston("right", to=50, speed=0.05).')
    if piston is not None and piston.wall not in names:
        raise ValueError(f'This system has no wall "{piston.wall}" to move; its walls are {names or "none (boundary is periodic)"}.')
    rows = []
    for name in names:
        index = WALLS.index(name)
        moving = piston is not None and piston.wall == name
        rows.append([index % 2, index//2, stiffness, piston.speed if moving else 0., piston.to if moving else 0.])
    return names, (np.ascontiguousarray(rows, dtype=float) if rows else None)
