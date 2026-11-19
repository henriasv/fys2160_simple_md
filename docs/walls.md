# Flat boxes, walls and pistons

Atoms in reduced Lennard-Jones units can move in **two dimensions**, be held in
by **walls** instead of periodic images, and be compressed by a **piston**. The
force on every wall is recorded, so the pressure can be read as the momentum
the atoms deliver to a wall.

```python
from fys2160_md import Random, MDSimulation, Piston

gas = Random(N=200, box=(100, 100), boundary="walls", min_distance=1.0)   # two box lengths: a flat box
sim = MDSimulation(pair_potential="lj", ensemble="nvt", temperature=2.0,
                   timestep=0.005, cutoff=2.5, friction=0.05)
live = sim.start(gas)
live.advance(20_000)                                   # settle with the thermostat: 100 time units
before = live.observables
live.advance(400_000)                                  # 2000 time units
after = live.observables
force = (after["wall_impulse_right"] - before["wall_impulse_right"]) / (after["time"] - before["time"])
print("force on the right wall", round(force, 2), "  N k T / Lx =", 200 * 2.0 / 100)

live.set(ensemble="nve", piston=Piston("right", to=50, speed=0.05))       # thermostat off; the wall moves in
energy, T = after["total_energy"], after["temperature"]
live.advance(200_000)                                  # the stroke: 50 / 0.05 = 1000 time units
after = live.observables
print("T", round(T, 2), "->", round(after["temperature"], 2), "  box", live.box)
print("work", round(after["work"], 1), "  change of energy", round(after["total_energy"] - energy, 1))
```

```text
force on the right wall 4.0   N k T / Lx = 4.0
T 1.76 -> 3.63   box [ 50. 100.]
work 368.0   change of energy 367.9
```

The same works with saved runs: `sim.run(gas, steps=..., piston=...)` writes the
wall columns to `thermo.csv`, and `run.state` is the compressed gas.

## A flat system

A box with two lengths is two-dimensional: `Random(N=200, box=(100, 100))`, or
`Random(N=200, rho=0.02, dimension=2)` for a square at that number of atoms per
area. With explicit arrays, `System(positions, box)` is flat when the positions
have two columns.

| | Three dimensions | Flat |
|---|---|---|
| `state.dimension` | 3 | 2 |
| `state.atoms.positions`, `.velocities`, `live.positions`, `live.velocities`, saved frames | (atoms, 3) | (atoms, 2) |
| `state.box.lengths`, `live.box` | three lengths | two lengths |
| `volume`, `state.box.volume` | Lx Ly Lz | the area Lx Ly |
| `density`, `state.rho` | atoms per volume | atoms per area |
| `pressure` | force per area | force per length |
| Kinetic energy per atom at temperature T | 3T/2 | T |

Flat systems take unbonded atoms with `pair_potential="lj"` or `"lj96"`, in
`nve` or `nvt`. `FCC` is a three-dimensional lattice and rejects a flat box.
`run.view()` and `show=True` draw a flat system face-on.

## Walls

`boundary="walls"` on `Random`, `FCC` or `System` puts a wall on every face of
the box in place of the periodic images. One word per axis chooses for each
pair of faces: `boundary=("walls", "periodic")` is a channel with walls across
x. The walls of a state are listed in `state.walls`:

| Wall | Where |
|---|---|
| `"left"`, `"right"` | x = 0 and x = Lx |
| `"bottom"`, `"top"` | y = 0 and y = Ly |
| `"back"`, `"front"` | z = 0 and z = Lz (three dimensions only) |

A wall is a short-range repulsive potential along its normal. An atom that has
passed the wall by the depth d is pushed back with the force
`wall_stiffness * d**2` and has the potential energy `wall_stiffness * d**3 / 3`,
which is part of `potential_energy`. A bounce therefore conserves energy, which
a mirror reflection in a finite time step does not. `wall_stiffness` is an
argument of `MDSimulation` (default 10 000 in reduced units); a thermal atom at
T* = 2 turns within about 0.08 of the wall, so the walls are at the stated box
edges for the atom centres to within a tenth of a diameter.

The force starts from zero without a kink at the wall. That matters for the
integration: with a spring wall (force proportional to d) of the same hardness,
every time step that straddles the wall costs energy accuracy, and the total
energy of 200 atoms wandered about ten times more.

For each wall the observables (and `thermo.csv`) have:

| Observable | Meaning |
|---|---|
| `wall_force_right` | Force of the atoms on that wall at this step, positive outwards. Zero most of the time in a dilute gas |
| `wall_impulse_right` | Momentum the atoms have delivered to that wall so far: the time integral of the force, along all the runs that led to this state |

The **mean force** on a wall between two readings is the difference in impulse
divided by the difference in time. Taken over a short interval it shows single
bounces; over a long one it is the pressure times the length (or area) of the
wall. The impulse is exact bookkeeping: in a run without a thermostat the
change in the gas's total momentum equals the impulses of the walls.

`pressure` is still computed from the momentum of the atoms and the forces
between them, as in a periodic box, and agrees with the force on the walls per
length on average.

With walls the total momentum along that axis is not conserved, so it is not
removed and not required to be zero, and the temperature counts every
direction of every atom: `T = 2K / (2N)` in a flat box with walls all round.
Along a periodic axis the total momentum stays zero and is left out of the
count, as before.

`state.with_box((200, 100))` moves the far wall and leaves every atom where it
is: the gas of a 100 x 100 box with twice the room, as when a partition is
removed. A state taken in the middle of a run may hold an atom in the middle of
a bounce, slightly beyond a wall; that is accepted.

## A piston

```python
Piston(wall, to, speed)
```

A piston is a wall with a schedule: `wall` moves from where it is to the
position `to` at constant `speed`, and stops there. It is a control of a run:
`sim.run(state, piston=...)`, `sim.start(state, piston=...)` or
`live.set(piston=...)`. The piston can be `"right"`, `"top"` or `"front"`, the
wall at the far end of its axis, whose position is the box length; the walls at
0 stay where they are. `to` may be smaller (compression) or larger (expansion).

| Observable | Meaning |
|---|---|
| `piston_position` | Position of the moving wall, equal to the box length along its axis |
| `volume`, `density` | Follow the piston |
| `work` | Work the moving walls have done **on the atoms** so far: minus the force on the wall times its displacement, summed over the steps. Positive for a compression |

With the thermostat off (`ensemble="nve"`) the change in `total_energy` equals
`work` to the accuracy of the integration, at any piston speed. The wall stays
where the piston left it: the next run starts from the smaller box, and a run
without `piston` keeps it at rest.

## What the thermostat does to a gas

The thermostat gives every atom a small random push each step and a friction
`friction` (default 1 in reduced units). At the default an atom forgets its
velocity within one time unit, which is right for holding a temperature but
turns the straight flights of a dilute gas into a random walk; the force on a
wall is then correct on average, but a mean over the same time scatters five
times more (see the validation below).
For a gas whose motion should look and behave like a gas, use a light
thermostat, `MDSimulation(..., friction=0.05)`: the temperature still settles
within about ten time units.

## Validation

`python scripts/validate_walls.py` prints the numbers below (a few minutes);
the fast checks are in `tests/test_walls.py`. All runs: 200 atoms, timestep
0.005, cutoff 2.5, eight independent starting states, native build on macOS.

**Energy.** Thermostat off, 1000 time units (200 000 steps): the total energy
changes by 0.05 (T\* = 2) and 0.14 (T\* = 4) between walls, root mean square over
the eight runs, which is 1 to 2 parts in 10 000 of the kinetic energy. The same
gas in a flat periodic box: 0.03 and 0.07.

**Force on a wall.** Thermostat on with friction 0.05, means over 2000 time
units. The dilute flat gas should have PV/NkT = 1 + B₂N/A, with the second
virial coefficient B₂ = 0.467 at T\* = 2 and 0.928 at T\* = 4 for this cut and
shifted potential.

| T\* | Box | Force on the right wall | N T / Lx | PV/NkT, right wall | PV/NkT, all four walls | Expected |
|---|---|---|---|---|---|---|
| 2 | 100 x 100 | 4.04 ± 0.06 | 4 | 1.008 | 1.015 | 1.009 |
| 4 | 100 x 100 | 8.10 ± 0.17 | 8 | 1.016 | 1.020 | 1.019 |
| 2 | 200 x 100 | 1.98 ± 0.07 | 2 | 0.986 | 1.005 | 1.005 |
| 2 | 100 x 200 | 3.97 ± 0.05 | 4 | 0.994 | 1.010 | 1.005 |

The ± is the standard deviation between the eight runs: one mean over 2000 time
units is good to 1.5 to 3.5 per cent, because the wall takes only about one hit
per time unit. The mean of eight runs has a third of that, and all four walls
together half of that again. The force doubles with the temperature (factor
2.01) and halves when the box is twice as long (0.49); in the box that is twice
as wide the force on the right wall is unchanged (0.98) and the pressure
halved. With the default friction 1 the first row reads 3.91 ± 0.28: the same
mean within its much larger scatter.

**Piston.** Thermostat off, the gas settled at T\* = 2, the right wall from
x = 100 to 50.

| Piston speed | Stroke | T after / T before | Work on the gas | Largest \|ΔE − W\| |
|---|---|---|---|---|
| 0.05 | 1000 time units | 2.03 ± 0.02 | 426 ± 19 | 0.13 |
| 0.1 | 500 | 2.05 ± 0.02 | 431 ± 23 | 0.09 |
| 0.2 | 250 | 2.06 ± 0.03 | 438 ± 17 | 0.37 |

The change in total energy equals the work to 3 parts in 10 000 at the slow
speed, and to better than 1 in 1000 at every speed. A flat ideal gas compressed
slowly to half its area doubles its temperature (T V = constant, γ = 2). The
Lennard-Jones gas at this density ends 1.5 to 2.5 per cent higher; the adiabat
computed with its second virial coefficient gives 2.044. The spread in the work
between runs is the spread in the starting temperature (2.06 ± 0.10 for 200
atoms). With atoms that do not collide (ε = 10⁻¹², σ = 10⁻⁹) the ratio is
2.46 ± 0.12: only the motion along x is heated, fourfold, which gives 2.5.

In the browser build (Pyodide) the trajectories differ from the native ones
after a short time, as rounding differs in the last digit and the motion is
chaotic; the numbers above agree within their scatter.
