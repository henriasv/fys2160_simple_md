# System & geometry

```python
from fys2160_md import System, FCC, Random
```

A `System` describes the physical particles and box. It has no timestep,
thermostat or output folder. `sim.run(system, ...)` evolves a copy; the state at
the end is `run.state`, a System that also carries the clock (`time`, `step`).

## FCC

```python
FCC(*, N=None, rho=None, box=None, molecule=None, bond=None, model=None, species=None, masses=1., seed=87287,
    boundary="periodic")
```

Returns a `System` with a complete periodic FCC lattice and velocities not yet
assigned. Default: 500 particles in a cubic box at rho = 0.001. Supply either rho or box; N counts
atoms or molecules according to molecule. Cubes require N=4n³; rectangular boxes
require commensurate whole conventional cells. `seed` controls molecular orientations.
`boundary="walls"` keeps the sites and replaces the periodic images by [walls](../walls.md).

## Random

```python
Random(*, N=None, rho=None, box=None, min_distance=.9,
       max_attempts=1000, molecule=None, bond=None, model=None, species=None, masses=1., seed=87287,
       dimension=None, boundary="periodic")
```

Returns a `System` with self-avoiding random positions and velocities not yet
assigned. Minimum-image distances between atoms of distinct particles are at
least `min_distance`. Partners in a molecule are exempt. `max_attempts` limits
trials **per particle**; unable-to-pack requests raise `ValueError`. The seed
controls placement and molecular orientations, independently of the simulation seed.

A box with two lengths, `box=(Lx, Ly)`, or `dimension=2` gives a flat,
two-dimensional system of atoms. `boundary="walls"` puts a wall on every face
of the box instead of periodic images, or one word per axis chooses, such as
`("walls", "periodic")`; atoms then start at least half of `min_distance` from a
wall. See [Flat boxes, walls and pistons](../walls.md).

## Explicit System

```python
System(positions, box, *, velocities=None, masses=None, molecule=None, bond=None, model=None, species=None, boundary="periodic")
```

Positions are a (number of atoms, 3) array; box is a scalar or length-3 array.
Positions with two columns and a box with two lengths make a flat system.
Optional velocities have the same shape. Masses are a positive per-atom array,
defaulting to one. Inputs are copied. With `boundary="periodic"` (the default), positions wrap into the box. With `boundary="open"`, coordinates are preserved and the box is only a viewing frame: no walls or periodic images. Open systems currently support unbonded particles and gravity. With `boundary="walls"`, or one of `"periodic"` and `"walls"` per axis, the faces across those axes are [walls](../walls.md): every atom must lie between them, and the total momentum along those axes need not be zero.
Use `molecule="diatomic"` with a `HarmonicBond` or `RigidBond`.
The default is unbonded atoms; diatomics default to `HarmonicBond()` if bond is omitted.

For diatomics, partners are consecutive. Different atom masses are supported. Rigid
bonds must have their specified length; supplied velocities must be tangent to those bonds.
The solver requires zero total momentum when explicit velocities are supplied.

```python
system = System([[4, 5, 5], [5.15, 5, 5]], box=10,
                velocities=[[.1, 0, 0], [-.1, 0, 0]])
```

If velocities are `None`, MDSimulation initializes them at its requested temperature.
If supplied, they are **preserved**; the simulation's temperature argument then
sets the thermostat target, not an immediate velocity reset.

### Read the state

| Attribute | Result |
|---|---|
| `system.species` | Tuple of distinct species names |
| `system.atoms.species` | Read-only per-atom species labels |
| `system.molecule` | None for atoms, `"diatomic"` or `"water"` for molecules |
| `system.bond` | Immutable bond settings, or None |
| `system.bonds` | Read-only (number of bonds, 2) array of connected atom indices |
| `system.model` | Label: `atomic`, `diatomic-flexible`, `diatomic-rigid` or `water-rigid` |
| `system.units` | `"reduced"` (default) or `"real"` (LAMMPS real units; water) |
| `system.N` | Number of atoms, or molecules for diatomics and water |
| `system.boundary` | `"periodic"`, `"walls"`, a tuple with one of the two per axis, or `"open"` |
| `system.dimension` | 3, or 2 for a flat system: the number of columns of positions and velocities |
| `system.walls` | Names of the walls, such as `("left", "right", "bottom", "top")`; empty for a periodic box |
| `system.rho` | Current N / volume (N / area for a flat system); NaN for an unconfined cluster |
| `system.box`, `.lengths`, `.volume` | Box lengths and volume (two lengths and the area for a flat system). The box can be changed in place, `system.box[0] *= 2` ([below](#change-the-box)) |
| `system.atoms.positions`, `.masses` | Read-only array snapshots |
| `system.atoms.velocities` | Read-only snapshot, or None before velocities are assigned |
| `system.atoms.molecule_ids` | IDs identifying molecular partners |
| `system.dof` | Kinetic degrees of freedom: one per atom and direction, less one for each periodic direction (total momentum zero) and for each rigid constraint |
| `system.time`, `.step` | Clock of the state: zero for a new geometry, where the run ended for `run.state` |
| `system.copy()` | Independent copy of the physical state |
| `system.with_velocities(velocities)` | Copy with other velocities and the same positions, box and clock ([below](#reverse-or-perturb-velocities)) |
| `system.with_box(box)` | Copy in a larger periodic box, contents centred and molecules whole; enlarging one side gives a slab with vacuum ([evaporation](../water.md#evaporation-a-slab-with-vacuum)). Between walls the atoms stay where they are and the far wall moves |

Array snapshots do not let you edit the running simulation. Construct a new
System to specify different coordinates; use `with_velocities` for different
velocities. The explicit-array convenience
`System.from_arrays(positions, velocities, box, ...)` is also available.

## Change the box

The box is the one thing in a state that can be changed directly:

```python
state = run.state
state.box[0] *= 2          # twice as long along x
state.box[:] = 20, 20, 20  # or set every length
next_run = sim.run(state, steps=1000)
```

The atoms keep their positions; the next run starts from the new box. Across a
periodic axis, atoms that are now outside the box are wrapped back in. A smaller
box can therefore put atoms almost on top of each other, which makes the run blow
up: shrink a dense system a little at a time, or let a barostat (`ensemble="npt"`)
do it. Between walls, the wall at 0 stays and the far wall moves to the new
length; every atom must still be inside, so squeeze atoms with a
[piston](../walls.md) instead. `run.state` is the same state each time you read
it, so `run.state.box[0] *= 2` is kept for `sim.run(run.state)`. `with_box`
instead centres the contents in a larger periodic box.

### Reverse or perturb velocities

```python
state = run.state
backwards = state.with_velocities(-state.atoms.velocities)

velocities = state.atoms.velocities.copy()
velocities[0] += [1e-6, 0, 0]              # nudge one atom
velocities -= velocities.mean(axis=0)      # keep the total momentum zero (equal masses)
nudged = state.with_velocities(velocities)
```

`with_velocities` returns a copy; the original state is unchanged. Positions,
box, masses, the clock and the thermostat's random sequence are kept, so two
runs from `state` and `nudged` differ only through the nudge. Velocities need
one finite row per atom. A run requires zero total momentum (subtract the
mass-weighted mean velocity when masses differ), and velocities along a rigid
bond are rejected.

With `ensemble="nve"`, running `backwards` for as many steps as led to `state`
retraces the motion to its starting positions, up to rounding errors that grow
with the length of the run.

Clock, forces, measured observables, experiment settings and checkpoint restart
belong to [MDSimulation](simulation.md).

The legacy `model=` argument is retained for old code. `diatomic-flexible` selects
the old `Class2Bond`, and `diatomic-rigid` selects `RigidBond()`. Do not combine
`model=` and `molecule=`. New code should name the bond explicitly.

## Species and masses

For FCC/Random, `species="A"` names a pure system. A mapping such as
`species={"A": 80, "B": 20}` specifies exact counts; omitted N is their sum.
If N is supplied too, it must match. With neither N nor species counts, N is 500.
Species are shuffled reproducibly among sites using the geometry seed.
`masses=1` applies one mass to every atom; `masses={"A": 1, "B": 4}` sets per-species
masses. Every named species must have a positive mass. Molecular counts are
numbers of homonuclear molecules, while masses remain **per atom**.

For an explicit System, `species` is a single string or one label per atom;
`masses` is a per-atom array. This also supports heteronuclear diatomics:
consecutive partners may have different species and masses.

## A separate bond for each molecular species

```python
bond = {
    "N2": RigidBond(length=0.65),
    "O2": RigidBond(length=0.75),
}
```

Pass this mapping to `Random`, `FCC` or `System` with `molecule="diatomic"`.
Keys must match every molecular species exactly. These lengths are examples,
not physical parameters for air. Each homonuclear molecule uses the bond for its
species, including during integration and after checkpoint restart.

For flexible molecules, use separate `HarmonicBond(length=..., stiffness=...)`
or `Class2Bond(...)` entries. All entries must be rigid or all must be flexible;
mixing constrained and flexible molecules in one simulation is not supported.
Per-species mappings require homonuclear partners. Explicit heteronuclear pairs
continue to support a shared bond description. `system.bond` returns a copy of
the mapping; editing it does not change the system or a running simulation.

See the [air-like example](../examples/index.md#an-air-like-mixture-with-different-bonds).


## Water: IceIh and Random

```python
from fys2160_md import IceIh, Random, RigidWater

spce = RigidWater(oh_length=1.0, hoh_angle=109.47)          # Å, degrees
ice = IceIh(N=768, density=0.92, molecule="water", bond=spce, seed=1)
liquid = Random(N=512, box=24.86, molecule="water", bond=spce, min_distance=1.8)
```

Both return systems in LAMMPS `real` units with species `"O"` and `"H"`, three
atoms per molecule (O first). `IceIh` needs N = 8·nx·ny·nz in a nearly cubic box
and rejects other counts with suggestions. `RigidWater(..., m_distance=...)`
adds a TIP4P-type M site. See [Water and electrostatics](../water.md).

## Plummer: an isolated gravitational cluster

```python
Plummer(*, N=256, scale_radius=1., mass=1., G=1., softening=.05, seed=87287)
```

Returns an open `System` with equal particle masses and **supplied velocities**.
The density is proportional to `(1 + r²/scale_radius²)^(-5/2)`. Isotropic speeds
are sampled from the Plummer distribution, COM motion is removed, and their
scale is adjusted so `2K = -W` for this finite sample and its chosen softened
force. This prepares an approximately stationary cluster; let it settle before
measuring. The initial virial ratio by itself is not proof of equilibrium.

Use the **same G and softening** in MDSimulation. `mass` is the mass of each
particle, so the total mass is `N*mass`. Supports 2–4096 particles. `scale_radius`
is not the gravitational smoothing length. Its viewing frame has side
`12*scale_radius`, centered on the COM. The spatial distribution has an
untruncated tail, and particles outside the view remain in the calculation.

See the [worked gravitational example](../examples/gravity.md).
