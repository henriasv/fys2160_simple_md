# System & geometry

```python
from fys2160_md import System, FCC, Random
```

## FCC

```python
FCC(*, N=500, rho=None, box=None, model="atomic", seed=87287, **settings)
```

Returns a `System` with a complete periodic FCC lattice. Default: cubic box at
rho = 0.001. Supply either rho or box; N counts atoms or molecules according to model.
Cubes require N=4n³. Rectangular boxes require commensurate whole conventional cells.

## Random

```python
Random(*, N=500, rho=None, box=None, min_distance=.9,
       max_attempts=1000, model="atomic", seed=87287, **settings)
```

Returns a `System` with self-avoiding random placement. Minimum-image distances
between atoms of distinct particles must be at least `min_distance`. The bond
inside a diatomic molecule is exempt. `max_attempts` limits trials **per particle**;
unable-to-pack requests raise `ValueError`.

## System

```python
System(positions, box, *, velocities=None, masses=None, model="atomic",
       temperature=2., ensemble="nvt", timestep=None, cutoff=2.5,
       skin=.4, friction=1., pressure=.01, pressure_time=5.,
       heat_rate=0., seed=87287, output_dir="runs", storage_name="simulation")
```

Explicit positions are an (number of atoms, 3) array; box is a scalar or length-3
array. Velocities have the same shape. Masses are a positive per-atom array.
Inputs are copied and positions are wrapped into the periodic box. If velocities
are omitted, they are generated with zero centre-of-mass momentum and scaled to
the initial temperature. Diatomic models currently require equal atom masses.

`FCC` and `Random` accept these initial physical/storage settings as keyword arguments.

| Setting | Meaning |
|---|---|
| `model` | `atomic`, `diatomic-flexible`, or `diatomic-rigid` |
| `temperature` | Initial temperature and initial thermostat target, in reduced units |
| `ensemble` | `nve`, `nvt`, or atomic-only `nph` / `npt` |
| `timestep` | Defaults to 0.005 for atoms, 0.002 for molecules |
| `cutoff`, `skin` | Interaction cutoff and Verlet neighbor-list skin, in σ |
| `friction` | Langevin friction, inverse reduced time |
| `pressure`, `pressure_time` | Target pressure and barostat time scale |
| `heat_rate` | Total energy added per reduced time; default zero |
| `seed` | Integer from 1 through 2⁶⁴−1 |
| `output_dir` | Parent folder for saved runs |
| `storage_name` | Default run folder name; reuse replaces previous contents |

### Read the state

| Attribute | Result |
|---|---|
| `system.N` | Number of atoms, or molecules for diatomics |
| `system.rho` | Current N / volume |
| `system.box.lengths`, `.volume` | Box lengths and volume |
| `system.atoms.positions`, `.velocities`, `.masses` | Read-only array snapshots |
| `system.atoms.molecule_ids` | IDs identifying molecular partners |
| `system.forces` | Read-only force array snapshot |
| `system.step`, `.time` | Cumulative step and reduced time |
| `system.settings` | Copy of current ensemble, dt, targets and heat rate |
| `system.observables` | Current temperature, pressure, energy and other recorded quantities |
| `system.dof` | Active kinetic degrees of freedom, with COM/rigid constraints removed |

Array snapshots are not handles for editing the running system. Construct a new
`System` to specify different coordinates or velocities.

### Load a checkpoint

```python
system = System.from_run("runs/gas")
# Or: System.from_arrays(positions, velocities, box, masses=None, model="atomic", ...)
```

`from_run` restores the final checkpoint of a completed or interrupted run,
including velocities, box, time, random state and barostat state. Failed runs
cannot be resumed. `output_dir=` optionally changes the destination directory.
