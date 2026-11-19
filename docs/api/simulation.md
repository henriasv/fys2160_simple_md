# MDSimulation

```python
from fys2160_md import FCC, MDSimulation

start = FCC(N=500, rho=.1)
sim = MDSimulation(temperature=2, ensemble="nvt", timestep=0.005, cutoff=2.5)
first = sim.run(start, steps=10_000, storage_name="gas")
second = sim.run(first.state, steps=10_000, storage_name="gas-continued")
```

An MDSimulation holds the **interactions and default run controls**, never
particles. Every `sim.run(state, ...)` takes a state (a [System](system.md)),
evolves a copy of it and saves the result. The final state is `run.state`: pass
it to the next run to continue. The input state is never changed, so one state
can start several independent runs, for example at two temperatures.

## Constructor

```python
MDSimulation(*, pair_potential="lj", epsilon=1., sigma=1.,
             mixing_rule="lorentz-berthelot", exclude_bonded_pairs=True,
             temperature=None, ensemble="nvt", timestep=None, cutoff=None,
             G=1., softening=.05, skin=None, friction=None, pressure=None,
             pressure_time=None, heat_rate=0., seed=87287, output_dir="runs",
             storage_name="simulation", units="reduced", charge=None,
             coulomb=None, coulomb_cutoff=None, kspace_accuracy=1e-5, lj_shift=True)
```

| Setting | Meaning |
|---|---|
| `pair_potential` | `"lj"` (12–6, default), `"lj96"` (9–6), or isolated `"gravity"` |
| `epsilon`, `sigma` | Positive scalar or per-species dictionary; both default to 1 |
| `mixing_rule` | `"lorentz-berthelot"`: arithmetic sigma and geometric epsilon; the currently supported rule |
| `exclude_bonded_pairs` | True by default: connected atom pairs feel the bond only |
| `temperature` | Default thermostat target, and the temperature for a state without velocities (default 2 in reduced units) |
| `ensemble` | Default ensemble: `nve`, `nvt`, or `nph` / `npt` for atoms and rigid water |
| `timestep` | Default timestep; None gives 0.005 for atoms and 0.002 for molecules (reduced units) |
| `cutoff`, `skin` | LJ cutoff (2.5 when None) and neighbour-list skin (0.4). Gravity requires cutoff=None |
| `G`, `softening` | Gravity coupling (default 1) and smoothing length (default 0.05); softening=0 gives exact inverse-square gravity |
| `friction` | Langevin friction, inverse time |
| `pressure`, `pressure_time` | Default target pressure and barostat time scale |
| `heat_rate` | Default total energy added per time; zero |
| `seed` | Seed for random velocities and for the thermostat of a state that has not run yet, integer from 1 through 2⁶⁴−1 |
| `output_dir` | Parent folder for saved runs |
| `storage_name` | Default run folder name; reuse replaces its previous contents |
| `units` | `"reduced"` (default) or `"real"` (LAMMPS real units, for water); must match the states it runs |
| `charge` | Per-species charges in e; the total must be zero. Requires `coulomb="pppm"` |
| `coulomb` | `"pppm"`: long-range electrostatics by particle-particle particle-mesh Ewald summation |
| `coulomb_cutoff` | Real-space Ewald cutoff; defaults to `cutoff` |
| `kspace_accuracy` | Relative PPPM force accuracy, default 1e-5 |
| `lj_shift` | True subtracts the Lennard-Jones energy at the cutoff; False is the plain truncated potential (water, charges and real units only) |

In real units (water), `temperature` (K), `timestep` (fs) and `cutoff` (Å) must be
given, and `pressure` (atm) for npt/nph. See [Water and electrostatics](../water.md).

Interaction parameters are the same for every run of a simulation and are saved
with each result. They do not change the reference units: for example,
`sigma=1.2` with `cutoff=3.0` gives a cutoff of 2.5 interaction diameters. Bond
lengths, box dimensions, timestep and temperature are not rescaled automatically.
See [Bonds and interactions](interactions.md) for equations and examples.
Per-species dictionaries must match the species of the state that is run.

## Run

```python
run = sim.run(state, steps=1000, *, timestep=None, ensemble=None,
              temperature=None, pressure=None, heat_rate=None,
              coulomb_cutoff=None, kspace_accuracy=None, deform=None, seed=None,
              thermo_every=100, trajectory_every=500, show=False, max_fps=None,
              frame_every=100, storage_name=None)
```

| Argument | Meaning |
|---|---|
| `state` | The System to start from: a geometry such as `FCC(...)`, or `run.state` from an earlier run |
| `steps` | Nonnegative integer number of MD steps |
| `timestep`, `ensemble`, `temperature`, `pressure`, `heat_rate` | Controls for this run only; None uses the constructor's value |
| `temperature` | Thermostat target from the first step of this run (a step change, not a ramp); velocities are not reset, so the thermostat relaxes the system towards it |
| `coulomb_cutoff`, `kspace_accuracy` | Electrostatics settings of a charged system for this run |
| `deform` | nve/nvt only: scale every box length linearly to this factor by the end of the run, moving whole molecules |
| `seed` | Restart the thermostat's random sequence (and set random velocities) for this run |
| `thermo_every` | Record thermodynamic quantities at this interval, in steps |
| `trajectory_every` | Save particle frames at this interval; None disables frames |
| `show` | Display one live notebook widget while running |
| `max_fps` | Optional positive frame-rate cap; paces MD when show=True |
| `frame_every` | MD steps between live snapshots, independent of disk sampling |
| `storage_name` | Override this run's folder name; None uses the simulation default |

**Nothing carries over between runs except the state you pass in.** A control
given to one run applies to that run only; the next run uses the constructor's
value again unless you repeat it. Adding heat with an active thermostat is
rejected. With heat input, `nve` means dynamics without a thermostat, not
constant energy. Periodic LJ keeps a fixed box; an isolated gravitational cluster
has no prescribed physical volume.

A state without velocities gets random ones at `temperature`. Explicit velocities
are kept. A state from a run carries its clock (time, step), the thermostat's
random state, the barostat velocity, the PPPM mesh choice and the heat added so
far, so continuing from `run.state` gives exactly the same trajectory as one
longer run. Running the same state twice with the same controls gives the same
result; pass `seed` for a different thermostat sequence.

### Thermodynamic measurements versus trajectory frames

`thermo_every` records temperature, pressure and energy in `thermo.csv`.
`trajectory_every` saves particle snapshots for replay and later particle analysis.
Both intervals count integration steps; neither changes the physical timestep.
Initial and final outputs are also saved. Use `trajectory_every=None` to omit
particle frames while keeping measurements and the final state.
`frame_every` is separate: it controls live display updates when `show=True`.

The earlier names `sample_every` and `save_every` remain accepted as aliases.
Do not pass both names for the same option. New examples use the descriptive names.

## Saved runs and states

```python
from fys2160_md import Run

saved = Run.load("runs/gas")
continued = sim.run(saved.state, steps=1000, storage_name="continued")
```

`run.state` reads the final checkpoint, so it also works in a later Python
session. Completed and interrupted runs have a final state; failed runs do not.
`sim.settings` returns the default controls; `sim.pair_potential`, `.epsilon`,
`.sigma`, `.charge` and similar read-only properties return the interactions.

## Advance a run yourself

`sim.run` takes all its steps and saves the result. `sim.start` gives the same
run **in memory**, for you to advance a few steps at a time: in a loop of your
own, behind a live plot, or from a web page. It writes no files.

```python
from fys2160_md import FCC, MDSimulation

sim = MDSimulation(
    pair_potential="lj",
    epsilon=1,
    sigma=1,
    ensemble="nvt",
    temperature=2,
    timestep=0.005,
    cutoff=2.5,
)
live = sim.start(FCC(N=500, rho=0.1))
while live.step < 2000:
    live.advance(10)                       # 10 MD steps
    if live.step == 1000:
        live.set(temperature=1)            # thermostat target from here on
    if live.step % 100 == 0:
        print(live.step, live.observables["temperature"])

speeds = (live.velocities**2).sum(axis=1)**0.5
cooled = sim.run(
    live.state,
    steps=5000,
    thermo_every=100,
    trajectory_every=500,
    storage_name="cooled",
)
```

```python
live = sim.start(state, *, timestep=None, ensemble=None, temperature=None,
                 pressure=None, heat_rate=None, coulomb_cutoff=None,
                 kspace_accuracy=None, seed=None)
```

The arguments are those of `run` that are not about output: the state is
copied, omitted controls take the constructor's values, and a state without
velocities gets random ones. The result is a `LiveRun`:

| Member | Meaning |
|---|---|
| `live.advance(steps)` | Take this many MD steps (a nonnegative integer) |
| `live.observables` | Dictionary of the current thermodynamic quantities: the columns `run` records in `thermo.csv` |
| `live.positions`, `live.velocities` | Read-only (atoms, 3) snapshots; positions are wrapped into a periodic box |
| `live.box` | Read-only snapshot of the three box lengths, as in a saved frame |
| `live.step`, `live.time` | Clock, continuing from the state's |
| `live.set(timestep=, ensemble=, temperature=, pressure=, heat_rate=)` | Change controls for the steps that follow; omitted ones keep their values |
| `live.settings` | The controls now in force |
| `live.state` | The current state, an independent `System`: pass it to `sim.run` or `sim.start` |

**How the steps are divided between calls does not change the trajectory.**
Advancing 1000 steps one at a time, in tens, or in one call gives the same
positions and velocities bit for bit, and so does `sim.run`, which is this loop
with files. Continuing from `live.state` matches going on with `live`. Only the
`time` and `heat_added` values are summed per call and may differ in the last digit.

`live.set(...)` has the effect of a new run started from `live.state` with those
controls: a new temperature or pressure is the target from the next step (a
step change, not a ramp) and velocities are not reset. The controls are checked
together, so `live.set(ensemble="nve", heat_rate=2)` is one valid change. The
simulation's own defaults (`sim.settings`) are not affected. Interactions,
cutoff, friction and `deform` cannot be changed on a live run; start another
from `live.state`.

Arrays are snapshots: they do not change when the run advances, and editing the
run through them is not possible. To change velocities, take the state and use
[`state.with_velocities`](system.md#reverse-or-perturb-velocities).

A keyboard interrupt stops `advance` after a whole step and is raised again, so
your loop ends and the run remains usable. After a numerical failure the run
cannot continue; start again from an earlier state.

The solver keeps its neighbour list and forces between calls, so short calls are
cheap: one step per call runs within about 15% of the rate of 500 steps per call
for a dilute gas, and within a few percent for a liquid.

## Live view

```python
sim.run(state, steps=10_000, show=True, max_fps=20, frame_every=100)
```

The cap slows wall time without changing the trajectory, physical timestep or
file sampling. Without a cap, the solver runs at full speed with display-only
updates limited to 20 fps. Slow calculations may run below a requested cap.
Only the latest snapshot is sent; at most 3000 atoms (whole molecules) are shown.

- One-finger click-drag rotates; two-finger click-drag / Shift-drag pans.
- Two-finger scrolling zooms. The control bar shows the run status, reset view, the projection and a **?** button for details; saved views add play/pause (also the space key), a frame slider and a playback speed menu from 0.25× to 8× of 10 frames per second (also the < and > keys).
- Perspective and orthographic views use sphere-surface depth for correct overlaps.
- LJ particle radius is 0.5σ; gravity uses display markers of radius 0.035. Reset restores pan, rotation and zoom.

Live and saved views are drawn with [OMOVI](https://github.com/andeplane/omovi)
by Anders Hafreager, the renderer used in Atomify, loaded from jsDelivr when the
view opens. It adds ambient occlusion in perspective view and draws bonds as
cylinders; drag rotates, right-drag pans and the mouse wheel zooms. Without a
network connection, or with `run.view(renderer="basic")`, the built-in renderer
is used. Views that are scrolled out of sight pause rendering.

## Storage and interruption

A run writes `output_dir/storage_name`. Reusing a name replaces that saved run,
including old frames. Use distinct names to keep comparisons. Unrelated folders
and symlinks are protected from replacement.

Keyboard interruption saves the completed state with status `interrupted`.
A numerical failure is marked `failed` and saves no final state.
Every completed/interrupted run has initial/final samples and a final state.

For mixtures, epsilon/sigma dictionaries must match `state.species`. See the
[mixture example](../examples/index.md#a-mixture-with-different-masses-and-interactions).


## Gravity observables

An open cluster has no bulk pressure, density, compressibility factor or enthalpy.
Its `thermo` output contains kinetic temperature, kinetic/potential/total energy,
`heat_added`, the force `virial` W, `virial_ratio=2K/(-W)` and `half_mass_radius`
(the radius enclosing half the mass around the COM). `pressure` is NaN for
compatibility with generic readers; no gas equation of state is reported.

Temperature is `2K/(3N-3)` in the COM frame, with kB=1. It measures particle
motion, not the internal temperatures of the stars that the particles represent.
G and softening are read-only as `sim.G` and `sim.softening`; the open boundary
condition belongs to the state.
