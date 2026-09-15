# MDSimulation

```python
from fys2160_md import FCC, MDSimulation

system = FCC(N=500, rho=.1)
sim = MDSimulation(system, temperature=2, ensemble="nvt")
result = sim.run(steps=10_000, storage_name="gas")
```

MDSimulation owns a **copy** of the supplied system. It initializes missing
velocities, controls how the system evolves and records results. Repeated
`sim.run(...)` calls continue its clock and state. The starting system remains
unchanged, so it can be reused for independent experiments.

## Constructor

```python
MDSimulation(system, *, pair_potential="lj", epsilon=1., sigma=1.,
           exclude_bonded_pairs=True, temperature=2., ensemble="nvt", timestep=None,
           cutoff=2.5, skin=.4, friction=1., pressure=.01, pressure_time=5.,
           heat_rate=0., seed=87287, output_dir="runs", storage_name="simulation")
```

| Setting | Meaning |
|---|---|
| `system` | A System returned by FCC, Random or explicit arrays; copied on construction |
| `pair_potential` | `"lj"` (12–6, default) or `"lj96"` (9–6); either works with atoms or molecules |
| `epsilon`, `sigma` | Positive LJ energy and length parameters; both default to 1 |
| `exclude_bonded_pairs` | True by default: connected atom pairs feel the bond only |
| `temperature` | Initial temperature for missing velocities and thermostat target |
| `ensemble` | `nve`, `nvt`, or atomic-only `nph` / `npt` |
| `timestep` | Defaults to 0.005 for atoms, 0.002 for molecules |
| `cutoff`, `skin` | Interaction cutoff and Verlet neighbor-list skin, in the fixed reference length unit |
| `friction` | Langevin friction, inverse reduced time |
| `pressure`, `pressure_time` | Target pressure and barostat time scale |
| `heat_rate` | Total energy added per reduced time; default zero |
| `seed` | Initial velocity / thermostat seed, integer from 1 through 2⁶⁴−1 |
| `output_dir` | Parent folder for saved runs |
| `storage_name` | Default run folder name; reuse replaces its previous contents |

Explicit velocities in the System are preserved. A geometry seed controls
placement and orientation; an MDSimulation seed controls velocities and stochastic
integration. Every box length must exceed twice the simulation's cutoff.

Interaction parameters are fixed when constructing the simulation and saved with
its results. They do not change the reference units. For example, `sigma=1.2`
with `cutoff=3.0` gives a cutoff of 2.5 interaction diameters. Bond lengths,
box dimensions, timestep and temperature are not rescaled automatically.
See [Bonds and interactions](interactions.md) for equations and examples.

`Simulation` remains an alias for `MDSimulation` for existing import statements.
The earlier static `Simulation.run(system, ...)` pattern is replaced by an instance.

## Run

```python
sim.run(timestep=None, steps=1000, *,
        ensemble=None, temperature=None, heat_rate=None, pressure=None,
        sample_every=100, save_every=500, show=False, max_fps=None,
        frame_every=100, storage_name=None)
```

| Argument | Meaning |
|---|---|
| `steps` | Nonnegative integer number of additional MD steps |
| `timestep` | Positive reduced timestep; None retains the current value |
| `ensemble` | `nve`, `nvt`, `nph`, `npt`; None retains the current value |
| `temperature` | New thermostat target; does not instantly reset velocities |
| `pressure` | Pressure-control target, atomic NPH/NPT only |
| `heat_rate` | Total energy per reduced time; positive heats, negative cools |
| `sample_every` | Record thermodynamic quantities at this interval, in steps |
| `save_every` | Save particle frames at this interval; None disables frames |
| `show` | Display one live notebook widget while running |
| `max_fps` | Optional positive frame-rate cap; paces MD when show=True |
| `frame_every` | MD steps between live snapshots, independent of disk sampling |
| `storage_name` | Override this call's name; None uses the simulation default |

Omitted physical settings persist. In particular, **set `heat_rate=0` to stop
heating**. Adding heat with an active thermostat is rejected. With heat input,
`nve` means fixed-volume dynamics without a thermostat, not constant energy.

## Current state and restart

| Attribute | Result |
|---|---|
| `sim.pair_potential`, `.epsilon`, `.sigma`, `.exclude_bonded_pairs` | Read-only interaction settings |
| `sim.system` | The owned, evolving physical System |
| `sim.step`, `.time` | Cumulative simulation step and reduced time |
| `sim.settings` | Copy of current ensemble, timestep, targets and heat rate |
| `sim.observables` | Current measured temperature, pressure, energy and other quantities |
| `sim.forces` | Read-only force snapshot |
| `sim.N`, `.rho`, `.atoms`, `.box` | Convenient access to the owned system's properties |

```python
sim = MDSimulation.from_run("runs/gas")
continued = sim.run(steps=1000, storage_name="continued")
```

Restart restores the checkpoint, clock, random state and barostat state. It
accepts a Run object or a directory, and optional `output_dir=`. Completed and
interrupted runs can be restarted; failed runs cannot. This also reads the
previous package's checkpoint format.

## Live view

```python
sim.run(steps=10_000, show=True, max_fps=20, frame_every=100)
```

The cap slows wall time without changing the trajectory, physical timestep or
file sampling. Without a cap, the solver runs at full speed with display-only
updates limited to 20 fps. Slow calculations may run below a requested cap.
Only the latest snapshot is sent; at most 1500 atoms are shown.

- One-finger click-drag rotates; two-finger click-drag / Shift-drag pans.
- Two-finger scrolling zooms; buttons and a slider are also available.
- Perspective and orthographic views use sphere-surface depth for correct overlaps.
- Particle radius is 0.5σ; Reset restores pan, rotation and zoom.

## Storage and interruption

A call writes `output_dir/storage_name`. Reusing a name replaces that saved run,
including old frames. Use distinct names to keep comparisons. Unrelated folders
and symlinks are protected from replacement.

Keyboard interruption saves the completed state with status `interrupted`.
A numerical failure is marked `failed` and blocks further use of that simulation.
Every completed/interrupted run has initial/final samples and a final checkpoint.
