# Simulation.run

```python
from fys2160_md import Simulation

result = Simulation.run(system, steps=10_000, storage_name="gas")
```

The runner advances the supplied `System` **in place** and returns a saved `Run`.
There is no separate simulation object to initialize. Different systems keep
independent state; successive calls on the same system continue from its current step.

```python
Simulation.run(system, timestep=None, steps=1000, *,
               ensemble=None, temperature=None, heat_rate=None, pressure=None,
               sample_every=100, save_every=500, show=False, max_fps=None,
               frame_every=100, storage_name=None)
```

| Argument | Meaning |
|---|---|
| `system` | A System returned by FCC, Random, explicit arrays or checkpoint loading |
| `steps` | Nonnegative integer number of additional MD steps |
| `timestep` | Positive reduced timestep; None retains the system's value |
| `ensemble` | `nve`, `nvt`, `nph`, `npt`; None retains the current value |
| `temperature` | Thermostat target; does not instantly reset velocities |
| `pressure` | Pressure-control target, atomic NPH/NPT only |
| `heat_rate` | Total energy per reduced time; positive heats, negative cools |
| `sample_every` | Record thermodynamic quantities every this many steps |
| `save_every` | Save particle frames at this interval; None disables frames |
| `show` | Display one live notebook widget while running |
| `max_fps` | Optional positive frame-rate cap; deliberately paces MD when show=True |
| `frame_every` | MD steps between live snapshots; independent of disk sampling |
| `storage_name` | Override this call's run name; None uses the system default |

Omitted physical settings persist. In particular, **set `heat_rate=0` to stop
heating**. Adding heat with an active thermostat is rejected. With heat input,
`nve` means fixed-volume dynamics without a thermostat, not a constant-energy ensemble.

## Live view

```python
Simulation.run(system, steps=10_000, show=True, max_fps=20, frame_every=100)
```

The cap slows wall time without changing the trajectory, physical timestep or
file sampling. Without a cap the solver runs at full speed, with display-only
updates limited to 20 fps. Slow calculations may run below a requested cap.
Only the latest snapshot is sent; at most 1500 atoms are shown.

- One-finger click-drag rotates; two-finger click-drag / Shift-drag pans.
- Two-finger scrolling zooms; buttons and a slider are also available.
- Choose perspective or orthographic projection; reset restores pan, rotation and zoom.
- Spheres have radius 0.5σ and use surface-depth testing for correct overlaps.

## Storage and interruption

A call writes `output_dir/storage_name`. Reusing a name replaces that saved run,
including old frames. Use distinct names to keep comparisons. The package refuses
to overwrite unrelated folders or symlinks.

Keyboard interruption saves the completed state and labels the run `interrupted`.
A numerical failure is marked `failed` and blocks further use of that in-memory system.
Every successful/interrupted run records its initial and final samples and a final checkpoint.
