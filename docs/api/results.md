# Saved Run

```python
from fys2160_md import Run, System, Simulation
saved = Run.load("runs/gas")
```

Loading a run reads files and does **not** start a simulation.

| Member | Meaning |
|---|---|
| `Run.load(path)` | Read one saved run folder |
| `Run.find(directory="runs")` | List recognized metadata folders |
| `run.path` | Resolved directory path |
| `run.metadata` | Configuration, units, completion status and output settings |
| `run.thermo` | Pandas DataFrame read from thermo.csv |
| `run.iter_frames()` | Stream saved frame dictionaries one at a time |
| `run.load_frames()` | Stack all frames in memory; raises if none were saved |
| `run.plot(*quantities, x="time")` | Matplotlib figure with one panel per quantity |
| `run.view(...)` | Self-contained HTML trajectory player for notebook output |

## Replay

```python
saved.view(max_frames=150, max_atoms=1500, projection="orthographic", zoom=1.0)
```

The default samples at most 150 saved frames and 1500 atoms; full recorded data
remain on disk. `projection` is `orthographic` or `perspective`; initial zoom
must be between 0.25 and 3. Camera controls and playback work in exported HTML.
This needs WebGL 2, but no Python kernel, remote scripts or simulation server.

## Recorded quantities

| Group | Columns |
|---|---|
| Clock | `step`, `time` |
| State | `temperature`, `pressure`, `volume`, `density`, `atom_density` |
| Energy | `kinetic_energy`, `potential_energy`, `total_energy` |
| Response | `compressibility_factor`, `enthalpy`, `heat_added` |
| Pressure control | `barostat_energy`, `extended_enthalpy` |

**Energies are totals.** Divide by N for energy per atom or molecule.
`density` uses that same particle count; `atom_density` always counts atoms.
For NPH/NPT, enthalpy uses target pressure; otherwise it uses instantaneous pressure.
See [physics conventions](../physics.md) for precise definitions.

## Files on disk

```text
runs/gas/
  metadata.json
  thermo.csv
  checkpoint.npz
  frames/0000000.npz
  frames/0000001.npz
  ...
```

Frame arrays include positions, unwrapped coordinates, velocities, box, step and
time. A final checkpoint is kept even when `save_every=None`.

```python
system = System.from_run(saved)
continued = Simulation.run(system, steps=1000, storage_name="continued")
```

A `Run` refers to a directory, not an immutable copy. Replacing that storage name
also changes what previously created Run objects read.
