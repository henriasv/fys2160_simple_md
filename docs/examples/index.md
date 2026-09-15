# Examples

The [executed notebook](notebook.md) shows a crystal, a gas, a mixture and bonded molecules,
with recorded plots and working particle replays. These are standalone examples
of the interface, not worked course exercises.

See also the [worked NPT temperature sweep](isobar.md): measurements of pressure,
volume and Z at one fixed pressure, with saved data and pyplot figures.

## Settings to change

Each complete example writes out the standard physical and recording settings.
They are the starting point for students to edit, even when they match the defaults.
All values use the same LJ reduced reference units.

| Setting | Meaning |
|---|---|
| `ensemble="nvt"`, `temperature=2` | Fixed volume, thermostat target T = 2 |
| `timestep=0.005` | Integration step for atomic examples; molecules use 0.002 |
| `cutoff=2.5` | Pair-interaction range in reference length units |
| `steps=10_000` | Number of integration steps; at timestep 0.005 this is time 50 |
| `thermo_every=100` | Record thermodynamic data every 100 MD steps |
| `trajectory_every=500` | Save particle coordinates every 500 MD steps for replay |
| `storage_name="gas"` | Saved run folder; using the name again replaces its contents |

`thermo_every` and `trajectory_every` control output, not the integration timestep.
Repeated `sim.run(...)` calls continue the same simulation and retain its physical
settings. Thermostat coupling (`friction=1`), random seed (`seed=87287`) and other
controls are described in the [API reference](../api/simulation.md).

## Pressure of a gas

```python
import matplotlib.pyplot as plt
from fys2160_md import FCC, Random, MDSimulation

system = Random(N=500, rho=.01)
sim = MDSimulation(
    system,
    pair_potential="lj",
    epsilon=1,
    sigma=1,
    ensemble="nvt",
    temperature=2,
    timestep=0.005,
    cutoff=2.5,
    storage_name="gas",
)
# Equilibrate, then record a measurement run under the same conditions.
sim.run(
    steps=10_000,
    thermo_every=100,
    trajectory_every=500,
)
result = sim.run(
    steps=10_000,
    thermo_every=100,
    trajectory_every=500,
)  # The same storage name replaces the equilibration output.
print(result.thermo.compressibility_factor.mean())
data = result.thermo
plt.figure(figsize=(8, 5))
plt.subplot(2, 1, 1)
plt.plot(data.time, data.temperature)
plt.ylabel("Temperature")
plt.subplot(2, 1, 2)
plt.plot(data.time, data.pressure)
plt.ylabel("Pressure")
plt.xlabel("Time")
plt.tight_layout()
plt.show()
```

The ideal-gas value of Z is approximately one. Compare densities and allow enough
time to equilibrate before interpreting a mean.

## Add energy without a thermostat

```python
system = FCC(N=256, rho=.1)
sim = MDSimulation(
    system,
    pair_potential="lj",
    epsilon=1,
    sigma=1,
    ensemble="nve",
    temperature=1,
    timestep=0.005,
    cutoff=2.5,
)
result = sim.run(
    steps=1000,
    heat_rate=1,
    thermo_every=100,
    trajectory_every=500,
    storage_name="warming",
)
data = result.thermo
plt.plot(data.time, data.temperature)
plt.xlabel("Time")
plt.ylabel("Temperature")
plt.show()
```

This selects fixed-volume dynamics without a thermostat, then adds energy at
the chosen rate. It demonstrates the controls; no response coefficient is fitted.

## Flexible and rigid molecules

```python
from fys2160_md import Random, MDSimulation, HarmonicBond, RigidBond

system = Random(N=100, rho=.01, molecule="diatomic",
                bond=HarmonicBond(length=.7, stiffness=100))
sim = MDSimulation(
    system,
    pair_potential="lj",
    epsilon=1,
    sigma=1,
    exclude_bonded_pairs=True,
    ensemble="nvt",
    temperature=2,
    timestep=0.002,
    cutoff=2.5,
)
result = sim.run(
    steps=5000,
    thermo_every=100,
    trajectory_every=500,
    storage_name="molecules",
)
result.view(projection="perspective")
```

N counts **molecules**: this system has 200 atoms. Each connected pair feels the
spring; atoms in different molecules feel LJ. Both atoms have the same type.
The exclusion applies to connected pairs, not all pairs of that atom type.

Replace the spring with `RigidBond(length=.7)` to fix each bond length while
allowing translation and rotation. The default molecular timestep is 0.002;
a stiffer spring may need a smaller timestep. See [Bonds and interactions](../api/interactions.md).

## Pressure control

```python
system = FCC(N=500, rho=.005)
sim = MDSimulation(
    system,
    pair_potential="lj",
    epsilon=1,
    sigma=1,
    ensemble="npt",
    temperature=2,
    timestep=0.005,
    cutoff=2.5,
    pressure=.01,
    storage_name="pressure-control",
)
result = sim.run(
    steps=10_000,
    thermo_every=100,
    trajectory_every=500,
)
data = result.thermo
plt.figure(figsize=(8, 5))
plt.subplot(2, 1, 1)
plt.plot(data.time, data.volume)
plt.ylabel("Volume")
plt.subplot(2, 1, 2)
plt.plot(data.time, data.pressure)
plt.ylabel("Pressure")
plt.xlabel("Time")
plt.tight_layout()
plt.show()
```

Pressure control changes the volume while the thermostat maintains the target
temperature. Compare the measured pressure trace with the requested target.

## Analyse later

```python
from fys2160_md import Run

saved = Run.load("runs/gas")
print(saved.thermo.temperature.mean())
saved.view()
```

No simulation runs when these files are loaded. Reusing a storage name intentionally
replaces its files; choose a new name when you want to retain both results.

## A mixture with different masses and interactions

All species share one set of LJ reduced reference units. Here B has four times
the mass, 1.2 times the size and half the well depth of A.

```python
from fys2160_md import Random, MDSimulation

system = Random(
    species={"A": 80, "B": 20},
    rho=.02,
    masses={"A": 1, "B": 4},
    min_distance=1.1,
)
sim = MDSimulation(
    system,
    pair_potential="lj",
    epsilon={"A": 1, "B": .5},
    sigma={"A": 1, "B": 1.2},
    mixing_rule="lorentz-berthelot",
    ensemble="nvt",
    temperature=2,
    timestep=0.005,
    cutoff=2.5,
)
mixture = sim.run(
    steps=3000,
    thermo_every=100,
    trajectory_every=30,
    storage_name="mixture",
)
mixture.view(max_frames=101, projection="perspective")
```

N is inferred from the counts (100 atoms). A scalar mass, epsilon or sigma applies
to every species; a dictionary must contain every species name exactly once.
For a mass-only comparison, retain `epsilon=1, sigma=1` and vary `masses`.
This changes the motion without changing the potential energy function.

Unlike interactions use sigma_AB = (sigma_A + sigma_B)/2 and
 epsilon_AB = sqrt(epsilon_A epsilon_B). The mixing rule is a stated model choice,
not a universal material law. The common cutoff is in reference length units;
it is not automatically multiplied by each species' sigma.

```python
import numpy as np

atoms = sim.system.atoms
for name in sim.system.species:
    selected = atoms.species == name
    mean_v2 = np.mean(np.sum(atoms.velocities[selected]**2, axis=1))
    print(name, selected.sum(), mean_v2)
```

This is a snapshot, not a time average. At equilibrium, mean squared speed scales
approximately as 1/mass at fixed temperature (with a finite-system correction
from removing total momentum). The [recorded notebook](notebook.md) includes
this mixture with working playback.

For a molecular mixture, counts refer to molecules and masses to **each atom**:

```python
system = Random(species={"A": 24, "B": 8}, rho=.01,
                masses={"A": 1, "B": 4}, molecule="diatomic",
                bond=RigidBond(length=.7))
```

This creates 24 A–A and 8 B–B molecules, with masses 2 and 8 per molecule.
This version shares one bond; the next example assigns a separate bond to each species. Use explicit System arrays for heteronuclear
molecules; consecutive atom pairs are bonded.

## An air-like mixture with different bonds

Species can have separate bond lengths as well as masses and LJ parameters.
Here 80 N2-like and 20 O2-like molecules contain 200 atoms. The composition and
interaction/bond values are illustrative teaching choices, not a calibrated
model of air. Names are labels; they do not load physical parameters.

Both species share LJ reference units. Masses are per atom: a nitrogen atom is
the mass reference, and an oxygen atom has approximately 16/14 times that mass.
Each molecule's bond remains at its own specified length. The LJ exclusion
applies to the two connected atoms, not to other molecules of the same species.

```python
from fys2160_md import Random, MDSimulation, RigidBond

# Illustrative reduced parameters, not a fitted force field for real air.
air = Random(
    species={"N2": 80, "O2": 20},
    masses={"N2": 1, "O2": 16/14},  # per atom, relative to nitrogen
    molecule="diatomic",
    bond={
        "N2": RigidBond(length=0.65),
        "O2": RigidBond(length=0.75),
    },
    rho=0.01,
    min_distance=1.1,
)
sim = MDSimulation(
    air,
    pair_potential="lj",
    epsilon={"N2": 1, "O2": 0.9},
    sigma={"N2": 1, "O2": 1.05},
    mixing_rule="lorentz-berthelot",
    exclude_bonded_pairs=True,
    ensemble="nvt",
    temperature=2,
    timestep=0.002,
    cutoff=2.5,
)
result = sim.run(
    steps=3000,
    thermo_every=100,
    trajectory_every=30,
    storage_name="air-like",
)
result.view(max_frames=101, projection="perspective")
```

To model vibrations, replace the bond mapping when constructing the system:

```python
from fys2160_md import HarmonicBond

bonds = {
    "N2": HarmonicBond(length=0.65, stiffness=120),
    "O2": HarmonicBond(length=0.75, stiffness=80),
}
```

Use `bond=bonds` in Random or FCC. These stiffnesses are illustrative too.
Each species can have its own spring or Class2Bond coefficients, but a simulation
must use either all rigid or all flexible bonds. A stiffer spring can require a
shorter timestep; check energy conservation in NVE before interpreting results.

For a quantitative model, choose a published force field suitable for the
properties of interest and convert all its parameters into the **same** LJ
reference units. This classical bonded-LJ solver does not automatically supply
real nitrogen/oxygen material parameters. It also does not include trace air components.
