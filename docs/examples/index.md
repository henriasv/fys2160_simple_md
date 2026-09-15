# Examples

The [executed notebook](notebook.md) shows a vibrating crystal and a dilute gas,
with recorded plots and working particle replays. These are standalone examples
of the interface, not worked course exercises.

## Pressure of a gas

```python
import matplotlib.pyplot as plt
from fys2160_md import FCC, Random, MDSimulation

system = Random(N=500, rho=.01)
sim = MDSimulation(system, temperature=2, storage_name="gas")
sim.run(steps=10_000)           # Equilibrate.
result = sim.run(steps=10_000)  # Replace with the measurement run.
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
sim = MDSimulation(system, temperature=1)
result = sim.run(steps=1000, ensemble="nve", heat_rate=1,
                        storage_name="warming")
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
sim = MDSimulation(system, pair_potential="lj", epsilon=1, sigma=1,
                   exclude_bonded_pairs=True, temperature=2)
result = sim.run(steps=5000, storage_name="molecules")
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
sim = MDSimulation(system, temperature=2, ensemble="npt", pressure=.01,
                 storage_name="pressure-control")
result = sim.run(steps=10_000)
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
