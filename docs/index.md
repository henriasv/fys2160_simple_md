<div class="hero" markdown>
<span class="eyebrow">FYS2160 · statistical physics in Python</span>
# A small engine. Physics you can see.

Build a system, let it evolve, and connect particle motion to temperature,
pressure and heat capacity. A native C solver does the numerical work;
a few Python lines define the experiment.

[Install & run](installation.md){ .md-button .md-button--primary }
[Explore a recorded notebook](examples/notebook.md){ .md-button }
</div>

## Three steps to an experiment

```python
from fys2160_md import FCC, MDSimulation

system = FCC(N=500, rho=0.1)
sim = MDSimulation(
    system,
    pair_potential="lj",
    epsilon=1,
    sigma=1,
    ensemble="nvt",
    temperature=2,
    timestep=0.005,
    cutoff=2.5,
)
result = sim.run(
    steps=10_000,
    sample_every=100,
    save_every=500,
    storage_name="gas",
)
result.view()
```

`System` holds the box and atoms. `MDSimulation(system, ...)` owns and evolves a copy. `Run` reads,
plots and replays the recorded result—even in a later Python session.

<div class="grid cards" markdown>

- **Choose the geometry**

    A periodic FCC crystal or self-avoiding random positions. Set density or box size.

    [Build a system →](geometry.md)

- **Ask a physics question**

    Measure pressure, heat a gas, compare molecular degrees of freedom.

    [Worked examples →](examples/index.md)

- **Watch a saved run**

    Play, rotate, pan and zoom recorded particle trajectories. No Python kernel needed.

    [Executed notebook →](examples/notebook.md)

- **Find an argument**

    Defaults, units, saved quantities and the full public interface in one place.

    [API reference →](api/system.md)

</div>

## Designed for the course

- Atomic Lennard–Jones and flexible or rigid diatomic models.
- Periodic 3D dynamics, thermostat and atomic pressure control.
- Explicit heating for simple heat-capacity measurements.
- Live notebook views and portable saved trajectories.

All numbers use LJ reduced reference units: σ_ref = ε_ref = mass_ref = kB = 1.
Species may have different masses and interaction parameters in those shared units.
The examples use ordinary NumPy, pandas and pyplot. Read the
[physics notes](physics.md) for the model conventions and limits.
