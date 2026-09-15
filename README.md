# FYS2160 Simple MD

A small native C molecular-dynamics engine with a Python interface for learning
statistical physics. Create a system, run it, then analyse or replay saved results.

**[Documentation & API](https://henriasv.github.io/fys2160_simple_md/)** ·
**[Executed notebook](https://henriasv.github.io/fys2160_simple_md/examples/notebook/)**

```sh
python -m pip install "fys2160-simple-md[notebook] @ git+https://github.com/henriasv/fys2160_simple_md.git"
```

Requires Python 3.10+, Git and a C compiler. On macOS install Xcode command-line
tools; on Linux install the usual build tools and Python development headers.
The Python import is `fys2160_md`.

```python
from fys2160_md import FCC, Random, MDSimulation, Run

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
    storage_name="argon",
)
print(result.thermo.compressibility_factor.mean())
result.view()  # Interactive replay in a notebook.
```

`Random(N=500, box=20, min_distance=.9)` gives self-avoiding random positions.
Use `show=True, max_fps=20` on a run for live notebook visualization.
Saved trajectories remain replayable after restarting Python. Reusing a
`storage_name` overwrites that result. All quantities use reduced LJ units.

The engine supports atomic LJ, flexible and rigid diatomics, NVE/NVT and atomic
NPH/NPT, plus controlled heating. It has no LAMMPS or Atomify runtime dependency.
See the [model details](https://henriasv.github.io/fys2160_simple_md/physics/) for
assumptions, limits and validation. This is a teaching tool, not a general MD engine.

## Development

```sh
python -m pip install -e '.[notebook,docs]'
python -m unittest discover -s tests
node tests/test_viewer.cjs
python scripts/build_docs.py --execute
python -m mkdocs serve
```

GitHub Actions runs the tests, executes the short example notebook and publishes
the documentation to GitHub Pages on pushes to `main`. The published notebook
is static HTML: particle trajectories play entirely in the browser, without a kernel.

Developed for FYS2160. Public examples are independent demonstrations; course
exercise sheets, lab solutions and their development history remain in the
private teaching-materials repository.

## Diatomic molecules

```python
from fys2160_md import Random, MDSimulation, HarmonicBond

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
    sample_every=100,
    save_every=500,
    storage_name="molecules",
)
```

N counts molecules here (200 atoms). Connected partners feel the spring;
atoms in different molecules feel Lennard–Jones. Use `RigidBond(length=.7)`
for fixed bond lengths. Geometry and interaction parameters use fixed reduced
reference units; changing sigma does not rescale the box or bond.
See the [interaction reference](https://henriasv.github.io/fys2160_simple_md/api/interactions/).

## Mixtures in shared LJ units

```python
system = Random(species={"A": 80, "B": 20}, rho=.02,
                masses={"A": 1, "B": 4}, min_distance=1.1)
sim = MDSimulation(
    system,
    pair_potential="lj",
    epsilon={"A": 1, "B": .5},
    sigma={"A": 1, "B": 1.2},
    ensemble="nvt",
    temperature=2,
    timestep=0.005,
    cutoff=2.5,
)
result = sim.run(
    steps=3000,
    sample_every=100,
    save_every=500,
    storage_name="mixture",
)
```

All species share fixed LJ reference units. Unlike pairs use the explicitly
named Lorentz–Berthelot mixing rule: arithmetic mean of sigma and geometric mean
of epsilon. Counts infer N, masses are per atom, and a scalar parameter applies
to every species. The examples page and executed notebook demonstrate mixtures.
