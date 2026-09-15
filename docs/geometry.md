# Build a system

Geometry belongs to a `System`. `MDSimulation(system, ...)` copies it and manages
the experiment; `sim.run(...)` continues that experiment.
All boxes are periodic in three dimensions. A scalar `box` means a cube side
length; a triple specifies rectangular side lengths.

## Density or box size

```python
from fys2160_md import FCC, Random, MDSimulation, RigidBond

crystal = FCC(N=500, rho=.8)
crystal = FCC(N=500, box=12)
rectangular = FCC(N=96, box=(6, 9, 12))
gas = Random(N=500, rho=.03, min_distance=.9)

sim = MDSimulation(
    gas,
    pair_potential="lj",
    epsilon=1,
    sigma=1,
    ensemble="nvt",
    temperature=2,
    timestep=0.005,
    cutoff=2.5,
)
result = sim.run(
    steps=5000,
    sample_every=100,
    save_every=500,
    storage_name="random-gas",
)
```

Specify **N and rho**, or **N and box**; specifying both rho and box is an error.
`system.N`, `system.rho` and `system.box.lengths` report the actual geometry.

## FCC: a complete periodic crystal

A cubic FCC box contains N = 4n³ atoms: 32, 108, 256, 500, 864, 1372, ….
N=1024 is rejected with nearby valid counts. The constructor never rounds N or
randomly removes lattice sites. Rectangular boxes must contain whole conventional
cells with the same lattice constant in all directions.

For diatomics, N counts **molecules**, rho is their number density, and FCC
sites locate the molecular centres. Partners start at the specified bond length with random orientations.

```python
molecules = FCC(N=500, rho=.01, molecule="diatomic", bond=RigidBond(length=.7))
```

## Random: no close starting pairs

```python
gas = Random(N=108, box=(10, 12, 14), min_distance=1.0, seed=12)
```

`min_distance` applies to atom separations, including across periodic edges.
The two bonded partners within a molecule are exempt. Placement uses rejection
sampling with at most 1000 trials per particle by default. An impossible or
very dense packing raises an error; it never silently weakens the separation.
Use FCC for a dense crystal. Reusing the seed reproduces the initial state.

## Explicit arrays

```python
from fys2160_md import System

pair = System([[4, 5, 5], [5.15, 5, 5]], box=10,
              velocities=[[.1, 0, 0], [-.1, 0, 0]])
sim = MDSimulation(
    pair,
    pair_potential="lj",
    epsilon=1,
    sigma=1,
    ensemble="nve",
    temperature=2,
    timestep=0.005,
    cutoff=2.5,
)
```

Inputs are copied. Explicit velocities must have zero total momentum. Diatomic
partners must be consecutive; rigid bonds must satisfy the length and tangency
constraints. The solver requires every box length to exceed twice its cutoff.

See the [System reference](api/system.md) for all defaults.

## Mixtures

```python
system = FCC(species={"A": 80, "B": 28}, rho=.2,
             masses={"A": 1, "B": 4})
```

The total (108 here) must be an exact FCC count. Individual species counts need
not be. Species are randomly assigned to lattice sites with a reproducible seed.
Random accepts the same species/masses arguments and checks all interparticle
separations regardless of species. Choose min_distance with the largest sigma
in mind; geometry is constructed before the interactions are chosen.
