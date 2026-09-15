# Bonds and interactions

A molecule specifies **which atoms are connected**. A bond describes their
internal motion. The pair potential describes interactions between other pairs.
No second atom type is needed for two identical atoms.

```python
from fys2160_md import Random, MDSimulation, HarmonicBond

system = Random(N=100, rho=.01, molecule="diatomic",
                bond=HarmonicBond(length=.7, stiffness=100))
sim = MDSimulation(system, pair_potential="lj", epsilon=1, sigma=1,
                   exclude_bonded_pairs=True, temperature=2)
```

This creates 100 molecules containing 200 atoms. `system.bonds` lists connected
atom indices: `(0, 1)`, `(2, 3)`, and so on.

| Pair in a two-molecule system | Interaction with exclusions enabled |
|---|---|
| 0–1, 2–3 | Bond |
| 0–2, 0–3, 1–2, 1–3 | Lennard–Jones |

## Pair potential

`pair_potential="lj"` is the default for both atoms and molecules:

`U(r) = 4ε [(σ/r)¹² − (σ/r)⁶]`

Epsilon is the unshifted well depth, and sigma is the unshifted zero crossing.
`epsilon=1, sigma=1` are defaults. Both must be positive.
The optional `"lj96"` uses

`U(r) = ε [2(σ/r)⁹ − 3(σ/r)⁶]`

For 9–6, sigma instead gives the minimum location. In both cases the pair energy
is shifted by subtracting U(cutoff) inside the cutoff, and is zero outside.

Lengths and energies use fixed course reference units. Changing sigma does not
rescale the geometry, bond or cutoff; changing epsilon does not rescale the
specified temperature. For example:

```python
sim = MDSimulation(system, epsilon=.8, sigma=1.2, cutoff=3.0)
```

Here cutoff is 2.5 times the chosen sigma. A stronger interaction or smaller
length scale may require a shorter timestep. Verify energy conservation in NVE.

## Bond choices

```python
from fys2160_md import HarmonicBond, RigidBond, Class2Bond

spring = HarmonicBond(length=.7, stiffness=100)
rigid = RigidBond(length=.7)
legacy = Class2Bond(length=.7, k2=35, k3=-40, k4=59)
```

| Bond | Physical meaning |
|---|---|
| HarmonicBond | Spring: U = stiffness/2 × (r − length)²; allows vibration |
| RigidBond | Constraint: r = length; allows translation and rotation, no stretching |
| Class2Bond | Legacy polynomial: U = k2 q² + k3 q³ + k4 q⁴, q = r − length |

Lengths and harmonic stiffness must be positive. Stiffness has units of
reference energy divided by reference length squared. Class2Bond requires positive
k2 and k4 and finite k3. Bonds are immutable; create a new system to change them.
A bond length must fit within half the shortest simulation box side.

## Why exclude bonded pairs?

With `exclude_bonded_pairs=True`, a flexible system's potential energy is the sum
of its bond energies and LJ energies **between molecules**. The bond potential
already supplies a restoring force for stretching and compression near its rest length.

With False, LJ also acts inside each molecule. The total bond + LJ potential then
determines equilibrium and stiffness. At r = 0.7 sigma the 12–6 repulsion is very
large, so simply turning exclusions off can make that starting geometry unstable.
For rigid bonds, the constraint still fixes the distance; the LJ energy is constant
within that bond. Molecular pressure includes intermolecular forces and COM motion.

Exclusion is determined by **connectivity**, not atom type. It never removes
interactions with matching atom types in other molecules.
