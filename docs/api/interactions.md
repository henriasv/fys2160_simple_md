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

## Mixtures

Species share fixed LJ reference units: sigma_ref = epsilon_ref = mass_ref = kB = 1.
Use a scalar epsilon/sigma for identical interactions or dictionaries for distinct
species. With `mixing_rule="lorentz-berthelot"` (the default and currently supported
rule), sigma_ij = (sigma_i + sigma_j)/2 and epsilon_ij = sqrt(epsilon_i epsilon_j).
Both 12–6 and 9–6 use this explicit rule. All pairs share the specified cutoff.
No cross-pair overrides are currently provided.

Masses belong to the System; they affect acceleration and thermal velocities,
not the potential function. Inspect `system.atoms.species` and `.masses`.
See [mixture examples](../examples/index.md#a-mixture-with-different-masses-and-interactions).

## A separate bond for each molecular species

```python
bond = {
    "N2": RigidBond(length=0.65),
    "O2": RigidBond(length=0.75),
}
```

Pass this mapping to `Random`, `FCC` or `System` with `molecule="diatomic"`.
Keys must match every molecular species exactly. These lengths are examples,
not physical parameters for air. Each homonuclear molecule uses the bond for its
species, including during integration and after checkpoint restart.

For flexible molecules, use separate `HarmonicBond(length=..., stiffness=...)`
or `Class2Bond(...)` entries. All entries must be rigid or all must be flexible;
mixing constrained and flexible molecules in one simulation is not supported.
Per-species mappings require homonuclear partners. Explicit heteronuclear pairs
continue to support a shared bond description. `system.bond` returns a copy of
the mapping; editing it does not change the system or a running simulation.

See the [air-like example](../examples/index.md#an-air-like-mixture-with-different-bonds).


## Isolated gravity

`pair_potential="gravity"` uses every pair, with no cutoff or periodic images:

- `U(r) = -G*m_i*m_j/sqrt(r² + a²)`
- `F_i = -G*m_i*m_j*(r_i-r_j)/(r² + a²)^(3/2)`

Here `a=softening`. With a=0 the force is exactly inverse square and potential
energy is exactly -G*m_i*m_j/r. Coincident particles are singular and rejected.
With a>0 close encounters are smoothed; no approximation truncates the far field.
Forces and energies use the same potential, integrated by velocity Verlet.

Gravity currently requires an unbonded open System and ensemble="nve". The C
solver directly sums all pairs, O(N²) per step, with a maximum of 4096 particles.
There is no Ewald solver: that would describe periodic copies, a different
physical problem. Particle masses enter both force and acceleration. G is
expressed in the same fixed reference units as the rest of the package; changing
it does not rescale coordinates or time. LJ epsilon/sigma are not applicable.

`heat_rate` provides an optional controlled energy source/sink by uniform velocity
rescaling. Its sign describes the **energy transfer**, not necessarily the final
temperature change. See [negative heat capacity](../examples/gravity.md).
