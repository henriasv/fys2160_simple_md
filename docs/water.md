# Water and electrostatics

The package simulates **rigid water** (three-site SPC/E, or four-site TIP4P-type
models) with **long-range electrostatics** by particle–particle particle–mesh
Ewald summation (PPPM). Water simulations use LAMMPS `real` units.

```python
from fys2160_md import IceIh, RigidWater, MDSimulation

spce = RigidWater(oh_length=1.0, hoh_angle=109.47)
system = IceIh(N=768, density=0.92, molecule="water", bond=spce, seed=1)
sim = MDSimulation(
    system,
    units="real",
    epsilon={"O": 0.15535, "H": 0.0},   # kcal/mol
    sigma={"O": 3.166, "H": 0.0},       # Å
    lj_shift=False,
    charge={"O": -0.8476, "H": 0.4238}, # e
    coulomb="pppm",
    kspace_accuracy=1e-5,
    cutoff=10.0,                        # Å, Lennard-Jones
    coulomb_cutoff=10.0,                # Å, real-space Ewald part
    ensemble="npt",
    temperature=251,                    # K
    pressure=1.0,                       # atm
    timestep=2.0,                       # fs
    friction=0.005,                     # 1/fs (Langevin thermostat)
    pressure_time=1000.0,               # fs (barostat inertia)
)
result = sim.run(steps=1000, thermo_every=10, trajectory_every=100, storage_name="ice")
```

## Units

| Quantity | `units="real"` |
|---|---|
| Length, time | Å, fs |
| Mass | g/mol |
| Energy, enthalpy | kcal/mol (totals for the whole box) |
| Temperature, pressure | K, atm |
| Charge | elementary charge e |

These are the units used in LAMMPS and in the Atomify water inputs. `thermo.csv`
records energies for the **whole simulated system**: divide by the number of
molecules, `system.N`, for kcal per mole of water (multiply by 4.184 for kJ/mol).
Real-unit runs also record `mass_density` in g/cm³. Charged runs record the split
into `lj_energy` and `coulomb_energy`.

The constants are those of LAMMPS: kB = 0.0019872067 kcal/(mol K), Coulomb
constant 332.06371 kcal Å/(mol e²), and 1 kcal/(mol Å³) = 68568.415 atm.
Reduced Lennard-Jones units remain the default for all other systems.

Real-unit simulations require explicit `temperature`, `timestep` and `cutoff`,
and `pressure` for npt/nph. Defaults that apply only to real units: `skin=2.0` Å,
`friction=0.005` 1/fs and `pressure_time=1000` fs.

## Rigid water geometry

`RigidWater(oh_length, hoh_angle, m_distance=0)` fixes both O–H distances and the
H–O–H angle. Each molecule is stored as three atoms, O then H then H, with
species `"O"` and `"H"`. Three distance constraints per molecule (O–H, O–H and
H–H) are enforced with RATTLE to 10⁻¹³ relative precision. A molecule has six
degrees of freedom: three translations and three rotations.

With `m_distance > 0` the model has a **massless M site** on the H–O–H bisector,
this far from O. As in LAMMPS `lj/cut/tip4p/long`, the charge given for O is
placed on M; O keeps the Lennard-Jones interaction. M is recomputed from the
atoms every step, and its force is shared between O and the two H atoms.

| Model | `oh_length` (Å) | `hoh_angle` (°) | `m_distance` (Å) | q(H) (e) | ε(O) (kcal/mol) | σ(O) (Å) | Melting point of the model |
|---|---|---|---|---|---|---|---|
| SPC/E | 1.0 | 109.47 | 0 | 0.4238 | 0.1553 | 3.166 | ≈ 215 K |
| TIP4P | 0.9572 | 104.52 | 0.15 | 0.52 | 0.1550 | 3.15365 | ≈ 232 K |
| TIP4P/2005 | 0.9572 | 104.52 | 0.1546 | 0.5564 | 0.1852 | 3.1589 | ≈ 250 K |
| TIP4P/Ice | 0.9572 | 104.52 | 0.1577 | 0.5897 | 0.2108 | 3.1668 | ≈ 270 K |

The O charge is −2 q(H). H has no Lennard-Jones interaction: give it `epsilon=0`
and `sigma=0`. The Atomify water inputs use SPC/E with ε(O) = 0.15535 kcal/mol.
Melting points are published estimates for each model (Vega and co-workers),
not properties computed here.

```python
tip4p2005 = RigidWater(oh_length=0.9572, hoh_angle=104.52, m_distance=0.1546)
system = IceIh(N=768, density=0.92, molecule="water", bond=tip4p2005)
sim = MDSimulation(system, units="real",
                   epsilon={"O": 0.1852, "H": 0.0}, sigma={"O": 3.1589, "H": 0.0},
                   lj_shift=False, charge={"O": -1.1128, "H": 0.5564},
                   coulomb="pppm", kspace_accuracy=1e-5, cutoff=10.0, coulomb_cutoff=10.0,
                   ensemble="npt", temperature=251, pressure=1.0, timestep=2.0)
```

`lj_shift=False` uses the plain truncated Lennard-Jones potential, as LAMMPS
`lj/cut/coul/long` does by default. `lj_shift=True` subtracts the energy at the
cutoff, which makes energy conservation easier to check. The two give the same
forces; energies differ by a constant for each pair inside the cutoff.

## Starting configurations

`IceIh(N, density=0.92, molecule="water", bond=..., masses=None, seed=...)` builds
proton-disordered hexagonal ice. N molecules fill nx·ny·nz orthorhombic cells of
eight molecules, chosen as close to cubic as possible: N=768 gives 6×4×4 cells
(27×31×29 Å at 0.92 g/cm³) and N=2880 gives 10×6×6, the Atomify box size.
Other counts are rejected with nearby valid suggestions. The hydrogens obey the
ice rules (each oxygen donates and accepts two hydrogen bonds), and loop flips
remove the net dipole. Default masses are O 15.9994 and H 1.00794 g/mol.

`Random(N, box=..., molecule="water", bond=..., min_distance=1.8)` places randomly
oriented molecules for a liquid start; `box` is in Å. Relax a random start with
a short timestep before production.

A perfect ice crystal can be **superheated**: without surfaces or defects it
needs time to start melting. 768 SPC/E molecules started from `IceIh` stay
crystalline for about 25 ps at 310 K, both here and in LAMMPS, and then melt.
Larger boxes and higher temperatures melt sooner.

## Electrostatics: PPPM

The Coulomb energy of the periodic system is split into a short-range part,
summed directly within `coulomb_cutoff`, and a smooth long-range part computed on
a mesh with fast Fourier transforms. The method follows LAMMPS
`kspace_style pppm`: fifth-order charge assignment, ik differentiation and the
optimal influence function of Hockney and Eastwood.

`kspace_accuracy` is the relative RMS force accuracy, compared with the force
between two unit charges 1 Å apart. From it and `coulomb_cutoff`, the Ewald
splitting parameter and mesh are chosen as in LAMMPS. They are selected again
at the start of every `run` and whenever the box changes by more than 5%; the
choice is recorded in `metadata["pppm"]`. The total Ewald energy does not depend
on the split, so changing `coulomb_cutoff` or `kspace_accuracy` between runs
changes only the numerical error, not the model.

The short-range part uses the Abramowitz–Stegun approximation of erfc (error
below 1.5×10⁻⁷, as LAMMPS without tables) together with its exact derivative, so
forces are gradients of the computed energy. Pairs within a molecule are excluded:
their mesh contribution is subtracted exactly. The total charge must be zero.

Mesh cost grows with the box volume. For dilute vapour in a greatly enlarged box,
lower `kspace_accuracy` (for example 10⁻⁴) or increase `coulomb_cutoff` in that run.

## Pressure, enthalpy and box expansion

Pressure uses the **molecular virial**: centre-of-mass kinetic energy and the
virial of the forces acting on whole molecules. Constraint forces then do not
enter. Under npt/nph, `enthalpy` is E + P_target·V. LAMMPS's `enthalpy` uses the
instantaneous pressure instead; the two agree on average, but the instantaneous
version fluctuates far more.

The npt/nph barostat scales centre-of-mass positions and velocities isotropically,
so molecules keep their shape. Its inertia is 3N_mol·kB·T·pressure_time².

`run(..., deform=s)` (nve or nvt only) changes every box length linearly in time,
to s times its starting value at the end of the run. Whole molecules move with
their centres of mass, so they are never stretched. This replaces LAMMPS
`fix deform ... scale s` with its default atom remapping, which stretched
constrained molecules and caused large temperature spikes.

## Validation

- **PPPM against a direct Ewald sum** for random charges: RMS force error below
  the requested accuracy; the NaCl Madelung constant is reproduced to 10⁻⁵.
- **Against LAMMPS (30 Mar 2026)** with identical mesh and splitting parameter:
  chosen meshes and splitting parameters agree. Total energies agree to 13
  significant digits for point charges. For 768 SPC/E and TIP4P/2005 molecules in
  ice, every energy component agrees to about 10⁻¹² relative, forces to about
  2×10⁻⁶, and the virial to about 10⁻⁵.
- **Dynamics:** constraints hold to 10⁻¹³ relative, NPH conserves the extended
  enthalpy, and trajectories are identical for any `thermo_every` and after
  resuming a saved run. Started from the same ice, LAMMPS and this package melt
  on the same timescale at 310 K.

The tests are `tests/test_pppm.py` and `tests/test_water.py`.

## Speed

One CPU core, 2 fs timestep, `cutoff=coulomb_cutoff=10` Å, `kspace_accuracy=1e-5`:

| System | ms per step |
|---|---|
| 768 SPC/E molecules, ice | 11 |
| 768 TIP4P/2005 molecules, ice | 12 |
| 2880 SPC/E molecules, ice | 43 |
| 768 molecules in a 10× larger box, `kspace_accuracy=1e-4` | 31 |

LAMMPS on one core takes about 9 ms per step for the first case.
