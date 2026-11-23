# Crystals and elasticity

An ionic crystal such as quartz needs three things the Lennard-Jones examples do
not: a pair potential given for each pair of species, long-range electrostatics in
a box that can be sheared, and pressure control that lets each box length find
its own size. This page shows them with alpha-quartz and the BKS force field.

```python
from fys2160_md import AlphaQuartz, Buckingham, MDSimulation

quartz = AlphaQuartz(cells=(4, 2, 3))      # 432 atoms, measured room-temperature structure
bks = {
    ("O", "O"): Buckingham(A=1388.7730, rho=0.362319, C=175.0),
    ("Si", "O"): Buckingham(A=18003.7572, rho=0.205205, C=133.5381),
}
relax = MDSimulation(
    units="metal",                  # Å, ps, eV, K, bar
    pair_potential=bks,
    charge={"Si": 2.4, "O": -1.2},
    coulomb="pppm",
    kspace_accuracy=1e-4,
    lj_shift=False,
    tail=True,                      # the -C/r^6 attraction beyond the cutoff
    ensemble="npt",
    barostat="anisotropic",         # x, y and z relax separately
    temperature=300,
    pressure=0.0,
    timestep=0.001,                 # 1 fs
    cutoff=8.0,
)
relaxed = relax.run(quartz, steps=20_000, thermo_every=100, trajectory_every=1000,
                    storage_name="quartz-relax")
```

The [executed quartz notebook](examples/quartz.md) continues from here to the
elastic constants, the speed of sound, the vibrational spectrum and the heat
capacity.

## Units

| Quantity | `units="metal"` |
|---|---|
| Length, time | Å, ps |
| Mass | g/mol |
| Energy | eV (totals for the whole box) |
| Temperature, pressure | K, bar (1 GPa = 10 000 bar) |
| Charge | elementary charge e |

These are LAMMPS `metal` units, with its constants: kB = 8.617343e-5 eV/K,
Coulomb constant 14.399645 eV Å/e², and 1 eV/Å³ = 1.6021765e6 bar. Metal-unit
simulations need explicit `temperature`, `timestep` and `cutoff`, and `pressure`
for npt/nph. Defaults that apply only to metal units: `skin=2.0` Å, `friction=5`
1/ps and `pressure_time=1` ps. Runs record `mass_density` in g/cm³.

## The BKS force field

BKS (van Beest, Kramer and van Santen, Phys. Rev. Lett. 64, 1955 (1990)) gives
silicon the charge +2.4 e and oxygen −1.2 e, and adds a
[Buckingham term](api/interactions.md#buckingham-terms-per-pair-of-species)
`A exp(−r/ρ) − C/r⁶` to the Si–O and O–O pairs; two Si atoms only repel through
their charges. The thermo table splits the potential energy into
`buckingham_energy` and `coulomb_energy`.

BKS was fitted with all pairs summed, however far apart. With a cutoff of 8 Å the
attraction `−C/r⁶` of the pairs beyond it is missing, about −1.1 GPa of pressure in
quartz, and the crystal swells: density 2.53 g/cm³ at 300 K and 0 bar, against
2.648 measured. `tail=True` (with `lj_shift=False`) adds that energy and pressure
back for a uniform density beyond the cutoff, as LAMMPS's `pair_modify tail yes`
does; the density is then 2.61 g/cm³ and the lattice constants a = 4.935 Å and
c = 5.441 Å are within 0.7% of the measured 4.913 and 5.405 Å.

Like every Buckingham potential, BKS falls without bound at short distances (inside
1.10 Å for Si–O). At room temperature the atoms never get near; a run that lets them
(a hot melt, a long timestep) stops with a message naming the two atoms.

## Anisotropic pressure control

`barostat="anisotropic"` gives each box length its own barostat, driven by its own
diagonal component of the pressure tensor (pxx for x, and so on), so that a crystal
whose c axis is stiffer than its a axis relaxes both to zero stress. The tilt of a
tilted box scales with the lengths along its row (xy and xz with x, yz with y), so
its shape stays a LAMMPS triclinic box. The equations are those of
Martyna, Tobias and Klein for each axis, with a third of the isotropic barostat's
inertia per axis; without a thermostat (`nph`) the recorded `extended_enthalpy`,
which includes `barostat_energy`, is conserved. The anisotropic barostat is
available for atoms with charges or Buckingham terms, or in real or metal units;
Lennard-Jones atoms in reduced units use the isotropic one.

## Strain, stress and elastic constants

```python
state = relaxed.state
e = 0.01
stretched, compressed = state.copy(), state.copy()
stretched.strain(xx=e)
compressed.strain(xx=-e)
fixed = MDSimulation(units="metal", pair_potential=bks, charge={"Si": 2.4, "O": -1.2},
                     coulomb="pppm", kspace_accuracy=1e-4, lj_shift=False, tail=True,
                     ensemble="nvt", temperature=300, timestep=0.001, cutoff=8.0)
plus = fixed.run(stretched, steps=10_000, thermo_every=10, trajectory_every=None, storage_name="plus")
minus = fixed.run(compressed, steps=10_000, thermo_every=10, trajectory_every=None, storage_name="minus")
C11 = -(plus.thermo["pxx"][100:].mean() - minus.thermo["pxx"][100:].mean())/(2*e)/1e4   # GPa
```

(The notebook first strains the relaxed state to its mean box, so that the
strains start from zero stress.) `state.strain` deforms the box and the crystal in it (see
[tilted boxes](api/system.md#tilted-boxes)); a shear such as `strain(yz=e)` tilts
the box. Charged atoms in a tilted box use the same PPPM electrostatics, with the
charges spread on a mesh in the box's own (fractional) coordinates. Every run of
atoms records the pressure tensor `pxx`, `pyy`, `pzz`, `pxy`, `pxz`, `pyz`; the
stress is its negative, so an elastic constant is minus the change in a pressure
component over the change in strain.

## Validation

`scripts/validate_quartz.py` compares a sheared, perturbed BKS quartz crystal of
1080 atoms (all three tilts nonzero) with LAMMPS (`buck/coul/long`, `pppm`, the
same mesh and splitting parameter, `pair_modify shift no tail yes table 0`):

| Quantity | Agreement |
|---|---|
| Potential energy | 3e-15 relative |
| Forces | 1.6e-5 eV/Å (largest force 8.5 eV/Å) |
| Off-diagonal pressure | 0.004 to 0.08 bar (of about 100 kbar) |
| Diagonal pressure | 3 bar (of 10 to 45 kbar) |

The 3 bar is understood: both codes approximate erfc the same way, but LAMMPS pairs
that approximation with the derivative of the exact erfc in its real-space force,
while ours is the exact derivative of the energy it reports. Summing the difference
over the pairs predicts the 3.0 to 3.2 bar seen. The Buckingham part alone agrees
with LAMMPS to rounding (forces 1e-12, pressures 1e-8 bar). The PPPM virial holds
the influence function fixed under strain, so in both codes it differs from the
strain derivative of the energy by a few bar; that is far below the stress
differences that give elastic constants (about 10⁴ bar for a 1% strain).

The tests also check PPPM in a tilted box against a direct Ewald sum (energy,
forces and the virial tensor), the pressure tensor against the strain derivative of
the energy, and that a tilted run is the same trajectory in pieces as in one go.
