# Physics & validation

All quantities use reduced units with reference mass, sigma, epsilon and kB
set to one. **All energies in `thermo.csv` are totals**, unlike LAMMPS's default
per-atom thermodynamic energies in LJ units. Divide by `system.N` to get energy
per atom or per molecule. `atom_density` is also recorded to expose that distinction.

- Default pair interaction for **atoms and molecules**: 12–6 LJ,
  `U(r) = 4*epsilon*((sigma/r)^12 - (sigma/r)^6)`.
- Optional `pair_potential="lj96"`: `U(r) = epsilon*(2*(sigma/r)^9 - 3*(sigma/r)^6)`.
  Here sigma is the minimum location; in 12–6 it is the zero crossing before shifting.
- Flexible default: `HarmonicBond(length=.7, stiffness=100)`,
  `U_b = stiffness/2 * (r - length)^2`.
- Rigid bond: fixed specified length, enforced with RATTLE position/velocity constraints.
- `exclude_bonded_pairs=True` omits LJ for the two connected partners. With False,
  their LJ energy and force are added to the bond; this changes the flexible
  bond's equilibrium length and stiffness. No exclusions apply between molecules.
- The original anharmonic bond remains available as `Class2Bond()`:
  `q = r - .7`, `U_b = 35*q^2 - 40*q^3 + 59*q^4`.

All inputs use **fixed reference units**. Changing epsilon or sigma adjusts the
interaction, not the unit system, bond length, temperature, box or cutoff.
The displayed particle radius follows 0.5 times the chosen sigma.

The pair potential is **shifted to zero at the cutoff**, with no tail correction;
forces below the cutoff match the stated potential. The original Atomify inputs
use an unshifted cutoff. Thus trajectories and absolute energies need not match,
and a changing number of pairs inside the cutoff also changes energy differences
slightly. The continuous potential makes energy-conservation checks meaningful.
The force still jumps at the cutoff; timestep convergence must be checked.
The recorded notebook uses the default shifted 12–6 model.
The older molecular Atomify model can be selected explicitly with `lj96` and
`Class2Bond()`; new harmonic/12–6 simulations need not reproduce its numbers.

Further deliberate differences: Langevin rather than Nosé–Hoover temperature
control, FCC molecular centres instead of random placement plus minimization,
stricter molecular timesteps, and a slower isotropic atomic barostat. These do
not constitute a bit-for-bit reproduction of the Atomify input files. Molecular models use a smaller default timestep to resolve their faster internal motion.

### Temperature, pressure and heating

`T = 2K/g`, with `g = 3N_atoms - 3` for flexible/atomic models and one fewer
DOF per rigid bond. Total centre-of-mass momentum is zero. Atomic pressure is
`(2K + sum(r_ij · F_ij))/(3V)`.

Molecular pressure uses centre-of-mass translational kinetic energy plus the
**intermolecular molecular virial**. This is an equilibrium-equivalent molecular
pressure estimator; it avoids explicit constraint-force estimates. Its instantaneous
trace differs from the atomic virial, particularly for flexible bonds. Internal
bond forces must not be omitted from an otherwise atomic pressure formula.
`Z = PV/(N_particles T)` always uses molecular count for diatomics. A finite ideal
atomic system with zero total momentum gives `Z = 1 - 1/N`, rather than exactly one.

`ensemble="nve"` selects fixed-volume dynamics without a thermostat. With
nonzero `heat_rate`, energy is deliberately changing, so this is not a strict
NVE ensemble. Similarly, NPH with heating does not conserve enthalpy.

Heating rescales the velocities each step to add `heat_rate * timestep` to the
kinetic energy, preserving zero momentum and rigid velocity constraints. For
constant volume, infer `C_V` from the slope of **total energy** versus temperature.
For constant pressure, infer `C_P` from **enthalpy** `E + P_target V`. The atomic
barostat has its own kinetic energy; the conserved quantity in unheated NPH is
`extended_enthalpy = E + P_target V + barostat_energy`. Under heating, its change
should equal the recorded `heat_added`. In NVE/NVT, `enthalpy` instead uses the
instantaneous pressure and is not a conserved quantity.

Dilute classical reference capacities per molecule / kB: atomic `3/2`, rigid
linear diatomic `5/2`, harmonic flexible diatomic **`7/2`**. The latter includes
three translational, two rotational, one vibrational kinetic and one vibrational
potential quadratic terms. The anharmonic, interacting model can depart
from 7/2. Finite system corrections, finite heating rate, and correlated samples
also matter when estimating response coefficients. The public notebook illustrates
trajectories and time series without working through a course lab.

### Integrators

The C kernel uses velocity Verlet for NVE, BAOAB Langevin splitting for NVT,
and RATTLE for the rigid bonds. Atomic pressure control uses a symmetric splitting
of isotropic MTK equations, with `a = d(log L)/dt`, `g = 3N - 3` and fixed
`W = (g+3) T_initial pressure_time²`:

```
dx/dt = v + a x
dv/dt = F/m - (1+3/g) a v
da/dt = [2K + virial - 3 P_target V + 6K/g] / W
```

NPT adds the particle Langevin thermostat; NPH removes it. Pressure control for
molecular models is deliberately rejected. The implementation checks energy,
pressure and kinetic-temperature behavior, but this initial validation does not
establish precision ensemble statistics near phase transitions.

References for the algorithms and model definitions:
[BAOAB](https://arxiv.org/abs/1203.5428),
[constrained Langevin integration](https://pmc.ncbi.nlm.nih.gov/articles/PMC4893190/),
[MTK](https://doi.org/10.1063/1.467468),
[LJ 9–6](https://docs.lammps.org/pair_class2.html),
[class2 bond](https://docs.lammps.org/bond_class2.html),
[heat-rate convention](https://docs.lammps.org/fix_heat.html).

## Saved runs and visualization

Each call replaces the result with the chosen name:

```
runs/<storage_name>/
    metadata.json       # model, units, effective settings, counts, status
    thermo.csv          # initial, sampled and final observables
    checkpoint.npz      # final arrays, RNG, step, time and barostat state
    frames/0000000.npz  # positions, unwrapped positions, velocities, box, time
```

Only existing MD run directories are replaced; arbitrary folders and symbolic
links are rejected. A `Run` refers to its folder, so overwriting that folder also
changes what an older `Run` object reads. `Run.find("runs")` remains available
for advanced comparisons; `MDSimulation.from_run(saved)` explicitly restarts the
final checkpoint.

Sampling/recording does not depend on visualization. Live display stores only
the latest snapshot and shows at most 1500 atoms. `Run.view()` samples at most
150 saved frames / 1500 atoms by default and reports those display limits.
Full-resolution frames remain on disk. The player offers orthographic and perspective
projection, with shaded sphere impostors of radius **0.5σ** in simulation
units. Their size relative to the box therefore reflects density; the radius is
an illustrative size, not a hard-sphere collision boundary. One-finger
click-and-drag rotates. Two-finger click-and-drag (secondary-button drag), or
Shift-drag, pans; two simultaneous touchscreen contacts also pan. Arrow keys
rotate a focused canvas. Two-finger scrolling zooms, with positive scroll
zooming in; buttons and the slider offer the same control. `view(zoom=1.5)` sets
the initial magnification. Reset view restores the angle, pan and zoom.
The renderer owns a separate 720 CSS-pixel-wide element inside the widget host,
so editor zoom scales it along with other notebook content. Narrow panes shrink
it to fit. Its canvas redraws at the current display resolution. After updating
the package in an open notebook, restart its kernel and rerun the visualization
cell to load the updated renderer.
Live and saved views share the same bundled WebGL 2 renderer. Ray–sphere
intersections write the actual surface depth for each pixel, so dense crystals
do not pop between whole-disc draw orders. Box edges and bonds share the depth
buffer. Graphics acceleration / WebGL 2 must be available in the notebook host. Camera controls run in
the browser and retain their state across updates while Python is busy. No
remote JavaScript is imported by the renderer. The saved player draws only while playing or
when its controls change; it has no continuous idle animation. No remote scripts
or notebook browser service are needed. For large trajectories, use
`iter_frames()` instead of loading all frames into memory. Unwrapped coordinates
include affine box motion under pressure control; account for this in diffusion
analyses. Time/steps continue across successive runs.

Interrupting a notebook run preserves the last completed state with status
`interrupted`. A numerical failure is marked `failed`, records its error and
blocks further use of that in-memory simulation; restart from a previous good
run. As with other programs, abrupt kernel termination cannot guarantee a final
checkpoint. `metadata.json` records completion so incomplete output is identifiable.

## Validation

```sh
python -m unittest discover -s tests
node tests/test_viewer.cjs
```

The 23 Python checks cover pair forces against analytical derivatives, periodic
boundaries, neighbor lists against brute-force sums, timestep convergence, rigid
constraints, thermodynamic sampling, energy addition, checkpoint restart, named
storage, pacing invariance, and both geometry constructors.

The viewer checks cover camera gestures, physical particle sizing, display scaling
and cleanup. `python tests/preview_viewer.py /tmp/depth.html` writes a browser
regression page: its GPU pixels are compared with independent ray–sphere
intersections at nearly coplanar angles, in both projections, and under reversed
draw order. This specifically guards against dense-crystal overlap artifacts.
