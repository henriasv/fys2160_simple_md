# FYS2160 Simple MD

This repository contains the standalone course MD package. Python import:
`fys2160_md`; distribution: `fys2160-simple-md`.

- Geometry/state: `src/fys2160_md/{geometry,system}.py`.
- Runner: `sim = MDSimulation(system, ...)`, then `sim.run(...)`; native solver: `_core.c`.
- Shared live/saved WebGL renderer: `viewer.js`.
- Docs: `docs/`, `mkdocs.yml`; static executed notebook built by
  `scripts/build_docs.py`. GitHub Actions publishes Pages on main.
- Surface ensemble, temperature, timestep, pair parameters, cutoff, step count and
  recording intervals in complete examples, even when they match defaults.
- Keep student examples short: no block averaging or run/fit wrapper functions.
- Reusing storage_name intentionally overwrites a recognized saved run.
- Exact cubic FCC counts N=4n³; reject invalid counts with suggestions.
- Never change geometry to satisfy a self-avoidance request silently.
- Tests: `python -m unittest discover -s tests` and `node tests/test_viewer.cjs`.
- Before publishing docs, build with `python scripts/build_docs.py --execute`
  then `python -m mkdocs build --strict` and check the static trajectory players.

The public examples must be independent demonstrations, never course lab
solutions. Private teaching materials and historical lab notebooks remain in
fys2160-groupsession-materials; do not copy them into this repository or site. Update this handoff when finishing.

## Handoff — 2026-09-15

MDSimulation now owns a copy of a geometry-only System. FCC/Random take
molecule="diatomic" and explicit HarmonicBond/RigidBond; Class2Bond retains the
old anharmonic model. Pair potential defaults to 12–6 for every system, with
optional lj96, epsilon/sigma, and explicit exclude_bonded_pairs. Checkpoints
record all parameters and old checkpoints retain the legacy force field.
32 solver/API tests pass, including numerical gradients and molecular pressure,
harmonic motion, rigid constraints and restart; viewer controls also pass.
The dilute free-rotation test now uses rho=1e-8: the former density allowed
transient collisions after velocity seeds were separated from geometry seeds.
Docs use pyplot directly, and the static notebook auto-sizes with the page.
The public notebook remains a separate crystal/gas showcase, not a lab solution.

Later the same day: species mappings specify exact counts (N is inferred), and
geometry accepts per-species masses. MDSimulation epsilon/sigma accept scalars
or complete species mappings, with explicit Lorentz–Berthelot mixing. Everything
uses shared LJ reduced reference units. Molecular geometry makes homonuclear
pairs; explicit System arrays also support heteronuclear pairs with unequal
masses. Per-species colors and radii work in live/saved viewers. The executed
notebook now includes mixture and bonded demonstrations plus direct pyplot.
37 tests pass, including mixture forces, mass acceleration, equipartition and
checkpoint restoration. Historical lab files were removed with targeted filtering;
private bundle backups are in the teaching repository's .git/md-history-backups.

2026-09-15: Standard simulation parameters are explicit in the recorded notebook,
recipes, installation, overview, geometry examples and README. Brief explanations
distinguish integration timestep, total duration and output intervals. Numerical
values and package defaults are unchanged.

2026-09-15: Per-species bond mappings are supported throughout geometry, C forces,
RATTLE lengths and checkpoint restart. All bonds must be rigid or all flexible;
per-species maps require homonuclear partners. Shared bonds still support
heteronuclear explicit arrays. The public notebook includes an explicitly
illustrative air-like mixture. Output options are now thermo_every and
trajectory_every; sample_every/save_every remain aliases, with conflicting names
rejected. Canonical names are used in metadata and all public examples.

The worked NPT example examples/isobar.ipynb varies T=1.5,2,2.5,3 at P=0.05,
with separate equilibration and production, disk reloading, pyplot P/V/Z curves,
production traces and playback. It estimates Z using target P,T and mean V;
measured P,T are separate control checks. Both notebooks execute in build_docs.py.
41 solver/API tests and viewer controls pass; per-species spring and constraint
checks include save/restart and recording-name compatibility.


2026-09-16: Added isolated all-pairs gravity in the C solver, selected with
pair_potential="gravity" on System(boundary="open"), ensemble="nve", cutoff=None.
G and softening are explicit, with softening=0 giving exact Newtonian gravity.
No periodic sums, walls or truncation. Open coordinates persist through copy,
saving and restart; pressure is NaN and no gas density/Z/enthalpy is reported.
Plummer generates approximate stationary cluster positions and velocities, with
G/softening matched to the simulation. It supplies a viewing frame, not a box.
Gravity supports at most 4096 particles and uses O(N²) direct pairs.

examples/gravity.ipynb removes 8 energy units from a 256-particle cluster and
compares relaxed averages, with matplotlib plots, energy/virial checks and a
saved player. The published protocol shows T 0.0971→0.1144 and half-mass radius
1.343→1.061. Half-dt, smaller-softening and no-sink comparisons are reproducible
with scripts/validate_gravity.py and documented in docs/physics.md. Interpretation
is a dynamical illustration, not a precision equilibrium heat-capacity curve.
46 tests pass, including exact two-body force/Kepler orbit, energy removal,
open checkpoint and native-buffer safety. Static build now executes three
notebooks. Gravity playback omits box edges; particle_radius controls display
markers only (default 0.035). Existing LJ view radii/default cutoff remain as before.
