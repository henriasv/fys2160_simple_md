# FYS2160 Simple MD

This repository contains the standalone course MD package. Python import:
`fys2160_md`; distribution: `fys2160-simple-md`.

- Geometry/state: `src/fys2160_md/{geometry,system}.py`.
- Runner: `sim = MDSimulation(...)` holds interactions and default controls only;
  `run = sim.run(state, ...)` evolves a copy of a System and `run.state` is the
  final state (clock, thermostat RNG, barostat rate, PPPM choice, heat added),
  so chained runs equal one long run. Run controls never carry over. Native solver: `_core.c`
  (LJ/diatomic/gravity, keep byte-identical) and `molecular.c` + `pppm.c` + `fft.c`
  (charged atoms and rigid water, LAMMPS `real` units).
- `tests/test_golden.py` freezes short trajectories of the old solvers; regenerate
  only for an intentional physics change.
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
- Raise `src/fys2160_md/_version.py` for every release: pip treats a git install
  with an unchanged version as already satisfied, even at a new pinned commit.
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

2026-10-09 (branch water-pppm): rigid water and PPPM for Lab 2. molecule="water"
with RigidWater(oh_length, hoh_angle, m_distance) gives SPC/E (m_distance=0) or
TIP4P-type models (O charge on the M site, as LAMMPS lj/cut/tip4p/long). Systems
carry units ("reduced" default, "real" = LAMMPS real); real units require explicit
temperature, timestep, cutoff (and pressure for npt). charge + coulomb="pppm"
select PPPM written from Hockney-Eastwood/Deserno-Holm, matching LAMMPS pppm
(ik, order 5, same mesh/g selection; energies agree to ~1e-12 with LAMMPS at the
same mesh, forces ~2e-6). The mesh choice is kept in the state and re-chosen
when kspace_accuracy/coulomb_cutoff change or the box changes >5%, so
trajectories are bitwise independent of thermo_every and of run splitting. Molecular MTK barostat (COM scaling), run(deform=s)
with COM remap (nve/nvt), IceIh(N=8nx ny nz) with ice rules and zero dipole,
Random(molecule="water"), lj_shift option, viewer bonds/element colours/real-unit
labels. Own radix-2/3/5 FFT instead of a vendored library. LAMMPS is a behavioural
reference only (GPL; nothing copied). Homebrew LAMMPS segfaults in PPPM; cross-checks
used a serial build from hsp_hymcp/lammps. 768 SPC/E waters: 11 ms/step on one core
(LAMMPS serial 9 ms). Perfect 768-molecule ice superheats ~25 ps at 310 K in both
codes, so the private Lab 2 notebook (fys2160-groupsession-materials/labs/lab2)
melts at 350 K; expansion/vapour use kspace_accuracy=1e-4. examples/water.ipynb is
an independent liquid-structure demo (SPC/E vs TIP4P/2005 density and g_OO).
61 Python tests and viewer controls pass.

2026-10-09: state is separated from the simulation. MDSimulation(...) takes no
System (passing one raises a TypeError explaining the new style); run(state, ...)
copies the state and run.state (from the checkpoint; also via Run.load) is the
output. System now carries time, step, thermostat RNG, barostat rate, PPPM mesh
choice and cumulative heat_added, so chained runs reproduce one long run bitwise.
Run controls apply to that run only (no carry-over); from_run and the Simulation
alias are removed. Tests use MDSimulation(...)._engine(state) for low-level
checks. Golden trajectories and all example outputs are unchanged by the
refactor. The private Lab 2 notebook passes each stage's .state explicitly.


2026-10-09: System.with_box(box) copies a state into a larger periodic box with
the contents centred and molecules whole (no stretching); enlarging one side
gives a liquid slab with vacuum for evaporation (direct coexistence), the
standard alternative to deform. 768 SPC/E waters in a 5x box: ~16 ms/step;
vapour holds ~3-6 molecules at 500 K, almost none at 400 K. The private Lab 2
notebook now uses ice -> melting (350 K) -> slab evaporation at 500 K.

2026-10-23: Live and saved views render with OMOVI (andeplane/omovi 0.25.0,
three 0.181.2, loaded from jsDelivr +esm at view time; Anders Hafreager allows any
use). viewer.js draws with the built-in renderer first, then switches; offline or
renderer="basic" keeps the built-in one. omovi looks colours/radii up by particle
id (not type), so they are set per atom; element-named species use omovi
AtomTypes colours. omovi 0.25 caches materials in a module-level map, so each view
imports omovi under its own URL (?view=N) to avoid shared colours; an upstream fix
would remove this. IntersectionObserver sets visualizer.idle off screen. Elongated
boxes (slabs) get a camera across the long axis. Default projection is now
perspective (SSAO needs it; orthographic disables SSAO).
omovi draws opaque (alpha:false in its Visualizer), so the scene background is set
to the first non-transparent ancestor background (re-read on light/dark switches);
the built-in canvas is transparent. True transparency would need an omovi option.

2026-10-23 (0.2.1): OMOVI camera and lighting. The camera distance is computed so
that all eight box corners are in view (5% margin) at zoom 1, with a 35° field of
view (omovi's 60° looked lopsided and stretched spheres at the edges). zoom < 1 steps
back; zoom > 1 narrows the field of view, so the camera never enters the box (the
old distance/zoom put it among the atoms). The quickstart crystal and gas views no
longer pass zoom=1.4, which was tuned for the built-in renderer. The orthographic frustum matches
the perspective scale at the box centre (omovi uses half the distance). The
glTexStorage2D zero-size warnings come from omovi's constructor and are harmless.

2026-10-24 (0.2.2): omovi places its point light half a box diagonal from the
centre (fine for Å-sized systems); in dense reduced-unit boxes that is ~1σ in front
of the nearest face, so spheres showed hard half-moon shading. viewer.js replaces
visualizer.updatePointLightPosition with a key light 0.3·distance above and left of
the camera, decay 0 and intensity 2.5, the same for every box size. This overrides
an omovi-internal method (fine for the pinned 0.25.0); an upstream light option
would be cleaner.

2026-10-26 (0.2.3): space plays and pauses a focused saved view. The OMOVI container
is focusable (tabIndex 0) and takes focus on pointerdown; the built-in canvas
already had tabindex. Key repeat is ignored and live views have no playback.
