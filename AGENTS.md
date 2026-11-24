# FYS2160 Simple MD

This repository contains the standalone course MD package. Python import:
`fys2160_md`; distribution: `fys2160-simple-md`.

- Geometry/state: `src/fys2160_md/{geometry,system}.py`.
- Runner: `sim = MDSimulation(...)` holds interactions and default controls only;
  `run = sim.run(state, ...)` evolves a copy of a System and `run.state` is the
  final state (clock, thermostat RNG, barostat rate, PPPM choice, heat added),
  so chained runs equal one long run. Run controls never carry over. Native solver: `_core.c`
  (LJ/diatomic/gravity, keep byte-identical) and `molecular.c` + `pppm.c` + `fft.c`
  (charged atoms, Buckingham terms and rigid water, LAMMPS `real` and `metal` units;
  tilted boxes and anisotropic pressure control for atoms). `tilt.h` is shared.
- A run in memory: `live = sim.start(state, ...)` returns a `LiveRun` with
  `advance(steps)`, `observables`, `positions`, `velocities`, `box`, `step`, `time`,
  `set(temperature=...)`, `settings` and `state`. `run` is that loop with files.
  New observables and controls (walls, pistons, external potentials) belong in
  `observables` and as keywords of `start`/`set`, not in new methods.
- `state.with_velocities(v)` and `state.with_box(box)` return changed copies.
  The box is the one editable part of a state (`state.box[0] *= 2`, for students);
  runs never change the state they are given, and `Run.state` is cached per Run so
  that an edit of `run.state` reaches the next run.
- Two dimensions and walls (atoms in reduced units only): a box with two lengths
  is flat; `boundary="walls"` or one word per axis. The solver always works on
  (N,3) arrays; a flat system has z = 0 and box[2] = 0 inside `System`/`_Engine`
  and two columns in everything public. Walls and pistons are rows of the
  solver's table of one-body terms (`walls.py`, `external_forces` in `_core.c`);
  a new external potential is a new row kind there, reported like a wall.
- The solver keeps a workspace per run (neighbour list, forces, the state it last
  returned) and reuses it only when its inputs are bit-for-bit unchanged. A new
  force input must be added to that comparison (`scalars`/tables in `_core.c`,
  `inputs` in `molecular.c`), or stale forces will be reused.
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
- Releases: `.github/workflows/wheels.yml` builds and tests wheels for Windows,
  macOS, Linux and Pyodide with cibuildwheel (settings under `[tool.cibuildwheel]`
  in `pyproject.toml`) and publishes to PyPI only for a pushed tag `vX.Y.Z` that
  matches `_version.py`.
- `import fys2160_md` must not import pandas, matplotlib, IPython or anywidget: the
  browser build has no pandas. `tests/test_imports.py` checks this; the test suite
  itself needs pandas and IPython.
- The C sources stay strict C99 without variable-length arrays, `complex.h`, `M_PI`
  or compiler extensions (MSVC), and must not assume a 64-bit `long` or `size_t`
  (Windows, WebAssembly).

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

2026-10-26 (0.2.4): redesigned view controls (VIEW_MARKUP/VIEW_STYLE in viewer.js).
The view and a footer bar form a square-cornered card: row 1 is playback (accent
play/pause, frame slider, "i / n" counter; hidden for live views), row 2 the
readout with species colour chips and tools (speed 0.5-4x, zoom buttons for the
built-in renderer only, reset, projection, "?" which reveals the note). Live views
show a status pill (running/completed/interrupted/failed). Colours derive from the
page text colour (color-mix with currentColor), so light and dark hosts both work;
root colour now inherits. The orthographic OMOVI frustum fits the box corners plus
the largest atom radius. Reset (and the initial view) is the angled view in
perspective and face-on to the xz plane (z up; z across for z-elongated boxes) in
orthographic. omovi's orbit controls are y-up, so the camera sits 1e-6 off the pole.

2026-10-26 (0.2.5): a stretched OMOVI view was reported on the desktop docs page
(not reproduced here). omovi resizes only inside its render loop, which the
IntersectionObserver pause stops, and it falls back to its canvas's 640x480 CSS size
when its element has no height. Now the container is in the page before the
Visualizer is built, the omovi canvas fills the container, and a ResizeObserver,
show() and pointer input wake the loop for two frames even while paused. Verified
with an observer that always reports "not visible".

2026-10-26 (0.2.6): playback speed is a menu (0.25x to 8x of 10 frames/s) instead
of a cycling button, and < > step it in a focused saved view. The projection menu
is found by [data-projection], since it is no longer the only select.

2026-10-27 (0.2.7): the view controls keep their layout during playback and live
runs. Readout numbers (step, time, frame counter) are padded with figure spaces to
the widest value of the run: saved views use the frame maxima, live views get
final_step/final_time from LiveView.start(sim, steps). Live T and P can widen only
(at most a few times). The status pill has a fixed minimum width; both menus
already size to their widest option.

2026-11-19 (0.2.8, branch claude/wheels-and-optional-pandas-5e75f6): packaging for
PyPI and the browser. pandas is imported inside `Run.thermo` (clear ImportError with
the numpy alternative when missing) and is a dependency everywhere except Emscripten
(`sys_platform != 'emscripten'`). `wheels.yml` has four jobs: native wheels on six
hosted runners (cp310-cp315; no win32, no free-threaded, no cp310 on Windows arm64,
which follow from where numpy and pandas publish wheels), the Pyodide wheel
(cp314, pyemscripten_2026_0_wasm32) plus `tests/pyodide_smoke.mjs` under Node with
Pyodide 314.0.8 and no pandas, an sdist, and trusted publishing for a version tag or a
manual run with "publish" ticked. Wheels no longer carry the C sources; the sdist
carries the tests and golden fixture. C audit for MSVC and wasm32: no source change
was needed except a guard in `pppm_init` against a mesh size that overflows a 32-bit
`size_t`; all six `PyMethodDef` entries match their flags. `viewer.js` is read as
UTF-8 (Windows would otherwise garble its non-ASCII labels).
Verified on this Mac only: native wheel, 64 tests (numpy 2.5.5, pandas 3.1.1,
Python 3.14); the Pyodide job run locally through cibuildwheel 4.3.5 (63 tests pass
in Pyodide, one skipped that needs a subprocess); the smoke script. The Windows and
Linux jobs and the publish job have never run. Golden trajectories differ from the
frozen ones by up to 2.4e-12 in Pyodide against a tolerance of 1e-11, so a platform
with another math library (MSVC, musl) may need a looser tolerance. The extension
already compiles unchanged against `Py_LIMITED_API=0x030B0000`, so `cp311-abi3`
wheels need only build settings; not done. The repository has no LICENSE file.

2026-11-19 (branch chunked-interface-neighbour-lists): a run can be advanced in
memory by its caller. `live = sim.start(state, ...)` returns a `LiveRun`:
`advance(steps)`, `observables` (the thermo.csv columns as a dict), read-only
snapshots `positions`, `velocities` and `box`, `step`, `time`,
`set(timestep=, ensemble=, temperature=, pressure=, heat_rate=)`, `settings`, and
`state` (a System that `sim.run` or `sim.start` continues). No files, threads or
Jupyter; it runs unchanged in Pyodide. `sim.run` starts a LiveRun and advances
it, so the two share one code path; saved output is unchanged. `set` gives the
controls a new run from the current state would have (barostat rate and inertia
included). `advance` re-raises a keyboard interrupt after whole steps; after any
other exception the run refuses to continue. `deform` stays a `run` option.
`state.with_velocities(v)` returns a copy with other velocities and the same
clock and random sequence (reverse the motion, or perturb one atom; total
momentum must stay zero). Saved frames already held velocities.

The solver keeps a workspace per run (a PyCapsule made by `_core.workspace()` or
`_core.molecular_workspace()`): the neighbour list with its reference positions
and box, the forces, and the positions, box and force parameters as the previous
call left them. A call continues from it only if those are bit-for-bit unchanged;
otherwise it rebuilds the list and recomputes the forces. The LJ list is built
with cells of reach 1 or 2 (chosen by a measured cost estimate in
`choose_cells`) or by checking all pairs in small boxes; all three give the same
pairs in the same order, so trajectories are bit-identical to 0.2.7 (checked on
18 cases, including cells, NPT, mixtures, diatomics, gravity, water, salt).
Only `time` and `heat_added` are summed per call and can differ in the last
digit between chunkings. The pressure of diatomic molecules is computed when the
observables are read (a zero-step call with `measure=1`). `molecular.c` keeps
its molecule list, PPPM mesh and forces the same way; its list build is still
all-pairs. Gravity keeps only its forces. `_core.neighbor_pairs` is a test helper.

Native steps per second at 1 and at 500 steps per call (Apple M5 Pro, NVE from an
equilibrated state, cutoff 2.5), before -> after: gas rho=0.1 T=2, 256 atoms
15k/67k -> 70k/85k; 500: 4.8k/26k -> 33k/39k; 2048: 400/2460 -> 7740/8590; 4000:
113/695 -> 3650/4030. Liquid rho=0.8 T=1, 256: 5.2k/14.0k -> 14.5k/15.1k; 500:
1670/6380 -> 6450/6650; 2048: 206/1031 -> 1326/1350; 4000: 74/414 -> 654/662.
In the liquid the force loop is now about 85% of a step; it was left untouched to
keep the summation order. SPC/E water, 768 molecules: 13.3 -> 6.2 ms/step at one
step per call. 72 Python tests and the viewer controls pass; test_chunks.py and
the golden trajectories also pass inside Pyodide 314.0.8 under Node. Not bumped:
`_version.py` (raise it with the release that contains this).

2026-11-19 (0.3.0, branch claude/2d-walls-piston-e91426): flat systems, walls and a
piston, for the week 34 session. Public surface: `Random(N=200, box=(100, 100),
boundary="walls")` (two lengths, or `dimension=2`, make a flat system; `boundary` is
"periodic", "walls", or one of the two per axis, also on `FCC` and `System`);
`state.dimension`, `state.walls` (names: left/right for x, bottom/top for y, back/front
for z, at 0 and at the box length); `MDSimulation(wall_stiffness=1e4)`;
`Piston(wall, to, speed)` as the keyword `piston` of `run`, `start` and `set`;
observables `wall_force_<name>` (now, positive outwards), `wall_impulse_<name>`
(momentum delivered so far), `work` (done on the atoms by moving walls) and
`piston_position`. `live.settings` has `piston` only for a state with walls.
`with_box` between walls keeps the atoms and moves the far wall. Docs: `docs/walls.md`.

How it is built. `_core.advance` takes three more optional arguments: the dimension,
a float64 (terms, 5) table of one-body terms (kind, axis, stiffness, speed, target) and
a (terms, 3) array the solver writes (force now, impulse, work). A direction without
periodic images (walls, or z of a flat system) gets the period `NO_IMAGE` (1e30), so
`forces` and the list build keep their unbranched minimum image, and a periodic 3D run
is bit-identical to 0.2.8 (golden trajectories and the two live scenes checked). `drift`
wraps only periodic directions; the thermostat adds noise in `dimension` directions and
removes the total momentum only along periodic ones; `System.dof` counts likewise, so
`T = K/N` exactly in a flat box with walls. Cells: a flat system has one layer; an atom
beyond a wall goes into the outermost cell. The wall is U = stiffness*d^3/3 for an atom a
depth d beyond it. A spring wall (d^2) was tried first: its kink in the force costs
energy accuracy at every step that straddles the wall (rms energy change 0.64 against
0.05 over 1000 time units for 200 atoms at T* = 2, at equal hardness). A piston is an
upper wall (kind 1) with speed > 0: `move_pistons` changes box[axis] each step, towards
the target, statelessly (the schedule is a function of the box and the controls), so
chunked and chained runs are bitwise equal. Impulse and work use the trapezoid of the
force before and after the step, which is exactly what the two half kicks give the atoms.
The external table and the wall forces are part of the workspace comparison. The
impulses and the work are carried in the state (`System._wall_sums`, checkpoint
`wall_sums`), like `heat_added`.

Not done, on purpose: pressure control, molecules, charges and real units in a flat box
or between walls (rejected with a message); a piston on a wall at 0 (the box origin is
fixed); several pistons at once (the solver could, `piston=` takes one); 2D lattices;
external potentials other than walls (gravity, a bias, a double well, a step: add a row
kind in `external_forces`, an object in `walls.py` that makes the row, and observables
under its name); a thermostat on a region. The `pressure` observable between walls is
the pair-virial pressure and leaves out the walls' own virial; it agrees with the wall
force per length on average. The Langevin thermostat at the default friction 1 turns a
dilute gas into a random walk (correct averages, slow wall-force statistics): the week 34
scenes use `friction=0.05`.

Viewer: a flat system has box thickness 0 and opens face-on to the xy plane in either
projection. Across an axis with walls the view is centred on the first frame's box, not
each frame's (`centre`, `reach` in `viewer.js`), so the wall at 0 stays put and a piston
moves; periodic boxes stay centred as before. Checked in a browser with OMOVI.

`examples/week34/` holds drafts of the two scenes for the session pages (`box.py`:
`run_box`, `wall_force`; `piston.py`: `compress`), in the style of the live scenes in
thermophysics-interactive-book, and `page_numbers.py`, which prints what they show.
`scripts/validate_walls.py` prints the validation numbers of `docs/walls.md`;
`scripts/pyodide_week34.mjs WHEEL` runs the scenes in Pyodide under Node.

Measured (Apple M5 Pro, 8 starting states, `scripts/validate_walls.py`): energy between
walls changes by 1 to 2 parts in 10 000 of the kinetic energy over 1000 time units;
PV/NkT from the force on all four walls 1.015, 1.020, 1.005, 1.010 against 1.009,
1.019, 1.005, 1.005 expected from the second virial coefficient; one 2000-time-unit
mean of the force on one wall scatters by 1.5 to 3.5% (about one hit per time unit);
piston 100 -> 50 at speed 0.05: T ratio 2.03 +- 0.02, |change of energy - work| at most
0.13 of 426. 200 atoms in 100 x 100 run at 147 000 steps per second with the thermostat
and 360 000 without (native); the week 34 scenes take 2.9 s (wall force, 420 000 steps)
and 0.6 s (piston) natively, 4.7 s and 0.9 s in Pyodide 314.0.8 under Node. Native and
browser trajectories differ (rounding in the last digit, chaotic motion); their
averages agree. 91 Python tests and the viewer checks pass; the Pyodide smoke test now
also runs a piston. Thermostat noise is one sequence per `seed`: runs that should be
independent need different seeds (`sim.start(state, seed=...)`), not only different
starting positions.

2026-11-23: `state.box` is the state's own length array (an ndarray subclass `Box`
with `.lengths` and `.volume`), so `state.box[0] *= 2` changes the box the next run
starts with. Positions are kept; the next run copies the state through the `System`
constructor, which wraps periodic axes and refuses atoms more than 1 beyond a wall
(the message points to `Piston`). A sudden squeeze of a dense periodic box overlaps
atoms and blows the run up (T* ~ 1e38 in a test; the solver does not stop it).
`Run.state` is cached on the Run object. `tests/test_box_edit.py`. Not released:
`_version.py` is unchanged. The dev `.venv` extension had to be rebuilt in place
(`python setup.py build_ext --inplace`) after the merged solver changes.

2026-11-23 (0.4.0, branch claude/measure-and-sampled-arrays): sampled arrays and
per-step measurements. `run.positions`, `.velocities`, `.unwrapped` (saved steps,
atoms, D), `.boxes`, `.steps`, `.times` stack the saved frames on first use and are
cached read-only (refuses above 2 GB, pointing to `iter_frames`). `measure={name:
function}` on `run` and `start`: the function takes one argument, the result,
and returns one value per step. During the run that argument is a `Frames` block
(simulation.py) with the attributes of a saved run's frames (positions,
velocities, unwrapped, boxes, steps, times; read-only), so `left_half(run)` on a
saved run and `measure={"N_left": left_half}` share one function. The solvers take an
optional `record` tuple (`record.h`: positions, velocities, unwrapped (steps,N,3),
box (steps,3)) and copy each step into it; `LiveRun._advance_measuring` advances in
blocks of at most 256 steps (each record array below 16 MB) and evaluates every
function once per block. Thermo rows get the mean since the previous row
(`LiveRun._row`), `measure.npz` the full series (`run.measured`,
`run.measured_steps`); the first value is the starting state. Recording does not
change trajectories (test). All four arrays are always recorded. Cost: 200 flat atoms, 113k -> 97k steps/s at
thermo_every=10 and 141k -> 133k at 100, with N_left. The same function applied to run.positions gives the in-run
values at the saved steps; the docs present the post-processing way first. The
Pyodide smoke test has a measurement, but was not run here (no wasm build).

2026-11-23 (0.4.0, same branch): still pictures. `state.draw()` and `run.draw(frame=-1)`
(`drawing.py`, matplotlib imported on use) return a matplotlib Axes: per atom a rim,
body and highlight ellipse in one EllipseCollection (units "xy"), interleaved and
sorted far to near, far atoms faded to the background; box edges split into behind
(under the atoms) and in front; bonds as lines with molecules made whole by minimum
image. Views, colours and radii follow the player (`_rotation` mirrors viewer.js:
angled, across a slab's long axis, or face-on xz in orthographic). 13 500 atoms draw
and save at 150 dpi in about 0.1 s. `tests/test_drawing.py`.
Views for `draw`: `view="angled"|"front"|"side"|"top"` or `(elevation, azimuth)` in
degrees (`VIEWS`, `_look`; matplotlib's convention, z up); without a view,
orthographic gives "front" and a z-long box lies across. x/y/z arrows on a soft
patch in the lower left (`_arrows`; a circle with a dot or cross for an axis along
the line of sight); `arrows=False` omits them. `scripts/build_docs.py` makes
`docs/images/draw_views.png` and `draw_systems.png` (gitignored) for the "Draw a
picture" section of docs/api/system.md, so the docs build must run before mkdocs.

2026-11-23 (0.5.0, branch claude/triclinic): tilted boxes, strain, the pressure tensor, and
editable states.
- Tilted (triclinic) boxes as in LAMMPS: `state.tilt` = (xy, xz, yz), edges a = (lx,0,0),
  b = (xy,ly,0), c = (xz,yz,lz); each tilt at most half the length it leans along;
  periodic 3D unbonded atoms only (walls, flat systems, molecules and the molecular
  engine refuse a tilt). `system.wrap` and C `tilted_wrap` wrap into the parallelepiped
  (fractional 0 <= s < 1). C `tilted_image` takes whole edges off z, y, x in turn; the
  result is the unique image inside the rectangle |d_i| <= l_i/2, which is the
  nearest image whenever it is shorter than half the smallest face distance (checked:
  face distances > 2*cutoff). A NULL/zero tilt keeps every per-axis code path, so
  rectangular trajectories are bit-identical (golden tests) and as fast as before
  (2048-atom liquid 1389 steps/s); tilted runs ~30% slower (972). Cells are laid out
  in fractional coordinates with counts from face distances; `needs_rebuild` maps the
  build positions by M = H H0^-1; the isotropic barostat scales the tilt; the tilt is
  part of the workspace comparison (`Workspace.tilt`, `Neighbors.tilt`). `_core.advance`
  takes optional trailing `tilt` (writable (3,)) and `stress` ((6,)) after `record`;
  `neighbor_pairs` takes an optional tilt.
- `state.strain(xx, yy, zz, yz, xz, xy)`: F = [[1+xx, xy, xz], [0, 1+yy, yz], [0, 0, 1+zz]],
  H' = F H, positions and unwrapped x' = F x (LAMMPS change_box ... remap), in place.
- Pressure tensor: `stress_numerator` (one pair pass, only when asked, so the force
  loop is unchanged) gives V P_ab = sum m v_a v_b + sum r_a f_b; observables `pxx, pyy,
  pzz, pxy, pxz, pyz` (pxx, pyy, pxy when flat) for unbonded atoms without walls, read
  by a zero-step call at each thermo row. Checked: the strain derivative of the energy
  (all six components) and LAMMPS `scripts/validate_triclinic.py` (energy 1e-15,
  forces 1e-11, tensor 1e-14 with Homebrew lmp_serial lj/cut).
- Frames and checkpoints store `tilt` when nonzero; `run.tilts`; `run.state` restores it;
  draw() and viewer.js (`cornerPoint`, `centre`, `cornersAbout`, `span`) draw the
  parallelepiped. Not yet: tilt in the molecular engine and PPPM (Part C).
- Editable states (the user wants everything editable, to see what happens):
  positions, velocities (None = random at the next run), masses, species labels
  (stored '<U32'), box, tilt, boundary (`_rebuild` through the constructor), time and
  step; `state.atoms` arrays are the state's own. The next run checks everything (it
  copies the state through the constructor); a run never changes its input. A nonzero
  total momentum is refused with the fix in the message. LiveRun getters stay
  read-only snapshots (so `before = live.positions` keeps its values) and take
  assignment for positions (unwrapped follow), velocities, masses, box and tilt.
- `examples/elastic.ipynb` (built by build_docs.py, page docs/examples/elastic.md): LJ FCC,
  N=500, NPT relax at P=0, T=0.1 (rho 1.059), then NVT strains of +-0.01 in xx and xy
  (20 000 steps each) give C11 81.5, C12 47.8, C44 49.5 (reduced); VRH moduli -> v_l 9.8,
  v_t 5.5, omega_D 24.4, theta_D 86 K in argon units (hbar* = 0.0296; argon ~92 K).
  Einstein omega_E 14.1 (theta_E 50 K) from <u^2> of `run.unwrapped`. The velocity power
  spectrum (4096 samples via `live.advance(2)`) is drawn against Debye omega^2 and the
  Einstein line; C_V(T) from it lies between Einstein and Debye and drops below Debye
  under a few K (finite box: no long waves). About 25 s to execute. Part C (alpha-quartz:
  metal units, BKS, triclinic PPPM, anisotropic barostat) is planned, not started.

2026-11-23 (0.6.0, branch claude/alpha-quartz): Part C, alpha-quartz. Public surface:
`AlphaQuartz(cells=(4, 2, 3), a=4.9134, c=5.4052)` (metal units, species Si/O, a
rectangular a x sqrt(3)a x c cell of 18 atoms, P3_2 21 positions of Levien et al. 1980,
checked: 4 O per Si at 1.604/1.613 Å, Si-O-Si 143.7 deg, 2.6487 g/cm^3);
`Buckingham(A, rho, C)` with `.energy(r)` and `.barrier` (inner maximum, 1.098 Å for
BKS Si-O); `MDSimulation(pair_potential={("Si","O"): Buckingham(...), ...})` (pairs left
out do not interact at short range; `sim.pair_potential == "buckingham"`, `sim.pairs`;
config `pairs`); `units="metal"` (LAMMPS constants; defaults skin 2 Å, friction 5/ps,
pressure_time 1 ps); `barostat="isotropic"|"anisotropic"` (constructor, run, start,
set; settings now always contain `barostat`); `tail=True` (needs lj_shift=False, as
LAMMPS refuses shift+tail) for the charged-atom solver, LJ and Buckingham. Observables
for atoms in the charged-atom solver: `pxx..pyz` from every step's tensor,
`buckingham_energy` (instead of lj_energy), `mass_density` in metal units.

How it is built. `tilt.h` holds the tilted-box helpers (moved out of _core.c) for
_core.c, molecular.c and pppm.c. molecular.c: Model.style/buckingham (types x types x 4:
A, rho, C, barrier; the barrier check raises a message naming the atoms), Model.tilt
(NULL when zero, so rectangular code paths are unchanged), `pair_shift` and
`build_list` use tilted_image, drift and scale_molecules wrap into the parallelepiped
(atoms only; tilt scales by row, F = diag(f): xy, xz with x, yz with y), face-distance
checks, `Terms.tensor` (pair + Coulomb + k-space + tail, atoms only), result "tensor"
= kinetic + virial tensor. Anisotropic MTK: `rates` (3,) buffer, inertia baromass/3
per axis, so equal rates reproduce the isotropic barostat; conserved
extended_enthalpy (checked in NPH). Tail: E V and W V precomputed from type counts,
divided by the current volume each force call. Workspace: ws->tilt compared like the
box; Buckingham tables and the style/tail flags are in the inputs (scalar slots 16, 17).
pppm.c: `pppm_select(width, lengths, ...)` (widths = face distances choose the mesh;
lengths keep the old volume arithmetic, which matters bitwise), `update_greens_tilted`
(non-separable Gaussian over aliases of k = 2 pi H^-T w; recomputed when box or tilt
change), charges spread in fractional coordinates, the field uses k with Nyquist
components zeroed (equivalent to LAMMPS taking the real part), k-space virial tensor.
The PPPM state stores widths as the reference for re-selection. Water NVT/NPT and
salt trajectories are bitwise identical to main (checked against a main worktree).

Validation: `scripts/validate_quartz.py PATH_TO_LMP` (serial LAMMPS with KSPACE built
from hsp_hymcp/lammps with `make yes-kspace; make serial` in a scratch copy; Homebrew's
segfaults in PPPM): 1080 atoms, all tilts, energy 4e-15, forces 1.6e-5 eV/Å, off-diagonal
pressure < 0.1 bar, diagonal 3 bar, fully explained by LAMMPS's real-space Coulomb force
using the exact erfc derivative with the approximate erfc (ours is the derivative of
our energy; the predicted offsets 3.17/3.07/3.04 bar match). Buckingham alone: forces
1e-12, pressure 1e-8 bar. The PPPM virial differs from the strain derivative by a few
bar in both codes. Tests: tests/test_quartz.py (structure, terms, barrier, tail against
numerical integrals for Buckingham and LJ, tilted PPPM vs a tilted direct Ewald sum
(tests/ewald_reference.py takes tilt and tensor), tensor vs strain derivative, tilted
chunks/chains, anisotropic relaxation, NPH conservation, rates through checkpoints).
Speed: 432 atoms, cutoff 8, 1.6 ms/step at kspace_accuracy 1e-4 (3.0 at 1e-5).

BKS at 300 K, 0 bar (tail on, cutoff 8, 432 atoms): a 4.935, c 5.442 Å, 2.608 g/cm^3.
Without the tail the -C/r^6 beyond 8 Å (about -1.1 GPa) is missing and the density is
2.53. Elastic constants (GPa, 10 ps per +-1% strain): C11 88.5, C12 7.6, C13 15.3,
|C14| 16.8-17.1, C33 106.7, C44 49.1, C66 40.4, against Ogi et al. 2006 (J. Appl. Phys.
100, 053511, Table I averages at 30 C): 86.76, 6.868, 11.85, -18.02, 105.46, 58.14,
39.95. Our C14 is positive in this setting; its sign depends on handedness and axes.
`examples/quartz.ipynb` (built by build_docs.py; generator kept out of the repo, edit
the notebook): NPT anisotropic relax 20 ps (a 4.935, c 5.441 Å, 2.609 g/cm^3), eight
NVT strains of 10 ps (C11 88.3, C12 7.4, C13 15.1, C14 16.9, C33 106.8, C44 49.1, C66 40.4
GPa; symmetry partners agree to 0.2 GPa), VRH from the full 6x6 matrix: K 39.6, G 41.3
GPa, v_l 6020, v_t 3976 m/s, theta_D 555 K (572 K from Ogi's constants). The
mass-weighted velocity spectrum has bands at 15-18, 22-27, 33 and 37 THz; C_V/3Nk at
300 K: spectrum 0.54, Debye 0.85, measured C_p 44.6 J/(mol K) -> 0.60. About 3 minutes.
`state.draw()` and `run.draw()` use element display sizes (O 1.2, Si 0.9 Å) in metal
units, as the player does.

2026-11-24 (0.6.1): `run` flushed thermo.csv after every row, which made `thermo_every=1`
slow on network and synced folders (a JupyterHub home, Google Drive, OneDrive). It now
flushes the header and first row at once, then at most once a second, and the file is
complete when the run ends. Locally a flush costs about 3 µs, so the local gain is small
(20 000 steps with a row per step: 1.44 -> 1.36 s); the rest of a row's cost (about 35 µs)
is writing the CSV row and reading the observables. `test_thermo_is_on_disk_during_the_run`.
