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
