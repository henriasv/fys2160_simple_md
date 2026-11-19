"""The numbers the week 34 page shows, from its two scenes, with the time they take.

    python examples/week34/page_numbers.py

The same file runs in the browser build: `node scripts/pyodide_week34.mjs WHEEL`
runs it in Pyodide under Node. Prints one JSON object.
"""
import hashlib
import json
import sys
import tempfile
import time
import numpy as np
import fys2160_md
from fys2160_md import MDSimulation, Random

import box
import piston


def digest(*arrays):
    """A short fingerprint of arrays, to see whether two builds agree bit for bit."""
    return hashlib.sha256(b''.join(np.ascontiguousarray(a).tobytes() for a in arrays)).hexdigest()[:12]


def wall_force_runs():
    """The four runs of the wall-force and wide-box tasks, as the page reads them off."""
    out = {}
    for T, Lx, Ly in ((2.0, 100.0, 100.0), (4.0, 100.0, 100.0), (2.0, 200.0, 100.0), (2.0, 100.0, 200.0)):
        start, frames, temperatures = time.perf_counter(), 0, []
        for live, extra in box.wall_force(T, Lx, Ly):
            frames += 1
            temperatures.append(live.observables["temperature"])
        seconds = time.perf_counter() - start
        F, T_mean = extra["mean_force"], float(np.mean(temperatures[len(temperatures)//2:]))
        out[f"T={T:g} {Lx:g}x{Ly:g}"] = dict(F=round(F, 3), P=round(F / Ly, 4), PV_NkT=round(F * Lx / (box.N * T), 3),
                                             T_late=round(T_mean, 3), frames=frames, steps=live.step,
                                             seconds=round(seconds, 2), fingerprint=digest(live.positions, live.velocities))
    return out


def piston_run():
    start, frames, rows = time.perf_counter(), 0, []
    for live, extra in piston.compress():
        frames += 1
        rows.append(extra)
    seconds = time.perf_counter() - start
    first, end = rows[0], rows[1000]                    # one frame per time unit; the stroke ends at t = 1000
    trace = {str(round(rows[k]["V"])): round(rows[k]["P"], 4) for k in (100, 500, 900, 1000)}
    return dict(T_before=round(first["T"], 3), T_after=round(end["T"], 3), factor=round(end["T"] / first["T"], 3),
                W=round(end["W"], 1), energy_change=round(end["E"] - first["E"], 1),
                NkdT=round(box.N * (end["T"] - first["T"]), 1), T_rest=round(rows[-1]["T"], 3), P_at_V=trace,
                frames=frames, steps=live.step, seconds=round(seconds, 2), fingerprint=digest(live.positions, live.velocities))


def code_blocks():
    """The two code blocks of the session, as written, read without pandas."""
    with tempfile.TemporaryDirectory() as folder:
        gas = Random(N=200, box=(100, 100), boundary="walls", min_distance=1.0)
        sim = MDSimulation(pair_potential="lj", ensemble="nvt", temperature=2.0,
                           timestep=0.005, cutoff=2.5, output_dir=folder)
        run = sim.run(gas, steps=4000, thermo_every=100, trajectory_every=None, storage_name="flat-box")
        v = run.state.atoms.velocities
        K = 0.5 * np.sum(v**2)
        thermo = np.genfromtxt(run.path/"thermo.csv", delimiter=",", names=True)
        flat = dict(K_per_atom=round(float(K / len(v)), 3), T_engine=round(float(thermo["temperature"][-1]), 3),
                    fingerprint=digest(v))
        gas = Random(N=500, rho=0.01, min_distance=0.9)
        run = sim.run(gas, steps=2000, thermo_every=50, trajectory_every=None, storage_name="install-check")
        thermo = np.genfromtxt(run.path/"thermo.csv", delimiter=",", names=True)
        late = thermo["step"] >= 1000
        ten = dict(temperature=round(float(thermo["temperature"][late].mean()), 3),
                   kinetic_energy_per_atom=round(float(thermo["kinetic_energy"][late].mean() / gas.N), 3),
                   PV_NkT=round(float(thermo["compressibility_factor"][late].mean()), 3),
                   fingerprint=digest(run.state.atoms.positions))
    return dict(flat_box=flat, ten_lines=ten)


if __name__ == "__main__":
    print(json.dumps(dict(engine=fys2160_md.__version__, platform=sys.platform, wall_force=wall_force_runs(),
                          piston=piston_run(), code_blocks=code_blocks()), indent=1))
