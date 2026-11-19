"""A flat box of atoms between four walls: temperature from the velocities, force on the right wall.

Scene for the sims `flat-box`, `wall-force` and `wide-box` of the session
"Energy, heat and temperature" (week 34). Plain use of the public API of
fys2160_md (0.3.0 or newer: two dimensions and walls). Reduced units: length
sigma, energy epsilon, mass m, k = 1.

Two ways in. The page runs `run_box` and `wall_force`: generators that advance
a run in memory and hand it over once per displayed frame. `saved_run` does the
same physics as a saved run, for your own Python; run this file to print the
table of the wall-force task.
"""
import numpy as np
from fys2160_md import MDSimulation, Random

N = 200                       # atoms
DT = 0.005                    # timestep; 200 steps per time unit
SETTLE = 100.0                # time with the thermostat before anything is shown or measured
# A light thermostat: an atom keeps its velocity for about 1/FRICTION = 20 time
# units, twice the time between two collisions, so it flies straight between
# collisions as in a gas. (The package default, 1, would turn the motion into a
# random walk, and the force on a wall would then take far longer to average.)
FRICTION = 0.05


def simulation(T=2.0):
    return MDSimulation(pair_potential="lj", ensemble="nvt", temperature=T,
                        timestep=DT, cutoff=2.5, friction=FRICTION)


def box(Lx=100.0, Ly=100.0, seed=34):
    """N atoms at random positions in a flat box with a wall on every side."""
    return Random(N=N, box=(Lx, Ly), boundary="walls", min_distance=1.0, seed=seed)


def settled(T=2.0, Lx=100.0, Ly=100.0, seed=34):
    """The gas after SETTLE time units with the thermostat at T, as a state.
    The seed sets the starting positions, the starting velocities and the
    random pushes of the thermostat, also in the run that continues from it."""
    live = simulation(T).start(box(Lx, Ly, seed), seed=seed)
    live.advance(round(SETTLE / DT))
    return live.state


def paced(live, until, slow, fast, ramp_from, ramp_to):
    """Advance `live` until `until` time units after its start, one displayed frame per yield.

    A frame covers `slow` steps until `ramp_from`, so the motion can be
    followed. From there the run speeds up smoothly, reaching `fast` steps a
    frame at `ramp_to`, and keeps that pace: a fast-forward, to collect
    statistics. Yields how many times faster than the slow pace the frame was.
    How the steps are divided does not change the trajectory."""
    per_time, first = round(1 / DT), live.step
    while live.step - first < round(until * per_time):
        t = (live.step - first) * DT
        if t < ramp_from - DT / 2:
            steps, stop = slow, ramp_from
        elif t < ramp_to - DT / 2:
            steps, stop = round(slow * (fast / slow) ** ((t - ramp_from) / (ramp_to - ramp_from))), ramp_to
        else:
            steps, stop = fast, until
        steps = max(1, min(steps, round(stop * per_time) - (live.step - first)))
        live.advance(steps)
        yield steps / slow


def run_box(T=2.0, Lx=100.0, Ly=100.0, until=60.0):
    """The page's `flat-box`: the atoms, and nothing else. Every yield is one frame.
    The thermostat holds T; the gas has settled before the first frame."""
    live = simulation(T).start(settled(T, Lx, Ly))
    first = live.step
    yield live, {"t": 0.0}
    while live.step - first < round(until / DT):
        live.advance(10)
        yield live, {"t": (live.step - first) * DT}


def wall_force(T=2.0, Lx=100.0, Ly=100.0, until=2000.0):
    """The page's `wall-force` and `wide-box`: the force on the right wall.

    "force" is the mean force over the frame, the momentum the wall received
    during the frame divided by the frame's duration: single kicks while the
    frames are short. "mean_force" is the mean since the start, the read-off.
    Slow for the first 10 time units, then a fast-forward that reaches 50 times
    the pace at t = 30."""
    live = simulation(T).start(settled(T, Lx, Ly))
    first, start = live.step, live.observables["wall_impulse_right"]
    before, t_before = start, 0.0
    yield live, {"t": 0.0, "force": 0.0, "mean_force": 0.0, "speedup": 1}
    for speedup in paced(live, until, slow=20, fast=1000, ramp_from=10.0, ramp_to=30.0):
        t = (live.step - first) * DT
        impulse = live.observables["wall_impulse_right"]
        yield live, {"t": t, "force": (impulse - before) / (t - t_before),
                     "mean_force": (impulse - start) / t, "speedup": speedup}
        before, t_before = impulse, t


def saved_run(T=2.0, Lx=100.0, Ly=100.0, time=2000.0, name="box"):
    """The same measurement as a saved run: the thermo table has the momentum each
    wall has received (wall_impulse_right, ...) at every row."""
    return simulation(T).run(settled(T, Lx, Ly), steps=round(time / DT), thermo_every=200,
                             trajectory_every=None, storage_name=name)


if __name__ == "__main__":
    # The table of the wall-force task, and the wider box of the task after it.
    for T, Lx, Ly in ((2, 100, 100), (4, 100, 100), (2, 200, 100), (2, 100, 200)):
        thermo = saved_run(T, Lx, Ly).thermo
        time, impulse = np.asarray(thermo.time), np.asarray(thermo.wall_impulse_right)
        F = (impulse[-1] - impulse[0]) / (time[-1] - time[0])
        T_measured = thermo.temperature.mean()
        print(f"T* = {T}, box {Lx} x {Ly}: F = {F:.3f}, P = F / Ly = {F / Ly:.4f}, "
              f"PV/NkT = {F * Lx / (N * T_measured):.3f}")
