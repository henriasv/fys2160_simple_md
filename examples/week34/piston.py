"""Compress the flat box with a piston, with the thermostat off.

Scene for the sim `piston` of the session "Energy, heat and temperature"
(week 34). Plain use of the public API of fys2160_md (0.3.0 or newer).

The page runs `compress`, a generator that hands the run over once per
displayed frame. `saved_run` is the same stroke as a saved run; run this file
to print the numbers of the two piston tasks.
"""
from collections import deque
import numpy as np
from fys2160_md import MDSimulation, Piston

from box import N, DT, settled

sim = MDSimulation(pair_potential="lj", ensemble="nve", timestep=DT, cutoff=2.5)


def compress(speed=0.05, x_end=50.0, rest=50.0, smooth=50):
    """The page's `piston`: the gas has settled at T* = 2 with the thermostat.
    Then the thermostat is off and the right wall moves from x = 100 to x_end
    at constant speed, stops, and the run goes on for `rest` time units.

    One frame per time unit. "T" is the temperature, "x" the position of the
    piston, "V" the volume (the area of the box), "W" the work the piston has
    done on the gas. "P" is the force on the piston per length of the piston,
    as a mean over the last `smooth` frames: the piston takes only a few hits
    per frame."""
    live = sim.start(settled(2.0), piston=Piston("right", to=x_end, speed=speed))
    first, Ly = live.step, float(live.box[1])
    stroke = (float(live.box[0]) - x_end) / speed
    recent = deque([(0.0, live.observables["wall_impulse_right"])], maxlen=smooth + 1)

    def columns():
        obs = live.observables
        t = (live.step - first) * DT
        if t > recent[-1][0]:
            recent.append((t, obs["wall_impulse_right"]))
        P = (recent[-1][1] - recent[0][1]) / (recent[-1][0] - recent[0][0]) / Ly if len(recent) > 1 else 0.0
        return {"t": t, "T": obs["temperature"], "x": obs["piston_position"], "V": obs["volume"],
                "P": P, "W": obs["work"], "E": obs["total_energy"]}

    yield live, columns()
    while live.step - first < round((stroke + rest) / DT):
        live.advance(round(1 / DT))
        yield live, columns()


def saved_run(speed=0.05, x_end=50.0, name="piston"):
    """The stroke as a saved run: the thermo table has temperature, volume,
    total_energy, work, piston_position and wall_impulse_right at every row."""
    start = settled(2.0)
    steps = round((start.box.lengths[0] - x_end) / speed / DT)
    return sim.run(start, steps=steps, piston=Piston("right", to=x_end, speed=speed),
                   thermo_every=200, trajectory_every=None, storage_name=name)


if __name__ == "__main__":
    thermo = saved_run().thermo
    T, E, W = np.asarray(thermo.temperature), np.asarray(thermo.total_energy), np.asarray(thermo.work)
    print(f"T: {T[0]:.3f} -> {T[-1]:.3f}   (2 T0 = {2 * T[0]:.3f})")
    print(f"change of total energy {E[-1] - E[0]:.1f};  N k (T1 - T0) = {N * (T[-1] - T[0]):.1f};  "
          f"work on the gas {W[-1]:.1f}")
    # The trace of the page: the mean force on the piston per length over a tenth of the stroke, against the volume.
    time, impulse, V = np.asarray(thermo.time), np.asarray(thermo.wall_impulse_right), np.asarray(thermo.volume)
    tenth = (len(time) - 1) // 10
    for a in range(0, len(time) - tenth, tenth):
        b = a + tenth
        print(f"   V = {V[a:b + 1].mean():6.0f}: P = {(impulse[b] - impulse[a]) / (time[b] - time[a]) / 100:.4f}")
