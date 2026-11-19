"""Validation of the flat box, the walls and the piston (slow: about five minutes).

Prints the numbers quoted in docs/physics.md and in the teacher notes of the
week 34 session: energy conservation between walls, the force on a wall
against N k T / L, and the piston's work against the change in energy.

    python scripts/validate_walls.py          # 8 seeds per case
    python scripts/validate_walls.py 2        # 2 seeds per case, a quick look
"""
import sys
import tempfile
from pathlib import Path
import numpy as np
import fys2160_md
from fys2160_md import MDSimulation, Piston, Random

sys.path.insert(0, str(Path(__file__).resolve().parent.parent/'examples'/'week34'))
import box as scene

N, DT = scene.N, scene.DT
WALLS = ('left', 'right', 'bottom', 'top')
SEEDS = range(34, 34+(int(sys.argv[1]) if len(sys.argv) > 1 else 8))


def second_virial(T, cutoff=2.5):
    """B2 of the flat Lennard-Jones gas (cut and shifted): -pi * integral of (exp(-u/T) - 1) r dr."""
    r = np.linspace(1e-4, cutoff, 20001)
    u = 4*(r**-12-r**-6)-4*(cutoff**-12-cutoff**-6)
    f = (np.exp(-u/T)-1)*r
    return float(-np.pi*np.sum(.5*(f[1:]+f[:-1])*np.diff(r)))


def adiabat_ratio(T0=2., rho0=.02, rho1=.04, steps=200):
    """T1/T0 along the adiabat of the flat gas with its second virial coefficient:
    P = rho T (1 + B2 rho) and U/N = T - T^2 rho dB2/dT, with dU = -P dV."""
    h = 1e-3
    d1 = lambda T: (second_virial(T+h)-second_virial(T-h))/(2*h)
    d2 = lambda T: (second_virial(T+h)-2*second_virial(T)+second_virial(T-h))/h**2
    slope = lambda rho, T: T/rho*(1+rho*(second_virial(T)+T*d1(T)))/(1-rho*(2*T*d1(T)+T*T*d2(T)))
    T, step = T0, (rho1-rho0)/steps
    for k in range(steps):                                      # midpoint rule
        rho = rho0+k*step
        T += step*slope(rho+.5*step, T+.5*step*slope(rho, T))
    return T/T0


def spread(values):
    values = np.asarray(values, dtype=float)
    return f'{values.mean():.4g} (sd {values.std(ddof=1):.2g} over {len(values)} seeds)' if len(values) > 1 else f'{values[0]:.4g}'


def energy_conservation():
    print('1. Energy, thermostat off, 200 atoms, 1000 time units (200 000 steps)')
    for boundary in ('walls', 'periodic'):
        for T in (2., 4.):
            drifts, kinetic = [], []
            for seed in SEEDS:
                gas = Random(N=N, box=(100, 100), boundary=boundary, min_distance=1.0, seed=seed)
                warm = MDSimulation(pair_potential='lj', ensemble='nvt', temperature=T, timestep=DT, cutoff=2.5, seed=seed).start(gas)
                warm.advance(20000)
                live = MDSimulation(pair_potential='lj', ensemble='nve', timestep=DT, cutoff=2.5).start(warm.state)
                first = live.observables
                live.advance(200000)
                drifts.append(live.observables['total_energy']-first['total_energy'])
                kinetic.append(first['kinetic_energy'])
            rms = float(np.sqrt(np.mean(np.square(drifts))))
            print(f'   {boundary:8s} T* = {T:g}: rms change of total energy {rms:.3f}, {rms/np.mean(kinetic):.1e} of the kinetic energy')


def wall_force(T, Lx, Ly, friction=scene.FRICTION, time=2000.):
    """Mean force on each wall over `time`, the mean temperature and the block scatter, for every seed."""
    rows = []
    for seed in SEEDS:
        sim = MDSimulation(pair_potential='lj', ensemble='nvt', temperature=T, timestep=DT, cutoff=2.5, friction=friction, seed=seed)
        live = sim.start(Random(N=N, box=(Lx, Ly), boundary='walls', min_distance=1.0, seed=seed))
        live.advance(round(scene.SETTLE/DT))
        impulses, temperatures = [], []
        for _ in range(round(time)):                       # one time unit at a time
            if not impulses:
                impulses.append([live.observables[f'wall_impulse_{w}'] for w in WALLS])
            live.advance(200)
            obs = live.observables
            impulses.append([obs[f'wall_impulse_{w}'] for w in WALLS])
            temperatures.append(obs['temperature'])
        impulses = np.array(impulses)
        force = (impulses[-1]-impulses[0])/time
        blocks = np.diff(impulses[::100, 1])/100.          # right wall, means over 100 time units
        rows.append(dict(F=force, T=np.mean(temperatures), block_sd=blocks.std(ddof=1)))
    return rows


def wall_forces():
    print('2. Force on a wall, thermostat on (friction %g), means over 2000 time units after %g' % (scene.FRICTION, scene.SETTLE))
    print('   second virial coefficient of the flat gas: B2(T*=2) = %.3f, B2(T*=4) = %.3f' % (second_virial(2.), second_virial(4.)))
    reference = None
    for run, (T, Lx, Ly) in enumerate(((2., 100., 100.), (4., 100., 100.), (2., 200., 100.), (2., 100., 200.)), start=1):
        rows = wall_force(T, Lx, Ly)
        right = np.array([row['F'][1] for row in rows])
        left = np.array([row['F'][0] for row in rows])
        Tm = np.array([row['T'] for row in rows])
        z_right = right*Lx/(N*Tm)
        # All four walls: the pressure as force per length of wall, times the area, over N k T.
        lengths = np.array([Ly, Ly, Lx, Lx])
        z_all = np.array([(row['F']/lengths).mean()*Lx*Ly/(N*row['T']) for row in rows])
        expected = 1+second_virial(T)*N/(Lx*Ly)
        print(f'   run {run}: T* = {T:g}, box {Lx:g} x {Ly:g}; ideal F = {N*T/Lx:.2f}, expected PV/NkT = {expected:.4f}')
        print(f'      F on the right wall {spread(right)}; P = F/Ly {right.mean()/Ly:.4f}; measured T {Tm.mean():.3f}')
        print(f'      F on the left wall {spread(left)}')
        print(f'      PV/NkT from the right wall {spread(z_right)}; from all four walls {spread(z_all)}')
        print(f'      scene seed {SEEDS[0]}: F = {right[0]:.3f}, P = {right[0]/Ly:.4f}, PV/NkT = {z_right[0]:.3f}; '
              f'sd of means over 100 time units {np.mean([row["block_sd"] for row in rows]):.2f}')
        if reference is None:
            reference = right.mean()
        else:
            print(f'      F / F(run 1) = {right.mean()/reference:.3f}')
    rows = wall_force(2., 100., 100., friction=1.)
    right = np.array([row['F'][1] for row in rows])
    print(f'   run 1 with the default friction 1: F on the right wall {spread(right)}; '
          f'sd of means over 100 time units {np.mean([row["block_sd"] for row in rows]):.2f}')


def stroke(seed, speed, epsilon=1., sigma=1., cutoff=2.5, x_end=50., count=False):
    """One compression from a state settled at T* = 2; returns the numbers of the piston tasks."""
    start = scene.settled(2., seed=seed)
    sim = MDSimulation(pair_potential='lj', epsilon=epsilon, sigma=sigma, ensemble='nve', timestep=DT, cutoff=cutoff)
    live = sim.start(start, piston=Piston('right', to=x_end, speed=speed))
    first = live.observables
    steps = round((100.-x_end)/speed/DT)
    bounces = 0
    if count:                                                  # a bounce: the force on the piston leaves zero
        touching = first['wall_force_right'] > 0
        for _ in range(steps):
            live.advance(1)
            now = live.observables['wall_force_right'] > 0
            bounces += now and not touching
            touching = now
    else:
        live.advance(steps)
    last = live.observables
    return dict(T0=first['temperature'], T1=last['temperature'], dE=last['total_energy']-first['total_energy'],
                W=last['work'], NkdT=N*(last['temperature']-first['temperature']),
                dU=last['potential_energy']-first['potential_energy'], bounces=bounces, x=last['piston_position'])


def pistons():
    print('3. Piston, thermostat off, right wall from x = 100 to 50, 200 atoms settled at T* = 2')
    print(f'   slow compression of the flat gas with its second virial coefficient, from T* = 2: T1/T0 = {adiabat_ratio():.3f}')
    for speed in (.05, .1, .2):
        rows = [stroke(seed, speed) for seed in SEEDS]
        get = lambda key: np.array([row[key] for row in rows])
        print(f'   speed {speed:g} ({50/speed:g} time units):')
        print(f'      T before {spread(get("T0"))}, after {spread(get("T1"))}, ratio {spread(get("T1")/get("T0"))}  (flat ideal gas, T V = const: 2)')
        print(f'      work {spread(get("W"))}; change of total energy - work: largest {np.abs(get("dE")-get("W")).max():.3f}, '
              f'{np.abs((get("dE")-get("W"))/get("W")).max():.1e} of the work')
        print(f'      N k (T1 - T0) {spread(get("NkdT"))}; change of potential energy {spread(get("dU"))}')
        row = rows[0]
        print(f'      scene seed {SEEDS[0]}: T {row["T0"]:.3f} -> {row["T1"]:.3f} (factor {row["T1"]/row["T0"]:.3f}), '
              f'N k dT = {row["NkdT"]:.1f}, change of energy {row["dE"]:.1f}, W = {row["W"]:.1f}')
    rows = [stroke(seed, .05, epsilon=1e-12, sigma=1e-9, cutoff=.01) for seed in SEEDS]
    ratio = np.array([row['T1']/row['T0'] for row in rows])
    print(f'   atoms that do not collide (epsilon 1e-12, sigma 1e-9, cutoff 0.01), speed 0.05: ratio {spread(ratio)}  (expected 2.5: only the motion along x heats, fourfold)')
    row = stroke(SEEDS[0], .05, count=True)
    print(f'   bounces off the piston during the stroke at speed 0.05 (scene seed): {row["bounces"]}')
    live = scene.simulation(2.).start(scene.settled(2., seed=SEEDS[0]))
    touching, bounces = False, 0
    for _ in range(100000):
        live.advance(1)
        now = live.observables['wall_force_right'] > 0
        bounces += now and not touching
        touching = now
    print(f'   bounces off the right wall at rest, T* = 2, 100 x 100: {bounces/500:.2f} per time unit '
          f'(ideal: N <|vx|> / (2 Lx) = {N*np.sqrt(2*2/np.pi)/200:.2f})')


def code_blocks():
    print('4. The code blocks of the session, as written (default thermostat friction)')
    with tempfile.TemporaryDirectory() as folder:
        gas = Random(N=200, box=(100, 100), boundary="walls", min_distance=1.0)
        sim = MDSimulation(pair_potential="lj", ensemble="nvt", temperature=2.0,
                           timestep=0.005, cutoff=2.5, output_dir=folder)
        run = sim.run(gas, steps=4000, thermo_every=100, trajectory_every=None, storage_name="flat-box")
        v = run.state.atoms.velocities
        K = 0.5*np.sum(v**2)
        print(f'   flat box: kinetic energy per atom {K/len(v):.3f}; T = K/N {K/len(v):.3f}; '
              f'temperature from the engine {float(np.asarray(run.thermo.temperature)[-1]):.3f}')
        T = np.asarray(run.thermo.temperature)[10:]
        print(f'      thermo rows after step 1000: mean T {T.mean():.3f}, sd of one row {T.std(ddof=1):.3f} (T/sqrt(N) = {2/np.sqrt(200):.3f})')
        gas = Random(N=500, rho=0.01, min_distance=0.9)
        run = sim.run(gas, steps=2000, thermo_every=50, trajectory_every=None, storage_name="install-check")
        late = run.thermo.step >= 1000
        print(f'   the ten lines: temperature {run.thermo.temperature[late].mean():.3f}, kinetic energy per atom '
              f'{run.thermo.kinetic_energy[late].mean()/gas.N:.3f}, PV/NkT {run.thermo.compressibility_factor[late].mean():.3f}')


if __name__ == '__main__':
    print(f'fys2160_md {fys2160_md.__version__}, {len(SEEDS)} seeds from {SEEDS[0]} (the scenes use {SEEDS[0]})')
    energy_conservation()
    wall_forces()
    pistons()
    code_blocks()
