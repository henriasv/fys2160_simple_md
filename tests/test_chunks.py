"""A run advanced in memory by its caller: chunks, continuation, controls and velocities."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from fys2160_md import FCC, Random, Plummer, IceIh, MDSimulation, LiveRun, RigidBond, RigidWater, Run
from fys2160_md import _core

LJ = dict(pair_potential='lj', epsilon=1., sigma=1., timestep=.005, cutoff=2.5)


def advance(live, total, size):
    """Take total steps in calls of at most size steps."""
    while total:
        live.advance(min(size, total))
        total -= min(size, total)


def pairs(x, box, radius, reach):
    """Listed pairs from checking all pairs (reach 0) or from cells that span radius/reach."""
    return np.frombuffer(_core.neighbor_pairs(np.ascontiguousarray(x), tuple(box), radius, reach), dtype=np.intc).reshape(-1, 2)


class ChunkTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)/'runs'

    def cases(self):
        water = RigidWater(oh_length=1.0, hoh_angle=109.47)
        yield 'gas, cells', FCC(N=500, rho=.1), dict(ensemble='nvt', temperature=2., **LJ)
        yield 'liquid, finer cells', FCC(N=864, rho=.8), dict(ensemble='nve', temperature=1., **LJ)
        yield 'small box, all pairs', FCC(N=108, rho=.8), dict(ensemble='nvt', temperature=1., **LJ)
        yield 'barostat', FCC(N=256, rho=.3), dict(ensemble='npt', temperature=2., pressure=.5, **LJ)
        yield 'heated', FCC(N=256, rho=.05), dict(ensemble='nve', temperature=1., heat_rate=3., **LJ)
        yield 'rigid molecules', Random(N=100, rho=.02, molecule='diatomic', bond=RigidBond(length=.7)), \
            dict(ensemble='nvt', temperature=1., **{**LJ, 'timestep': .002})
        yield 'gravity', Plummer(N=64, G=1/64, softening=.05), \
            dict(pair_potential='gravity', G=1/64, softening=.05, ensemble='nve', timestep=.005, cutoff=None)
        yield 'water', IceIh(N=288, density=.93, bond=water), \
            dict(units='real', epsilon={'O': .15535, 'H': 0.}, sigma={'O': 3.166, 'H': 0.}, charge={'O': -.8476, 'H': .4238},
                 coulomb='pppm', kspace_accuracy=1e-4, cutoff=6., temperature=250., timestep=1., ensemble='nvt')

    def test_chunk_size_does_not_change_the_trajectory(self):
        for name, state, settings in self.cases():
            steps = 21 if name == 'water' else 500
            with self.subTest(name):
                sim = MDSimulation(output_dir=self.root, **settings)
                whole = sim.start(state)
                whole.advance(steps)
                for size in (1, 7) if name == 'water' else (1, 7, 500):
                    live = sim.start(state)
                    advance(live, steps, size)
                    np.testing.assert_array_equal(live.positions, whole.positions)
                    np.testing.assert_array_equal(live.velocities, whole.velocities)
                    np.testing.assert_array_equal(live.box, whole.box)
                    self.assertEqual(live.step, steps)
                    a, b = live.observables, whole.observables
                    self.assertEqual(list(a), list(b))
                    for key in a:
                        # The clock and the heat are summed per call; everything else is exact.
                        if key in ('time', 'heat_added'):
                            self.assertAlmostEqual(a[key], b[key], delta=1e-12*max(1., abs(b[key])), msg=key)
                        else:
                            np.testing.assert_array_equal(a[key], b[key], err_msg=key)
                # sim.run is the same loop with files.
                run = sim.run(state, steps=steps, thermo_every=13, trajectory_every=None, storage_name='saved')
                np.testing.assert_array_equal(run.state.atoms.positions, whole.positions)
                np.testing.assert_array_equal(run.state.atoms.velocities, whole.velocities)

    def test_start_writes_nothing_and_leaves_the_state_alone(self):
        state = FCC(N=108, rho=.1)
        sim = MDSimulation(output_dir=self.root, ensemble='nvt', temperature=2., **LJ)
        live = sim.start(state)
        self.assertIsInstance(live, LiveRun)
        live.advance(50)
        self.assertFalse(self.root.exists())
        self.assertIsNone(state.atoms.velocities)
        self.assertEqual((state.step, live.step), (0, 50))
        self.assertAlmostEqual(live.time, .25)
        self.assertEqual(live.settings, dict(ensemble='nvt', timestep=.005, temperature=2., pressure=.01, heat_rate=0., barostat='isotropic'))
        observables = live.observables
        self.assertEqual((observables['step'], observables['time']), (live.step, live.time))
        self.assertAlmostEqual(observables['temperature'], np.sum(live.velocities**2)/(3*108-3))
        for array, shape in ((live.positions, (108, 3)), (live.velocities, (108, 3)), (live.box, (3,))):
            self.assertEqual(array.shape, shape)
            with self.assertRaises(ValueError):
                array[0] = 0
        before = live.positions
        live.advance(1)
        self.assertFalse(np.array_equal(before, live.positions))  # a snapshot, not a view
        for bad in (-1, 1.5, True):
            with self.assertRaises(ValueError):
                live.advance(bad)
        live.advance(0)
        self.assertEqual(live.step, 51)

    def test_continuing_from_a_returned_state(self):
        for name, state, settings in self.cases():
            first, second = (8, 6) if name == 'water' else (300, 200)
            with self.subTest(name):
                sim = MDSimulation(output_dir=self.root, **settings)
                live = sim.start(state)
                live.advance(first)
                middle = live.state
                self.assertEqual(middle.step, first)
                live.advance(second)
                # The returned state is independent of the run it came from.
                self.assertEqual((middle.step, live.step), (first, first + second))
                continued = sim.start(middle)
                self.assertEqual(continued.step, first)
                self.assertAlmostEqual(continued.time, middle.time)
                continued.advance(second)
                saved = sim.run(middle, steps=second, thermo_every=50, trajectory_every=None, storage_name='continued')
                for other in (continued.state, saved.state):
                    np.testing.assert_array_equal(other.atoms.positions, live.positions)
                    np.testing.assert_array_equal(other.atoms.velocities, live.velocities)
                    np.testing.assert_array_equal(other.box.lengths, live.box)
                    self.assertEqual(other.step, first + second)
                # And a saved run continues in memory.
                again = sim.start(Run.load(saved.path).state)
                again.advance(5)
                live.advance(5)
                np.testing.assert_array_equal(again.positions, live.positions)

    def test_a_control_changed_between_chunks(self):
        sim = MDSimulation(output_dir=self.root, ensemble='nvt', temperature=2., **LJ)
        live = sim.start(FCC(N=500, rho=.1))
        live.advance(400)
        middle = live.state
        live.set(temperature=.8)
        self.assertEqual(live.settings['temperature'], .8)
        self.assertEqual(sim.settings['temperature'], 2.)
        # As a new run from the same state with that control.
        other = sim.start(middle, temperature=.8)
        temperatures = []
        for _ in range(40):
            live.advance(50)
            temperatures.append(live.observables['temperature'])
        other.advance(2000)
        np.testing.assert_array_equal(live.positions, other.positions)
        np.testing.assert_array_equal(live.velocities, other.velocities)
        self.assertLess(abs(np.mean(temperatures[20:]) - .8), .05)
        # Several controls at once, checked together; a rejected change changes nothing.
        with self.assertRaises(ValueError):
            live.set(heat_rate=1.)
        with self.assertRaises(ValueError):
            live.set(ensemble='npt', pressure=-1.)
        with self.assertRaises(ValueError):
            live.set(timestep=0.)
        self.assertEqual(live.settings, dict(ensemble='nvt', timestep=.005, temperature=.8, pressure=.01, heat_rate=0., barostat='isotropic'))
        live.set(ensemble='nve', heat_rate=2., timestep=.002)
        energy = live.observables['total_energy']
        live.advance(100)
        self.assertAlmostEqual(live.observables['heat_added'], 2.*100*.002, places=12)
        self.assertAlmostEqual(live.observables['total_energy'] - energy, 2.*100*.002, delta=.05)
        self.assertAlmostEqual(live.time, 2400*.005 + 100*.002)

    def test_pressure_control_switched_on_between_chunks(self):
        sim = MDSimulation(output_dir=self.root, ensemble='nvt', temperature=2., **LJ)
        live = sim.start(FCC(N=256, rho=.3))
        live.advance(200)
        middle = live.state
        live.set(ensemble='npt', pressure=.4)
        other = sim.start(middle, ensemble='npt', pressure=.4)
        advance(live, 300, 11)
        other.advance(300)
        np.testing.assert_array_equal(live.positions, other.positions)
        np.testing.assert_array_equal(live.box, other.box)
        self.assertNotEqual(live.box[0], middle.box.lengths[0])
        # Back at constant volume the box stays, and the barostat starts from rest when switched on again.
        live.set(ensemble='nvt')
        box = live.box
        live.advance(50)
        np.testing.assert_array_equal(live.box, box)
        self.assertEqual(live.state._rate, 0.)

    def test_velocity_reversal_returns_to_the_start(self):
        sim = MDSimulation(output_dir=self.root, ensemble='nve', temperature=1., **LJ)
        melted = sim.start(FCC(N=108, rho=.3))
        melted.advance(300)
        start = melted.state
        forward = sim.start(start)
        forward.advance(200)
        turned = forward.state
        box = start.box.lengths
        moved = turned.atoms.positions - start.atoms.positions
        self.assertGreater(np.abs(moved - box*np.rint(moved/box)).max(), .5)
        backward = sim.start(turned.with_velocities(-turned.atoms.velocities))
        self.assertEqual(backward.step, turned.step)
        backward.advance(200)
        returned = backward.positions - start.atoms.positions
        self.assertLess(np.abs(returned - box*np.rint(returned/box)).max(), 1e-9)
        np.testing.assert_allclose(backward.velocities, -start.atoms.velocities, atol=1e-9)
        # The same through saved runs.
        there = sim.run(start, steps=200, thermo_every=100, trajectory_every=100, storage_name='there')
        back = sim.run(there.state.with_velocities(-there.state.atoms.velocities), steps=200,
                       thermo_every=100, trajectory_every=100, storage_name='back')
        np.testing.assert_array_equal(back.state.atoms.positions, backward.positions)
        # Saved trajectories hold the velocities of every frame.
        frames = back.load_frames()
        self.assertEqual(frames['velocities'].shape, (3, 108, 3))
        np.testing.assert_array_equal(frames['velocities'][0], -there.state.atoms.velocities)
        np.testing.assert_array_equal(frames['velocities'][-1], back.state.atoms.velocities)
        self.assertEqual([sorted(frame) for frame in back.iter_frames()][0],
                         ['box', 'positions', 'step', 'time', 'unwrapped', 'velocities'])

    def test_velocities_set_on_a_state(self):
        sim = MDSimulation(output_dir=self.root, ensemble='nvt', temperature=1., **LJ)
        live = sim.start(FCC(N=108, rho=.3))
        live.advance(100)
        state = live.state
        velocities = state.atoms.velocities.copy()
        velocities[0] += [1e-3, 0, 0]
        velocities -= velocities.mean(axis=0)
        perturbed = state.with_velocities(velocities)
        np.testing.assert_array_equal(perturbed.atoms.velocities, velocities)
        np.testing.assert_array_equal(perturbed.atoms.positions, state.atoms.positions)
        np.testing.assert_array_equal(state.atoms.velocities, live.velocities)  # the original is untouched
        # Clock and random sequence are kept, so only the perturbation differs.
        self.assertEqual((perturbed.step, perturbed.time, perturbed._rng), (state.step, state.time, state._rng))
        a, b = sim.start(state), sim.start(perturbed)
        a.advance(1); b.advance(1)
        self.assertLess(np.abs(a.positions - b.positions).max(), 1e-5)
        a.advance(2000); b.advance(2000)
        self.assertGreater(np.abs(a.velocities - b.velocities).max(), .1)
        for bad in (None, velocities[:5], np.full_like(velocities, np.nan)):
            with self.assertRaises(ValueError):
                state.with_velocities(bad)
        velocities[0] += [1., 0, 0]
        with self.assertRaisesRegex(ValueError, 'total momentum'):
            sim.start(state.with_velocities(velocities))
        rigid = Random(N=20, rho=.02, molecule='diatomic', bond=RigidBond(length=.7))
        moving = MDSimulation(ensemble='nve', temperature=1.).start(rigid).state
        along = moving.atoms.velocities.copy()
        along[0] += moving.atoms.positions[0] - moving.atoms.positions[1]
        with self.assertRaises(ValueError):
            moving.with_velocities(along)
        moving.with_velocities(-moving.atoms.velocities)

    def test_failure_and_interrupt(self):
        sim = MDSimulation(output_dir=self.root, ensemble='nve', temperature=1., **LJ)
        live = sim.start(FCC(N=108, rho=.8))
        native = _core.advance

        def interrupted(*args):
            result = list(native(*args[:6], min(args[6], 3), *args[7:]))
            result[4] = 1
            return tuple(result)
        with patch.object(_core, 'advance', interrupted):
            with self.assertRaises(KeyboardInterrupt):
                live.advance(10)
        # Whole steps were taken and counted; the run goes on.
        self.assertEqual(live.step, 3)
        live.advance(7)
        whole = sim.start(FCC(N=108, rho=.8))
        whole.advance(10)
        np.testing.assert_array_equal(live.positions, whole.positions)
        def failing(*args):
            raise FloatingPointError('Non-finite position; simulation is unstable.')
        with patch.object(_core, 'advance', failing):
            with self.assertRaises(FloatingPointError):
                live.advance(100)
        for use in (lambda: live.advance(1), lambda: live.state, lambda: live.set(timestep=.005)):
            with self.assertRaisesRegex(RuntimeError, 'cannot continue'):
                use()


class NeighbourListTests(unittest.TestCase):
    def test_cells_and_all_pairs_give_the_same_list(self):
        rng = np.random.default_rng(5)
        for atoms, box, radius, reaches in ((400, (12., 12., 12.), 2.9, (1, 2)), (700, (9., 15., 31.), 2.9, (1, 2)),
                                            (900, (8., 7.3, 8.6), 2.9, (2,)), (60, (500., 400., 300.), 2.9, (1, 2)),
                                            (2000, (14., 14., 14.), 2.9, (1, 2)), (300, (7., 7.3, 40.), 1.4, (1, 2))):
            with self.subTest(atoms=atoms, box=box):
                x = rng.uniform(0, 1, (atoms, 3))*box
                x[::7] += np.array(box)*rng.integers(-2, 3, (len(x[::7]), 3))  # also atoms outside the box
                x[1] = [0., 0., 0.]
                x[2] = np.nextafter(box, 0)
                x[3] = -1e-18
                every = pairs(x, box, radius, 0)
                for reach in reaches:
                    np.testing.assert_array_equal(pairs(x, box, radius, reach), every)
                dr = x[:, None, :] - x[None, :, :]
                dr -= box*np.rint(dr/box)
                i, j = np.nonzero(np.triu((dr**2).sum(axis=2) < radius**2, 1))
                self.assertEqual(len(every), len(i))
                np.testing.assert_array_equal(every, np.column_stack((i, j)))
        small = rng.uniform(0, 8, (50, 3))
        for box, reach in (((8., 8., 8.), 1), ((7., 8., 8.), 2)):
            with self.assertRaisesRegex(ValueError, 'too small for cells'):
                pairs(small, box, 2.9, reach)

    def test_kept_list_and_forces_follow_every_change(self):
        """One workspace through changes of positions, box, cutoff and atom count."""
        def arrays(atoms, rho, seed):
            state = Random(N=atoms, rho=rho, min_distance=.95, seed=seed)
            x = state.atoms.positions.copy()
            return [x, np.zeros_like(x), x.copy(), np.zeros_like(x), state.box.lengths.copy(), np.ones(atoms)]

        def forces(a, cutoff, workspace=None, steps=0):
            parameters, bonds = np.ones((len(a[0]), 2)), np.empty((0, 4))
            result = _core.advance(*a, steps, .002, cutoff, .4, 1., 0., 1, 0, 0., -1., 1., 0.,
                                   0, 1, .7, 0., 0., 0., 1., 1., parameters, bonds, workspace)
            return a[3].copy(), result[0]

        workspace = _core.workspace()
        a = arrays(600, .3, 1)
        forces(a, 2.5, workspace, steps=20)
        kept = forces(a, 2.5, workspace)
        fresh = forces(a, 2.5)
        np.testing.assert_array_equal(kept[0], fresh[0])
        self.assertEqual(kept[1], fresh[1])
        # An overwritten force array is restored from what was kept.
        a[3][:] = 7.
        np.testing.assert_array_equal(forces(a, 2.5, workspace)[0], fresh[0])
        changes = [lambda a: a[0].__setitem__(5, a[0][17] + [.97, 0, 0]),        # one atom moved far, next to another
                   lambda a: [a[k].__imul__(1.3) for k in (0, 2, 4)],            # box and positions scaled
                   lambda a: a[4].__imul__([1., 1.5, 1.])]                       # one box length alone
        for change in changes:
            change(a)
            kept, fresh = forces(a, 2.5, workspace), forces(a, 2.5)
            np.testing.assert_array_equal(kept[0], fresh[0])
            self.assertEqual(kept[1], fresh[1])
        for cutoff in (3.4, 1.5, 2.5):
            kept, fresh = forces(a, cutoff, workspace), forces(a, cutoff)
            np.testing.assert_array_equal(kept[0], fresh[0])
            self.assertEqual(kept[1], fresh[1])
        for atoms in (900, 150, 600):
            b = arrays(atoms, .3, atoms)
            forces(b, 2.5, workspace, steps=5)
            kept, fresh = forces(b, 2.5, workspace), forces(b, 2.5)
            np.testing.assert_array_equal(kept[0], fresh[0])
            self.assertEqual(kept[1], fresh[1])
        with self.assertRaises(ValueError):
            forces(a, 2.5, object())


if __name__ == '__main__':
    unittest.main()
