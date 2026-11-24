"""Per-step measurements (measure=) and the sampled arrays of a saved run."""
import tempfile
import unittest
import numpy as np
from fys2160_md import FCC, MDSimulation, Plummer, Random, RigidWater


def left_half(result):
    return np.sum(result.positions[:, :, 0] < result.boxes[0, 0]/2, axis=1)


def flat_gas():
    return Random(N=60, box=(30, 30), boundary='walls', min_distance=1.0, seed=3)


class Measure(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.sim = MDSimulation(pair_potential='lj', ensemble='nvt', temperature=2.0, timestep=0.005,
                                cutoff=2.5, friction=0.05, output_dir=self.directory.name)

    def tearDown(self):
        self.directory.cleanup()

    def test_the_same_function_on_the_saved_run_and_during_it(self):
        run = self.sim.run(flat_gas(), steps=60, thermo_every=10, trajectory_every=1,
                           measure={'N_left': left_half})
        np.testing.assert_array_equal(run.positions, run.load_frames()['positions'])
        self.assertEqual(run.positions.shape, (61, 60, 2))
        self.assertEqual(run.boxes.shape, (61, 2))
        np.testing.assert_array_equal(run.steps, np.arange(61))
        np.testing.assert_array_equal(run.measured_steps, np.arange(61))
        np.testing.assert_array_equal(run.measured['N_left'], left_half(run))
        # Each thermo row after the first is the mean over the steps since the previous row.
        thermo = run.thermo
        series = run.measured['N_left']
        self.assertEqual(thermo['N_left'][0], series[0])
        for row in range(1, len(thermo)):
            self.assertAlmostEqual(thermo['N_left'][row], series[10*row-9:10*row+1].mean())
        self.assertEqual(run.metadata['measure'], ['N_left'])

    def test_measuring_does_not_change_the_trajectory(self):
        plain = self.sim.run(flat_gas(), steps=300, thermo_every=10, trajectory_every=None, storage_name='plain')
        measured = self.sim.run(flat_gas(), steps=300, thermo_every=10, trajectory_every=None, storage_name='measured',
                                measure={'N_left': left_half})
        np.testing.assert_array_equal(plain.state.atoms.positions, measured.state.atoms.positions)
        np.testing.assert_array_equal(plain.state.atoms.velocities, measured.state.atoms.velocities)

    def test_start_gives_the_series_of_run(self):
        run = self.sim.run(flat_gas(), steps=100, thermo_every=10, trajectory_every=None,
                           measure={'N_left': left_half})
        live = self.sim.start(flat_gas(), measure={'N_left': left_half})
        for _ in range(4):
            live.advance(25)
        np.testing.assert_array_equal(live.measured['N_left'], run.measured['N_left'])
        np.testing.assert_array_equal(live.measured_steps, run.measured_steps)
        self.assertEqual(live.observables['N_left'], np.sum(live.positions[:, 0] < live.box[0]/2))

    def test_thermo_is_on_disk_during_the_run(self):
        # Rows are flushed at most once a second, but the header and the first row at once.
        path = f'{self.directory.name}/readable/thermo.csv'
        seen = []

        def lines_on_disk(result):
            try:
                with open(path) as stream:
                    seen.append(len(stream.readlines()))
            except FileNotFoundError:  # the starting state is measured before the file is made
                seen.append(None)
            return np.zeros(len(result.steps))

        run = self.sim.run(flat_gas(), steps=50, thermo_every=1, trajectory_every=None,
                           storage_name='readable', measure={'lines': lines_on_disk})
        self.assertEqual(seen[1:3], [2, 2])
        self.assertEqual(len(run.thermo), 51)

    def test_every_attribute_and_a_changing_box(self):
        def length(result):
            self.assertEqual(result.positions.shape, result.velocities.shape)
            self.assertEqual(result.unwrapped.shape, result.positions.shape)
            self.assertEqual(result.boxes.shape, (len(result.steps), 3))
            np.testing.assert_allclose(result.times, result.steps*0.005)
            self.assertFalse(result.positions.flags.writeable)
            return result.boxes[:, 0]
        sim = MDSimulation(pair_potential='lj', ensemble='npt', temperature=1.0, pressure=1.0, timestep=0.005,
                           cutoff=2.5, output_dir=self.directory.name)
        run = sim.run(FCC(N=108, rho=.7), steps=50, thermo_every=10, trajectory_every=10,
                      measure={'length': length})
        lengths = run.measured['length']
        self.assertGreater(np.ptp(lengths), 0)  # the barostat moves the box
        np.testing.assert_array_equal(lengths[run.steps], run.boxes[:, 0])

    def test_water_and_gravity_can_be_measured(self):
        spce = RigidWater(oh_length=1.0, hoh_angle=109.47)
        water = Random(N=64, box=16.0, molecule='water', bond=spce, min_distance=2.4, seed=2)
        sim = MDSimulation(units='real', epsilon={'O': .15535, 'H': 0.}, sigma={'O': 3.166, 'H': 0.},
                           lj_shift=False, ensemble='nvt', temperature=300, timestep=1.0, cutoff=5.0,
                           output_dir=self.directory.name)
        live = sim.start(water, measure={'O_left': lambda r: np.sum(r.positions[:, ::3, 0] < r.boxes[:, None, 0]/2, axis=1)})
        live.advance(20)
        self.assertEqual(len(live.measured['O_left']), 21)
        cluster = Plummer(N=64, seed=4)
        gravity = MDSimulation(pair_potential='gravity', ensemble='nve', timestep=0.01, cutoff=None,
                               output_dir=self.directory.name)
        spread = gravity.start(cluster, measure={'spread': lambda r: r.unwrapped[:, :, 0].std(axis=1)})
        spread.advance(10)
        self.assertEqual(len(spread.measured['spread']), 11)

    def test_mistakes_are_explained(self):
        with self.assertRaisesRegex(ValueError, 'one value per step.*axis=1'):
            self.sim.start(flat_gas(), measure={'total': lambda r: np.sum(r.positions[:, :, 0] < 15)})
        with self.assertRaisesRegex(ValueError, 'not finite'):
            self.sim.start(flat_gas(), measure={'bad': lambda r: np.full(len(r.steps), np.nan)})
        with self.assertRaisesRegex(TypeError, 'one argument, the result'):
            self.sim.start(flat_gas(), measure={'bad': lambda positions, box: positions})
        with self.assertRaisesRegex(ValueError, 'is taken'):
            self.sim.start(flat_gas(), measure={'temperature': left_half})
        with self.assertRaisesRegex(TypeError, 'dictionary'):
            self.sim.start(flat_gas(), measure=left_half)


if __name__ == '__main__':
    unittest.main()
