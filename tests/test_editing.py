"""Everything in a state can be changed, and a running simulation can be nudged."""
import unittest
import numpy as np
from fys2160_md import FCC, MDSimulation, Random


def sim(**options):
    settings = dict(pair_potential='lj', ensemble='nve', temperature=1.0, timestep=.005, cutoff=2.5)
    settings.update(options)
    return MDSimulation(**settings)


class StateEdits(unittest.TestCase):
    def test_positions_and_masses_are_the_states_own(self):
        state = FCC(N=108, rho=.5)
        state.positions[0, 0] += .25
        self.assertEqual(state.atoms.positions[0, 0], state.positions[0, 0])
        state.masses[:10] = 4.
        live = sim().start(state)
        np.testing.assert_array_equal(live.masses[:10], 4.)
        self.assertAlmostEqual(live.positions[0, 0], state.positions[0, 0])
        state.positions = state.positions + 100  # outside the box: wrapped back in by the next run
        self.assertTrue(np.all((sim().start(state).positions >= 0) & (sim().start(state).positions < state.box)))
        with self.assertRaisesRegex(ValueError, 'one row per atom'):
            state.positions = np.zeros((3, 3))

    def test_velocities_can_be_set_scaled_and_kicked(self):
        state = FCC(N=108, rho=.5)
        self.assertIsNone(state.velocities)
        state.velocities = sim().start(state, seed=1).velocities  # give a fresh state velocities
        before = sim().start(state).observables['temperature']
        state.velocities *= 2
        self.assertAlmostEqual(sim().start(state).observables['temperature'], 4*before)
        state.velocities[0] += (1., 0, 0)
        with self.assertRaisesRegex(ValueError, 'total momentum.*np.average'):
            sim().start(state)
        state.velocities -= np.average(state.velocities, axis=0, weights=state.masses)
        sim().start(state).advance(10)
        state.velocities = None  # back to random velocities at the run's temperature
        self.assertAlmostEqual(sim(temperature=3.).start(state).observables['temperature'], 3.)

    def test_species_boundary_and_clock(self):
        state = Random(N=100, rho=.05, min_distance=1.0, seed=1)
        state.atoms.species[:40] = 'Heavy'
        self.assertEqual(state.species, ('Heavy', 'A'))
        mixed = sim(epsilon={'Heavy': 1.5, 'A': 1.}, sigma={'Heavy': 1., 'A': 1.})
        mixed.start(state).advance(10)
        with self.assertRaisesRegex(ValueError, 'Heavy'):
            sim(epsilon={'A': 1.}, sigma={'A': 1.}).start(state)
        state.boundary = 'walls'
        self.assertEqual(len(state.walls), 6)
        sim(ensemble='nvt').start(state).advance(10)
        state.time, state.step = 5., 1000
        self.assertEqual(sim().start(state).step, 1000)
        with self.assertRaises(ValueError):
            state.step = -1

    def test_a_run_never_changes_the_state_it_is_given(self):
        state = FCC(N=108, rho=.5)
        state.velocities = sim().start(state, seed=1).velocities
        x, v = state.positions.copy(), state.velocities.copy()
        run = sim().start(state)
        run.advance(100)
        np.testing.assert_array_equal(state.positions, x)
        np.testing.assert_array_equal(state.velocities, v)


class LiveEdits(unittest.TestCase):
    def test_snapshots_stay_and_assignment_changes_the_run(self):
        live = sim(ensemble='nvt').start(FCC(N=108, rho=.5), seed=2)
        before = live.positions
        live.advance(50)
        self.assertFalse(np.array_equal(before, live.positions))  # a snapshot does not follow the run
        with self.assertRaises(ValueError):
            live.velocities[0] = 0  # read-only: assign instead
        hot = 4*live.velocities
        live.velocities = hot
        np.testing.assert_array_equal(live.velocities, hot)
        x = live.positions.copy()
        x[0] += .1
        unwrapped = live.state._u[0].copy()
        live.positions = x
        np.testing.assert_allclose(live.state._u[0]-unwrapped, .1)  # displacements stay continuous
        live.box = live.box*1.1
        live.advance(20)
        self.assertAlmostEqual(live.box[0], 1.1*FCC(N=108, rho=.5).box.lengths[0])
        live.masses = np.full(108, 2.)
        live.advance(10)
        with self.assertRaisesRegex(ValueError, 'shape'):
            live.positions = x[:5]


if __name__ == '__main__':
    unittest.main()
