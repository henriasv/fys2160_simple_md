"""Changing a state's box in place: state.box[0] *= 2."""
import tempfile
import unittest
import numpy as np
from fys2160_md import FCC, MDSimulation, Random


def sim():
    return MDSimulation(pair_potential='lj', ensemble='nvt', temperature=1.0, timestep=0.005, cutoff=2.5)


class BoxEdit(unittest.TestCase):
    def test_a_periodic_box_can_grow_and_shrink(self):
        state = FCC(N=256, rho=.5)
        L = state.box.lengths[0]
        x = state.atoms.positions.copy()
        state.box[0] *= 2
        self.assertEqual(state.box.lengths.tolist(), [2*L, L, L])
        self.assertAlmostEqual(state.box.volume, 2*L**3)
        np.testing.assert_array_equal(state.atoms.positions, x)  # atoms stay where they are
        live = sim().start(state)
        live.advance(10)
        self.assertEqual(live.box.tolist(), [2*L, L, L])
        state.box[:] = .8*L
        squeezed = sim().start(state)  # atoms outside the smaller box are wrapped back in
        self.assertTrue(np.all(squeezed.positions >= 0) and np.all(squeezed.positions < .8*L))

    def test_an_edited_run_state_is_kept_for_the_next_run(self):
        with tempfile.TemporaryDirectory() as root:
            simulation = MDSimulation(pair_potential='lj', ensemble='nvt', temperature=1.0, timestep=0.005,
                                      cutoff=2.5, output_dir=root)
            first = simulation.run(FCC(N=108, rho=.2), steps=20, thermo_every=10, trajectory_every=None)
            L = first.state.box.lengths[0]
            first.state.box[0] *= 2
            second = simulation.run(first.state, steps=20, thermo_every=10, trajectory_every=None,
                                    storage_name='second')
            self.assertEqual(second.state.box.lengths.tolist(), [2*L, L, L])

    def test_between_walls_the_far_wall_moves(self):
        state = Random(N=50, box=(20, 20), boundary='walls', min_distance=1.0, seed=2)
        state.box[0] *= 2
        live = sim().start(state)
        live.advance(10)
        self.assertEqual(live.box.tolist(), [40., 20.])
        state.box[0] = 5
        with self.assertRaisesRegex(ValueError, 'Piston'):
            sim().start(state)


if __name__ == '__main__':
    unittest.main()
