"""Still pictures in matplotlib: state.draw() and run.draw()."""
import importlib.util
import os
import tempfile
import unittest
import numpy as np

# matplotlib comes with the notebook extras; the wheel builds test without it.
if importlib.util.find_spec('matplotlib') is None:
    raise unittest.SkipTest('drawing needs matplotlib')
os.environ.setdefault('MPLBACKEND', 'Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import EllipseCollection, LineCollection
from fys2160_md import FCC, IceIh, MDSimulation, Plummer, Random, RigidWater


def balls(ax):
    collection, = [c for c in ax.collections if isinstance(c, EllipseCollection)]
    return collection


class Draw(unittest.TestCase):
    def tearDown(self):
        plt.close('all')

    def test_a_crystal_is_drawn_with_its_box(self):
        crystal = FCC(N=108, rho=.8)
        ax = crystal.draw()
        self.assertEqual(len(balls(ax).get_offsets()), 3*108)  # rim, body and highlight per atom
        edges = sum(len(c.get_segments()) for c in ax.collections if isinstance(c, LineCollection))
        self.assertEqual(edges, 12)
        self.assertFalse(ax.axison)
        ortho = crystal.draw(projection='orthographic')
        x0, x1 = ortho.get_xlim()
        self.assertGreaterEqual(x1-x0, crystal.box.lengths[0])  # face-on: the whole box width shows
        with self.assertRaisesRegex(ValueError, 'projection'):
            crystal.draw(projection='fisheye')

    def test_views_and_the_arrows_that_show_them(self):
        crystal = FCC(N=108, rho=.8)
        L = crystal.box.lengths[0]
        front = crystal.draw(view='front', projection='orthographic')
        labels = sorted(text.get_text() for text in front.texts if text.get_text())
        self.assertEqual(labels, ['x', 'y', 'z'])
        for ax, limits in ((front, front.get_xlim()), (front, front.get_ylim())):
            self.assertTrue(L <= np.ptp(limits) <= L+1.2)  # face-on: the box, and at most the atoms' overhang
        top = crystal.draw(view='top', projection='orthographic')
        self.assertTrue(L <= np.ptp(top.get_ylim()) <= L+1.2)
        crystal.draw(view=(30, -30))
        self.assertEqual(len(crystal.draw(arrows=False).texts), 0)
        for bad in ('below', (1, 2, 3), 'xy'):
            with self.assertRaisesRegex(ValueError, 'view'):
                crystal.draw(view=bad)

    def test_a_flat_box_with_walls_has_axes_and_solid_walls(self):
        gas = Random(N=50, box=(30, 20), boundary=('walls', 'periodic'), min_distance=1.0, seed=1)
        ax = gas.draw()
        self.assertTrue(ax.axison)
        x0, x1 = ax.get_xlim()
        self.assertLessEqual(x0, 0)
        self.assertGreaterEqual(x1, 30)
        styles = sorted(line.get_linestyle() for line in ax.lines)
        self.assertEqual(styles.count('-'), 2)  # the two walls; the periodic sides are dashed

    def test_water_molecules_are_drawn_whole(self):
        ice = IceIh(N=96, density=.92, molecule='water', bond=RigidWater(oh_length=1.0, hoh_angle=109.47), seed=1)
        ax = ice.draw()
        bonds, = [c for c in ax.collections if isinstance(c, LineCollection) and len(c.get_segments()) == 192]
        lengths = [np.linalg.norm(np.subtract(*segment)) for segment in bonds.get_segments()]
        self.assertLess(max(lengths), 3.)  # no bond drawn across the box
        sizes = balls(ax).get_offsets()
        self.assertEqual(len(sizes), 3*288)

    def test_radius_choices(self):
        mixture = Random(species={'A': 20, 'B': 10}, rho=.05, min_distance=1.0, seed=2)
        mixture.draw(radius={'A': .5, 'B': .7})
        with self.assertRaisesRegex(ValueError, "species 'B'"):
            mixture.draw(radius={'A': .5})
        with self.assertRaisesRegex(ValueError, 'positive'):
            mixture.draw(radius=-1)
        Plummer(N=64, seed=1).draw()  # an open cluster: no box

    def test_a_saved_run_draws_any_frame(self):
        with tempfile.TemporaryDirectory() as folder:
            sim = MDSimulation(pair_potential='lj', ensemble='nvt', temperature=1.0, timestep=.005,
                               cutoff=2.5, output_dir=folder)
            run = sim.run(FCC(N=108, rho=.8), steps=20, thermo_every=10, trajectory_every=10)
            fig, axes = plt.subplots(1, 2)
            self.assertIs(run.draw(0, ax=axes[0]), axes[0])
            run.draw(ax=axes[1], projection='orthographic')
            self.assertEqual(len(balls(axes[0]).get_offsets()), 3*108)


if __name__ == '__main__':
    unittest.main()
