"""Tilted (triclinic) boxes, strain and the pressure tensor for unbonded atoms."""
import itertools
import unittest
import numpy as np
from fys2160_md import FCC, MDSimulation, Random, System, _core


def sim(**options):
    settings = dict(pair_potential='lj', ensemble='nve', temperature=1.0, timestep=.005, cutoff=2.5)
    settings.update(options)
    return MDSimulation(**settings)


def edges(state):
    (lx, ly, lz), (xy, xz, yz) = state.box.lengths, state.tilt
    return np.array([[lx, 0, 0], [xy, ly, 0], [xz, yz, lz]])


def brute_force(state, cutoff=2.5):
    """Shifted 12-6 energy and forces summed over periodic images, the slow way."""
    x = state.atoms.positions
    a, b, c = edges(state)
    shift = 4*(cutoff**-12-cutoff**-6)
    energy, forces = 0., np.zeros_like(x)
    for i, j in itertools.combinations(range(len(x)), 2):
        for n in itertools.product((-1, 0, 1), repeat=3):
            r = x[i]-x[j]-(n[0]*a+n[1]*b+n[2]*c)
            r2 = r@r
            if r2 < cutoff**2:
                inv6 = r2**-3
                energy += 4*inv6*(inv6-1)-shift
                f = 24*inv6*(2*inv6-1)/r2*r
                forces[i] += f
                forces[j] -= f
    return energy, forces


def tilted_gas(seed=1, scale=1., count=40):
    box = scale*np.array([7., 6.5, 6.])
    tilt = scale*np.array([2.1, -1.4, 1.2])
    rng = np.random.default_rng(seed)
    x = rng.random((count, 3))@np.array([[box[0], 0, 0], [tilt[0], box[1], 0], [tilt[1], tilt[2], box[2]]])
    # Keep atoms apart, so the forces stay moderate.
    keep = [0]
    for i in range(1, len(x)):
        if all(np.linalg.norm(x[i]-x[k]) > 1.0 for k in keep):
            keep.append(i)
    # At rest, so that the pressure tensor is the pair virial alone.
    return System(x[keep], box, tilt=tilt, velocities=np.zeros((len(keep), 3)))


class Triclinic(unittest.TestCase):
    def test_forces_match_a_sum_over_images(self):
        state = tilted_gas()
        engine = sim()._engine(state)
        engine._advance(0)
        energy, forces = brute_force(engine.state())
        self.assertAlmostEqual(engine._potential, energy, places=10)
        np.testing.assert_allclose(engine._f, forces, atol=1e-9)

    def test_a_whole_lattice_vector_of_tilt_is_the_same_crystal(self):
        crystal = FCC(N=500, rho=1.0)
        a = crystal.box.lengths[0]/5  # one unit cell
        leaning = crystal.copy()
        leaning.tilt[0] = a          # the y edge leans by one cell along x
        leaning = leaning.copy()      # through the constructor: positions wrapped into the tilted box
        energies = []
        for state in (crystal, leaning):
            engine = sim()._engine(state)
            engine._advance(0)
            energies.append(engine._potential)
            self.assertLess(np.abs(engine._f).max(), 1e-9)  # still a perfect lattice
        self.assertAlmostEqual(energies[0], energies[1], places=9)

    def test_the_pressure_tensor_is_the_strain_derivative_of_the_energy(self):
        state = tilted_gas(seed=2)
        V = np.prod(state.box.lengths)
        engine = sim()._engine(state)
        P = engine.observables
        h = 1e-6
        for component, key in (('xx', 'pxx'), ('yy', 'pyy'), ('zz', 'pzz'), ('xy', 'pxy'), ('xz', 'pxz'), ('yz', 'pyz')):
            energy = []
            for sign in (1, -1):
                strained = state.copy()
                strained.strain(**{component: sign*h})
                e = sim()._engine(strained)
                e._advance(0)
                energy.append(e._potential)
            derivative = (energy[0]-energy[1])/(2*h)
            # Without velocities, V P_ab is the pair virial, -dU/d(strain).
            self.assertAlmostEqual(P[key]*V, -derivative, delta=1e-5*max(1, abs(derivative)), msg=component)

    def test_cells_find_the_pairs_of_all_pairs(self):
        state = tilted_gas(seed=3, scale=2., count=300)
        x = np.ascontiguousarray(state.atoms.positions)
        box, tilt = tuple(state.box.lengths), tuple(state.tilt)
        pairs = lambda reach: sorted(np.frombuffer(_core.neighbor_pairs(x, box, 2.9, reach, 3, (1, 1, 1), tilt),
                                                   dtype=np.int32).reshape(-1, 2).tolist())
        everything = pairs(0)
        self.assertGreater(len(everything), 20)
        for reach in (1, 2):
            self.assertEqual(pairs(reach), everything)

    def test_chunks_and_chains_are_one_trajectory(self):
        crystal = FCC(N=500, rho=1.0)
        crystal.strain(xy=.05, yz=-.03)
        md = sim(ensemble='nvt', temperature=.5)
        whole = md.start(crystal, seed=4)
        whole.advance(200)
        pieces = md.start(crystal, seed=4)
        for _ in range(8):
            pieces.advance(25)
        np.testing.assert_array_equal(whole.positions, pieces.positions)
        chained = md.start(md.start(crystal, seed=4).state, seed=4)  # a state carries its tilt
        self.assertEqual(chained.state.tilt.tolist(), crystal.tilt.tolist())

    def test_the_barostat_scales_the_tilt_with_the_box(self):
        crystal = FCC(N=500, rho=1.0)
        crystal.strain(xy=.05)
        live = sim(ensemble='npt', temperature=.5, pressure=2.0).start(crystal, seed=5)
        live.advance(300)
        state = live.state
        self.assertNotAlmostEqual(state.box.lengths[0], crystal.box.lengths[0], places=3)
        self.assertAlmostEqual(state.tilt[0]/state.box.lengths[1], crystal.tilt[0]/crystal.box.lengths[1], places=12)

    def test_limits_and_refusals(self):
        crystal = FCC(N=108, rho=1.0)
        with self.assertRaisesRegex(ValueError, 'at most half'):
            System(crystal.atoms.positions, crystal.box.lengths, tilt=(3., 0, 0))
        with self.assertRaisesRegex(ValueError, 'at most half'):
            crystal.copy().strain(xy=.6)
        with self.assertRaisesRegex(ValueError, 'periodic three-dimensional box'):
            Random(N=20, box=(20, 20), seed=1, min_distance=1.0).copy().strain(xy=.1)
        leaning = crystal.copy()
        leaning.strain(xy=.1)
        with self.assertRaisesRegex(ValueError, 'rectangular'):
            leaning.with_box(crystal.box.lengths*2)
        small = FCC(N=108, rho=1.0)  # 4.76 wide: a strong tilt brings the faces within twice the cutoff
        small.strain(xy=.45)
        with self.assertRaisesRegex(ValueError, 'faces'):
            sim(cutoff=2.3).start(small)


if __name__ == '__main__':
    unittest.main()
