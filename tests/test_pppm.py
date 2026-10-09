"""PPPM electrostatics against direct Ewald sums and the NaCl Madelung constant."""
import tempfile
import unittest
import numpy as np
from fys2160_md import System, MDSimulation
from fys2160_md import _core
from ewald_reference import ewald

MADELUNG_NACL = 1.747564594633


def rock_salt(cells, spacing=1.):
    sites = np.indices((2*cells,)*3).reshape(3, -1).T
    charges = np.where(sites.sum(axis=1) % 2, -1., 1.)
    return sites*spacing + .25, charges, np.full(3, 2*cells*spacing)


class PPPMTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)

    def charged(self, x, q, box, *, units='reduced', accuracy=1e-6, cutoff=None, **kwargs):
        labels = np.where(q > 0, 'Na', 'Cl')
        system = System(x, box, species=labels, units=units, masses=np.ones(len(x)))
        cutoff = cutoff or .45*min(box)
        return MDSimulation(units=units, epsilon={'Na': 0., 'Cl': 0.}, sigma={'Na': 0., 'Cl': 0.},
                            charge={'Na': 1., 'Cl': -1.}, coulomb='pppm', kspace_accuracy=accuracy,
                            cutoff=cutoff, temperature=1. if units == 'reduced' else 300., timestep=.001 if units == 'reduced' else 1.,
                            ensemble='nve', output_dir=self.temp.name, **kwargs)._engine(system)

    def test_fft_matches_numpy(self):
        rng = np.random.default_rng(4)
        for shape in [(2, 3, 5), (8, 9, 10), (12, 15, 16), (27, 25, 32)]:
            data = rng.normal(size=shape) + 1j*rng.normal(size=shape)
            forward, inverse = data.copy(), data.copy()
            _core.fft_transform(forward, -1)
            _core.fft_transform(inverse, 1)
            np.testing.assert_allclose(forward, np.fft.fftn(data), atol=1e-11)
            np.testing.assert_allclose(inverse, np.fft.ifftn(data)*data.size, atol=1e-11)
        with self.assertRaises(ValueError):
            _core.fft_transform(np.zeros((7, 4, 4), complex), -1)

    def test_parameter_choice_matches_lammps_and_is_factorable(self):
        # Reference meshes and g from LAMMPS (30 Mar 2026) for 120 unit charges in this box.
        for accuracy, mesh, g in [(1e-4, (15, 15, 15), .30850154), (1e-5, (20, 24, 24), .34180908), (1e-6, (36, 36, 36), .37894025)]:
            chosen, g_ewald, estimate = _core.pppm_parameters((21., 22., 23.), 120, 120., 332.06371, 9., accuracy, 5)
            self.assertEqual(chosen, mesh)
            self.assertAlmostEqual(g_ewald, g, places=7)
            self.assertLess(estimate, 2*accuracy)

    def test_nacl_madelung_constant(self):
        x, q, box = rock_salt(2)
        # Energy per ion pair is -M e^2/r0 for nearest-neighbour distance r0 = 1.
        reference, _, _ = ewald(x, q, box, alpha=1.5, real_cutoff=4., kmax=12)
        self.assertAlmostEqual(reference/(len(q)/2), -MADELUNG_NACL, places=9)
        sim = self.charged(x, q, box, accuracy=1e-7, cutoff=1.9)
        energy = sim.observables['coulomb_energy']
        self.assertAlmostEqual(energy/(len(q)/2), -MADELUNG_NACL, places=5)
        np.testing.assert_allclose(sim.forces, 0, atol=1e-5)

    def test_random_charges_match_direct_ewald(self):
        rng = np.random.default_rng(3)
        box = np.array([14., 15., 16.])
        x = []
        while len(x) < 60:
            p = rng.uniform(0, box)
            if all(np.linalg.norm((p-y)-box*np.rint((p-y)/box)) > 2. for y in x):
                x.append(p)
        x = np.array(x); q = np.where(np.arange(60) % 2, 1., -1.)
        energy, forces, virial = ewald(x, q, box, alpha=.5, coulomb=332.06371)
        for accuracy in (1e-4, 1e-6):
            sim = self.charged(x, q, box, units='real', accuracy=accuracy, cutoff=6.5)
            rms = np.sqrt(np.mean((sim.forces-forces)**2))/332.06371
            self.assertLess(rms, 2*accuracy)
            self.assertLess(abs(sim.observables['coulomb_energy']-energy)/abs(energy), 20*accuracy)
            self.assertLess(abs(sim._virial-virial)/abs(virial), 20*accuracy)
            np.testing.assert_allclose(sim.forces.sum(axis=0), 0, atol=1e-9)
            self.assertLess(sim._pppm_state[7], 2*accuracy)

    def test_neutrality_and_explicit_settings_are_required(self):
        x, q, box = rock_salt(1)
        system = System(x, box, species=np.where(q > 0, 'Na', 'Cl'))
        common = dict(epsilon={'Na': 0., 'Cl': 0.}, sigma={'Na': 0., 'Cl': 0.}, cutoff=.9, output_dir=self.temp.name)
        with self.assertRaises(ValueError):
            MDSimulation(charge={'Na': 1., 'Cl': -.5}, coulomb='pppm', **common)._engine(system)
        with self.assertRaises(ValueError):
            MDSimulation(charge={'Na': 1., 'Cl': -1.}, **common)
        with self.assertRaises(ValueError):
            MDSimulation(coulomb='pppm', **common)
        with self.assertRaises(ValueError):
            MDSimulation(charge={'Na': 1., 'Cl': -1.}, coulomb='pppm', **common)._engine(
                System(x, box, species=np.where(q > 0, 'Na', 'Cl'), units='real'))


if __name__ == '__main__':
    unittest.main()
