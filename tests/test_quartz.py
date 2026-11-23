"""Alpha-quartz with Buckingham terms, charges in tilted boxes, and anisotropic pressure control."""
import tempfile
import unittest
import numpy as np
from fys2160_md import AlphaQuartz, Buckingham, MDSimulation, Run, System
from ewald_reference import ewald

BKS = {('O', 'O'): Buckingham(A=1388.7730, rho=0.362319, C=175.0),
       ('Si', 'O'): Buckingham(A=18003.7572, rho=0.205205, C=133.5381)}
CHARGE = {'Si': 2.4, 'O': -1.2}
BAR = 1.6021765e6  # bar per eV/Å^3


def bks(**options):
    settings = dict(units='metal', pair_potential=BKS, charge=CHARGE, coulomb='pppm', ensemble='nvt',
                    temperature=300, timestep=0.001, cutoff=8.0)
    settings.update(options)
    return MDSimulation(**settings)


def neighbours(state, cutoff=1.7):
    x, box = state.positions, state.box.lengths
    d = x[:, None]-x[None]
    d -= box*np.rint(d/box)
    r = np.linalg.norm(d, axis=2)
    np.fill_diagonal(r, np.inf)
    return d, r


class Structure(unittest.TestCase):
    def test_quartz_is_built_of_corner_sharing_tetrahedra(self):
        quartz = AlphaQuartz()
        self.assertEqual(len(quartz.atoms), 432)
        self.assertEqual(quartz.species, ('Si', 'O'))
        self.assertEqual(quartz.units, 'metal')
        si = quartz.atoms.species == 'Si'
        self.assertEqual(si.sum(), 144)
        d, r = neighbours(quartz)
        self.assertTrue(np.all((r[si][:, ~si] < 1.7).sum(axis=1) == 4))  # four O around each Si
        self.assertTrue(np.all((r[~si][:, si] < 1.7).sum(axis=1) == 2))  # each O joins two Si
        self.assertAlmostEqual(r.min(), 1.604, places=3)
        density = quartz.masses.sum()/(0.602214076*quartz.box.volume)
        self.assertAlmostEqual(density, 2.6487, places=3)  # g/cm^3, measured 2.648
        o = np.flatnonzero(~si)[0]
        a, b = np.flatnonzero(si & (r[o] < 1.7))
        angle = np.degrees(np.arccos(d[a, o]@d[b, o]/r[a, o]/r[b, o]))
        self.assertAlmostEqual(angle, 143.7, places=1)
        self.assertEqual(len(AlphaQuartz(cells=(1, 1, 1)).atoms), 18)
        with self.assertRaises(ValueError):
            AlphaQuartz(cells=(2, 2))


class Terms(unittest.TestCase):
    def test_buckingham_barrier_and_refusals(self):
        si_o = BKS[('Si', 'O')]
        self.assertAlmostEqual(si_o.barrier, 1.0983, places=4)
        slope = (si_o.energy(si_o.barrier+1e-6)-si_o.energy(si_o.barrier-1e-6))/2e-6
        self.assertAlmostEqual(slope, 0, places=3)  # the inner maximum
        self.assertEqual(Buckingham(A=1., rho=.3, C=0.).barrier, 0.)
        with self.assertRaisesRegex(ValueError, 'attractive at every distance'):
            Buckingham(A=1e-3, rho=.1, C=100.)
        with self.assertRaises(ValueError):
            Buckingham(A=1., rho=0., C=1.)
        with self.assertRaisesRegex(ValueError, 'twice'):
            bks(pair_potential={('Si', 'O'): si_o, ('O', 'Si'): si_o})
        with self.assertRaisesRegex(ValueError, 'pair of species'):
            bks(pair_potential={'Si': si_o})
        with self.assertRaisesRegex(ValueError, 'leave out epsilon'):
            bks(epsilon={'Si': 1., 'O': 1.})
        with self.assertRaisesRegex(ValueError, 'Ge'):
            bks(pair_potential={**BKS, ('Ge', 'O'): si_o}).start(AlphaQuartz(cells=(4, 2, 3)))
        with self.assertRaisesRegex(ValueError, 'explicit'):
            MDSimulation(units='metal', pair_potential=BKS, temperature=300)

    def test_atoms_inside_the_barrier_stop_the_run(self):
        pair = System([[5., 5., 5.], [6.05, 5., 5.]], (20., 20., 20.), species=['Si', 'O'], masses=[28.0855, 15.9994],
                      units='metal', velocities=np.zeros((2, 3)))
        with self.assertRaisesRegex(ValueError, 'barrier'):
            MDSimulation(units='metal', pair_potential=BKS, ensemble='nve', temperature=300, timestep=.001,
                         cutoff=8.).start(pair)


class Tail(unittest.TestCase):
    def test_the_tail_adds_the_pairs_beyond_the_cutoff(self):
        quartz = AlphaQuartz(cells=(4, 2, 3))
        quartz.velocities = np.zeros_like(quartz.positions)
        cut, tail = (bks(ensemble='nve', lj_shift=False, tail=flag)._engine(quartz).observables for flag in (False, True))
        V, rc = quartz.box.volume, 8.0
        r = np.linspace(rc, 80., 400_001)
        counts = {'Si': 144, 'O': 288}
        energy = virial = 0.
        for (a, b), term in BKS.items():
            pairs = counts[a]*counts[b]*(2 if a != b else 1)  # ordered pairs of atoms
            U = term.energy(r)
            dU = -term.A/term.rho*np.exp(-r/term.rho) + 6*term.C/r**7
            # Beyond the grid only -C/r^6 is left: its integrals there are -C/(3R^3) and 2C/R^3.
            energy += 2*np.pi/V*pairs*(np.trapezoid(r**2*U, r) - term.C/(3*r[-1]**3))
            virial -= 2*np.pi/V*pairs*(np.trapezoid(r**3*dU, r) + 2*term.C/r[-1]**3)
        self.assertAlmostEqual(tail['potential_energy']-cut['potential_energy'], energy, delta=1e-6*abs(energy))
        for name in ('pxx', 'pyy', 'pzz'):
            self.assertAlmostEqual((tail[name]-cut[name])/BAR, virial/(3*V), delta=1e-6*abs(virial/V))
        self.assertAlmostEqual(tail['pxy'], cut['pxy'], places=9)

    def test_lennard_jones_tail_and_refusals(self):
        x = np.indices((4, 4, 4)).reshape(3, -1).T*3.2 + .5
        salt = System(x, (12.8, 12.8, 12.8), species=['Na', 'Cl']*32, masses=[22.99, 35.45]*32, units='real',
                      velocities=np.zeros((64, 3)))
        sim = lambda tail: MDSimulation(units='real', epsilon={'Na': .1, 'Cl': .15}, sigma={'Na': 2.5, 'Cl': 4.},
                                        charge={'Na': 1, 'Cl': -1}, coulomb='pppm', temperature=300, timestep=1,
                                        cutoff=6., ensemble='nve', lj_shift=False, tail=tail)._engine(salt)
        difference = sim(True)._potential-sim(False)._potential
        V, rc = 12.8**3, 6.
        expected = 0.
        for a, na in (('Na', 32), ('Cl', 32)):
            for b, nb in (('Na', 32), ('Cl', 32)):
                eps = np.sqrt({'Na': .1, 'Cl': .15}[a]*{'Na': .1, 'Cl': .15}[b])
                sig = (({'Na': 2.5, 'Cl': 4.}[a])+({'Na': 2.5, 'Cl': 4.}[b]))/2
                expected += 2*np.pi/V*na*nb*4*eps*sig**3*((sig/rc)**9/9-(sig/rc)**3/3)
        self.assertAlmostEqual(difference, expected, places=10)
        with self.assertRaisesRegex(ValueError, 'lj_shift=False'):
            bks(tail=True)
        from fys2160_md import FCC
        with self.assertRaisesRegex(ValueError, 'tail=True'):
            MDSimulation(lj_shift=False, tail=True).start(FCC(N=500, rho=1.))


class TiltedCharges(unittest.TestCase):
    def test_pppm_in_a_tilted_box_matches_an_ewald_sum(self):
        rng = np.random.default_rng(5)
        box, tilt = np.array([7., 6.5, 6.2]), np.array([2.1, -1.3, 1.6])
        edges = np.array([[box[0], 0, 0], [tilt[0], box[1], 0], [tilt[1], tilt[2], box[2]]])
        x = (np.indices((2, 2, 3)).reshape(3, -1).T + .5 + rng.normal(scale=.12, size=(12, 3)))/[2, 2, 3]@edges
        q = np.where(np.arange(12) % 2, -1., 1.)
        state = System(x, box, tilt=tilt, species=np.where(q > 0, 'Na', 'Cl'), masses=np.ones(12),
                       velocities=np.zeros((12, 3)))
        sim = MDSimulation(epsilon={'Na': 0., 'Cl': 0.}, sigma={'Na': 0., 'Cl': 0.}, charge={'Na': 1., 'Cl': -1.},
                           coulomb='pppm', kspace_accuracy=1e-7, cutoff=2.4, temperature=1., timestep=.001, ensemble='nve')
        engine = sim._engine(state)
        observed = engine.observables
        energy, forces, _, tensor = ewald(engine._x, q, box, tilt=tilt, alpha=1.4, real_cutoff=6., tensor=True)
        self.assertAlmostEqual(observed['coulomb_energy'], energy, delta=2e-5*abs(energy))
        np.testing.assert_allclose(engine._f, forces, atol=2e-4*np.abs(forces).max())
        V = np.prod(box)
        ours = np.array([observed[f'p{k}'] for k in ('xx', 'yy', 'zz', 'xy', 'xz', 'yz')])*V
        np.testing.assert_allclose(ours, tensor, atol=1e-3*np.abs(tensor).max())

    def test_the_pressure_tensor_is_the_strain_derivative_of_the_energy(self):
        quartz = AlphaQuartz(cells=(4, 2, 3))
        quartz.strain(xy=.04, yz=-.03)
        quartz.positions += np.random.default_rng(2).normal(scale=.02, size=quartz.positions.shape)
        quartz.velocities = np.zeros_like(quartz.positions)
        sim = bks(ensemble='nve', kspace_accuracy=1e-6)
        engine = sim._engine(quartz)
        P = engine.observables
        state, V, h = engine.state(), quartz.box.volume, 1e-4
        for component in ('xx', 'zz', 'xy', 'yz'):
            energy = []
            for sign in (1, -1):
                strained = state.copy()
                strained._pppm = state._pppm  # the same mesh, as in a run
                strained.strain(**{component: sign*h})
                energy.append(sim._engine(strained)._potential)
            derivative = -(energy[0]-energy[1])/(2*h)/V*BAR
            # The PPPM virial holds the influence function fixed: a few bar of tens of thousands.
            self.assertAlmostEqual(P[f'p{component}'], derivative, delta=2e-4*max(abs(derivative), 3e4), msg=component)

    def test_chunks_and_chains_in_a_tilted_box_are_one_trajectory(self):
        quartz = AlphaQuartz(cells=(4, 2, 3))
        quartz.strain(yz=.02)
        sim = bks()
        whole = sim.start(quartz, seed=3)
        whole.advance(40)
        pieces = sim.start(quartz, seed=3)
        for _ in range(4):
            pieces.advance(10)
        np.testing.assert_array_equal(whole.positions, pieces.positions)
        half = sim.start(quartz, seed=3)
        half.advance(20)
        rest = sim.start(half.state)
        rest.advance(20)
        np.testing.assert_array_equal(whole.positions, rest.positions)
        np.testing.assert_array_equal(rest.state.tilt, quartz.tilt)


class Anisotropic(unittest.TestCase):
    def test_each_axis_relaxes_on_its_own(self):
        quartz = AlphaQuartz(cells=(4, 2, 3))  # the measured cell is under tension in BKS
        ratio = quartz.box.lengths[2]/quartz.box.lengths[0]
        iso = bks(ensemble='npt', pressure=0., pressure_time=.2).start(quartz, seed=1)
        iso.advance(200)
        self.assertAlmostEqual(iso.box[2]/iso.box[0], ratio, places=12)  # the shape is kept
        aniso = bks(ensemble='npt', pressure=0., pressure_time=.2, barostat='anisotropic').start(quartz, seed=1)
        aniso.advance(200)
        self.assertGreater(abs(aniso.box[2]/aniso.box[0]-ratio), 1e-3)
        self.assertEqual(aniso.settings['barostat'], 'anisotropic')

    def test_without_a_thermostat_the_extended_enthalpy_is_conserved(self):
        quartz = AlphaQuartz(cells=(4, 2, 3))
        live = bks(ensemble='nph', pressure=0., pressure_time=.3, barostat='anisotropic',
                   kspace_accuracy=1e-6).start(quartz, seed=2)
        values = []
        for _ in range(10):
            live.advance(20)
            values.append(live.observables['extended_enthalpy'])
        drift = np.ptp(values)
        self.assertLess(drift, 0.02*live.observables['kinetic_energy'])
        self.assertGreater(live.observables['barostat_energy'], 0)

    def test_rates_are_kept_in_the_state_and_checkpoint(self):
        quartz = AlphaQuartz(cells=(4, 2, 3))
        sim = bks(ensemble='npt', pressure=0., pressure_time=.2, barostat='anisotropic')
        whole = sim.start(quartz, seed=4)
        whole.advance(30)
        with tempfile.TemporaryDirectory() as folder:
            sim = bks(ensemble='npt', pressure=0., pressure_time=.2, barostat='anisotropic', output_dir=folder)
            first = sim.run(quartz, steps=15, seed=4, thermo_every=5, trajectory_every=15, storage_name='first')
            self.assertEqual(first.state._rate.shape, (3,))
            second = sim.run(Run(first.path).state, steps=15, thermo_every=5, trajectory_every=15, storage_name='second')
            np.testing.assert_array_equal(second.state.positions, whole.positions)
            thermo = second.thermo
            for name in ('pxx', 'pyz', 'buckingham_energy', 'coulomb_energy', 'mass_density'):
                self.assertIn(name, thermo)
            self.assertEqual(second.metadata['configuration']['pairs'][0][:2], ['O', 'O'])

    def test_refusals(self):
        from fys2160_md import FCC
        with self.assertRaisesRegex(ValueError, 'isotropic'):
            MDSimulation(ensemble='npt', pressure=1., barostat='anisotropic').start(FCC(N=500, rho=1.))
        with self.assertRaisesRegex(ValueError, 'barostat must be'):
            bks(barostat='sideways')


if __name__ == '__main__':
    unittest.main()
