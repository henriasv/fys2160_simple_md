"""Rigid SPC/E and TIP4P-type water: geometry, forces, pressure, integrators and runs."""
import tempfile
import unittest
import numpy as np
from fys2160_md import IceIh, RigidWater, MDSimulation, System, Run
from fys2160_md.geometry import _ice_cells

SPCE = RigidWater(oh_length=1.0, hoh_angle=109.47)
TIP4P2005 = RigidWater(oh_length=0.9572, hoh_angle=104.52, m_distance=0.1546)
PARAMETERS = {'spce': (SPCE, 0.15535, 3.166, 0.4238), 'tip4p': (TIP4P2005, 0.1852, 3.1589, 0.5564)}


class WaterTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)

    def simulation(self, model='spce', **kwargs):
        water, epsilon, sigma, q = PARAMETERS[model]
        settings = dict(units='real', epsilon={'O': epsilon, 'H': 0.}, sigma={'O': sigma, 'H': 0.},
                        charge={'O': -2*q, 'H': q}, coulomb='pppm', kspace_accuracy=1e-5, cutoff=6.,
                        temperature=250., timestep=1., ensemble='nve', output_dir=self.temp.name)
        settings.update(kwargs)
        return MDSimulation(**settings)

    def ice(self, model='spce'):
        return IceIh(N=288, density=.93, bond=PARAMETERS[model][0])

    def sim(self, model='spce', system=None, **kwargs):
        """Working engine of one run, for direct checks of forces and dynamics."""
        return self.simulation(model, **kwargs)._engine(system or self.ice(model))

    def test_rigid_water_geometry(self):
        self.assertAlmostEqual(SPCE.hh_length, 2*np.sin(np.radians(109.47/2)), places=14)
        self.assertEqual(SPCE.alpha, 0)
        self.assertAlmostEqual(TIP4P2005.alpha, .1546/(2*.9572*np.cos(np.radians(52.26))), places=14)
        for bad in (dict(oh_length=0), dict(hoh_angle=180), dict(m_distance=-.1), dict(m_distance=.7)):
            with self.assertRaises(ValueError):
                RigidWater(**{**dict(oh_length=1., hoh_angle=109.47), **bad})

    def test_ice_counts_rules_and_dipole(self):
        self.assertEqual(_ice_cells(768)[1], (6, 4, 4))
        self.assertEqual(_ice_cells(2880)[1], (10, 6, 6))
        with self.assertRaisesRegex(ValueError, 'Nearby counts'):
            IceIh(N=700)
        for water in (SPCE, TIP4P2005):
            system = IceIh(N=288, density=.92, bond=water)
            self.assertEqual((system.N, system.units, system.molecule, system.dof), (288, 'real', 'water', 6*288-3))
            atoms = system.atoms
            np.testing.assert_array_equal(atoms.species[:3], ['O', 'H', 'H'])
            np.testing.assert_allclose(atoms.masses[:3], [15.9994, 1.00794, 1.00794])
            self.assertAlmostEqual(atoms.masses.sum()/(.602214076*system.box.volume), .92, places=12)
            x, box = atoms.positions, system.box.lengths
            oxygen, hydrogen = x[0::3], np.concatenate((x[1::3], x[2::3]))
            dr = hydrogen[:, None, :] - oxygen[None, :, :]
            dr -= box*np.rint(dr/box)
            distance = np.linalg.norm(dr, axis=2)
            # Every H is covalently bound to one O and hydrogen-bonded to one other O.
            np.testing.assert_array_equal((distance < 1.1).sum(axis=1), 1)
            np.testing.assert_array_equal((distance < 2.2).sum(axis=1), 2)
            # Every O accepts exactly two hydrogen bonds.
            np.testing.assert_array_equal(((distance > 1.1) & (distance < 2.2)).sum(axis=0), 2)
            bisector = (x[1::3] + x[2::3] - 2*x[0::3])
            bisector -= box*np.rint(bisector/box)
            self.assertLess(np.linalg.norm(bisector.sum(axis=0))/system.N, 1e-6)

    def test_water_system_validation(self):
        good = IceIh(N=288, bond=SPCE)
        x, box, labels = good.atoms.positions.copy(), good.box.lengths, good.atoms.species
        with self.assertRaisesRegex(ValueError, 'RigidWater'):
            System(x, box, molecule='water', bond=TIP4P2005, species=labels, units='real')
        with self.assertRaisesRegex(ValueError, 'O first'):
            System(x, box, molecule='water', bond=SPCE, species=np.roll(labels, 1), units='real')
        v = np.zeros_like(x); v[1] = x[1] - x[0]
        with self.assertRaisesRegex(ValueError, 'perpendicular'):
            System(x, box, velocities=v, molecule='water', bond=SPCE, species=labels, units='real')
        with self.assertRaisesRegex(ValueError, 'units'):
            MDSimulation()._engine(good)
        with self.assertRaisesRegex(ValueError, 'explicit temperature'):
            MDSimulation(units='real', timestep=1., cutoff=6.)

    def test_forces_are_energy_gradients_including_m_site(self):
        for model in ('spce', 'tip4p'):
            sim = self.sim(model, kspace_accuracy=1e-6)
            sim._advance(10)
            x0, forces = sim._x.copy(), sim.forces.copy()
            def energy(x):
                sim._x[:] = x
                sim._advance(0)
                return sim._potential
            h = 1e-4
            for atom in (0, 1, 2, 31):
                for d in range(3):
                    xp, xm = x0.copy(), x0.copy()
                    xp[atom, d] += h; xm[atom, d] -= h
                    # ik-PPPM forces are not exact gradients; errors scale with the accuracy.
                    self.assertAlmostEqual(-(energy(xp)-energy(xm))/(2*h), forces[atom, d], delta=5e-3)
            sim._x[:] = x0

    def test_molecular_pressure_is_volume_derivative(self):
        for model in ('spce', 'tip4p'):
            water = PARAMETERS[model][0]
            sim = self.sim(model, kspace_accuracy=1e-7)
            sim._advance(10)
            state = sim._pppm_state.copy()
            x = sim._x.reshape(-1, 3, 3)
            L, m = sim._box.copy(), sim._mass.reshape(-1, 3)
            d = x - x[:, :1]; d -= L*np.rint(d/L)
            centre = (m[:, :, None]*d).sum(axis=1)/m.sum(axis=1)[:, None]
            def energy(factor):
                y = ((x[:, :1] + centre[:, None, :])*factor + d - centre[:, None, :]).reshape(-1, 3)
                other = self.sim(model, System(y % (L*factor), L*factor, molecule='water', bond=water,
                                               species=sim.system._species, masses=sim._mass, units='real'),
                                 kspace_accuracy=1e-7)
                other._pppm_state[:] = state
                other._advance(0)
                return other._potential
            e, volume = 1e-5, np.prod(L)
            dUdV = (energy(1+e) - energy(1-e))/(volume*((1+e)**3 - (1-e)**3))
            W = sim._terms['molecular_virial']
            self.assertLess(abs(W + 3*volume*dUdV), 2e-3*np.abs(sim._virial))

    def test_constraints_temperature_and_degrees_of_freedom(self):
        sim = self.sim(ensemble='nvt', temperature=300., friction=.01)
        self.assertAlmostEqual(sim.observables['temperature'], 300., places=9)
        sim._advance(200)
        x, L = sim._x, sim._box
        for (a, b), length in zip(((0, 1), (0, 2), (1, 2)), SPCE.constraints):
            dr = x[a::3]-x[b::3]; dr -= L*np.rint(dr/L)
            np.testing.assert_allclose(np.linalg.norm(dr, axis=1), length, atol=1e-11)
            np.testing.assert_allclose(((sim._v[a::3]-sim._v[b::3])*dr).sum(axis=1), 0, atol=1e-12)
        np.testing.assert_allclose((sim._v*sim._mass[:, None]).sum(axis=0), 0, atol=1e-9)

    def test_energy_and_enthalpy_conservation(self):
        sim = self.sim(lj_shift=True)
        energies = []
        for _ in range(5):
            sim._advance(20); energies.append(sim.observables['total_energy'])
        self.assertLess(np.ptp(energies), .5)
        sim = self.sim(lj_shift=True, ensemble='nph', pressure=1.)
        values = []
        for _ in range(5):
            sim._advance(20); o = sim.observables; values.append((o['extended_enthalpy'], o['enthalpy']))
        values = np.array(values)
        self.assertLess(np.ptp(values[:, 0]), .25*np.ptp(values[:, 1]))
        heated = self.sim(lj_shift=True, ensemble='nph', pressure=1., heat_rate=.5)
        start = heated.observables['extended_enthalpy']
        heated._advance(100)
        self.assertAlmostEqual(heated.observables['extended_enthalpy']-start, heated.heat_added, delta=.5)

    def test_deform_expands_box_and_keeps_molecules_rigid(self):
        sim, ice = self.simulation(ensemble='nvt'), self.ice()
        initial = ice.box.lengths
        result = sim.run(ice, steps=40, thermo_every=10, trajectory_every=None, deform=1.2, storage_name='expand')
        state = result.state
        np.testing.assert_allclose(state.box.lengths, 1.2*initial, rtol=1e-13)
        np.testing.assert_array_equal(ice.box.lengths, initial)
        self.assertEqual(result.metadata['deform'], 1.2)
        x, L = state.atoms.positions, state.box.lengths
        for (a, b), length in zip(((0, 1), (0, 2), (1, 2)), SPCE.constraints):
            dr = x[a::3]-x[b::3]; dr -= L*np.rint(dr/L)
            np.testing.assert_allclose(np.linalg.norm(dr, axis=1), length, atol=1e-11)
        self.assertGreater(result.thermo.volume.iloc[-1], 1.7*result.thermo.volume.iloc[0])
        with self.assertRaisesRegex(ValueError, 'nve'):
            sim.run(state, steps=10, ensemble='npt', pressure=1., deform=2.)

    def test_runs_are_independent_of_recording_and_resume_exactly(self):
        sim, ice = self.simulation(ensemble='npt', pressure=1.), self.ice()
        a = sim.run(ice, steps=30, thermo_every=1, trajectory_every=None, storage_name='a')
        b = sim.run(ice, steps=30, thermo_every=7, trajectory_every=10, storage_name='b')
        np.testing.assert_array_equal(a.state.atoms.positions, b.state.atoms.positions)
        np.testing.assert_array_equal(a.state.atoms.velocities, b.state.atoms.velocities)
        whole = sim._engine(ice)
        whole._advance(30)
        self.assertEqual(a.metadata['pppm']['mesh'], [int(n) for n in whole._pppm_state[:3]])
        state = Run(a.path).state
        self.assertEqual((state.units, state.molecule, state.step), ('real', 'water', 30))
        resumed = sim._engine(state)
        whole._advance(20); resumed._advance(20)
        np.testing.assert_array_equal(whole._x, resumed._x)
        np.testing.assert_array_equal(whole._v, resumed._v)
        # A different accuracy is a different numerical choice: the mesh is chosen again.
        coarse = sim.run(state, steps=1, kspace_accuracy=1e-4, trajectory_every=None, storage_name='coarse')
        self.assertLess(np.prod(coarse.metadata['pppm']['mesh']), np.prod(a.metadata['pppm']['mesh']))


if __name__ == '__main__':
    unittest.main()
