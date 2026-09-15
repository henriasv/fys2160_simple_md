"""Mixture forces, mass-dependent dynamics, geometry and persistence."""
import tempfile
import unittest
import numpy as np
from fys2160_md import FCC, Random, System, MDSimulation, RigidBond


class SpeciesTests(unittest.TestCase):
    def test_composition_and_mass_assignment(self):
        for maker in (FCC,Random):
            system=maker(species={'light':24,'heavy':8},rho=.01,masses={'light':1,'heavy':4})
            self.assertEqual(system.N,32)
            labels=system.atoms.species
            self.assertEqual(sum(labels=='light'),24);self.assertEqual(sum(labels=='heavy'),8)
            np.testing.assert_array_equal(system.atoms.masses,np.where(labels=='light',1,4))
            self.assertIsNone(system.atoms.velocities)
            np.testing.assert_array_equal(system.copy().atoms.species,labels)
            with self.assertRaises(ValueError):maker(N=108,species={'A':24,'B':8})
            with self.assertRaises(ValueError):maker(species={'A':24,'B':8},masses={'A':1})
            with self.assertRaises(ValueError):maker(species={'A':24,'B':8},masses={'A':1,'B':-1})
            with self.assertRaises(ValueError):maker(species={'A':0})

    def test_unlike_pair_force_mixing_and_mass_acceleration(self):
        for potential in ('lj','lj96'):
            x=np.array([[3.,3,3],[4.6,3,3]])
            system=System(x,10,velocities=np.zeros_like(x),species=['A','B'],masses=[1,4])
            sim=MDSimulation(system,epsilon={'A':.5,'B':2},sigma={'A':.8,'B':1.2},
                             pair_potential=potential,ensemble='nve',timestep=1e-5)
            r=1.6;eps=1.;sig=1.
            def U(d):return eps*(4*((sig/d)**12-(sig/d)**6) if potential=='lj' else 2*(sig/d)**9-3*(sig/d)**6)
            f=-(U(r+1e-6)-U(r-1e-6))/2e-6
            self.assertAlmostEqual(sim.forces[1,0],f,places=8)
            self.assertAlmostEqual(sim._potential,U(r)-U(2.5),places=11)
            sim._advance(1)
            np.testing.assert_allclose(sim.atoms.velocities[:,0],[-f*1e-5,f*1e-5/4],rtol=1e-8)
            np.testing.assert_allclose((sim.atoms.velocities*sim.atoms.masses[:,None]).sum(axis=0),0,atol=1e-14)

    def test_species_thermostat_equipartition_and_restart(self):
        system=Random(species={'A':40,'B':40},rho=1e-8,masses={'A':1,'B':4})
        with tempfile.TemporaryDirectory() as folder:
            sim=MDSimulation(system,epsilon={'A':1,'B':.6},sigma={'A':1,'B':1.1},output_dir=folder)
            sim._advance(5000)
            samples={name:[] for name in system.species}
            for _ in range(800):
                sim._advance(100)
                for name in samples:
                    mask=sim.atoms.species==name
                    samples[name].append(np.mean(np.sum(sim._v[mask]**2,axis=1)))
            for name,mass in [('A',1),('B',4)]:
                expected=3*2/mass*(1-mass/200)  # zero total momentum constraint
                self.assertLess(abs(np.mean(samples[name])/expected-1),.04)
            run=sim.run(steps=100)
            resumed=MDSimulation.from_run(run)
            self.assertEqual(resumed.epsilon,sim.epsilon);self.assertEqual(resumed.sigma,sim.sigma)
            np.testing.assert_array_equal(resumed.atoms.species,sim.atoms.species)
            np.testing.assert_array_equal(resumed.atoms.masses,sim.atoms.masses)
            sim._advance(100);resumed._advance(100)
            np.testing.assert_array_equal(sim.atoms.positions,resumed.atoms.positions)
            self.assertIn('"species"',run.view().data)
            parameters=sim.epsilon;parameters['A']=999
            self.assertEqual(sim.epsilon['A'],1)

    def test_molecular_mixtures_and_unequal_bonded_masses(self):
        system=Random(species={'A':10,'B':6},molecule='diatomic',bond=RigidBond(length=.9),rho=.001,
                      masses={'A':1,'B':3})
        self.assertEqual(system.N,16);self.assertEqual(len(system.atoms),32)
        np.testing.assert_array_equal(system.atoms.species[::2],system.atoms.species[1::2])
        # Explicit arrays also support different masses within a rigid molecule.
        x=np.array([[3.,3,3],[3.9,3,3],[13.,13,13],[13.9,13,13]])
        sim=MDSimulation(System(x,30,molecule='diatomic',bond=RigidBond(length=.9),
                               masses=[1,3,1,3],species=['A','B','A','B']),ensemble='nve')
        energy=sim.observables['total_energy'];sim._advance(1000)
        np.testing.assert_allclose(np.linalg.norm(sim._bond_vectors(),axis=1),.9,atol=1e-12)
        np.testing.assert_allclose((sim._v*sim._mass[:,None]).sum(axis=0),0,atol=1e-12)
        self.assertAlmostEqual(sim.observables['total_energy'],energy,places=9)

    def test_parameter_maps_are_complete(self):
        system=Random(species={'A':10,'B':6},rho=.001)
        for kw in ({'epsilon':{'A':1}},{'sigma':{'A':1,'B':1,'C':1}},{'mixing_rule':'unknown'}):
            with self.assertRaises(ValueError):MDSimulation(system,**kw)
