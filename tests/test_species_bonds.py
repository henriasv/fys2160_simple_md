"""Species-dependent bond geometry, forces, constraints and output names."""
import tempfile
import unittest
import numpy as np
from fys2160_md import FCC,Random,System,MDSimulation,HarmonicBond,RigidBond,Class2Bond


class SpeciesBondTests(unittest.TestCase):
    def test_species_lengths_in_fcc_and_periodic_random_geometry(self):
        for maker in (FCC,Random):
            settings={'N2':RigidBond(length=.6),'O2':RigidBond(length=.9)}
            system=maker(species={'N2':24,'O2':8},rho=.03,molecule='diatomic',bond=settings)
            expected=np.where(system.atoms.species[::2]=='N2',.6,.9)
            dr=system.atoms.positions[::2]-system.atoms.positions[1::2]
            dr-=system.box.lengths*np.rint(dr/system.box.lengths)
            np.testing.assert_allclose(np.linalg.norm(dr,axis=1),expected,atol=1e-14)
            settings['N2']=RigidBond(length=1.1)
            copy=system.bond;copy['O2']=RigidBond(length=1.2)
            self.assertEqual(system.bond['N2'].length,.6)
            self.assertEqual(system.bond['O2'].length,.9)
            sim=MDSimulation(system,ensemble='nve',temperature=.3)
            for _ in range(10):sim._advance(100)
            np.testing.assert_allclose(np.linalg.norm(sim._bond_vectors(),axis=1),expected,atol=1e-12)
            np.testing.assert_allclose(np.sum(sim._bond_vectors()*(sim._v[::2]-sim._v[1::2]),axis=1),0,atol=1e-12)

    def test_independent_species_spring_energies_forces_and_periods(self):
        x=np.array([[3.,3,3],[3.8,3,3],[13.,13,13],[14.,13,13]])
        bonds={'A':HarmonicBond(length=.7,stiffness=40),'B':HarmonicBond(length=.9,stiffness=80)}
        sim=MDSimulation(System(x,30,velocities=np.zeros_like(x),masses=[1,1,2,2],
            species=['A','A','B','B'],molecule='diatomic',bond=bonds),ensemble='nve',timestep=.0002)
        np.testing.assert_allclose(sim.forces[:,0],[4,-4,8,-8],atol=1e-12)
        self.assertAlmostEqual(sim._potential,.6,places=12)
        initial=sim.observables['total_energy'];sim._advance(1000)
        # Frequencies sqrt(k / reduced_mass) both equal sqrt(80).
        expected=np.array([.7,.9])+.1*np.cos(np.sqrt(80)*.2)
        np.testing.assert_allclose(np.linalg.norm(sim._bond_vectors(),axis=1),expected,atol=1e-7)
        self.assertLess(abs(sim.observables['total_energy']-initial),6e-7)

    def test_bond_maps_survive_checkpoint_and_recording_intervals(self):
        with tempfile.TemporaryDirectory() as folder:
            for bonds in ({'A':RigidBond(.6),'B':RigidBond(.9)},
                          {'A':HarmonicBond(.7,60),'B':Class2Bond(.8,40,-30,50)}):
                sim=MDSimulation(Random(species={'A':12,'B':4},rho=.005,molecule='diatomic',bond=bonds),output_dir=folder)
                result=sim.run(steps=111,thermo_every=17,trajectory_every=31)
                self.assertEqual(result.thermo.step.tolist(),[0,17,34,51,68,85,102,111])
                self.assertEqual([int(f['step']) for f in result.iter_frames()],[0,31,62,93,111])
                resumed=MDSimulation.from_run(result)
                self.assertEqual(resumed.system.bond,bonds)
                sim._advance(100);resumed._advance(100)
                np.testing.assert_array_equal(sim.atoms.positions,resumed.atoms.positions)
                np.testing.assert_array_equal(sim.atoms.velocities,resumed.atoms.velocities)
            legacy=sim.run(steps=10,sample_every=2,save_every=None)
            self.assertEqual(legacy.metadata['thermo_every'],2)
            self.assertIsNone(legacy.metadata['trajectory_every'])
            self.assertFalse(list((legacy.path/'frames').glob('*.npz')))
            with self.assertRaises(TypeError):sim.run(thermo_every=100,sample_every=5)
            with self.assertRaises(TypeError):sim.run(trajectory_every=500,save_every=10)

    def test_mapping_errors_are_explicit(self):
        for bonds in ({'A':RigidBond(.7)}, {'A':RigidBond(.7),'B':HarmonicBond()},
                      {'A':RigidBond(.7),'B':RigidBond(.8),'C':RigidBond(.9)}):
            with self.assertRaises(ValueError):Random(species={'A':4,'B':4},molecule='diatomic',bond=bonds)
        with self.assertRaises(TypeError):Random(species={'A':4},molecule='diatomic',bond={'A':.7})
        with self.assertRaises(ValueError):Random(species={'A':4},box=3,molecule='diatomic',bond={'A':RigidBond(2)})
        with self.assertRaises(ValueError):
            System([[0,0,0],[.7,0,0],[3,3,3],[3.7,3,3]],10,molecule='diatomic',
                   species=['A','B','A','B'],bond={'A':RigidBond(.7),'B':RigidBond(.7)})
