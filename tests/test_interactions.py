"""Independent checks of configurable pair and bond physics."""
import json
import tempfile
import unittest
import numpy as np
from fys2160_md import MDSimulation, System, Random, FCC, HarmonicBond, RigidBond, Class2Bond


class InteractionTests(unittest.TestCase):
    def test_lj_parameters_force_energy_and_pressure(self):
        for potential in ('lj', 'lj96'):
            for epsilon, sigma in ((1., 1.), (.4, 1.3), (2., .8)):
                r = 1.6
                sim = MDSimulation(System([[3, 3, 3], [3+r, 3, 3]], 12,
                    velocities=np.zeros((2, 3))), pair_potential=potential,
                    epsilon=epsilon, sigma=sigma)
                def U(d):
                    z = sigma/d
                    return epsilon*(4*(z**12-z**6) if potential == 'lj' else 2*z**9-3*z**6)
                force = -(U(r+1e-6)-U(r-1e-6))/2e-6
                self.assertAlmostEqual(sim.forces[1, 0], force, places=8)
                self.assertAlmostEqual(sim.observables['potential_energy'], U(r)-U(2.5), places=11)
                self.assertAlmostEqual(sim.observables['pressure'], r*force/(3*12**3), places=10)

    def test_molecular_gradient_exclusions_and_com_pressure(self):
        x = np.array([[3., 3, 3], [4.1, 3.1, 3], [5.4, 3.2, 3], [6.5, 3.3, 3]])
        bond = HarmonicBond(length=1.05, stiffness=70)
        for potential in ('lj', 'lj96'):
            for exclude in (True, False):
                def U(p):
                    value = 0.
                    for i, j in ((0, 1), (2, 3)):
                        value += .5*bond.stiffness*(np.linalg.norm(p[i]-p[j])-bond.length)**2
                    for i in range(4):
                        for j in range(i+1, 4):
                            if exclude and i//2 == j//2: continue
                            r = np.linalg.norm(p[i]-p[j])
                            if r >= 2.5: continue
                            def pair(d):
                                z=.9/d
                                return .6*(4*(z**12-z**6) if potential=='lj' else 2*z**9-3*z**6)
                            value += pair(r)-pair(2.5)
                    return value
                sim = MDSimulation(System(x, 12, velocities=np.zeros_like(x), molecule='diatomic', bond=bond),
                    pair_potential=potential, exclude_bonded_pairs=exclude, epsilon=.6, sigma=.9)
                gradient=np.zeros_like(x)
                for i in range(4):
                    for d in range(3):
                        a=x.copy(); b=x.copy(); a[i,d]+=1e-6; b[i,d]-=1e-6
                        gradient[i,d]=-(U(a)-U(b))/2e-6
                np.testing.assert_allclose(sim.forces, gradient, rtol=2e-8, atol=1e-8)
                self.assertAlmostEqual(sim.observables['potential_energy'], U(x), places=11)
                # Scale molecule centres while holding internal coordinates fixed.
                centres=np.repeat(x.reshape(2,2,3).mean(axis=1),2,axis=0)
                virial=-(U(x+1e-6*centres)-U(x-1e-6*centres))/2e-6
                self.assertAlmostEqual(sim.observables['pressure'], virial/(3*12**3), places=9)

    def test_harmonic_oscillator_matches_analytic_motion(self):
        bond=HarmonicBond(length=.8,stiffness=40)
        x=np.array([[4.,4,4],[4.9,4,4],[14.,14,14],[14.8,14,14]])
        sim=MDSimulation(System(x,30,velocities=np.zeros_like(x),molecule='diatomic',bond=bond),
                         ensemble='nve',timestep=.0002)
        initial=sim.observables['total_energy']
        sim._advance(1000)
        separation=np.linalg.norm(sim._bond_vectors()[0])
        expected=.8+.1*np.cos(np.sqrt(80)*.2)  # reduced mass = 1/2
        self.assertAlmostEqual(separation,expected,places=7)
        self.assertLess(abs(sim.observables['total_energy']-initial),2e-7)

    def test_custom_rigid_length_and_saved_physics(self):
        with tempfile.TemporaryDirectory() as folder:
            for bond in (HarmonicBond(length=.85,stiffness=60), RigidBond(length=.85), Class2Bond()):
                sim=MDSimulation(Random(N=16,rho=.001,molecule='diatomic',bond=bond),
                                 pair_potential='lj96',epsilon=.3,sigma=.8,output_dir=folder)
                run=sim.run(steps=100,storage_name=type(bond).__name__)
                resumed=MDSimulation.from_run(run)
                self.assertEqual(resumed.system.bond,bond)
                self.assertEqual(resumed.pair_potential,'lj96')
                self.assertEqual(resumed.epsilon,.3);self.assertEqual(resumed.sigma,.8)
                self.assertTrue(resumed.exclude_bonded_pairs)
                sim._advance(200);resumed._advance(200)
                np.testing.assert_array_equal(sim.atoms.positions,resumed.atoms.positions)
                np.testing.assert_array_equal(sim.atoms.velocities,resumed.atoms.velocities)
                if isinstance(bond,RigidBond):
                    np.testing.assert_allclose(np.linalg.norm(sim._bond_vectors(),axis=1),.85,atol=1e-12)
                    np.testing.assert_allclose(np.sum(sim._bond_vectors()*(sim._v[::2]-sim._v[1::2]),axis=1),0,atol=1e-12)

    def test_legacy_checkpoint_preserves_old_force_field(self):
        with tempfile.TemporaryDirectory() as folder:
            sim=MDSimulation(FCC(N=32,rho=.001,model='diatomic-flexible'),pair_potential='lj96',output_dir=folder)
            run=sim.run(steps=100)
            path=run.path/'metadata.json';metadata=json.loads(path.read_text())
            for name in ('bond','pair_potential','exclude_bonded_pairs','epsilon','sigma'):
                metadata['configuration'].pop(name)
            path.write_text(json.dumps(metadata))
            resumed=MDSimulation.from_run(run.path)
            self.assertEqual(resumed.system.bond,Class2Bond())
            self.assertEqual(resumed.pair_potential,'lj96')
            sim._advance(100);resumed._advance(100)
            np.testing.assert_array_equal(sim.atoms.positions,resumed.atoms.positions)

    def test_validation_and_geometry(self):
        for constructor in (Random,FCC):
            s=constructor(N=32,rho=.01,molecule='diatomic',bond=RigidBond(length=1.1))
            self.assertEqual(s.N,32);self.assertEqual(len(s.atoms),64)
            np.testing.assert_array_equal(s.bonds,np.arange(64).reshape(-1,2))
            dr=s.atoms.positions[::2]-s.atoms.positions[1::2]
            dr-=s.box.lengths*np.rint(dr/s.box.lengths)
            np.testing.assert_allclose(np.linalg.norm(dr,axis=1),1.1)
            with self.assertRaises(ValueError):constructor(bond=HarmonicBond())
            with self.assertRaises(ValueError):constructor(molecule='diatomic',model='atomic')
        system=FCC(N=32,rho=.01)
        for kw in ({'epsilon':0},{'sigma':-1},{'sigma':float('nan')},{'pair_potential':'wrong'}):
            with self.assertRaises(ValueError):MDSimulation(system,**kw)
        with self.assertRaises(TypeError):MDSimulation(system,exclude_bonded_pairs='yes')
        with self.assertRaises(ValueError):HarmonicBond(stiffness=-1)
        with self.assertRaises(ValueError):RigidBond(length=0)
