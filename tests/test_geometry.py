"""Independent periodic geometry checks and explicit-system runner semantics."""
import tempfile
import unittest
import numpy as np
from fys2160_md import FCC, Random, System, MDSimulation


def distances(system):
    x=system.atoms.positions;box=system.box.lengths
    dr=x[:,None,:]-x[None,:,:];dr-=box*np.rint(dr/box)
    d=np.sqrt(np.sum(dr*dr,axis=-1));np.fill_diagonal(d,np.inf)
    return d


class GeometryTests(unittest.TestCase):
    def test_cubic_fcc_has_exact_density_and_twelve_neighbors(self):
        s=FCC(N=108,rho=.5)
        self.assertEqual(s.N,108)
        self.assertAlmostEqual(s.rho,.5)
        np.testing.assert_allclose(s.box.lengths,s.box.lengths[0])
        nearest=(4/.5)**(1/3)/np.sqrt(2)
        self.assertTrue(np.all(np.isclose(distances(s),nearest).sum(axis=1)==12))
        other=FCC(N=108,box=s.box.lengths[0])
        np.testing.assert_allclose(other.atoms.positions,s.atoms.positions)

    def test_rectangular_fcc_and_incompatible_counts(self):
        s=FCC(N=96,box=[6,9,12])
        np.testing.assert_array_equal(s.box.lengths,[6,9,12])
        self.assertAlmostEqual(s.rho,96/(6*9*12))
        self.assertTrue(np.all(np.isclose(distances(s),3/np.sqrt(2)).sum(axis=1)==12))
        with self.assertRaisesRegex(ValueError,'864 or 1372'):FCC(N=1024,rho=.1)
        with self.assertRaisesRegex(ValueError,'aspect ratio'):FCC(N=108,box=[6,7,8])
        for constructor in (FCC,Random):
            with self.assertRaises(ValueError):constructor(N=108,rho=.1,box=10)
            with self.assertRaises(ValueError):constructor(N=108,box=[10,-1,10])

    def test_random_placement_checks_periodic_separation_and_seed(self):
        a=Random(N=108,box=[8,9,10],min_distance=1.1,seed=1)
        b=Random(N=108,box=[8,9,10],min_distance=1.1,seed=1)
        c=Random(N=108,box=[8,9,10],min_distance=1.1,seed=2)
        self.assertGreaterEqual(distances(a).min(),1.1)
        np.testing.assert_array_equal(a.atoms.positions,b.atoms.positions)
        self.assertFalse(np.array_equal(a.atoms.positions,c.atoms.positions))
        self.assertTrue(np.all(a.atoms.positions>=0))
        self.assertTrue(np.all(a.atoms.positions<a.box.lengths))
        with self.assertRaisesRegex(ValueError,'Could only place'):
            Random(N=32,box=6,min_distance=10,max_attempts=2)

    def test_random_molecules_avoid_other_atoms_and_preserve_bonds(self):
        s=Random(N=32,rho=.03,min_distance=.9,model='diatomic-rigid')
        d=distances(s)
        for i in range(0,64,2):
            self.assertAlmostEqual(d[i,i+1],.7)
            d[i,i+1]=d[i+1,i]=np.inf
        self.assertGreaterEqual(d.min(),.9)
        self.assertEqual(s.N,32);self.assertEqual(len(s.atoms),64)

    def test_runs_take_a_state_and_return_a_new_one(self):
        with tempfile.TemporaryDirectory() as folder:
            system=Random(N=32,box=[8,9,10])
            initial=system.atoms.positions.copy()
            sim=MDSimulation(temperature=1.,output_dir=folder)
            first=sim.run(system,steps=20,storage_name='first')
            hot=sim.run(system,steps=20,temperature=3.,storage_name='hot')  # branch from the same state
            np.testing.assert_array_equal(system.atoms.positions,initial)
            self.assertIsNone(system.atoms.velocities)
            self.assertEqual((system.step,system.time),(0,0.))
            self.assertEqual((first.state.step,first.state.time),(20,20*.005))
            self.assertAlmostEqual(hot.thermo.temperature.iloc[0],3.)
            self.assertEqual(sim.settings['temperature'],1.)  # run arguments do not carry over
            second=sim.run(first.state,steps=10,storage_name='second')
            whole=sim.run(system,steps=30,storage_name='whole')
            self.assertEqual(second.thermo.step.tolist(),[20,30])
            np.testing.assert_array_equal(second.state.atoms.positions,whole.state.atoms.positions)
            np.testing.assert_array_equal(second.state.atoms.velocities,whole.state.atoms.velocities)
            np.testing.assert_array_equal(second.state.box.lengths,[8,9,10])
            with self.assertRaisesRegex(TypeError,'Pass the state to each run'):MDSimulation(system)
            with self.assertRaises(TypeError):sim.run(None)

    def test_geometry_has_no_protocol_and_does_not_require_solver_box_size(self):
        tiny=FCC(N=4,rho=1.)
        self.assertIsNone(tiny.atoms.velocities)
        self.assertFalse(hasattr(tiny,'run'))
        self.assertFalse(hasattr(tiny,'settings'))
        self.assertEqual((tiny.step,tiny.time),(0,0.))
        with self.assertRaises(ValueError):MDSimulation()._engine(tiny)
        MDSimulation(cutoff=.5)._engine(tiny)  # Integration, not geometry, imposes a cutoff.
        for constructor in (FCC,Random):
            with self.assertRaises(TypeError):constructor(temperature=2)
            with self.assertRaises(TypeError):constructor(storage_name='gas')

    def test_explicit_state_is_copied_and_velocities_are_not_reinitialized(self):
        x=np.array([[4.,5,5],[5.2,5,5]])
        v=np.array([[.1,0,0],[-.1,0,0]])
        system=System(x,10,velocities=v)
        x[:]=0;v[:]=0
        sim=MDSimulation(temperature=99,ensemble='nve')._engine(system)
        np.testing.assert_array_equal(sim.atoms.velocities,[[.1,0,0],[-.1,0,0]])
        sim._advance(5)
        np.testing.assert_array_equal(system.atoms.positions,[[4,5,5],[5.2,5,5]])
        np.testing.assert_array_equal(system.atoms.velocities,[[.1,0,0],[-.1,0,0]])
        snapshot=sim.system.atoms.positions
        self.assertFalse(snapshot.flags.writeable)
        snapshot.setflags(write=True);snapshot[:]=0
        self.assertFalse(np.all(sim.system.atoms.positions==0))

    def test_pressure_control_changes_only_the_output_box(self):
        system=FCC(N=108,rho=.05)
        original=system.box.lengths.copy()
        sim=MDSimulation(ensemble='npt',temperature=2,pressure=.2)._engine(system)
        sim._advance(100)
        state=sim.state()
        self.assertFalse(np.array_equal(state.box.lengths,original))
        np.testing.assert_array_equal(system.box.lengths,original)
        self.assertAlmostEqual(state.rho,sim.observables['density'])
