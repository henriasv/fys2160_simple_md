"""Independent periodic geometry checks and explicit-system runner semantics."""
import tempfile
import unittest
import numpy as np
from fys2160_md import FCC, Random, System, Simulation


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

    def test_simulations_copy_the_same_starting_system_and_resume(self):
        with tempfile.TemporaryDirectory() as folder:
            system=Random(N=32,box=[8,9,10])
            initial=system.atoms.positions.copy()
            a=Simulation(system,temperature=1.,output_dir=folder)
            b=Simulation(system,temperature=3.,output_dir=folder)
            untouched=b.atoms.positions.copy()
            run=a.run(steps=20,storage_name='random')
            a.run(steps=10,storage_name='continued')
            self.assertEqual(a.step,30);self.assertEqual(b.step,0)
            np.testing.assert_array_equal(b.atoms.positions,untouched)
            np.testing.assert_array_equal(system.atoms.positions,initial)
            self.assertIsNone(system.atoms.velocities)
            self.assertIsNot(a.system,system)
            np.testing.assert_array_equal(a.system.atoms.positions,a.atoms.positions)
            self.assertAlmostEqual(b.observables['temperature'],3.)
            resumed=Simulation.from_run(run)
            self.assertEqual(resumed.step,20)
            np.testing.assert_array_equal(resumed.system.box.lengths,[8,9,10])
            resumed.run(steps=10,storage_name='resumed')
            np.testing.assert_array_equal(resumed.system.atoms.positions,a.system.atoms.positions)
            np.testing.assert_array_equal(resumed.system.atoms.velocities,a.system.atoms.velocities)
            with self.assertRaises(TypeError):Simulation(None)

    def test_geometry_has_no_protocol_and_does_not_require_solver_box_size(self):
        tiny=FCC(N=4,rho=1.)
        self.assertIsNone(tiny.atoms.velocities)
        self.assertFalse(hasattr(tiny,'run'))
        self.assertFalse(hasattr(tiny,'settings'))
        self.assertFalse(hasattr(tiny,'step'))
        with self.assertRaises(ValueError):Simulation(tiny)
        Simulation(tiny,cutoff=.5)  # Integration, not geometry, imposes a cutoff.
        for constructor in (FCC,Random):
            with self.assertRaises(TypeError):constructor(temperature=2)
            with self.assertRaises(TypeError):constructor(storage_name='gas')

    def test_explicit_state_is_copied_and_velocities_are_not_reinitialized(self):
        x=np.array([[4.,5,5],[5.2,5,5]])
        v=np.array([[.1,0,0],[-.1,0,0]])
        system=System(x,10,velocities=v)
        x[:]=0;v[:]=0
        sim=Simulation(system,temperature=99,ensemble='nve')
        np.testing.assert_array_equal(sim.atoms.velocities,[[.1,0,0],[-.1,0,0]])
        sim._advance(5)
        np.testing.assert_array_equal(system.atoms.positions,[[4,5,5],[5.2,5,5]])
        np.testing.assert_array_equal(system.atoms.velocities,[[.1,0,0],[-.1,0,0]])
        snapshot=sim.system.atoms.positions
        self.assertFalse(snapshot.flags.writeable)
        snapshot.setflags(write=True);snapshot[:]=0
        self.assertFalse(np.all(sim.system.atoms.positions==0))

    def test_pressure_control_updates_owned_system_box(self):
        system=FCC(N=108,rho=.05)
        original=system.box.lengths.copy()
        sim=Simulation(system,ensemble='npt',temperature=2,pressure=.2)
        sim._advance(100)
        self.assertFalse(np.array_equal(sim.system.box.lengths,original))
        np.testing.assert_array_equal(system.box.lengths,original)
        self.assertAlmostEqual(sim.system.rho,sim.observables['density'])
