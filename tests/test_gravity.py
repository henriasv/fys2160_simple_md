"""Gravity physics, open boundaries, native buffer safety and saved state."""
import tempfile
import unittest
import numpy as np
from fys2160_md import System, Plummer, FCC, MDSimulation, Run
from fys2160_md import _core


class GravityTests(unittest.TestCase):
    def test_pair_force_energy_and_virial(self):
        x=np.array([[-2.,0,0],[1.,0,0]])
        for a in (0.,.2):
            sim=MDSimulation(pair_potential='gravity',G=.7,softening=a,ensemble='nve',timestep=1e-5)._engine(
                System(x,1.,masses=[2.,3.],velocities=np.zeros_like(x),boundary='open'))
            r=3.;expected=.7*2*3*r/(r*r+a*a)**1.5
            np.testing.assert_allclose(sim.forces,[[expected,0,0],[-expected,0,0]])
            self.assertAlmostEqual(sim._potential,-.7*2*3/np.sqrt(r*r+a*a))
            self.assertAlmostEqual(sim._virial,-expected*r)
            np.testing.assert_array_equal(sim.atoms.positions,x)  # outside viewing frame: no wrapping
            eps=1e-5
            derivative=(-.7*6/np.sqrt((r+eps)**2+a*a)+.7*6/np.sqrt((r-eps)**2+a*a))/(2*eps)
            self.assertAlmostEqual(sim.forces[0,0],derivative,places=9)
            sim._advance(1)
            np.testing.assert_allclose(sim.atoms.velocities[:,0],[expected*1e-5/2,-expected*1e-5/3],rtol=1e-9)

    def test_kepler_orbit_and_time_averaged_negative_response(self):
        # Exact equal-mass circular binary: analytic period, energy and E=-K.
        measurements=[]
        for separation in (1.,2.):
            speed=np.sqrt(1/(2*separation));period=2*np.pi*separation/(2*speed)
            x=np.array([[-separation/2,0,0],[separation/2,0,0]])
            v=np.array([[0,-speed,0],[0,speed,0]])
            sim=MDSimulation(pair_potential='gravity',G=1.,softening=0.,ensemble='nve',
                             timestep=period/4000)._engine(System(x,1.,velocities=v,boundary='open'))
            initial=sim.observables['total_energy'];energies=[];kinetic=[]
            angular=np.cross(x,v).sum(axis=0)
            for _ in range(40):
                sim._advance(100);energies.append(sim.observables['total_energy']);kinetic.append(sim.observables['kinetic_energy'])
            self.assertLess(max(abs(np.array(energies)-initial)),1e-10)
            np.testing.assert_allclose(sim.atoms.positions,x,atol=1e-5)
            np.testing.assert_allclose(np.cross(sim._x,sim._v).sum(axis=0),angular,atol=1e-12)
            self.assertAlmostEqual(np.mean(kinetic),-initial,places=5)
            measurements.append((initial,2*np.mean(kinetic)/sim.dof))
        slope=(measurements[1][0]-measurements[0][0])/(measurements[1][1]-measurements[0][1])
        self.assertAlmostEqual(slope,-1.5,places=4)  # kinetic temperature, not a thermalized binary gas

    def test_open_checkpoint_recording_and_energy_removal(self):
        with tempfile.TemporaryDirectory() as folder:
            system=Plummer(N=32,G=1/32,softening=.1)
            original=system.atoms.positions.copy()
            sim=MDSimulation(pair_potential='gravity',G=1/32,softening=.1,
                             ensemble='nve',timestep=.001,output_dir=folder)
            run=sim.run(system,steps=1000,heat_rate=-.1,thermo_every=50,trajectory_every=200)
            data=run.thermo;balance=data.total_energy-data.heat_added
            self.assertLess(balance.max()-balance.min(),1e-5)
            self.assertAlmostEqual(data.heat_added.iloc[-1],-.1)
            np.testing.assert_array_equal(system.atoms.positions,original)
            self.assertTrue(data.pressure.isna().all())
            self.assertNotIn('compressibility_factor',data)
            state=Run.load(run.path).state
            self.assertEqual(state.boundary,'open');self.assertEqual(run.metadata['configuration']['softening'],.1)
            self.assertIsNone(run.metadata['configuration']['cutoff'])
            whole=sim._engine(system,heat_rate=-.1);whole._advance(1000)
            whole.settings['heat_rate']=0.;whole._advance(100)
            resumed=sim._engine(state);resumed._advance(100)
            np.testing.assert_array_equal(whole._x,resumed._x)
            self.assertIn('"boundary": "open"',run.view(particle_radius=.03).data)
            with self.assertRaises(ValueError):run.view(particle_radius=-1)
            frames=run.load_frames()
            np.testing.assert_array_equal(frames['positions'],frames['unwrapped'])

    def test_plummer_initial_state_and_rejections(self):
        cluster=Plummer(N=100,G=.01,softening=.05)
        twin=Plummer(N=100,G=.01,softening=.05)
        np.testing.assert_array_equal(cluster.atoms.positions,twin.atoms.positions)
        np.testing.assert_allclose(cluster.atoms.velocities.mean(axis=0),0,atol=1e-15)
        gravity=MDSimulation(pair_potential='gravity',G=.01,softening=.05,ensemble='nve')
        self.assertAlmostEqual(gravity._engine(cluster).observables['virial_ratio'],1.)
        with self.assertRaises(ValueError):MDSimulation(pair_potential='gravity')
        with self.assertRaises(ValueError):MDSimulation(pair_potential='gravity',ensemble='nve',cutoff=20)
        with self.assertRaises(ValueError):MDSimulation()._engine(cluster)
        with self.assertRaises(ValueError):gravity._engine(FCC(N=32,rho=.01))
        with self.assertRaises(ValueError):gravity.run(cluster,steps=1,ensemble='npt')
        with self.assertRaises(ValueError):Plummer(N=4097)

    def test_native_validation_and_collision(self):
        x=np.zeros((2,3));v=np.zeros_like(x);f=np.zeros_like(x);mass=np.ones(2)
        with self.assertRaises(ValueError):_core.advance_gravity(x,v,f,mass,0,.001,1.,0.,0.)
        _core.advance_gravity(x,v,f,mass,0,.001,1.,.1,0.)  # finite overlapping softened particles
        with self.assertRaises(ValueError):_core.advance_gravity(x,x,f,mass,0,.001,1.,.1,0.)
        with self.assertRaises((ValueError,TypeError)):_core.advance_gravity(x.astype('float32'),v,f,mass,0,.001,1.,.1,0.)
        with self.assertRaises(ValueError):_core.advance_gravity(x,v,f,mass,0,.001,1.,-.1,0.)
        with self.assertRaises(ValueError):_core.advance_gravity(x,v,f,np.ones(3),0,.001,1.,.1,0.)
        with self.assertRaises(ValueError):_core.advance_gravity(x,v,f,mass,1,.001,1.,.1,-1.)
