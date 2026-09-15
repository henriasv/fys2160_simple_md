"""Independent analytical and numerical checks of the native solver and persistence."""
import tempfile
import unittest
from pathlib import Path
import numpy as np
from fys2160_md import Simulation, System, FCC, Run
from fys2160_md import _core


class SolverTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = self.temp.name

    def sim(self, **kwargs):
        return Simulation(FCC(N=kwargs.pop('N',500),rho=kwargs.pop('rho',.001),
                              model=kwargs.pop('model','atomic')),output_dir=self.root,**kwargs)

    def test_pair_force_energy_and_periodicity(self):
        for separation in (1.1, 1.8, 2.49, 2.6):
            s = Simulation(System.from_arrays([[0.2, 1, 1], [10.2 - separation, 1, 1]], np.zeros((2, 3)), 10), ensemble='nve', output_dir=self.root)
            r=separation
            expected=24*(2*r**-13-r**-7) if r<2.5 else 0
            energy=4*(r**-12-r**-6)-4*(2.5**-12-2.5**-6) if r<2.5 else 0
            np.testing.assert_allclose(s.forces[0], [expected,0,0], atol=2e-11)
            np.testing.assert_allclose(s.forces.sum(axis=0), 0, atol=1e-12)
            self.assertAlmostEqual(s.observables['potential_energy'],energy,places=11)
            self.assertAlmostEqual(s._virial, r*expected, places=10)

    def test_molecular_force_matches_independent_potential_gradient(self):
        x=np.array([[4.,4,4],[4.8,4,4],[6.,4.2,4],[6.65,4.1,4]])
        def potential(p):
            U=0.
            for a in (0,2):
                q=np.linalg.norm(p[a]-p[a+1])-.7
                U+=35*q*q-40*q**3+59*q**4
            for i in (0,1):
                for j in (2,3):
                    r=np.linalg.norm(p[i]-p[j])
                    if r<2.5:U+=2*r**-9-3*r**-6-(2*2.5**-9-3*2.5**-6)
            return U
        s=Simulation(System.from_arrays(x, np.zeros_like(x), 12, model='diatomic-flexible'), pair_potential='lj96', ensemble='nve', output_dir=self.root)
        gradient=np.zeros_like(x)
        for i in range(4):
            for d in range(3):
                a=x.copy();b=x.copy();a[i,d]+=1e-6;b[i,d]-=1e-6
                gradient[i,d]=-(potential(a)-potential(b))/2e-6
        np.testing.assert_allclose(s.forces,gradient,rtol=2e-8,atol=1e-8)
        self.assertAlmostEqual(s.observables['potential_energy'],potential(x),places=11)

    def test_neighbor_list_matches_brute_force_after_motion(self):
        s=self.sim(N=108,rho=.25)
        for _ in range(8):
            s._advance(200)
            x=s.atoms.positions; f=np.zeros_like(x);u=0.
            for i in range(len(x)):
                for j in range(i+1,len(x)):
                    dr=x[i]-x[j];dr-=s._box*np.rint(dr/s._box);r2=dr@dr
                    if r2<2.5**2:
                        inv6=r2**-3; force=24*inv6*(2*inv6-1)/r2*dr
                        f[i]+=force;f[j]-=force
                        u+=4*inv6*(inv6-1)-4*(2.5**-12-2.5**-6)
            np.testing.assert_allclose(s.forces,f,atol=5e-11,rtol=1e-12)
            self.assertAlmostEqual(s._potential,u,places=10)

    def test_verlet_energy_converges_with_timestep(self):
        errors=[]
        for dt in (.004,.002):
            s=Simulation(System.from_arrays([[4, 5, 5], [5.15, 5, 5]], [[0.1, 0, 0], [-0.1, 0, 0]], 10), ensemble='nve', timestep=dt, output_dir=self.root)
            initial=s.observables['total_energy'];e=[]
            for _ in range(200):
                s._advance(round(.02/dt));e.append(s.observables['total_energy'])
            errors.append(max(abs(np.array(e)-initial)))
        self.assertLess(errors[1],errors[0]*.3)
        self.assertLess(errors[1],1e-4)

    def test_rigid_bonds_rotation_and_heat(self):
        s=self.sim(model='diatomic-rigid',N=32,rho=1e-8,ensemble='nve')
        initial=s.observables['total_energy']
        s._advance(10000)
        np.testing.assert_allclose(np.linalg.norm(s._bond_vectors(),axis=1),.7,atol=1e-12)
        np.testing.assert_allclose((s._bond_vectors()*(s._v[::2]-s._v[1::2])).sum(axis=1),0,atol=1e-12)
        self.assertAlmostEqual(s.observables['total_energy'],initial,places=7)
        s._settings['heat_rate']=2.
        s._advance(1000)
        self.assertAlmostEqual(s.observables['total_energy']-initial,4.,places=7)

    def test_heat_is_energy_per_time_and_atomic_ideal_cv(self):
        s=self.sim(N=32,rho=1e-7,ensemble='nve',heat_rate=5.)
        a=s.observables;s._advance(2000);b=s.observables
        self.assertAlmostEqual(b['total_energy']-a['total_energy'],50.,places=8)
        self.assertAlmostEqual((b['total_energy']-a['total_energy'])/(b['temperature']-a['temperature'])/32,1.5*(1-1/32),places=10)
        self.assertAlmostEqual(b['compressibility_factor'],1-1/32,places=12)

    def test_langevin_temperature_mean_and_fluctuations(self):
        for model in ('atomic','diatomic-flexible','diatomic-rigid'):
            s=self.sim(model=model,N=108,rho=.001)
            s._advance(5000)
            T=[]
            for i in range(1000):
                s._advance(100);T.append(s.observables['temperature'])
            self.assertLess(abs(np.mean(T)-2),.07,model)
            self.assertLess(abs(np.var(T)/(8/s.dof)-1),.25,model)
            np.testing.assert_allclose((s._v*s._mass[:,None]).sum(axis=0),0,atol=2e-11)

    def test_barostat_extended_enthalpy_and_pressure(self):
        s=self.sim(N=108,rho=.005,ensemble='npt')
        s._advance(15000)
        s._settings['ensemble']='nph';e=[];p=[]
        for i in range(300):
            s._advance(100);o=s.observables;e.append(o['extended_enthalpy']);p.append(o['pressure'])
        self.assertLess(np.ptp(e)/np.mean(e),.001)
        self.assertLess(abs(np.mean(p)-.01),.0005)

    def test_ideal_npt_volume_statistics(self):
        s=self.sim(N=108,rho=1e-7,pressure=2e-7,ensemble='npt')
        s._advance(100000)
        volumes=[];pressures=[]
        for _ in range(10000):
            s._advance(100)
            volumes.append(s.observables['volume'])
            pressures.append(s.observables['pressure'])
        v=np.array(volumes)
        # Fixed scaled COM and zero momentum: gamma(N) volume distribution.
        self.assertLess(abs(v.mean()/(108*2/2e-7)-1),.03)
        self.assertLess(abs(v.var()/v.mean()**2*108-1),.25)
        self.assertLess(abs(np.mean(pressures)/2e-7-1),.03)

    def test_saved_run_reload_resume_and_sampling_independence(self):
        s=self.sim(N=32,model='diatomic-rigid')
        a=s.run(storage_name="first",steps=251,thermo_every=43,trajectory_every=97)
        self.assertEqual(a.metadata['completed_steps'],251)
        self.assertEqual(a.thermo.step.tolist(),[0,43,86,129,172,215,251])
        self.assertEqual([int(f['step']) for f in a.iter_frames()],[0,97,194,251])
        resumed=Simulation.from_run(a)
        s.run(storage_name="continued",steps=333,thermo_every=100,trajectory_every=None)
        resumed.run(storage_name="restarted",steps=333,thermo_every=17,trajectory_every=71)
        np.testing.assert_allclose(s.atoms.positions,resumed.atoms.positions,atol=1e-12)
        np.testing.assert_allclose(s.atoms.velocities,resumed.atoms.velocities,atol=1e-12)
        self.assertEqual(len(Run.find(self.root)),3)
        self.assertEqual(Run.load(a.path).thermo.shape,a.thermo.shape)
        self.assertEqual(len(a.load_frames()['positions']),4)

    def test_interrupt_and_numerical_failure_are_recorded(self):
        from unittest.mock import patch
        import json
        s=self.sim(N=32)
        advance=_core.advance
        def interrupted_advance(*args):
            result=list(advance(*args))
            if args[6]:result[4]=1
            return tuple(result)
        with patch.object(_core, 'advance', interrupted_advance):
            run=s.run(steps=1000,thermo_every=100)
        self.assertEqual(run.metadata['status'],'interrupted')
        self.assertEqual(run.metadata['completed_steps'],100)
        self.assertEqual(int(run.thermo.step.iloc[-1]),100)
        restart=Simulation.from_run(run)
        np.testing.assert_array_equal(restart.atoms.velocities,s.atoms.velocities)
        broken=self.sim(model='diatomic-rigid',N=32)
        with self.assertRaises(ValueError):
            broken.run(timestep=10,steps=100,storage_name='unstable')
        failed=[r for r in Run.find(self.root) if r.metadata['label']=='unstable'][0]
        self.assertEqual(failed.metadata['status'],'failed')
        self.assertFalse((failed.path/'checkpoint.npz').exists())
        with self.assertRaises(RuntimeError):broken.run(steps=10)
        with self.assertRaises(ValueError):Simulation.from_run(failed)

    def test_input_errors_and_readonly_views(self):
        for kwargs in (dict(N=1),dict(rho=-1),dict(temperature=-1),dict(ensemble='foo'),dict(heat_rate=5),dict(model='diatomic-rigid',ensemble='npt')):
            with self.assertRaises(ValueError):self.sim(**kwargs)
        s=self.sim(N=32)
        with self.assertRaises(ValueError):s.atoms.positions[0,0]=0
        with self.assertRaises(ValueError):s.run(steps=1.5)
        with self.assertRaises(ValueError):s.run(heat_rate=1)
        with self.assertRaises(ValueError):
            _core.advance(s._x,s._v,s._x,s._f,s._box,s._mass,0,.005,2.5,.4,2.,1.,1,0,0.,-1.,1.,0.)
        self.assertFalse(list(Path(self.root).iterdir()))


if __name__=='__main__':unittest.main()
