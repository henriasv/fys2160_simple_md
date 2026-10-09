"""Frozen short trajectories guarding the existing atomic and diatomic solvers.

Regenerate only for an intentional physics change:
    python tests/test_golden.py --regenerate
"""
import sys
import unittest
from pathlib import Path
import numpy as np
from fys2160_md import FCC, MDSimulation, HarmonicBond, RigidBond, Class2Bond, Plummer

FIXTURE = Path(__file__).with_name('fixtures')/'golden_trajectories.npz'


def cases():
    yield 'atomic_nve', MDSimulation(FCC(N=108, rho=.8), ensemble='nve', temperature=1., timestep=.005, cutoff=2.5), {}
    yield 'atomic_nvt', MDSimulation(FCC(N=108, rho=.8), ensemble='nvt', temperature=1., timestep=.005, cutoff=2.5), {}
    yield 'atomic_npt', MDSimulation(FCC(N=256, rho=.8), ensemble='npt', temperature=1., pressure=1., timestep=.005, cutoff=2.5), {}
    yield 'atomic_nph_heat', MDSimulation(FCC(N=256, rho=.8), ensemble='nph', temperature=1., pressure=1., heat_rate=2., timestep=.005, cutoff=2.5), {}
    yield 'lj96_nvt', MDSimulation(FCC(N=108, rho=.5), pair_potential='lj96', ensemble='nvt', temperature=1., timestep=.005, cutoff=2.5), {}
    yield 'mixture_nvt', MDSimulation(FCC(species={'A': 80, 'B': 28}, rho=.6, masses={'A': 1., 'B': 2.}),
                                      epsilon={'A': 1., 'B': .5}, sigma={'A': 1., 'B': 1.2}, ensemble='nvt', temperature=1.2, timestep=.005, cutoff=2.5), {}
    yield 'harmonic_nvt', MDSimulation(FCC(N=32, rho=.2, molecule='diatomic', bond=HarmonicBond(length=.7, stiffness=100)),
                                       ensemble='nvt', temperature=1., timestep=.002, cutoff=2.5), {}
    yield 'class2_nve', MDSimulation(FCC(N=32, rho=.2, molecule='diatomic', bond=Class2Bond()),
                                     ensemble='nve', temperature=1., timestep=.002, cutoff=2.5), {}
    yield 'rigid_nve', MDSimulation(FCC(N=32, rho=.2, molecule='diatomic', bond=RigidBond(length=.7)),
                                    ensemble='nve', temperature=1., timestep=.002, cutoff=2.5), {}
    yield 'rigid_nvt', MDSimulation(FCC(N=32, rho=.2, molecule='diatomic', bond=RigidBond(length=.7)),
                                    ensemble='nvt', temperature=1., timestep=.002, cutoff=2.5), {}
    yield 'gravity_nve', MDSimulation(Plummer(N=64), pair_potential='gravity', ensemble='nve', temperature=0.,
                                      timestep=.001, cutoff=None, G=1., softening=.05), {}


def trajectory(sim):
    for _ in range(4):
        sim._advance(50)
    obs = sim.observables
    return {'positions': sim._x.copy(), 'velocities': sim._v.copy(), 'box': sim._box.copy(),
            'energies': np.array([obs['potential_energy'], obs['kinetic_energy'], obs['pressure']])}


class GoldenTrajectoryTests(unittest.TestCase):
    def test_existing_solvers_reproduce_frozen_trajectories(self):
        with np.load(FIXTURE) as data:
            for name, sim, _ in cases():
                with self.subTest(name):
                    for key, value in trajectory(sim).items():
                        np.testing.assert_allclose(value, data[f'{name}/{key}'], rtol=1e-11, atol=1e-12,
                                                   equal_nan=True, err_msg=f'{name} {key}')


if __name__ == '__main__' and '--regenerate' in sys.argv:
    FIXTURE.parent.mkdir(exist_ok=True)
    np.savez_compressed(FIXTURE, **{f'{name}/{key}': value for name, sim, _ in cases()
                                    for key, value in trajectory(sim).items()})
    print('Wrote', FIXTURE)
elif __name__ == '__main__':
    unittest.main()
