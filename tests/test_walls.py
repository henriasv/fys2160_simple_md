"""Flat (two-dimensional) systems, walls with the force on them, and pistons."""
import tempfile
import unittest
import numpy as np
from fys2160_md import FCC, MDSimulation, Piston, Random, Run, System, _core

WALLS_2D = ('left', 'right', 'bottom', 'top')


def sim(**options):
    settings = dict(pair_potential='lj', ensemble='nvt', temperature=2.0, timestep=0.005, cutoff=2.5)
    settings.update(options)
    return MDSimulation(**settings)


def gas(N=100, box=(40, 40), boundary='walls', seed=5):
    return Random(N=N, box=box, boundary=boundary, min_distance=1.0, seed=seed)


def settled(state=None, steps=2000, **options):
    """A state after a short run with the thermostat."""
    live = sim(**options).start(gas() if state is None else state)
    live.advance(steps)
    return live.state


class Surface(unittest.TestCase):
    def test_a_box_with_two_lengths_is_flat(self):
        state = gas(N=50, box=(20, 30))
        self.assertEqual(state.dimension, 2)
        self.assertEqual(state.atoms.positions.shape, (50, 2))
        self.assertEqual(state.box.lengths.tolist(), [20., 30.])
        self.assertEqual(state.box.volume, 600.)
        self.assertAlmostEqual(state.rho, 50/600)
        self.assertEqual(state.boundary, 'walls')
        self.assertEqual(state.walls, WALLS_2D)
        self.assertIn('dimension=2', repr(state))
        x = state.atoms.positions
        self.assertTrue(np.all(x >= .5) and np.all(x <= np.array([20, 30])-.5))  # half of min_distance from each wall
        self.assertEqual(Random(N=50, rho=.05, dimension=2).box.lengths.tolist(), [np.sqrt(1000)]*2)
        self.assertEqual(Random(N=50, box=20., dimension=2).box.lengths.tolist(), [20., 20.])

    def test_degrees_of_freedom_count_momentum_only_where_it_is_conserved(self):
        self.assertEqual(gas().dof, 200)                                   # walls all round: 2N
        self.assertEqual(gas(boundary='periodic').dof, 198)                # flat and periodic: 2N - 2
        mixed = gas(boundary=('walls', 'periodic'))
        self.assertEqual((mixed.walls, mixed.dof, mixed.boundary), (('left', 'right'), 199, ('walls', 'periodic')))
        self.assertEqual(FCC(N=32, rho=.5).dof, 93)                        # unchanged: 3N - 3
        slab = FCC(N=32, rho=.5, boundary=('periodic', 'periodic', 'walls'))
        self.assertEqual((slab.walls, slab.dof), (('back', 'front'), 94))

    def test_what_is_refused(self):
        with self.assertRaisesRegex(ValueError, 'three-dimensional lattice'):
            FCC(N=32, box=(10, 10))
        with self.assertRaisesRegex(ValueError, 'unbonded atoms'):
            Random(N=10, box=(20, 20), molecule='diatomic')
        with self.assertRaisesRegex(ValueError, 'unbonded atoms'):
            Random(N=10, box=20., molecule='diatomic', boundary='walls')
        with self.assertRaisesRegex(ValueError, 'between the walls'):
            System([[1., 1.], [25., 1.]], (20, 20), boundary='walls')
        with self.assertRaisesRegex(ValueError, 'boundary must be'):
            System([[1., 1.], [5., 1.]], (20, 20), boundary='wall')
        with self.assertRaisesRegex(ValueError, 'boundary must be'):
            System([[1., 1.], [5., 1.]], (20, 20), boundary=('walls', 'periodic', 'periodic'))
        with self.assertRaisesRegex(ValueError, 'Pressure control'):
            sim(ensemble='npt', pressure=.1).start(gas())
        with self.assertRaisesRegex(ValueError, 'Pressure control'):
            sim().start(gas(boundary='periodic')).set(ensemble='nph')
        with self.assertRaisesRegex(ValueError, 'far end'):
            Piston('left', to=10, speed=.1)
        with self.assertRaisesRegex(ValueError, 'no wall "top"'):
            sim().start(gas(boundary=('walls', 'periodic')), piston=Piston('top', to=10, speed=.1))
        with self.assertRaisesRegex(ValueError, 'no wall "right"'):
            sim().start(FCC(N=32, rho=.1), piston=Piston('right', to=5, speed=.1))
        with self.assertRaisesRegex(TypeError, 'must be a Piston'):
            sim().start(gas(), piston=50)
        with self.assertRaisesRegex(ValueError, 'reduced units'):
            MDSimulation(units='real', temperature=300, timestep=1, cutoff=8).start(
                System([[1., 1., 1.], [5., 5., 5.]], 30., boundary='walls', units='real'))

    def test_live_run_and_saved_run_have_two_columns(self):
        live = sim().start(gas())
        self.assertEqual((live.positions.shape, live.velocities.shape, live.box.tolist()), ((100, 2), (100, 2), [40., 40.]))
        self.assertEqual(live.settings['piston'], None)
        self.assertEqual(live.state.atoms.velocities.shape, (100, 2))
        with tempfile.TemporaryDirectory() as folder:
            run = sim(output_dir=folder).run(gas(), steps=40, thermo_every=20, trajectory_every=20, storage_name='flat')
            frames = run.load_frames()
            self.assertEqual((frames['positions'].shape, frames['velocities'].shape, frames['box'].shape), ((3, 100, 2), (3, 100, 2), (3, 2)))
            configuration = run.metadata['configuration']
            self.assertEqual((configuration['dimension'], configuration['boundary'], configuration['walls']), (2, 'walls', list(WALLS_2D)))
            self.assertEqual((run.state.dimension, run.state.walls, run.state.step), (2, WALLS_2D, 40))
            self.assertIn('"boundary": "walls"', run.view().data)

    def test_with_box_keeps_the_atoms_where_they_are_between_walls(self):
        state = settled()
        wider = state.with_box((80, 40))
        self.assertEqual(wider.box.lengths.tolist(), [80., 40.])
        np.testing.assert_array_equal(wider.atoms.positions, state.atoms.positions)
        np.testing.assert_array_equal(wider.atoms.velocities, state.atoms.velocities)
        with self.assertRaisesRegex(ValueError, 'every atom inside'):
            state.with_box((20, 40))
        periodic = gas(boundary='periodic').with_box((80, 40))            # periodic: centred, as in three dimensions
        self.assertAlmostEqual(float(periodic.atoms.positions[:, 0].mean()), float(gas(boundary='periodic').atoms.positions[:, 0].mean())+20)


class NeighbourCells(unittest.TestCase):
    def pairs(self, x, box, radius, reach, dimension, periodic):
        found = np.frombuffer(_core.neighbor_pairs(x, tuple(box), radius, reach, dimension, tuple(periodic)), dtype=np.intc)
        return found.reshape(-1, 2)

    def reference(self, x, box, radius, periodic):
        dr = x[:, None, :]-x[None, :, :]
        for d in range(3):
            if periodic[d]:
                dr[..., d] -= box[d]*np.rint(dr[..., d]/box[d])
        i, j = np.nonzero(np.triu((dr**2).sum(axis=2) < radius**2, k=1))
        return np.stack((i, j), axis=1)

    def test_cells_find_the_pairs_of_a_direct_search(self):
        rng = np.random.default_rng(3)
        cases = [(2, (30., 40., 0.), (0, 0, 0)), (2, (30., 40., 0.), (1, 1, 0)), (2, (30., 40., 0.), (0, 1, 0)),
                 (3, (12., 14., 16.), (0, 0, 0)), (3, (12., 14., 16.), (0, 1, 1)), (3, (12., 14., 16.), (1, 1, 0))]
        for dimension, box, periodic in cases:
            with self.subTest(dimension=dimension, periodic=periodic):
                x = rng.uniform(0, 1, size=(600, 3))*np.array(box)
                # Atoms caught in a bounce lie a little beyond a wall.
                for d in range(dimension):
                    if not periodic[d]:
                        x[:20, d] = -rng.uniform(0, .3, 20)
                        x[20:40, d] = box[d]+rng.uniform(0, .3, 20)
                x = np.ascontiguousarray(x)
                expected = self.reference(x, box, 2.9, periodic)
                self.assertGreater(len(expected), 300)
                for reach in (0, 1, 2):
                    np.testing.assert_array_equal(self.pairs(x, box, 2.9, reach, dimension, periodic), expected)


class Dynamics(unittest.TestCase):
    def test_energy_is_conserved_between_walls(self):
        live = sim().start(settled(), ensemble='nve')
        before = live.observables
        momentum = live.velocities.sum(axis=0)
        live.advance(20000)                                                # 100 time units, about 400 bounces
        after = live.observables
        self.assertLess(abs(after['total_energy']-before['total_energy']), 2e-4*before['kinetic_energy'])
        self.assertEqual(after['work'], 0.)
        # Each wall has taken up exactly the momentum the gas lost to it.
        change = live.velocities.sum(axis=0)-momentum
        impulse = {name: after[f'wall_impulse_{name}']-before[f'wall_impulse_{name}'] for name in WALLS_2D}
        self.assertTrue(all(value > 100 for value in impulse.values()))
        np.testing.assert_allclose(change, [impulse['left']-impulse['right'], impulse['bottom']-impulse['top']], rtol=0, atol=1e-9)

    def test_energy_is_conserved_in_a_flat_periodic_box(self):
        state = settled(gas(boundary='periodic'))
        np.testing.assert_allclose(state.atoms.velocities.sum(axis=0), 0, atol=1e-12)
        live = sim().start(state, ensemble='nve')
        before = live.observables
        live.advance(10000)
        after = live.observables
        self.assertLess(abs(after['total_energy']-before['total_energy']), 2e-4*before['kinetic_energy'])
        self.assertNotIn('work', after)
        self.assertTrue(np.all(live.positions >= 0) and np.all(live.positions < 40))
        # P = (2K + virial)/(2A): close to N T / A in a dilute gas.
        self.assertAlmostEqual(after['compressibility_factor'], 1., delta=.15)

    def test_temperature_counts_two_terms_per_atom(self):
        live = sim().start(gas())
        self.assertAlmostEqual(live.observables['temperature'], 2., places=12)   # starting velocities are scaled to it
        temperatures = []
        for _ in range(100):
            live.advance(100)
            v = live.velocities
            K = .5*np.sum(v**2)
            self.assertAlmostEqual(live.observables['temperature'], K/len(v), places=12)
            self.assertAlmostEqual(live.observables['kinetic_energy'], K, places=9)
            temperatures.append(K/len(v))
        self.assertAlmostEqual(np.mean(temperatures[20:]), 2., delta=.12)
        self.assertGreater(abs(live.velocities.sum(axis=0)).max(), .5)            # the walls take up momentum

    def test_thermostat_holds_a_changed_temperature(self):
        live = sim().start(settled())
        live.set(temperature=4.)
        live.advance(2000)
        temperatures = []
        for _ in range(300):
            live.advance(100)
            temperatures.append(live.observables['temperature'])
        self.assertAlmostEqual(np.mean(temperatures), 4., delta=.15)       # canonical scatter of one reading: 0.4
        self.assertAlmostEqual(np.std(temperatures), 4/np.sqrt(100), delta=.1)

    def test_wall_force_is_the_ideal_gas_pressure_times_the_wall_length(self):
        # 100 atoms in 60 x 40 at T = 2, light thermostat. A wall across x feels N T / Lx.
        live = sim(friction=.05).start(settled(gas(box=(60, 40)), friction=.05))
        first = live.observables
        temperatures = []
        for _ in range(200):
            live.advance(1000)
            temperatures.append(live.observables['temperature'])
        last = live.observables
        T, elapsed = np.mean(temperatures), last['time']-first['time']
        force = {name: (last[f'wall_impulse_{name}']-first[f'wall_impulse_{name}'])/elapsed for name in WALLS_2D}
        self.assertAlmostEqual(.5*(force['left']+force['right'])*60/(100*T), 1., delta=.08)
        self.assertAlmostEqual(.5*(force['bottom']+force['top'])*40/(100*T), 1., delta=.08)
        self.assertGreaterEqual(min(last[f'wall_force_{name}'] for name in WALLS_2D), 0.)

    def test_motion_reverses_between_walls(self):
        start = settled()
        live = sim().start(start, ensemble='nve')
        live.advance(400)
        self.assertGreater(abs(live.velocities.sum(axis=0)).max(), 1e-3)
        back = sim().start(System.from_arrays(live.positions, -live.velocities, live.box, boundary='walls'), ensemble='nve')
        back.advance(400)
        np.testing.assert_allclose(back.positions, start.atoms.positions, atol=1e-8)

    def test_walls_across_one_axis_in_three_dimensions(self):
        state = Random(N=108, box=12., boundary=('walls', 'periodic', 'periodic'), min_distance=1.0, seed=2)
        live = sim().start(state)
        live.advance(1000)
        self.assertAlmostEqual(live.observables['temperature'], np.sum(live.velocities**2)/(3*108-2), places=12)
        np.testing.assert_allclose(live.velocities.sum(axis=0)[1:], 0, atol=1e-11)
        nve = sim().start(live.state, ensemble='nve')
        before = nve.observables
        nve.advance(4000)
        after = nve.observables
        self.assertLess(abs(after['total_energy']-before['total_energy']), 5e-4*before['kinetic_energy'])
        self.assertEqual([key for key in after if key.startswith('wall_force')], ['wall_force_left', 'wall_force_right'])
        self.assertGreater(after['wall_impulse_right']-before['wall_impulse_right'], 10)
        x = nve.positions
        self.assertTrue(np.all(x[:, 0] > -.5) and np.all(x[:, 0] < 12.5) and np.all(x[:, 1:] >= 0) and np.all(x[:, 1:] < 12))


class PistonRuns(unittest.TestCase):
    piston = Piston('right', to=20, speed=.2)

    def test_work_of_the_piston_is_the_energy_the_gas_gains(self):
        live = sim().start(settled(), ensemble='nve', piston=self.piston)
        before = live.observables
        self.assertEqual((before['piston_position'], before['work'], before['volume']), (40., 0., 1600.))
        live.advance(10000)                                                # half way: x = 30
        middle = live.observables
        self.assertAlmostEqual(middle['piston_position'], 30., places=9)
        self.assertAlmostEqual(middle['volume'], 1200., places=6)
        self.assertAlmostEqual(live.box[0], 30., places=9)
        live.advance(12000)                                                # arrives after 20 000 steps and stays
        after = live.observables
        self.assertEqual((after['piston_position'], after['volume'], live.box.tolist()), (20., 800., [20., 40.]))
        gained = after['total_energy']-before['total_energy']
        self.assertGreater(after['work'], 100)
        self.assertAlmostEqual(gained, after['work'], delta=1e-3*after['work'])
        # Slow enough to follow the adiabat of a flat ideal gas, T V = constant, within a few per cent.
        self.assertAlmostEqual(after['temperature']/before['temperature'], 2., delta=.15)
        self.assertTrue(np.all(live.positions[:, 0] < 20.5))
        resting = live.observables['work']
        live.advance(2000)
        self.assertEqual(live.observables['work'], resting)                # a wall at rest does no work
        self.assertEqual(live.state.box.lengths.tolist(), [20., 40.])

    def test_expansion_takes_energy_out(self):
        live = sim().start(settled(), ensemble='nve', piston=Piston('top', to=60, speed=.2))
        before = live.observables
        live.advance(20000)
        after = live.observables
        self.assertEqual(live.box.tolist(), [40., 60.])
        self.assertLess(after['work'], -20)
        self.assertAlmostEqual(after['total_energy']-before['total_energy'], after['work'], delta=1e-3*abs(after['work']))

    def test_pieces_and_a_restart_give_one_long_run(self):
        start = settled()
        whole = sim().start(start, ensemble='nve', piston=self.piston)
        whole.advance(3000)
        pieces = sim().start(start, ensemble='nve', piston=self.piston)
        for steps in (1, 7, 992, 500):
            pieces.advance(steps)
        continued = sim().start(pieces.state, ensemble='nve', piston=self.piston)
        continued.advance(1500)
        np.testing.assert_array_equal(continued.positions, whole.positions)
        np.testing.assert_array_equal(continued.velocities, whole.velocities)
        np.testing.assert_array_equal(continued.box, whole.box)
        a, b = whole.observables, continued.observables
        for key in a:
            if key != 'time':                                              # summed per call
                self.assertEqual(a[key], b[key], key)
        # The piston can be set going in a run that is under way.
        later = sim().start(start, ensemble='nve')
        later.set(piston=self.piston)
        later.advance(3000)
        np.testing.assert_array_equal(later.positions, whole.positions)

    def test_saved_run_records_the_wall_and_continues(self):
        with tempfile.TemporaryDirectory() as folder:
            simulation = sim(output_dir=folder)
            start = settled()
            run = simulation.run(start, steps=2000, ensemble='nve', piston=self.piston,
                                 thermo_every=500, trajectory_every=1000, storage_name='stroke')
            thermo = run.thermo
            for column in ('wall_force_right', 'wall_impulse_right', 'wall_impulse_left', 'work', 'piston_position', 'volume'):
                self.assertIn(column, thermo.columns)
            np.testing.assert_allclose(thermo.piston_position, 40-.2*.005*(thermo.step-start.step), atol=1e-9)
            np.testing.assert_allclose(thermo.volume, 40*thermo.piston_position)
            self.assertTrue(np.all(np.diff(thermo.work) > 0) and np.all(np.diff(thermo.wall_impulse_right) > 0))
            self.assertEqual(run.metadata['configuration']['piston'], dict(wall='right', to=20., speed=.2))
            np.testing.assert_allclose(run.load_frames()['box'][:, 0], [40., 39., 38.], atol=1e-9)
            state = Run.load(run.path).state
            self.assertAlmostEqual(state.box.lengths[0], 38.)
            second = simulation.run(state, steps=1000, ensemble='nve', piston=self.piston,
                                    thermo_every=500, trajectory_every=None, storage_name='stroke-2')
            whole = simulation.start(start, ensemble='nve', piston=self.piston)
            whole.advance(3000)
            np.testing.assert_array_equal(second.state.atoms.positions, whole.positions)
            self.assertEqual(float(second.thermo.work.iloc[-1]), whole.observables['work'])
            self.assertEqual(float(second.thermo.wall_impulse_left.iloc[-1]), whole.observables['wall_impulse_left'])
            # Without a piston the next run leaves the wall where it is.
            rest = simulation.run(second.state, steps=200, ensemble='nve', thermo_every=100, trajectory_every=None, storage_name='rest')
            self.assertNotIn('piston_position', rest.thermo.columns)
            self.assertEqual(rest.state.box.lengths.tolist(), second.state.box.lengths.tolist())


if __name__ == '__main__':
    unittest.main()
