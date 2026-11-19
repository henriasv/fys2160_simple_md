"""Importing the package must stay cheap: pandas and the notebook stack load on use."""
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from fys2160_md import FCC, MDSimulation


class ImportTests(unittest.TestCase):
    @unittest.skipIf(sys.platform == 'emscripten', 'needs a fresh interpreter')
    def test_import_loads_no_optional_packages(self):
        code = ('import sys, fys2160_md; '
                'print([m for m in ("pandas", "matplotlib", "IPython", "anywidget") if m in sys.modules])')
        result = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, check=True)
        self.assertEqual(result.stdout.strip(), '[]')

    def test_run_works_without_pandas_and_thermo_says_what_is_missing(self):
        with tempfile.TemporaryDirectory() as folder:
            sim = MDSimulation(pair_potential='lj', epsilon=1, sigma=1, ensemble='nvt', temperature=2,
                               timestep=.005, cutoff=2.5, output_dir=folder)
            with patch.dict(sys.modules, {'pandas': None}):
                run = sim.run(FCC(N=108, rho=.1), steps=20, thermo_every=10, trajectory_every=None,
                              storage_name='no-pandas')
                with self.assertRaisesRegex(ImportError, 'pandas is not installed'):
                    run.thermo
                table = np.genfromtxt(run.path/'thermo.csv', delimiter=',', names=True)
            self.assertEqual(table['step'].tolist(), [0, 10, 20])
            self.assertEqual(run.state.step, 20)


if __name__ == '__main__':
    unittest.main()
