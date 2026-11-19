// Install the browser wheel in Pyodide under Node, the way a web page does, and run a
// short Lennard-Jones simulation with pandas absent.
//   npm install --no-save pyodide@314.0.8
//   node tests/pyodide_smoke.mjs dist/fys2160_simple_md-*-pyemscripten_2026_0_wasm32.whl
import { readFileSync } from 'node:fs';
import { basename } from 'node:path';
import { loadPyodide } from 'pyodide';

const wheel = process.argv[2];
if (!wheel) {
  console.error('usage: node tests/pyodide_smoke.mjs WHEEL');
  process.exit(2);
}
const pyodide = await loadPyodide();
await pyodide.loadPackage('micropip', { messageCallback() {} });
pyodide.FS.writeFile(`/tmp/${basename(wheel)}`, readFileSync(wheel));
pyodide.globals.set('wheel', `emfs:/tmp/${basename(wheel)}`);
console.log(await pyodide.runPythonAsync(`
import importlib.util, json, sys
import micropip
await micropip.install(wheel)  # checks the platform tag and installs numpy only
assert importlib.util.find_spec('pandas') is None, 'pandas was installed as a dependency'

import numpy as np
import fys2160_md
from fys2160_md import FCC, MDSimulation
loaded = {'pandas', 'matplotlib', 'IPython', 'anywidget'} & set(sys.modules)
assert not loaded, f'import fys2160_md loaded {loaded}'

system = FCC(N=108, rho=0.8)
sim = MDSimulation(pair_potential='lj', epsilon=1, sigma=1, ensemble='nve',
                   temperature=1, timestep=0.005, cutoff=2.5)
run = sim.run(system, steps=200, thermo_every=50, trajectory_every=100, storage_name='smoke')
table = np.genfromtxt(run.path/'thermo.csv', delimiter=',', names=True)
energy = table['total_energy']
drift = float(abs(energy[-1]-energy[0])/abs(energy[0]))
assert table['step'].tolist() == [0, 50, 100, 150, 200]
assert drift < 1e-3, drift
assert run.state.step == 200 and run.load_frames()['positions'].shape == (3, 108, 3)
try:
    run.thermo
except ImportError as exc:
    assert 'pandas is not installed' in str(exc)
else:
    raise AssertionError('Run.thermo should need pandas')
json.dumps(dict(pyodide=sys.modules['pyodide'].__version__, python=sys.version.split()[0],
                numpy=np.__version__, fys2160_md=fys2160_md.__version__,
                temperature=float(table['temperature'][-1]), relative_energy_drift=drift))
`));
