// Runs the week 34 scenes (examples/week34) in Pyodide under Node, with a browser
// wheel, and prints the numbers and the time they take. Compare with the native
// build: python examples/week34/page_numbers.py
//   npm install --no-save pyodide@314.0.8
//   node scripts/pyodide_week34.mjs dist/fys2160_simple_md-*-pyemscripten_2026_0_wasm32.whl
import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { loadPyodide } from 'pyodide';

const wheel = process.argv[2];
if (!wheel) {
  console.error('usage: node scripts/pyodide_week34.mjs WHEEL');
  process.exit(2);
}
const scenes = join(dirname(fileURLToPath(import.meta.url)), '..', 'examples', 'week34');
const pyodide = await loadPyodide();
await pyodide.loadPackage(['numpy'], { messageCallback() {} });
pyodide.unpackArchive(new Uint8Array(readFileSync(wheel)), 'wheel');
for (const name of ['box.py', 'piston.py', 'page_numbers.py']) {
  pyodide.FS.writeFile(`/home/pyodide/${name}`, new Uint8Array(readFileSync(join(scenes, name))));
}
console.log(pyodide.runPython(`
import contextlib, io, runpy
printed = io.StringIO()
with contextlib.redirect_stdout(printed):
    runpy.run_path('/home/pyodide/page_numbers.py', run_name='__main__')
printed.getvalue()
`));
