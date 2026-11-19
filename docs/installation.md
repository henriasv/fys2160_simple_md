# Installation

## Install with pip

```sh
python -m pip install "fys2160-simple-md[notebook]"
```

You need **Python 3.10 or newer**. pip downloads a ready-built package (a wheel)
for Windows (x64 and arm64), macOS (Intel and Apple Silicon) and Linux (x86_64 and
aarch64), so no C compiler is needed. There is no Atomify or LAMMPS dependency.

!!! warning "Before the first release on PyPI"
    The command above works once a release has been published on
    [PyPI](https://pypi.org/project/fys2160-simple-md/). If pip answers that it
    cannot find `fys2160-simple-md`, no release is out yet: use
    [Install from GitHub](#install-from-github) below, which needs Git and a C compiler.
    The same holds for the browser. Windows is first built and tested when that
    release is made; building from source on Windows has not been validated.

Open a notebook in Jupyter or VS Code and select the Python environment where
you installed the package. The distribution is called `fys2160-simple-md`;
its Python import is **`fys2160_md`**.

```python
import matplotlib.pyplot as plt
from fys2160_md import FCC, MDSimulation

system = FCC(N=108, rho=.1)
sim = MDSimulation(
    pair_potential="lj",
    epsilon=1,
    sigma=1,
    ensemble="nvt",
    temperature=2,
    timestep=0.005,
    cutoff=2.5,
)
result = sim.run(
    system,
    steps=1000,
    thermo_every=100,
    trajectory_every=500,
    storage_name="first-run",
)
data = result.thermo
plt.figure(figsize=(8, 5))
plt.subplot(2, 1, 1)
plt.plot(data.time, data.temperature)
plt.ylabel("Temperature")
plt.subplot(2, 1, 2)
plt.plot(data.time, data.pressure)
plt.ylabel("Pressure")
plt.xlabel("Time")
plt.tight_layout()
plt.show()
result.view()
```

`timestep` sets the integration step; `thermo_every` and `trajectory_every` set output
intervals in MD steps. All main settings are explicit so you can edit them.

`show=True` needs the notebook extras and a working widget renderer. The 3D
viewer needs WebGL 2 / graphics acceleration. Restart the kernel after updating
the package so its Python code and embedded viewer are reloaded.

## Reproducible installation

Name the release when sharing a fixed course version:

```sh
python -m pip install "fys2160-simple-md[notebook]==VERSION"
```

Replace `VERSION` with a release number from
[PyPI](https://pypi.org/project/fys2160-simple-md/#history).

## In the browser

Each release also carries a WebAssembly wheel for [Pyodide](https://pyodide.org) 314
(Python 3.14), so the same package and the same Python code run in a web page:

```python
import micropip
await micropip.install("fys2160-simple-md==VERSION")

from fys2160_md import FCC, MDSimulation
```

pandas is left out in the browser because it is large and only `Run.thermo` uses it.
Either install it as well with `await micropip.install("pandas")`, or read the same
table with numpy:

```python
import numpy as np
data = np.genfromtxt(run.path / "thermo.csv", delimiter=",", names=True)
data["temperature"]
```

Saved runs are written to Pyodide's in-memory file system and disappear when the
page is reloaded. A wheel is built for one Pyodide generation; a page that moves to
a newer Pyodide needs a release built for it.

## Install from GitHub

The development version, or any commit, can be installed from source:

```sh
python -m pip install "fys2160-simple-md[notebook] @ git+https://github.com/henriasv/fys2160_simple_md.git"
```

This builds the small native solver during installation, so it needs Git and a C compiler:

- **macOS:** install Xcode command-line tools with `xcode-select --install`.
- **Ubuntu/Debian:** install `build-essential`, `python3-dev` and `git`.
- **Windows:** install the Visual Studio Build Tools with the C++ workload
  (Visual Studio 2019 version 16.8 or newer).

Add `@COMMIT` after `.git` (an actual Git commit hash) to fix the version.

## Develop locally

```sh
git clone https://github.com/henriasv/fys2160_simple_md.git
cd fys2160_simple_md
python -m pip install -e '.[notebook,docs]'
python -m unittest discover -s tests
python scripts/build_docs.py --execute
python -m mkdocs serve
```

The [executed example](examples/notebook.md) can be explored without installing anything.
