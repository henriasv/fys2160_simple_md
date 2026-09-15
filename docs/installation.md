# Installation

## Install directly from GitHub

```sh
python -m pip install "fys2160-simple-md[notebook] @ git+https://github.com/henriasv/fys2160_simple_md.git"
```

You need **Python 3.10 or newer**, Git and a C compiler. The package builds its
small native solver during installation. There is no Atomify, LAMMPS or Pyodide dependency.

- **macOS:** install Xcode command-line tools with `xcode-select --install`.
- **Ubuntu/Debian:** install `build-essential`, `python3-dev` and `git`.
- **Windows:** a compatible Visual Studio C/C++ build environment is required;
  Windows has not yet been validated.

Open a notebook in Jupyter or VS Code and select the Python environment where
you installed the package. The distribution is called `fys2160-simple-md`;
its Python import is **`fys2160_md`**.

```python
import matplotlib.pyplot as plt
from fys2160_md import FCC, MDSimulation

system = FCC(N=108, rho=.1)
sim = MDSimulation(
    system,
    pair_potential="lj",
    epsilon=1,
    sigma=1,
    ensemble="nvt",
    temperature=2,
    timestep=0.005,
    cutoff=2.5,
)
result = sim.run(
    steps=1000,
    sample_every=100,
    save_every=500,
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

`timestep` sets the integration step; `sample_every` and `save_every` set output
intervals in MD steps. All main settings are explicit so you can edit them.

`show=True` needs the notebook extras and a working widget renderer. The 3D
viewer needs WebGL 2 / graphics acceleration. Restart the kernel after updating
the package so its Python code and embedded viewer are reloaded.

## Reproducible installation

Use a commit instead of the moving main branch when sharing a fixed course version:

```sh
python -m pip install "fys2160-simple-md[notebook] @ git+https://github.com/henriasv/fys2160_simple_md.git@COMMIT"
```

Replace `COMMIT` with an actual Git commit hash.

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
