import os
from setuptools import Extension, setup

sources = ["_core.c", "molecular.c", "pppm.c", "fft.c"]
setup(ext_modules=[Extension(
    "fys2160_md._core", [f"src/fys2160_md/{name}" for name in sources],
    depends=[f"src/fys2160_md/{name}" for name in ("molecular.h", "pppm.h", "fft.h")],
    extra_compile_args=["/O2", "/std:c11"] if os.name == "nt" else ["-O3", "-std=c99"],
    libraries=[] if os.name == "nt" else ["m"],
)])
