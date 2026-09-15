import os
from setuptools import Extension, setup

setup(ext_modules=[Extension(
    "fys2160_md._core", ["src/fys2160_md/_core.c"],
    extra_compile_args=["/O2", "/std:c11"] if os.name == "nt" else ["-O3", "-std=c99"],
    libraries=[] if os.name == "nt" else ["m"],
)])
