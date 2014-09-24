from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import os
import shutil

# set C and C++ compiler
os.environ["CC"] = "g++"
os.environ["CXX"] = "g++"

src = ["ppi_networkit.pyx"]  # list of source files

modules = [Extension("ppi_networkit",
           src,
           language="c++",
           extra_compile_args=["-fopenmp", "-std=c++11", "-O3", "-DNOGTEST"],
           include_dirs=['../src/', '../ext/SQLiteCpp/include',
                         '../ext/NetworKit/networkit/cpp',
                         '../ext/NetworKit/networkit'],
           extra_link_args=["-fopenmp", "-std=c++11"],
           libraries=["ppi_networkit", "SQLiteCpp", "sqlite3",
                      "NetworKit-Core-Opt"],
           library_dirs=["../build/src/", "../build/ext/SQLiteCpp",
                         "../ext/NetworKit"])]


setup(name="ppi_networkit",
      author="Patrick Flick",
      cmdclass={"build_ext": build_ext},
      #py_modules = ["NetworKit.py"],
      ext_modules=modules)
