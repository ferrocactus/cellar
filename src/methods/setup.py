from setuptools import setup
from Cython.Build import cythonize
import numpy

setup(
    ext_modules = cythonize(["cellar/src/methods/_k_means_fast.pyx",
                            "cellar/src/methods/_k_means_elkan.pyx"]),
    include_dirs=[numpy.get_include()]
)