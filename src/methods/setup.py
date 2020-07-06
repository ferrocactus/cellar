# Run if want to compile only cython
# From root folder run
#       python src/methods/setup.py build_ext --inplace
# then
#       python src/methods/setup.py clean
from Cython.Build import cythonize
from setuptools import Command
from setuptools import setup
from setuptools import Extension

import numpy
import os

class CleanCommand(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')

ext_modules = [
    Extension("src.methods._k_means_fast", [
              "src/methods/_k_means_fast.pyx"]),
    Extension("src.methods._k_means_elkan", [
              "src/methods/_k_means_elkan.pyx"]),
]

cmdclass = {'clean': CleanCommand}

setup(
    cmdclass=cmdclass,
    ext_modules = cythonize(ext_modules),
    include_dirs=[numpy.get_include()]
)

