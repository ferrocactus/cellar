import os
from setuptools import setup
from setuptools import find_packages
from setuptools import Command
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from setuptools import Extension
import numpy

install_requires = [
    'matplotlib',
    'numpy',
    'pandas',
    'scikit-learn',
    'sklearn',
    'umap-learn',
    'kneed',
    'gtfparse',
    'statsmodels',
    'joblib',
    'anndata',
    'leidenalg',
    'scanpy',
    'psutil',
    'plotly',
    'bidict'
]

ext_modules = [
    Extension("cellar.methods._k_means_fast", [
              "src/methods/_k_means_fast.pyx"]),
    Extension("cellar.methods._k_means_elkan", [
              "src/methods/_k_means_elkan.pyx"]),
]


class CleanCommand(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')


cmdclass = {'build_ext': build_ext, 'clean': CleanCommand}

setup(
    name='cellar',
    package_dir={'cellar': 'src'},
    packages=['cellar', 'cellar.units', 'cellar.utils', 'cellar.methods'],
    data_files=[('datasets', ['datasets/RNAseqTPM.csv'])],
    version=1.0,
    include_package_data=True,
    description='Cell Annotation Robot.',
    author='Euxhen Hasanaj',
    author_email='ehasanaj@cs.cmu.edu',
    url='https://github.com/ferrocactus/cellar',
    install_requires=install_requires,
    python_requires='>3.7',
    license='MIT',
    cmdclass=cmdclass,
    ext_modules=cythonize(ext_modules),
    include_dirs=[numpy.get_include()]
)
