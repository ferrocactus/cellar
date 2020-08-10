import os
from setuptools import setup
from setuptools import find_packages
from setuptools import Command
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
    'statsmodels',
    'joblib',
    'anndata',
    'leidenalg',
    'scanpy',
    'psutil',
    'plotly',
    'bidict',
    'pydiffmap',
    'diffxpy',
    'scipy'
]


class CleanCommand(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')


cmdclass = {'clean': CleanCommand}

setup(
    name='cellar',
    package_dir={'cellar': 'src'},
    packages=['cellar', 'cellar.units', 'cellar.utils', 'cellar.methods'],
    version=1.0,
    include_package_data=True,
    description='Cell Annotation Robot.',
    author='Systems Biology Group, Carnegie Mellon Universityj',
    author_email='ehasanaj@cs.cmu.edu',
    url='https://github.com/ferrocactus/cellar',
    install_requires=install_requires,
    python_requires='>3.7',
    license='MIT',
    cmdclass=cmdclass,
    include_dirs=[numpy.get_include()]
)
