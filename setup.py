from setuptools import setup, find_packages
from pip.req import parse_requirements

#install_reqs = parse_requirements('requirements.txt')
#reqs = [str(ir.req) for ir in install_reqs]

install_requires = ['anndata',
                    'h5py',
                    'jupyter',
                    'matplotlib',
                    'numpy',
                    'pandas',
                    'scikit-learn',
                    'seaborn',
                    'sklearn',
                    'tqdm',
                    'umap-learn',
                    'git+https://github.com/idc9/jackstraw']

setup(
    name='acip',
    version=1.0,
    description='Workflow for identifying cell types.',
    author='Euxhen Hasanaj',
    author_email='ehasanaj@cs.cmu.edu',
    url='https://github.com/ferrocactus/cell_identification',
    packages=find_packages(),
    install_requires=install_requires,
    python_requires='>3.7',
    license='MIT'
)