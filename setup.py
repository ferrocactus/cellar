from setuptools import setup, find_packages

setup(
    name='scRNA-Identification-Project',
    version=1.0,
    description='Workflow for identifying cell types.',
    author='Euxhen Hasanaj',
    author_email='ehasanaj@cs.cmu.edu',
    url='https://github.com/ferrocactus/cell_identification',
    packages=find_packages(),
    python_requires='3.7',
    license='MIT'
)