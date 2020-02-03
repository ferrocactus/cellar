from setuptools import setup, find_packages

install_requires = [
    'matplotlib',
    'numpy',
    'pandas',
    'scikit-learn',
    'seaborn',
    'sklearn',
    'umap-learn',
    'tqdm',
    'kneed'
]

setup(
    name='acip',
    version=1.0,
    description='Automatic Cell Identification Platform.',
    author='Euxhen Hasanaj',
    author_email='ehasanaj@cs.cmu.edu',
    url='https://github.com/ferrocactus/cell_identification',
    packages=find_packages(),
    install_requires=install_requires,
    python_requires='>3.7',
    license='MIT'
)