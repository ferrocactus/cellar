from setuptools import setup, find_packages

install_requires = [
    'matplotlib',
    'numpy',
    'pandas',
    'scikit-learn',
    'seaborn',
    'sklearn',
    'umap-learn',
    'kneed',
    'gtfparse',
    'statsmodels',
    'joblib',
    'gseapy',
    'anndata',
    'leidenalg',
    'scanpy',
    'Cluster_Ensembles',
    'psutil',
    'plotly'
]

setup(
    name='cellar',
    package_dir={'': 'src'},
    version=1.0,
    description='Cell Identification Pipeline.',
    author='Euxhen Hasanaj',
    author_email='ehasanaj@cs.cmu.edu',
    url='https://github.com/ferrocactus/cell_identification',
    install_requires=install_requires,
    python_requires='>3.7',
    license='MIT'
)
