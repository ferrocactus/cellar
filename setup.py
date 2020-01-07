from setuptools import setup, find_packages

install_requires = ['numpy', 'pandas', 'matplotlib', 'seaborn', 'umap',
                'scikit-learn', 'anndata', 'jupyter', 'jackstraw', 'tqdm']

setup(
    name='src',
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