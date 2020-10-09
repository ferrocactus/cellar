from src import BinToGene
import anndata

a = anndata.read_h5ad('resources/cell_by_bin.h5ad')
btg = BinToGene(n_jobs=4)
counts, ids = btg.convert(a.X, a.var_names)

