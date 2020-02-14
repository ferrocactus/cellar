import anndata
import gtfparse

def read_h5ad(file):
    ann = anndata.read_h5ad(file)
    return ann

def read_gtf(file):
    gtf = gtfparse.read_gtf(file)
    return gtf