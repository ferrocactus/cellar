import anndata
from anndata import AnnData
import numpy as np
from bidict import bidict

def has_key(adata, attr, key):
    if key in getattr(adata, attr):
        return True
    return False

def is_active(adata):
    return isinstance(adata, AnnData)

def get_n_obs(adata):
    return int(adata.n_obs)

def get_labels(adata):
    return adata.obs.labels.to_numpy()

def get_subsets(adata):
    return list(adata.uns['subsets'].keys())

def get_label_names(adata):
    names = np.zeros_like(adata.obs.labels).astype('U200')
    for i in adata.uns['cluster_names'].values():
        names[np.where(adata.obs.labels == adata.uns['cluster_names'].inverse[i])] = i
    return names

def get_emb_2d(adata):
    return adata.obsm.x_emb_2d.to_numpy()

def get_cluster_label_list(adata):
    return list(adata.uns['cluster_names'].keys())

def get_cluster_name_list(adata):
    return list(adata.uns['cluster_names'].values())

def get_unique_labels(adata):
    return adata.uns['cluster_info']['unique_labels']

def get_var_names(adata):
    return adata.var_names.to_numpy().astype('U')

def get_all_gene_ids(adata):
    return adata.var['parsed_ids'].to_numpy().astype('U')

def get_all_gene_names(adata):
    return adata.var['parsed_names'].to_numpy().astype('U')

def get_gene_names(adata, indices, from_r=False):
    if from_r:
        indices -= 1
    return adata.var['parsed_names'].to_numpy()[indices].astype('U')

def get_gene_names_de(adata):
    return adata.var['parsed_names'].to_numpy()[adata.uns['de']['indices']].astype('U')

def get_gene_pvals_de(adata):
    return adata.uns['de']['pvals']

def get_gene_logFC_de(adata):
    return adata.uns['de']['diffs']

def write_h5ad(adata, path, compression=9):
    adata.write_h5ad(path, compression=compression)

def read_h5ad(path):
    adata = anndata.read_h5ad(path)
    if 'cluster_names' in adata.uns:
        adata.uns['cluster_names'] = bidict(adata.uns['cluster_names'])
    return adata