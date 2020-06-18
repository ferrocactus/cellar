from anndata import AnnData
import numpy as np

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