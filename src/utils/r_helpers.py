import anndata
from anndata import AnnData
import numpy as np
import pandas as pd
from bidict import bidict
from .. import name_genes
from scipy.sparse import issparse
from .validation import InappropriateArgument


def has_key(adata, attr, key):
    if key in getattr(adata, attr):
        return True
    return False


def get_key(adata, attr, key):
    return adata[attr][key]


def has_key_tri(adata, key1, key2, key3):
    if key2 not in getattr(adata, key1):
        return False
    if key3 not in getattr(adata, key1)[key2]:
        return False
    return True


def get_key_tri(adata, key1, key2, key3):
    return getattr(adata, key1)[key2][key3]


def is_active(adata):
    return isinstance(adata, AnnData)


def is_str(obj):
    return isinstance(obj, str)


def get_n_obs(adata):
    return int(adata.n_obs)


def get_labels(adata):
    if 'labels' not in adata.obs:
        return "No labels found"
    return adata.obs['labels'].to_numpy()


def get_subsets(adata):
    return list(adata.uns['subsets'].keys())


def get_uncertain_subsets(adata):
    return list(adata.uns['uncertain_subsets'].keys())


def get_label_names(adata):
    names = np.zeros_like(adata.obs.labels).astype('U200')
    if isinstance(adata.uns['cluster_names'], bidict) == False:
        adata.uns['cluster_names'] = bidict(adata.uns['cluster_names'])
        for key in list(adata.uns['cluster_names'].keys()):
            adata.uns['cluster_names'][int(key)] = \
                adata.uns['cluster_names'].pop(key, None)
        #name_genes(adata, inplace=True)

    for i in adata.uns['cluster_names'].values():
        names[np.where(adata.obs.labels ==
                       adata.uns['cluster_names'].inverse[i])] = i
    return names


def get_emb_2d(adata):
    return adata.obsm.x_emb_2d.to_numpy()


def get_cluster_label_list(adata):
    if 'cluster_names' in adata.uns:
        return list(adata.uns['cluster_names'].keys())
    else:
        return []


def get_cluster_name_list(adata):
    if 'cluster_names' in adata.uns:
        return list(adata.uns['cluster_names'].values())
    else:
        return []


def get_unique_labels(adata):
    return adata.uns['cluster_info']['unique_labels']


def get_var_names(adata):
    return adata.var_names.to_numpy().astype('U')


def get_obs_names(adata):
    return adata.obs_names.to_numpy().astype('U')


def get_all_gene_ids(adata):
    if issparse(adata.X):
        return []
    return adata.var['parsed_ids'].to_numpy().astype('U')


def get_all_gene_names(adata):
    if issparse(adata.X):
        return []
    return adata.var['parsed_names'].to_numpy().astype('U')


def get_gene_names(adata, indices, from_r=False):
    if from_r:
        indices -= 1
    return adata.var['parsed_names'].to_numpy()[indices].astype('U')


def get_gene_names_de(adata):
    return adata.var['parsed_names'].to_numpy()[adata.uns['de'].index].astype('U')


def get_gene_pvals_de(adata):
    if 'de' not in adata.uns:
        return None
    return adata.uns['de']['qvals']


def get_gene_logFC_de(adata):
    if 'de' not in adata.uns:
        return None
    return adata.uns['de']['log2fc']


def get_de_table(adata):
    if 'de' not in adata.uns:
        return None
    return pd.DataFrame(adata.uns['de'])


def write_h5ad(adata, path, compression=9):
    adata.var['parsed_names'] = None
    adata.var['parsed_ids'] = None
    adata.write_h5ad(path, compression=compression)


def read_h5ad(path):
    try:
        adata = anndata.read_h5ad(path)
    except:
        print('Error reading file')
        return "file_error"

    if 'cluster_names' in adata.uns:
        adata.uns['cluster_names'] = bidict(adata.uns['cluster_names'])
        for key in list(adata.uns['cluster_names'].keys()):
            adata.uns['cluster_names'][int(key)] = \
                adata.uns['cluster_names'].pop(key, None)
    name_genes(adata, inplace=True)
    return adata


def write_key(adata, keyname, value):
    adata.uns[keyname] = np.array(value).reshape(-1)


def get_col(adata, i, from_r=True):
    if from_r:
        i -= 1
    if issparse(adata.X):
        return adata.X[:, i].todense()
    return adata.X[:, i]


def get_x(adata):
    if issparse(adata.X):
        raise InappropriateArgument("Sparse matrix detected. "
                                    "Cannot run chosen method.")
    if adata.isbacked:
        return adata.X[()]
    return adata.X
