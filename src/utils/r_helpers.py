import anndata
from anndata import AnnData
import numpy as np
import pandas as pd
from bidict import bidict
from .. import name_genes
from scipy.sparse import issparse
from .validation import InappropriateArgument
import os
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
matplotlib.use('agg')


def has_key(adata, attr, key):
    if key in getattr(adata, attr):
        return True
    return False


def get_key(adata, attr, key):
    return adata[attr][key]


def has_labels(adata):
    if 'labels' in adata.obs:
        return True
    return False


def has_emb(adata):
    if 'x_emb' in adata.obsm:
        return True
    return False


def has_key_tri(adata, key1, key2, key3):
    if key2 not in getattr(adata, key1):
        return False
    if key3 not in getattr(adata, key1)[key2]:
        return False
    return True


def get_key_tri(adata, key1, key2, key3):
    return getattr(adata, key1)[key2][key3]


def more_than_one_subset(adata):
    return len(adata.uns['cluster_names']) > 1


def is_active(adata):
    return isinstance(adata, AnnData)


def is_str(obj):
    return isinstance(obj, str)


def get_n_obs(adata):
    return int(adata.n_obs)


def get_n_clusters(adata):
    if 'cluster_info' in adata.uns:
        if 'n_clusters' in adata.uns['cluster_info']:
            return adata.uns['cluster_info']['n_clusters']
    else:
        return -1


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
    return []


def get_cluster_name_list(adata):
    if 'cluster_names' in adata.uns:
        return list(adata.uns['cluster_names'].values())
    return []


def get_cluster_names(adata):
    if 'cluster_names' in adata.uns:
        return adata.uns['cluster_names']
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
    adata = adata.copy()
    # adata.var.pop('parsed_names')
    if 'parsed_ids' in adata.var:
        adata.var.pop('parsed_ids')
        # adata.var.pop('parsed_names')
    if compression == "None":
        compression = None
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


def is_sparse(adata):
    return issparse(adata.X)


def has_x_emb_sparse(adata, method, n_components):
    if n_components == 'knee':
        n_components = 50

    n_components = int(n_components)

    if 'dim_reduction_info' not in adata.uns:
        return False

    if adata.uns['dim_reduction_info']['method'] != method:
        return False

    if adata.uns['dim_reduction_info']['n_components'] != n_components:
        return False

    return True


def set_adata(x, transpose=False):
    a = AnnData(x)
    if transpose:
        a = a.T
    return a


def read_10x(file):
    # read 10x data
    # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.settings.verbosity = 3
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=80, facecolor='white')

    # the file that will store the analysis results
    results_file = 'write/pbmc3k.h5ad'

    # os.chdir(file)
    dir_name = ''
    path10x = file
    for i in os.listdir(file):
        if (i[:6] == 'filter'):
            dir_name = i
            if (os.listdir(file+'/'+dir_name)[0] == 'hg19'):
                path10x = file+'/'+dir_name+'/hg19/'  # the directory with the `.mtx` file
            else:
                path10x = file+'/'+dir_name+'/'  # the directory with the `.mtx` file

    adata = sc.read_10x_mtx(
        path10x,
        # use gene symbols for the variable names (variables-axis index)
        var_names='gene_symbols',
        cache=True)
    adata.var_names_make_unique()
    # os.chdir('../../..')
    adata.uns['dataset'] = 'data10x'
    return AnnData(adata)


def generate_violin(adata, labels, genename='S100A6', v1=-1, v2=10):
    matplotlib.use('agg')
    x = np.array(adata)
    x = x.reshape(-1, 1)
    try:
        labels = np.array(labels)
        labels = labels.reshape((x.shape[0], 1))

        x = np.hstack((x, labels))
        genes = [genename, 'cluster']

        df = pd.DataFrame(x, columns=genes)
        df.drop(df[df[genename] < v1].index, inplace=True)
        df.drop(df[df[genename] > v2].index, inplace=True)
        sns.set_theme(style="whitegrid")
        df = df.astype({'cluster': 'int32'})
        plt.figure()
        ax = sns.violinplot(x="cluster", y=genename, data=df)
        fig = ax.get_figure()
        fig.savefig('violin'+genename+'.png')
    except:
        return -1

    # return ax
    return v2


def get_color_by(adata, key):
    if key in adata.obs:
        return np.array(getattr(adata.obs, key))
    else:
        return "error"


def get_obs_keys(adata):
    temp = list(adata.obs.keys())
    if 'labels' in temp:
        temp.remove('labels')
    return np.array(temp)
