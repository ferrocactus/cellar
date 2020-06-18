import itertools
import json
import re
from ast import literal_eval
from functools import reduce

import anndata
import numpy as np
import pandas as pd


def parse(x):
    if x.size == 0:
        return x
    parsed_x = np.char.split(x.flatten(), sep='.', maxsplit=1)
    parsed_x = np.array([i[0] for i in parsed_x])
    parsed_x = np.char.replace(parsed_x, '-', '')
    parsed_x = np.char.replace(parsed_x, ' ', '')
    # Remove everything that goes inside brackets
    parsed_x = [re.sub(r" ?\([^)]+\)", "", item) for item in parsed_x]
    parsed_x = np.char.upper(parsed_x)
    parsed_x = parsed_x.reshape(x.shape)
    return parsed_x


def _emb_exists_in_adata(adata, method, n_components):
    if 'x_emb' not in adata.obsm:
        return None
    if adata.uns['dim_reduction_info']['method'] != method:
        return None
    # The following will pass also if n_components_used = n_components = knee
    if n_components == adata.uns['dim_reduction_info']['n_components_used']:
        return adata.obsm['x_emb']
    if isinstance(n_components, str):
        return None
    if n_components <= adata.uns['dim_reduction_info']['n_components']:
        return adata.obsm['x_emb'][:, :n_components]
    return None


def _2d_emb_exists_in_adata(adata, method, dim, use_emb):
    if f'x_emb_{dim}d' not in adata.obsm:
        return False
    if adata.uns[f'visualization_info_{dim}d']['method'] != method:
        return False
    if adata.uns[f'visualization_info_{dim}d']['used_emb'] != use_emb:
        return False
    return True


def store_subset(adata, name, indices, from_r=False):
    # ARRAYS START FROM 0
    if from_r:
        indices = np.array(indices).astype(int) - 1
    if 'subsets' not in adata.uns:
        adata.uns['subsets'] = {}
    adata.uns['subsets']['name'] = indices


def update_subset_label(adata, name):
    name = str(name)
    if name not in adata.uns['subsets']:
        raise ValueError("Subset not found in adata.")

    indices = adata.uns['subsets'][name]

    # adata.uns['cluster_names] is a bidict
    # with keys=cluster ids and values=cluster names
    if name in adata.uns['cluster_names'].values():
        # name exists, so simply update cell labels
        label = adata.uns['cluster_names'].inverse[name]
        adata.obs['labels'][indices] = label
    else:
        # name does not exist, so add it
        # careful not to mess up existing label configuration
        labels_cp = adata.obs['labels'].copy()
        temp_label = np.max(labels_cp) + 10
        adata.uns['cluster_names'][temp_label] = name
        labels_cp[indices] = temp_label

        unq_new_labels = np.unique(labels_cp)
        for label in adata.uns['cluster_names']:
            # remove label if it dissapeared after update
            if label not in unq_new_labels:
                adata.uns['cluster_names'].pop(label, None)

        # create empty arrays
        new_cluster_names = bidict({})
        new_labels = np.zeros_like(labels_cp) * (-1)
        # assign labels starting from 0
        for i, label in enumerate(adata.uns['cluster_names']):
            new_labels[labels_cp == label] = i
            new_cluster_names['i'] = adata.uns['cluster_names'][label]

        adata.obs['labels'] = new_labels
        adata.uns['cluster_names'] = new_cluster_names

        print(adata.uns['cluster_names'])