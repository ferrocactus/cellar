import itertools
import json
import re
from ast import literal_eval
from functools import reduce

import anndata
import numpy as np
import pandas as pd
from bidict import bidict

from .exceptions import InvalidArgument
from .exceptions import InappropriateArgument
from .validation import _validate_cluster_list


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
    if 'x_emb' not in adata.obsm and 'X_pca' not in adata.obsm:
        return None
    if 'dim_reduction_info' in adata.uns:
        if adata.uns['dim_reduction_info']['method'] != method:
            return None
    if isinstance(n_components, str):
        n_components = 10
    if 'x_emb' in adata.obsm:
        if n_components <= adata.uns['dim_reduction_info']['n_components']:
            return adata.obsm['x_emb'][:, :n_components]
    if 'X_pca' in adata.obsm:
        if n_components <= adata.obsm['X_pca'].shape[1]:
            return adata.obsm['X_pca'][:, :n_components]
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
    adata.uns['subsets'][name] = indices


def store_labels(adata, labels, method):
    adata.obs['labels'] = np.array(labels).astype(int)
    # Update entries
    adata.uns['cluster_info'] = {}
    unq_labels = np.unique(adata.obs['labels'])
    adata.uns['cluster_info']['unique_labels'] = unq_labels
    adata.uns['cluster_info']['n_clusters'] = len(unq_labels)
    adata.uns['cluster_info']['method'] = method
    adata.uns['cluster_names'] = bidict(
        {i: str(i) for i in unq_labels})
    populate_subsets(adata)


def update_subset_label(adata, subset_name, name):
    subset_name = str(subset_name)
    name = str(name)

    if subset_name not in adata.uns['subsets']:
        raise InvalidArgument(f"Subset {subset_name} not found in adata.")

    indices = adata.uns['subsets'][subset_name]

    # adata.uns['cluster_names] is a bidict
    # with keys=cluster ids and values=cluster names
    if name in list(adata.uns['cluster_names'].values()):
        label = adata.uns['cluster_names'].inverse[name]
    else:  # name does not exist, so add it
        label = np.max(adata.obs['labels']) + 1
        adata.uns['cluster_names'][label] = name

    adata.obs['labels'][indices] = label
    unq_new_labels = np.unique(adata.obs['labels'])

    for old_label in list(adata.uns['cluster_names'].keys()):
        if old_label not in unq_new_labels:
            adata.uns['cluster_names'].pop(old_label, None)

    # Update entries
    adata.uns['cluster_info'] = {}
    unq_labels = np.unique(adata.obs['labels'])
    adata.uns['cluster_info']['unique_labels'] = unq_labels
    adata.uns['cluster_info']['n_clusters'] = len(unq_labels)
    adata.uns['cluster_info']['method'] = "Modified"
    populate_subsets(adata)


def populate_subsets(adata):
    # Populate subsets with the found clusters and
    # any don't remove user defined subsets
    if 'subsets' not in adata.uns:
        adata.uns['subsets'] = {}

    # Remove old clusters
    for subset in list(adata.uns['subsets'].keys()):
        if subset[:7] == 'Cluster':
            adata.uns['subsets'].pop(subset, None)
    # Update with new clusters
    unq_labels = np.unique(adata.obs['labels'])
    for label in unq_labels:
        adata.uns['subsets'][f'Cluster_{label}'] = \
            np.squeeze(np.where(adata.obs['labels'] == label))


def merge_clusters(adata, clusters):
    # Merge all clusters to the first one
    clusters = _validate_cluster_list(adata.obs['labels'], clusters)
    if len(clusters) < 2:
        raise InvalidArgument("Not enough clusters found to merge")

    c0 = clusters[0]
    c0_name = adata.uns['cluster_names'][c0]

    for cluster in clusters[1:]:
        update_subset_label(adata, f'Cluster_{cluster}', c0_name)

def get_dic():
    f=open("markers/cl-simple.json","rb")
    onto_data=json.load(f)
    cells=onto_data['graphs'][0]['nodes']    # dictionary of length 2632
    tissue=[]
    cell_type=[]
    tissues=['kidney', 'thymus', 'spleen', "liver",'lymph', 'stomach', 'heart',
         'small intestine', 'large intestine','blood',"brain","thyroid",
         'placenta','eye','other', "embryo", "muscle"]

    identifiers=[['kidney'], ['thymocyte'], ['splenic'], ["hepat","liver"],
            ['T','B'], ['stomach','gastric'], ['heart','cardiac'],
         ['small intestine'], ['large intestine'],['blood'],["brain"],["thyroid"],
         ['placenta'],['eye','retina'], ['other'],["embryo"], ["muscle",'myo']]

    for i in cells:
        f=0
        if ('lbl' in i.keys()):
            for j in range(len(identifiers)):
                for k in identifiers[j]:
                    if ((k in i['lbl']) and f==0):
                        tissue.append(tissues[j])
                        f=1
            if (f==0):
                tissue.append('other')
            cell_type.append(i['id'].split('/')[-1]+': '+i['lbl'])
        else:
            tissue.append('other')
            cell_type.append(i['id'].split('/')[-1])
    dic={}
    for i in tissue:
        dic[i]=[]
    for i in range(len(cell_type)):
        dic[tissue[i]].append(cell_type[i])
    return dic

