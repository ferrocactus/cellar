import itertools
import json
import re
import os


from typing import Optional, Union
from anndata import AnnData
from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import kneighbors_graph

from ast import literal_eval
from functools import reduce

import anndata
import numpy as np
from bidict import bidict

from .exceptions import InvalidArgument
from .exceptions import InappropriateArgument
from .validation import _validate_cluster_list


this_dir = os.path.dirname(__file__)


def join_root(path):
    return os.path.abspath(os.path.join(this_dir, path))


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
        elif adata.obsm['x_emb'].shape[1] == n_components:
            return adata.obsm['x_emb']
    if method != 'PCA':
        return None
    if isinstance(n_components, str):
        n_components = 10
    if 'x_emb' in adata.obsm:
        if n_components <= adata.obsm['x_emb'].shape[1]:
            return adata.obsm['x_emb'][:, :n_components]
    if 'X_pca' in adata.obsm:
        if n_components <= adata.obsm['X_pca'].shape[1]:
            return adata.obsm['X_pca'][:, :n_components]
    return None


def _labels_exist_in_adata(adata, method, n_clusters):
    if 'labels' not in adata.obs:
        return False
    # Manually labelled file, used for debugging
    # if adata.uns['cluster_info']['method'] == 'F':
    #    return True
    if method != adata.uns['cluster_info']['method']:
        return False
    if np.array_equal(n_clusters, adata.uns['cluster_info']['n_clusters_used']):
        return True
    return False


def _2d_emb_exists_in_adata(adata, method, dim, use_emb):
    if f'x_emb_{dim}d' not in adata.obsm:
        return False
    if adata.uns[f'visualization_info_{dim}d']['method'] != method:
        return False
    if adata.uns[f'visualization_info_{dim}d']['used_emb'] != use_emb:
        return False
    return True


def merge_cluster_names(adata, ref):
    # Currently used for alignment only
    if 'cluster_names' not in ref.uns or 'cluster_names' not in adata.uns:
        return

    for i in adata.uns['cluster_names']:
        if i in ref.uns['cluster_names']:
            adata.uns['cluster_names'][i] = ref.uns['cluster_names'][i]


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
        {int(i): str(i) for i in unq_labels})
    populate_subsets(adata)


def store_x_emb(adata, x_emb, method):
    adata.obsm['x_emb'] = np.array(x_emb).astype(float)
    adata.uns['dim_reduction_info'] = {}
    adata.uns['dim_reduction_info']['method'] = method
    adata.uns['dim_reduction_info']['n_components'] = x_emb.shape[1]
    adata.uns['dim_reduction_info']['n_components_used'] = x_emb.shape[1]
    adata.uns['dim_reduction_info']['kwargs'] = {}


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
            np.reshape(np.where(adata.obs['labels'] == label), (-1))


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
    f = open(join_root("../markers/cl-simple.json"), "rb")
    onto_data = json.load(f)
    cells = onto_data['graphs'][0]['nodes']    # dictionary of length 2632
    tissue = []
    cell_type = []
    tissues = ['kidney', 'thymus', 'spleen', "liver", 'lymph', 'stomach', 'heart',
               'small intestine', 'large intestine', 'blood', "brain", "thyroid",
               'placenta', 'eye', 'other', "embryo", "muscle"]

    identifiers = [['kidney'], ['thymocyte'], ['splenic'], ["hepat", "liver"],
                   ['T', 'B'], ['stomach', 'gastric'], ['heart', 'cardiac'],
                   ['small intestine'], ['large intestine'], [
                       'blood'], ["brain"], ["thyroid"],
                   ['placenta'], ['eye', 'retina'], ['other'], ["embryo"], ["muscle", 'myo']]

    for i in cells:
        f = 0
        if ('lbl' in i.keys()):
            for j in range(len(identifiers)):
                for k in identifiers[j]:
                    if ((k in i['lbl']) and f == 0):
                        tissue.append(tissues[j])
                        f = 1
            if (f == 0):
                tissue.append('other')
            cell_type.append(i['id'].split('/')[-1]+': '+i['lbl'])
        else:
            tissue.append('other')
            cell_type.append(i['id'].split('/')[-1])
    dic = {}
    for i in tissue:
        dic[i] = []
    for i in range(len(cell_type)):
        dic[tissue[i]].append(cell_type[i])
    return dic


def get_neighbors(
        x: Union[AnnData, np.ndarray, list],
        n_neighbors: int = 5,
        metric: str = 'euclidean',
        inplace: Optional[bool] = True,
        **kwargs):
    """
    Calculate the neighbors based on the dimension reduced data

    Parameters
    __________
    x: AnnData object containing the data.

    metric: The kind of metric used to calculate neighbors

    n_neighbors: the number of closest points as neighbors

    inplace: If set to true will update x.obs['uncertainty'], otherwise,
        return a new AnnData object.

    **kwargs: Additional parameters that will get passed to the
        clustering object as specified in method. For a full list
        see the documentation of the corresponding method.

    Returns
    _______
    If x is an AnnData object, will either return an AnnData object
    or None depending on the value of inplace.
    """
    # Validations
    is_AnnData = isinstance(x, AnnData)
    if not is_AnnData:
        raise InvalidArgument("x is not in AnnData format.")
    adata = x.copy() if not inplace else x

    n_graph = kneighbors_graph(adata.obsm['x_emb'], n_neighbors=n_neighbors,
                               include_self=False, metric=metric,
                               mode='distance').toarray()

    nn = n_graph.nonzero()
    indices = nn[1].reshape(-1, n_neighbors)
    adata.obsm['neighbor_labels'] = adata.obs['labels'][indices]

    distances = [n_graph[nn[0][i], nn[1][i]]
                 for i in range(len(x)*n_neighbors)]
    distances = np.array(distances).reshape(-1, n_neighbors)
    adata.obsm['neighbor_dist'] = distances

    if not inplace:
        return adata


def uncertainty(
        x: Union[AnnData, np.ndarray, list],
        method: str = 'conformal',
        n_neighbors: int = 5,
        inplace: Optional[bool] = True,
        **kwargs):
    """
    Calculate the uncertainty for each point

    Parameters
    __________
    x: AnnData object containing the data.

    method: The method used to calculate uncertainty

    n_neighbor: the number of closest points as neighbors

    inplace: If set to true will update x.obs['uncertainty'], otherwise,
        return a new AnnData object.

    **kwargs: Additional parameters that will get passed to the
        clustering object as specified in method. For a full list
        see the documentation of the corresponding method.

    Returns
    _______
    If x is an AnnData object, will either return an AnnData object
    or None depending on the value of inplace.
    """

    # Validations
    is_AnnData = isinstance(x, AnnData)
    if not is_AnnData:
        raise InvalidArgument("x is not in AnnData format.")

    adata = x.copy() if not inplace else x

    if 'labels' not in x.obs:
        raise InvalidArgument(
            "Cannot calculate uncertainty for dataset without labels")

    if 'neighbor_labels' not in x.obsm:
        get_neighbors(x, n_neighbors=n_neighbors)

    # uncertainty calculation:
    uncertainty = []

    if (method == "conformal"):
        threshold = 1  # make sure 'b' is larger than 0
        # x.obsm['neighbor_dist'].sort()
        dist = x.obsm['neighbor_dist']
        n_labels = x.obsm['neighbor_labels']
        labels = np.array(x.obs['labels'])
        nonconformity = []
        for i in range(len(np.unique(labels))):
            # same_label=(n_labels.transpose()==labels.transpose()).transpose()
            same_label = (n_labels == i)
            diff_label = (n_labels != i)
            a = (dist*same_label).sum(axis=1)
            b = (dist*diff_label).sum(axis=1)
            nonconformity.append(a/(b+threshold))
        nonconformity = np.array(nonconformity).transpose()
        # argsort twice to get the ranking
        order = nonconformity.argsort(axis=1)
        # the reversed ranking of nonconformity score
        ranking = order.argsort(axis=1)
        p = ranking/(len(np.unique(labels))+1)  # p_value for each label
        sp = p
        sp = sp[range(len(sp)), labels]  # selected p
        uncertainty = 1-sp
        uncertainty = np.around(uncertainty, 3)
        # calculate possible labels
        mask = (np.ones((len(labels), len(np.unique(labels)))) == 1)
        mask[range(len(labels)), labels] = False
        p = p*mask  # p value expect the selected ones

        p = p*(p > (np.average(p, axis=1)+np.std(p, axis=1)
                    ).reshape(len(p), -1))  # only leave entry>ave+std
        p = p/(p.sum(axis=1).reshape(len(p), -1))  # normalization
        order = p.argsort(axis=1)
        rank = order.argsort(axis=1)
        rank = rank[:, ::-1]
        p.sort(axis=1)
        p = p[:, ::-1]
        p = np.around(p, 3)
        # adata.uns['pos_labels']=rank
        # adata.uns['label_prob']=p

    elif (method == 'confidence'):
        for i in range(len(adata.obs.labels)):
            u, indices = np.unique(adata.obsm['neighbor_labels'][i],
                                   return_inverse=True)
            knn_label_freq = np.bincount(indices).max()
            confidence = knn_label_freq / n_neighbors
            uncertain_score = (1 - confidence) * 100
            uncertainty.append(uncertain_score)
    else:
        raise InvalidArgument("Invalid uncertainty method")

    adata.obs['uncertainty'] = uncertainty

    # make different labels for certain and uncertain points
    # threshhold for defining uncertain points
    t = np.average(uncertainty)+np.std(uncertainty)
    uncertain_matrix = (uncertainty > t)
    certainty = []
    t = np.average(uncertainty)+np.std(uncertainty)
    for i in range(len(labels)):
        if (uncertainty[i] > t):
            pos_label = ""
            for j in range(len(np.unique(labels))):
                prob = p[i, j]
                if (prob > 0):
                    lb = str(rank[i, j])
                    pos_label = pos_label+lb+": "+str(prob)+"\n"
            certainty.append(
                "Uncertain point, other possible label(s):\n"+str(pos_label))
    adata.uns['uncertainty_text'] = certainty

    # adata.obs['uncertain_matrix']=uncertain_matrix
    labels = np.array(adata.obs['labels'])
    uncertain_labels = labels*2+uncertain_matrix  # certain labels are even numbers,
    # uncertain labels are odd numbers
    # origin labels = uncertain labels -1 / 2
    #               = certain labels / 2
    # adata.obs['uncertain_labels']=uncertain_labels

    # make different subsets for certain and uncertain points
    clusters = list(adata.uns['subsets'].keys())
    uncertain_clusters = {}
    for i in clusters:
        uncertain_clusters[i+'_certain'] = np.array([], dtype='int64')
        uncertain_clusters[i+'_uncertain'] = np.array([], dtype='int64')
    for i in range(len(labels)):  # put points into certain/uncertain clusters
        if uncertain_labels[i] % 2 == 0:  # even, certain points
            cluster = list(adata.uns['subsets'].keys())[
                int(uncertain_labels[i]/2)]  # the cluster it belongs to
            uncertain_clusters[cluster+'_certain'] = np.append(
                uncertain_clusters[cluster+'_certain'], i)
        else:  # odd, uncertain points
            cluster = list(adata.uns['subsets'].keys())[
                int((uncertain_labels[i]-1)/2)]  # the cluster it belongs to
            uncertain_clusters[cluster+'_uncertain'] = np.append(
                uncertain_clusters[cluster+'_uncertain'], i)

    adata.uns['uncertain_subsets'] = uncertain_clusters

    if not inplace:
        return adata
