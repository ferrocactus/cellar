from ast import literal_eval
import pandas as pd
import numpy as np
import json
from scipy.stats import hypergeom
import re
import itertools
from functools import reduce

def read_config(dataset):
    with open("configs/" + dataset + ".json", "r") as f:
        config = literal_eval(f.read())
    return config

def load_data(dataset):
    # return X, Y
    if dataset == 'spleen':
        import anndata
        ann = anndata.read_h5ad('datasets/spleen/dim_reduced_clustered.h5ad')
        return ann.X, ann.obs['leiden'].to_numpy().astype(np.int), ann.var.index.to_numpy().astype('U')
    if dataset == 'spellman':
        import pandas as pd
        return pd.read_csv('datasets/Spellman.csv', index_col=0).to_numpy(), None
    elif dataset == 'iris':
        from sklearn.datasets import load_iris
        X = load_iris()
        return X.data, X.target
    elif dataset == 'digits':
        from sklearn.datasets import load_digits
        X = load_digits()
        return X.data, X.target, np.array([str(i) for i in range(1, 11)])
    elif dataset == 'mnist':
        from mnist import MNIST
        mndata = MNIST('datasets/mnist')
        images, labels = mndata.load_training()
        return np.asarray(images), np.asarray(labels)
    else:
        raise FileNotFoundError('Dataset not supported.')

def gene_id_to_name(gene_ids, path="gene_id_name.csv"):
    """
    Given an array of gene ids, convert them to gene names using the provided file
    and return an array of the same shape with the converted names.
    """
    gene_dict = pd.read_csv('markers/' + path, index_col=0, squeeze=True).to_dict()
    # Split gene ids by dot char if there is any and take the first part
    splitted = np.char.split(gene_ids.flatten(), sep='.')
    return np.char.upper(np.vectorize(gene_dict.get)(
        np.array([i[0] for i in splitted]).reshape(gene_ids.shape))
    )

def read_markers(path="CellTypeMarker.json"):
    with open("markers/" + path) as f:
        markers = json.load(f)

    marker_organized = {}
    for key in markers:
        organ = {}
        for cell in markers[key]:
            group = []
            for tp in markers[key][cell]:
                pos = list(filter(lambda x: x!='', markers[key][cell][tp]['positive']))
                neg = list(filter(lambda x: x!='', markers[key][cell][tp]['negative']))
                if len(pos) > 0:
                    group += pos
                if len(neg) > 0:
                    group += neg
            group = [re.sub(r'[^A-Za-z0-9]+', '', i.upper()) for i in group]
            organ[cell] = group
        marker_organized[key] = organ
    return marker_organized

def gene_name_to_cell(gene_names, path="CellTypeMarker.json"):
    with open("markers/" + path, "r") as f:
        markers_dict = json.load(f)
    gene_names = [[re.sub(r'[^A-Za-z0-9]+', '', i.upper()) for i in j] for j in gene_names]
    return identify_genes(gene_names, markers_dict)

def identify_genes(gene_lists, markers_dict):
    """
    Given an array of arrays of significant genes in "genes_list"
    for every such array of genes, find the best matching marker population
    and subpopulation in "markers_dict". "markers_dict" is assumed to be
    a dictionary of dictionaries of a single list.
    params:
        genes_list: 2D array
        markers_dict: {organ : {cell1: [], cell2: []}}
    """
    # Group subpopulations
    pops = [reduce(lambda a, b: a+b, markers_dict[markers].values())
                    for markers in markers_dict]
    gene_lists = [[gene.upper() for gene in gene_list] for gene_list in gene_lists]
    h, svs, intersec = find_populations(gene_lists, pops)
    names = np.array(list(markers_dict.keys()))
    pop_names = np.where(h >= 0, names[h], "None")

    # Find subpopulation if any
    sub_pop_names = []
    sub_svs = []
    sub_intersec = []
    for i, gene_list in enumerate(gene_lists):
        if h[i] >= 0:
            subpops = list(markers_dict[list(markers_dict.keys())[h[i]]].values())
            h_i, sv_i, intersec_i = find_population(gene_list, subpops)
            subnames = list(markers_dict[list(markers_dict.keys())[h[i]]].keys())
            sub_pop_names.append(subnames[h_i]) # we are guaranteed that h_i >= 0
            sub_svs.append(sv_i)
            sub_intersec.append(intersec_i)
        else:
            sub_pop_names.append("None")
            sub_svs.append(1)
            sub_intersec.append([])
    return pop_names, svs, intersec, sub_pop_names, sub_svs, sub_intersec

def find_populations(X, pops):
    """
    Compute for every x in X, the pop in pops where x is most likely
    to have been drawn from by using a hypergeometric test.
    params:
        X: a list of 1D arrays
        pops: a list of 1D arrays
        return_intersec: to return the intersecion or no
    returns:
        (list, list, [list]): where the first list
                contains the indices of the best pop for x
                and the second list contains the survival
                value (1 - cdf) for the returned pop.
                If return_intersec is set to True, the returned
                list will contain an extra list which contains
                a list of the intersection between x and pop
    """
    h, svs, intersec = zip(*[find_population(x, pops) for x in X])
    return np.array(h), np.array(svs), np.array(intersec)

def find_population(x, pops):
    """
    See find_populations. Assumes x is a single list.
    """
    M = sum(map(lambda pop: len(pop), pops))
    N = len(x)
    svs = []
    for pop in pops:
        n = len(pop)
        k = len(np.intersect1d(x, pop))
        sv = hypergeom.sf(k-1, M=M, n=n, N=N) if k > 0 else 1
        svs.append(sv)
    h = np.argmin(svs)
    intersec = np.intersect1d(x, pops[h])
    if len(intersec) == 0: # in case of no intersection, return -1
        h = -1
    return h, svs[h], intersec