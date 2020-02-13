from ast import literal_eval
import pandas as pd
import numpy as np
import json

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
    return np.char.upper(np.vectorize(gene_dict.get)(np.array([i[0] for i in splitted]).reshape(gene_ids.shape)))

def gene_name_to_cell(gene_names, path="markers.json"):
    """
    Given an array of arrays of gene names, find the biggest overlap with
    gene names given in file, and determine cell type. Returns an array
    of size len(gene_names).
    """
    with open("markers/" + path) as f:
        markers = json.load(f)
    level1 = ["None"]
    organ_names = ["None"]
    for i, organ in enumerate(markers):
        organ_genes = []
        for cell in markers[organ]:
            cell_genes = []
            for specialty in markers[organ][cell]:
                cell_genes += markers[organ][cell][specialty].get('positive')
                cell_genes += markers[organ][cell][specialty].get('negative')
            cell_genes = list(filter(None, cell_genes)) # filter empty strings
            if len(cell_genes) > 0:
                organ_genes += cell_genes
        if len(organ_genes) > 0:
            level1.append(np.char.upper(np.array(organ_genes, dtype='U')))
            organ_names.append(organ)
    organ_names = np.array(organ_names, dtype='U')
    top_classes, common_genes = match(gene_names, level1)
    return organ_names[top_classes], common_genes

def match(X, Y):
    """
    Given two arrays of arrays X and Y, find for each array x of X
    the array y of Y which has most elements in common. Return top
    indices found in Y.
    """
    top_classes = np.full(X.shape[0], 0)
    common_genes = (X.shape[0]) * [[]]
    for i, x in enumerate(X):
        best_len_cap = 0
        for j, y in enumerate(Y):
            cap = np.intersect1d(x, y)
            if len(cap) > best_len_cap:
                best_len_cap = len(cap)
                top_classes[i] = j
                common_genes[i] = cap
    return top_classes, common_genes