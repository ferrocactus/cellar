from ast import literal_eval
import numpy as np

def read_config(dataset):
    with open("configs/" + dataset + ".json", "r") as f:
        config = literal_eval(f.read())
    return config

def load_data(dataset):
    # return X, Y
    if dataset == 'spleen':
        import anndata
        ann = anndata.read_h5ad('datasets/spleen/dim_reduced_clustered.h5ad')
        return ann.X, ann.obs['leiden'].to_numpy().astype(np.int), ann.var.index.to_numpy()
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
        return X.data, X.target
    elif dataset == 'mnist':
        from mnist import MNIST
        mndata = MNIST('datasets/mnist')
        images, labels = mndata.load_training()
        return np.asarray(images), np.asarray(labels)
    else:
        raise FileNotFoundError('Dataset not supported.')