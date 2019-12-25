from ast import literal_eval
import numpy as np
import pandas as pd
from mnist import MNIST
import anndata
from sklearn.datasets import load_iris, load_digits

def read_config(dataset):
    with open("configs/" + dataset + ".json", "r") as f:
        config = literal_eval(f.read())
    return config

def load_data(dataset):
    if dataset == 'spleen':
        ann = anndata.read_h5ad('datasets/spleen/dim_reduced_clustered.h5ad')
        return ann.X, ann.obs['leiden'].to_numpy().astype(np.int)
    if dataset == 'spellman':
        return pd.read_csv('datasets/Spellman.csv', index_col=0).to_numpy(), None
    elif dataset == 'iris':
        X = load_iris()
        return X.data, X.target
    elif dataset == 'digits':
        X = load_digits()
        return X.data, X.target
    elif dataset == 'mnist':
        mndata = MNIST('datasets/mnist')
        images, labels = mndata.load_training()
        return np.asarray(images), np.asarray(labels)
    else:
        raise FileNotFoundError('Dataset not supported.')