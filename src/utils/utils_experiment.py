from ast import literal_eval
import pandas as pd
import numpy as np
import json
import re
import itertools
from functools import reduce

def read_config(dataset):
    with open("configs/" + dataset + ".json", "r") as f:
        config = literal_eval(f.read())
    return config

def parse(x):
    parsed_x = np.char.split(x.flatten(), sep='.', maxsplit=1)
    parsed_x = np.array([i[0] for i in parsed_x])
    parsed_x = np.char.replace(parsed_x, '-', '')
    parsed_x = np.char.replace(parsed_x, ' ', '')
    # Remove everything that goes inside brackets
    parsed_x = [re.sub(r" ?\([^)]+\)", "", item) for item in parsed_x]
    parsed_x = np.char.upper(parsed_x)
    parsed_x = parsed_x.reshape(x.shape)
    return parsed_x

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