import numpy as np
import pandas as pd
from mnist import MNIST
from sklearn.datasets import load_iris, load_digits

def load_data(dataset):
    if dataset == 'spellman':
        return pd.read_csv('datasets/Spellman.csv', index_col=0).to_numpy(), []
    if dataset == 'iris':
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