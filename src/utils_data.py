import numpy as np
from sklearn.datasets import load_iris, load_digits

def read_iris():
    X = load_iris()
    return X.data, X.target

def read_data():
    return np.random.random((20, 32, 32))