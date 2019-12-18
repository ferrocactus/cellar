import numpy as np
from mnist import MNIST
from sklearn.datasets import load_iris, load_digits

def load_data(dataset):
    if dataset == 'iris':
        X = load_iris()
        s = int(2/3 * X.data.shape[0])
        return X.data[:s], X.target[:s], X.data[s:], X.target[s:]
    elif dataset == 'digits':
        X = load_digits()
        s = int(2/3 * X.data.shape[0])
        return X.data[:s], X.target[:s], X.data[s:], X.target[s:]
    elif dataset == 'mnist':
        mndata = MNIST('datasets/mnist')
        images, labels = mndata.load_training()
        images, labels = np.asarray(images), np.asarray(labels)
        return images[:2000], labels[:2000], images[5000:7000], labels[5000:7000]
    else:
        raise FileNotFoundError('Dataset not supported.')