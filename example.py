from src.workflow import Workflow
from src.utils_experiment import read_config, load_data

import numpy as np
from matplotlib import pyplot as plt

import warnings
warnings.filterwarnings('ignore')


if __name__ == '__main__':
    config = read_config('digits')
    X, Y = load_data(config['dataset']['name'])
    plt.ion()

    w = Workflow(X, Y)
    w.boxplot()
    w.pca_plot_var_ratio()
    w.reduce_dim(method='pca', **config['reduce_dim']['pca'])
    w.cluster(method='spectral', **config['cluster']['spectral'])
    w.reduce_plot(labels=w.y_train_pred, method='umap', **config['reduce_plot']['umap'])

    plt.ioff()
    plt.show()
