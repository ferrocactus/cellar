import numpy as np
import pandas as pd
from matplotlib import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import tqdm
import json

# Visualization and Dimensionality Reduction
from umap import UMAP
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
# Clustering
from sklearn.cluster import KMeans
from sklearn.cluster import SpectralClustering
from sklearn.neighbors import kneighbors_graph
from src.k_medoids import KMedoids

from sklearn.metrics import mean_squared_error as mse

from src.utils_visualization import reduce_and_plot
from src.utils_experiment import read_config

class Workflow:
    def __init__(self, x, y=None, config=None, verbose=False):
        self.set_train_data(x, y)
        # If use_config is set for a function, all its kwargs will be ignored
        self.config         = config
        self.verbose        = verbose
    
    def set_train_data(self, x, y=None, row_names=None, col_names=None):
        assert len(x.shape) == 2, "Data needs to be of shape (n x d)."
        self.x_train        = x
        self.n_train        = x.shape[0]
        self.dims           = x.shape[1]
        self.y_train        = y
        self.row_names      = row_names
        self.col_names      = col_names
    
    def set_row_names(self, row_names):
        self.row_names      = row_names
    
    def set_col_names(self, col_names):
        self.col_names      = col_names

    def set_config(self, config):
        self.config         = config
    
    def normalize(self):
        raise NotImplementedError()
    
    def boxplot(self, data='train'):
        """
        Boxplot of the columns of the data
        """
        if self.y_train is not None and self.verbose:
            print("WARNING: This is a boxplot of the columns of the data. " +
                    "Results may not be what expected.")
        sns.set_style("whitegrid")
        fig, ax = plt.subplots()
        if data == 'train':
            sns.boxplot(data=self.x_train)
        elif data == 'test' and self.has_test_data:
            sns.boxplot(data=self.x_test)
        else:
            raise NameError("Dataset not found.")
        sns.despine(left=True)
        if self.col_names is not None and len(self.col_names) < 40:
            plt.xticks(np.arange(self.dims), self.col_names, rotation=90)
        plt.xlabel('Columns')
        plt.ylabel('Rows')
        plt.title('Boxplot')
        fig.set_size_inches(10, 5)
    
    def reduce_plot(self, labels=None, method='umap', use_emb=True, **kwargs):
        if labels is None:      # By default use predicted labels
            labels = self.y_train_pred
        if self.config is not None:
            kwargs = read_config(self.config)['reduce_plot'][method]
            if self.verbose:
                print('Using', method, 'with the following params:')
                pretty_print_dict(kwargs)
        self.x_train_2d_emb = reduce_and_plot(
                        x=self.x_train_emb if use_emb else self.x_train,
                        y=labels,
                        method=method,
                        **kwargs)

    def pca_plot_var_ratio(self, n_components=None, **kwargs):
        """
        Plots the percentage of variance explained by each of the PCA components.
        """
        kwargs['n_components'] = n_components or min(self.dims, self.n_train)
        pca = PCA(**kwargs)
        pca.fit(self.x_train)

        # Plotting
        sns.set_style('whitegrid')

        fig, (ax1, ax2) = plt.subplots(1, 2)
        ax1.plot(list(range(1, kwargs['n_components'] + 1)),
                pca.explained_variance_ratio_)
        ax1.set_xlabel("Number of components")
        ax1.set_ylabel("Percentage of Explained Variance per Component")
        ax1.set_ylim(0)
        # Cumulative plot
        ax2.plot(list(range(1, kwargs['n_components'] + 1)),
                np.cumsum(pca.explained_variance_ratio_))
        ax2.set_xlabel("Number of components")
        ax2.set_ylabel("Cumulative Percentage of Explained Variance")
        ax2.set_ylim(0)

        sns.despine()
        fig.set_size_inches(10, 5)
    
    def reduce_dim(self, method='pca', **kwargs):
        """
        Reduces the dimensionality of the data and stores it in self.emb.
        """
        if self.config is not None:  # Use the provided config file instead
            kwargs = read_config(self.config)['reduce_dim'][method]
            if self.verbose:
                print('Using', method, 'with the following params:')
                pretty_print_dict(kwargs)

        if method == 'pca':
            pca = PCA(**kwargs)
            pca.fit(self.x_train)
            # Print scores
            self.x_train_emb = pca.transform(self.x_train)
            print("Embedding created. Train MSE:",
                    mse(self.x_train, pca.inverse_transform(self.x_train_emb)))
            #print("Train Average Log Likelihood:", pca.score(self.x_train))
        elif method == 'umap':
            umap = UMAP(**kwargs)
            umap.fit(self.x_train)
            self.x_train_emb = umap.transform(self.x_train)
        else:
            raise NotImplementedError()
    
    def cluster(self, method='kmeans', **kwargs):
        """
        Clusters the embeddings as constructed by reduce_dim.
        """
        if self.config is not None:
            kwargs = read_config(self.config)['cluster'][method]
            if self.verbose:
                print("Using", method, "clustering with the following params:")
                pretty_print_dict(kwargs)

        if method == 'kmeans':
            kmeans = KMeans(**kwargs)
            kmeans.fit(self.x_train_emb)
            self.y_train_pred = kmeans.predict(self.x_train_emb)
        elif method == 'spectral':
            spectral = SpectralClustering(**kwargs)
            spectral.fit(self.x_train_emb)
            self.y_train_pred = spectral.labels_
        elif method == 'kmedoids':
            kmedoids = KMedoids(**kwargs)
            kmedoids.fit(self.x_train_emb)
            self.y_train_pred = kmedoids.predict(self.x_train_emb)
        else:
            raise NotImplementedError()
        
        if self.verbose:
            print("Clustering complete.")
    
    def get_markers(self):
        df = pd.DataFrame(self.x_train)
        df['y_train_pred'] = self.y_train_pred
        self.markers = df.groupby('y_train_pred').median().to_numpy()
        return self.markers

def pretty_print_dict(mydict):
    print(json.dumps(mydict, sort_keys=True, indent=4))