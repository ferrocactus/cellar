import numpy as np
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

from sklearn.metrics import mean_squared_error as mse

from src.utils_visualization import reduce_and_plot

class Workflow:
    def __init__(self, x, y=None, config=None, verbose=False):
        self.set_train_data(x, y)
        self.has_test_data = False
        # If use_config is set for a function, all its kwargs will be ignored
        self.config = config
        self.verbose = verbose
    
    def set_train_data(self, x, y=None, row_names=None, col_names=None):
        assert len(x.shape) == 2, "Data needs to be of shape (n x d)."
        self.x_train        = x
        self.n_train        = x.shape[0]
        self.dims           = x.shape[1]
        self.y_train        = y
        self.row_names      = row_names
        self.col_names      = col_names
    
    def set_test_data(self, x, y=None):
        assert len(x.shape) == 2, "Data needs to be of shape (n x d)."
        assert x.shape[1]   == self.dims, "Dimension mismatch between test & train."
        self.has_test_data  = True
        self.x_test         = x
        self.n_test         = x.shape[0]
        self.has_y_test     = False if y is None else True
        self.y_test         = y
    
    def set_row_names(self, row_names):
        self.row_names = row_names
    
    def set_col_names(self, col_names):
        self.col_names = col_names

    def set_config(self, config):
        self.config = config
    
    def assert_config(self):
        if self.config is None:
            assert False, "Config file not provided."
    
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
    
    def reduce_plot(self, labels=None, method='umap', use_config=True, **kwargs):
        if labels is None:      # By default use predicted labels
            labels = self.y_train_pred
        if use_config:
            self.assert_config()
            kwargs = self.config['reduce_plot'][method]
            if self.verbose:
                print('Using', method, 'with the following params:')
                pretty_print_dict(kwargs)
        reduce_and_plot(x=self.x_train, y=labels, method=method, **kwargs)

    def pca_plot_var_ratio(self, n_components=None, **kwargs):
        """
        Plots the percentage of variance explained by each of the PCA components.
        """
        kwargs['n_components'] = n_components or self.dims
        pca = PCA(**kwargs)
        pca.fit(self.x_train)

        # Plotting
        sns.set_style("whitegrid")

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
    
    def reduce_dim(self, method='pca', use_config=True, **kwargs):
        """
        Reduces the dimensionality of the data and stores it in self.emb.
        """
        if use_config:  # Use the provided config file instead
            self.assert_config()
            kwargs = self.config['reduce_dim'][method]
            if self.verbose:
                print('Using', method, 'with the following params:')
                pretty_print_dict(kwargs)

        if method == 'pca':
            pca = PCA(**kwargs)
            pca.fit(self.x_train)

            # Print scores
            self.x_train_emb        = pca.transform(self.x_train)
            self.x_train_emb_mse    = mse(self.x_train,
                                        pca.inverse_transform(self.x_train_emb))
            self.x_train_emb_score  = pca.score(self.x_train)
            print("Embedding created. Train MSE:", self.x_train_emb_mse)
            print("Train Average Log Likelihood:", self.x_train_emb_score)

            if self.has_test_data:
                self.x_test_emb         = pca.transform(self.x_test)
                self.x_test_emb_mse     = mse(self.x_test,
                                            pca.inverse_transform(self.x_test_emb))
                self.x_test_emb_score   = pca.score(self.x_test)
                print("Embedding created. Test MSE:", self.x_test_emb_mse)
                print("Test Average Log Likelihood:", self.x_test_emb_score)
        else:
            raise NotImplementedError()
    
    def cluster(self, method='kmeans', use_config=True, **kwargs):
        """
        Clusters the embeddings as constructed by reduce_dim.
        """
        if use_config:
            self.assert_config()
            kwargs = self.config['cluster'][method]
            if self.verbose:
                print('Using', method, 'with the following params:')
                pretty_print_dict(kwargs)

        if method == 'kmeans':
            kmeans = KMeans(**kwargs)
            kmeans.fit(self.x_train_emb)

            # Print scores
            self.y_train_pred      = kmeans.predict(self.x_train_emb)
            x_train_cluster_score  = kmeans.score(self.x_train_emb)
            print("Clustering complete. Train Score:", x_train_cluster_score)

            if self.has_test_data:
                self.y_test_pred       = kmeans.predict(self.x_test_emb)
                x_test_cluster_score   = kmeans.score(self.x_test_emb)
                print("Test Score:", x_test_cluster_score)
        elif method == 'spectral':
            spectral = SpectralClustering(**kwargs)
            spectral.fit(self.x_train_emb)
            self.y_train_pred = spectral.labels_
            print('Clustering complete.')
        else:
            raise NotImplementedError()

def pretty_print_dict(mydict):
    print(json.dumps(mydict, sort_keys=True, indent=4))