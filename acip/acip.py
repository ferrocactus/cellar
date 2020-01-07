import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import json
import tqdm
import time

# Visualization and Dimensionality Reduction
from umap import UMAP
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
# Clustering
from sklearn.cluster import KMeans
from sklearn.cluster import SpectralClustering
from acip.k_medoids import KMedoids
# Metrics
from sklearn.metrics import mean_squared_error as mse, silhouette_score
# Utils
from utils.utils_visualization import reduce_and_plot
from utils.utils_experiment import read_config

class ACIP:
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
        print_header("Visualizing")
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
        Reduces the dimensionality of the data and stores it in self.x_train_emb.
        """
        print_header("Reducing dimensionality")
        print("Using " + method + ".")
        if self.config is not None:  # Use the provided config file instead
            kwargs = read_config(self.config)['reduce_dim'][method]
            if self.verbose:
                print('Using the following params:')
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
    
    def cluster(self, method='kmedoids', n_clusters_list=list(range(2, 17, 2)), **kwargs):
        """
        Clusters the embeddings as constructed by reduce_dim.
        """
        print_header("Clustering")
        print("Using " + method + ".")
        if self.config is not None:
            kwargs = read_config(self.config)['cluster'][method]
            if self.verbose:
                print("Using the following params:")
                pretty_print_dict(kwargs)

        silhouette_score_list = []
        max_sscore = -2
        
        time.sleep(0.2) # To avoid tqdm bar from appearing before prints
        pbar = tqdm.tqdm(n_clusters_list)
        for n_clusters in pbar:
            pbar.set_description("Trying n_clusters=" + str(n_clusters))
            if method == 'kmeans':
                kmeans = KMeans(n_clusters=n_clusters, **kwargs)
                y_train_pred = kmeans.fit_predict(self.x_train_emb)
            elif method == 'spectral':
                spectral = SpectralClustering(n_clusters=n_clusters, **kwargs)
                y_train_pred = spectral.fit_predict(self.x_train_emb)
            elif method == 'kmedoids':
                kmedoids = KMedoids(n_clusters=n_clusters, **kwargs)
                y_train_pred = kmedoids.fit_predict(self.x_train_emb)
            else:
                raise NotImplementedError()

            sscore = silhouette_score(self.x_train_emb, y_train_pred)
            silhouette_score_list.append(sscore)
            if max_sscore < sscore:
                max_sscore = sscore
                self.n_clusters = n_clusters
                self.y_train_pred = y_train_pred
        
        print("Clustering complete.")
        print("Highest silhouette score is achieved for n_clusters =", self.n_clusters)
        
        if self.verbose:
            sns.set_style('whitegrid')
            fig, ax = plt.subplots()
            ax.plot(n_clusters_list, silhouette_score_list)
            ax.set_xlabel('Number of Clusters')
            ax.set_ylabel('Silhouette Score')
            sns.despine()
            fig.set_size_inches(10, 5)
    
    def flow(self, reduce_dim='pca', cluster='kmedoids', reduce_plot='umap'):
        self.reduce_dim(method=reduce_dim)
        self.cluster(method=cluster)
        self.reduce_plot(method=reduce_plot)
    
    def get_markers(self):
        df = pd.DataFrame(self.x_train)
        df['y_train_pred'] = self.y_train_pred
        self.markers = df.groupby('y_train_pred').median().to_numpy()
        return self.markers

def pretty_print_dict(mydict):
    print(json.dumps(mydict, sort_keys=True, indent=4))

def print_header(text, max_length=50, prepend_newline=True):
    if len(text) >= max_length:
        if prepend_newline:
            print()
        print("*** " + text + " ***")
    else:
        rems = int((max_length - len(text)) / 2) - 2
        if prepend_newline:
            print()
        print('*' * rems, text, '*' * (max_length - rems - len(text) - 2))