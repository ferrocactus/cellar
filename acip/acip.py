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
from sklearn.cluster import AgglomerativeClustering
from acip.k_medoids import KMedoids
# Differential Expression
from scipy import stats
# Metrics
from sklearn.metrics import mean_squared_error as mse, silhouette_score
# Utils
from utils.utils_visualization import reduce_and_plot
from utils.utils_experiment import read_config
# Constrained clustering
from copkmeans.cop_kmeans import cop_kmeans
from active_semi_clustering.semi_supervised.pairwise_constraints import PCKMeans

class ACIP:
    def __init__(self, x, y=None, config=None, verbose=False):
        self.set_train_data(x, y)
        # If use_config is set for a function, all its kwargs will be ignored
        self.config = config
        self.verbose = verbose
    
    def set_train_data(self, x, y=None, row_names=None, col_names=None):
        assert len(x.shape) == 2, "Data needs to be of shape (n x d)."
        self.x_train = x
        self.n_train = x.shape[0]
        self.dims = x.shape[1]
        self.y_train = y
        self.row_names = row_names
        self.col_names = col_names
    
    def set_row_names(self, row_names):
        self.row_names = row_names
    
    def set_col_names(self, col_names):
        self.col_names = col_names

    def set_config(self, config):
        self.config = config

    def set_verbose(self, verbose):
        self.verbose = verbose
    
    def normalize(self):
        raise NotImplementedError()
    
    def boxplot(self):
        """
        Boxplot of the columns of the data
        """
        if self.y_train is not None and self.verbose:
            print("WARNING: This is a boxplot of the columns of the data. " +
                    "Results may not be what expected.")
        sns.set_style("whitegrid")
        fig, ax = plt.subplots()
        sns.boxplot(data=self.x_train)

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

    def filter_genes(self):
        """
        Computes the variance for every gene (column) and picks the highest 20%
        variable genes. Stores the new matrix into self.x_train_filtered
        """
        print_header("Filtering genes")
        gene_variances = np.var(self.x_train, axis=0)
        highest_var_gene_indices = gene_variances.argsort()[-int(0.2 * self.x_train.shape[1]):]
        self.x_train_filtered = self.x_train[:, highest_var_gene_indices]

        # Plotting
        if self.verbose:
            fig, ax = plt.subplots()
            ax.hist(gene_variances)
            ax.set_xlabel("Variance")
            ax.set_ylabel("Gene Count")
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
            pca.fit(self.x_train_filtered)
            # Print scores
            self.x_train_emb = pca.transform(self.x_train_filtered)
            print("Embedding created. MSE:",
                    mse(self.x_train_filtered, pca.inverse_transform(self.x_train_emb)))
            print("Used", pca.n_components_, "components.")
            #print("Train Average Log Likelihood:", pca.score(self.x_train_filtered))
        else:
            raise NotImplementedError()
    
    def cluster(self, method='kmedoids', n_clusters_list=list(range(2, 15, 2)), **kwargs):
        """
        Clusters the embeddings as constructed by reduce_dim.
        """
        print_header("Clustering")
        print("Using " + method + ".")
        if self.config is not None:
            kwargs = read_config(self.config)['cluster'][method]
            if 'n_clusters' in kwargs:
                n_clusters_list=[kwargs['n_clusters']]
                kwargs.pop('n_clusters', None)
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

            sscore = silhouette_score(self.x_train_emb, y_train_pred, metric='correlation')
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
    
    def get_markers(self):
        self.markers = []

        fig, ax = plt.subplots(self.n_clusters, 2)

        for cluster_id in tqdm.tqdm(range(self.n_clusters), desc="Cluster id"):
            x_in = self.x_train[np.where(self.y_train_pred == cluster_id)]
            x_not_in = self.x_train[np.where(self.y_train_pred != cluster_id)]

            pvals = []
            mads = []
            for gene in range(self.dims):
                gene_vector_in = x_in[:, gene]
                gene_vector_not_in = x_not_in[:, gene]

                t, pval = stats.ttest_ind(gene_vector_in, gene_vector_not_in, equal_var=False)
                pvals.append(pval)
                mads.append(np.mean(gene_vector_in) - np.mean(gene_vector_not_in))
            marker = self.col_names[np.where((np.array(pvals) < (0.05 / self.dims)) & (np.array(mads) > 1))]
            self.markers.append(marker)

            if self.verbose:
                ax[cluster_id][0].hist(pvals)
                ax[cluster_id][0].set_xlabel("p-values")
                ax[cluster_id][0].set_ylabel("gene count")
                ax[cluster_id][0].set_title("Cluster:" + str(cluster_id))
                ax[cluster_id][1].hist(mads, color='r')
                ax[cluster_id][1].set_xlabel("absolute difference")
                ax[cluster_id][1].set_ylabel("gene count")
                ax[cluster_id][1].set_title("Cluster:" + str(cluster_id))

        self.markers = np.asarray(self.markers)
        
        if self.verbose:
            sns.despine()
            fig.set_size_inches(10, self.n_clusters * 5)
    
    def cluster_constraints(self, method='agglomerative', n_clusters_list=list(range(2, 65, 2)), constraints=None, **kwargs):
        assert constraints is not None, "No constraints provided."
        print_header("Clustering with constraints")
        print("Using " + method + ".")

        silhouette_score_list = []
        max_sscore = -2
        
        time.sleep(0.2) # To avoid tqdm bar from appearing before prints
        pbar = tqdm.tqdm(n_clusters_list)
        for n_clusters in pbar:
            pbar.set_description("Trying n_clusters=" + str(n_clusters))
            if method == 'cop_kmeans':
                y_train_pred, centers = cop_kmeans(dataset=self.x_train_emb, k=n_clusters,
                                                    ml=constraints['ml'], cl=constraints['cl'],
                                                    **kwargs)
            elif method == 'agglomerative':
                ac = AgglomerativeClustering(n_clusters=n_clusters, connectivity=constraints, **kwargs)
                y_train_pred = ac.fit_predict(self.x_train_emb)
            elif method == 'pckmeans':
                pck = PCKMeans(n_clusters=self.n_clusters)
                pck.fit(self.x_train_emb, ml=constraints['ml'], cl=constraints['cl'])
                y_train_pred = pck.labels_
            else:
                raise NotImplementedError()
            
            sscore = silhouette_score(self.x_train_emb, y_train_pred, metric='correlation')
            silhouette_score_list.append(sscore)
            if max_sscore < sscore:
                max_sscore = sscore
                self.n_clusters = n_clusters
                self.y_train_pred = y_train_pred
        
        print("Clustering with constraints complete.")
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
        self.filter_genes()
        self.reduce_dim(method=reduce_dim)
        self.cluster(method=cluster)
        self.reduce_plot(method=reduce_plot)
        self.get_markers()

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