import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import json
import tqdm
import time
import sys
import configparser
from ast import literal_eval

# Visualization and Dimensionality Reduction
from umap import UMAP
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from kneed import KneeLocator
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
from utils.utils_visualization import ( plot_2d,
                                        plot_explained_variance,
                                        plot_gene_variances,
                                        plot_scores,
                                        plot_marker_hist,
                                        plot_top_markers)
from utils.utils_experiment import (read_config,
                                    gene_id_to_name,
                                    gene_name_to_cell,
                                    dict_literal_eval)
# Constrained clustering
#from copkmeans.cop_kmeans import cop_kmeans
#from active_semi_clustering.semi_supervised.pairwise_constraints import PCKMeans

class ACIP:
    def __init__(self, x, config="defaults", verbose=False, row_ids=None, col_ids=None):
        self.set_train_data(x)
        self.parse_config(config)
        self.verbose = verbose
        self.set_ids(row_ids, col_ids)
        self.x_emb = None # Embeddings after running PCA
        self.y_pred = None # Labels after running clustering
        self.marker_ids = None # Markers for cluster
        self.x_2d_emb = None # Embeddings used for visualization

    def set_train_data(self, x):
        assert len(x.shape) == 2, "Data needs to be of shape (n x d)."
        self.x = x
        self.rownum = x.shape[0]
        self.colnum = x.shape[1]

    def set_ids(self, row_ids=None, col_ids=None):
        self.row_ids = row_ids.astype('U') if row_ids is not None else None
        self.col_ids = col_ids.astype('U') if col_ids is not None else None

    def parse_config(self, filename='configs/config.ini'):
        configp = configparser.ConfigParser()
        configp.read(filename)

        self.config = configp._sections
        for key in self.config:
            self.config[key] = dict_literal_eval(self.config[key])

        # Read from config file
        self.method_reduce_dim = self.config['methods'].get('reduce_dim', 'PCA')
        assert(self.method_reduce_dim in ['PCA'], self.method_reduce_dim + " not supported.")
        if self.method_reduce_dim not in self.config:
            self.config[self.method_reduce_dim] = {}

        self.method_cluster = self.config['methods'].get('cluster', 'KMedoids')
        assert(self.method_cluster in ['KMedoids', "KMeans", "SpectralClustering"], self.method_cluster + " not supported.")
        if self.method_cluster not in self.config:
            self.config[self.method_cluster] = {}

        self.method_cluster_eval = self.config['methods'].get('cluster_eval', 'silhouette_score')
        assert(self.method_cluster_eval in ['silhouette_score'], self.method_cluster_eval + " not supported.")
        if self.method_cluster_eval not in self.config:
            self.config[self.method_cluster_eval] = {}

        self.method_visualization = self.config['methods'].get('visualization', 'UMAP')
        assert(self.method_visualization in ['UMAP', 'TSNE'], self.method_visualization + " not supported.")
        if self.method_visualization not in self.config:
            self.config[self.method_visualization] = {}

        self.markers = 'markers'

    def normalize(self):
        raise NotImplementedError()

    def print_verbose(self, method):
        # Prints configs for given processing step. Reads from self.config
        if self.verbose == True:
            print("Step:", method)
            pretty_print_dict(self.config[method])

    def get_visual_emb(self):
        """
        Constructs 2d embeddings for visualization. Methods available: UMAP, TSNE.
        """
        if self.x_2d_emb is None:
            if self.x_emb is None:
                self.cluster()
            print_header("Constructing visual embeddings.")
            method = self.method_visualization
            self.print_verbose(method)
            # Use the method specified in the config file
            method_obj = globals()[method](**self.config[method])
            self.x_2d_emb = method_obj.fit_transform(self.x_emb, y=self.y_pred)

    def find_knee(self):
        """
        Runs PCA with evr_components number of components
        and uses KneeLocator to find the "knee" in the graph, which is the
        point where the graph turns the fastest. Use that value for running
        PCA when calling reduce_dim.
        """
        pca = PCA(self.config['PCA']['evr_components'])
        pca.fit(self.x)
        kn = KneeLocator(list(range(1, pca.n_components_ + 1)),
                        pca.explained_variance_ratio_,
                        curve='convex',  # not sctrict, but we can assume it
                        direction='decreasing')  # by virtue of PCA
        self.knee = kn.knee
        print("Ankle found at", self.knee, "components.")
        self.explained_variance = pca.explained_variance_ratio_
        return self.knee

    def reduce_dim(self):
        """
        Reduces the dimensionality of the data and stores it in self.x_emb.
        Requires filter_genes to be run first.
        """
        print_header("Reducing dimensionality")
        self.print_verbose(self.method_reduce_dim)

        if self.method_reduce_dim == 'PCA':
            n_components = self.config['PCA'].get("n_components", "auto")
            if n_components == "auto":
                n_components = self.find_knee()
            pca = PCA(n_components)
            pca.fit(self.x)
            self.x_emb = pca.transform(self.x)

            if self.verbose:
                print("Embedding created. MSE:", mse(self.x, pca.inverse_transform(self.x_emb)))
            print("Used", pca.n_components_, "components.")
        else:
            raise NotImplementedError()

    def cluster(self):
        """
        Clusters the embeddings created by reduce_dim. Tries various n_clusters
        as specified by the list (low, high, step) and stores into y_pred the
        highest scoring labels according to eval_method.
        """
        if self.x_emb is None:
            self.reduce_dim()
        print_header("Clustering")
        self.print_verbose(self.method_cluster)

        # Construct range of clusters to try
        self.n_clusters_list = list(range(*self.config[self.method_cluster].get('n_clusters', (2, 15, 2))))
        self.config[self.method_cluster].pop('n_clusters', None)
        self.score_list = []
        max_score = -np.Inf

        time.sleep(0.2) # To avoid tqdm bar from appearing before prints
        pbar = tqdm.tqdm(self.n_clusters_list)
        for n_clusters in pbar:
            pbar.set_description("Trying n_clusters=" + str(n_clusters))
            method_obj = globals()[self.method_cluster](
                n_clusters=n_clusters, **self.config[self.method_cluster]
            )
            y_pred = method_obj.fit_predict(self.x_emb)
            score = globals()[self.method_cluster_eval](
                self.x_emb, y_pred, **self.config[self.method_cluster_eval]
            )
            self.score_list.append(score)
            if max_score < score: # Update best score
                max_score = score
                self.n_clusters = n_clusters
                self.y_pred = y_pred

        print("\nClustering complete.")
        print("Highest score is achieved for n_clusters =", self.n_clusters)

    def find_markers(self):
        """
        After having constructed the clusters, look at all the genes and try
        to determine which genes are significant for a given cluster. Collect
        all such gene ids into "self.marker_ids" list.
        """
        if self.y_pred is None:
            self.cluster()
        print_header("Finding Markers")
        self.print_verbose('markers')
        # Read params from config
        alpha = self.config[self.markers].get('alpha', 0.05)
        difference = self.config[self.markers].get('difference', 1)
        top_k = self.config[self.markers].get('top_k', 100)
        # Initialize
        self.marker_indices = []
        self.pvals = np.zeros((self.n_clusters, self.colnum))
        self.mds = np.zeros((self.n_clusters, self.colnum))

        time.sleep(0.2) # To avoid tqdm bar from appearing before prints
        for cluster_id in tqdm.tqdm(range(self.n_clusters), desc="Completed clusters:"):
            # Current cluster vs the rest
            x_in = self.x[np.where(self.y_pred == cluster_id)]
            x_not_in = self.x[np.where(self.y_pred != cluster_id)]

            # Consider all genes and run T-test to determine if that gene
            # is significant for the current cluster
            for gene in range(self.colnum):
                gene_vector_in = x_in[:, gene]
                gene_vector_not_in = x_not_in[:, gene]
                # T-test + mean difference (md)
                t, pval = stats.ttest_ind(gene_vector_in, gene_vector_not_in, equal_var=False)
                self.pvals[cluster_id][gene] = pval
                # No absolute value needed for the following difference
                # While large negative values would also imply significance,
                # they are not considered markers for the given cluster
                self.mds[cluster_id][gene] = np.mean(gene_vector_in) - np.mean(gene_vector_not_in)

            # Sort the mean differences and use the indices to sort p-values
            sorted_indices = np.flip(np.argsort(self.mds[cluster_id]))
            pvals_sorted_by_mad = self.pvals[cluster_id][sorted_indices]
            # Apply Bonferroni correction to the p-values (i.e., divide alpha by cols)
            marker_index = sorted_indices[np.where(pvals_sorted_by_mad < (alpha / self.colnum))]
            # Choose top_k genes
            self.marker_indices.append(marker_index[:top_k])
        self.marker_indices = np.asarray(self.marker_indices)

        # Convert found indices to gene IDs
        # Also keep track of p-values and md's
        self.marker_ids = self.col_ids[self.marker_indices]
        self.marker_pvals = np.array([self.pvals[i][self.marker_indices[i]] for i in range(len(self.marker_indices))])
        self.marker_mds = np.array([self.mds[i][self.marker_indices[i]] for i in range(len(self.marker_indices))])

    def convert_markers(self):
        """
        Converts gene IDs to gene names.
        """
        if self.marker_ids is None:
            self.find_markers()
        self.marker_names = gene_id_to_name(self.marker_ids)
        self.pop_names, self.svs, self.intersec, self.sub_pop_names, self.sub_svs, self.sub_intersec = gene_name_to_cell(self.marker_names)

    def get_clusters_csv(self):
        self.filter_genes()
        self.reduce_dim()
        self.cluster()
        emb = self.get_visual_emb()
        result = np.hstack([emb, np.expand_dims(self.y_pred, axis=1)])
        np.savetxt('csv/clusters.csv', result, delimiter=',')

    def new_hard_cluster(self, labels, indices=None, all_points=False):
        """
        Update the label of points whose indices are given in indices,
        or it all_points flag is set, then assume all labels are given.
        """
        if all_points == False:
            assert(indices is not None)
            self.y_pred[indices] = labels
        else:
            assert(len(labels) == len(y_pred))
            self.y_pred = labels

    def new_soft_cluster(self, point_index, k=None):
        """
        Given a single point, find the cluster where that point belongs
        determined by using elbow heuristics and update the labels.
        """
        if self.x_2d_emb is None:
            self.get_visual_emb()
        distances = np.linalg.norm(self.x_2d_emb[point_index] - self.x_2d_emb, axis=1)
        sorted_indices = np.argsort(distances)
        plt.plot(range(0, len(distances)), distances[sorted_indices])
        knee = KneeLocator(range(0, len(distances)), distances[sorted_indices],
                            curve='convex',
                            direction='increasing',
                            S=2)
        print('knee at', knee.knee)
        if k is not None:
            self.y_pred[sorted_indices[:k]] = self.n_clusters
            self.n_clusters += 1
            self.find_markers()
            self.convert_markers()

    # def cluster_constraints(self, method='agglomerative', n_clusters_list=list(range(2, 65, 2)), constraints=None, **kwargs):
    #     assert constraints is not None, "No constraints provided."
    #     print_header("Clustering with constraints")
    #     print("Using " + method + ".")

    #     silhouette_score_list = []
    #     max_sscore = -2

    #     time.sleep(0.2) # To avoid tqdm bar from appearing before prints
    #     pbar = tqdm.tqdm(n_clusters_list)
    #     for n_clusters in pbar:
    #         pbar.set_description("Trying n_clusters=" + str(n_clusters))
    #         if method == 'cop_kmeans':
    #             y_pred, centers = cop_kmeans(dataset=self.x_emb, k=n_clusters,
    #                                                 ml=constraints['ml'], cl=constraints['cl'],
    #                                                 **kwargs)
    #         elif method == 'agglomerative':
    #             ac = AgglomerativeClustering(n_clusters=n_clusters, connectivity=constraints, **kwargs)
    #             y_pred = ac.fit_predict(self.x_emb)
    #         elif method == 'pckmeans':
    #             pck = PCKMeans(n_clusters=self.n_clusters)
    #             pck.fit(self.x_emb, ml=constraints['ml'], cl=constraints['cl'])
    #             y_pred = pck.labels_
    #         else:
    #             raise NotImplementedError()

    #         sscore = silhouette_score(self.x_emb, y_pred, metric='correlation')
    #         silhouette_score_list.append(sscore)
    #         if max_sscore < sscore:
    #             max_sscore = sscore
    #             self.n_clusters = n_clusters
    #             self.y_pred = y_pred

    #     print("Clustering with constraints complete.")
    #     print("Highest silhouette score is achieved for n_clusters =", self.n_clusters)

    #     if self.verbose:
    #         sns.set_style('whitegrid')
    #         fig, ax = plt.subplots()
    #         ax.plot(n_clusters_list, silhouette_score_list)
    #         ax.set_xlabel('Number of Clusters')
    #         ax.set_ylabel('Silhouette Score')
    #         sns.despine()
    #         fig.set_size_inches(10, 5)

    def flow(self):
        self.reduce_dim()
        self.cluster()
        self.find_markers()
        self.convert_markers()

    def plot(self, what):
        if what == "explained_variance":
            plot_explained_variance(list(range(1, len(self.explained_variance) + 1)),
                                    self.explained_variance)
        elif what == "gene_variances":
            plot_gene_variances(self.gene_variances)
        elif what == "cluster_scores":
            plot_scores(self.n_clusters_list, self.score_list)
        elif what == "marker_hist":
            plot_marker_hist(self.n_clusters, self.pvals, self.mds)
        elif what == "top_markers":
            plot_top_markers(self.marker_ids, self.marker_pvals, self.marker_mds)
        elif what == "2d":
            self.get_visual_emb()
            labels = None
            if hasattr(self, 'pop_names'):
                labels = ["{0} ({1} common genes; sv={2:.2f})\n{3} ({4} common genes; sv={5:.2f})".format(
                                                        x,
                                                        len(self.intersec[i]),
                                                        self.svs[i],
                                                        self.sub_pop_names[i],
                                                        len(self.sub_intersec[i]),
                                                        self.sub_svs[i])
                                for i, x in enumerate(self.pop_names)]
            plot_2d(self.x_2d_emb, self.y_pred, labels=labels)
        else:
            raise ValueError()

def pretty_print_dict(mydict):
    print(json.dumps(mydict, sort_keys=True, indent=4))

def print_header(text, max_length=80, prepend_newline=True):
    if len(text) >= max_length:
        if prepend_newline:
            print()
        print("*** " + text + " ***")
    else:
        rems = int((max_length - len(text)) / 2) - 2
        if prepend_newline:
            print()
        print('*' * rems, text, '*' * (max_length - rems - len(text) - 2))