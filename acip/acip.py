import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import json
import tqdm
import time
import sys

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
from utils.utils_visualization import (plot_2d,
                                    plot_explained_variance,
                                    plot_gene_variances,
                                    plot_scores,
                                    plot_marker_hist,
                                    plot_top_markers)
from utils.utils_experiment import read_config
# Constrained clustering
from copkmeans.cop_kmeans import cop_kmeans
from active_semi_clustering.semi_supervised.pairwise_constraints import PCKMeans

class ACIP:
    def __init__(self, x, config="defaults", verbose=False, row_ids=None, col_ids=None):
        assert type(x).__module__ == np.__name__
        assert type(row_ids).__module__ == np.__name__ or row_ids is None
        assert type(col_ids).__module__ == np.__name__ or col_ids is None
        self.set_train_data(x)
        # If use_config is set for a function, all its kwargs will be ignored
        self.config = config
        self.read_config()
        self.verbose = verbose
        self.set_ids(row_ids, col_ids)
    
    def set_train_data(self, x):
        assert len(x.shape) == 2, "Data needs to be of shape (n x d)."
        self.x = x
        self.rows = x.shape[0]
        self.cols = x.shape[1]
    
    def set_ids(self, row_ids=None, col_ids=None):
        self.row_ids = row_ids # Cell IDs
        self.col_ids = col_ids # Gene IDs
    
    def read_config(self):
        self.params = read_config(self.config)
    
    def normalize(self):
        raise NotImplementedError()

    def print_verbose(self, step):
        # Prints configs for given processing step. Reads from self.params
        if self.verbose == False:
            return
        
        if "method" in self.params[step]:
            method = self.params[step]['method']
            print("Using", method, "with params:")
            pretty_print_dict(self.params[step][method])
        else:
            print("Using params:")
            pretty_print_dict(self.params[step])
    
    def get_visual_emb(self):
        """
        Constructs 2d/3d embeddings for visualization.
        If use_emb=1, then use embeddings x_emb otherwise use x.
        Method should be one of UMAP, PCA, TSNE.
        """
        use_x_emb = self.params['visual_emb']['use_x_emb']
        if use_x_emb and not hasattr(self, 'x_emb'):
            sys.exit("No embeddings found. Consider setting use_x_emb=0.")
        print_header("Constructing visual embeddings.")
        self.print_verbose('visual_emb')

        x = self.x_emb if use_x_emb else self.x

        # Use the method specified in the config file
        method = self.params['visual_emb']['method']
        assert method in ["UMAP", "PCA", "TSNE"]
        method_obj = globals()[method](**self.params['visual_emb'][method])
        return method_obj.fit_transform(x)

    def get_explained_variance(self):
        """
        Returns the percentage of variance explained by each of the PCA components.
        """
        print_header("Getting explained variance using PCA.")
        self.print_verbose('explained_variance')

        pca = PCA(**self.params['explained_variance'])
        pca.fit(self.x)
        return list(range(1, pca.n_components_ + 1)), pca.explained_variance_ratio_

    def filter_genes(self):
        """
        Computes the variance for every gene (column). If 0 < n_genes <= 1, then
        use the top (n_genes*cols) genes, otherwise use the top (n_genes). It is
        possible to set n_genes="all" in which case no filtering is done.
        """
        print_header("Filtering genes")
        self.print_verbose("filter_genes")

        n_genes = self.params['filter_genes']['n_genes']
        assert n_genes > 0

        self.gene_variances = np.var(self.x, axis=0)
        if n_genes <= 1: # Fraction of genes to keep
            highest_var_gene_indices = self.gene_variances.argsort()[-int(n_genes * self.cols):]
        else:            # Count of genes to keep
            highest_var_gene_indices = self.gene_variances.argsort()[-n_genes:]
        
        print("Keeping", len(highest_var_gene_indices), "genes.")

        self.x_filtered = self.x[:, highest_var_gene_indices]

    def reduce_dim(self):
        """
        Reduces the dimensionality of the data and stores it in self.x_emb.
        Requires filter_genes to be run first.
        """
        if not hasattr(self, 'x_filtered'):
            sys.exit("Genes not filtered. Run filter_genes first.")
        print_header("Reducing dimensionality")
        self.print_verbose('reduce_dim')

        method = self.params['reduce_dim']['method']
        assert method in ["PCA"]
        method_obj = globals()[method](**self.params['reduce_dim'][method])
        
        if method == 'PCA':
            method_obj.fit(self.x_filtered)
            self.x_emb = method_obj.transform(self.x_filtered)
            if self.verbose:
                print("Embedding created. MSE:",
                        mse(self.x_filtered, method_obj.inverse_transform(self.x_emb)))
            print("Used", method_obj.n_components_, "components.")
        else:
            raise NotImplementedError()
    
    def cluster(self):
        """
        Clusters the embeddings created by reduce_dim. Tries various n_clusters
        as specified by the list (low, high, step) and stores into y_pred the
        highest scoring labels according to eval_method.
        """
        if not hasattr(self, 'x_emb'):
            sys.exit("Embeddings not created. Run reduce_dim first.")
        print_header("Clustering")
        self.print_verbose('cluster')

        method = self.params['cluster']['method']
        assert method in ["KMeans", "KMedoids", "SpectralClustering"]
        eval_method = self.params['cluster']['eval_method']
        assert eval_method in ["silhouette_score"]

        # Construct list of n_clusters to try
        low = self.params['cluster']['low']
        high = self.params['cluster']['high']
        step = self.params['cluster']['step']
        self.n_clusters_list = list(range(low, high, step))
        self.score_list = []
        max_score = -np.Inf
        
        time.sleep(0.2) # To avoid tqdm bar from appearing before prints
        pbar = tqdm.tqdm(self.n_clusters_list)
        for n_clusters in pbar:
            pbar.set_description("Trying n_clusters=" + str(n_clusters))
            method_obj = globals()[method](n_clusters=n_clusters, **self.params['cluster'][method])
            y_pred = method_obj.fit_predict(self.x_emb)
            score = globals()[eval_method](self.x_emb, y_pred, **self.params['cluster'][eval_method])
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
        all such gene indices into "self.markers" list.
        """
        if not hasattr(self, 'y_pred'):
            sys.exit("Clustering not performed. Run cluster first.")
        print_header("Finding Markers")
        self.print_verbose('markers')
        # Read params from config
        alpha = self.params['markers']['alpha']
        difference = self.params['markers']['difference']
        top_k = self.params['markers']['top_k']
        # Initialize
        self.marker_indices = []
        self.pvals = np.zeros((self.n_clusters, self.cols))
        self.mds = np.zeros((self.n_clusters, self.cols))

        time.sleep(0.2) # To avoid tqdm bar from appearing before prints
        for cluster_id in tqdm.tqdm(range(self.n_clusters), desc="Completed clusters:"):
            # Current cluster vs the rest
            x_in = self.x[np.where(self.y_pred == cluster_id)]
            x_not_in = self.x[np.where(self.y_pred != cluster_id)]

            # Consider all genes and run T-test to determine if that gene
            # is significant for the current cluster
            for gene in range(self.cols):
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
            marker_index = sorted_indices[np.where(pvals_sorted_by_mad < (alpha / self.cols))]
            # Choose top_k genes
            self.marker_indices.append(marker_index[:top_k])
        self.marker_indices = np.asarray(self.marker_indices)
    
    def get_markers(self):
        """
        Returns two arrays where the first array consists of tuples (p-value, md)
        and the second array consists of the corresponding gene IDs
        """
        if not hasattr(self, 'marker_indices'):
            sys.exit("Markers not found. Run find_markers first.")
        
        marker_ids = []
        marker_pvals = []
        marker_mds = []

        for i, marker_index in enumerate(self.marker_indices):
            marker_ids.append(self.col_ids[marker_index]) # Get gene names
            marker_pvals.append(np.array([self.pvals[i][m] for m in marker_index]))
            marker_mds.append(np.array([self.mds[i][m] for m in marker_index]))

        return np.array(marker_ids), np.array(marker_pvals), np.array(marker_mds)

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
                y_pred, centers = cop_kmeans(dataset=self.x_emb, k=n_clusters,
                                                    ml=constraints['ml'], cl=constraints['cl'],
                                                    **kwargs)
            elif method == 'agglomerative':
                ac = AgglomerativeClustering(n_clusters=n_clusters, connectivity=constraints, **kwargs)
                y_pred = ac.fit_predict(self.x_emb)
            elif method == 'pckmeans':
                pck = PCKMeans(n_clusters=self.n_clusters)
                pck.fit(self.x_emb, ml=constraints['ml'], cl=constraints['cl'])
                y_pred = pck.labels_
            else:
                raise NotImplementedError()
            
            sscore = silhouette_score(self.x_emb, y_pred, metric='correlation')
            silhouette_score_list.append(sscore)
            if max_sscore < sscore:
                max_sscore = sscore
                self.n_clusters = n_clusters
                self.y_pred = y_pred
        
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
        
    def flow(self, reduce_dim='pca', cluster='kmedoids', visualize='umap'):
        self.filter_genes()
        self.reduce_dim()
        self.cluster()
        self.find_markers()

    def plot(self, what):
        if what == "explained_variance":
            plot_explained_variance(*self.get_explained_variance())
        elif what == "gene_variances":
            if not hasattr(self, 'gene_variances'):
                sys.exit("No gene variances found. Run filter_genes.")
            plot_gene_variances(self.gene_variances)
        elif what == "cluster_scores":
            if not hasattr(self, 'score_list'):
                sys.exit("No scores found. Run cluster.")
            plot_scores(self.n_clusters_list, self.score_list)
        elif what == "marker_hist":
            if not hasattr(self, 'marker_indices'):
                sys.exit("No markers found. Run find_markers.")
            plot_marker_hist(self.n_clusters, self.pvals, self.mds)
        elif what == "top_markers":
            if not hasattr(self, 'marker_indices'):
                sys.exit("No markers found. Run find_markers.")
            plot_top_markers(*self.get_markers())
        elif what == "2d":
            emb = self.get_visual_emb()
            plot_2d(emb, self.y_pred)
        else:
            raise ValueError()

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