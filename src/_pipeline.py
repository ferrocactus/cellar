import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import json
import time
import tqdm

from .wrapper import wrap

# Differential Expression
from scipy import stats
# Utils
from .utils.utils_visualization import (
    plot_2d,
    plot_explained_variance,
    plot_gene_variances,
    plot_scores,
    plot_marker_hist,
    plot_top_markers
)
from .utils.utils_experiment import (
    read_config,
    gene_id_to_name,
    gene_name_to_cell
)
from .utils.utils_read import parse_config
# Constrained clustering
#from copkmeans.cop_kmeans import cop_kmeans
#from active_semi_clustering.semi_supervised.pairwise_constraints import PCKMeans


class Pipeline:
    def __init__(self, x, config='configs/config.ini', verbose=False, row_ids=None, col_ids=None):
        assert len(x.shape) == 2, "Data needs to be of shape (n x d)."
        self.x = x
        self.row_ids = row_ids.astype('U') if row_ids is not None else None
        self.col_ids = col_ids.astype('U') if col_ids is not None else None
        self.config = parse_config(config)
        self.verbose = verbose
        self.create_objects()
        self.x_emb = None # Embeddings after running PCA
        self.labels = None # Labels after running clustering
        self.marker_ids = None # Markers for cluster

    def create_objects(self, methods=None):
        try:
            dim_method = self.config["methods"]["dim_reduction"]
            clu_method = self.config["methods"]["cluster"]
            eval_method = self.config["methods"]["cluster_eval"]
            vis_method = self.config["methods"]["visualization"]
        except:
            raise "Config file malformed or method missing."

        self.dim_obj = wrap("dim_reduction", dim_method)(self.verbose, **self.config[dim_method])
        self.clu_obj = wrap("cluster", clu_method)(self.verbose, **self.config[clu_method])
        self.eval_obj = wrap("cluster_eval", eval_method)(self.verbose, **self.config[eval_method])
        self.vis_obj = wrap("dim_reduction", vis_method)(self.verbose, **self.config[vis_method])

    def run(self):
        self.x_emb = self.dim_obj.get(self.x)
        self.labels = self.clu_obj.get(self.x_emb, self.eval_obj)
        self.n_clusters = self.clu_obj._n_clusters
        self.find_markers()
        self.convert_markers()

    def find_markers(self):
        """
        After having constructed the clusters, look at all the genes and try
        to determine which genes are significant for a given cluster. Collect
        all such gene ids into "self.marker_ids" list.
        """
        print_header("Finding Markers")
        # Read params from config
        alpha = self.config['markers'].get('alpha', 0.05)
        difference = self.config['markers'].get('difference', 1)
        top_k = self.config['markers'].get('top_k', 100)
        # Initialize
        self.marker_indices = []
        self.pvals = np.zeros((self.n_clusters, self.x.shape[1]))
        self.mds = np.zeros((self.n_clusters, self.x.shape[1]))

        time.sleep(0.2) # To avoid tqdm bar from appearing before prints
        for cluster_id in tqdm.tqdm(range(self.n_clusters), desc="Completed clusters:"):
            # Current cluster vs the rest
            x_in = self.x[np.where(self.labels == cluster_id)]
            x_not_in = self.x[np.where(self.labels != cluster_id)]

            # Consider all genes and run T-test to determine if that gene
            # is significant for the current cluster
            for gene in range(self.x.shape[1]):
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
            marker_index = sorted_indices[np.where(pvals_sorted_by_mad < (alpha / self.x.shape[1]))]
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

    def new_hard_cluster(self, labels, indices=None, all_points=False):
        """
        Update the label of points whose indices are given in indices,
        or it all_points flag is set, then assume all labels are given.
        """
        if all_points == False:
            assert(indices is not None)
            self.labels[indices] = labels
        else:
            assert(len(labels) == len(labels))
            self.labels = labels

    def new_soft_cluster(self, point_index, k=None):
        """
        Given a single point, find the cluster where that point belongs
        determined by using elbow heuristics and update the labels.
        """
        distances = np.linalg.norm(self.x_2d_emb[point_index] - self.x_2d_emb, axis=1)
        sorted_indices = np.argsort(distances)
        plt.plot(range(0, len(distances)), distances[sorted_indices])
        knee = KneeLocator(range(0, len(distances)), distances[sorted_indices],
                            curve='convex',
                            direction='increasing',
                            S=2)
        print('knee at', knee.knee)
        if k is not None:
            self.labels[sorted_indices[:k]] = self.n_clusters
            self.n_clusters += 1
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
            plot_top_markers(self.marker_ids[:10], self.marker_pvals[:10], self.marker_mds[:10])
        elif what == "2d":
            emb_2d = self.vis_obj.get(self.x_emb, self.labels)
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
            plot_2d(emb_2d, self.labels, labels=labels)
        else:
            raise ValueError()


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
