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
        assert isinstance(config, str)
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
            mark_method = self.config["methods"]["markers"]
        except:
            raise "Config file malformed or method missing."

        self.dim_obj = wrap("dim_reduction", dim_method)(self.verbose, **self.config["dim_reduction"])
        self.clu_obj = wrap("cluster", clu_method)(self.verbose, **self.config["cluster"])
        self.eval_obj = wrap("cluster_eval", eval_method)(self.verbose, **self.config["cluster_eval"])
        self.vis_obj = wrap("dim_reduction", vis_method)(self.verbose, **self.config["visualization"])
        self.mark_obj = wrap("markers", mark_method)(self.verbose, **self.config["markers"])

    def run(self):
        self.x_emb = self.dim_obj.get(self.x)
        self.labels = self.clu_obj.get(self.x_emb, self.eval_obj)
        self.n_clusters = self.clu_obj._n_clusters
        self.mark_results = self.mark_obj.get(self.x, self.labels)
        self.marker_indices = []
        for i in self.mark_results.keys():
            self.marker_indices.append(self.mark_results[i]['indices'])
        self.marker_indices = np.array(self.marker_indices)
        self.marker_ids = self.col_ids[self.marker_indices]
        self.convert_markers()

    def convert_markers(self):
        """
        Converts gene IDs to gene names.
        """
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
