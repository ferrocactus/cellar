from ._wrapper import wrap

import numpy as np

# Utils
from .utils.utils_experiment import read_config
from .utils.utils_read import parse_config
# Constrained clustering
#from copkmeans.cop_kmeans import cop_kmeans
#from active_semi_clustering.semi_supervised.pairwise_constraints import PCKMeans


class Pipeline:
    def __init__(self, x, config, verbose=False,
                row_ids=None, col_ids=None):
        assert len(x.shape) == 2, "Pipe: Data needs to be of shape (n x d)."
        assert isinstance(config, str)
        self.x = x
        self.config = parse_config(config)
        self.verbose = verbose
        self.row_ids = row_ids.astype('U') if row_ids is not None else None
        self.col_ids = col_ids.astype('U') if col_ids is not None else None
        self.create_objects()

    def create_objects(self, methods=None):
        try:
            dim_method = self.config["methods"]["dim_reduction"]
            clu_method = self.config["methods"]["cluster"]
            eval_method = self.config["methods"]["cluster_eval"]
            vis_method = self.config["methods"]["visualization"]
            mark_method = self.config["methods"]["markers"]
            con_method = self.config["methods"]["conversion"]
            ide_method = self.config["methods"]["identification"]
        except:
            raise ValueError("Pipe: Config file malformed or method missing.")

        self.dim = wrap("dim_reduction", dim_method)(
            self.verbose, **self.config["dim_reduction"]
        )
        self.clu = wrap("cluster", clu_method)(
            self.verbose, **self.config["cluster"]
        )
        self.eval = wrap("cluster_eval", eval_method)(
            self.verbose, **self.config["cluster_eval"]
        )
        self.vis = wrap("dim_reduction", vis_method)(
            self.verbose, **self.config["visualization"]
        )
        self.mark = wrap("markers", mark_method)(
            self.verbose, **self.config["markers"]
        )
        self.con = wrap("conversion", con_method)(
            self.verbose, **self.config["conversion"]
        )
        self.ide = wrap("identification", ide_method)(
            self.verbose, **self.config["identification"]
        )

    def run(self):
        self.x_emb = self.dim.get(self.x)
        self.labels = self.clu.get(self.x_emb, self.eval)
        self.unq_labels = np.unique(self.labels)
        self.markers = self.mark.get(self.x, self.labels, self.unq_labels)
        for marker in self.markers:
            self.markers[marker]['ids'] = self.col_ids[
                self.markers[marker]['indices']
            ]
            self.markers[marker]['names'] = self.con.get(
                self.markers[marker]['ids']
            )
        self.markers = self.ide.get(self.markers)

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