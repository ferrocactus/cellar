from ._wrapper import wrap
from .utils.utils_read import parse_config

import numpy as np
import pandas as pd
import pickle
import datetime


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
        self.updated = False

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
        # 1. Create embedding
        self.x_emb = self.dim.get(self.x)
        # 2. Cluster
        self.labels = self.clu.get(self.x_emb, self.eval)
        self.unq_labels = np.unique(self.labels)
        # 3. Differential expression
        self.markers = self.mark.get(self.x, self.labels, self.unq_labels)
        for marker in self.markers:
            self.markers[marker]['ids'] = self.col_ids[
                self.markers[marker]['indices']
            ]
            self.markers[marker]['names'] = self.con.get(
                self.markers[marker]['ids']
            )
        # 4. Perform identification
        self.markers = self.ide.get(self.markers)

    def save_plot_info(self, path=None):
        """
        Saves x, y coordinates and label for every point.
        """
        df = pd.DataFrame()
        if not hasattr(self, 'x_emb_2d'):
            self.x_emb_2d = self.vis.get(self.x_emb)

        df['x'] = self.x_emb_2d[:, 0]
        df['y'] = self.x_emb_2d[:, 1]
        df['label'] = self.labels

        if path is None:
            fn = datetime.datetime.now().strftime("%y%m%d-%H-%M-%S")
            path = "states/plot-info-" + fn + ".csv"
        df.to_csv(path)

    def save_marker_info(self, path=None):
        """
        Saves all the information obtained (indices, p-values, differences
        significant gene names, lvl1/2 cell type, lvl1/2 survival values,
        lvl1/2 common gene names, lvl1/2 total gene names for type)
        for every label.
        """
        df = pd.DataFrame.from_dict(self.markers, orient='index')

        if path is None:
            fn = datetime.datetime.now().strftime("%y%m%d-%H-%M-%S")
            path = "states/marker-info-" + fn + ".csv"
        df.to_csv(path)

    def save(self, path=None):
        """
        Saves current object state into a pickle file.
        """
        if path is None:
            fn = datetime.datetime.now().strftime("%y%m%d-%H-%M-%S")
            path = "states/pipe-" + fn + ".pkl"
        with open(path, "wb") as f:
            pickle.dump(self, f)
        return path