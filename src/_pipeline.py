from ._wrapper import wrap
from .utils.read import parse_config
from .units._unit import Unit

import numpy as np
import pandas as pd
import pickle
import datetime


class Pipeline(Unit):
    def __init__(self, x, config, verbose=False,
                 row_ids=None, col_ids=None):
        super().__init__(verbose)
        assert len(x.shape) == 2, "Pipe: Data needs to be of shape (n x d)."
        assert isinstance(config, str)
        self.x = x
        self.config = parse_config(config)
        self.verbose = verbose
        self.row_ids = row_ids.astype('U') if row_ids is not None else None
        self.col_ids = col_ids.astype('U') if col_ids is not None else None
        self.create_objects()
        self.updated = False
        self.name = 'Pipe'

    def create_objects(self, methods=None):
        try:
            dim_method = self.config["methods"]["dim_reduction"]
            clu_method = self.config["methods"]["cluster"]
            eval_method = self.config["methods"]["cluster_eval"]
            vis_method = self.config["methods"]["visualization"]
            mark_method = self.config["methods"]["markers"]
            con_method = self.config["methods"]["conversion"]
            ide_method = self.config["methods"]["identification"]
            ssclu_method = self.config["methods"]["ss_cluster"]
        except:
            raise ValueError("Pipe: Config file malformed or method missing.")

        self.dim = wrap("dim_reduction", dim_method)(
            self.verbose, **self.config["dim_reduction"]
        )
        self.eval = wrap("cluster_eval", eval_method)(
            self.verbose, **self.config["cluster_eval"]
        )
        self.clu = wrap("cluster", clu_method)(
            eval_obj=self.eval, **self.config["cluster"]
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
        self.ssclu = wrap("ss_cluster", ssclu_method)(
            self.verbose, **self.config["ss_cluster"]
        )

    def run(self):
        self.emb()
        self.cluster()
        self.get_markers()
        self.convert()
        self.identify()

    def emb(self):
        self.x_emb = self.dim.get(self.x)

    def cluster(self):
        self.labels = self.clu.get(self.x_emb)

    def get_markers(self):
        self.unq_labels = np.unique(self.labels)
        # 3. Differential expression
        self.markers = self.mark.get(
            self.x,
            self.labels,
            self.unq_labels
        )

    def get_markers_subset(self, indices):
        markers = self.mark.get_subset(self.x, indices)
        # Convert
        for marker in markers: # should be only 1
            markers[marker]['inp_names'] = self.col_ids[
                markers[marker]['indices']
            ]
            markers[marker]['outp_names'] = self.con.get(
                markers[marker]['inp_names']
            )
        markers = self.ide.get(markers)
        return markers

    def convert(self):
        for marker in self.markers:
            self.markers[marker]['inp_names'] = self.col_ids[
                self.markers[marker]['indices']
            ]
            self.markers[marker]['outp_names'] = self.con.get(
                self.markers[marker]['inp_names']
            )

    def identify(self):
        self.markers = self.ide.get(self.markers)

    def get_emb_2d(self):
        self.x_emb_2d = self.vis.get(self.x_emb, self.labels)
        return self.x_emb_2d

    def save_plot_info(self, path=None):
        """
        Saves x, y coordinates and label for every point.
        """
        df = pd.DataFrame()
        if not hasattr(self, 'x_emb_2d'):
            self.x_emb_2d = self.vis.get(self.x_emb, self.labels)

        df['x'] = self.x_emb_2d[:, 0]
        df['y'] = self.x_emb_2d[:, 1]
        df['label'] = self.labels

        if path is None:
            fn = datetime.datetime.now().strftime("%y%m%d-%H-%M-%S")
            path = "states/plot-info-" + fn + ".csv"
        df.to_csv(path)
        return path

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
            path = "csv/marker-info-" + fn + ".csv"

        df.to_csv(path)
        return path

    def save(self, path=None):
        """
        Saves current object state into a pickle file.
        """
        if path is None:
            fn = datetime.datetime.now().strftime("%y%m%d-%H-%M-%S")
            path = "csv/pipe-" + fn + ".pkl"
        with open(path, "wb") as f:
            pickle.dump(self, f)
        return path

    def update(self, new_labels=None, code=100):
        """
        Given new_labels (1D, same size as labels), update according to code;
        Codes:
            100: Hard cluster. Update the given labels. Don't run clustering.
            200: Soft cluster. Run constrained clustering using new labels.
        """
        # labels assumed to be 1D and same dimension as original labels
        if code == 100:
            self.labels = new_labels
        elif code == 200:
            self.labels = self.ssclu.get(self.x, new_labels)
        else:
            raise ValueError("Invalid code.")

        self.get_markers()
        self.convert()
        self.identify()

    def get(self):
        return self.x_emb_2d, self.labels, self.markers
