import datetime
import pickle
from ast import literal_eval

import numpy as np
import pandas as pd
import gseapy

from ._wrapper import wrap
from .log import setup_logger
from .units._unit import Unit
from .utils.read import parse_config
from .utils.read import load_data
from .utils.validation import _effective_numerical_value as env


def is_float(val):
    try:
        num = float(val)
    except ValueError:
        return False
    return True


class Pipeline(Unit):
    def __init__(self, x='default', config='configs/config.ini', col_ids=None):
        if x == 'default':
            x, col_ids = load_data('default')
            print("Loaded data.")
        if type(x) == str:
            x, col_ids = load_data(x)
            print("Loaded data.")
        if type(x) != np.ndarray:
            x = np.array(x)
        if len(x.shape) != 2:
            raise ValueError(
                "Data needs to be of shape (n_samples, n_features).")
        assert isinstance(config, str)
        self.x = x
        self.config = parse_config(config)
        self.col_ids = np.array(col_ids).astype('U').reshape(-1)
        self.create_objects()
        self.updated = False

    def create_objects(self, methods=None):
        try:
            clu_method = self.config["methods"]["cluster"]
            eval_method = self.config["methods"]["cluster_eval"]
            vis_method = self.config["methods"]["visualization"]
            mark_method = self.config["methods"]["markers"]
            con_method = self.config["methods"]["conversion"]
            ide_method = self.config["methods"]["identification"]
            ssclu_method = self.config["methods"]["ss_cluster"]
        except:
            raise ValueError("Pipe: Config file malformed or method missing.")

        self.ssclu = wrap("ss_cluster", ssclu_method)(
            **self.config["ss_cluster"]
        )

    def run(self, dim_method='PCA', dim_n_components='knee',
            clu_method='KMeans', eval_method='Silhouette',
            clu_n_clusters='(3, 5, 1)', clu_n_jobs=-1,
            mark_method='TTest', mark_alpha=0.05, mark_markers_n=200,
            mark_correction='hs', mark_n_jobs=1, con_method="Converter",
            con_convention='id-to-name', con_path="markers/gene_id_name.csv",
            ide_method='HyperGeom', ide_path='markers/cell_type_marker.json',
            ide_tissue='all'
            ):
        self.x_emb = self.get_emb(
            self.x, method=dim_method, n_components=dim_n_components)

        self.labels = self.get_labels(
            self.x_emb, method=clu_method,
            eval_method=eval_method, n_clusters=clu_n_clusters,
            n_jobs=clu_n_jobs)

        self.markers = self.get_markers(
            self.x, self.labels, method=mark_method,
            alpha=mark_alpha, markers_n=mark_markers_n,
            correction=mark_correction, n_jobs=mark_n_jobs)

        self.markers = self.convert(
            self.markers, col_ids=self.col_ids, method=con_method,
            convention=con_convention, path=con_path)

        self.markers = self.identify(
            self.markers, method="HyperGeom",
            path=ide_path, tissue=ide_tissue)

    def get_emb(self, x=None, method="PCA", n_components='knee', **kwargs):
        if x is None:
            x = self.x

        if is_float(n_components):
            n_components = literal_eval(n_components)

        self.dim = wrap("dim_reduction", method)(
            n_components=n_components,
            **kwargs
        )
        self.x_emb = self.dim.get(x)
        return self.x_emb

    def get_labels(self, x=None, method="KMeans", eval_method="Silhouette",
                   n_clusters='(3, 5, 1)', n_jobs=-1, **kwargs):
        if x is None:
            x = self.x_emb

        n_clusters = env(n_clusters)
        n_jobs = env(n_jobs)

        self.eval = wrap("cluster_eval", eval_method)()
        self.clu = wrap("cluster", method)(
            eval_obj=self.eval, n_clusters=n_clusters, n_jobs=n_jobs, **kwargs
        )
        self.labels = self.clu.get(x)
        return self.labels

    def get_markers(self, x=None, y=None, method="TTest", alpha=0.05,
                    markers_n=200, correction='hs', n_jobs=1):
        if x is None:
            x = self.x
        if y is None:
            y = self.labels

        alpha = env(alpha)
        markers_n = env(markers_n)
        n_jobs = env(n_jobs)

        self.mark = wrap("markers", method)(
            alpha=alpha, markers_n=markers_n,
            correction=correction, n_jobs=n_jobs
        )
        unq_labels = np.unique(y)
        # 3. Differential expression
        self.markers = self.mark.get(x, y, unq_labels)
        return self.markers

    def convert(self, markers=None, col_ids=None, method="Converter",
                convention='id-to-name', path='markers/gene_id_name.csv'):
        if markers is None:
            markers = self.markers
        if col_ids is None:
            col_ids = self.col_ids

        self.con = wrap("conversion", method)(
            convention=convention, path=path)

        markers = markers.copy()
        for marker in markers:
            markers[marker]['inp_names'] = col_ids[
                markers[marker]['indices']
            ]
            markers[marker]['outp_names'] = self.con.get(
                markers[marker]['inp_names']
            )
        self.markers = markers
        return self.markers

    def identify(self, markers=None, method="HyperGeom",
                 path='markers/cell_type_marker.json', tissue='all'):
        if markers is None:
            markers = self.markers

        self.ide = wrap("identification", method)(
            path=path, tissue=tissue,
        )
        markers = markers.copy()
        self.markers = self.ide.get(markers)
        return self.markers

    def update(self, x=None, new_labels=None, method='SeededKMeans', **kwargs):
        if x is None:
            x = self.x_emb

        self.ssclu = wrap("ss_cluster", method)()

        self.labels = self.ssclu.get(x, new_labels, **kwargs)

        self.get_markers()
        self.convert()
        self.identify()

    def get_emb_2d(self, x=None, y=None, method='UMAP'):
        if x is None:
            x = self.x_emb

        if method == 'UMAP':
            self.vis = wrap("dim_reduction", method)(
                n_components=2,
                n_neighbors=15,
                min_dist=0,
                metric="euclidean"
            )
        else:
            self.vis = wrap("dim_reduction", method)()

        return self.vis.get(x, y)

    def get_markers_subset(self, indices1, indices2=None):
        if indices2 is None:
            markers = self.mark.get_subset(self.x, indices1)
            # Convert
            for marker in markers:  # should be only 1
                markers[marker]['inp_names'] = self.col_ids[
                    markers[marker]['indices']
                ]
                markers[marker]['outp_names'] = self.con.get(
                    markers[marker]['inp_names']
                )
            markers = self.ide.get(markers)
            return markers
        else:
            x1 = self.x[indices1]
            labels1 = np.zeros((x1.shape[0],), dtype=int)
            x2 = self.x[indices2]
            labels2 = np.ones((x2.shape[0],), dtype=int)

            x = np.concatenate([x1, x2])
            labels = np.concatenate([labels1, labels2])

            markers = self.mark.get(x, labels)
            for marker in markers:
                markers[marker]['inp_names'] = self.col_ids[
                    markers[marker]['indices']
                ]
                markers[marker]['outp_names'] = self.con.get(
                    markers[marker]['inp_names']
                )
            markers = self.ide.get(markers)
            return markers

    def enrich(self, indices1, indices2=None):
        if type(indices1) != np.ndarray:
            indices1 = np.array(indices1).astype(int).reshape(-1)
        if indices2 is not None and type(indices2) != np.ndarray:
            indices2 = np.array(indices2).astype(int).reshape(-1)

        markers = self.get_markers_subset(indices1=indices1, indices2=indices2)
        for key in markers:
            gseapy.enrichr(gene_list=markers[key]['outp_names'].tolist(),
                           description=key,
                           format='png',
                           gene_sets='Human_Gene_Atlas')

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

    def get(self):
        return self.x_emb_2d, self.labels, self.markers
