import datetime
import pickle
from ast import literal_eval
import traceback
import sys

import numpy as np
import pandas as pd
import gseapy
import anndata
import scanpy

from .units import wrap
from .units import _method_exists
from .log import setup_logger
from .units._unit import Unit
from .utils.read import parse_config
from .utils.read import load_data
from .utils.validation import _validate_dim_n_components
from .utils.validation import _validate_clu_n_clusters
from .utils.validation import _validate_n_jobs
from .utils.validation import _validate_mark_alpha
from .utils.validation import _validate_mark_markers_n
from .utils.validation import _validate_mark_correction
from .utils.validation import _validate_con_convention
from .utils.validation import _validate_ensemble_methods


class Pipeline(Unit):
    def __init__(self, x='default', dataset_source=None,
                    config='configs/config.ini', col_ids=None):
        if isinstance(x, str):
            self.dataset = x
            if x == 'default':
                x, col_ids = load_data('default')
                print("Loaded data.")
            else:
                x, col_ids = load_data(x, dataset_source=dataset_source)
                print("Loaded data.")
        if not isinstance(x, np.ndarray):
            try:
                x = np.array(x)
            except:
                raise ValueError("Incorrect format for x.")
        if len(x.shape) != 2:
            raise ValueError(
                "Data needs to be of shape (n_samples, n_features).")
        assert isinstance(config, str)
        self.x = x
        self.config = parse_config(config)
        self.col_ids = np.array(col_ids).astype('U').reshape(-1)

    def validate_params(self, dim_method='PCA', dim_n_components='knee',
            clu_method='KMeans', eval_method='Silhouette',
            clu_n_clusters='(3, 5, 1)', clu_n_jobs=1,
            mark_method='TTest', mark_alpha=0.05, mark_markers_n=200,
            mark_correction='hs', mark_n_jobs=1, con_method="Converter",
            con_convention='id-to-name', con_path="markers/gene_id_name.csv",
            ide_method='HyperGeom', ide_path='markers/cell_type_marker.json',
            ide_tissue='all', vis_method='UMAP', ssc_method=None,
            ensemble_methods=None):

        try:
            _method_exists('dim_reduction', dim_method)
            _validate_dim_n_components(dim_n_components, *self.x.shape)
            _method_exists('cluster', clu_method)
            _method_exists('cluster_eval', eval_method)
            _validate_clu_n_clusters(clu_n_clusters, self.x.shape[0])
            _validate_n_jobs(clu_n_jobs)
            _method_exists('markers', mark_method)
            _validate_mark_alpha(mark_alpha)
            _validate_mark_markers_n(mark_markers_n, self.x.shape[1])
            _validate_mark_correction(mark_correction)
            _validate_n_jobs(mark_n_jobs)
            _method_exists('conversion', con_method)
            _validate_con_convention(con_convention)
            _method_exists('identification', ide_method)
            _method_exists('visualization', vis_method)
            if clu_method == 'Ensemble':
                _validate_ensemble_methods(ensemble_methods)
            if ssc_method is not None:
                _method_exists('ss_cluster', ssc_method)
        except NotImplementedError as nie:
            print(str(nie))
            return str(nie)
        except ValueError as ve:
            print(str(ve))
            return str(ve)
        except Exception as e:
            print(str(e))
            traceback.print_exc(file=sys.stdout)
            return "A problem occurred."

        return "good"


    def run(self, dim_method='PCA', dim_n_components='knee',
            clu_method='KMeans', eval_method='Silhouette',
            clu_n_clusters='(3, 5, 1)', clu_n_jobs=1,
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

        n_components = _validate_dim_n_components(n_components, *self.x.shape)
        if hasattr(self, 'x_emb') and self.n_components == n_components:
            return self.x_emb

        self.n_components = n_components

        self.dim = wrap("dim_reduction", method)(
            n_components=n_components,
            **kwargs
        )
        self.x_emb = self.dim.get(x)
        return self.x_emb

    def get_labels(self, x=None, method="KMeans", eval_method="Silhouette",
                   n_clusters='(3, 5, 1)', n_jobs=1, **kwargs):
        if x is None:
            x = self.x_emb

        n_clusters = _validate_clu_n_clusters(n_clusters, self.x.shape[0])
        n_jobs = _validate_n_jobs(n_jobs)

        self.clu_method = method
        self.eval = wrap("cluster_eval", eval_method)()
        self.clu = wrap("cluster", method)(
            eval_obj=self.eval, n_clusters=n_clusters, n_jobs=n_jobs, **kwargs
        )
        self.labels = self.clu.get(x)
        self.n_clusters = np.unique(self.labels)
        return self.labels

    def get_markers(self, x=None, y=None, method="TTest", alpha=0.05,
                    markers_n=200, correction='holm-sidak', n_jobs=1):
        if x is None:
            x = self.x
        if y is None:
            y = self.labels

        alpha = _validate_mark_alpha(alpha)
        markers_n = _validate_mark_markers_n(markers_n, self.x.shape[1])
        correction = _validate_mark_correction(correction)
        n_jobs = _validate_n_jobs(n_jobs)

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

        convention = _validate_con_convention(convention)

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
        elif type(x) != np.ndarray:
            x = np.asarray(x)

        if type(new_labels) != np.ndarray:
            new_labels = np.asarray(new_labels)

        self.ssclu = wrap("ss_cluster", method)()

        self.labels = self.ssclu.get(x, new_labels, **kwargs)
        return self.labels

    def get_emb_2d(self, x=None, y=None, method='UMAP'):
        if hasattr(self, 'x_emb_2d'):
            return self.x_emb_2d

        if x is None:
            x = self.x_emb

        if self.clu_method == "Scanpy":
            ann = anndata.AnnData(X=self.x_emb)
            scanpy.pp.neighbors(ann, n_neighbors=10, n_pcs=40)
            scanpy.tl.umap(ann)
            self.x_emb_2d = ann.obsm['X_umap']
            return self.x_emb_2d

        if method == 'UMAP':
            self.vis = wrap("visualization", method)(
                n_components=2,
                n_neighbors=15,
                min_dist=0,
                metric="euclidean"
            )
        else:
            self.vis = wrap("visualization", method)()

        self.x_emb_2d = self.vis.get(x)
        return self.x_emb_2d

    def get_markers_subset(self, indices1, indices2=None, alpha=0.05,
                           correction='holm-sidak', markers_n=200, n_jobs=1,
                           method='TTest'):

        alpha = _validate_mark_alpha(alpha)
        markers_n = _validate_mark_markers_n(markers_n, self.x.shape[1])
        correction = _validate_mark_correction(correction)
        n_jobs = _validate_n_jobs(n_jobs)

        mark = wrap("markers", method)(
            alpha=alpha, markers_n=markers_n,
            correction=correction, n_jobs=n_jobs
        )

        indices1 = np.array(indices1).astype(np.int)
        if indices2 is not None:
            indices2 = np.array(indices2).astype(np.int)

        if indices2 is None:
            markers = mark.get_subset(self.x, indices1)
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

            markers = mark.get(x, labels)
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

    def get_details(self):
        details = {}
        details['w'] = self.x.shape[0]
        details['h'] = self.x.shape[1]
        details['min'] = np.min(self.x)
        details['max'] = np.max(self.x)
        details['range'] = details['max'] - details['min']
        details['means'] = np.mean(self.x, axis=0)
        details['max_mean'] = np.max(details['means'])
        details['min_mean'] = np.min(details['means'])
        details['variances'] = np.var(self.x, axis=0)
        details['max_var'] = np.max(details['variances'])
        details['min_var'] = np.min(details['variances'])

        details_txt =   f"Number of cells: {details['w']}\n" + \
                        f"Number of genes: {details['h']}\n" + \
                        f"Minimum value: {details['min']}\n" + \
                        f"Maximum value: {details['max']}\n" + \
                        f"Range: {details['range']}\n" + \
                        f"Maximum gene mean: {details['max_mean']}\n" + \
                        f"Minimum gene mean: {details['min_mean']}\n" + \
                        f"Maximum gene variance: {details['max_var']}\n" + \
                        f"Minimum gene variance: {details['min_var']}\n"
        return details_txt


    def get(self):
        return self.x_emb_2d, self.labels, self.markers
