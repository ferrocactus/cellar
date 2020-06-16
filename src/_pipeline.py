import traceback
import sys

import numpy as np
import scanpy
import anndata

from .units import wrap
from .units import _method_exists
from .log import setup_logger
from .units._unit import Unit
from .utils.read import load_data
from .utils.validation import _validate_clu_n_clusters
from .utils.validation import _validate_con_convention
from .utils.validation import _validate_dim_n_components
from .utils.validation import _validate_ensemble_methods
from .utils.validation import _validate_mark_alpha
from .utils.validation import _validate_mark_correction
from .utils.validation import _validate_mark_markers_n
from .utils.validation import _validate_n_clusters
from .utils.validation import _validate_n_jobs
from .utils.validation import _validate_subset
from .utils.validation import _validate_new_labels
from .utils.validation import _validate_cluster_list
from .utils.validation import _categorify
from .utils.exceptions import InappropriateArgument
from .utils.exceptions import InvalidArgument
from .utils.exceptions import MethodNotImplementedError

OK = 'good'


class Pipeline():
    def __init__(self, x='default', col_ids=None, row_ids=None):
        self.loadmsg = self.load(x, col_ids, row_ids)

    def restate(self, x='default', col_ids=None, row_ids=None):
        keys = list(self.__dict__.keys())
        for var in keys:
            delattr(self, var)
        self.loadmsg = self.load(x, col_ids, row_ids)

    def load(self, x, col_ids, row_ids):
        if isinstance(x, str):
            try:
                self.dataset = x
                df = load_data(x)
                self.x = df['x']
                self.col_ids = df['col_ids']
                self.row_ids = df['row_ids']
                self.x_emb_precomp = df['precomputed_pca']
                self.col_ids = np.array(self.col_ids).reshape(-1)
                if self.row_ids is not None:
                    self.row_ids = np.array(self.row_ids).reshape(-1)
            except Exception as e:
                traceback.print_exc(file=sys.stdout)
                return "Incorrect data format."
        else:
            try:
                self.dataset = "Noname"
                self.x = np.array(x)
                self.col_ids = np.array(col_ids).astype('U').reshape(-1)
                self.row_ids = np.array(row_ids).astype('U').reshape(-1)
            except Exception as e:
                traceback.print_exc(file=sys.stdout)
                return "Incorrect data format."

        self.col_ids = np.char.split(self.col_ids.astype('U').flatten(),
                                     sep='.', maxsplit=1)
        self.col_ids = np.array([i[0] for i in self.col_ids])
        self.col_names = wrap("conversion", "Converter")().get(self.col_ids)

        if len(self.x.shape) != 2:
            return "Data needs to be of shape (n_samples, n_features)."
        print("Loaded data.")
        return OK

    def validate_all(self, dim_method='PCA', dim_n_components='knee',
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
            _validate_dim_n_components(
                dim_n_components, dim_method, *self.x.shape)
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

        return OK

    def run_all(self, dim_method='PCA', dim_n_components='knee',
                clu_method='KMeans', eval_method='Silhouette',
                clu_n_clusters='(3, 5, 1)', clu_n_jobs=1,
                de_method='TTest', de_alpha=0.05, de_n_genes=200,
                de_correction='holm-sidak', de_n_jobs=1, con_method="Converter",
                con_convention='id-to-name', con_path="markers/gene_id_name.csv",
                ide_method='HyperGeom', ide_path='markers/cell_type_marker.json',
                ide_tissue='all'):
        try:
            self.run_step('dim', dim_method=dim_method,
                          dim_n_components=dim_n_components)
            self.run_step('clu', clu_method=clu_method, eval_method=eval_method,
                          clu_n_clusters=clu_n_clusters, clu_n_jobs=clu_n_jobs)
            self.run_step('de', de_method=de_method, de_alpha=de_alpha,
                          de_n_genes=de_n_genes, de_correction=de_correction,
                          de_n_jobs=de_n_jobs)
            self.run_step('con', con_method=con_method,
                          con_convention=con_convention, con_path=con_path)
            self.run_step('ide', ide_method=ide_method,
                          ide_path=ide_path, ide_tissue=ide_tissue)
        except MethodNotImplementedError as e:
            print(str(e))
            return str(e)
        except InappropriateArgument as e:
            print(str(e))
            return str(e)
        except InvalidArgument as e:
            print(str(e))
            return str(e)
        except Exception as e:
            print(str(e))
            traceback.print_exc(file=sys.stdout)
            return "An error occurred."

    def run_step(self, step, **kwargs):
        function_name = "run_" + step
        if hasattr(self, function_name):
            try:
                getattr(self, function_name)(**kwargs)
                return OK
            except MethodNotImplementedError as e:
                print(str(e))
                return str(e)
            except InappropriateArgument as e:
                print(str(e))
                return str(e)
            except InvalidArgument as e:
                print(str(e))
                return str(e)
            except Exception as e:
                print(str(e))
                traceback.print_exc(file=sys.stdout)
                return "An error occurred."
        else:
            raise MethodNotImplementedError("Step not found.")

    def run_dim(self, dim_method="PCA", dim_n_components='knee', **kwargs):
        dim_n_components = _validate_dim_n_components(
            dim_n_components, dim_method, *self.x.shape)

        if dim_method == 'Precomputed PCA':
            if self.x_emb_precomp is not None:
                if dim_n_components > 60:
                    raise InappropriateArgument(
                        "A max of 60 components is precomputed.")
                self.x_emb = self.x_emb_precomp[:, :dim_n_components]
                return
            else:
                raise InappropriateArgument("No precomputed PCA found.")

        _method_exists('dim_reduction', dim_method)

        if hasattr(self, 'dim_method') and hasattr(self, 'dim_n_components_inp')\
                and hasattr(self, 'x_emb'):
            if self.dim_method == dim_method\
                    and self.dim_n_components_inp == dim_n_components:
                return

        self.dim_method = dim_method
        self.dim_n_components_inp = dim_n_components
        self.x_emb = wrap('dim_reduction', dim_method)(
            n_components=dim_n_components, **kwargs).get(self.x)
        self.n_components = self.x_emb.shape[1]

    def run_clu(self, clu_method="KMeans", eval_method="Silhouette",
                clu_n_clusters='(3, 6, 1)', clu_n_jobs=1, **kwargs):
        _method_exists('cluster', clu_method)
        _method_exists('cluster_eval', eval_method)
        clu_n_clusters = _validate_clu_n_clusters(
            clu_n_clusters, self.x.shape[0])
        clu_n_jobs = _validate_n_jobs(clu_n_jobs)
        if clu_method != 'Ensemble':
            kwargs.pop('ensemble_methods', None)
        if "ensemble_methods" in kwargs:
            kwargs['ensemble_methods'] = _validate_ensemble_methods(
                kwargs['ensemble_methods'])

        self.clu_method = clu_method
        self.eval_method = eval_method
        self.clu_n_clusters_inp = clu_n_clusters
        self.labels, scores = wrap("cluster", clu_method)(
            eval_obj=wrap("cluster_eval", eval_method)(),
            n_clusters=clu_n_clusters,
            n_jobs=clu_n_jobs, **kwargs).get(self.x_emb)
        self.n_clusters = np.unique(self.labels)
        self.key_maps = {str(i): i for i in np.unique(self.labels)}

    def set_labels(self, labels):
        self.labels = labels
        self.n_clusters = np.unique(self.labels)
        self.key_maps = {str(i): i for i in np.unique(self.labels)}

    def get_cluster_id(self, name):
        return self.key_maps[name]

    def get_cluster_name(self, id):
        all_values = list(self.key_maps.values())
        all_keys = list(self.key_maps.keys())
        name = all_keys[all_values.index(id)]
        return name

    def run_merge_clusters(self, clusters):
        clusters = _validate_cluster_list(self.labels, clusters)
        if len(clusters) < 2:
            raise InvalidArgument("Not enough clusters found to merge")

        for cluster in clusters[1:]:
            name = self.get_cluster_name(cluster)
            self.key_maps.pop(name, None)
            self.labels[self.labels == cluster] = clusters[0]

    def update_labels(self, name, indices):
        indices = np.array(indices).astype(np.int)
        name = str(name)

        if name in self.key_maps:
            cid = self.get_cluster_id(name)
            self.labels[indices] = cid
        else:
            unq = np.unique(self.labels)
            for x in range(1000):
                if x not in unq:
                    self.key_maps[name] = x
                    self.labels[indices] = x
                    newunq = np.unique(self.labels)
                    oldlabels = list(self.key_maps.values())
                    for oldlabel in oldlabels:
                        if oldlabel not in newunq:
                            oldname = self.get_cluster_name(oldlabel)
                            self.key_maps.pop(oldname, None)
                    break

        self.n_clusters = np.unique(self.labels)

    def run_de(self, de_method="TTest", de_alpha=0.05, de_n_genes=200,
               de_correction='holm-sidak', de_n_jobs=1):
        _method_exists('de', de_method)
        de_alpha = _validate_mark_alpha(de_alpha)
        de_n_genes = _validate_mark_markers_n(de_n_genes, self.x.shape[1])
        de_correction = _validate_mark_correction(de_correction)
        de_n_jobs = _validate_n_jobs(de_n_jobs)

        self.de_method = de_method
        self.de_alpha = de_alpha
        self.de_n_genes = de_n_genes
        self.de_correction = de_correction
        self.markers = wrap("de", de_method)(
            alpha=de_alpha, markers_n=de_n_genes,
            correction=de_correction, n_jobs=de_n_jobs
        ).get(self.x, self.labels, self.n_clusters)

    def run_de_subset(self, subset1, subset2=None, de_method='TTest',
                      de_alpha=0.05, de_n_genes=200,
                      de_correction='holm-sidak'):
        _method_exists('de', de_method)
        subset1 = _validate_subset(subset1)
        subset2 = _validate_subset(subset2)
        de_alpha = _validate_mark_alpha(de_alpha)
        de_n_genes = _validate_mark_markers_n(de_n_genes, self.x.shape[1])
        de_correction = _validate_mark_correction(de_correction)

        self.de_method = de_method
        self.de_alpha = de_alpha
        self.de_n_genes = de_n_genes
        self.de_correction = de_correction

        if subset2 is None:
            self.markers = wrap("de", de_method)(
                alpha=de_alpha, markers_n=de_n_genes,
                correction=de_correction).get_subset(self.x, subset1)
        else:
            # Artificially create dataset and labels
            x1 = self.x[subset1]
            x2 = self.x[subset2]
            x = np.concatenate([x1, x2])

            labels1 = np.zeros((x1.shape[0],), dtype=int)
            labels2 = np.ones((x2.shape[0],), dtype=int)
            labels = np.concatenate([labels1, labels2])

            self.markers = wrap("de", de_method)(
                alpha=de_alpha, markers_n=de_n_genes,
                correction=de_correction).get(x, labels)

    def run_con(self, con_method="Converter", con_convention='id-to-name',
                con_path='markers/gene_id_name.csv'):
        _method_exists('conversion', con_method)
        con_convention = _validate_con_convention(con_convention)

        self.con_method = con_method
        self.con_convention = con_convention
        con = wrap("conversion", con_method)(
            convention=con_convention, path=con_path)
        for marker in self.markers:
            self.markers[marker]['inp_names'] = self.col_ids[
                self.markers[marker]['indices']
            ]
            self.markers[marker]['outp_names'] = con.get(
                self.markers[marker]['inp_names']
            )

    def run_ide(self, ide_method="HyperGeom",
                ide_path='markers/cell_type_marker.json', ide_tissue='all'):
        _method_exists('identification', ide_method)

        self.ide_method = ide_method
        self.markers = wrap("identification", ide_method)(
            path=ide_path, tissue=ide_tissue,
        ).get(self.markers.copy())

    def run_ssclu(self, ssclu_method='SeededKMeans', **kwargs):
        _method_exists('ss_cluster', ssclu_method)
        if 'saved_clusters' in kwargs:
            kwargs['saved_clusters'] = _validate_cluster_list(self.labels,
                                                              kwargs['saved_clusters'])
            saved = {}
            for saved_cluster in kwargs['saved_clusters']:
                name = self.get_cluster_name(saved_cluster)
                if not name.isnumeric():
                    saved[saved_cluster] = (
                        name, np.where(self.labels == saved_cluster))

            self.labels, kwargs['saved_clusters'], repl = _categorify(
                self.labels, kwargs['saved_clusters']
            )

        self.ssclu_method = ssclu_method
        self.labels = wrap(
            "ss_cluster",
            ssclu_method)().get(self.x_emb, self.labels, **kwargs)
        self.n_clusters = np.unique(self.labels)

        self.key_maps = {str(i): i for i in np.unique(self.labels)}

        if 'saved_clusters' in kwargs:
            for i in saved:
                new_lbl = self.labels[saved[i][1]][0]
                old_name = self.get_cluster_name(new_lbl)
                if old_name in self.key_maps:
                    self.key_maps.pop(old_name, None)
                self.key_maps[saved[i][0]] = new_lbl

    def run_vis(self, vis_method='UMAP', **kwargs):
        _method_exists('visualization', vis_method)

        if hasattr(self, 'vis_method') and hasattr(self, 'x_emb_2d'):
            if self.vis_method == vis_method:
                return

        if not hasattr(self, 'x_emb'):
            self.run_step('dim')

        if hasattr(self, 'clu_method'):
            if self.clu_method == "Scanpy":
                ann = anndata.AnnData(X=self.x_emb)
                scanpy.pp.neighbors(ann, n_neighbors=10, n_pcs=40)
                scanpy.tl.umap(ann)
                self.x_emb_2d = ann.obsm['X_umap']
                return
        if hasattr(self, 'align_method'):
            if self.align_method == "Scanpy":
                ann = anndata.AnnData(X=self.x_emb)
                scanpy.pp.neighbors(ann, n_neighbors=10, n_pcs=40)
                scanpy.tl.umap(ann)
                self.x_emb_2d = ann.obsm['X_umap']
                return

        self.vis_method = vis_method
        self.x_emb_2d = wrap("visualization", vis_method)().get(self.x_emb)

    def run_align(self, align_method, x_ref, col_ids_ref, labels_ref,
                  key_maps, **kwargs):
        _method_exists('align', align_method)
        self.align_method = align_method

        self.labels = wrap("align", align_method)().get(
            self.x, self.col_ids, x_ref, col_ids_ref, labels_ref
        ).astype(np.int)
        self.n_clusters = np.unique(self.labels)
        self.key_maps = key_maps.copy()
        for name in key_maps:
            if self.get_cluster_id(name) not in self.n_clusters:
                self.key_maps.pop(name, None)

    def has(self, attr):
        # Needed to make calls from R easier
        if hasattr(self, attr):
            return True
        return False

    def get_cluster_names(self):
        sd = {k: v for k, v in sorted(
            self.key_maps.items(), key=lambda item: item[1])}
        return {'labels': list(sd.values()), 'names': list(sd.keys())}

    def get_label_names(self):
        names = np.zeros_like(self.labels).astype('U200')
        for i in self.key_maps:
            names[np.where(self.labels == self.key_maps[i])] = i
        return names

    def save_session(self):
        sess = {}
        to_save = [
            'dataset', 'dim_method', 'dim_n_components_inp', 'n_components',
            'clu_method', 'eval_method', 'clu_n_clusters_inp', 'labels',
            'n_clusters', 'de_alpha', 'de_n_genes', 'de_correction',
            'markers', 'vis_method', 'key_maps'
        ]

        for attr in to_save:
            if hasattr(self, attr):
                sess[attr] = getattr(self, attr)
        if hasattr(self, 'x_emb'):
            for comp in range(self.x_emb.shape[1]):
                sess['comp' + str(comp)] = self.x_emb[:, comp]
        if hasattr(self, 'x_emb_2d'):
            sess['x1'] = self.x_emb_2d[:, 0]
            sess['x2'] = self.x_emb_2d[:, 1]

        return sess

    def load_session(self, sess):
        to_load = [
            'dataset', 'dim_method', 'dim_n_components_inp', 'n_components',
            'clu_method', 'eval_method', 'clu_n_clusters_inp', 'labels',
            'n_clusters', 'de_alpha', 'de_n_genes', 'de_correction',
            'markers', 'vis_method', 'key_maps'
        ]

        for attr in to_load:
            if attr in sess:
                try:
                    setattr(self, attr, sess[attr])
                except:
                    raise KeyError("Error loading session.")

        if 'x1' in sess:
            self.x_emb_2d = np.vstack([sess['x1'], sess['x2']]).T

        if 'comp0' in sess:
            comps = []
            for comp in range(1000):
                if 'comp' + str(comp) in sess:
                    comps.append(sess['comp' + str(comp)])
                else:
                    break
            self.x_emb = np.vstack(comps).T
