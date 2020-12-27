import logging

import numpy as np
from numpy.testing import assert_array_equal as aae
from sklearn.decomposition import IncrementalPCA
import scanpy as sc
from anndata import AnnData, read_h5ad

from ..utils.tools import merge_cluster_names
from ..log import setup_logger
from ..utils.exceptions import InvalidArgument
from ._converter import convert
from ._unit import Unit


class Ali_Scanpy_Ingest(Unit):
    """
    See https://scanpy-tutorials.readthedocs.io/en/latest/integrating-data-using-ingest.html
    """

    def __init__(self, **kwargs):
        self.logger = setup_logger("Scanpy")

    def get(self, adata_main, ref, pca_comps=50):
        self.logger.info("Initializing Scanpy ingest method.")
        is_str_ref = isinstance(ref, str)
        is_anndata_ref = isinstance(ref, AnnData)

        col_ids = adata_main.var_names.to_numpy().astype('U')

        if not is_anndata_ref and not is_str_ref:
            raise InvalidArgument("Invalid reference dataset format.")

        if is_str_ref:
            #adata_ref = read_h5ad(ref, backed='r')
            adata_ref = read_h5ad(ref)
        else:
            adata_ref = ref.copy()

        if 'labels' not in adata_ref.obs:
            raise InvalidArgument("No labels found in reference dataset.")

        d_main = convert(col_ids)
        d_ref = convert(adata_ref.var_names.to_numpy().astype('U'))

        adata_main.var_names = d_main['names']
        adata_main.var_names_make_unique()
        adata_ref.var_names = d_ref['names']
        adata_ref.var_names_make_unique()

        common_genes, comm1, comm2 = np.intersect1d(
            d_main['names'], d_ref['names'], return_indices=True)

        adata_main = adata_main[:, common_genes]
        adata_ref = adata_ref[:, common_genes]

        #common_genes_mask = np.isin(adata_ref.var_names.to_numpy().astype('U'), common_genes)
        #assert common_genes_mask.shape[0] == adata_ref.shape[1]

        # aae(adata_main.var_names.to_numpy().astype('U'),
        #     adata_ref.var_names.to_numpy().astype('U'))

        self.logger.info(f'Found {len(common_genes)} genes in common.')

        self.logger.info('Running PCA on main data.')
        sc.pp.pca(adata_main, n_comps=pca_comps)
        sc.pp.neighbors(adata_main)
        self.logger.info('Running UMAP on main data.')
        sc.tl.umap(adata_main)

        # ipca = IncrementalPCA(pca_comps)
        # i = 0
        # self.logger.info('Running Incremental PCA on reference data.')
        # print(adata_ref.shape)
        # while i < adata_ref.shape[0]:
        #     upp_lim = min(i+1000, adata_ref.shape[0])
        #     ipca.partial_fit(adata_ref.X[i:upp_lim][:, comm2])
        #     i += 1000
        #     self.logger.info(f'Incremental PCA :: Finished Epoch {int(i/1000)}')

        # i = 0
        # x_ref_pca = np.zeros((adata_ref.shape[0], pca_comps))

        # self.logger.info('Transforming data.')
        # while i < adata_ref.shape[0]:
        #     upp_lim = min(i+1000, adata_ref.shape[0])
        #     x_ref_pca[i:upp_lim] = ipca.transform(adata_ref.X[i:upp_lim][:, comm2])
        #     i += 1000

        # adata_ref = adata_ref[:, common_genes]
        # adata_ref.obsm['X_pca'] = x_ref_pca

        self.logger.info('Running PCA on reference data.')
        sc.pp.pca(adata_ref, n_comps=pca_comps)
        sc.pp.neighbors(adata_ref)
        self.logger.info('Running UMAP on reference data.')
        sc.tl.umap(adata_ref)

        self.logger.info('Running Scanpy Ingest.')
        sc.tl.ingest(adata_main, adata_ref, obs='labels')

        merge_cluster_names(adata_main, adata_ref)

        return np.array(adata_main.obs['labels']).astype(np.int).reshape(-1), adata_main.uns['cluster_names']


# class Ali_Scanpy_Ingest2(Unit):
#     """
#     See https://scanpy-tutorials.readthedocs.io/en/latest/integrating-data-using-ingest.html
#     """

#     def __init__(self, **kwargs):
#         self.logger = setup_logger("Scanpy")

#     def get(self, x, col_ids, x_ref, col_ids_ref, labels_ref):
#         self.logger.info("Initializing Scanpy ingest method.")
#         # Consider only genes in common
#         common_col_ids = np.intersect1d(col_ids, col_ids_ref)
#         self.logger.info(f"Found {len(common_col_ids)} genes in common.")

#         x = x[:, np.isin(col_ids, common_col_ids)]
#         x_ref = x_ref[:, np.isin(col_ids_ref, common_col_ids)]

#         a1 = AnnData(x)
#         a2 = AnnData(x_ref)
#         sc.pp.pca(a1)
#         sc.pp.pca(a2)
#         sc.pp.neighbors(a1)
#         sc.pp.neighbors(a2)
#         sc.tl.umap(a1)
#         sc.tl.umap(a2)
#         a2.obs['labels'] = labels_ref

#         sc.tl.ingest(a1, a2, obs='labels')
#         return np.array(a1.obs['labels']).astype(np.int).reshape(-1)
