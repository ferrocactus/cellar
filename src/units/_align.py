import logging

import anndata
import scanpy as sc
from ..log import setup_logger
from ._unit import Unit
import numpy as np

class Ali_Scanpy_Ingest(Unit):
    """
    See https://scanpy-tutorials.readthedocs.io/en/latest/integrating-data-using-ingest.html
    """
    def __init__(self, **kwargs):
        self.logger = setup_logger("Scanpy")

    def get(self, x, col_ids, x_ref, col_ids_ref, labels_ref):
        self.logger.info("Initializing Scanpy ingest method.")
        # Consider only genes in common
        common_col_ids = np.intersect1d(col_ids, col_ids_ref)
        self.logger.info(f"Found {len(common_col_ids)} genes in common.")

        x = x[:, np.isin(col_ids, common_col_ids)]
        x_ref = x_ref[:, np.isin(col_ids_ref, common_col_ids)]
        a1 = anndata.AnnData(x)
        a2 = anndata.AnnData(x_ref)
        sc.pp.pca(a1)
        sc.pp.pca(a2)
        sc.pp.neighbors(a1)
        sc.pp.neighbors(a2)
        sc.tl.umap(a1)
        sc.tl.umap(a2)
        a2.obs['labels'] = labels_ref

        sc.tl.ingest(a1, a2, obs='labels')
        return a1.obs['labels']
