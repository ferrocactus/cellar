import logging

import scanpy as sc
from anndata import AnnData

from ..log import setup_logger
from ._unit import Unit
from ..utils.exceptions import InvalidArgument


DEFAULTS_SCANPY = {
    'filter_cells': {
        'run1': {
            'min_genes': 200,
        },
        'run2': {
            'max_genes': 3000
        }
    },
    'filter_genes': {
        'run1': {
            'min_cells': 3
        }
    },
    'normalize_total': {
        'target_sum': 1e4
    },
    'highly_variable_genes': {
        'min_mean': 0.0125,
        'max_mean': 3,
        'min_disp': 0.5
    },
    'scale': {
        'max_value': 10
    }
}


class Pre_Scanpy(Unit):
    """
    See https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html
    """

    def __init__(self, filter_cells=DEFAULTS_SCANPY['filter_cells'],
                 filter_genes=DEFAULTS_SCANPY['filter_genes'],
                 pct_counts_mt=40,
                 normalize_total=DEFAULTS_SCANPY['normalize_total'],
                 apply_log1p=True,
                 highly_variable_genes=DEFAULTS_SCANPY['highly_variable_genes'],
                 scale=DEFAULTS_SCANPY['scale'], **kwargs):
        """
        Parameters
        __________

        filter_cells: dict of dicts
            Will apply sc.filter_cells for every key in filter_cells.

        filter_genes: dict of dicts
            Will apply sc.filter_genes for every key in filter_genes.

        normalize_total: dict
            Will apply sc.normalize_total using these arguments.

        apply_log1p: bool
            Whether to apply sc.log1p or not.

        highly_variable_genes: dict
            Will apply sc.highly_variable_genes using these arguments.

        scale: dict
            Will apply sc.scale using these arguments.

        **kwargs: dict
            Additional Parameters.
        """

        self.logger = setup_logger('Scanpy')
        self.filter_cells = filter_cells
        self.filter_genes = filter_genes
        self.pct_counts_mt = pct_counts_mt
        self.normalize_total = normalize_total
        self.apply_log1p = apply_log1p
        self.highly_variable_genes = highly_variable_genes
        self.scale = scale
        self.kwargs = kwargs

    def get(self, x):
        """
        Parameters
        __________
        x: AnnData object, shape (n_samples, n_features) or np.ndarray
            Contains the raw data.

        Returns
        _______
        adata: AnnData object, shape (n_samples_p, n_features_p)
            The preprocessed data matrix.
        """

        self.logger.info("Running Scanpy Preprocessing.")

        if not isinstance(x, AnnData):
            adata = AnnData(x)
        else:
            adata = x.copy()

        for key in self.filter_cells:
            sc.pp.filter_cells(
                adata, **self.filter_cells[key], inplace=True)
        for key in self.filter_genes:
            sc.pp.filter_genes(
                adata, **self.filter_genes[key], inplace=True)

        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=['mt'], percent_top=None, inplace=True)

        adata = adata[adata.obs.pct_counts_mt < self.pct_counts_mt, :]

        sc.pp.normalize_total(
            adata, **self.normalize_total, inplace=True)
        if self.apply_log1p:
            sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(
            adata, **self.highly_variable_genes, inplace=True)
        adata = adata[:, adata.var.highly_variable]
        sc.pp.scale(adata, **self.scale)

        return adata
