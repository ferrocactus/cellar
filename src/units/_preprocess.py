import logging
from typing import Optional, Union

import scanpy as sc
import numpy as np
from anndata import AnnData

from ..log import setup_logger
from ._unit import Unit
from ..utils.exceptions import InvalidArgument
from ..utils.validation import _validate_unprocessed_x
from ..utils.validation import _validate_atac_operation
from ..utils.validation import _validate_interval_extension
from ..methods import BinToGene


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

        _validate_unprocessed_x(adata.X)

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


class Pre_ATAC(Unit):
    """
    Converts a cell-by-bin matrix to cell-by-gene and performs basic
    normalization.

    See https://github.com/ferrocactus/BinToGene
    """

    def __init__(self,
                 gencode_path='gencode_v34_genes_protein_coding.csv',
                 operation: str = 'sum',
                 extend: Optional[Union[str, int]] = '5x',
                 max_extend: Optional[Union[str, int]] = 50000,
                 stream_direction: bool = True,
                 op_extend: Optional[Union[str, int]] = '1x',
                 max_op_extend: Optional[Union[str, int]] = 10000,
                 n_jobs: Optional[int] = 4,
                 normalize_total=DEFAULTS_SCANPY['normalize_total'],
                 apply_log1p=True,
                 highly_variable_genes=DEFAULTS_SCANPY['highly_variable_genes'],
                 scale=DEFAULTS_SCANPY['scale'], **kwargs):
        """
        Parameters:
        ___________
        gencode_path: string specifying the path of the gencode file

        operation: Can be 'sum' or 'mean'. Specifies whether to add bins
            that intersect with the gene or take their mean instead.

        extend: String ending in x specifying a multiplier of the gene length
            to add to its extremes, or int specifying exact number of base
            pairs to add.

        max_extend: Same format as extend, but will be used as a threshold.

        stream_direction: Specifies whether to take stream direction
            into consideration. If set to True, then op_stream_extend specifies
            the value that will be added to the opposite stream end.

        op_extend: Same format as extend, but applied to the opposite
            stream end only. Ignored if stream_direction set to False.

        max_op_extend: Same as max_extend applied to op_stream_extend.

        n_jobs: If None, 0, or 1 will use 1 job. Otherwise use multithreading.

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
        self.logger = setup_logger('BinToGene')

        self.gencode_path = gencode_path
        self.operation = _validate_atac_operation(operation)
        self.extend = _validate_interval_extension(extend)
        self.max_extend = _validate_interval_extension(max_extend)
        self.stream_direction = stream_direction
        self.op_extend = _validate_interval_extension(op_extend)
        self.max_op_extend = _validate_interval_extension(max_op_extend)
        self.n_jobs = n_jobs
        self.normalize_total = normalize_total
        self.apply_log1p = apply_log1p
        self.highly_variable_genes = highly_variable_genes
        self.scale = scale
        self.kwargs = kwargs

    def get(self, adata: AnnData):
        """
        Given an AnnData object whose X is a cell-by-bin matrix, convert it
        to cell-by-gene and return new anndata object.
        """
        if not isinstance(adata, AnnData):
            raise InvalidArgument("Object passed is not AnnData.")

        btg = BinToGene(
            gencode_path=self.gencode_path,
            operation=self.operation,
            extend=self.extend,
            max_extend=self.max_extend,
            stream_direction=self.stream_direction,
            op_extend=self.op_extend,
            max_op_extend=self.max_op_extend,
            n_jobs=self.n_jobs)

        counts, ids = btg.convert(adata.X, adata.var_names)

        cbg = AnnData(counts)
        cbg.var_names = ids
        cbg.obs_names = adata.obs_names.copy()

        self.logger.info(f"Bin to gene matrix has shape {adata.shape}.")

        sc.pp.normalize_total(cbg, **self.normalize_total, inplace=True)
        if self.apply_log1p:
            sc.pp.log1p(cbg)
        sc.pp.highly_variable_genes(
            cbg, **self.highly_variable_genes, inplace=True)
        cbg = cbg[:, cbg.var.highly_variable]
        sc.pp.scale(cbg, **self.scale)

        return cbg
