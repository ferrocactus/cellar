from typing import Optional, Union

import numpy as np
import plotly.express as px
from plotly.graph_objects import Figure
from anndata import AnnData

from .units import wrap
from . import reduce_dim_vis
from . import name_genes


COLORS = [
    '#cc5151', '#51cccc', '#337f7f', '#8ecc51', '#7f3333', '#597f33', '#8e51cc',
    '#59337f', '#ccad51', '#7f6c33', '#51cc70', '#337f46', '#5170cc', '#33467f',
    '#cc51ad', '#7f336c', '#cc7f51', '#7f4f33', '#bccc51', '#757f33', '#60cc51',
    '#3c7f33', '#51cc9e', '#337f62', '#519ecc', '#33627f', '#6051cc', '#3c337f'
]


def _find_gene_index(adata, gene):
    if isinstance(gene, int):
        if gene > adata.var.index.to_numpy().shape[0]:
            raise ValueError("Index out of bounds.")
        return gene
    if not isinstance(gene, str):
        raise ValueError("Incorrect gene format.")
    if gene in adata.var.index.to_numpy():
        return np.where(adata.var.index.to_numpy() == gene)[0]
    if 'parsed_names' not in adata.var:
        name_genes(adata)
    if gene in adata.var['parsed_names'].to_numpy():
        return np.where(adata.var['parsed_names'] == gene)[0]
    if gene in adata.var['parsed_ids'].to_numpy():
        return np.where(adata.var['parsed_ids'] == gene)[0]

    return -1


def _plot_labels(
    adata: AnnData,
    show_title: Optional[bool] = False,
    return_fig: Optional[bool] = False
) -> Figure:
    """
    Helper function for plot.
    """
    has_labels = True
    if 'labels' not in adata.obs:
        has_labels = False
        print("Labels not found. Plotting 2d embeddings.")
        #raise ValueError("Labels not found in object.")
    if 'x_emb_2d' not in adata.obsm:
        print("2d embeddings not found.")
        print("Running default visualization method.")
        reduce_dim_vis(adata)
    if has_labels:
        color = adata.obs['labels'].to_numpy().astype(str)
    method = adata.uns['visualization_info_2d']['method']
    fig = px.scatter(
        x=adata.obsm['x_emb_2d'][:, 0],
        y=adata.obsm['x_emb_2d'][:, 1],
        color=color if has_labels else None,
        hover_data={'Cell': adata.obs.index.to_numpy()},
        labels={
            'x': f'{method}1',
            'y': f'{method}2'
        },
        title=adata.uns['dataset'] if show_title else None,
        template='none'
    )
    if return_fig:
        return fig
    fig.show()


def _plot_gene(
    adata: AnnData,
    gene: Optional[Union[str, int]] = None,
    show_title: Optional[bool] = False,
    return_fig: Optional[bool] = False
) -> Figure:
    """
    Helper function for plot.
    """
    if gene is None:
        raise ValueError("Please specify gene to plot.")
    index = _find_gene_index(adata, gene)
    if index == -1:
        print("Gene not found.")
        return
    color = adata.X[:, index]
    method = adata.uns['visualization_info_2d']['method']
    fig = px.scatter(
        x=adata.obsm['x_emb_2d'][:, 0],
        y=adata.obsm['x_emb_2d'][:, 1],
        color=color,
        hover_data={'Cell': adata.obs.index.to_numpy()},
        labels={
            'x': f'{method}1',
            'y': f'{method}2'
        },
        title=adata.uns['dataset'] if show_title else None,
        template='none'
    )
    if return_fig:
        return fig
    fig.show()


def _plot_scores(
    adata: AnnData,
    show_title: Optional[bool] = False,
    return_fig: Optional[bool] = False
) -> Figure:
    """
    Helper function for plot.
    """
    if 'scores' not in adata.uns['cluster_info']:
        raise ValueError("Scores not found in object.")

    eval_method = adata.uns['cluster_info']['eval_method']
    fig = px.line(
        x=adata.uns['cluster_info']['n_clusters_used'],
        y=adata.uns['cluster_info']['scores'],
        labels={
            'x': 'n_clusters',
            'y': f'{eval_method} score'
        },
        title=adata.uns['dataset'] if show_title else None,
        template='none'
    )
    if return_fig:
        return_fig
    fig.show()


def plot(
    x: AnnData,
    by: Optional[str] = None,
    gene: Optional[Union[str, int]] = None,
    show_title: Optional[bool] = False,
    return_fig: Optional[bool] = False
) -> None:
    """
    Plotting functionality.

    Parameters
    __________
    x: AnnData object containing the data matrix and the plot keys.

    by: String specifying what to plot.

    gene: Will be used only if by is None or by == 'gene'.
        Specify the name of the gene for which to plot expression for.
        Can be in ensembl format, gene name, or index. If name is
        specified and names are not found in adata, then will run
        name_genes, but will not save the names in the adata object.

    show_title: Boolean specifying whether to show the name of the
        dataset in the plot.

    return_fig: Boolean specifying whether to return a fig object if
        set to True, otherwise will plot immediately.
    """
    # Validations
    is_AnnData = isinstance(x, AnnData)
    if not is_AnnData:
        raise ValueError("Object not in AnnData format.")

    if by == 'labels' or (by is None and gene is None):
        return _plot_labels(x, show_title, return_fig)
    elif by is None or by == 'gene':
        return _plot_gene(x, gene, show_title, return_fig)
    elif by == 'scores':
        return _plot_scores(x, show_title, return_fig)
