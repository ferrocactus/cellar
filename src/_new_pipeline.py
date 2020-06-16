from typing import Optional, Union

import numpy as np
from anndata import AnnData
from bidict import bidict

from .units import wrap
from .units import _method_exists
from .utils.validation import _validate_x
from .utils.validation import _validate_clu_n_clusters
from .utils.validation import _validate_n_jobs
from .utils.exceptions import InappropriateArgument
from .utils.exceptions import InvalidArgument


def cluster(
        x: Union[AnnData, np.ndarray, list],
        method: str = 'Leiden',
        eval_method: str = 'Silhouette',
        n_clusters: Union[str, tuple, int, list, np.ndarray] = (3, 6, 1),
        use_emb: bool = True,
        inplace: bool = True,
        n_jobs_multiple: Optional[int] = 1,
        **kwargs) -> Optional[Union[AnnData, np.ndarray]]:
    """
    Parameters
    __________
    x: AnnData object or np array containing the data matrix.

    method: String specifying the clustering method to use. See
        https://github.com/ferrocactus/cellar/tree/master/doc for
        a full list of methods available.

    eval_method: String specifying the clustering method to use.
        See https://github.com/ferrocactus/cellar/tree/master/doc
        for a full list of methods available.

    n_clusters: Can be a single integer of a list/tuple specifying
        multiple number of clusters to use. The evaluation method
        in eval_method will be used to determine the number of
        clusters with the best score. In the case of a tuple,
        a list is formed by treating the tuple as a pythonic range;
        (a, b, c) will start at a, finish at b-1, in increments of c.

    use_emb: If True, will attempt to cluster on an embedding of x
        as specified in x.obsm['x_emb] if x is AnnData object. If x
        is not AnnData, this argument is ignored.

    inplace: Only valid when x is an AnnData object.
        If set to true will update x.obs['labels'], otherwise,
        return a new AnnData object.

    n_jobs_multiple: Only valid when n_clusters results is a list
        of more than one integer. If n_jobs > 1, will run clustering
        for each n_cluster in parallel.

    **kwargs: Additional parameters that will get passed to the
        clustering object as specified in method. For a full list
        see the documentation of the corresponding method.

    Returns
    _______
    If x is an AnnData object, will either return an AnnData object
    or None depending on the value of inplace.
    If x is not AnnData, will return np.ndarray of shape
    (n_observations,) containing the labels.
    """

    # Validations
    is_AnnData = isinstance(x, AnnData)

    if is_AnnData:
        adata = x.copy() if not inplace else x
    else:
        adata = _validate_x(x)
    _method_exists('cluster', method)
    _method_exists('cluster_eval', eval_method)
    n_clusters = _validate_clu_n_clusters(n_clusters, adata.X.shape[0])
    n_jobs_multiple = _validate_n_jobs(n_jobs_multiple)

    # Determine if embeddings should be used or the full data matrix
    if is_AnnData and use_emb:
        if 'x_emb' not in adata.obsm:
            raise ValueError("x_emb not found in AnnData object.")
        x_to_use = adata.obsm['x_emb']
    else:
        x_to_use = adata.X

    # Create clustering object and get labels
    labels = wrap("cluster", method)(
        eval_obj=wrap("cluster_eval", eval_method)(),
        n_clusters=n_clusters,
        n_jobs=n_jobs_multiple, **kwargs).get(x_to_use)

    # If x was a list or numpy array
    if not is_AnnData:
        return labels

    adata.obs['labels'] = labels
    adata.uns['unq_labels'] = np.unique(labels)
    adata.uns['n_clusters'] = len(adata.uns['unq_labels'])
    adata.uns['cluster_method'] = method
    adata.uns['eval_method'] = eval_method
    adata.uns['n_clusters_tested'] = n_clusters
    adata.uns['cluster_names'] = bidict(
        {i: str(i) for i in adata.uns['unq_labels']})

    if not inplace:
        return adata