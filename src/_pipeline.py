from typing import Optional, Union

import numpy as np
from anndata import AnnData
from bidict import bidict

from .units import wrap
from .units import _method_exists
from .units import convert

from .utils.tools import _emb_exists_in_adata
from .utils.tools import _labels_exist_in_adata
from .utils.tools import _2d_emb_exists_in_adata
from .utils.tools import populate_subsets


from .utils.validation import _validate_clu_n_clusters
from .utils.validation import _validate_cluster_list
from .utils.validation import _validate_dim_n_components
from .utils.validation import _validate_ensemble_methods
from .utils.validation import _validate_mark_alpha
from .utils.validation import _validate_mark_correction
from .utils.validation import _validate_mark_markers_n
from .utils.validation import _validate_n_jobs
from .utils.validation import _validate_subset
from .utils.validation import _validate_uncertain_subset
from .utils.validation import _validate_x

from .utils.exceptions import InappropriateArgument
from .utils.exceptions import InvalidArgument


def reduce_dim(
        x: Union[AnnData, np.ndarray, list],
        method: str = 'PCA',
        n_components: Union[str, int, float] = 'knee',
        inplace: Optional[bool] = True,
        check_if_exists: Optional[bool] = False,
        clear_2d_emb: Optional[bool] = True,
        **kwargs) -> Optional[Union[AnnData, np.ndarray]]:
    """
    Reduce dimensionality of the data.

    Parameters
    __________
    x: AnnData object or np array containing the data matrix.

    method: String specifying the dimensionality reduction method
        to use. See https://github.com/ferrocactus/cellar/tree/master/doc
        for a full list of methods available.

    n_components: Number of components to use. If method is 'PCA', this
        parameters can also be 'knee' in which case the number of
        components will be determined automatically based on the
        ankle heuristic.

    inplace: Only used when x is an AnnData object.
        If set to true will update x.obsm['x_emb'], otherwise,
        return a new AnnData object.

    check_if_exists: Only used when x is an AnnData object.
        If set to true, will check x if it contains embeddings
        from the same method. If found, and the number of components
        in x.obsm['x_emb'] is greater than or equal to n_components,
        then the first n_components columns of x.obsm['x_emb']
        will be returned. No computation takes place.

    clear_2d_emb: If set to true and x_emb changes, then will also
        clear f'x_emb_{dim}d' if any exist.

    **kwargs: Additional parameters that will get passed to the
        clustering object as specified in method. For a full list
        see the documentation of the corresponding method.

    Returns
    _______
    If x is an AnnData object, will either return an AnnData object
    or None depending on the value of inplace.
    If x is not AnnData, will return np.ndarray of shape
    (n_observations, n_components) containing the embedding.
    """
    # Validations
    is_AnnData = isinstance(x, AnnData)

    if is_AnnData:
        adata = x.copy() if not inplace else x
    else:
        adata = _validate_x(x)
    _method_exists('dim_reduction', method)
    n_components_used = n_components
    n_components = _validate_dim_n_components(n_components,
                                              method, *adata.X.shape)

    x_emb = None
    # Check if embeddings exist in x to save computation
    if is_AnnData and check_if_exists:
        x_emb = _emb_exists_in_adata(adata, method, n_components)

    # If x_emb was not found in adata, it is None
    # Create dimensionality reduction object and find embedding
    if x_emb is None:
        x_emb = wrap('dim_reduction', method)(
            n_components=n_components, **kwargs).get(adata.X)
        # Clear visualization if it used previous embeddings
        for dim in [2, 3]:
            if f'x_emb_{dim}d' in adata.obsm and clear_2d_emb:
                if adata.uns[f'visualization_info_{dim}d']['used_emb'] == True:
                    print('Clearing 2d embeddings...')
                    adata.obsm.pop(f'x_emb_{dim}d', None)
                    adata.uns.pop(f'visualization_info_{dim}d', None)

    if not is_AnnData:
        return x_emb

    # Populate
    adata.obsm['x_emb'] = x_emb
    adata.uns['dim_reduction_info'] = {}
    adata.uns['dim_reduction_info']['method'] = method
    adata.uns['dim_reduction_info']['n_components'] = x_emb.shape[1]
    if n_components_used == 'knee':
        n_components_used = 'Automatic'
    adata.uns['dim_reduction_info']['n_components_used'] = n_components_used
    adata.uns['dim_reduction_info']['kwargs'] = kwargs

    if not inplace:
        return adata


def cluster(
        x: Union[AnnData, np.ndarray, list],
        method: str = 'Leiden',
        eval_method: str = 'Silhouette',
        n_clusters: Union[str, tuple, int, list, np.ndarray] = (3, 6, 1),
        use_emb: bool = True,
        inplace: Optional[bool] = True,
        check_if_exists: Optional[bool] = False,
        n_jobs_multiple: Optional[int] = 1,
        **kwargs) -> Optional[Union[AnnData, np.ndarray]]:
    """
    Cluster the data.

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
        Ignored if method is 'Leiden' or 'Scanpy'.

    use_emb: If True, will attempt to cluster on an embedding of x
        as specified in x.obsm['x_emb'] if x is AnnData object. If x
        is not AnnData, this argument is ignored.

    inplace: Only used when x is an AnnData object.
        If set to true will update x.obs['labels'], otherwise,
        return a new AnnData object.

    check_if_exists: Only used when x is an AnnData object.
        If set to true, will check if x contains labels
        from the same method.

    n_jobs_multiple: Only used when n_clusters results is a list
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
    if method not in ['Leiden', 'Scanpy']:
        n_clusters = _validate_clu_n_clusters(n_clusters, adata.X.shape[0])
    else:
        n_clusters = 'NA'

    # Check if labels exist and return if they do
    if _labels_exist_in_adata(adata, method, n_clusters):
        if not inplace:
            return adata
        return

    n_jobs_multiple = _validate_n_jobs(n_jobs_multiple)

    # Determine if embeddings should be used or the full data matrix
    if is_AnnData and use_emb:
        if 'x_emb' not in adata.obsm:
            raise InvalidArgument("x_emb not found in AnnData object.")
        x_to_use = adata.obsm['x_emb']
    else:
        x_to_use = adata.X

    # Create clustering object and get labels
    labels, scores = wrap("cluster", method)(
        eval_obj=wrap("cluster_eval", eval_method)(),
        n_clusters=n_clusters,
        n_jobs=n_jobs_multiple, **kwargs).get(x_to_use)

    labels = np.array(labels).astype(int)
    scores = np.array(scores).astype(float).reshape(-1)

    # If x was a list or numpy array
    if not is_AnnData:
        return labels

    # Populate entries
    adata.obs['labels'] = labels
    adata.uns['cluster_info'] = {}
    unq_labels = np.unique(labels)
    adata.uns['cluster_info']['unique_labels'] = unq_labels
    adata.uns['cluster_info']['n_clusters'] = len(unq_labels)
    adata.uns['cluster_info']['method'] = method
    adata.uns['cluster_info']['n_clusters_used'] = n_clusters
    if method in ['Leiden', 'Scanpy']:
        eval_method = "NA"
        scores = "NA"
    adata.uns['cluster_info']['eval_method'] = eval_method
    adata.uns['cluster_info']['scores'] = scores
    adata.uns['cluster_info']['used_emb'] = use_emb
    adata.uns['cluster_info']['kwargs'] = kwargs
    adata.uns['cluster_names'] = bidict(
        {i: str(i) for i in unq_labels})

    populate_subsets(adata)

    if not inplace:
        return adata


def reduce_dim_vis(
        x: Union[AnnData, np.ndarray, list],
        method: str = 'UMAP',
        dim: int = 2,
        use_emb: bool = True,
        inplace: Optional[bool] = True,
        check_if_exists: Optional[bool] = False,
        **kwargs) -> Optional[Union[AnnData, np.ndarray]]:
    """
    Reduce dimensionality of the data for visualization.

    Parameters
    __________
    x: AnnData object or np array containing the data matrix.

    method: String specifying the dimensionality reduction method
        to use. See https://github.com/ferrocactus/cellar/tree/master/doc
        for a full list of methods available.

    dim: Dimensionality, can be either 2 or 3.

    use_emb: If True, will attempt to run on an embedding of x
        as specified in x.obsm['x_emb] if x is AnnData object. If x
        is not AnnData, this argument is ignored.

    inplace: Only used when x is an AnnData object.
        If set to true will update x.obsm['x_emb_2d'] if dim=2,
        otherwise, x.obsm['x_emb_3d'] if dim=3, otherwise
        return a new AnnData object with the appropriate
        key added.

    check_if_exists: Only used when x is an AnnData object.
        If set to true, will check x if it contains 2d/3d embeddings
        from the same method. If found, return those instead.

    **kwargs: Additional parameters that will get passed to the
        clustering object as specified in method. For a full list
        see the documentation of the corresponding method.

    Returns
    _______
    If x is an AnnData object, will either return an AnnData object
    or None depending on the value of inplace.
    If x is not AnnData, will return np.ndarray of shape
    (n_observations, dim) containing the visualization embedding.
    """
    # Validations
    is_AnnData = isinstance(x, AnnData)

    if is_AnnData:
        adata = x.copy() if not inplace else x
    else:
        adata = _validate_x(x)
    _method_exists('visualization', method)

    dim = int(dim)
    if dim != 2 and dim != 3:
        raise InvalidArgument("Incorrect number of dimensions specified.")

    # Check if embedding already exists in adata
    if is_AnnData and check_if_exists:
        if _2d_emb_exists_in_adata(adata, method, dim, use_emb):
            if not inplace:
                return adata
            return

    # Determine if embeddings should be used or the full data matrix
    if is_AnnData and use_emb:
        if 'x_emb' not in adata.obsm:
            print("x_emb not found. Running PCA with default parameters.")
            print("If you don't want to use embeddings, pass use_emb=False.")
            reduce_dim(adata, inplace=True)
        x_to_use = adata.obsm['x_emb']
    else:
        x_to_use = adata.X

    # Create dimensionality reduction object and find embedding
    x_emb_nd = wrap("visualization", method)(
        n_components=dim, **kwargs).get(x_to_use)

    if not is_AnnData:
        return x_emb_nd

    # Update the 2d or 3d keys
    adata.obsm[f'x_emb_{dim}d'] = x_emb_nd
    adata.uns[f'visualization_info_{dim}d'] = {}
    adata.uns[f'visualization_info_{dim}d']['method'] = method
    adata.uns[f'visualization_info_{dim}d']['used_emb'] = use_emb
    adata.uns[f'visualization_info_{dim}d']['kwargs'] = kwargs

    if not inplace:
        return adata


def name_genes(
    x: Union[AnnData, np.ndarray, list],
    inplace: Optional[bool] = True
) -> Optional[Union[AnnData, np.ndarray]]:
    """
    Find synonyms for gene IDs/names.

    Parameters
    __________
    x: AnnData object or np array containing names to be converted.

    inplace: Only used when x is an AnnData object.
        If set to true will update x.var['names'] if var.index is
        determined to consist of gene ensembl ids, or will update
        x.var['ids'] if var.index is determined to consist of gene names.
        In either case, a new column with the other name is crated for
        consistency.

    Returns
    _______
    If x is an AnnData object, will either return an AnnData object
    or None depending on the value of inplace.
    If x is not AnnData, will return a dictionary with two keys
    'ids' and 'names'.
    """
    # Validations
    is_AnnData = isinstance(x, AnnData)

    if is_AnnData:
        adata = x.copy() if not inplace else x
    else:
        d = convert(np.array(x))
        return d

    d = convert(adata.var.index.to_numpy().astype('U'))
    for key in d:
        adata.var['parsed_' + key] = d[key]

    if not inplace:
        return adata


def de(
        x: Union[AnnData, np.ndarray, list],
        subset1: Union[str, np.ndarray, list],
        subset2: Optional[Union[str, np.ndarray, list]] = None,
        method: str = 'TTest',
        alpha: float = 0.05,
        max_n_genes: int = 200,
        correction: str = 'holm-sidak',
        inplace: Optional[bool] = True,
        uncertain: Optional[bool] = False,
        **kwargs) -> Optional[Union[AnnData, dict]]:
    """
    Run differential expression.

    Parameters
    __________
    x: AnnData object or np array containing the data matrix.

    method: String specifying the de method to use. See
        https://github.com/ferrocactus/cellar/tree/master/doc
        for a full list of methods available.

    subset1: An array consisting of the indices of the cells
        for which de genes need to be found, or a string
        specifying the subset as stored in adata.uns['subsets].
        If string, x has to be AnnData object.

    subset2: Same format at subset1. If set to None, then will
        run differential expression of cells in subset1 vs all,
        otherwise, will run de of subset1 vs subset2.

    alpha: Cutoff for p values. Default to 0.05.

    max_n_genes: Maximum number of differential genes to look for.
        The returned number may be smaller if not enough p values
        are found to be less then alpha.

    correction: String specifying the correction method to
        use for p values. For a full list see
        https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html
        Can also be set to 'None' (string), in which case no correction
        is performed.

    inplace: Only used when x is an AnnData object.
        If set to true will update x.uns['de'].

    **kwargs: Additional parameters.

    Returns
    _______
    If x is an AnnData object, will either return an AnnData object
    or None depending on the value of inplace.
    If x is not AnnData, will return a dictionary containing all the
    de information.
    """
    # Validations
    is_AnnData = isinstance(x, AnnData)

    if is_AnnData:
        adata = x.copy() if not inplace else x
    else:
        adata = _validate_x(x)
    _method_exists('de', method)
    if (uncertain):
        subset1 = _validate_uncertain_subset(subset1, adata)
        subset2 = _validate_uncertain_subset(subset2, adata)
    else:
        subset1 = _validate_subset(subset1, adata)
        subset2 = _validate_subset(subset2, adata)
    alpha = _validate_mark_alpha(alpha)
    max_n_genes = _validate_mark_markers_n(max_n_genes, adata.X.shape[1])
    correction = _validate_mark_correction(correction)

    if subset2 is None:
        # Run subset1 vs all
        de_genes_d = wrap("de", method)(
            alpha=alpha,
            markers_n=max_n_genes,
            correction=correction).get_subset(adata.X, subset1)
    else:
        # Run subset1 vs subset2
        # First, artificially create dataset and labels
        x1 = adata.X[subset1]
        x2 = adata.X[subset2]
        x = np.concatenate([x1, x2])
        labels1 = np.zeros((x1.shape[0],), dtype=int)
        labels2 = np.ones((x2.shape[0],), dtype=int)
        labels = np.concatenate([labels1, labels2])

        de_genes_d = wrap("de", method)(
            alpha=alpha,
            markers_n=max_n_genes,
            correction=correction).get(x, labels)

    # Assign first dict key to anndata object
    # This is the only key if subset2 is None
    adata.uns['de'] = de_genes_d[next(iter(de_genes_d))]
    adata.uns['de_info'] = {}
    adata.uns['de_info']['method'] = method
    adata.uns['de_info']['alpha'] = alpha
    adata.uns['de_info']['max_n_genes'] = max_n_genes
    adata.uns['de_info']['correction'] = correction
    adata.uns['de_info']['subset1'] = subset1
    adata.uns['de_info']['subset2'] = subset2
    adata.uns['de_info']['kwargs'] = kwargs

    if not inplace:
        return adata



def ss_cluster(
        x: AnnData,
        method: str = 'SeededKMeans',
        use_emb: bool = True,
        inplace: Optional[bool] = True,
        preserved_labels: Optional[Union[np.ndarray, list, int]] = None,
        **kwargs) -> Optional[Union[AnnData, np.ndarray]]:
    """
    Semi-supervised clustering of the data.

    Parameters
    __________
    x: AnnData object with x.var['labels'] populated.

    method: String specifying the ss clustering method to use. See
        https://github.com/ferrocactus/cellar/tree/master/doc for
        a full list of methods available.

    use_emb: If True, will attempt to cluster on an embedding of x
        as specified in x.obsm['x_emb'].

    inplace: If set to true will update x.obs['labels'], otherwise,
        return a new AnnData object.

    preserved_labels: A list of labels that should be preserved
        by the algorithm.

    **kwargs: Additional parameters that will get passed to the
        clustering object as specified in method. For a full list
        see the documentation of the corresponding method.

    Returns
    _______
    Will either return an AnnData object
    or None depending on the value of inplace.
    """
    # Validations
    is_AnnData = isinstance(x, AnnData)
    if not is_AnnData:
        raise InvalidArgument("x is not in AnnData format.")
    if 'labels' not in x.obs:
        raise InvalidArgument("'labels' not found in object.")
    adata = x.copy() if not inplace else x
    _method_exists('ss_cluster', method)
    if method != 'SeededKMeans':
        preserved_labels = _validate_cluster_list(
            adata.obs['labels'],
            preserved_labels)
    else:
        preserved_labels = np.array([])

    # Determine if embeddings should be used or the full data matrix
    if use_emb:
        if 'x_emb' not in adata.obsm:
            raise InvalidArgument("x_emb not found in AnnData object.")
        x_to_use = adata.obsm['x_emb']

    labels_cp = adata.obs['labels'].to_numpy().copy()
    cluster_names = bidict({})
    unq_labels = np.unique(labels_cp)
    n_clusters = len(unq_labels)
    preserved_labels_upd = []

    for i, label in enumerate(unq_labels):
        labels_cp[adata.obs['labels'] == label] = i
        if label in preserved_labels:
            cluster_names[i] = adata.uns['cluster_names'][label]
            preserved_labels_upd.append(i)

    preserved_labels_upd = np.array(preserved_labels_upd)
    # Create clustering object and get labels
    labels = wrap("ss_cluster", method)(**kwargs).get(
        x_to_use, labels_cp, preserved_labels=preserved_labels_upd)

    labels = np.array(labels).astype(int)
    unq_labels = np.unique(labels)

    for label in unq_labels:
        if label not in list(cluster_names.keys()):
            cluster_names[label] = str(label)

    # Populate entries
    adata.obs['labels'] = labels
    adata.uns['cluster_info'] = {}
    adata.uns['cluster_info']['unique_labels'] = unq_labels
    adata.uns['cluster_info']['n_clusters'] = len(unq_labels)
    adata.uns['cluster_info']['method'] = method
    adata.uns['cluster_info']['n_clusters_used'] = n_clusters
    adata.uns['cluster_info']['used_emb'] = use_emb
    adata.uns['cluster_info']['kwargs'] = kwargs
    adata.uns['cluster_names'] = cluster_names

    populate_subsets(adata)

    if not inplace:
        return adata


def transfer_labels(
        x: AnnData,
        ref: AnnData,
        method: str = 'Scanpy Ingest',
        inplace: Optional[bool] = True,
        **kwargs) -> Optional[Union[AnnData, np.ndarray]]:
    """
    Transfer labels using a reference dataset.

    Parameters
    __________
    x: AnnData object containing the data.

    ref: AnnData object with x.var['labels'] populated.

    method: String specifying the label transfer method to use. See
        https://github.com/ferrocactus/cellar/tree/master/doc for
        a full list of methods available.

    inplace: If set to true will update x.obs['labels'], otherwise,
        return a new AnnData object.

    **kwargs: Additional parameters that will get passed to the
        clustering object as specified in method. For a full list
        see the documentation of the corresponding method.

    Returns
    _______
    If x is an AnnData object, will either return an AnnData object
    or None depending on the value of inplace.
    """
    # Validations
    is_AnnData = isinstance(x, AnnData) and isinstance(ref, AnnData)
    if not is_AnnData:
        raise InvalidArgument("x is not in AnnData format.")
    adata = x.copy() if not inplace else x

    if 'labels' not in ref.obs:
        raise InvalidArgument("labels not found in reference dataset.")

    _method_exists('align', method)

    # Create alignment object and get labels
    labels = wrap("align", method)().get(
        adata.X, adata.var.index.to_numpy().astype('U'),
        ref.X, ref.var.index.to_numpy().astype('U'),
        ref.obs['labels'].to_numpy().astype(np.int)
    ).astype(np.int)

    # Populate entries
    adata.obs['labels'] = labels
    adata.uns['cluster_info'] = {}
    unq_labels = np.unique(labels)
    adata.uns['cluster_info']['unique_labels'] = unq_labels
    adata.uns['cluster_info']['n_clusters'] = len(unq_labels)
    adata.uns['cluster_info']['method'] = method
    adata.uns['cluster_info']['kwargs'] = kwargs
    adata.uns['cluster_names'] = bidict(
        {i: str(i) for i in unq_labels})

    populate_subsets(adata)

    if not inplace:
        return adata


