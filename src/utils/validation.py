from ast import literal_eval

import numpy as np


def _effective_n_jobs(n_jobs):
    if isinstance(n_jobs, str) and n_jobs is not None:
        n_jobs = literal_eval(n_jobs)

    if n_jobs is None or n_jobs == 1:
        return 1
    else:
        return n_jobs


def _effective_numerical_value(x):
    if isinstance(x, str) and x is not None:
        return literal_eval(x)
    else:
        return x


def _effective_n_clusters(n_clusters):
    if isinstance(n_clusters, str):
        n_clusters = literal_eval(n_clusters)

    if isinstance(n_clusters, int):
        return n_clusters, 'int'
    elif isinstance(n_clusters, list) or isinstance(n_clusters, np.ndarray):
        if len(n_clusters) == 1:
            return n_clusters[0], 'int'
        elif len(n_clusters) == 0:
            raise ValueError("Incorrect list of clusters. Got empty list.")
        else:
            return n_clusters, 'list'
    elif isinstance(n_clusters, tuple):
        k_list = list(range(*n_clusters))
        if len(k_list) == 1:
            return k_list[0], 'int'
        elif len(k_list) == 0:
            raise ValueError("Incorrect range of clusters. Got empty range.")
        else:
            return k_list, 'list'
    else:
        raise ValueError("Incorrect format for clusters specified.")


def _validate_dim_n_components(dim_n_components, h, w):
    if isinstance(dim_n_components, str):
        if dim_n_components != 'knee':
            try:
                dim_n_components = literal_eval(dim_n_components)
            except:
                raise ValueError("Incorrect number of components specified.")
        else:
            return dim_n_components

    if isinstance(dim_n_components, int):
        if dim_n_components >= h or dim_n_components >= w:
            raise ValueError("Number of components needs to be less than "
                             f"min(width, height)={min(w, h)}.")
        if dim_n_components < 1:
            raise ValueError(
                "Number of components needs to be greater than 0.")
    else:
        raise ValueError("Incorrect number of components specified.")
    return dim_n_components


def _validate_n_clusters(n_clusters, h):
    if n_clusters >= h:
        raise ValueError("Number of clusters needs to be less than "
                         "the number of samples.")
    if n_clusters < 2:
        raise ValueError("Number of clusters needs to be greater "
                         "than 2.")


def _validate_clu_n_clusters(clu_n_clusters, h):
    if isinstance(clu_n_clusters, str):
        try:
            clu_n_clusters = literal_eval(clu_n_clusters)
        except:
            raise ValueError("Incorrect format for the number of clusters.")


    if isinstance(clu_n_clusters, float):
        clu_n_clusters = int(clu_n_clusters)

    if isinstance(clu_n_clusters, int):
        _validate_n_clusters(clu_n_clusters, h)
        return clu_n_clusters
    elif isinstance(clu_n_clusters, tuple):
        try:
            clu_n_clusters = list(range(*clu_n_clusters))
        except:
            raise ValueError("Incorrect tuple specified for the number of "
                             "clusters.")
        if len(clu_n_clusters) < 1:
            raise ValueError("Empty tuple encountered for the number of "
                             "clusters.")
    elif isinstance(clu_n_clusters, list):
        if len(clu_n_clusters) < 1:
            raise ValueError("Empty list encountered for the number of "
                             "clusters.")
    else:
        raise ValueError("Incorrect format for the number of clusters.")

    _validate_n_clusters(clu_n_clusters[0], h)
    _validate_n_clusters(clu_n_clusters[-1], h)

    return clu_n_clusters


def _validate_n_jobs(n_jobs):
    if isinstance(n_jobs, str):
        try:
            n_jobs = literal_eval(n_jobs)
        except:
            raise ValueError("Incorrect number of jobs specified.")

    if isinstance(n_jobs, float):
        n_jobs = int(n_jobs)

    if isinstance(n_jobs, int):
        if n_jobs < -1:
            raise ValueError("Incorrect number of jobs specified.")
        elif n_jobs > 8:
            raise ValueError("Number of jobs is too high.")
        elif n_jobs == 0:
            raise ValueError("Number of jobs is 0.")
    elif n_jobs is not None:
        raise ValueError("Incorrect number of jobs specified.")

    return n_jobs


def _validate_mark_alpha(mark_alpha):
    if isinstance(mark_alpha, str):
        try:
            mark_alpha = literal_eval(mark_alpha)
        except:
            raise ValueError("Incorrect significance alpha set.")

    if isinstance(mark_alpha, int):
        mark_alpha = float(mark_alpha)

    if isinstance(mark_alpha, float):
        if mark_alpha > 0.15 or mark_alpha < 0.001:
            raise ValueError("Significance alpha needs to be in the interval "
                             "(0.001, 0.15)")
    else:
        raise ValueError("Incorrect significance alpha set.")

    return mark_alpha


def _validate_mark_markers_n(mark_markers_n, h):
    if isinstance(mark_markers_n, str):
        try:
            mark_markers_n = literal_eval(mark_markers_n)
        except:
            raise ValueError("Incorrect number of markers set.")

    if isinstance(mark_markers_n, float):
        mark_markers_n = int(mark_markers_n)

    if isinstance(mark_markers_n, int):
        if mark_markers_n < 1:
            raise ValueError("Number of markers needs to be greater than 1.")
        if mark_markers_n > h:
            raise ValueError("Number of markers needs to be less "
                             "than number of genes.")
    else:
        raise ValueError("Incorrect number of markers set.")

    return mark_markers_n


def _validate_mark_correction(mark_correction):
    corrections_methods = ["bonferroni", "sidak", "holm-sidak", "holm",
                           "simes-hochberg", "hommel", "fdr_bh", "fdr_by",
                           "fdr_tsbh", "fdr_tsbky"]

    if isinstance(mark_correction, str):
        if mark_correction not in corrections_methods:
            raise ValueError("Incorrect correction method.")
    else:
        raise ValueError("Incorrect correction method.")

    return mark_correction


def _validate_con_convention(con_convention):
    conventions = ["id-to-name", "name-to-id"]

    if isinstance(con_convention, str):
        if con_convention not in conventions:
            raise ValueError("Incorrect convention specified.")
    else:
        raise ValueError("Incorrect convention specified.")

    return con_convention
