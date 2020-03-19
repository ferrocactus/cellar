import numpy as np


def _effective_n_jobs(n_jobs):
    if n_jobs is None or n_jobs == 1:
        return 1
    else:
        return n_jobs


def _effective_n_clusters(n_clusters):
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