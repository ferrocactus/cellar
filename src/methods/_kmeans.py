"""K-means clustering"""

# Authors: Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Thomas Rueckstiess <ruecksti@in.tum.de>
#          James Bergstra <james.bergstra@umontreal.ca>
#          Jan Schlueter <scikit-learn@jan-schlueter.de>
#          Nelle Varoquaux
#          Peter Prettenhofer <peter.prettenhofer@gmail.com>
#          Olivier Grisel <olivier.grisel@ensta.org>
#          Mathieu Blondel <mathieu@mblondel.org>
#          Robert Layton <robertlayton@gmail.com>
# License: BSD 3 clause

import warnings

import numpy as np
import scipy.sparse as sp
from joblib import Parallel, delayed, effective_n_jobs

from sklearn.base import BaseEstimator, ClusterMixin, TransformerMixin
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.metrics.pairwise import pairwise_distances_argmin_min
from sklearn.utils.extmath import row_norms, squared_norm, stable_cumsum
from sklearn.utils.sparsefuncs_fast import assign_rows_csr
from sklearn.utils.sparsefuncs import mean_variance_axis
from sklearn.utils.validation import _num_samples
from sklearn.utils import check_array
from sklearn.utils import gen_batches
from sklearn.utils import check_random_state
from sklearn.utils.validation import check_is_fitted, _check_sample_weight
from sklearn.utils.validation import FLOAT_DTYPES
from sklearn.exceptions import ConvergenceWarning
from . import _k_means_fast as _k_means
from ._k_means_elkan import k_means_elkan


###############################################################################
# Initialization heuristic


def _k_init(X, n_clusters, x_squared_norms, random_state, n_local_trials=None):
    """Init n_clusters seeds according to k-means++

    Parameters
    ----------
    X : array or sparse matrix, shape (n_samples, n_features)
        The data to pick seeds for. To avoid memory copy, the input data
        should be double precision (dtype=np.float64).

    n_clusters : integer
        The number of seeds to choose

    x_squared_norms : array, shape (n_samples,)
        Squared Euclidean norm of each data point.

    random_state : int, RandomState instance
        The generator used to initialize the centers. Use an int to make the
        randomness deterministic.
        See :term:`Glossary <random_state>`.

    n_local_trials : integer, optional
        The number of seeding trials for each center (except the first),
        of which the one reducing inertia the most is greedily chosen.
        Set to None to make the number of trials depend logarithmically
        on the number of seeds (2+log(k)); this is the default.

    Notes
    -----
    Selects initial cluster centers for k-mean clustering in a smart way
    to speed up convergence. see: Arthur, D. and Vassilvitskii, S.
    "k-means++: the advantages of careful seeding". ACM-SIAM symposium
    on Discrete algorithms. 2007

    Version ported from http://www.stanford.edu/~darthur/kMeansppTest.zip,
    which is the implementation used in the aforementioned paper.
    """
    n_samples, n_features = X.shape

    centers = np.empty((n_clusters, n_features), dtype=X.dtype)

    assert x_squared_norms is not None, 'x_squared_norms None in _k_init'

    # Set the number of local seeding trials if none is given
    if n_local_trials is None:
        # This is what Arthur/Vassilvitskii tried, but did not report
        # specific results for other than mentioning in the conclusion
        # that it helped.
        n_local_trials = 2 + int(np.log(n_clusters))

    # Pick first center randomly
    center_id = random_state.randint(n_samples)
    if sp.issparse(X):
        centers[0] = X[center_id].toarray()
    else:
        centers[0] = X[center_id]

    # Initialize list of closest distances and calculate current potential
    closest_dist_sq = euclidean_distances(
        centers[0, np.newaxis], X, Y_norm_squared=x_squared_norms,
        squared=True)
    current_pot = closest_dist_sq.sum()

    # Pick the remaining n_clusters-1 points
    for c in range(1, n_clusters):
        # Choose center candidates by sampling with probability proportional
        # to the squared distance to the closest existing center
        rand_vals = random_state.random_sample(n_local_trials) * current_pot
        candidate_ids = np.searchsorted(stable_cumsum(closest_dist_sq),
                                        rand_vals)
        # XXX: numerical imprecision can result in a candidate_id out of range
        np.clip(candidate_ids, None, closest_dist_sq.size - 1,
                out=candidate_ids)

        # Compute distances to center candidates
        distance_to_candidates = euclidean_distances(
            X[candidate_ids], X, Y_norm_squared=x_squared_norms, squared=True)

        # update closest distances squared and potential for each candidate
        np.minimum(closest_dist_sq, distance_to_candidates,
                   out=distance_to_candidates)
        candidates_pot = distance_to_candidates.sum(axis=1)

        # Decide which candidate is the best
        best_candidate = np.argmin(candidates_pot)
        current_pot = candidates_pot[best_candidate]
        closest_dist_sq = distance_to_candidates[best_candidate]
        best_candidate = candidate_ids[best_candidate]

        # Permanently add best center candidate found in local tries
        if sp.issparse(X):
            centers[c] = X[best_candidate].toarray()
        else:
            centers[c] = X[best_candidate]

    return centers


###############################################################################
# K-means batch estimation by EM (expectation maximization)

def _validate_center_shape(X, n_centers, centers):
    """Check if centers is compatible with X and n_centers"""
    if len(centers) != n_centers:
        raise ValueError('The shape of the initial centers (%s) '
                         'does not match the number of clusters %i'
                         % (centers.shape, n_centers))
    if centers.shape[1] != X.shape[1]:
        raise ValueError(
            "The number of features of the initial centers %s "
            "does not match the number of features of the data %s."
            % (centers.shape[1], X.shape[1]))


def _tolerance(X, tol):
    """Return a tolerance which is independent of the dataset"""
    if sp.issparse(X):
        variances = mean_variance_axis(X, axis=0)[1]
    else:
        variances = np.var(X, axis=0)
    return np.mean(variances) * tol


def _check_normalize_sample_weight(sample_weight, X):
    """Set sample_weight if None, and check for correct dtype"""

    sample_weight_was_none = sample_weight is None

    sample_weight = _check_sample_weight(sample_weight, X, dtype=X.dtype)
    if not sample_weight_was_none:
        # normalize the weights to sum up to n_samples
        # an array of 1 (i.e. samples_weight is None) is already normalized
        n_samples = len(sample_weight)
        scale = n_samples / sample_weight.sum()
        sample_weight *= scale
    return sample_weight


def k_means(X, n_clusters, sample_weight=None, init='k-means++',
            precompute_distances='auto', n_init=10, max_iter=300,
            verbose=False, tol=1e-4, random_state=None, copy_x=True,
            n_jobs=None, algorithm="auto", return_n_iter=False, mask=None,
            fixed_labels=None):
    """K-means clustering algorithm.

    Read more in the :ref:`User Guide <k_means>`.

    Parameters
    ----------
    X : array-like or sparse matrix, shape (n_samples, n_features)
        The observations to cluster. It must be noted that the data
        will be converted to C ordering, which will cause a memory copy
        if the given data is not C-contiguous.

    n_clusters : int
        The number of clusters to form as well as the number of
        centroids to generate.

    sample_weight : array-like, shape (n_samples,), optional
        The weights for each observation in X. If None, all observations
        are assigned equal weight (default: None)

    init : {'k-means++', 'random', or ndarray, or a callable}, optional
        Method for initialization, default to 'k-means++':

        'k-means++' : selects initial cluster centers for k-mean
        clustering in a smart way to speed up convergence. See section
        Notes in k_init for more details.

        'random': choose k observations (rows) at random from data for
        the initial centroids.

        If an ndarray is passed, it should be of shape (n_clusters, n_features)
        and gives the initial centers.

        If a callable is passed, it should take arguments X, k and
        and a random state and return an initialization.

    precompute_distances : {'auto', True, False}
        Precompute distances (faster but takes more memory).

        'auto' : do not precompute distances if n_samples * n_clusters > 12
        million. This corresponds to about 100MB overhead per job using
        double precision.

        True : always precompute distances

        False : never precompute distances

    n_init : int, optional, default: 10
        Number of time the k-means algorithm will be run with different
        centroid seeds. The final results will be the best output of
        n_init consecutive runs in terms of inertia.

    max_iter : int, optional, default 300
        Maximum number of iterations of the k-means algorithm to run.

    verbose : boolean, optional
        Verbosity mode.

    tol : float, optional
        The relative increment in the results before declaring convergence.

    random_state : int, RandomState instance or None (default)
        Determines random number generation for centroid initialization. Use
        an int to make the randomness deterministic.
        See :term:`Glossary <random_state>`.

    copy_x : bool, optional
        When pre-computing distances it is more numerically accurate to center
        the data first.  If copy_x is True (default), then the original data is
        not modified, ensuring X is C-contiguous.  If False, the original data
        is modified, and put back before the function returns, but small
        numerical differences may be introduced by subtracting and then adding
        the data mean, in this case it will also not ensure that data is
        C-contiguous which may cause a significant slowdown.

    n_jobs : int or None, optional (default=None)
        The number of jobs to use for the computation. This works by computing
        each of the n_init runs in parallel.

        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    algorithm : "auto", "full" or "elkan", default="auto"
        K-means algorithm to use. The classical EM-style algorithm is "full".
        The "elkan" variation is more efficient by using the triangle
        inequality, but currently doesn't support sparse data. "auto" chooses
        "elkan" for dense data and "full" for sparse data.

    return_n_iter : bool, optional
        Whether or not to return the number of iterations.

    Returns
    -------
    centroid : float ndarray with shape (k, n_features)
        Centroids found at the last iteration of k-means.

    label : integer ndarray with shape (n_samples,)
        label[i] is the code or index of the centroid the
        i'th observation is closest to.

    inertia : float
        The final value of the inertia criterion (sum of squared distances to
        the closest centroid for all observations in the training set).

    best_n_iter : int
        Number of iterations corresponding to the best results.
        Returned only if `return_n_iter` is set to True.
    """

    est = KMeans(
        n_clusters=n_clusters, init=init, n_init=n_init, max_iter=max_iter,
        verbose=verbose, precompute_distances=precompute_distances, tol=tol,
        random_state=random_state, copy_x=copy_x, n_jobs=n_jobs,
        algorithm=algorithm, mask=mask, fixed_labels=fixed_labels
    ).fit(X, sample_weight=sample_weight)
    if return_n_iter:
        return est.cluster_centers_, est.labels_, est.inertia_, est.n_iter_
    else:
        return est.cluster_centers_, est.labels_, est.inertia_


def _kmeans_single_elkan(X, sample_weight, n_clusters, max_iter=300,
                         init='k-means++', verbose=False, x_squared_norms=None,
                         random_state=None, tol=1e-4,
                         precompute_distances=True, mask=None,
                         fixed_labels=None):
    if sp.issparse(X):
        raise TypeError("algorithm='elkan' not supported for sparse input X")
    random_state = check_random_state(random_state)
    if x_squared_norms is None:
        x_squared_norms = row_norms(X, squared=True)
    # init
    centers = _init_centroids(X, n_clusters, init, random_state=random_state,
                              x_squared_norms=x_squared_norms)
    centers = np.ascontiguousarray(centers)
    if verbose:
        print('Initialization complete')

    checked_sample_weight = _check_normalize_sample_weight(sample_weight, X)
    centers, labels, n_iter = k_means_elkan(X, checked_sample_weight,
                                            n_clusters, centers, mask,
                                            fixed_labels, tol=tol,
                                            max_iter=max_iter, verbose=verbose)
    if sample_weight is None:
        inertia = np.sum((X - centers[labels]) ** 2, dtype=np.float64)
    else:
        sq_distances = np.sum((X - centers[labels]) ** 2, axis=1,
                              dtype=np.float64) * checked_sample_weight
        inertia = np.sum(sq_distances, dtype=np.float64)
    return labels, inertia, centers, n_iter


def _kmeans_single_lloyd(X, sample_weight, n_clusters, max_iter=300,
                         init='k-means++', verbose=False, x_squared_norms=None,
                         random_state=None, tol=1e-4,
                         precompute_distances=True, mask=None,
                         fixed_labels=None):
    """A single run of k-means, assumes preparation completed prior.

    Parameters
    ----------
    X : array-like of floats, shape (n_samples, n_features)
        The observations to cluster.

    n_clusters : int
        The number of clusters to form as well as the number of
        centroids to generate.

    sample_weight : array-like, shape (n_samples,)
        The weights for each observation in X.

    max_iter : int, optional, default 300
        Maximum number of iterations of the k-means algorithm to run.

    init : {'k-means++', 'random', or ndarray, or a callable}, optional
        Method for initialization, default to 'k-means++':

        'k-means++' : selects initial cluster centers for k-mean
        clustering in a smart way to speed up convergence. See section
        Notes in k_init for more details.

        'random': choose k observations (rows) at random from data for
        the initial centroids.

        If an ndarray is passed, it should be of shape (k, p) and gives
        the initial centers.

        If a callable is passed, it should take arguments X, k and
        and a random state and return an initialization.

    tol : float, optional
        The relative increment in the results before declaring convergence.

    verbose : boolean, optional
        Verbosity mode

    x_squared_norms : array
        Precomputed x_squared_norms.

    precompute_distances : boolean, default: True
        Precompute distances (faster but takes more memory).

    random_state : int, RandomState instance or None (default)
        Determines random number generation for centroid initialization. Use
        an int to make the randomness deterministic.
        See :term:`Glossary <random_state>`.

    mask : None or ndarray
        If set to 0 do not change the label for that point.

    Returns
    -------
    centroid : float ndarray with shape (k, n_features)
        Centroids found at the last iteration of k-means.

    label : integer ndarray with shape (n_samples,)
        label[i] is the code or index of the centroid the
        i'th observation is closest to.

    inertia : float
        The final value of the inertia criterion (sum of squared distances to
        the closest centroid for all observations in the training set).

    n_iter : int
        Number of iterations run.
    """
    random_state = check_random_state(random_state)

    sample_weight = _check_normalize_sample_weight(sample_weight, X)

    best_labels, best_inertia, best_centers = None, None, None
    # init
    centers = _init_centroids(X, n_clusters, init, random_state=random_state,
                              x_squared_norms=x_squared_norms)
    if verbose:
        print("Initialization complete")

    # Allocate memory to store the distances for each sample to its
    # closer center for reallocation in case of ties
    distances = np.zeros(shape=(X.shape[0],), dtype=X.dtype)

    # iterations
    for i in range(max_iter):
        centers_old = centers.copy()
        # labels assignment is also called the E-step of EM
        labels, inertia = \
            _labels_inertia(X, sample_weight, x_squared_norms, centers,
                            precompute_distances=precompute_distances,
                            distances=distances,
                            mask=mask, fixed_labels=fixed_labels)

        # computation of the means is also called the M-step of EM
        if sp.issparse(X):
            centers = _k_means._centers_sparse(X, sample_weight, labels,
                                               n_clusters, distances)
        else:
            centers = _k_means._centers_dense(X, sample_weight, labels,
                                              n_clusters, distances)

        if verbose:
            print("Iteration %2d, inertia %.3f" % (i, inertia))

        if best_inertia is None or inertia < best_inertia:
            best_labels = labels.copy()
            best_centers = centers.copy()
            best_inertia = inertia

        center_shift_total = squared_norm(centers_old - centers)
        if center_shift_total <= tol:
            if verbose:
                print("Converged at iteration %d: "
                      "center shift %e within tolerance %e"
                      % (i, center_shift_total, tol))
            break

    if center_shift_total > 0:
        # rerun E-step in case of non-convergence so that predicted labels
        # match cluster centers
        best_labels, best_inertia = \
            _labels_inertia(X, sample_weight, x_squared_norms, best_centers,
                            precompute_distances=precompute_distances,
                            distances=distances, mask=mask,
                            fixed_labels=fixed_labels)

    return best_labels, best_inertia, best_centers, i + 1


def _labels_inertia_precompute_dense(X, sample_weight, x_squared_norms,
                                     centers, distances):
    """Compute labels and inertia using a full distance matrix.

    This will overwrite the 'distances' array in-place.

    Parameters
    ----------
    X : numpy array, shape (n_sample, n_features)
        Input data.

    sample_weight : array-like, shape (n_samples,)
        The weights for each observation in X.

    x_squared_norms : numpy array, shape (n_samples,)
        Precomputed squared norms of X.

    centers : numpy array, shape (n_clusters, n_features)
        Cluster centers which data is assigned to.

    distances : numpy array, shape (n_samples,)
        Pre-allocated array in which distances are stored.

    Returns
    -------
    labels : numpy array, dtype=np.int, shape (n_samples,)
        Indices of clusters that samples are assigned to.

    inertia : float
        Sum of squared distances of samples to their closest cluster center.

    """
    n_samples = X.shape[0]

    # Breakup nearest neighbor distance computation into batches to prevent
    # memory blowup in the case of a large number of samples and clusters.
    # TODO: Once PR #7383 is merged use check_inputs=False in metric_kwargs.
    labels, mindist = pairwise_distances_argmin_min(
        X=X, Y=centers, metric='euclidean', metric_kwargs={'squared': True})
    # cython k-means code assumes int32 inputs
    labels = labels.astype(np.int32, copy=False)
    if n_samples == distances.shape[0]:
        # distances will be changed in-place
        distances[:] = mindist
    inertia = (mindist * sample_weight).sum()
    return labels, inertia


def _labels_inertia(X, sample_weight, x_squared_norms, centers,
                    precompute_distances=True, distances=None,
                    mask=None, fixed_labels=None):
    """E step of the K-means EM algorithm.

    Compute the labels and the inertia of the given samples and centers.
    This will compute the distances in-place.

    Parameters
    ----------
    X : float64 array-like or CSR sparse matrix, shape (n_samples, n_features)
        The input samples to assign to the labels.

    sample_weight : array-like, shape (n_samples,)
        The weights for each observation in X.

    x_squared_norms : array, shape (n_samples,)
        Precomputed squared euclidean norm of each data point, to speed up
        computations.

    centers : float array, shape (k, n_features)
        The cluster centers.

    precompute_distances : boolean, default: True
        Precompute distances (faster but takes more memory).

    distances : float array, shape (n_samples,)
        Pre-allocated array to be filled in with each sample's distance
        to the closest center.

    Returns
    -------
    labels : int array of shape(n)
        The resulting assignment

    inertia : float
        Sum of squared distances of samples to their closest cluster center.
    """
    n_samples = X.shape[0]
    sample_weight = _check_normalize_sample_weight(sample_weight, X)
    # set the default value of centers to -1 to be able to detect any anomaly
    # easily
    labels = np.full(n_samples, -1, np.int32)
    if distances is None:
        distances = np.zeros(shape=(0,), dtype=X.dtype)
    # distances will be changed in-place
    if sp.issparse(X):
        inertia = _k_means._assign_labels_csr(
            X, sample_weight, x_squared_norms, centers, labels,
            distances=distances, mask=mask, fixed_labels=fixed_labels)
    else:
        if precompute_distances:
            return _labels_inertia_precompute_dense(X, sample_weight,
                                                    x_squared_norms, centers,
                                                    distances)
        inertia = _k_means._assign_labels_array(
            X, sample_weight, x_squared_norms, centers, labels,
            distances=distances, mask=mask, fixed_labels=fixed_labels)
    return labels, inertia


def _init_centroids(X, k, init, random_state=None, x_squared_norms=None,
                    init_size=None):
    """Compute the initial centroids

    Parameters
    ----------

    X : array, shape (n_samples, n_features)

    k : int
        number of centroids

    init : {'k-means++', 'random' or ndarray or callable} optional
        Method for initialization

    random_state : int, RandomState instance or None (default)
        Determines random number generation for centroid initialization. Use
        an int to make the randomness deterministic.
        See :term:`Glossary <random_state>`.

    x_squared_norms : array, shape (n_samples,), optional
        Squared euclidean norm of each data point. Pass it if you have it at
        hands already to avoid it being recomputed here. Default: None

    init_size : int, optional
        Number of samples to randomly sample for speeding up the
        initialization (sometimes at the expense of accuracy): the
        only algorithm is initialized by running a batch KMeans on a
        random subset of the data. This needs to be larger than k.

    Returns
    -------
    centers : array, shape(k, n_features)
    """
    random_state = check_random_state(random_state)
    n_samples = X.shape[0]

    if x_squared_norms is None:
        x_squared_norms = row_norms(X, squared=True)

    if init_size is not None and init_size < n_samples:
        if init_size < k:
            warnings.warn(
                "init_size=%d should be larger than k=%d. "
                "Setting it to 3*k" % (init_size, k),
                RuntimeWarning, stacklevel=2)
            init_size = 3 * k
        init_indices = random_state.randint(0, n_samples, init_size)
        X = X[init_indices]
        x_squared_norms = x_squared_norms[init_indices]
        n_samples = X.shape[0]
    elif n_samples < k:
        raise ValueError(
            "n_samples=%d should be larger than k=%d" % (n_samples, k))

    if isinstance(init, str) and init == 'k-means++':
        centers = _k_init(X, k, random_state=random_state,
                          x_squared_norms=x_squared_norms)
    elif isinstance(init, str) and init == 'random':
        seeds = random_state.permutation(n_samples)[:k]
        centers = X[seeds]
    elif hasattr(init, '__array__'):
        # ensure that the centers have the same dtype as X
        # this is a requirement of fused types of cython
        centers = np.array(init, dtype=X.dtype)
    elif callable(init):
        centers = init(X, k, random_state=random_state)
        centers = np.asarray(centers, dtype=X.dtype)
    else:
        raise ValueError("the init parameter for the k-means should "
                         "be 'k-means++' or 'random' or an ndarray, "
                         "'%s' (type '%s') was passed." % (init, type(init)))

    if sp.issparse(centers):
        centers = centers.toarray()

    _validate_center_shape(X, k, centers)
    return centers


class KMeans(TransformerMixin, ClusterMixin, BaseEstimator):
    """K-Means clustering.

    Read more in the :ref:`User Guide <k_means>`.

    Parameters
    ----------

    n_clusters : int, default=8
        The number of clusters to form as well as the number of
        centroids to generate.

    init : {'k-means++', 'random'} or ndarray of shape \
            (n_clusters, n_features), default='k-means++'
        Method for initialization, defaults to 'k-means++':

        'k-means++' : selects initial cluster centers for k-mean
        clustering in a smart way to speed up convergence. See section
        Notes in k_init for more details.

        'random': choose k observations (rows) at random from data for
        the initial centroids.

        If an ndarray is passed, it should be of shape (n_clusters, n_features)
        and gives the initial centers.

    n_init : int, default=10
        Number of time the k-means algorithm will be run with different
        centroid seeds. The final results will be the best output of
        n_init consecutive runs in terms of inertia.

    max_iter : int, default=300
        Maximum number of iterations of the k-means algorithm for a
        single run.

    tol : float, default=1e-4
        Relative tolerance with regards to inertia to declare convergence.

    precompute_distances : 'auto' or bool, default='auto'
        Precompute distances (faster but takes more memory).

        'auto' : do not precompute distances if n_samples * n_clusters > 12
        million. This corresponds to about 100MB overhead per job using
        double precision.

        True : always precompute distances.

        False : never precompute distances.

    verbose : int, default=0
        Verbosity mode.

    random_state : int, RandomState instance, default=None
        Determines random number generation for centroid initialization. Use
        an int to make the randomness deterministic.
        See :term:`Glossary <random_state>`.

    copy_x : bool, default=True
        When pre-computing distances it is more numerically accurate to center
        the data first.  If copy_x is True (default), then the original data is
        not modified, ensuring X is C-contiguous.  If False, the original data
        is modified, and put back before the function returns, but small
        numerical differences may be introduced by subtracting and then adding
        the data mean, in this case it will also not ensure that data is
        C-contiguous which may cause a significant slowdown.

    n_jobs : int, default=None
        The number of jobs to use for the computation. This works by computing
        each of the n_init runs in parallel.

        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    algorithm : {"auto", "full", "elkan"}, default="auto"
        K-means algorithm to use. The classical EM-style algorithm is "full".
        The "elkan" variation is more efficient by using the triangle
        inequality, but currently doesn't support sparse data. "auto" chooses
        "elkan" for dense data and "full" for sparse data.

    mask : None or ndarray of shape (n_samples,) Default None
        If set to 0, then that point's label cannot change.
        If set to 1, run the algorithm for that point normally.

    Attributes
    ----------
    cluster_centers_ : ndarray of shape (n_clusters, n_features)
        Coordinates of cluster centers. If the algorithm stops before fully
        converging (see ``tol`` and ``max_iter``), these will not be
        consistent with ``labels_``.

    labels_ : ndarray of shape (n_samples,)
        Labels of each point

    inertia_ : float
        Sum of squared distances of samples to their closest cluster center.

    n_iter_ : int
        Number of iterations run.

    See Also
    --------

    MiniBatchKMeans
        Alternative online implementation that does incremental updates
        of the centers positions using mini-batches.
        For large scale learning (say n_samples > 10k) MiniBatchKMeans is
        probably much faster than the default batch implementation.

    Notes
    -----
    The k-means problem is solved using either Lloyd's or Elkan's algorithm.

    The average complexity is given by O(k n T), were n is the number of
    samples and T is the number of iteration.

    The worst case complexity is given by O(n^(k+2/p)) with
    n = n_samples, p = n_features. (D. Arthur and S. Vassilvitskii,
    'How slow is the k-means method?' SoCG2006)

    In practice, the k-means algorithm is very fast (one of the fastest
    clustering algorithms available), but it falls in local minima. That's why
    it can be useful to restart it several times.

    If the algorithm stops before fully converging (because of ``tol`` or
    ``max_iter``), ``labels_`` and ``cluster_centers_`` will not be consistent,
    i.e. the ``cluster_centers_`` will not be the means of the points in each
    cluster. Also, the estimator will reassign ``labels_`` after the last
    iteration to make ``labels_`` consistent with ``predict`` on the training
    set.

    Examples
    --------

    >>> from sklearn.cluster import KMeans
    >>> import numpy as np
    >>> X = np.array([[1, 2], [1, 4], [1, 0],
    ...               [10, 2], [10, 4], [10, 0]])
    >>> kmeans = KMeans(n_clusters=2, random_state=0).fit(X)
    >>> kmeans.labels_
    array([1, 1, 1, 0, 0, 0], dtype=int32)
    >>> kmeans.predict([[0, 0], [12, 3]])
    array([1, 0], dtype=int32)
    >>> kmeans.cluster_centers_
    array([[10.,  2.],
           [ 1.,  2.]])
    """

    def __init__(self, n_clusters=8, init='k-means++', n_init=10,
                 max_iter=300, tol=1e-4, precompute_distances='auto',
                 verbose=0, random_state=None, copy_x=True,
                 n_jobs=None, algorithm='auto', mask=None,
                 fixed_labels=None):

        self.n_clusters = n_clusters
        self.init = init
        self.max_iter = max_iter
        self.tol = tol
        self.precompute_distances = precompute_distances
        self.n_init = n_init
        self.verbose = verbose
        self.random_state = random_state
        self.copy_x = copy_x
        self.n_jobs = n_jobs
        self.algorithm = algorithm
        self.mask = mask
        self.fixed_labels = fixed_labels

    def _check_test_data(self, X):
        X = check_array(X, accept_sparse='csr', dtype=FLOAT_DTYPES)
        n_samples, n_features = X.shape
        expected_n_features = self.cluster_centers_.shape[1]
        if not n_features == expected_n_features:
            raise ValueError("Incorrect number of features. "
                             "Got %d features, expected %d" % (
                                 n_features, expected_n_features))

        return X

    def fit(self, X, y=None, sample_weight=None):
        """Compute k-means clustering.

        Parameters
        ----------
        X : array-like or sparse matrix, shape=(n_samples, n_features)
            Training instances to cluster. It must be noted that the data
            will be converted to C ordering, which will cause a memory
            copy if the given data is not C-contiguous.

        y : Ignored
            Not used, present here for API consistency by convention.

        sample_weight : array-like, shape (n_samples,), optional
            The weights for each observation in X. If None, all observations
            are assigned equal weight (default: None).

        Returns
        -------
        self
            Fitted estimator.
        """
        random_state = check_random_state(self.random_state)

        n_init = self.n_init
        if n_init <= 0:
            raise ValueError("Invalid number of initializations."
                             " n_init=%d must be bigger than zero." % n_init)

        if self.max_iter <= 0:
            raise ValueError(
                'Number of iterations should be a positive number,'
                ' got %d instead' % self.max_iter
            )

        # avoid forcing order when copy_x=False
        order = "C" if self.copy_x else None
        X = check_array(X, accept_sparse='csr', dtype=[np.float64, np.float32],
                        order=order, copy=self.copy_x)
        # verify that the number of samples given is larger than k
        if _num_samples(X) < self.n_clusters:
            raise ValueError("n_samples=%d should be >= n_clusters=%d" % (
                _num_samples(X), self.n_clusters))

        tol = _tolerance(X, self.tol)

        # If the distances are precomputed every job will create a matrix of
        # shape (n_clusters, n_samples). To stop KMeans from eating up memory
        # we only activate this if the created matrix is guaranteed to be
        # under 100MB. 12 million entries consume a little under 100MB if they
        # are of type double.
        precompute_distances = self.precompute_distances
        if precompute_distances == 'auto':
            n_samples = X.shape[0]
            precompute_distances = (self.n_clusters * n_samples) < 12e6
        elif isinstance(precompute_distances, bool):
            pass
        else:
            raise ValueError(
                "precompute_distances should be 'auto' or True/False"
                ", but a value of %r was passed" %
                precompute_distances
            )

        # Validate init array
        init = self.init
        if hasattr(init, '__array__'):
            init = check_array(init, dtype=X.dtype.type, copy=True)
            _validate_center_shape(X, self.n_clusters, init)

            if n_init != 1:
                warnings.warn(
                    'Explicit initial center position passed: '
                    'performing only one init in k-means instead of n_init=%d'
                    % n_init, RuntimeWarning, stacklevel=2)
                n_init = 1

        # subtract of mean of x for more accurate distance computations
        if not sp.issparse(X):
            X_mean = X.mean(axis=0)
            # The copy was already done above
            X -= X_mean

            if hasattr(init, '__array__'):
                init -= X_mean

        # precompute squared norms of data points
        x_squared_norms = row_norms(X, squared=True)

        best_labels, best_inertia, best_centers = None, None, None
        algorithm = self.algorithm
        if self.n_clusters == 1:
            # elkan doesn't make sense for a single cluster, full will produce
            # the right result.
            algorithm = "full"
        if algorithm == "auto":
            algorithm = "full" if sp.issparse(X) else 'elkan'
        if algorithm == "full":
            kmeans_single = _kmeans_single_lloyd
        elif algorithm == "elkan":
            kmeans_single = _kmeans_single_elkan
        else:
            raise ValueError("Algorithm must be 'auto', 'full' or 'elkan', got"
                             " %s" % str(algorithm))

        seeds = random_state.randint(np.iinfo(np.int32).max, size=n_init)
        if effective_n_jobs(self.n_jobs) == 1:
            # For a single thread, less memory is needed if we just store one
            # set of the best results (as opposed to one set per run per
            # thread).
            for seed in seeds:
                # run a k-means once
                labels, inertia, centers, n_iter_ = kmeans_single(
                    X, sample_weight, self.n_clusters,
                    max_iter=self.max_iter, init=init, verbose=self.verbose,
                    precompute_distances=precompute_distances, tol=tol,
                    x_squared_norms=x_squared_norms, random_state=seed,
                    mask=self.mask, fixed_labels=self.fixed_labels)
                # determine if these results are the best so far
                if best_inertia is None or inertia < best_inertia:
                    best_labels = labels.copy()
                    best_centers = centers.copy()
                    best_inertia = inertia
                    best_n_iter = n_iter_
        else:
            # parallelisation of k-means runs
            results = Parallel(n_jobs=self.n_jobs, verbose=0)(
                delayed(kmeans_single)(
                    X, sample_weight, self.n_clusters,
                    max_iter=self.max_iter, init=init,
                    verbose=self.verbose, tol=tol,
                    precompute_distances=precompute_distances,
                    x_squared_norms=x_squared_norms,
                    # Change seed to ensure variety
                    random_state=seed,
                    mask=self.mask,
                    fixed_labels=self.fixed_labels
                )
                for seed in seeds)
            # Get results with the lowest inertia
            labels, inertia, centers, n_iters = zip(*results)
            best = np.argmin(inertia)
            best_labels = labels[best]
            best_inertia = inertia[best]
            best_centers = centers[best]
            best_n_iter = n_iters[best]

        if not sp.issparse(X):
            if not self.copy_x:
                X += X_mean
            best_centers += X_mean

        distinct_clusters = len(set(best_labels))
        if distinct_clusters < self.n_clusters:
            warnings.warn(
                "Number of distinct clusters ({}) found smaller than "
                "n_clusters ({}). Possibly due to duplicate points "
                "in X.".format(distinct_clusters, self.n_clusters),
                ConvergenceWarning, stacklevel=2
            )

        self.cluster_centers_ = best_centers
        self.labels_ = best_labels
        self.inertia_ = best_inertia
        self.n_iter_ = best_n_iter
        return self

    def fit_predict(self, X, y=None, sample_weight=None):
        """Compute cluster centers and predict cluster index for each sample.

        Convenience method; equivalent to calling fit(X) followed by
        predict(X).

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            New data to transform.

        y : Ignored
            Not used, present here for API consistency by convention.

        sample_weight : array-like, shape (n_samples,), optional
            The weights for each observation in X. If None, all observations
            are assigned equal weight (default: None).

        Returns
        -------
        labels : array, shape [n_samples,]
            Index of the cluster each sample belongs to.
        """
        return self.fit(X, sample_weight=sample_weight).labels_

    def fit_transform(self, X, y=None, sample_weight=None):
        """Compute clustering and transform X to cluster-distance space.

        Equivalent to fit(X).transform(X), but more efficiently implemented.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            New data to transform.

        y : Ignored
            Not used, present here for API consistency by convention.

        sample_weight : array-like, shape (n_samples,), optional
            The weights for each observation in X. If None, all observations
            are assigned equal weight (default: None).

        Returns
        -------
        X_new : array, shape [n_samples, k]
            X transformed in the new space.
        """
        # Currently, this just skips a copy of the data if it is not in
        # np.array or CSR format already.
        # XXX This skips _check_test_data, which may change the dtype;
        # we should refactor the input validation.
        return self.fit(X, sample_weight=sample_weight)._transform(X)

    def transform(self, X):
        """Transform X to a cluster-distance space.

        In the new space, each dimension is the distance to the cluster
        centers.  Note that even if X is sparse, the array returned by
        `transform` will typically be dense.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            New data to transform.

        Returns
        -------
        X_new : array, shape [n_samples, k]
            X transformed in the new space.
        """
        check_is_fitted(self)

        X = self._check_test_data(X)
        return self._transform(X)

    def _transform(self, X):
        """guts of transform method; no input validation"""
        return euclidean_distances(X, self.cluster_centers_)

    def predict(self, X, sample_weight=None):
        """Predict the closest cluster each sample in X belongs to.

        In the vector quantization literature, `cluster_centers_` is called
        the code book and each value returned by `predict` is the index of
        the closest code in the code book.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            New data to predict.

        sample_weight : array-like, shape (n_samples,), optional
            The weights for each observation in X. If None, all observations
            are assigned equal weight (default: None).

        Returns
        -------
        labels : array, shape [n_samples,]
            Index of the cluster each sample belongs to.
        """
        check_is_fitted(self)

        X = self._check_test_data(X)
        x_squared_norms = row_norms(X, squared=True)
        return _labels_inertia(X, sample_weight, x_squared_norms,
                               self.cluster_centers_, self.mask, self.fixed_labels)[0]

    def score(self, X, y=None, sample_weight=None):
        """Opposite of the value of X on the K-means objective.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            New data.

        y : Ignored
            Not used, present here for API consistency by convention.

        sample_weight : array-like, shape (n_samples,), optional
            The weights for each observation in X. If None, all observations
            are assigned equal weight (default: None).

        Returns
        -------
        score : float
            Opposite of the value of X on the K-means objective.
        """
        check_is_fitted(self)

        X = self._check_test_data(X)
        x_squared_norms = row_norms(X, squared=True)
        return -_labels_inertia(X, sample_weight, x_squared_norms,
                                self.cluster_centers_, self.mask,
                                self.fixed_labels)[1]
