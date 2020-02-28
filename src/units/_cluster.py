from ._unit import Unit
from ._k_medoids import KMedoids

from abc import abstractmethod
import numpy as np
from sklearn.cluster import KMeans
from sklearn.cluster import SpectralClustering


class Cluster(Unit):
    """
    Base class for Clustering methods. A child class needs to instantiate
    a self._obj object with a fit_predict member function.
    """
    @abstractmethod
    def __init__(self, verbose=False, **kwargs):
        """
        Args:
            verbose (bool): Printing flag.
            **kwargs: Argument list.
        Raises:
            ValueError: If number of clusters to use is not provided.
        """
        if 'n_clusters' not in kwargs:
            raise ValueError("n_clusters not provided.")
        super().__init__(verbose, **kwargs)
        self.name = 'Clu'
        self._labels = None
        self._n_clusters = None

    def get(self, x, eval_obj=None):
        """
        Runs the clustering algorithm and returns labels for each point.
        If is possible to define a single k to use for clustering or pass
        a range tuple instead. In the latter case, an evaluation Eval object
        needs to be provided to choose the best k.

        WARNING! This method assumes a higher score to be a better one.

        Args:
            x (np.ndarray): Data in matrix (n x d) form.
            eval_obj (Sco): Sco object to use for evaluating clusters.
                            Must be set if kwargs['n_clusters'] is a tuple.
        Returns:
            clusters (np.ndarray): The labels for each x (must be integers).
        Raise:
            ValueError: If a tuple passed for k and Eval object not provided.
            ValueError: If range of k passed is invalid or empty.
            ValueError: If anything other than scalar or tuple is passed for k.
        """
        self._labels = None
        self._n_clusters = None
        if isinstance(self.kwargs['n_clusters'], int): # single cluster num
            self._labels = self.fit_predict(self._obj(**self.kwargs), x)
            self._n_clusters = self.kwargs['n_clusters']
            self.score_list = np.array([eval_obj.get(x, self._labels)])
        elif isinstance(self.kwargs['n_clusters'], tuple): # range of clusters
            temp_kwargs = self.kwargs.copy() # avoid editting the original dict
            k_list = range(*temp_kwargs['n_clusters'])

            if eval_obj is None: # Need evaluation method if using range
                raise ValueError("Evaluation object not provided.")
            if len(k_list) < 1:
                raise ValueError("Invalid k list encountered in clustering.")

            self.score_list = np.zeros(len(k_list))
            best_score = -np.Inf
            for i, k in enumerate(k_list): # Iterate over k
                temp_kwargs['n_clusters'] = k
                labels = self.fit_predict(self._obj(**temp_kwargs), x)
                score = eval_obj.get(x, labels)
                self.score_list[i] = score
                if best_score < score: # Update if best score found
                    best_score, self._labels, self._n_clusters = score,labels,k
                self.vprint("Finished clustering for k={0}. Score={1:.2f}.".format(
                    k, score
                ))

            self.vprint("Best score achieved for k={0} at {1:.2f}.".format(
                self._n_clusters, best_score
            ))
        else:
            raise ValueError("Incorrect number of clusters used.")
        return self._labels

    def fit_predict(self, obj, x):
        """
        Wrapper around fit predict if the corresponding clustering
        method uses a different interface. Reimplement in children
        classes if the fit predict method is different.

        Args:
            obj (Cluster): Cluster object.
            x (np.ndarray): Data to use for clustering.
        """
        return obj.fit_predict(x) # By default for most methods


class Clu_KMedoids(Cluster):
    def __init__(self, verbose=False, **kwargs):
        super().__init__(verbose, **kwargs)
        self._obj = KMedoids


class Clu_KMeans(Cluster):
    def __init__(self, verbose=False, **kwargs):
        super().__init__(verbose, **kwargs)
        self._obj = KMeans


class Clu_SpectralClustering(Cluster):
    def __init__(self, verbose=False, **kwargs):
        super().__init__(verbose, **kwargs)
        if affinity not in self.kwargs: # For consistency.
            self.kwargs['affinity'] = 'nearest_neighbors'
        self._obj = SpectralClustering