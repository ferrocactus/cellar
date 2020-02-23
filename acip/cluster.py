from acip.unit import Unit

from abc import abstractmethod
import numpy as np
from acip.k_medoids import KMedoids
from sklearn.cluster import KMeans
from sklearn.cluster import SpectralClustering

class Cluster(Unit):
    @abstractmethod
    def __init__(self, verbose=False, **args):
        """
        Base class for Clustering methods. A child class needs to instantiate
        a self._obj object with a fit_predict member function.

        Args:
            verbose (bool): Printing flag.
            **args: Argument list.
        """
        super().__init__(verbose, **args)
        if 'n_clusters' not in args:
            raise ValueError("n_clusters not provided.")
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
                            Must be set if args['n_clusters'] is a tuple.
        Returns:
            clusters (np.ndarray): The labels for each x.
        """
        self._labels = None
        self._n_clusters = None
        if isinstance(self.args['n_clusters'], int): # single cluster num
            self._labels = self.fit_predict(self._obj(**self.args), x)
            self._n_clusters = self.args['n_clusters']
        elif isinstance(self.args['n_clusters'], tuple): # range of cluster nums
            temp = self.args['n_clusters'] # cache for later
            if eval_obj is None: # Need evaluation method if using range
                raise ValueError("eval_obj not provided.")
            k_list = range(*self.args['n_clusters'])
            if len(k_list) < 1:
                raise ValueError()
            self.score_list = [0] * len(k_list)
            best_score = -np.Inf
            for i, k in enumerate(k_list): # Iterate over k
                self.args['n_clusters'] = k
                labels = self.fit_predict(self._obj(**self.args), x)
                score = eval_obj.get(x, labels)
                self.score_list[i] = score
                if best_score < score: # Update if best score found
                    best_score, self._labels, self._n_clusters = score, labels, k
            self.vprint("Best score achieved for n={0} at {1:.2f}.".format(
                self._n_clusters, best_score
            ))
            self.args['n_clusters'] = temp # reset to original arg
        else:
            raise ValueError("Incorrect number of clusters used.")
        return self._labels

    def fit_predict(self, obj, x):
        """
        Wrapper around fit predict if the corresponding clustering
        method uses a different interface. Reimplement in children
        classes if the fit predict method is different.

        Args:
            obj (Clu): Clu object.
            x (np.ndarray): Data to use for clustering.
        """
        return obj.fit_predict(x) # By default for most methods

class Clu_KMedoids(Cluster):
    def __init__(self, verbose=False, **args):
        super().__init__(verbose, **args)
        self._obj = KMedoids

class Clu_KMeans(Cluster):
    def __init__(self, verbose=False, **args):
        super().__init__(verbose, **args)
        self._obj = KMeans

class Clu_SpectralClustering(Cluster):
    def __init__(self, verbose=False, **args):
        super().__init__(verbose, **args)
        if affinity not in self.args:
            self.args['affinity'] = 'nearest_neighbors'
        self._obj = SpectralClustering