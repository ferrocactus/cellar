from acip.unit import Unit

import numpy as np
from acip.k_medoids import KMedoids
from sklearn.cluster import KMeans
from sklearn.cluster import SpectralClustering

class Clu(Unit):
    def __init__(self, verbose=False, **args):
        """
        Base class for Clustering methods.

        Args:
            verbose (bool): Printing flag.
            **args: Argument list.
        """
        super().__init__(verbose, **args)
        self._labels = None
        self._n_clusters = None

    def get(self, x, eval_obj=None):
        """
        Args:
            x (np.ndarray): Data in matrix (n x d) form.
            eval_obj (Sco): Sco object to use for evaluating clusters.
                            Must be set if args['n_clusters'] is a tuple.
        Returns:
            clusters (np.ndarray): The labels for each x.
        """
        if isinstance(self.args['n_clusters'], int): # single cluster num
            self._labels = self.fit_predict(self._obj(**self.args), x)
            self._n_clusters = self.args['n_clusters']
        elif isinstance(self.args['n_clusters'], tuple): # range of cluster nums
            k_list = range(*self.args['n_clusters'])
            if len(k_list) < 1:
                raise ValueError()
            self.s_list = [0] * len(k_list)
            best_score = -np.Inf
            for i, k in enumerate(k_list): # Iterate over k
                self.args['n_clusters'] = k
                labels = self.fit_predict(self._obj(**self.args), x)
                score = eval_obj.score(x, labels)
                s_list[i] = score
                if best_score < score: # Update if best score found
                    best_score, self._labels, self._n_clusters = score, labels, k
        else:
            raise ValueError("Incorrect number of clusters used.")
        return self._labels

    def fit_predict(self, obj, x):
        """
        Wrapper around fit predict if the corresponding clustering
        method uses a different interface.

        Args:
            obj (Clu): Clu object.
            x (np.ndarray): Data to use for clustering.
        """
        return obj.fit_predict(x) # By default for most methods

class Clu_KMedoids(Clu):
    def __init__(self, verbose=False, **args):
        super().__init__(verbose, **args)
        self._obj = KMedoids

class Clu_KMeans(Clu):
    def __init__(self, verbose=False, **args):
        super().__init__(verbose, **args)
        self._obj = KMeans


class Clu_SpectralClustering(Clu):
    def __init__(self, verbose=False, **args):
        super().__init__(verbose, **args)
        if affinity not in self.args:
            self.args['affinity'] = 'nearest_neighbors'
        self._obj = SpectralClustering