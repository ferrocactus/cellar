from acip.unit import Unit

from abc import abstractmethod
from sklearn.metrics import silhouette_score
from sklearn.metrics import davies_bouldin_score

class Eval(Unit):
    def __init__(self, verbose=False, **args):
        """
        Base class for Evaluation methods.

        Args:
            verbose (bool): Printing flag.
            **args: Argument list.
        """
        super().__init__(verbose, **args)
        self._score = None

    @abstractmethod
    def get(self, x, labels):
        """
        Args:
            x (np.ndarray): Data in matrix (n x d) form.
            labels (np.ndarray): Labels corresponding to x
        Returns:
            score (int): The cluster score.
        """
        return self._score

class Eval_SilhouetteScore(Eval):
    def __init__(self, verbose=False, **args):
        super().__init__(verbose, **args)

    def get(self, x, labels):
        self._score = silhouette_score(x, labels, **self.args)
        self.vprint("Silhouette Score: {0:.2f}.".format(self._score))
        return self._score

class Eval_DaviesBouldinScore(Eval):
    def __init__(self, verbose=False, **args):
        super().__init__(verbose, **args)

    def get(self, x, labels):
        self._score = davies_bouldin_score(x, labels, **self.args)
        self.vprint("Davies Bouldin Score: {0:.2f}.".format(self._score))
        # Return negative the result, because the db score
        # assumes a better clustering if the score is lower
        return -self._score