from acip.unit import Unit

from abc import abstractmethod
from sklearn.decomposition import PCA
from kneed import KneeLocator

PCA_EVR_MAX_N = 100

class Dim(Unit):
    def __init__(self, verbose=False, **args):
        """
        Base class for Dimensionality Reduction methods.

        Args:
            verbose (bool): Printing flag.
            **args: Argument list.
        """
        super().__init__(verbose, **args)

    @abstractmethod
    def get(self, x):
        """
        Returns the embedding of x.

        Args:
            x (np.ndarray): Data in matrix (n x d) form.
        Returns:
            emb (np.ndarray): The embedding of x.
        """
        pass

class Dim_PCA(Dim):
    def __init__(self, verbose=False, **args):
        super().__init__(verbose, **args)
        if 'n_components' not in args:
            raise ValueError("n_components not provided.")

    def get(self, x):
        if self.args['n_components'] == 'knee':
            self.args['n_components'] = min(PCA_EVR_MAX_N, x.shape[1])
            self.pca = PCA(**self.args)
            self.pca.fit(x)
            # Construct axis for KneeLocator
            x_axis = list(range(1, self.pca.n_components_ + 1))
            y_axis = self.pca.explained_variance_ratio_
            # Find knee
            self.args['n_components'] = self.knee = KneeLocator(
                x_axis,
                y_axis,
                curve='convex', # approximately
                direction='decreasing' # sklearn PCA eigenvalues are sorted
            ).knee
            self.args['n_components'] = max(self.knee, 2)
            self.vprint("Knee found at {0} components. Using n={1}.".format(
                self.knee, self.args['n_components'])
            )
            self.pca = PCA(**self.args)
            self.args['n_components'] = 'knee' # Reset to original arg
        else:
            self.pca = PCA(**self.args)
        return self.pca.fit_transform(x)