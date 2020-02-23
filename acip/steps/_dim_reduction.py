from ._unit import Unit

from abc import abstractmethod
from sklearn.decomposition import PCA
from kneed import KneeLocator
from umap import UMAP
from sklearn.manifold import TSNE

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
            (np.ndarray): The embedding of x.
        """
        return x


class Dim_PCA(Dim):
    def __init__(self, verbose=False, **args):
        super().__init__(verbose, **args)
        if 'n_components' not in args:
            raise ValueError("n_components not provided.")

    def get(self, x):
        if self.args['n_components'] == 'knee':
            temp_args = self.args.copy()
            temp_args['n_components'] = min(PCA_EVR_MAX_N, x.shape[1])
            self.ev_pca = PCA(**temp_args)
            self.ev_pca.fit(x)
            # Construct axis for KneeLocator
            x_axis = list(range(1, self.ev_pca.n_components_ + 1))
            y_axis = self.ev_pca.explained_variance_ratio_
            # Find knee
            temp_args['n_components'] = self.knee = KneeLocator(
                x_axis,
                y_axis,
                curve='convex', # approximately
                direction='decreasing' # sklearn PCA eigenvalues are sorted
            ).knee
            temp_args['n_components'] = max(self.knee, 2)
            self.vprint("Knee found at {0} components. Using n={1}.".format(
                self.knee, temp_args['n_components'])
            )
            self.pca = PCA(**temp_args)
        else:
            self.pca = PCA(**self.args)
        return self.pca.fit_transform(x)

    def __getattr__(self, attr):
        """
        Allows user to return member values of pca object without
        having to remember or know any variable names.
        """
        return getattr(self.pca, attr)


class Dim_UMAP(Dim):
    def __init__(self, verbose=False, **args):
        super().__init__(verbose, **args)
        if 'n_components' not in args:
            raise ValueError("n_components not provided.")

    def get(self, x, y=None):
        self.umap = UMAP(**self.args)
        return self.umap.fit_transform(x, y=y)

    def __getattr__(self, attr):
        return getattr(self.umap, attr)


class Dim_TSNE(Dim):
    def __init__(self, verbose=False, **args):
        super().__init__(verbose, **args)
        if 'n_components' not in args:
            raise ValueError("n_components not provided.")

    def get(self, x, y=None):
        self.tsne = TSNE(**self.args)
        return self.tsne.fit_transform(x) # y is ignored

    def __getattr__(self, attr):
        return getattr(self.tsne, attr)