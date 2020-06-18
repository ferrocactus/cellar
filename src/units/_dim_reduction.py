import logging

#import torch
import numpy as np
from kneed import KneeLocator
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from umap import UMAP

from ..log import setup_logger
from ._unit import Unit

#from src.methods import Autoencoder
# # Autoencoder
# N_COMPONENTS_AE = 15
# EPOCHS = 5
# BATCH = 64
# FACTOR = 0.1
# ACTIVATION = 'ReLU'
# LR = 1e-3
# WEIGHT_DECAY = 1e-5


class Dim_PCA(Unit):
    """
    Reduces the dimensionality of the data.
    See https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
    """

    def __init__(self, n_components='knee', n_components_max=100, **kwargs):
        """
        Parameters
        __________

        n_components: int, float or string, default 'knee'
            Number of components to use for PCA. If int, use that exact
            number of components. If float between 0 and 1, use that fraction
            of n_features of x. If 'knee', use Knee Detector algorithm to
            automatically figure out the n_components based on the curvature
            of the explained variance ratio graph.

        n_components_max: int, default 100
            If is ignored if n_components is different than 'knee'.
            Specifies the number of components to use for constructing the
            initial graph of explained variance ratio.

        **kwargs: dictionary
            Dictionary of parameters that will get passed to obj
            when instantiating it.
        """
        self.logger = setup_logger('PCA')
        self.n_components = n_components
        self.n_components_max = n_components_max
        self.kwargs = kwargs

    def get(self, x, return_evr=False):
        """
        Runs clustering for multiple n_clusters.

        Parameters
        __________

        x: array, shape (n_samples, n_features)
            The data array.

        return_evr: Bool
            If set, function will also return an array of
            explained variance ratios for every component.

        Returns
        _______

        x_emb: array, shape (n_samples, n_components)
            Data in the reduced dimensionality.

        [y_axis]: array, if return_evr==True, shape (n_components_max,)
            Explained Variance Ratio array.

        """
        self.logger.info("Initializing PCA.")

        if self.n_components == 'knee':
            n_components = min(self.n_components_max, np.min(x.shape))
            obj = PCA(n_components=n_components,
                      svd_solver='arpack', **self.kwargs).fit(x)

            # Construct axis for KneeLocator
            x_axis = list(range(1, obj.n_components_ + 1))
            y_axis = obj.explained_variance_ratio_
            # Find knee
            n_components = KneeLocator(x_axis, y_axis,
                                       curve='convex',  # Approximately
                                       direction='decreasing'
                                       ).knee
            n_components = max(n_components, 10)

            self.logger.info("Knee found at {0}.".format(n_components))
            x_emb = PCA(n_components=n_components,
                        svd_solver='arpack',
                        **self.kwargs).fit_transform(x)
        else:
            x_emb = PCA(n_components=self.n_components,
                        svd_solver='arpack',
                        **self.kwargs).fit_transform(x)

        if return_evr == True:
            return x_emb, y_axis
        else:
            return x_emb


class Dim_UMAP(Unit):
    def __init__(self, n_components=2, **kwargs):
        """
        Parameters
        __________

        use_y: Bool, default True
            If set, will also pass labels to UMAP for constrained
            dimensionality reduction.

        **kwargs: dictionary
            Dictionary of parameters that will get passed to obj
            when instantiating it.

        """
        import warnings
        from numba.errors import NumbaPerformanceWarning
        warnings.filterwarnings("ignore", category=NumbaPerformanceWarning)

        self.logger = setup_logger('UMAP')
        self.n_components = n_components
        self.kwargs = kwargs

    def get(self, x):
        """
        Reduces the dimensionality of the data.
        See https://umap-learn.readthedocs.io/en/latest/

        Parameters
        __________

        x: array, shape (n_samples, n_features)
            The data array.

        Returns
        _______

        x_emb: array, shape (n_samples, n_components)
            Data in the reduced dimensionality.

        """
        self.logger.info("Initializing UMAP.")
        return UMAP(n_components=self.n_components,
                    **self.kwargs).fit_transform(x)


class Dim_TSNE(Unit):
    def __init__(self, n_components=2, **kwargs):
        """
        Parameters
        __________

        **kwargs: dictionary
            Dictionary of parameters that will get passed to obj
            when instantiating it.

        """
        self.logger = setup_logger('TSNE')
        self.n_components = n_components
        self.kwargs = kwargs

    def get(self, x):
        """
        Reduces the dimensionality of the data.
        See https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html

        Parameters
        __________

        x: array, shape (n_samples, n_features)
            The data array.

        Returns
        _______

        x_emb: array, shape (n_samples, n_components)
            Data in the reduced dimensionality.
        """
        self.logger.info("Initializing TSNE.")
        return TSNE(**self.kwargs).fit_transform(x)


# class Dim_AE(Unit):
#     """
#     Use autoencoder to reduce dimensionality.
#     """

#     def __init__(self, verbose=False, name='AE', **kwargs):
#         super().__init__(verbose, name, **kwargs)
#         self.epochs = kwargs.get('epochs', EPOCHS)
#         self.batch = kwargs.get('batch', BATCH)
#         self.n_components = kwargs.get('n_components', N_COMPONENTS_AE)
#         self.factor = kwargs.get('factor', FACTOR)
#         self.activation = kwargs.get('activation', ACTIVATION)
#         self.lr = kwargs.get('lr', LR)
#         self.weight_decay = kwargs.get('weight_decay', WEIGHT_DECAY)

#     def get(self, x):
#         model = Autoencoder(
#             x.shape[1],
#             factor=self.factor,
#             output_dim=self.n_components,
#             non_linearity=self.activation
#         )

#         loss_f = torch.nn.MSELoss()
#         optimizer = torch.optim.Adam(
#             model.parameters(),
#             lr=self.lr,
#             weight_decay=self.weight_decay
#         )

#         self.vprint("Training Autoencoder.")
#         for epoch in range(self.epochs):
#             indices = np.arange(0, x.shape[0])
#             np.random.shuffle(indices)
#             for b in range(0, x.shape[0] // self.batch):
#                 data = x[indices[b * self.batch: (b + 1) * self.batch]]
#                 #device = torch.device("cpu")
#                 data = torch.tensor(data, device='cpu').float()
#                 # ===================forward=====================
#                 output = model(data)
#                 loss = loss_f(output, data)
#                 # ===================backward====================
#                 optimizer.zero_grad()
#                 loss.backward()
#                 optimizer.step()
#             # ===================log========================
#             self.vprint(
#                 f'epoch [{epoch+1}/{self.epochs}], loss:{loss.item():.4f}')
#         return model.encoder(
#             torch.tensor(x, device='cpu').float()
#         ).detach().numpy()
