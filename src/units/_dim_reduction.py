from ._unit import Unit
#from src.methods import Autoencoder

from abc import abstractmethod
from sklearn.decomposition import PCA
from kneed import KneeLocator
from umap import UMAP
from sklearn.manifold import TSNE
#import torch
import numpy as np

PCA_EVR_MAX_N = 100
N_COMPONENTS = 'knee'

# Autoencoder
N_COMPONENTS_AE = 15
EPOCHS = 5
BATCH = 64
FACTOR = 0.1
ACTIVATION = 'ReLU'
LR = 1e-3
WEIGHT_DECAY = 1e-5


class Dim(Unit):
    """
    Base class for Dimensionality Reduction methods.
    """

    def __init__(self, verbose=False, name='Dim', **kwargs):
        """
        Args:
            verbose (bool): Printing flag.
            **kwargs: Argument dict.
        """
        super().__init__(verbose, name, **kwargs)

    @abstractmethod
    def get(self, x):
        """
        Returns the embedding of x.

        Args:
            x (np.ndarray): Data in matrix (n x d) form.
        Returns:
            (np.ndarray): The embedding of x.
        """
        pass


class Dim_PCA(Dim):
    def __init__(self, verbose=False, name='PCA', **kwargs):
        super().__init__(verbose, name, **kwargs)
        if 'n_components' not in kwargs:
            raise ValueError("n_components not provided.")

    def get(self, x):
        if self.kwargs['n_components'] == 'knee':
            temp_kwargs = self.kwargs.copy()
            temp_kwargs['n_components'] = min(PCA_EVR_MAX_N, min(x.shape))
            self.ev_pca = PCA(**temp_kwargs)
            self.ev_pca.fit(x)
            # Construct axis for KneeLocator
            x_axis = list(range(1, self.ev_pca.n_components_ + 1))
            y_axis = self.ev_pca.explained_variance_ratio_
            # Find knee
            temp_kwargs['n_components'] = self.knee = KneeLocator(
                x_axis,
                y_axis,
                curve='convex',  # approximately
                direction='decreasing'  # sklearn PCA eigenvalues are sorted
            ).knee
            temp_kwargs['n_components'] = max(self.knee, 2)
            self.vprint("Knee found at {0} components. Using n={1}.".format(
                self.knee, temp_kwargs['n_components'])
            )
            self.pca = PCA(**temp_kwargs)
        else:
            self.pca = PCA(**self.kwargs)
        return self.pca.fit_transform(x)


class Dim_UMAP(Dim):
    def __init__(self, verbose=False, name='UMAP', **kwargs):
        super().__init__(verbose, name, **kwargs)
        if 'n_components' not in kwargs:
            raise ValueError("n_components not provided.")

    def get(self, x, y=None):
        self.umap = UMAP(**self.kwargs)
        return self.umap.fit_transform(x, y=y)


class Dim_TSNE(Dim):
    def __init__(self, verbose=False, name='TSNE', **kwargs):
        super().__init__(verbose, name, **kwargs)
        if 'n_components' not in kwargs:
            raise ValueError("n_components not provided.")

    def get(self, x, y=None):
        self.tsne = TSNE(**self.kwargs)
        return self.tsne.fit_transform(x)  # y is ignored


# class Dim_AE(Dim):
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
