from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import seaborn as sns

from umap import UMAP
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

def plot_images(imlist, cols=3, titlelist=None, cmap='gray'):
    """
    Plots numpy images in a grid with "cols" columns.
    """
    n = imlist.shape[0]
    rows = int((n + cols - 1) / cols)
    for i, im in enumerate(imlist):
        plt.subplot(rows, cols, i + 1)
        plt.axis('off')
        plt.imshow(im, cmap=cmap)
        if titlelist is not None:
            plt.title(titlelist[i])
    plt.show()

def reduce_and_plot(x=None, y=None, method='umap', dims=2, **kwargs):
    """
    Reduce the dimensionality of the data to 2D using the specified method
    and plot.
    params:
        method: Method to use for dimensionality reduction;
                choose between: umap, pca, t-sne
        dims:   Plot dimension; Choose between 2 and 3
    """
    assert x is not None or emb is not None, "Data not provided."

    fig = plt.figure()
    if   dims == 2:     ax = fig.add_subplot(111)
    elif dims == 3:     ax = fig.add_subplot(111, projection='3d')
    else:               raise NotImplementedError('Can only visualize 2 or 3 dims.')

    if method == 'umap':
        """
        Possible arguments for UMAP and their default values
            n_neighbors     = 15
            min_dist        = 0.1
            n_components    = 2
            metric          = 'euclidean'
        """
        print("Reducing dimensionality using UMAP.")
        emb = UMAP(n_components=dims, **kwargs).fit_transform(x)
    elif method == 'pca':
        print("Reducing dimensionality using PCA.")
        emb = PCA(n_components=dims, **kwargs).fit_transform(x)
    elif method == 'tsne':
        print("Reducing dimensionality using T-SNE.")
        emb = TSNE(n_components=dims, **kwargs).fit_transform(x)
    else:
        raise NotImplementedError('Method not found')
    
    if dims == 2:
        if y is not None:
            sns.scatterplot(x=emb[:, 0], y=emb[:, 1],
                            hue=y,
                            palette='Dark2',
                            linewidth=0,
                            s=10,
                            legend='full')
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        else:
            sns.scatterplot(x=emb[:, 0], y=emb[:, 1], linewidth=0, s=10)
    else:
        if y is not None:
            ax.scatter(emb[:, 0], emb[:, 1], emb[:, 2],
                        c=[sns.color_palette()[i] for i in y])
        else:
            ax.scatter(emb[:, 0], emb[:, 1], emb[:, 2])
    sns.despine(left=True, bottom=True)
    plt.xticks([])
    plt.yticks([])
    fig.set_size_inches(10, 5)