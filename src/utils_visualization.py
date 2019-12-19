from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import seaborn as sns
import umap
from sklearn.decomposition import PCA

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

def reduce_and_plot(x, y=None, ax=None,
                                method='umap',
                                emb=None,
                                dims=2,
                                n_neighbors=15,
                                min_dist=0.1):
    """
    Gives a "dims" dimensional plot of "x" using UMAP or PCA
    params:
        x:      nunmpy array of size (d, dims); ignored if "emb" is set
        y:      integer labels for data points in x (will be used for coloring);
                    if None, no coloring applied
        ax:     pyplot axis; if None, new axis will be created
        method: method to use for dimensionality reduction; 'umap' or 'pca'
        dims:   2 or 3 dimensions; ignored if "emb" is set
        emb:    reduced embedding to plot
    """
    assert dims == 2 or dims == 3, "Can only visualize 2 or 3 dimensions"
    axis_is_set = False if ax is None else True

    if emb is None:
        if method == 'umap':
            emb = umap.UMAP(n_components=dims,
                            n_neighbors=n_neighbors,
                            min_dist=min_dist).fit_transform(x)
        elif method == 'pca':
            emb = PCA(n_components=dims).fit_transform(x)
    else:
        assert emb.shape[1] == 2 or emb.shape[1] == 3, "Can only visualize 2 or 3 dimensions"

    if dims == 2:
        if not axis_is_set:
            fig, ax = plt.subplots()
        if y is not None and len(y) > 0:   # colors
            ax.scatter(emb[:, 0], emb[:, 1], c=[sns.color_palette()[i] for i in y])
        else:
            ax.scatter(emb[:, 0], emb[:, 1])
    elif dims == 3:
        if not axis_is_set:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        if y is not None and len(y) > 0:   # colors
            ax.scatter(emb[:, 0], emb[:, 1], emb[:, 2], c=[sns.color_palette()[i] for i in y])
        else:
            ax.scatter(emb[:, 0], emb[:, 1], emb[:, 2])
    sns.despine()
    if not axis_is_set:
        plt.show()
    else:
        return ax

def plot_pca_var_ratio(x, ax=None, n_components=None, cumulative=False):
    """
    Plot n_components vs variance_ratio
    """
    n_components = n_components or x.shape[1]
    print(n_components)
    pca = PCA(n_components=n_components)
    pca.fit(x)
    var = pca.explained_variance_ratio_
    
    axis_is_set = False if ax is None else True
    if not axis_is_set:
        fig, ax = plt.subplots()
    
    if cumulative:
        cumulative = np.cumsum(var)
        ax.plot(range(1, n_components+1), cumulative)
    else:
        ax.plot(range(1, n_components+1), var)

    sns.despine()
    if not axis_is_set:
        plt.show()
    else:
        return ax