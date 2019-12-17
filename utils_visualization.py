from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import umap
from sklearn.decomposition import PCA

def plot_images(imlist, cols=3, titlelist=None):
    """
    Plots numpy images in a grid with "cols" columns.
    """
    n = imlist.shape[0]
    rows = int((n + cols - 1) / cols)
    for i, im in enumerate(imlist):
        plt.subplot(rows, cols, i + 1)
        plt.imshow(im)
        if titlelist is not None:
            plt.title(titlelist[i])
    plt.show()

def reduce_and_plot(x, y=None, ax=None, method='umap', dims=2, emb=None):
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

    if emb is None:
        if method == 'umap':
            emb = umap.UMAP(n_components=dims).fit_transform(x)
        elif method == 'pca':
            emb = PCA(n_components=dims).fit_transform(x)
    else:
        assert emb.shape[1] == 2 or emb.shape[1] == 3, "Can only visualize 2 or 3 dimensions"
    
    if ax == None:
        fig, ax = plt.subplots()

    if dims == 2:
        if y is not None:   # colors
            ax.scatter(emb[:, 0], emb[:, 1], c=[sns.color_palette()[i] for i in y])
        else:
            ax.scatter(emb[:, 0], emb[:, 1])
    elif dims == 3:
        if y is not None:   # colors
            ax.scatter(emb[:, 0], emb[:, 1], emb[:, 2], c=[sns.color_palette()[i] for i in y])
        else:
            ax.scatter(emb[:, 0], emb[:, 1], emb[:, 2])
    sns.despine()
    return ax

def plot_pca_var_ratio(x, ax=None):
    """
    Plot n_components vs variance_ratio
    """
    pca = PCA(n_components=x.shape[1])
    pca.fit(x)
    var = pca.explained_variance_ratio_
    
    if ax == None:
        fig, ax = plt.subplots()
    
    ax.plot(range(1, x.shape[1]+1), var)
    sns.despine()
    return ax