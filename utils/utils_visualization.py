from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import seaborn as sns
from umap import UMAP
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

sns.set_style('whitegrid')

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

def plot_explained_variance(x, y, cumulative=True):
    fig, ax = plt.subplots(1, cumulative+1, squeeze=False)

    ax[0][0].plot(x, y*100)
    ax[0][0].set_xlabel("Number of components")
    ax[0][0].set_ylabel("Percentage of Explained Variance per Component")
    ax[0][0].set_ylim(0)
    if cumulative:
        ax[0][1].plot(x, np.cumsum(y*100))
        ax[0][1].set_xlabel("Number of components")
        ax[0][1].set_ylabel("Cumulative Percentage of Explained Variance")
        ax[0][1].set_ylim(0)

    sns.despine()
    fig.set_size_inches(5*(cumulative + 1), 5)

def plot_gene_variances(gene_variances):
    fig, ax = plt.subplots()
    ax.hist(gene_variances)
    ax.set_xlabel("Variance")
    ax.set_ylabel("Gene Count")
    sns.despine()
    fig.set_size_inches(10, 5)

def plot_scores(x, y):
    fig, ax = plt.subplots()
    ax.plot(x, y)
    ax.set_xlabel('Number of Clusters')
    ax.set_ylabel('Score')
    sns.despine()
    fig.set_size_inches(10, 5)

def plot_markers(n_clusters, pvals, mads):
    fig, ax = plt.subplots(n_clusters, 2)

    for cluster_id in range(n_clusters):
        ax[cluster_id][0].hist(pvals[cluster_id])
        ax[cluster_id][0].set_xlabel("p-values")
        ax[cluster_id][0].set_ylabel("gene count")
        ax[cluster_id][0].set_title("Cluster:" + str(cluster_id))
        ax[cluster_id][1].hist(mads[cluster_id], color='r')
        ax[cluster_id][1].set_xlabel("absolute difference")
        ax[cluster_id][1].set_ylabel("gene count")
        ax[cluster_id][1].set_title("Cluster:" + str(cluster_id))
    
    sns.despine()
    fig.set_size_inches(10, n_clusters * 5)

def plot_2d(x, y=None, dims=2):
    fig = plt.figure()
    if dims == 2: ax = fig.add_subplot(111)
    elif dims == 3: ax = fig.add_subplot(111, projection='3d')
    else: raise NotImplementedError('Can only visualize 2 or 3 dims.')

    if dims == 2:
        if y is not None:
            pal = []
            for i in range(len(np.unique(y))):
                pal.append("#" + '%06x' % np.random.randint(16**6))
            sns.scatterplot(x=x[:, 0], y=x[:, 1],
                            hue=y,
                            palette=pal,
                            linewidth=0,
                            s=10,
                            legend='full')
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        else:
            sns.scatterplot(x=x[:, 0], y=x[:, 1], linewidth=0, s=10)
    else:
        if y is not None:
            ax.scatter(x[:, 0], x[:, 1], x[:, 2],
                        c=[sns.color_palette()[i] for i in y])
        else:
            ax.scatter(x[:, 0], x[:, 1], x[:, 2])
    sns.despine(left=True, bottom=True)
    plt.xticks([])
    plt.yticks([])
    fig.set_size_inches(10, 5)