import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import seaborn as sns
from umap import UMAP
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

sns.set_style('whitegrid')

def plot_marker_hist(n_clusters, pvals, mads):
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

def plot_top_markers(marker_ids, marker_pvals, marker_mds):
    clusters = marker_ids.shape[0]
    fig, ax = plt.subplots(int((clusters + 1)/2), 2)
    for i in range(clusters):
        x = np.arange(len(marker_pvals[i]))
        y = marker_mds[i]
        ax[int(i/2)][i%2].scatter(x, y, c=np.arange(len(marker_pvals[i])))
        ax[int(i/2)][i%2].set_xticklabels(np.round(marker_pvals[i], 3))
        ax[int(i/2)][i%2].set_title("Cluster " + str(i+1))
        ax[int(i/2)][i%2].set_xlabel("p-value")
        ax[int(i/2)][i%2].set_ylabel("mean difference")
        for j, txt in enumerate(marker_ids[i]):
            ax[int(i/2)][i%2].text(x[j], y[j]+0.02, marker_ids[i][j], fontsize=10, rotation=90)

    sns.despine()
    fig.set_size_inches(10, int((clusters + 1)/2)*8)

