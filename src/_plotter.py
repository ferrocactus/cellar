from ._wrapper import wrap

from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

COLORS = [
    '#cc5151', '#7f3333', '#51cccc', '#337f7f', '#8ecc51', '#597f33', '#8e51cc',
    '#59337f', '#ccad51', '#7f6c33', '#51cc70', '#337f46', '#5170cc', '#33467f',
    '#cc51ad', '#7f336c', '#cc7f51', '#7f4f33', '#bccc51', '#757f33', '#60cc51',
    '#3c7f33', '#51cc9e', '#337f62', '#519ecc', '#33627f', '#6051cc', '#3c337f'
]

class Plotter:
    """
    Class used for plotting intermediary steps of a Pipeline object.
    """
    def __init__(self, pipe):
        """
        Args:
            pipe (Pipeline): Pipeline object.
        """
        self.pipe = pipe

    def plot_dim(self, cumulative=True):
        """
        Plots explained variance ratio per component.
        Args:
            cumulative (Boolean): If set, show extra plot of cumulative var.
        """
        fig, ax = plt.subplots(1, cumulative+1, squeeze=False)

        y = self.pipe.dim_obj.ev_pca.explained_variance_ratio_
        x = list(range(1, y+1))

        ax[0][0].plot(x, y*100)
        ax[0][0].set_xlabel("Number of components")
        ax[0][0].set_ylabel("Percentage of Explained Variance per Component")
        ax[0][0].set_ylim(0)

        if cumulative:
            ax[0][1].plot(x, np.cumsum(y*100))
            ax[0][1].set_xlabel("Number of components")
            ax[0][1].set_ylabel("Cumulative Percentage of Explained Variance")
            ax[0][1].set_ylim(0)

        ax[0][0].set_title('PCA')

        sns.despine()
        fig.set_size_inches(5*(cumulative + 1), 5)
        plt.show()

    def plot_clu(self):
        """
        Plots the clusters in 2D for a single k.
        """
        hue = self.pipe.labels
        x = self.pipe.vis_obj.get(self.pipe.x_emb, hue)
        if x.shape[1] != 2:
            raise ValueError('Incorrect number of dimensions.')

        labels = "sdf"

        unq = len(np.unique(hue))
        pal = COLORS[:unq]
        for i in range(len(COLORS), unq): # In case more colors are needed.
            pal.append("#" + '%06x' % np.random.randint(16**6))

        plt.figure(figsize=(10,5))
        sns.scatterplot(x=x[:, 0], y=x[:, 1], hue=hue, palette=pal, linewidth=0,
                                                            s=10, legend='full')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0, labels=labels)

        sns.despine(left=True, bottom=True)
        plt.xticks([])
        plt.yticks([])
        plt.show()

    def plot_clu_all(self, cols=3):
        """
        Plots the clusters in 2D for different k.
        Args:
            cols (int): Number of columns to use in the plot.
        """
        clu_method = self.pipe.config["methods"]["cluster"]
        kwargs = self.pipe.clu_obj.kwargs.copy() # Copy original kwargs
        k = kwargs['n_clusters']

        if not isinstance(k, tuple): # Single cluster_n
            self.plot_clu()
            return

        # Get the embedding using vis_obj
        emb = self.pipe.vis_obj.get(self.pipe.x_emb)
        if emb.shape[1] != 2:
            raise ValueError("Incorrect number of dimensions.")

        ks = list(range(*k))
        rows = (len(ks) + cols - 1) // cols
        fig, ax = plt.subplots(rows, cols, squeeze=False)

        # Iterate over all k
        for i, kk in enumerate(ks):
            kwargs['n_clusters'] = kk
            clu_obj = wrap("cluster", clu_method)(**kwargs)
            labels = clu_obj.get(self.pipe.x_emb, self.pipe.eval_obj)

            ax[i // cols][i % cols].scatter(emb[:, 0], emb[:, 1], s=.5, c=labels)
            ax[i // cols][i %cols].set_title('k={0}, score={1:.2f}'.format(
                                                kk, clu_obj.score_list))
            ax[i // cols][i %cols].set_xticks([])
            ax[i // cols][i %cols].set_yticks([])

        # Remove unused boxes
        for i in range(len(ks), cols*rows):
            ax[i // cols][i % cols].spines['right'].set_visible(False)
            ax[i // cols][i % cols].spines['top'].set_visible(False)
            ax[i // cols][i % cols].spines['left'].set_visible(False)
            ax[i // cols][i % cols].spines['bottom'].set_visible(False)
            ax[i // cols][i % cols].set_xticks([])
            ax[i // cols][i % cols].set_yticks([])

        fig.set_size_inches(5 * cols, 5 * rows)
        plt.show()

    def plot_eval(self):
        """
        Plots the scores of the clusters for different k. Branch if k is int.
        """
        fig, ax = plt.subplots()

        k = self.pipe.clu_obj.kwargs['n_clusters']
        if isinstance(k, tuple):
            x = list(range(*k))
            ax.set_xticks(x)
        else:
            x = k
            ax.set_xticks([x-1, x, x+1])

        y = self.pipe.clu_obj.score_list

        ax.bar(x, y)
        ax.set_xlabel('Number of Clusters')
        ax.set_ylabel('Score')
        sns.despine()
        fig.set_size_inches(10, 5)
        plt.show()