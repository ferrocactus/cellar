from ._wrapper import wrap

from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

COLORS = [
    '#cc5151', '#51cccc', '#337f7f', '#8ecc51', '#7f3333', '#597f33', '#8e51cc',
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

        y = self.pipe.dim.ev_pca.explained_variance_ratio_
        x = list(range(1, len(y)+1))

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

    def plot_clu(self, use_markers=True):
        """
        Plots the clusters in 2D for a single k.
        """
        hue = self.pipe.labels
        x = self.pipe.vis.get(self.pipe.x_emb, hue)
        if x.shape[1] != 2:
            raise ValueError('Incorrect number of dimensions.')

        labels = []
        if use_markers == True:
            for label in self.pipe.markers:
                lvl1_type = self.pipe.markers[label]['lvl1_type']
                lvl1_sv = self.pipe.markers[label]['lvl1_sv']
                lvl1_intersec_n =len(self.pipe.markers[label]['lvl1_intersec'])
                lvl1_total = self.pipe.markers[label]['lvl1_total']
                lvl2_type = self.pipe.markers[label]['lvl2_type']
                lvl2_sv = self.pipe.markers[label]['lvl2_sv']
                lvl2_total = self.pipe.markers[label]['lvl2_total']
                lvl2_intersec_n =len(self.pipe.markers[label]['lvl2_intersec'])

                labels.append("Lvl1: {0}, sv={1:.2f}, IoU={2}/{3}".format(
                        lvl1_type, lvl1_sv, lvl1_intersec_n, lvl1_total
                    ) + "\nLvl2: {0}, sv={1:.2f}, IoU={2}/{3}".format(
                        lvl2_type, lvl2_sv, lvl2_intersec_n, lvl2_total
                    )
                )

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
        kwargs = self.pipe.clu.kwargs.copy() # Copy original kwargs
        k = kwargs['n_clusters']

        if not isinstance(k, tuple): # Single cluster_n
            self.plot_clu()
            return

        # Get the embedding using vis
        if not hasattr(self.pipe, 'x_emb_2d'):
            emb = self.pipe.vis.get(self.pipe.x_emb)
        else:
            emb = self.pipe.x_emb_2d

        ks = list(range(*k))
        rows = (len(ks) + cols - 1) // cols
        fig, ax = plt.subplots(rows, cols, squeeze=False)

        # Iterate over all k
        for i, kk in enumerate(ks):
            kwargs['n_clusters'] = kk
            clu = wrap("cluster", clu_method)(**kwargs)
            labels = clu.get(self.pipe.x_emb, self.pipe.eval)

            ax[i // cols][i % cols].scatter(emb[:, 0], emb[:, 1], s=.5, c=labels)
            ax[i // cols][i %cols].set_title('k={0}, score={1:.2f}'.format(
                                                kk, clu.score_list[0]))
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

        k = self.pipe.clu.kwargs['n_clusters']
        if isinstance(k, tuple):
            x = list(range(*k))
            ax.set_xticks(x)
        else:
            x = k
            ax.set_xticks([x-1, x, x+1])

        y = self.pipe.clu.score_list

        ax.bar(x, y)
        ax.set_xlabel('Number of Clusters')
        ax.set_ylabel('Score')
        sns.despine()
        fig.set_size_inches(10, 5)
        plt.show()

    def plot_mark(self, convention="names"):
        """
        Plots marker information.
        """
        n_clusters = len(self.pipe.unq_labels)

        fig, ax = plt.subplots(n_clusters, 2)

        for i, label in enumerate(self.pipe.unq_labels):
            pvals = self.pipe.markers[label]['pvals']
            diffs = self.pipe.markers[label]['diffs']
            names = self.pipe.markers[label]['outp_names']

            ax[i][0].hist(pvals, color='b', alpha=.4, label="pvals")
            ax[i][0].set_xlabel("p-values")
            ax[i][0].set_ylabel("gene count")
            ax[i][0].legend(loc=1)
            ax[i][0].spines['bottom'].set_color('blue')

            ax2 = ax[i][0].twiny()
            ax2.hist(diffs, color='r', alpha=.4, label="diffs")
            ax2.set_xlabel("absolute difference")
            ax2.xaxis.tick_top()
            ax2.xaxis.set_label_position('top')
            ax2.legend(loc=5)
            ax2.spines['top'].set_color('red')

            k = 10
            x = np.arange(len(pvals[:k]))
            ax[i][1].scatter(x, diffs[:k], c=np.arange(len(pvals[:k])))
            ax[i][1].set_xticklabels(np.round(pvals[:k], 3))
            ax[i][1].set_title(f"Cluster: {label}")
            ax[i][1].set_xlabel("p-value")
            ax[i][1].set_ylabel("mean difference")
            for j, txt in enumerate(names[:k]):
                ax[i][1].text(x[j], diffs[j]+0.02, names[j], fontsize=10, rotation=90)

        fig.set_size_inches(10, n_clusters * 6)
        plt.show()