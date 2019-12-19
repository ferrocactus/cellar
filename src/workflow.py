import numpy as np
from matplotlib import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import tqdm

from umap import UMAP
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

from sklearn.metrics import mean_squared_error as mse

class Workflow:
    def __init__(self, x, y=None):
        self.set_train_data(x, y)
        self.has_test_data = False
    
    def set_train_data(self, x, y=None, row_names=None, col_names=None):
        assert len(x.shape) == 2, "Data needs to be of shape (n x d)."
        self.x_train        = x
        self.n_train        = x.shape[0]
        self.dims           = x.shape[1]
        self.has_y_train    = False if y is None else True
        self.y_train        = y
        self.row_names      = row_names
        self.col_names      = col_names
    
    def set_test_data(self, x, y=None):
        assert len(x.shape) == 2, "Data needs to be of shape (n x d)."
        assert x.shape[1]   == self.dims, "Dimension mismatch between test & train."
        self.has_test_data  = True
        self.x_test         = x
        self.n_test         = x.shape[0]
        self.has_y_test     = False if y is None else True
        self.y_test         = y
    
    def set_row_names(self, row_names):
        self.row_names = row_names
    
    def set_col_names(self, col_names):
        self.col_names = col_names

    def normalize(self):
        raise NotImplementedError()
    
    def boxplot(self, data='train'):
        """
        Boxplot of the columns of the data
        """
        if self.has_y_train:
            print("WARNING: This is a boxplot of the columns of the data. " +
                    "Results may not be what expected.")
        sns.set_style("whitegrid")
        fig, ax = plt.subplots()
        if data == 'train':
            sns.boxplot(data=self.x_train)
        elif data == 'test' and self.has_test_data:
            sns.boxplot(data=self.x_test)
        else:
            raise NameError("Dataset not found.")
        sns.despine(left=True)
        if self.col_names is not None and len(self.col_names) < 40:
            plt.xticks(np.arange(self.dims), self.col_names, rotation=90)
        fig.set_size_inches(10, 5)
        sns.reset_defaults()
        plt.show()
    
    def reduce_and_plot(self, method='umap', dims=2, **kwargs):
        """
        Reduce the dimensionality of the data to 2D using the specified method
        and plot.
        params:
            method: Method to use for dimensionality reduction;
                    choose between: umap, pca, t-sne
            dims:   Plot dimension; Choose between 2 and 3
        """
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
            emb = UMAP(n_components=dims, **kwargs).fit_transform(self.x_train)
        elif method == 'pca':
            print("Reducing dimensionality using PCA.")
            emb = PCA(n_components=dims, **kwargs).fit_transform(self.x_train)
        elif method == 'tsne':
            print("Reducing dimensionality using T-SNE.")
            emb = TSNE(n_components=dims, **kwargs).fit_transform(self.x_train)
        else:
            raise NotImplementedError('Method not found')
        
        if dims == 2:
            if self.has_y_train:
                sns.scatterplot(x=emb[:, 0], y=emb[:, 1],
                                hue=self.y_train,
                                palette='Dark2',
                                linewidth=0,
                                s=10,
                                legend='full')
                plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            else:
                sns.scatterplot(x=emb[:, 0], y=emb[:, 1], linewidth=0, s=10)
        else:
            if self.has_y_train:
                ax.scatter(emb[:, 0], emb[:, 1], emb[:, 2],
                            c=[sns.color_palette()[i] for i in self.y_train])
            else:
                ax.scatter(emb[:, 0], emb[:, 1], emb[:, 2])
        sns.despine(left=True, bottom=True)
        plt.xticks([])
        plt.yticks([])
        plt.show()
    
    def reduce_dim(self, method='pca', **kwargs):
        """
        Reduces the dimensionality of the data and stores it in self.emb.
        """
        if method == 'pca':
            print("Reducing dimensionality using PCA.")
            pca = PCA(**kwargs)
            pca.fit(self.x_train)
            self.x_train_emb = pca.transform(self.x_train)

            # Print scores
            x_train_mserror = mse(self.x_train, pca.inverse_transform(self.x_train_emb))
            x_train_score = pca.score(self.x_train)
            print("Embedding created. Train MSE:", x_train_mserror)
            print("Train Average Log Likelihood:", x_train_score)

            if self.has_test_data:
                self.x_test_emb = pca.transform(self.x_test)
                x_test_mserror = mse(self.x_test, pca.inverse_transform(self.x_test_emb))
                x_test_score = pca.score(self.x_test)
                print("Embedding created. Test MSE:", x_test_mserror)
                print("Test Average Log Likelihood:", x_test_score)
        else:
            raise NotImplementedError()
    
    def pca_plot_var_ratio(self, n_components=None, **kwargs):
        """
        Plots the percentage of variance explained by each of the PCA components.
        """
        kwargs['n_components'] = n_components or self.dims
        pca = PCA(**kwargs)
        pca.fit(self.x_train)

        # Plotting
        sns.set_style("whitegrid")

        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.set_size_inches(10, 5)

        ax1.plot(list(range(1, kwargs['n_components'] + 1)),
                pca.explained_variance_ratio_)
        ax1.set_xlabel("Number of components")
        ax1.set_ylabel("Percentage of Explained Variance per Component")
        ax1.set_ylim(0)
        # Cumulative plot
        ax2.plot(list(range(1, kwargs['n_components'] + 1)),
                np.cumsum(pca.explained_variance_ratio_))
        ax2.set_xlabel("Number of components")
        ax2.set_ylabel("Cumulative Percentage of Explained Variance")
        ax2.set_ylim(0)

        sns.despine()
        plt.show()