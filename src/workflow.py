import numpy as np
from matplotlib import pylab as plt
import seaborn as sns

class Workflow:    
    def __init__(self, x, y=None):
        self.set_train_data(x, y)
    
    def set_train_data(self, x, y=None, row_names=None, col_names=None):
        assert len(x.shape)   == 2, "Data needs to be of shape (n x d)."
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
        sns.set_style("whitegrid")
        fig, ax = plt.subplots()
        if data == 'train':  sns.boxplot(data=self.x_train)
        elif data == 'test': sns.boxplot(data=self.x_test)
        else:                raise NameError()
        sns.despine(left=True)
        if self.col_names is not None and len(self.col_names) < 40:
            plt.xticks(np.arange(self.dims), self.col_names, rotation=90)
        fig.set_size_inches(10, 5)
        sns.reset_defaults()

    def heatmap(self):
        sns.heatmap(self.x_train)