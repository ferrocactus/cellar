from src.workflow import Workflow
from src.utils_experiment import load_data
from matplotlib import pyplot as plt
import warnings
warnings.filterwarnings('ignore')

if __name__ == '__main__':
    plt.ion()

    dataset = 'spleen'
    X, Y = load_data(dataset)
    w = Workflow(X, Y, config=dataset, verbose=True)
    w.pca_plot_var_ratio()
    w.reduce_dim(method='pca')
    w.cluster(method='kmedoids')
    w.reduce_plot(labels=w.y_train_pred, method='umap')

    plt.ioff()
    plt.show()