from acip.acip import ACIP
from utils.utils_experiment import load_data
from matplotlib import pyplot as plt
import warnings
warnings.filterwarnings('ignore')

if __name__ == '__main__':
    plt.ion()

    #plt.ioff()
    dataset = 'spleen'
    X, Y, gene_ids = load_data(dataset)
    w = ACIP(X, Y, config=dataset, verbose=True)
    w.set_col_names(gene_ids)
    w.filter_genes()
    w.reduce_dim(method='pca')
    w.cluster(method='kmedoids')
    w.reduce_plot(method='umap')
    w.get_markers()
    print(w.markers)

    plt.ioff()
    plt.show()