from src import Pipeline
from src.utils.utils_experiment import load_data
import warnings
warnings.filterwarnings('ignore')

if __name__ == '__main__':
    dataset = 'spleen'
    X, ids = load_data(dataset)
    pipe = Pipeline(X, config='configs/'+dataset+".ini", verbose=True, col_ids=ids)
    pipe.run()

    plotter = Plotter(pipe)
    plotter.plot_clu()