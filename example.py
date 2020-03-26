from src import Pipeline, Plotter
from src.utils.read import load_data
import warnings
warnings.filterwarnings("ignore")

if __name__ == '__main__':
    X, ids = load_data('spleen')
    pipe = Pipeline(X, config='configs/config.ini', col_ids=ids)
    pipe.run()

    plotter = Plotter(pipe)
    plotter.plot_clu()
