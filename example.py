from src import Pipeline, Plotter
from src.utils.read import load_data
import warnings
warnings.filterwarnings("ignore")

if __name__ == '__main__':
    pipe = Pipeline("hubmap/testdataset.h5ad")
    pipe.run_all()

    plotter = Plotter(pipe)
    plotter.plot_clu()
