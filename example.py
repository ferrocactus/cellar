from src import Pipeline
from src.utils.utils_experiment import load_data
import warnings
warnings.filterwarnings('ignore')

if __name__ == '__main__':
    dataset = 'spleen'
    X, Y, gene_ids = load_data(dataset)
    pipe = Pipeline(X, config='configs/'+dataset+".ini", verbose=True, col_ids=gene_ids)
    pipe.run()
    pipe.plot('2d')