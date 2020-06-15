from ast import literal_eval
from configparser import ConfigParser
import anndata
import gtfparse
import pandas as pd
import os
import shutil


def read_h5ad(filename):
    ann = anndata.read_h5ad(filename)
    return ann


def read_gtf(filename):
    gtf = gtfparse.read_gtf(filename)
    return gtf


def dict_literal_eval(d):
    return {key: literal_eval(d[key]) for key in d}


def parse_config(filename):
    configp = ConfigParser()
    configp.read(filename)

    config = configp._sections
    for key in config:
        config[key] = dict_literal_eval(config[key])
    return config


def read_config(dataset):
    with open("configs/" + dataset + ".json", "r") as f:
        config = literal_eval(f.read())
    return config


def load_data(dataset, check_precomputed_PCA=False):
    # return dict with keys x, y, precomputed_pca, col_ids, row_ids
    df = {}
    if (dataset == 'default' or dataset == "brain"):
        rnaseqtpm = pd.read_csv(
            'datasets/brain/RNAseqTPM.csv', index_col=0, header=None).T
        df['x'] = rnaseqtpm.to_numpy()
        df['col_ids'] = rnaseqtpm.columns.to_numpy()
        df['row_ids'] = None
        df['precomputed_pca'] = None
    elif dataset == 'spleen':
        ann = anndata.read_h5ad('datasets/spleen/dim_reduced_clustered.h5ad')
        df['x'] = ann.X
        df['col_ids'] = ann.var.index.to_numpy().astype('U')
        df['row_ids'] = ann.obs.index.to_numpy().astype('U')
        df['precomputed_pca'] = ann.obsm['X_pca']
    else:
        if dataset[-4:] == 'h5ad':
            ann = anndata.read_h5ad("datasets/"  + dataset)
            df['x'] = ann.X
            df['col_ids'] = ann.var.index.to_numpy().astype('U')
            df['row_ids'] = ann.obs.index.to_numpy().astype('U')
            df['precomputed_pca'] = ann.obsm['X_pca']
        elif dataset[-3:] == 'csv':
            csvdf = pd.read_csv("datasets/"  + dataset,
                             index_col=0, header=None).T
            df['x'] = csvdf.to_numpy()
            df['col_ids'] = csvdf.columns.to_numpy()
            df['row_ids'] = None
            df['precomputed_pca'] = None
        else:
            raise ValueError("Dataset format not implemented.")
    return df


def upload_file(dataset, path):
    try:
        filename, file_extension = os.path.splitext(path)
        shutil.move(path, "datasets/user_uploaded/" + dataset + file_extension)
    except Exception as e:
        print(str(e))
        raise OSError("A problem occured when reading dataset.")
