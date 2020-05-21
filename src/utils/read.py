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


def load_data(dataset, dataset_source=None):
    # return X, Y
    if (dataset == 'default' or dataset == "brain"):
        rnaseqtpm = pd.read_csv(
            'datasets/brain/RNAseqTPM.csv', index_col=0, header=None).T
        return rnaseqtpm.to_numpy(), rnaseqtpm.columns.to_numpy()
    elif dataset == 'spleen':
        ann = anndata.read_h5ad('datasets/spleen/dim_reduced_clustered.h5ad')
        return [ann.X, ann.var.index.to_numpy().astype('U')]
    else:
        if dataset_source is None:
            raise ValueError("Dataset source not specified.")

        if dataset_source == 'hubmap':
            if len(dataset) <= 5 or dataset[-4:] != 'h5ad':
                dataset = dataset + ".h5ad"
            ann = anndata.read_h5ad("datasets/" + dataset_source + "/" + dataset)
            return [ann.X, ann.var.index.to_numpy().astype('U')]
        elif dataset_source == 'user_uploaded':
            if dataset[-4:] == 'h5ad':
                ann = anndata.read_h5ad(
                    "datasets/" + dataset_source + "/" + dataset)
                return [ann.X, ann.var.index.to_numpy().astype('U')]
            elif dataset[-3:] == 'csv':
                df = pd.read_csv("datasets/" + dataset_source + "/" + dataset,
                                 index_col=0, header=None).T
                return df.to_numpy(), df.columns.to_numpy()
            else:
                raise ValueError("Dataset format not implemented.")
        else:
            raise ValueError("Dataset source not found.")


def upload_file(dataset, path):
    try:
        filename, file_extension = os.path.splitext(path)
        shutil.move(path, "datasets/user_uploaded/" + dataset + file_extension)
    except Exception as e:
        print(str(e))
        raise OSError("A problem occured when reading dataset.")
