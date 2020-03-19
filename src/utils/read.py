from ast import literal_eval
from configparser import ConfigParser

import anndata
import gtfparse


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


def load_data(dataset):
    # return X, Y
    if dataset == 'spleen':
        ann = anndata.read_h5ad('datasets/spleen/dim_reduced_clustered.h5ad')
        return ann.X, ann.var.index.to_numpy().astype('U')
    elif dataset == 'brain':
        rnaseqtpm = pd.read_csv(
            'datasets/brain/RNAseqTPM.csv', index_col=0, header=None).T
        return rnaseqtpm.to_numpy(), rnaseqtpm.columns.to_numpy()
    elif dataset == 'brain_micro':
        brain_micro = pd.read_csv(
            'datasets/microarray_brain/curated_microarray.csv')
        return brain_micro.to_numpy()[:, 1:], brain_micro.columns.to_numpy()[1:]
    elif dataset == 'spellman':
        return pd.read_csv('datasets/Spellman.csv', index_col=0).to_numpy(), None
    else:
        raise FileNotFoundError('Dataset not supported.')
