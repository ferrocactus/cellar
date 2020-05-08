from ast import literal_eval
from configparser import ConfigParser
import os
import anndata
import gtfparse
import pandas as pd

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
    if (dataset == 'default' or dataset == "brain"):
        rnaseqtpm = pd.read_csv(
            'datasets/brain/RNAseqTPM.csv', index_col=0, header=None).T
        return rnaseqtpm.to_numpy(), rnaseqtpm.columns.to_numpy()
    elif dataset == 'spleen':
        ann = anndata.read_h5ad('datasets/spleen/dim_reduced_clustered.h5ad')
        return [ann.X, ann.var.index.to_numpy().astype('U')]
    elif dataset[-4:] == 'h5ad':
        filename=dataset[0:-5]
        ann = anndata.read_h5ad(str("datasets/"+filename+"/"+filename+".h5ad"))
        return [ann.X, ann.var.index.to_numpy().astype('U')]
    elif dataset[-3:] == 'csv':
        filename=dataset[0:-4]
        df = pd.read_csv(str("datasets/"+filename+"/"+filename+".csv"), index_col=0, header=None).T
        return df.to_numpy(), df.columns.to_numpy()
    else:
        return "error"
'''
def upload(path):
    if dataset[-4:] == 'h5ad':
        filename=dataset[0:-5]
        ann = anndata.read_h5ad(str("datasets/"+filename+"/"+filename+".h5ad"))
        ann.write_h5ad()
'''

def write_data(dataset,path):

    if path[-4:] == 'h5ad':
        filename=dataset[0:-5]
        ann = anndata.read_h5ad(str(path))
        os.mkdir("datasets/"+filename)
        ann.write_h5ad(str("datasets/"+filename+"/"+filename+".h5ad"))
    elif path[-3:] == 'csv':
        filename=dataset[0:-4]
        os.mkdir("datasets/"+filename)
        df = pd.read_csv(str(path)).T
        pd.write_csv(df,str("datasets/"+filename+"/"+filename+".csv"), index_col=0, header=None)
    else:
        return "error"

    return
