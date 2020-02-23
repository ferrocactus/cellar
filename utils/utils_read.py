import anndata
import gtfparse
from ast import literal_eval
from configparser import ConfigParser

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