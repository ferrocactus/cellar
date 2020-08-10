from ast import literal_eval
from configparser import ConfigParser
import anndata
import os
import shutil
import traceback
import sys
from .exceptions import IncorrectFileFormat
from bidict import bidict

this_dir = os.path.dirname(__file__)
def join_root(path):
    return os.path.abspath(os.path.join(this_dir, path))


def safe_load_file(filepath):
    try:
        return load_file(filepath)
    except IncorrectFileFormat as e:
        print("Incorrect file format")
        return "file_error"
    except Exception as e:
        traceback.print_exc(file=sys.stdout)
        return "file_error"


def load_file(filepath, backed='r+'):
    t_flag = False
    if filepath == 'default' or filepath == 'datasets/user_uploaded/default':
        filepath = join_root("../datasets/default.csv")
        t_flag = True
    elif filepath == 'test':
        filepath = join_root('../../datasets/server/testdataset.h5ad')

    dataset = os.path.basename(filepath)
    dataset = os.path.splitext(dataset)[0]

    try:
        if filepath[-4:] == 'h5ad':
            adata = anndata.read_h5ad(filepath, backed=backed)
        if filepath[-3:] == 'csv':
            adata = anndata.read_csv(filepath)
            if t_flag:
                adata = adata.T
        if filepath[-4:] == 'xlsx':
            adata = anndata.read_excel(filepath)
        if filepath[-3:] == 'mtx':
            adata = anndata.read_mtx(filepath)
        if filepath[-3:] == 'txt' or filepath[-3:] == 'tab' or filepath[-4:] == 'data':
            adata = anndata.read_text(filepath)
        if filepath[-2:] == 'h5':
            adata = anndata.read_hdf(filepath)
        if filepath[-4:] == 'loom':
            adata = anndata.read_loom(filepath)
    except Exception as e:
        print(str(e))
        raise IncorrectFileFormat(
            "File does not exist or file format is incorrect.")

    # Make sure cluster names are in proper format
    if 'cluster_names' in adata.uns:
        adata.uns['cluster_names'] = bidict(adata.uns['cluster_names'])
        for key in list(adata.uns['cluster_names'].keys()):
            adata.uns['cluster_names'][int(key)] = \
                adata.uns['cluster_names'].pop(key, None)

    adata.uns['dataset'] = dataset
    return adata


def upload_file(dataset, path):
    try:
        filename, file_extension = os.path.splitext(path)
        shutil.move(path, join_root("../../datasets/user_uploaded") + "/" + dataset + file_extension)
    except Exception as e:
        traceback.print_exc(file=sys.stdout)
        raise OSError("A problem occured when reading dataset.")
