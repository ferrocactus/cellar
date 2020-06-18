from ast import literal_eval
from configparser import ConfigParser
import anndata
import gtfparse
import pandas as pd
import os
import shutil


def load_file(filepath):
    if filepath == 'default':
        filepath = 'datasets/brain/RNAseqTPM.csv'
    elif filepath == 'test':
        filepath = 'datasets/hubmap/testdataset.h5ad'

    if filepath[-4:] == 'h5ad':
        return anndata.read_h5ad(filepath)
    if filepath[-3:] == 'csv':
        # TODO remove transpose
        return anndata.read_csv(filepath).T
    if filepath[-4:] == 'xlsx':
        return anndata.read_excel(filepath)
    if filepath[-3:] == 'mtx':
        return anndata.read_mtx(filepath)
    if filepath[-3:] == 'txt' or filepath[-3:] == 'tab' or filepath[-4:] == 'data':
        return anndata.read_text(filepath)
    if filepath[-2:] == 'h5':
        return anndata.read_hdf(filepath)
    if filepath[-4:] == 'loom':
        return anndata.read_loom(filepath)


def load_data(dataset):
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
            ann = anndata.read_h5ad("datasets/" + dataset)
            df['x'] = ann.X
            df['col_ids'] = ann.var.index.to_numpy().astype('U')
            df['row_ids'] = ann.obs.index.to_numpy().astype('U')
            if hasattr(ann, 'obsm') and 'X_pca' in ann.obsm:
                df['precomputed_pca'] = ann.obsm['X_pca']
            else:
                df['precomputed_pca'] = None
        elif dataset[-3:] == 'csv':
            csvdf = pd.read_csv("datasets/" + dataset,
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
