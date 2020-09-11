#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 01:12:37 2020

@author: jingtao
"""


import os 
from typing import Optional, Union
from ast import literal_eval
from configparser import ConfigParser
import anndata
import os
import shutil
import traceback
import sys
from scipy.sparse import issparse

#from bidict import bidict
import numpy as np
from anndata import AnnData
from bidict import bidict
from scipy.sparse import issparse



from anndata import AnnData



def reduce_dim(
        x: Union[AnnData, np.ndarray, list],
        method: str = "Truncated SVD",
        n_components: Union[str, int, float] = 'knee',
        inplace: Optional[bool] = True,
        check_if_exists: Optional[bool] = False,
        clear_dependents: Optional[bool] = True,
        **kwargs) -> Optional[Union[AnnData, np.ndarray]]:
    """
    Reduce dimensionality of the data.

    Parameters
    __________
    x: AnnData object or np array containing the data matrix.

    method: String specifying the dimensionality reduction method
        to use. See https://github.com/ferrocactus/cellar/tree/master/doc
        for a full list of methods available.

    n_components: Number of components to use. If method is 'PCA', this
        parameters can also be 'knee' in which case the number of
        components will be determined automatically based on the
        ankle heuristic.

    inplace: Only used when x is an AnnData object.
        If set to true will update x.obsm['x_emb'], otherwise,
        return a new AnnData object.

    check_if_exists: Only used when x is an AnnData object.
        If set to true, will check x if it contains embeddings
        from the same method. If found, and the number of components
        in x.obsm['x_emb'] is greater than or equal to n_components,
        then the first n_components columns of x.obsm['x_emb']
        will be returned. No computation takes place.

    clear_dependents: If set to true and x_emb changes, then will also
        clear labels and 2d embeddings if any exist.

    **kwargs: Additional parameters that will get passed to the
        clustering object as specified in method. For a full list
        see the documentation of the corresponding method.

    Returns
    _______
    If x is an AnnData object, will either return an AnnData object
    or None depending on the value of inplace.
    If x is not AnnData, will return np.ndarray of shape
    (n_observations, n_components) containing the embedding.
    """
    # Validations
    is_AnnData = isinstance(x, AnnData)
    if issparse(adata.X):
        allowed = ['Truncated SVD', 'Spectral Embedding',
                   'UMAP', 'Diffusion Map']


    n_components_used = n_components
    n_components = _validate_dim_n_components(n_components,
                                              method, *adata.X.shape)

    x_emb = None

    clear_x_emb_dependents(adata)
    # Create dimensionality reduction object and find embedding
    if x_emb is None:
        x_emb = wrap('dim_reduction', method)(
            n_components=n_components, **kwargs).get(adata.X)
    _clear_x_emb_dependents(adata)

    if not is_AnnData:
        return x_emb

    # Populate
    adata.obsm['x_emb'] = x_emb
    adata.uns['dim_reduction_info'] = {}
    adata.uns['dim_reduction_info']['method'] = method
    adata.uns['dim_reduction_info']['n_components'] = x_emb.shape[1]
    if n_components_used == 'knee':
        n_components_used = 'Automatic'
    adata.uns['dim_reduction_info']['n_components_used'] = n_components_used
    adata.uns['dim_reduction_info']['kwargs'] = kwargs

    if not inplace:
        return adata

def _clear_x_emb_dependents(adata):
    # Clear labels that used old x_emb
    if 'labels' in adata.obs:
        adata.obs.pop('labels')
        print('Clearing labels...')
    # Clear visualization if it used previous embeddings
    for dim in [2, 3]:
        if f'x_emb_{dim}d' in adata.obsm:
            if adata.uns[f'visualization_info_{dim}d']['used_emb'] == True:
                print('Clearing 2d embeddings...')
                adata.obsm.pop(f'x_emb_{dim}d', None)
                adata.uns.pop(f'visualization_info_{dim}d', None)


#os.chdir('Desktop/cellar')

datapath='../datasets/user_uploaded/HBMP_spleen_merged_cellar.h5ad'
adata = anndata.read_h5ad(datapath)
import timeit
import scanpy as sc
start = timeit.default_timer()
sc.pp.pca(adata,n_comps=40)
end = timeit.default_timer()
time = end - start
