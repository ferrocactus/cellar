import itertools
import json
import re
from ast import literal_eval
from functools import reduce

import anndata
import numpy as np
import pandas as pd


def parse(x):
    if x.size == 0:
        return x
    parsed_x = np.char.split(x.flatten(), sep='.', maxsplit=1)
    parsed_x = np.array([i[0] for i in parsed_x])
    parsed_x = np.char.replace(parsed_x, '-', '')
    parsed_x = np.char.replace(parsed_x, ' ', '')
    # Remove everything that goes inside brackets
    parsed_x = [re.sub(r" ?\([^)]+\)", "", item) for item in parsed_x]
    parsed_x = np.char.upper(parsed_x)
    parsed_x = parsed_x.reshape(x.shape)
    return parsed_x


def _emb_exists_in_adata(adata, method, n_components):
    if 'x_emb' not in adata.obsm:
        return None
    if adata.uns['dim_reduction_info']['method'] != method:
        return None
    # The following will pass also if n_components_used = n_components = knee
    if n_components == adata.uns['dim_reduction_info']['n_components_used']:
        return adata.obsm['x_emb']
    if isinstance(n_components, str):
        return None
    if n_components <= adata.uns['dim_reduction_info']['n_components']:
        return adata.obsm['x_emb'][:, :n_components]
    return None


def _2d_emb_exists_in_adata(adata, method, dim, use_emb):
    if f'x_emb_{dim}d' not in adata.obsm:
        return False
    if adata.uns[f'visualization_info_{dim}d']['method'] != method:
        return False
    if adata.uns[f'visualization_info_{dim}d']['used_emb'] != use_emb:
        return False
    return True
