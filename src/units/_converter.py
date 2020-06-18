from abc import abstractmethod

import numpy as np
import pandas as pd

from ..log import setup_logger
from ..utils.tools import parse
from ._unit import Unit


path = 'markers/gene_id_name.csv'


def convert(str_list):
    convention = find_naming_convention(str_list)
    d = {}
    if convention == 'ids':
        d['names'], d['ids'] = id_to_name(str_list)
    else:
        d['ids'], d['names'] = name_to_id(str_list)
    return d


def find_naming_convention(str_list):
    if str_list[0][:3] == 'ENS':
        return 'ids'
    else:
        return 'names'


def id_to_name(ids):
    """
    Given an array of gene ids, convert them to gene names.
    Ids: ENSG* format.
    Args:
        ids (np.ndarray): Array of strings.
    """
    gene_dict = pd.read_csv(path, index_col=0, squeeze=True).to_dict()
    parsed_ids = parse(ids)

    if parsed_ids.size == 0:
        return np.array([])

    # Leave unchanged if not found
    return np.char.upper([gene_dict.get(i, i) for i in parsed_ids]), parsed_ids


def name_to_id(names):
    """
    Does the opposite of id_to_name.
    """
    gene_dict = pd.read_csv(path, squeeze=True)
    col1, col2 = gene_dict.columns

    # Revert the order of columns
    gene_dict = gene_dict[gene_dict.columns[::-1]
                          ].set_index(col2)[col1].to_dict()
    parsed_names = parse(names)

    if parsed_names.size == 0:
        return np.array([])

    # Leave unchanged if not found
    return np.char.upper([gene_dict.get(i, i) for i in parsed_names]), parsed_names


class Con(Unit):
    """
    Base class for converting marker names.
    """

    def __init__(self, convention='id-to-name',
                 path='markers/gene_id_name.csv'):
        """
        Parameters
        __________

        convention: string, default 'id-to-name'
            Convention to use for converting names. Currently
            accepts 'id-to-name', 'name-to-id', or None.
            If set to None, will not convert names.

        path: string or list of strings
            Path to markers file. Note that the columns should be
            ordered according to convention. If list of strings
            provided instead, will merge the marker files in a
            single dictionary. Care should be taken that the files
            have the same format and hierarchy.

        """
        self.logger = setup_logger('Converter')
        self.convention = convention
        self.path = path

    def get(self, x):
        """
        Returns the names of the markers according to convention.

        Parameters
        __________
        x: array: Names of the genes.

        Returns
        _______
        array: Converted names.

        """
        if self.convention is None:
            return x
        elif self.convention == 'id-to-name':
            return self.id_to_name(x)
        elif self.convention == 'name-to-id':
            return self.name_to_id(x)
        else:
            raise NotImplementedError("Convention not found.")

    def id_to_name(self, ids):
        """
        Given an array of gene ids, convert them to gene names.
        Ids: ENSG* format. Path is assumed to be well formatted.
        Args:
            ids (np.ndarray): List of strings.
        """
        gene_dict = pd.read_csv(self.path, index_col=0, squeeze=True).to_dict()
        parsed_ids = parse(ids)

        if parsed_ids.size == 0:
            return np.array([])

        # Leave unchanged if not found
        return np.char.upper([gene_dict.get(i, i) for i in parsed_ids])

    def name_to_id(self, names):
        """
        Does the opposite of id_to_name.
        """
        gene_dict = pd.read_csv(self.path, squeeze=True)
        col1, col2 = gene_dict.columns

        # Revert the order of columns
        gene_dict = gene_dict[gene_dict.columns[::-1]
                              ].set_index(col2)[col1].to_dict()
        parsed_names = parse(names)

        if parsed_names.size == 0:
            return np.array([])

        # Leave unchanged if not found
        return np.char.upper([gene_dict.get(i, i) for i in parsed_names])
