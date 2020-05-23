from abc import abstractmethod

import numpy as np
import pandas as pd

from ..log import setup_logger
from ..utils.experiment import parse
from ._unit import Unit


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
    
    def check_name(self, names):
        """
        Check if the name exists
        """
        if (type(names)!="numpy.ndarray"):
            names=np.array([names])
        
        gene_dict = pd.read_csv(self.path, squeeze=True)
        col1, col2 = gene_dict.columns

        # Revert the order of columns
        gene_dict = gene_dict[gene_dict.columns[::-1]
                              ].set_index(col2)[col1].to_dict()
        parsed_names = parse(names)

        if parsed_names.size == 0:
            return np.array([])

        # Leave unchanged if not found
        ids = np.char.upper([gene_dict.get(i, i) for i in parsed_names])
        
        if ids==np.char.upper(names):
            return 0 # not found
        else:
            return 1 # found
