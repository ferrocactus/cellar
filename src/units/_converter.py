from ._unit import Unit
from src.utils.utils_experiment import parse

from abc import abstractmethod
import pandas as pd
import numpy as np

CONVENTION = 'id-to-name'
PATH = 'markers/gene_id_name.csv'


class Con(Unit):
    """
    Base class for converting marker names.
    """
    def __init__(self, verbose=False, **kwargs):
        """
        Args:
            verbose (bool): Printing flag.
            **kwargs: Argument dict.
        """
        super().__init__(verbose, **kwargs)
        self.name = 'Con'
        self.convention = kwargs.get('convention', CONVENTION)
        self.path = kwargs.get('path', PATH)

    def get(self, x):
        """
        Returns the names of the markers according to convention.

        Args:
            x (np.ndarray): Marker names.
        Returns:
            (np.ndarray): Names.
        """
        if self.convention == 'id-to-name':
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
        # Leave unchanged if not found
        return np.char.upper([gene_dict.get(i, i) for i in parsed_ids])

    def name_to_id(self, names):
        """
        Does the opposite of id_to_name.
        """
        gene_dict = pd.read_csv(self.path, squeeze=True)
        col1, col2 = gene_dict.columns
        # Revert the order of columns
        gene_dict = gene_dict[gene_dict.columns[::-1]].set_index(col2)[col1].to_dict()
        parsed_names = parse(names)
        # Leave unchanged if not found
        return np.char.upper([gene_dict.get(i, i) for i in parsed_names])