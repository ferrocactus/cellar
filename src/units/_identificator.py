from ._unit import Unit
from src.utils.utils_experiment import parse

from abc import abstractmethod
import json
from functools import reduce
import numpy as np
from scipy.stats import hypergeom

PATH = "markers/cell_type_marker.json"


class Ide(Unit):
    """
    Base class for gene identification methods.
    """
    def __init__(self, verbose=False, **kwargs):
        """
        Args:
            verbose (bool): Printing flag.
            **kwargs: Argument dict.
        """
        super().__init__(verbose, **kwargs)
        self.name = 'Ide'
        self.path = kwargs.get('path', PATH)

    @abstractmethod
    def get(self, x):
        """
        Returns the types of cells in x.

        Args:
            x (dict): x = {
                label_1: {
                    outp_names: [name_1, ...],
                    ...
                },
                ...
            }
        Returns:
            (dict): Extends x with new keys (returns copy).
        """
        pass


class Ide_HyperGeom(Ide):
    """
    Runs hypergeom to find matching populations. Compute for every label
    in x, the pop in pops where x is most likely to have been drawn from.
    It is assumed that the dictionary that is passed has two levels of
    hierarchy of types. First determine the lvl1 type, then the lvl2 subtype.
    """
    def __init__(self, verbose=False, **kwargs):
        super().__init__(verbose, **kwargs)

    def get(self, x):
        """
        Extended keys are: lvl1_type, lvl1_sv, lvl1_intersec, lvl1_total,
                           lvl2_type, lvl2_sv, lvl2_intersec, lvl2_total
                    type (string): identified type
                    sv (float): survival value from Hypergeometric Test
                    intersec (np.ndarray): array of names that overlap
                    total (int): total number of names in dict[type]
        """
        x = x.copy()
        lvl2 = self.get_dict() # Assumed to have two nested levels

        # Construct lvl1 dict by merging all lvl2 dicts
        lvl1 = {}
        for pop in lvl2:
            lvl1[pop] = parse(
                np.array(reduce(lambda a, b: a+b, lvl2[pop].values()))
            )

        # Level 1 in the hierarchy identification loop
        for key in x:
            tp, sv, intersec, total = self.find_population(
                x[key]['outp_names'],
                lvl1
            )
            x[key]['lvl1_type'] = tp
            x[key]['lvl1_sv'] = sv
            x[key]['lvl1_intersec'] = intersec
            x[key]['lvl1_total'] = total
        self.vprint("Finished finding lvl1 types.")

        # Level 2 in the hierarchy identification loop
        for key in x:
            if x[key]['lvl1_type'] == 'None':
                tp, sv, intersec, total = "None", 1, np.array([]), 0
            else:
                tp, sv, intersec, total = self.find_population(
                    x[key]['outp_names'],
                    lvl2[x[key]['lvl1_type']]
                )
            x[key]['lvl2_type'] = tp
            x[key]['lvl2_sv'] = sv
            x[key]['lvl2_intersec'] = intersec
            x[key]['lvl2_total'] = total
        self.vprint("Finished finding lvl2 types.")

        return x

    def get_dict(self):
        """
        Reads json file and converts to dict. In case a list of paths
        is provided instead, read them all and merge then into a single
        dict.

        Returns: dict.
        """
        if isinstance(self.path, str):
            with open(self.path, "r") as f:
                return json.load(f)
        else:
            d = {}
            for path in self.path:
                with open(path, "r") as f:
                    d = {**d, **json.load(f)}
            return d

    def find_population(self, x, pops):
        """
        See find_populations. Assumes x is a single list.

        Args:
            x (np.ndarray): 1D list of names.
            pops (dict): Dictionary of populations: pops = {
                type: [name_1, name_2, ...],
                ...
            }
        Returns:
            (string): population name
            (float): survival value
            (np.ndarray): common names
            (int): total number of names in matched population
        """
        M = sum([len(pops[pop]) for pop in pops])
        N = len(x)

        rsv, rpop, rk = 2, -1, 0

        for pop in pops:
            n = len(pops[pop])
            k = len(np.intersect1d(x, pops[pop]))
            sv = hypergeom.sf(k-1, M=M, n=n, N=N) if k > 0 else 1
            if sv < rsv:
                rsv, rpop, rk = sv, pop, k
        if rk == 0: # in case of no intersection, return -1
            return "None", 1, np.array([]), 0
        else:
            return rpop, rsv, np.intersect1d(x, pops[rpop]), len(pops[rpop])