import json
from abc import abstractmethod
from functools import reduce

import numpy as np
from scipy.stats import hypergeom

from ..log import setup_logger
from ..utils.experiment import parse
from ._unit import Unit


def _get_dict(path):
    """
    Parameters
    __________

    path: string or array
        Path to json file. In case a list of paths is provided instead,
        read them all and merge then into a single dict. Assumes depth two.

    Returns
    _______
    d: dict
        Dictionary containing marker information.
        d = {
            key: {
                subkey: [...],
                ...
            },
            ...
        }
    """
    # TODO straighten up the spaghetti

    if isinstance(path, str):
        with open(path, "r") as f:
            return json.load(f)
    else:
        d = {}
        for path in path:
            with open(path, "r") as f:
                d_part = json.load(f)
                for key in d_part:
                    if key in d:
                        for subkey in d_part[key]:
                            if subkey in d[key]:
                                # to remove duplicates
                                d[key][subkey] = list(set().union(
                                    d[key][subkey], d_part[key][subkey]))
                            else:
                                d[key][subkey] = d_part[key][subkey]
                    else:
                        d[key] = d_part[key]
        return d


class Ide_HyperGeom(Unit):
    """
    Runs hypergeom to find matching populations. Compute for every label
    in x, the pop in pops where x is most likely to have been drawn from.
    It is assumed that the dictionary that is passed has two levels of
    hierarchy of types. First determine the lvl1 type, then the lvl2 subtype.
    """

    def __init__(self, path='markers/cell_type_marker.json', tissue='all'):
        self.logger = setup_logger("HyperGeom")
        self.path = path
        self.tissue = tissue

    def get(self, x):
        """
        Extended keys are: lvl1_type, lvl1_sv, lvl1_intersec, lvl1_total,
                           lvl2_type, lvl2_sv, lvl2_intersec, lvl2_total
                    type (string): identified type
                    sv (float): survival value from Hypergeometric Test
                    intersec (np.ndarray): array of names that overlap
                    total (int): total number of names in dict[type]

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
        x = x.copy()
        lvl2 = _get_dict(self.path)

        # Construct lvl1 dict by merging all lvl2 dicts
        lvl1 = {}
        for pop in lvl2:
            lvl1[pop] = parse(
                np.array(reduce(lambda a, b: a+b, lvl2[pop].values()))
            )

        if self.tissue == 'all':
            self.process_level(x, lvl1, level=1)
            self.process_level(x, lvl2, level=2)
        else:
            self.logger.info(f"Running HyperGeom for {self.tissue} only.")
            self.process_tissue(x, tissue=self.tissue, level_dict=lvl2)

        return x

    def process_level(self, x, level_dict, level):
        for key in x:
            if level > 1 and x[key][f'lvl{level-1}_type'] == 'None':
                tp, sv, intersec, total = "None", 1, np.array([]), 0
            else:
                if level > 1:
                    tp, sv, intersec, total, all_pops = self.find_population(
                        # x[key]['outp_names'],
                        x[key][f'lvl{level-1}_intersec'],
                        level_dict[x[key][f'lvl{level-1}_type']]
                    )
                else:
                    tp, sv, intersec, total, all_pops = self.find_population(
                        x[key]['outp_names'],
                        level_dict
                    )
            x[key][f'lvl{level}_type'] = tp
            x[key][f'lvl{level}_sv'] = sv
            x[key][f'lvl{level}_intersec'] = intersec
            x[key][f'lvl{level}_total'] = total
            x[key][f'lvl{level}_all'] = all_pops
        self.logger.info(f"Finished finding lvl{level} types.")

    def process_tissue(self, x, tissue, level_dict):
        for key in x:
            tp, sv, intersec, total, all_pops = self.find_population(
                x[key]['outp_names'],
                level_dict[tissue]
            )
            x[key]['lvl1_type'] = "User Defined"
            x[key]['lvl1_sv'] = 1
            x[key]['lvl1_intersec'] = np.array([])
            x[key]['lvl1_total'] = 0
            x[key]['lvl1_all'] = {}

            x[key]['lvl2_type'] = tp
            x[key]['lvl2_sv'] = sv
            x[key]['lvl2_intersec'] = intersec
            x[key]['lvl2_total'] = total
            x[key]['lvl2_all'] = all_pops
        self.logger.info("Finished finding lvl2 types.")

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

        survival_values = []
        intersections = []
        lens = []

        rsv, rpop, rk = 2, -1, 0

        for pop in pops:
            n = len(pops[pop])
            intersec = np.intersect1d(x, pops[pop])
            k = len(intersec)
            sv = hypergeom.sf(k-1, M=M, n=n, N=N) if k > 0 else 1

            survival_values.append(sv)
            intersections.append(intersec)
            lens.append(len(pops[pop]))

            if sv <= rsv or (rsv == 2 and k > 0):
                rsv, rpop, rk = sv, pop, k

        all_pops = {'svs': np.array(survival_values),
                    'intersecs': np.array(intersections),
                    'lens': np.array(lens)}

        if rk == 0:  # in case of no intersection, return -1
            return "None", 1, np.array([]), 0, all_pops
        else:
            return rpop, rsv, np.intersect1d(x, pops[rpop]), len(pops[rpop]), all_pops
