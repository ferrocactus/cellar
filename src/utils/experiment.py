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
