from ._pipeline import preprocess
from ._pipeline import reduce_dim
from ._pipeline import cluster
from ._pipeline import reduce_dim_vis
from ._pipeline import name_genes
from ._pipeline import de
from ._pipeline import ss_cluster
from ._pipeline import transfer_labels


from ._plotter import plot

from .utils.tools import parse
from .utils.tools import store_subset
from .utils.tools import store_labels
from .utils.tools import match_labels
from .utils.tools import merge_cluster_names
from .utils.tools import update_subset_label
from .utils.tools import populate_subsets
from .utils.tools import merge_clusters
from .utils.tools import get_neighbors
from .utils.tools import uncertainty


from .utils.read import load_file

from .utils.codex_tile_generator import generate_tile

from .utils.exceptions import InvalidArgument
from .utils.exceptions import InappropriateArgument
from .utils.exceptions import MethodNotImplementedError

import traceback
import sys


OK = 'good'


def safe(f, **kwargs):
    try:
        if f == preprocess:
            adata = f(**kwargs)
            return adata
        else:
            f(**kwargs)
            return OK
    except (InappropriateArgument,
            InvalidArgument,
            MethodNotImplementedError) as e:
        print(str(e))
        return str(e)
    except Exception as e:
        traceback.print_exc(file=sys.stdout)
        return "Something went wrong."
