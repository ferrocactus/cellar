from ._pipeline import reduce_dim
from ._pipeline import cluster
from ._pipeline import reduce_dim_vis
from ._pipeline import name_genes
from ._pipeline import de
from ._pipeline import ss_cluster
from ._pipeline import transfer_labels

from ._plotter import plot

from .utils.read import load_file

from .utils.exceptions import InvalidArgument
from .utils.exceptions import InappropriateArgument
from .utils.exceptions import MethodNotImplementedError
from .utils.exceptions import IncorrectFileFormat

import traceback
import sys

OK = 'good'

units = ['reduce_dim', 'cluster', 'reduce_dim_vis',
         'name_genes', 'de', 'ss_cluster', 'transfer_labels',
         'load_file']


def safe(f, **kwargs):
    try:
        f(**kwargs)
        return OK
    except (InappropriateArgument,
            InvalidArgument,
            MethodNotImplementedError,
            IncorrectFileFormat) as e:
        print(str(e))
        return str(e)
    except Exception as e:
        traceback.print_exc(file=sys.stdout)
        return "An error occurred."
