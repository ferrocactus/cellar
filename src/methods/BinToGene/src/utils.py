import os
from os import path

import pandas as pd


def load_gene_file(gencode_path, logger=None):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    paths = [gencode_path,
             path.join("resources", gencode_path),
             path.join("src/resources", gencode_path),
             path.join(dir_path, 'resources', gencode_path),
             path.join(dir_path, '../resources', gencode_path)]

    for g_path in paths:
        if path.exists(g_path):
            path_to_use = g_path
            if logger is not None:
                logger.info("Using the gene file " + path_to_use)
            break
    else:
        raise FileNotFoundError("No gencode file found in path.")

    gencode = pd.read_csv(path_to_use)

    gencode.drop([
        'feature', 'ccdsid', 'protein_id', 'ont', 'exon_id', 'exon_number',
        'havana_transcript', 'tag', 'transcript_support_level',
        'transcript_name', 'transcript_type', 'transcript_id', 'score', 'frame'],
        axis=1, inplace=True, errors='ignore')

    return gencode


def str_to_int(extend):
    v = ValueError("Invalid extension value found.")

    if extend[-1] != 'x':
        raise v
    extend = extend[:-1]
    try:
        extend = int(extend)
    except:
        raise v

    return extend


def extend_interval(start, end, extend=None, max_extend=None, stream=None,
                    op_extend=None, max_op_extend=None):
    """
    If stream = "+"
        start = start - min(extend, max_extend)
        end = end + min(op_extend, max_op_extend)

    This functions resolves all cases when either is None
    """
    if extend is None:
        return start, end
    if max_extend is None:
        max_extend = extend

    gene_len = int(end - start)
    assert gene_len > 0

    if isinstance(extend, str):
        extend = str_to_int(extend) * gene_len
    if isinstance(max_extend, str):
        max_extend = str_to_int(max_extend) * gene_len

    assert extend > 0
    assert max_extend > 0

    if stream is None:
        start = int(start - min(extend, max_extend))
        end = int(end + min(extend, max_extend))
        return start, end

    if op_extend is None:
        raise ValueError("Please specify opposite strand extension "
                         "or set stream to None.")
    if max_op_extend is None:
        max_op_extend = op_extend

    if isinstance(op_extend, str):
        op_extend = str_to_int(op_extend) * gene_len
    if isinstance(max_op_extend, str):
        max_op_extend = str_to_int(max_op_extend) * gene_len

    assert op_extend > 0
    assert max_op_extend > 0

    if stream == '+':  # upstream
        start = int(start - min(extend, max_extend))
        end = int(end + min(op_extend, max_op_extend))
        return start, end
    else:  # downstream
        start = int(start - min(op_extend, max_op_extend))
        end = int(end + min(extend, max_extend))
        return start, end


def validate_n_jobs(n_jobs):
    if n_jobs is None:
        return 1
    if n_jobs == 0:
        return 1
    return int(n_jobs)
