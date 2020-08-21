import pandas as pd
import numpy as np
from gtfparse import read_gtf
from tqdm import tqdm

gene_ensembl_path = "genes_ensembl.csv"


def get_gene_df(path=gene_ensembl_path):
    return pd.read_csv(path)


class Gene:
    def __init__(self, start, end, gene_index, seqname):
        self.start = int(start)
        self.end = int(end)
        self.gene_index = gene_index
        if seqname[:3] == 'chr':
            self.seqname = seqname[3:]
        else:
            self.seqname = seqname
        assert self.start < self.end

    def contains(self, val):
        return val >= self.start and val <= self.end

    def precedes(self, val):
        return val > self.end

    def succeeds(self, val):
        return val < self.start

    def __str__(self):
        return f'[{self.start}, {self.end}] at {self.gene_index}'


def parse_bin_id(bin_id):
    interval, chromo = bin_id[::-1].split(":", maxsplit=1)
    start, end = interval[::-1].split("-")
    start, end = int(start), int(end)
    return chromo[::-1], start, end


class Bin:
    def __init__(self, start, end, seqname, index):
        self.start = int(start)
        self.end = int(end)
        self.seqname = seqname
        self.index = index
        assert self.start < self.end

    def contains(self, val):
        return val >= self.start and val <= self.end

    def precedes(self, val):
        return val > self.end

    def succeeds(self, val):
        return val < self.start

    def __str__(self):
        return f'[{self.start}, {self.end}] on {self.seqname} at {self.index}.'


def get_bin_dict(cell_by_bin):
    bin_dict = {}

    for index, bin_id in enumerate(cell_by_bin.var_names):
        chromo, start, end = parse_bin_id(bin_id)
        b = Bin(start, end, chromo, index)
        if chromo not in bin_dict:
            bin_dict[chromo] = []
        bin_dict[chromo].append(b)

    for chromo in bin_dict:
        bin_dict[chromo].sort(key=lambda x: x.start)

    return bin_dict


def get_chromos_dict():
    gene_ensembl = get_gene_df(gene_ensembl_path)
    chromos = np.unique(gene_ensembl['seqname'])
    chromos_dict = {}

    for chromo in tqdm(chromos):
        this_chromo = gene_ensembl[gene_ensembl['seqname'] == chromo]
        chromos_dict[chromo] = [Gene(this_chromo.loc[i]['start'],
                                     this_chromo.loc[i]['end'], i)
                                for i in this_chromo.index]
        chromos_dict[chromo].sort(key=lambda x: x.start)
    return chromos_dict


def bin_search(bin_dict, gene):
    chromo = gene.seqname

    if chromo not in bin_dict:
        return -1, -1
    this_chromo = bin_dict[chromo]

    start_loc, end_loc = -1, -1

    i, j = 0, len(this_chromo) - 1
    while i <= j:
        m = i + (j - i) // 2
        if (m == 0 or this_chromo[m-1].precedes(gene.start))\
                and this_chromo[m].contains(gene.start):
            start_loc = m
            break
        elif this_chromo[m].precedes(gene.start):
            i = m + 1
        else:
            j = m - 1

    i, j = 0, len(this_chromo) - 1
    while i <= j:
        m = i + (j - i) // 2
        if (m == len(this_chromo) - 1 or this_chromo[m+1].succeeds(gene.end))\
                and this_chromo[m].contains(gene.end):
            end_loc = m
            break
        elif this_chromo[m].precedes(gene.end):
            i = m + 1
        else:
            j = m - 1

    if end_loc == -1 and start_loc == -1:
        pass
    elif end_loc == -1:
        last_end = start_loc
        for index in range(start_loc, len(this_chromo)):
            if this_chromo[index].precedes(gene.end):
                last_end = index
            else:
                break
        end_loc = last_end
    elif start_loc == -1:
        first_start = end_loc
        for index in range(end_loc, 0, -1):
            if this_chromo[index].succeeds(gene.start):
                first_start = index
            else:
                break
        start_loc = first_start

    return start_loc, end_loc


def sum_bins(bin_dict, gene, cell_by_bin):
    start_loc, end_loc = bin_search(bin_dict, gene)
    if start_loc == end_loc and start_loc == -1:
        return None
    assert end_loc >= start_loc

    this_chromo = bin_dict[gene.seqname]

    gene_expr = np.zeros(cell_by_bin.X[:, 0].shape)

    for i in range(start_loc, end_loc+1):
        gene_expr += cell_by_bin.X[:, this_chromo[i].index].todense()

    #gene_expr /= (end_loc + 1 - start_loc)

    return np.squeeze(gene_expr)


def get_expr_mat(cell_by_bin):
    bin_dict = get_bin_dict(cell_by_bin)
    gene_df = get_gene_df()

    df = pd.DataFrame()

    for i, gene in tqdm(enumerate(gene_df.index)):
        g = Gene(gene_df.loc[gene].start, gene_df.loc[gene].end,
                 i, gene_df.loc[gene].seqname[3:])
        s = sum_bins(bin_dict, g, cell_by_bin)
        if s is not None:
            df[gene_df.loc[gene].gene_name] = s

    return df
