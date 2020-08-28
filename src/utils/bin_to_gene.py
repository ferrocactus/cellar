import numpy as np
import pandas as pd
import anndata
import gtfparse
from matplotlib import pyplot as plt
from tqdm import tqdm

#g = gtfparse.read_gtf('gencode.v34.chr_patch_hapl_scaff.annotation.gtf')
genes = pd.read_csv('genes_only.csv', index_col=0)
genes.drop([
    'feature', 'ccdsid', 'protein_id', 'ont', 'exon_id', 'exon_number',
    'havana_transcript', 'tag', 'transcript_support_level',
    'transcript_name', 'transcript_type', 'transcript_id', 'score', 'frame'],
    axis=1, inplace=True)

a = anndata.read_h5ad('cell_by_bin.h5ad')


class Bin:
    def __init__(self, index, start, end, seqname):
        self.index = index
        self.start = start
        self.end = end
        assert start < end
        self.seqname = seqname

    def __repr__(self):
        return f'Bin({self.index}, {self.start}, {self.end}, {self.seqname})'

    def intersects(self, start, end):
        assert start < end, "Invalid interval"
        return not ((end < self.start) or (self.end < start))

    def precedes(self, start, end):
        return self.end < start

    def suceeds(self, start, end):
        return end < self.start


bins = {}

chromos = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13',
           '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

for i, bin_name in enumerate(a.var_names):
    seqname, interval = bin_name.split(':')
    if seqname not in chromos:
        continue
    seqname = 'chr' + seqname
    start, end = interval.split('-')
    start, end = int(start), int(end)
    b = Bin(index=i, start=start, end=end, seqname=seqname)
    if seqname not in bins:
        bins[seqname] = []
    bins[seqname].append(b)

for seqname in bins:
    bins[seqname].sort(key=lambda x: x.start)
    bins[seqname] = np.array(bins[seqname])


counts = pd.DataFrame()


def bin_search(gene, bins):
    # bin bin search; Get it? Get it?
    bin_list = bins[gene.seqname]
    assert gene.start < gene.end, "Wrong gene interval"

    if gene.strand == '+':  # upstream gene
        s, e = int(gene.start - 5e4), int(gene.end + 1e4)
    elif gene.strand == '-':  # downstream gene
        s, e = int(gene.start - 1e4), int(gene.end + 5e4)

    # find starting bin
    i, j, first = 0, len(bin_list) - 1, -1
    while i <= j:
        m = (j + i) // 2
        if bin_list[m].suceeds(s, e):
            j = m - 1
        elif bin_list[m].precedes(s, e):
            i = m + 1
        else:
            first = m
            j = m - 1

    # find ending bin
    i, j, last = 0, len(bin_list) - 1, -1
    while i <= j:
        m = (j + i) // 2
        if bin_list[m].suceeds(s, e):
            j = m - 1
        elif bin_list[m].precedes(s, e):
            i = m + 1
        else:
            last = m
            i = m + 1

    assert last >= first, f"last less than first, {last} < {first}"

    if last != -1 and first != -1:  # gene is expressed in some bin
        counts[gene.gene_name] = np.sum(
            a.X[:, first:(last+1)].toarray(), axis=1)


for gene_index in tqdm(genes.index):
    gene = genes.loc[gene_index]
    bin_search(gene, bins)
