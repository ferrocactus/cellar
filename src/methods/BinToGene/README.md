<h2>BinToGene</h2>
Library to convert a cell-by-bin matrix to cell-by-gene.
This can be useful, for e.g., when one wishes to use existing gene expression
tools to a sc-ATAC-seq cell-by-bin matrix.</br></br>

The package uses a gencode file (v34) to determine start and end locations
for protein coding genes. This interval is then extended by fixed or gene-length
dependent values to get a new interval I=(start-v1, end+v2). Then the bins in
the provided data matrix are searched to see if they intersect I. An efficient
binary search implementation is used to determine the first and last bin for
which such intersection is found. Then the values of the bins are added or
averaged as specified by the user. The resulting vector is a count vector for
the gene in consideration. Note, the counts here do not represent expression
of the gene, but chromatin accessibility counts in the case of sc-ATAC-seq.
In the end, these vectors are stacked to form a cell-by-gene matrix.
Tools that are used for analyzing gene expression data can be useful in
analyzing the cell-by-gene matrix formed this way.

Example

```python
import numpy as np

x = np.random.randint(0, 2, (200, 4000)) # your cell by bin matrix

bin_names = ['chr1:0-1000', 'chr1:1001-2000', ...] # your list of bin names

btg = BinToGene()  # Use default interval extension parameters
counts, ids = btg.convert(x, bin_names, prefix='chr', delim1=':', delim2='-')
```

Example of a cell-by-gene matrix obtained via UMAP after running some basic preprocessing
<img src="https://github.com/ferrocactus/BinToGene/blob/master/images/example.png" style="zoom:82%;" />
