import unittest
import numpy as np
from numpy.testing import assert_array_almost_equal as aae
from numpy.testing import assert_array_equal as ae

from src.bin_to_gene import BinToGene as BTG


class bin_test(unittest.TestCase):
    def test1(self):
        btg = BTG(gencode_path='src/tests/gencode_example.csv', extend=None)

        bin_names = np.array([
            'chr1:0-100',  # 0
            'chr1:101-200',  # 1
            'chr1:201-300',  # 2

            'chr2:0-100',  # 3
            'chr2:101-200',  # 4
            'chr2:201-300',  # 5
            'chr2:301-400',  # 6

            'chr3:0-100',  # 7
            'chr3:101-200',  # 8
            'chr3:201-300',  # 9
            'chr3:301-400',  # 10
            'chr3:401-500',  # 11
            'chr3:501-600'  # 12
        ])

        x = np.random.random((1000, len(bin_names)))

        counts, ids = btg.convert(x, bin_names, prefix='chr')
        ae(ids, ["gene" + str(i) for i in range(1, 8)])

        ground_counts = np.zeros((1000, 7))
        ground_counts[:, 0] += x[:, 0]
        ground_counts[:, 1] += x[:, 1]
        ground_counts[:, 2] += x[:, 4] + x[:, 5]
        ground_counts[:, 3] += x[:, 3] + x[:, 4] + x[:, 5]
        ground_counts[:, 4] += x[:, 8] + x[:, 9]
        ground_counts[:, 5] += x[:, 9:13].sum(axis=1)
        ground_counts[:, 6] += x[:, 9:13].sum(axis=1)

        aae(ground_counts, counts)

    def test2(self):
        btg = BTG(
            gencode_path='src/tests/gencode_example.csv',
            operation='mean', extend=None)

        bin_names = np.array([
            'chr1:0-100',  # 0
            'chr1:101-200',  # 1
            'chr1:201-300',  # 2

            'chr2:0-100',  # 3
            'chr2:101-200',  # 4
            'chr2:201-300',  # 5
            'chr2:301-400',  # 6

            'chr3:0-100',  # 7
            'chr3:101-200',  # 8
            'chr3:201-300',  # 9
            'chr3:301-400',  # 10
            'chr3:401-500',  # 11
            'chr3:501-600'  # 12
        ])

        x = np.random.random((1000, len(bin_names)))

        counts, ids = btg.convert(x, bin_names, prefix='chr')
        ae(ids, ["gene" + str(i) for i in range(1, 8)])

        ground_counts = np.zeros((1000, 7))
        ground_counts[:, 0] += x[:, 0]
        ground_counts[:, 1] += x[:, 1]
        ground_counts[:, 2] += (x[:, 4] + x[:, 5]) / 2
        ground_counts[:, 3] += (x[:, 3] + x[:, 4] + x[:, 5]) / 3
        ground_counts[:, 4] += (x[:, 8] + x[:, 9]) / 2
        ground_counts[:, 5] += x[:, 9:13].mean(axis=1)
        ground_counts[:, 6] += x[:, 9:13].mean(axis=1)

        aae(ground_counts, counts)

    def test3(self):
        btg = BTG(
            gencode_path='src/tests/gencode_example.csv', extend=None)

        bin_names = np.array([
            'chr1:0-100',  # 0 -> 0
            'chr2:201-300',  # 5 -> 1
            'chr2:301-400',  # 6 -> 2
            'chr1:101-200',  # 1 -> 3
            'chr3:0-100',  # 7 -> 4
            'chr3:101-200',  # 8 -> 5
            'chr2:0-100',  # 3 -> 6
            'chr2:101-200',  # 4 -> 7
            'chr3:201-300',  # 9 -> 8
            'chr3:301-400',  # 10 -> 9
            'chr1:201-300',  # 2 -> 10
            'chr3:401-500',  # 11 -> 11
            'chr3:501-600'  # 12 -> 12
        ])

        x = np.random.random((1000, len(bin_names)))

        counts, ids = btg.convert(x, bin_names, prefix='chr')
        ae(ids, ["gene" + str(i) for i in range(1, 8)])

        ground_counts = np.zeros((1000, 7))
        ground_counts[:, 0] += x[:, 0]
        ground_counts[:, 1] += x[:, 3]
        ground_counts[:, 2] += x[:, 7] + x[:, 1]
        ground_counts[:, 3] += x[:, 6] + x[:, 7] + x[:, 1]
        ground_counts[:, 4] += x[:, 5] + x[:, 8]
        ground_counts[:, 5] += x[:, 8] + x[:, 9] + x[:, 11] + x[:, 12]
        ground_counts[:, 6] += x[:, 8] + x[:, 9] + x[:, 11] + x[:, 12]

        aae(ground_counts, counts)

    def test4(self):
        btg = BTG(
            gencode_path='src/tests/gencode_example.csv', extend=None)

        bin_names = np.array([
            'chr1:0-100',  # 0
            'chr3:501-600',  # 1
            'chr2:201-300',  # 2
            'chr2:301-400',  # 3
            'chr1:101-200',  # 4
            'chr2:601-700',  # 5
            'chr3:700-900'  # 6
        ])

        x = np.random.random((1000, len(bin_names)))

        counts, ids = btg.convert(x, bin_names, prefix='chr')
        ae(ids, ["gene1", "gene2", "gene3", "gene4", "gene6", "gene7"])

        ground_counts = np.zeros((1000, 6))
        ground_counts[:, 0] += x[:, 0]
        ground_counts[:, 1] += x[:, 4]
        ground_counts[:, 2] += x[:, 2]
        ground_counts[:, 3] += x[:, 2]
        ground_counts[:, 4] += x[:, 1]
        ground_counts[:, 5] += x[:, 1] + x[:, 6]

        aae(ground_counts, counts)

    def test5(self):
        btg = BTG(
            gencode_path='src/tests/gencode_example.csv',
            extend=None, n_jobs=4)

        bin_names = np.array([
            'chr1:0-100',  # 0
            'chr1:101-200',  # 1
            'chr1:201-300',  # 2

            'chr2:0-100',  # 3
            'chr2:101-200',  # 4
            'chr2:201-300',  # 5
            'chr2:301-400',  # 6

            'chr3:0-100',  # 7
            'chr3:101-200',  # 8
            'chr3:201-300',  # 9
            'chr3:301-400',  # 10
            'chr3:401-500',  # 11
            'chr3:501-600'  # 12
        ])

        x = np.random.random((1000, len(bin_names)))

        counts, ids = btg.convert(x, bin_names, prefix='chr')
        ae(ids, ["gene" + str(i) for i in range(1, 8)])

        ground_counts = np.zeros((1000, 7))
        ground_counts[:, 0] += x[:, 0]
        ground_counts[:, 1] += x[:, 1]
        ground_counts[:, 2] += x[:, 4] + x[:, 5]
        ground_counts[:, 3] += x[:, 3] + x[:, 4] + x[:, 5]
        ground_counts[:, 4] += x[:, 8] + x[:, 9]
        ground_counts[:, 5] += x[:, 9:13].sum(axis=1)
        ground_counts[:, 6] += x[:, 9:13].sum(axis=1)

        aae(ground_counts, counts)

    def test6(self):
        btg = BTG(
            gencode_path='src/tests/gencode_example.csv',
            operation='mean', extend=None, n_jobs=3)

        bin_names = np.array([
            'chr1:0-100',  # 0
            'chr1:101-200',  # 1
            'chr1:201-300',  # 2

            'chr2:0-100',  # 3
            'chr2:101-200',  # 4
            'chr2:201-300',  # 5
            'chr2:301-400',  # 6

            'chr3:0-100',  # 7
            'chr3:101-200',  # 8
            'chr3:201-300',  # 9
            'chr3:301-400',  # 10
            'chr3:401-500',  # 11
            'chr3:501-600'  # 12
        ])

        x = np.random.random((1000, len(bin_names)))

        counts, ids = btg.convert(x, bin_names, prefix='chr')
        ae(ids, ["gene" + str(i) for i in range(1, 8)])

        ground_counts = np.zeros((1000, 7))
        ground_counts[:, 0] += x[:, 0]
        ground_counts[:, 1] += x[:, 1]
        ground_counts[:, 2] += (x[:, 4] + x[:, 5]) / 2
        ground_counts[:, 3] += (x[:, 3] + x[:, 4] + x[:, 5]) / 3
        ground_counts[:, 4] += (x[:, 8] + x[:, 9]) / 2
        ground_counts[:, 5] += x[:, 9:13].mean(axis=1)
        ground_counts[:, 6] += x[:, 9:13].mean(axis=1)

        aae(ground_counts, counts)

    def test7(self):
        btg = BTG(
            gencode_path='src/tests/gencode_example.csv',
            extend=None, n_jobs=4)

        bin_names = np.array([
            'chr1:0-100',  # 0 -> 0
            'chr2:201-300',  # 5 -> 1
            'chr2:301-400',  # 6 -> 2
            'chr1:101-200',  # 1 -> 3
            'chr3:0-100',  # 7 -> 4
            'chr3:101-200',  # 8 -> 5
            'chr2:0-100',  # 3 -> 6
            'chr2:101-200',  # 4 -> 7
            'chr3:201-300',  # 9 -> 8
            'chr3:301-400',  # 10 -> 9
            'chr1:201-300',  # 2 -> 10
            'chr3:401-500',  # 11 -> 11
            'chr3:501-600'  # 12 -> 12
        ])

        x = np.random.random((1000, len(bin_names)))

        counts, ids = btg.convert(x, bin_names, prefix='chr')
        ae(ids, ["gene" + str(i) for i in range(1, 8)])

        ground_counts = np.zeros((1000, 7))
        ground_counts[:, 0] += x[:, 0]
        ground_counts[:, 1] += x[:, 3]
        ground_counts[:, 2] += x[:, 7] + x[:, 1]
        ground_counts[:, 3] += x[:, 6] + x[:, 7] + x[:, 1]
        ground_counts[:, 4] += x[:, 5] + x[:, 8]
        ground_counts[:, 5] += x[:, 8] + x[:, 9] + x[:, 11] + x[:, 12]
        ground_counts[:, 6] += x[:, 8] + x[:, 9] + x[:, 11] + x[:, 12]

        aae(ground_counts, counts)

    def test8(self):
        btg = BTG(
            gencode_path='src/tests/gencode_example.csv',
            extend=None, n_jobs=4)

        bin_names = np.array([
            'chr1:0-100',  # 0
            'chr3:501-600',  # 1
            'chr2:201-300',  # 2
            'chr2:301-400',  # 3
            'chr1:101-200',  # 4
            'chr2:601-700',  # 5
            'chr3:700-900'  # 6
        ])

        x = np.random.random((1000, len(bin_names)))

        counts, ids = btg.convert(x, bin_names, prefix='chr')
        ae(ids, ["gene1", "gene2", "gene3", "gene4", "gene6", "gene7"])

        ground_counts = np.zeros((1000, 6))
        ground_counts[:, 0] += x[:, 0]
        ground_counts[:, 1] += x[:, 4]
        ground_counts[:, 2] += x[:, 2]
        ground_counts[:, 3] += x[:, 2]
        ground_counts[:, 4] += x[:, 1]
        ground_counts[:, 5] += x[:, 1] + x[:, 6]

        aae(ground_counts, counts)

    def test9(self):
        btg = BTG(
            gencode_path='src/tests/gencode_example.csv',
            extend=None, n_jobs=4)

        bin_names = np.array([
            'chr1:0-100',  # 0
        ])

        x = np.random.random((1000, len(bin_names)))

        counts, ids = btg.convert(x, bin_names, prefix='chr')
        ae(ids, ["gene1"])

        aae(x[:, 0].reshape(-1, 1), counts)

    def test10(self):
        btg = BTG(
            gencode_path='src/tests/gencode_example.csv',
            extend=None, n_jobs=4)

        bin_names = np.array([
            'chr1:1000-2000',  # 0
        ])

        x = np.random.random((1000, len(bin_names)))

        counts, ids = btg.convert(x, bin_names, prefix='chr')
        ae(ids, [])

        aae([], counts)


if __name__ == '__main__':
    unittest.main()
