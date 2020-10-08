import unittest

from src.bin_to_gene import BinToGene
from src.bin import Bin

class bin_test(unittest.TestCase):
    def test1(self):
        b = BinToGene()
        bin_list = [Bin("s", 0, 10, 0),
                    Bin("s", 11, 20, 0),
                    Bin("s", 21, 30, 0),
                    Bin("s", 31, 40, 0),
                    Bin("s", 41, 50, 0)]

        result = b.binary_search(bin_list, 2, 21, True)
        assert result == 0
        result = b.binary_search(bin_list, 12, 22, True)
        assert result == 1
        result = b.binary_search(bin_list, 22, 32, True)
        assert result == 2
        result = b.binary_search(bin_list, 32, 42, True)
        assert result == 3
        result = b.binary_search(bin_list, 42, 52, True)
        assert result == 4
        result = b.binary_search(bin_list, 12, 71, True)
        assert result == 1
        result = b.binary_search(bin_list, -10, 21, True)
        assert result == 0
        result = b.binary_search(bin_list, -10, 71, True)
        assert result == 0
        result = b.binary_search(bin_list, -20, -10, True)
        assert result == -1
        result = b.binary_search(bin_list, 61, 81, True)
        assert result == -1

    def test2(self):
        b = BinToGene()
        bin_list = [Bin("s", 0, 10, 0),
                    Bin("s", 11, 20, 0),
                    Bin("s", 21, 30, 0),
                    Bin("s", 31, 40, 0)]

        result = b.binary_search(bin_list, 0, 10, True)
        assert result == 0
        result = b.binary_search(bin_list, 11, 20, True)
        assert result == 1
        result = b.binary_search(bin_list, 21, 31, True)
        assert result == 2
        result = b.binary_search(bin_list, 31, 41, True)
        assert result == 3
        result = b.binary_search(bin_list, 12, 71, True)
        assert result == 1
        result = b.binary_search(bin_list, -10, 21, True)
        assert result == 0
        result = b.binary_search(bin_list, -10, 71, True)
        assert result == 0
        result = b.binary_search(bin_list, -20, -10, True)
        assert result == -1
        result = b.binary_search(bin_list, 61, 81, True)
        assert result == -1

    def test3(self):
        b = BinToGene()
        bin_list = [Bin("s", 0, 10, 0),
                    Bin("s", 11, 20, 0),
                    Bin("s", 21, 30, 0),
                    Bin("s", 31, 40, 0),
                    Bin("s", 41, 50, 0)]

        result = b.binary_search(bin_list, 2, 21, False)
        assert result == 2
        result = b.binary_search(bin_list, 12, 22, False)
        assert result == 2
        result = b.binary_search(bin_list, 22, 32, False)
        assert result == 3
        result = b.binary_search(bin_list, 32, 42, False)
        assert result == 4
        result = b.binary_search(bin_list, 42, 52, False)
        assert result == 4
        result = b.binary_search(bin_list, 12, 71, False)
        assert result == 4
        result = b.binary_search(bin_list, -10, 21, False)
        assert result == 2
        result = b.binary_search(bin_list, -10, 71, False)
        assert result == 4
        result = b.binary_search(bin_list, -10, 1, False)
        assert result == 0
        result = b.binary_search(bin_list, -20, -10, False)
        assert result == -1
        result = b.binary_search(bin_list, 61, 81, False)
        assert result == -1

    def test4(self):
        b = BinToGene()
        bin_list = [Bin("s", 0, 10, 0),
                    Bin("s", 21, 30, 0),
                    Bin("s", 31, 40, 0),
                    Bin("s", 41, 50, 0),
                    Bin("s", 71, 90, 0)]

        result = b.binary_search(bin_list, 2, 13, True)
        assert result == 0
        result = b.binary_search(bin_list, 2, 13, False)
        assert result == 0

        result = b.binary_search(bin_list, 12, 22, True)
        assert result == 1
        result = b.binary_search(bin_list, 12, 22, False)
        assert result == 1

        result = b.binary_search(bin_list, 41, 62, True)
        assert result == 3
        result = b.binary_search(bin_list, 41, 62, False)
        assert result == 3
        result = b.binary_search(bin_list, 41, 92, False)
        assert result == 4

        result = b.binary_search(bin_list, 11, 20, True)
        assert result == -1
        result = b.binary_search(bin_list, 11, 20, False)
        assert result == -1
        result = b.binary_search(bin_list, 51, 60, True)
        assert result == -1
        result = b.binary_search(bin_list, 51, 60, False)
        assert result == -1

if __name__ == '__main__':
    unittest.main()

