import unittest

from src.bin import Bin

class bin_test(unittest.TestCase):
    def test1(self):
        b = Bin("s", 20, 60, 0)
        assert b.intersects(-10, 30)
        assert b.intersects(30, 80)
        assert b.intersects(30, 50)
        assert b.intersects(10, 90)

        assert not b.intersects(-10, 0)
        assert not b.intersects(80, 100)

    def test2(self):
        b = Bin("s", 20, 60, 0)
        assert b.precedes(70, 90)
        assert not b.precedes(50, 90)
        assert not b.precedes(30, 40)
        assert not b.precedes(0, 30)
        assert not b.precedes(0, 10)

    def test3(self):
        b = Bin("s", 20, 60, 0)
        assert not b.suceeds(70, 90)
        assert not b.suceeds(50, 90)
        assert not b.suceeds(30, 40)
        assert not b.suceeds(0, 30)
        assert b.suceeds(0, 10)

if __name__ == '__main__':
    unittest.main()

