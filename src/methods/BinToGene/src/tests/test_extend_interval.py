import unittest

from src.utils import extend_interval

class extend_interval_test(unittest.TestCase):
    def test1(self):
        start, end = 20, 80
        s, e = extend_interval(start, end)

        assert s == start
        assert e == end

    def test2(self):
        start, end = 20, 80
        s, e = extend_interval(start, end, 20)

        assert s == 0
        assert e == 100

    def test3(self):
        start, end = 20, 80
        s, e = extend_interval(start, end, '2x')

        assert s == -100
        assert e == 200

    def test4(self):
        start, end = 20, 80
        s, e = extend_interval(start, end, '2x', max_extend=40)

        assert s == -20
        assert e == 120

    def test5(self):
        start, end = 20, 80
        s, e = extend_interval(start, end, '2x', max_extend=150)

        assert s == -100
        assert e == 200

    def test6(self):
        start, end = 20, 80
        s, e = extend_interval(start, end, '2x', max_extend=40,
                               stream='+', op_extend=300)

        assert s == -20
        assert e == 380

    def test7(self):
        start, end = 20, 80
        s, e = extend_interval(start, end, '2x', max_extend=40,
                               stream='+', op_extend='1x')

        assert s == -20
        assert e == 140

    def test8(self):
        start, end = 20, 80
        s, e = extend_interval(start, end, '2x', max_extend=40,
                               stream='+', op_extend=300, max_op_extend='2x')

        assert s == -20
        assert e == 200

    def test9(self):
        start, end = 20, 80
        s, e = extend_interval(start, end, '3x',
                               stream='+', op_extend='1x', max_op_extend=100)

        assert s == -160
        assert e == 140

    def test10(self):
        start, end = 20, 80
        s, e = extend_interval(start, end, '3x',
                               stream='-', op_extend='1x', max_op_extend=100)

        assert s == -40
        assert e == 260

if __name__ == '__main__':
    unittest.main()

