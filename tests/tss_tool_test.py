import unittest

from tss_tool import min_tss_dist


class TssToolTestCase(unittest.TestCase):
    def test_min_tss_dist(self):
        self.assertEqual(min_tss_dist([1, 2, 3]), 1)
        self.assertEqual(min_tss_dist([-1, -2, -3]), -1)
        self.assertEqual(min_tss_dist([1, -2, 3]), 1)
        self.assertEqual(min_tss_dist([1, 1, 2]), 1)
        self.assertEqual(min_tss_dist([-1, 1, 1]), -1)


if __name__ == '__main__':
    unittest.main()
