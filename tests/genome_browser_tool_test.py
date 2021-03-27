import unittest
from genome_browser_tool import bin_from_range


class GenomeBrowserToolTestCase(unittest.TestCase):
    def test_bin_from_range(self):
        _bin = bin_from_range(10019, 10020)
        self.assertEqual(_bin, 585)


if __name__ == '__main__':
    unittest.main()
