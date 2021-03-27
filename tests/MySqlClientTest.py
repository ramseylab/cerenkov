import unittest
from MySqlClient import in_protein_coding_exon


class MySqlClientTestCase(unittest.TestCase):
    def test_in_protein_coding_exon(self):
        is_in = in_protein_coding_exon(chrom='chr7', chrom_start=117306990)
        self.assertTrue(is_in)


if __name__ == '__main__':
    unittest.main()
