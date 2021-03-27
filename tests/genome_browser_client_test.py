import unittest
from genome_browser_client import GenomeBrowserClient


class MySqlClientTestCase(unittest.TestCase):
    def test_in_protein_coding_exon(self):
        with GenomeBrowserClient('local_hg19') as gb_client:
            is_in = gb_client.in_protein_coding_exon(chrom='chr7', chrom_start=117306990)
            self.assertTrue(is_in)


if __name__ == '__main__':
    unittest.main()
