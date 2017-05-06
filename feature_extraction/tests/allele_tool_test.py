"""
Created on Jul 11, 2016

@author: Yao
"""

import unittest
from allele_tool import parse_bm_alleles_col, parse_gb_alleles_col, parse_gb_allele_freqs_col, flip_allele, build_allele_freq_map


class AlleleToolTest(unittest.TestCase):
    def test_parse_bm_alleles_col(self):
        _str_1 = 'A/G'
        _str_2 = ''
        self.assertEqual(parse_bm_alleles_col(_str_1), ['A', 'G'])
        self.assertEqual(parse_bm_alleles_col(_str_2), [])

    def test_parse_gb_alleles_col(self):
        _bytes_1 = b'A,'
        _bytes_2 = b'A,C,'
        _bytes_3 = b''

        self.assertEqual(parse_gb_alleles_col(_bytes_1), ['A'])
        self.assertEqual(parse_gb_alleles_col(_bytes_2), ['A', 'C'])
        self.assertEqual(parse_gb_alleles_col(_bytes_3), [])

    def test_parse_gb_allele_freqs_col(self):
        _bytes_1 = b'0.2,'
        _bytes_2 = b'0.2,0.8,'
        _bytes_3 = b''

        self.assertEqual(parse_gb_allele_freqs_col(_bytes_1), [0.2])
        self.assertEqual(parse_gb_allele_freqs_col(_bytes_2), [0.2, 0.8])
        self.assertEqual(parse_gb_allele_freqs_col(_bytes_3), [])

    def test_flip_allele(self):
        self.assertEqual(flip_allele(['T']), ['A'], 'T => A')
        self.assertEqual(flip_allele(['A']), ['T'], 'A => T')
        self.assertEqual(flip_allele(['G']), ['C'], 'G => C')
        self.assertEqual(flip_allele(['C']), ['G'], 'C => G')
        
        self.assertEqual(flip_allele(['A', 'C']), ['T', 'G'], 'A/C => T/G')
        self.assertEqual(flip_allele(['A', 'C', 'T']), ['T', 'G', 'A'], 'A/C/T => T/G/A')

        self.assertEqual(flip_allele('T'), 'A', 'T => A')
        self.assertEqual(flip_allele('A'), 'T', 'A => T')
        self.assertEqual(flip_allele('G'), 'C', 'G => C')
        self.assertEqual(flip_allele('C'), 'G', 'C => G')

        self.assertEqual(flip_allele('AC'), 'TG', 'AC => TG')
        self.assertEqual(flip_allele('ACT'), 'TGA', 'ACT => TGA')
        
    def test_build_allele_freq_map(self):
        m_1 = build_allele_freq_map(['G', 'A', 'T'], [0.3, 0.5, 0.2])
        m_2 = build_allele_freq_map(['C', 'T', 'A'], [0.2, 0.5, 0.3])
        m_3 = build_allele_freq_map([], [])
        m_4 = build_allele_freq_map(['G', 'A', '0'], [0.3, 0.7, 0])
        m_5 = build_allele_freq_map(['0', 'T', 'A'], [0, 0.7, 0.3])
        
        self.assertEqual(m_1, [('T', 0.2), ('G', 0.3), ('A', 0.5)])
        self.assertEqual(m_2, [('C', 0.2), ('A', 0.3), ('T', 0.5)])
        self.assertEqual(m_3, [])
        self.assertEqual(m_4, [('G', 0.3), ('A', 0.7)])
        self.assertEqual(m_5, [('A', 0.3), ('T', 0.7)])


if __name__ == '__main__':
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
