"""
Created on Aug 9, 2016

@author: ramseylab
"""

import unittest
import pandas
import numpy
from tss_util import TssDistUtil


class TssUtil2Test(unittest.TestCase):

    def testMinTssDist(self):
        dfm = pandas.DataFrame({
            'name': ["rs1", "rs1", "rs2", "rs2", "rs3", "rs3", "rs3"], 
            'chrom': ["chr1", "chr1", "chr2", "chr2", "chr3", "chr3", "chrY"], 
            "tssDistance": [1, -1, 3, 4, 5, numpy.NaN, 6]
        })
        min_td = TssDistUtil._select_min(dfm)
        min_td = min_td.set_index(["name", "chrom"])

        self.assertEqual(min_td.shape[0], 4)
        self.assertEqual(min_td.loc[("rs1", "chr1"), "tssDistance"], -1.0)
        self.assertEqual(min_td.loc[("rs2", "chr2"), "tssDistance"], 3.0)
        self.assertEqual(min_td.loc[("rs3", "chr3"), "tssDistance"], 5.0)
        self.assertEqual(min_td.loc[("rs3", "chrY"), "tssDistance"], 6.0)

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
