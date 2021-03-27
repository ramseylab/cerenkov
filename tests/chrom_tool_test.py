"""
Created on Aug 9, 2016

@author: ramseylab
"""
import unittest
import pandas
from chrom_tool import remove_dup_on_chrY


class ChromToolTest(unittest.TestCase):

    def testRemoveDupOnChrY(self):
        dfm = pandas.DataFrame({
            'name': ["rs1", "rs2", "rs3", "rs3"], 
            'chrom': ["chr1", "chr2", "chr3", "chrY"], 
            "tssDistance": [1, 2, 3, 4]
        })
        
        dedup = remove_dup_on_chrY(dfm)
        
        self.assertEqual(dedup.shape[0], 3)
        self.assertFalse((dedup["chrom"] == "chrY").any())
        self.assertEqual((dedup["name"] == "rs3").sum(), 1)
        
    def testRemoveDupOnChrYError(self):
        dfm = pandas.DataFrame({
            'name': ["rs1", "rs2", "rs2", "rs3", "rs3"], 
            'chrom': ["chr1", "chr2", "chr3", "chr3", "chrY"], 
            "tssDistance": [1, 2, 3, 4, 5]
        })
        
        self.assertRaises(AssertionError, remove_dup_on_chrY, dfm)

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
