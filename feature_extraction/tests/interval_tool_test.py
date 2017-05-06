"""
Created on Aug 1, 2016

@author: ramseylab
"""
import unittest
from io import StringIO
import pandas
from interval_tool import match, group_match


class IntervalToolTest(unittest.TestCase):
    
    def testMatchWhenYEmpty(self):
        x_str = """name    chromStart
            rs0    50
            rs1    150
            rs2    250
            rs3    350
            rs4    450
            rs5    550
            rs6    650     
            rs7    750"""
            
        x = pandas.read_csv(StringIO(x_str), encoding='utf8', header=0,
                            delim_whitespace=True, skipinitialspace=True)
        
        y = pandas.DataFrame()
        
        result = match(x, "name", "chromStart", y, "chromStart", "chromEnd", "score", True)
        
        self.assertEqual(result.loc['rs0', 'score'], [])
        self.assertEqual(result.loc['rs1', 'score'], [])
        self.assertEqual(result.loc['rs2', 'score'], [])
        self.assertEqual(result.loc['rs3', 'score'], [])
        self.assertEqual(result.loc['rs4', 'score'], [])
        self.assertEqual(result.loc['rs5', 'score'], [])
        self.assertEqual(result.loc['rs6', 'score'], [])
        self.assertEqual(result.loc['rs7', 'score'], [])

    def testMatchSorted(self):
        x_str = """name    chromStart
            rs0    50
            rs1    150
            rs2    250
            rs3    350
            rs4    450
            rs5    550
            rs6    650     
            rs7    750"""
            
        x = pandas.read_csv(StringIO(x_str), encoding='utf8', header=0,
                            delim_whitespace=True, skipinitialspace=True)
        
        y_str = """chromStart    chromEnd    score
            20    100    0
            40    60    1
            51    99    2
            101    151    3
            131    181    4
            200    260    5
            300    360    6"""
            
        y = pandas.read_csv(StringIO(y_str), encoding='utf8', header=0,
                            delim_whitespace=True, skipinitialspace=True)
        
        result = match(x, "name", "chromStart", y, "chromStart", "chromEnd", "score", True)
        
        self.assertEqual(result.loc['rs0', 'score'], [0, 1])
        self.assertEqual(result.loc['rs1', 'score'], [3, 4])
        self.assertEqual(result.loc['rs2', 'score'], [5])
        self.assertEqual(result.loc['rs3', 'score'], [6])
        self.assertEqual(result.loc['rs4', 'score'], [])
        self.assertEqual(result.loc['rs5', 'score'], [])
        self.assertEqual(result.loc['rs6', 'score'], [])
        self.assertEqual(result.loc['rs7', 'score'], [])
        
    def testMatchUnsorted(self):
        x_str = """name    chromStart
            rs2    250
            rs5    550
            rs4    450
            rs1    150
            rs3    350
            rs0    50
            rs7    750
            rs6    650     
            """
            
        x = pandas.read_csv(StringIO(x_str), encoding='utf8', header=0,
                            delim_whitespace=True, skipinitialspace=True)
        
        y_str = """chromStart    chromEnd    score
            51    99    2
            20    100    0
            131    181    4
            40    60    1
            300    360    6
            101    151    3
            200    260    5
            """
            
        y = pandas.read_csv(StringIO(y_str), encoding='utf8', header=0,
                            delim_whitespace=True, skipinitialspace=True)
        
        result = match(x, "name", "chromStart", y, "chromStart", "chromEnd", "score", False)
        
        self.assertEqual(result.loc['rs0', 'score'], [0, 1])
        self.assertEqual(result.loc['rs1', 'score'], [3, 4])
        self.assertEqual(result.loc['rs2', 'score'], [5])
        self.assertEqual(result.loc['rs3', 'score'], [6])
        self.assertEqual(result.loc['rs4', 'score'], [])
        self.assertEqual(result.loc['rs5', 'score'], [])
        self.assertEqual(result.loc['rs6', 'score'], [])
        self.assertEqual(result.loc['rs7', 'score'], [])
        
    def testGroupMatch(self):
        x_str = """name    chromStart    chrom
            rs2    250    chr1
            rs5    550    chr2
            rs4    450    chr2
            rs1    150    chr1
            rs3    350    chr3
            rs0    50    chr3
            rs7    750    chr4
            rs6    650     chr4
            """
            
        x = pandas.read_csv(StringIO(x_str), encoding='utf8', header=0,
                            delim_whitespace=True, skipinitialspace=True)
        
        y_str = """chromStart    chromEnd    score    chrom
            51    99    2    chr1
            20    100    0    chr3
            131    181    4    chr1
            40    60    1    chr3
            300    360    6    chr3
            101    151    3    chr2
            200    260    5    chr2
            """
            
        y = pandas.read_csv(StringIO(y_str), encoding='utf8', header=0,
                            delim_whitespace=True, skipinitialspace=True)
        
        result = group_match(x, "chrom", "name", "chromStart", y, "chrom", "chromStart", "chromEnd", "score")
        
        self.assertEqual(result.loc['rs0', 'score'], [0, 1])
        self.assertEqual(result.loc['rs1', 'score'], [4])
        self.assertEqual(result.loc['rs2', 'score'], [])
        self.assertEqual(result.loc['rs3', 'score'], [6])
        self.assertEqual(result.loc['rs4', 'score'], [])
        self.assertEqual(result.loc['rs5', 'score'], [])
        self.assertEqual(result.loc['rs6', 'score'], [])
        self.assertEqual(result.loc['rs7', 'score'], [])

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
