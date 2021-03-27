"""
Created on Jul 29, 2016

@author: ramseylab
"""

import unittest
import pandas.util.testing as pdt
from io import StringIO
import pandas
import ast
from allele_tool import parse_gb_alleles_col, parse_gb_allele_freqs_col
from metadata_util import __remove_non_regular_chrom as remove_non_regular_chrom
from metadata_util import __remove_non_single_class as remove_non_single_class
from metadata_util import __normalize_allele_strand as normalize_allele_strand
from metadata_util import __build_allele_freq_map as build_allele_freq_map
from metadata_util import __identify_major_minor_alleles as identify_major_minor_alleles
from metadata_util import __revise_alleles_with_equal_freqs as revise_alleles_with_equal_freqs
from metadata_util import __normalize_chrom_coord as normalize_chrom_coord


class MetadataUtilTest(unittest.TestCase):

    def test_remove_non_regular_chrom(self):
        df = pandas.DataFrame(columns=['name', 'chrom'])
        data = [
            ['rs1', 'chr1'],
            ['rs2', 'chr22'],
            ['rs3', 'chr23'],
            ['rs4', 'chrX'],
            ['rs5', 'chrY'],
            ['rs6', 'chrM'],
            ['rs7', 'chr6_cox_hap2'],
        ]
        for i in range(0, len(data)):
            df.loc[i] = data[i]
        
        result = remove_non_regular_chrom(df)
        
        self.assertEqual(result.shape[0], 4)
        pdt.assert_series_equal(result.loc[:, "chrom"], 
                                pandas.Series(data=["chr1", "chr22", "chrX", "chrY"], index=[0, 1, 3, 4]),
                                check_names=False)

    def test_remove_non_single_class(self):
        df = pandas.DataFrame(columns=['name', 'class'])
        data = [
            ['rs1', 'single'],
            ['rs2', 'insertion'],
            ['rs3', 'deletion'],
            ['rs4', 'indel'],
        ]
        for i in range(0, len(data)):
            df.loc[i] = data[i]

        result = remove_non_single_class(df)
        
        self.assertEqual(result.shape[0], 1)
        self.assertEqual(result.loc[0, "name"], "rs1")
        
    def test_normalize_allele_strand(self):
        df = pandas.DataFrame(columns=['name', 'strand', 'alleles'])
        data = [
            ['rs1', '+', b'A,C,'],
            ['rs2', '-', b'A,C,'],
            ['rs3', '+', b''],
            ['rs4', '-', b''],
        ]
        for i in range(0, len(data)):
            df.loc[i] = data[i]

        df.loc[:, "alleles"] = df["alleles"].apply(parse_gb_alleles_col)
        
        result = normalize_allele_strand(df)
        
        self.assertEqual(result.loc[0, "alleles"], ['A', 'C'])
        self.assertEqual(result.loc[1, "alleles"], ['T', 'G'])
        self.assertEqual(len(result.loc[2, "alleles"]), 0)
        self.assertEqual(len(result.loc[3, "alleles"]), 0)
    
    def test_build_allele_freq_map(self):
        df = pandas.DataFrame(columns=['name', 'alleles', 'alleleFreqs'])
        data = [
            ['rs0', b'', b''],
            ['rs1', b'0,A,C,', b'0.1,0.2,0.7,'],
            ['rs2', b'A,', b'0.0,'],
            ['rs3', b'C,', b'1.0,'],
            ['rs4', b'A,C,', b'0.2,0.8,'],
            ['rs5', b'A,C,', b'0.8,0.2,'],
            ['rs6', b'T,G,', b''],
            ['rs7', b'', b'0.5,0.5,'],
        ]
        for i in range(0, len(data)):
            df.loc[i] = data[i]

        df.loc[:, "alleles"] = df["alleles"].apply(parse_gb_alleles_col)
        df.loc[:, "alleleFreqs"] = df["alleleFreqs"].apply(parse_gb_allele_freqs_col)

        result = build_allele_freq_map(df)
        
        self.assertEqual(result.loc[0, "afMap"], [])
        self.assertEqual(result.loc[1, "afMap"], [("A", 0.2), ("C", 0.7)])
        self.assertEqual(result.loc[2, "afMap"], [("A", 0.0)])
        self.assertEqual(result.loc[3, "afMap"], [("C", 1.0)])
        self.assertEqual(result.loc[4, "afMap"], [("A", 0.2), ("C", 0.8)])
        self.assertEqual(result.loc[5, "afMap"], [("C", 0.2), ("A", 0.8)])
        self.assertEqual(result.loc[6, "afMap"], [])
        self.assertEqual(result.loc[7, "afMap"], [])
        
    def test_identify_major_minor_alleles(self):
        df_str = """name    afMap
            rs1002076    []
            rs1002076    [('A',1.0)]
            rs3    [('A',0.2),('C',0.3),('T',0.5)]
            rs4    [('A',0.02),('C',0.3),('T',0.68)]
            rs5    [('A',0.2),('C',0.2),('T',0.3),('G',0.3)]
            rs6    [('A',0.02),('C',0.2),('T',0.3),('G',0.48)]
            rs7    [('A',0.02),('C',0.02),('T',0.38),('G',0.58)]
            rs1002076    [('A',0.0),('C',0.99)]
            rs1002076    [('A',0.01),('C',1.0)]
            rs10    [('A',0.3),('C',0.7)]
            rs11    [('A',0.03),('C',0.97)]"""
        
        # `ast.literal_eval` converts string representation of a list into a real list
        # However, `ast.literal_eval` raises "SyntaxError: unexpected EOF while parsing"
        # if you pass an empty string to it
        df = pandas.read_csv(StringIO(df_str), encoding='utf8', header=0,
                             delim_whitespace=True, skipinitialspace=True, converters={'afMap': ast.literal_eval})
        
        result = identify_major_minor_alleles(df, verbose=True)
        
        self.assertEqual(sum(result.loc[:, "name"] == "rs10"), 1)  # 1 count
        self.assertEqual(sum(result.loc[:, "name"] == "rs1002076"), 4)  # 4 counts
        
        rs4 = (result.loc[:, "name"] == "rs4")
        self.assertEqual(sum(rs4), 1)
        
        # Because `rs4` is a Series and `result.loc[rs4, majorAllele]` cannot be treated as a single string
        # but a Series of stings, so you cannot compare `result.loc[rs4, majorAllele]` and "T" directly
        self.assertEqual(result.loc[rs4, "majorAllele"].item(), "T")
        self.assertAlmostEqual(result.loc[rs4, "majorAlleleFreq"].item(), 0.68)
        self.assertEqual(result.loc[rs4, "minorAllele"].item(), "C")
        self.assertAlmostEqual(result.loc[rs4, "minorAlleleFreq"].item(), 0.3)
        
        rs7 = (result.loc[:, "name"] == "rs7")
        self.assertEqual(sum(rs7), 1)
        self.assertEqual(result.loc[rs7, "majorAllele"].item(), "G")
        self.assertAlmostEqual(result.loc[rs7, "majorAlleleFreq"].item(), 0.58)
        self.assertEqual(result.loc[rs7, "minorAllele"].item(), "T")
        self.assertAlmostEqual(result.loc[rs7, "minorAlleleFreq"].item(), 0.38)
        
    def test_normalize_chrom_coord(self):
        df_str = """name    chrom    chromStart
            rs0    chr1    24895642"""
            
        df = pandas.read_csv(StringIO(df_str), encoding='utf8', header=0,
                             delim_whitespace=True, skipinitialspace=True)
        
        result = normalize_chrom_coord(df)
        
        self.assertAlmostEqual(result.loc[0, "normChromCoord"], 0.1)
    
    def test_revise_alleles_with_equal_freqs(self):
        df_str = """name    chrom    chromStart    chromEnd    majorAllele    minorAllele    majorAlleleFreq    minorAlleleFreq
            rs2024758    chr1    16308411    16308412    T    A    0.5    0.5
            rs74136818    chr1    186571524    186571525    T    C    0.5    0.5
            rs1775145    chr1    205756142    205756143    C    A    0.5    0.5
            rs73262817    chr12    7011123    7011124    G    C    0.5    0.5"""
            
        df = pandas.read_csv(StringIO(df_str), encoding='utf8', header=0,
                             delim_whitespace=True, skipinitialspace=True)
        
        result = revise_alleles_with_equal_freqs(df).set_index("name")
        
        self.assertEqual(result.shape[0], 4)
        
        self.assertEqual(result.loc["rs2024758", "majorAllele"], "A")
        self.assertEqual(result.loc["rs2024758", "minorAllele"], "T")
        
        self.assertEqual(result.loc["rs74136818", "majorAllele"], "T")
        self.assertEqual(result.loc["rs74136818", "minorAllele"], "C")
        
        self.assertEqual(result.loc["rs1775145", "majorAllele"], "C")
        self.assertEqual(result.loc["rs1775145", "minorAllele"], "A")
        
        self.assertEqual(result.loc["rs73262817", "majorAllele"], "C")
        self.assertEqual(result.loc["rs73262817", "minorAllele"], "G")

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
