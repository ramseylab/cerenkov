"""
Created on Jul 29, 2016

@author: ramseylab
"""
import unittest
from biomart_client import BiomartClient


class BiomartClientTest(unittest.TestCase):

    # def test_filter(self):
    #     df_str = """name,chrom,strand,alleles,minorAllele,minorAlleleFreq
    #         rs1,1,1,G/A,A,0.333666
    #         rs2,1,1,C/G/T,G,0.019169
    #         rs3,1,1,C/T,C,NaN
    #         rs4,1,1,C/T,G,0
    #         rs5,HSCHR6_MHC_APD,1,T/C,T,0.444"""
    #
    #     df = pandas.read_csv(StringIO(df_str), encoding='utf8', header=0, sep=",", skipinitialspace=True)
    #     df.loc[:, "alleles"] = df["alleles"].apply(parse_bm_alleles_col)
    #
    #     result = BiomartClient.filter(df)
    #
    #     self.assertEqual(result.shape, (1, 6))
    #     self.assertEqual(result.loc[0, "name"], "rs1")
    #
    # def test_transform(self):
    #     df_str = """name,chrom,strand,alleles,minorAllele,minorAlleleFreq
    #         rs1,1,1,G/A,A,0.333666
    #         rs2,1,-1,C/G,G,0.019169"""
    #
    #     df = pandas.read_csv(StringIO(df_str), encoding='utf8', header=0, sep=",", skipinitialspace=True)
    #     df.loc[:, "alleles"] = df["alleles"].apply(parse_bm_alleles_col)
    #
    #     result = BiomartClient.transform(df)
    #
    #     self.assertEqual(result.shape, (2, 8))
    #
    #     self.assertEqual(result.loc[0, "name"], "rs1")
    #     self.assertEqual(result.loc[0, "majorAllele"], "G")
    #     self.assertEqual(result.loc[0, "majorAlleleFreq"], 1 - 0.333666)
    #
    #     self.assertEqual(result.loc[1, "name"], "rs2")
    #     self.assertEqual(result.loc[1, "alleles"], ["G", "C"])
    #     self.assertEqual(result.loc[1, "minorAllele"], "C")
    #     self.assertEqual(result.loc[1, "majorAllele"], "G")
    #     self.assertEqual(result.loc[1, "majorAlleleFreq"], 1 - 0.019169)
        
    def test_query_snp(self):
        rsid_list = ["rs1048238", "rs1002076", "rs2232945"]

        with BiomartClient() as bm_client:
            result = bm_client.query_snp(rsid_list)

            self.assertEqual(result.shape, (1, 6))

            snp = result.loc[0, :]
            self.assertEqual(snp["name"], "rs1002076")
            self.assertEqual(snp["alleles"], ['G', 'A'])
            self.assertEqual(snp["minorAllele"], "A")
            self.assertAlmostEqual(snp["minorAlleleFreq"], 0.333666)
            self.assertEqual(snp["majorAllele"], "G")
            self.assertAlmostEqual(snp["majorAlleleFreq"], 0.666334)

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
