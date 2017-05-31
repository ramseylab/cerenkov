import unittest
import pandas as pd
import numpy as np
from allele_util import AlleleUtil
from allele_tool import parse_gb_allele_freqs_col, parse_gb_alleles_col


class AlleleUtilCase(unittest.TestCase):
    def test_transform_cols(self):
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

        df = pd.DataFrame(data=data, columns=['name', 'alleles', 'alleleFreqs'])

        # for i in range(0, len(data)):
        #     df.loc[i] = data[i]

        df.loc[:, "alleles"] = df["alleles"].apply(parse_gb_alleles_col)
        df.loc[:, "alleleFreqs"] = df["alleleFreqs"].apply(parse_gb_allele_freqs_col)

        result = AlleleUtil.transform_cols(df)

        self.assertEqual(result.loc[0, 'A'], 0)
        self.assertEqual(result.loc[0, 'T'], 0)
        self.assertEqual(result.loc[0, 'C'], 0)
        self.assertEqual(result.loc[0, 'G'], 0)

        self.assertEqual(result.loc[1, 'A'], 0.2)
        self.assertEqual(result.loc[1, 'T'], 0)
        self.assertEqual(result.loc[1, 'C'], 0.7)
        self.assertEqual(result.loc[1, 'G'], 0)

        self.assertEqual(result.loc[2, 'A'], 0)
        self.assertEqual(result.loc[2, 'T'], 0)
        self.assertEqual(result.loc[2, 'C'], 0)
        self.assertEqual(result.loc[2, 'G'], 0)

        self.assertEqual(result.loc[3, 'A'], 0)
        self.assertEqual(result.loc[3, 'T'], 0)
        self.assertEqual(result.loc[3, 'C'], 1)
        self.assertEqual(result.loc[3, 'G'], 0)

        self.assertEqual(result.loc[4, 'A'], 0.2)
        self.assertEqual(result.loc[4, 'T'], 0)
        self.assertEqual(result.loc[4, 'C'], 0.8)
        self.assertEqual(result.loc[4, 'G'], 0)

        self.assertEqual(result.loc[5, 'A'], 0.8)
        self.assertEqual(result.loc[5, 'T'], 0)
        self.assertEqual(result.loc[5, 'C'], 0.2)
        self.assertEqual(result.loc[5, 'G'], 0)

        self.assertEqual(result.loc[6, 'A'], 0)
        self.assertEqual(result.loc[6, 'T'], 0)
        self.assertEqual(result.loc[6, 'C'], 0)
        self.assertEqual(result.loc[6, 'G'], 0)

        self.assertEqual(result.loc[7, 'A'], 0)
        self.assertEqual(result.loc[7, 'T'], 0)
        self.assertEqual(result.loc[7, 'C'], 0)
        self.assertEqual(result.loc[7, 'G'], 0)

    def test_maf_filter(self):
        data = [
            # ['name',   'A',  'T',  'C',  'G']
            ['rs1002076', 0,    0,    0,    0],     # biomart; keep
            ['rs1002076', 1,    0,    0,    0],     # biomart; keep
            ['rs3',       0.2,  0.5,  0.3,  0],     # biomart; discard
            ['rs4',       0.02, 0.68, 0.3,  0],     # keep
            ['rs5',       0.2,  0.2,  0.3,  0.3],   # biomart; discard
            ['rs6',       0.02, 0.3,  0.2,  0.48],  # biomart; discard
            ['rs7',       0.02, 0.38, 0.02, 0.58],  # keep
            ['rs1002076', 0,    0,    0.99, 0],     # biomart; keep
            ['rs1002076', 0.01, 0,    1,    0],     # biomart; keep
            ['rs10',      0.3,  0,    0.7,  0],     # keep
            ['rs11',      0.03, 0,    0.97, 0],     # biomart; discard
        ]

        df = pd.DataFrame(data=data, columns=['name', 'A', 'T', 'C', 'G'])

        # is_eligible = df.loc[:, list('ATCG')].apply(lambda x: sum(x >= 0.05) == 2, axis=1)
        # self.assertEqual(len(is_eligible), 11)
        # self.assertEqual(sum(is_eligible), 3)

        fltr_df = AlleleUtil.maf_filter(df, maf_threshold=0.05, use_biomart=False)
        fltr_df = fltr_df.set_index('name')
        self.assertEqual(fltr_df.shape[0], 3)
        self.assertTrue('rs4' in fltr_df.index)
        self.assertTrue('rs7' in fltr_df.index)
        self.assertTrue('rs10' in fltr_df.index)


if __name__ == '__main__':
    unittest.main()
