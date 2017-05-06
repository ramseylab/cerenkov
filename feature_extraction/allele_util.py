import pandas as pd
import numpy as np
from abstract_feature_util import AbstractFeatureUtil
from genome_browser_client import GenomeBrowserClient
from biomart_client import BiomartClient2
from chrom_tool import remove_dup_on_chrY


class AlleleUtil(AbstractFeatureUtil):
    def __init__(self, use_biomart=True, maf_threshold=0.05):
        super(self.__class__, self).__init__()
        self.use_biomart = use_biomart
        self.maf_thld = maf_threshold

    @staticmethod
    def transform_cols(df):
        df = df.assign(A=np.nan, T=np.nan, C=np.nan, G=np.nan)

        # N.B. it is possible that `alleles` column contains a '0' for uncertain nucleobase.
        # It will be omitted in this step
        df = df.apply(lambda x: x.fillna(dict(zip(x['alleles'], x['alleleFreqs']))), axis=1)
        df = df.fillna(0)
        df = df.drop(['alleles', "alleleFreqs"], axis=1)

        return df

    # @staticmethod
    # def __print_allele_summary(df):
    #     # Identify mono-, bi-, tri-, quad-allelic
    #     n_allele = df.loc[:, list("ATCG")].apply(lambda x: sum(x > 0), axis=1)
    #
    #     allele_types = ["Mono-allelic", "Bi-allelic", "Tri-allelic", "Quad-allelic", 'Total']
    #
    #     cnt_mono = df.loc[n_allele == 1].shape[0]
    #     cnt_bi = df.loc[n_allele == 2].shape[0]
    #     cnt_tri = df.loc[n_allele == 3].shape[0]
    #     cnt_quad = df.loc[n_allele == 4].shape[0]
    #     cnt_total = cnt_mono + cnt_bi + cnt_tri + cnt_quad
    #     counts = [cnt_mono, cnt_bi, cnt_tri, cnt_quad, cnt_total]
    #
    #     summary = pd.DataFrame(counts, allele_types, ["count"])
    #     print(summary)

    def _freq_filter(self, freq):
        if self.maf_thld > 0:
            return freq >= self.maf_thld
        else:
            return freq > 0  # filter `freq >= 0` is meaningless

    def get_feat(self, _input):
        rsid = _input

        with GenomeBrowserClient(self.db_config_key) as gb_client:
            gb_df = gb_client.fetch_alleles(rsid)
            gb_df = remove_dup_on_chrY(gb_df)

            gb_df = AlleleUtil.transform_cols(gb_df)

            # self.__print_allele_summary(gb_df)

            # We call a SNP not valid if:
            #   CASE 1: it has only one allele;
            #   CASE 2: it has two allele but maf < 0.05 (`self.maf_thld`);
            #   CASE 3: it has three or four allele but you cannot tell which is minor because there are at least 3
            #       alleles with freq >= 0.05. (If there is one allele with freq < 0.05, we simply discard it)
            # SNPs of CASE 1 and CASE 3 will be queried in Biomart for a second chance if you set `self.use_biomart`
            gb_valid = gb_df.loc[:, list('ATCG')].apply(lambda x: sum(self._freq_filter(x)) == 2, axis=1)
            gb_biallelic = gb_df.loc[:, list('ATCG')].apply(lambda x: sum(x > 0) == 2, axis=1)

            if self.use_biomart:
                with BiomartClient2() as bm_client:
                    bm_df = bm_client.query_snp(rsid_list=gb_df.loc[~gb_valid & ~gb_biallelic, 'name'].tolist(),
                                                verbose=True)
                    bm_df = remove_dup_on_chrY(bm_df)

                    bm_valid = bm_df.loc[:, list('ATCG')].apply(lambda x: sum(self._freq_filter(x)) == 2, axis=1)

                    df = pd.concat([gb_df.loc[gb_valid, :], bm_df.loc[bm_valid, :]], axis=0)
            else:
                df = gb_df.loc[gb_valid, :]

            return df

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', index=False, header=True)


if __name__ == '__main__':
    allele_util = AlleleUtil()
    allele_util.db_config_key = 'local_hg19'

    allele_dfm = allele_util.extract(_input=['rs115493313', 'rs2297233'])
    print(allele_dfm)