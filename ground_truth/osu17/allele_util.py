import pandas as pd
import numpy as np
from abstract_feature_util import AbstractFeatureUtil
from genome_browser_client import GenomeBrowserClient
from biomart_client import BiomartClient2
from chrom_tool import remove_dup_on_chrY
from allele_tool import flip_allele


class AlleleUtil(AbstractFeatureUtil):
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

    @staticmethod
    def maf_filter(snp_dfm, maf_threshold=0.05, use_biomart=False, verbose=False):
        def __is_freq_valid(freq):
            return (freq >= maf_threshold) & (freq > 0)

        def __is_allele(freq):
            return freq > 0

        # We call a SNP not valid if:
        #   CASE 1: it has only one allele;
        #   CASE 2: it has two allele but maf < 0.05 (`self.maf_thld`);
        #   CASE 3: it has three or four allele but you cannot tell which is minor because there are at least 3
        #       alleles with freq >= 0.05. (If there is one allele with freq < 0.05, we simply discard it)
        # SNPs of CASE 1 and CASE 3 will be queried in Biomart for a second chance if you set `self.use_biomart`
        snp_valid = snp_dfm.loc[:, list('ATCG')].apply(lambda x: sum(__is_freq_valid(x)) == 2, axis=1, reduce=True)

        if use_biomart:
            snp_biallelic = snp_dfm.loc[:, list('ATCG')].apply(lambda x: sum(__is_allele(x)) == 2, axis=1, reduce=True)

            with BiomartClient2() as bm_client:
                bm_df = bm_client.query_snp(rsid_list=snp_dfm.loc[~snp_valid & ~snp_biallelic, 'name'].tolist(),
                                            verbose=verbose)
                bm_df = remove_dup_on_chrY(bm_df)

                bm_valid = bm_df.loc[:, list('ATCG')].apply(lambda x: sum(__is_freq_valid(x)) == 2, axis=1)

                df = pd.concat([snp_dfm.loc[snp_valid, :], bm_df.loc[bm_valid, :]], axis=0)
        else:
            df = snp_dfm.loc[snp_valid, :]

        return df

    @staticmethod
    def identify_major_vs_minor(snp_dfm, maf_threshold=0.05, remove_acgt_cols=False, revise_when_freqs_equal=False, verbose=False):
        def identify(freq_series):
            valid_freqs = freq_series[(freq_series >= maf_threshold) & (freq_series > 0)]

            minor_freq = min(valid_freqs)
            major_freq = max(valid_freqs)

            minor_allele = valid_freqs.idxmin()
            major_allele = valid_freqs.idxmax()

            return pd.Series({"minorAlleleFreq": minor_freq,
                              "majorAlleleFreq": major_freq,
                              "minorAllele": minor_allele,
                              "majorAllele": major_allele
                              })

        mm_dfm = snp_dfm.loc[:, list('ATCG')].apply(lambda x: identify(x), axis=1, reduce=False)
        snp_dfm = pd.concat([snp_dfm, mm_dfm], axis=0)

        if remove_acgt_cols:
            if verbose:
                print("[identify_major_vs_minor] going to drop columns 'A', 'C', 'G' and 'T'")
            snp_dfm = snp_dfm.drop(list('ATCG'), axis=1)

        if revise_when_freqs_equal:
            """
            There exist SNPs with majorAlleleFreq == minorAlleleFreq.

            Steve:

            > For those 57 SNPs, I would maybe use the following approach:
            > if the last (least significant) digit of the chromosomal coordinate is even, select the first (out of two)
                alleles listed as the minor allele.
            > If the last (least significant) digit of the chromosomal coordinate is odd, select the second (out of two)
                alleles listed as the minor allele.
            > This approach is deterministic but should more or less "randomly" and with equal probability assign the
                first or second allele to be the 'minor allele'.
            """
            if 'chromStart' not in snp_dfm.columns:
                raise KeyError("'chromStart' not in the data frame's columns")

            has_equal_freq = (snp_dfm.loc[:, "minorAlleleFreq"] == snp_dfm.loc[:, "majorAlleleFreq"])

            has_odd_chrom_start = (snp_dfm.loc[:, "chromStart"] % 2 == 1)

            idx = has_equal_freq & has_odd_chrom_start

            if verbose:
                print("[identify_major_vs_minor] alleles with equal freqs detected: \r\n{}".format(
                      snp_dfm.loc[has_equal_freq, ["name", "chrom", "chromStart", "chromEnd",
                                                   "majorAllele", "minorAllele",
                                                   "majorAlleleFreq", "minorAlleleFreq"]]))

                print("[identify_major_vs_minor] entries to swap alleles: \r\n{}".format(
                      snp_dfm.loc[idx, ["name", "chrom", "chromStart", "chromEnd",
                                        "majorAllele", "minorAllele",
                                        "majorAlleleFreq", "minorAlleleFreq"]]))

            # Swap
            snp_dfm.loc[idx, ["majorAllele", "minorAllele"]] = snp_dfm.loc[idx, ["minorAllele", "majorAllele"]].values

        return snp_dfm

    @staticmethod
    def flip_alleles_on_rev_strand(snp_dfm, remove_strand_col=False, verbose=False):
        on_rev = (snp_dfm.loc[:, "strand"] == "-")

        if not snp_dfm.loc[on_rev, :].empty:
            snp_dfm.loc[on_rev, "majorAllele"] = snp_dfm.loc[on_rev, "majorAllele"].apply(flip_allele)
            snp_dfm.loc[on_rev, "minorAllele"] = snp_dfm.loc[on_rev, "minorAllele"].apply(flip_allele)

        if remove_strand_col:
            if verbose:
                print("[flip_alleles_on_rev_strand] going to drop columns 'strand'")

            snp_dfm = snp_dfm.drop('strand', axis=1)

        return snp_dfm

    def get_feat(self, _input):
        rsid = _input

        with GenomeBrowserClient(self.db_config_key) as gb_client:
            gb_df = gb_client.fetch_alleles(rsid)
            gb_df = remove_dup_on_chrY(gb_df)

            gb_df = AlleleUtil.transform_cols(gb_df)

            return gb_df

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', index=False, header=True)


if __name__ == '__main__':
    allele_util = AlleleUtil()
    allele_util.db_config_key = 'local_hg19'

    allele_dfm = allele_util.extract(_input=['rs115493313', 'rs2297233'])
    print(allele_dfm)