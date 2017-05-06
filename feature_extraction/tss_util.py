import pandas
from genome_browser_client import GenomeBrowserClient
import chrom_tool as ct
from tss_tool import min_tss_dist
from abstract_feature_util import AbstractFeatureUtil


# def _get_candidate_dfm(rsid):
#     with GenomeBrowserClient('local_hg19') as gb_client:
#         first_run_result = gb_client.compute_tss_dist(rsid, adjacent_bins=1)
#
#         remainder = rsid[~rsid.isin(first_run_result.loc[:, "name"])]
#         if not remainder.empty:
#             print("[tss_util.__compute_tss_dist]: No distance found for {n} SNP(s) after 1st run: \r\n{snps}".
#                   format(n=remainder.shape[0], snps=remainder.tolist()))
#
#         second_run_result = gb_client.compute_tss_dist(remainder, adjacent_bins=-1)
#
#         remainder = remainder[~remainder.isin(second_run_result.loc[:, "name"])]
#
#         if not remainder.empty:
#             print("[tss_util.__compute_tss_dist]: No distance found for {n} SNP(s) after 2st run: \r\n{snps}".
#                   format(n=remainder.shape[0], snps=remainder.tolist()))
#
#         # IMPORTANT: do apply `reset_index()` or set `ignore_index = True` here
#         # if not, multiple entries would share a same index
#         #   which would be a disaster for `idxmin()` operation in `__minTssDist`
#         # return pandas.concat([firstRunResult, secondRunResult]).reset_index()
#         return pandas.concat([first_run_result, second_run_result], ignore_index=True)
#
#
# def _select_min(candidate_dfm):
#     """
#     Select the closest distances from TSS candidates.
#
#     N.B. Negative distance to TSS (which means the SNP is located to the upstream of the TSS)
#         is given higher priority if there existed 2 TSS's for a single SNP and the absolute values
#         of these 2 distances are equal.
#
#     E.g. for SNP rs123, if the distance to TSS-1 is -5 and 5 to TSS-2, we pick -5 as the closest one.
#
#     We prefer negative TSS distances because we prefer upstream SNPs.
#
#     The preference is performed by adding 0.1 to the distances and then rank their absolute value ascendingly.
#
#     E.g. |-5 + 0.1| < |5 + 0.1|, so -5 is picked.
#     """
#
#     # add a new column named `amenedAbsTssDistance`
#     # add 0.1 here to provide negative distances with higher priorities
#     # candidate_dfm.loc[:, "amenedAbsTssDistance"] = candidate_dfm.loc[:, "tssDistance"].add(0.1).abs()
#     #
#     # min_index = candidate_dfm.groupby(["name", "chrom"]).apply(
#     #     lambda group: group["amenedAbsTssDistance"].idxmin(skipna=True))
#     #
#     # return candidate_dfm.loc[min_index, :]
#
#     # We apply to a single column of each group. Therefore `apply` would return a MultiIndex Series
#     # `reset_index()` to a Series would return a DataFrame
#     # `reset_index(name='xxx')` specifies the name of the column corresponding to the Series values
#     #   in the output DataFrame
#     min_tss_dist_dfm = candidate_dfm.groupby(["name", "chrom"]).\
#         apply(lambda x: min_tss_dist(x['tssDistance'].tolist())).reset_index(name='tssDistance')
#
#     return min_tss_dist_dfm


# def extract_min_tss_dist(src, dest):
#     """
#     Exact TSS Distances (actually the distances to TSS, Transcription Start Site) in 3 steps:
#
#         (1) Run the sql to extract TSS candidates by querying UCSC database.
#         (2) Select the closest distances from TSS candidates.
#         (3) Remove duplication on the chrY entries
#     """
#
#     rsid = pandas.read_csv(src, sep='\t').loc[:, "name"]
#
#     candidate_dfm = _get_candidate_dfm(rsid)
#
#     min_tss_dist_dfm = _select_min(candidate_dfm)
#     min_tss_dist_dfm = ct.remove_dup_on_chrY(min_tss_dist_dfm)
#
#     min_tss_dist_dfm.loc[:, ['name', 'tssDistance']].to_csv(dest, header=True, index=False, sep='\t')


class TssDistUtil(AbstractFeatureUtil):
    def __get_candidate_dfm(self, rsid):
        with GenomeBrowserClient(self.db_config_key) as gb_client:
            first_run_result = gb_client.compute_tss_dist(rsid, adjacent_bins=1)

            remainder = rsid[~rsid.isin(first_run_result.loc[:, "name"])]
            if not remainder.empty:
                print("[tss_util.__compute_tss_dist]: No distance found for {n} SNP(s) after 1st run: \r\n{snps}".
                      format(n=remainder.shape[0], snps=remainder.tolist()))

            second_run_result = gb_client.compute_tss_dist(remainder, adjacent_bins=-1)

            remainder = remainder[~remainder.isin(second_run_result.loc[:, "name"])]

            if not remainder.empty:
                print("[tss_util.__compute_tss_dist]: No distance found for {n} SNP(s) after 2st run: \r\n{snps}".
                      format(n=remainder.shape[0], snps=remainder.tolist()))

            # IMPORTANT: do apply `reset_index()` or set `ignore_index = True` here
            # if not, multiple entries would share a same index
            #   which would be a disaster for `idxmin()` operation in `__minTssDist`
            # return pandas.concat([firstRunResult, secondRunResult]).reset_index()
            return pandas.concat([first_run_result, second_run_result], ignore_index=True)

    @staticmethod
    def _select_min(candidate_dfm):
        """
        Select the closest distances from TSS candidates.

        N.B. Negative distance to TSS (which means the SNP is located to the upstream of the TSS)
            is given higher priority if there existed 2 TSS's for a single SNP and the absolute values
            of these 2 distances are equal.

        E.g. for SNP rs123, if the distance to TSS-1 is -5 and 5 to TSS-2, we pick -5 as the closest one.

        We prefer negative TSS distances because we prefer upstream SNPs.

        The preference is performed by adding 0.1 to the distances and then rank their absolute value ascendingly.

        E.g. |-5 + 0.1| < |5 + 0.1|, so -5 is picked.
        """

        # We apply to a single column of each group. Therefore `apply` would return a MultiIndex Series
        # `reset_index()` to a Series would return a DataFrame
        # `reset_index(name='xxx')` specifies the name of the column corresponding to the Series values
        #   in the output DataFrame
        min_tss_dist_dfm = candidate_dfm.groupby(["name", "chrom"]). \
            apply(lambda x: min_tss_dist(x['tssDistance'].tolist())).reset_index(name='tssDistance')

        return min_tss_dist_dfm

    def get_feat(self, _input):
        """
        Exact TSS Distances (actually the distances to TSS, Transcription Start Site) in 3 steps:

        (1) Run the sql to extract TSS candidates by querying UCSC database.
        (2) Select the closest distances from TSS candidates.
        (3) Remove duplication on the chrY entries

        :param _input: a sequence of RSID
        :return:
        """

        # rsid = pandas.read_csv(src, sep='\t').loc[:, "name"]
        rsid = _input

        candidate_dfm = self.__get_candidate_dfm(rsid)

        min_tss_dist_dfm = self._select_min(candidate_dfm)
        min_tss_dist_dfm = ct.remove_dup_on_chrY(min_tss_dist_dfm)

        return min_tss_dist_dfm.loc[:, ['name', 'tssDistance']]

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, header=True, index=False, sep='\t')
