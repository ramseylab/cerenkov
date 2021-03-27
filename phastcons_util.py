from abstract_feature_util import AbstractFeatureUtil
from genome_browser_client import GenomeBrowserClient


class PhastconsUtil(AbstractFeatureUtil):
    def get_feat(self, _input):
        """
        :param _input: the SNP data frame
        :return:
        """
        snp_dfm = _input

        with GenomeBrowserClient(self.db_config_key) as gb_client:
            phastcons_df = gb_client.fetch_phastcons(snp_dfm.loc[:, 'name'].tolist())

        snp_dfm = snp_dfm.merge(phastcons_df, how='left', on=['name', 'chrom'], copy=True)

        return snp_dfm.loc[:, ["name", "phastCons"]].fillna(0)

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', header=True, index=False,
                       columns=["name", "phastCons"])
