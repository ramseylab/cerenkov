from genome_browser_client import GenomeBrowserClient
from abstract_feature_util import AbstractFeatureUtil
import pandas as pd


class NkiLadUtil(AbstractFeatureUtil):
    def get_feat(self, _input):
        """

        :param _input: the SNP data frame in BED format
        :return:
        """
        snps = _input

        with GenomeBrowserClient(self.db_config_key) as gb_client:
            in_lad = snps.apply(lambda x: gb_client.in_nki_lad(x['chrom'], x['chromStart']), axis=1)

            overlap = pd.DataFrame(data=dict(name=snps['name'], NkiLad=in_lad))
            # Rearrange column order
            overlap = overlap[['name', 'NkiLad']]
            return overlap

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', index=False, header=True)
