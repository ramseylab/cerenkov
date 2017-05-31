import pandas as pd
from genome_browser_client import GenomeBrowserClient
from abstract_feature_util import AbstractFeatureUtil


class VistaEnhancerUtil(AbstractFeatureUtil):
    def get_feat(self, _input):
        snps = _input.loc[:, ['chrom', 'chromStart', 'chromEnd', 'name']]

        with GenomeBrowserClient(self.db_config_key) as gb_client:
            vh_list_series = snps.apply(lambda x: gb_client.select_vista_enhancer(x['chrom'], x['chromStart']), axis=1)
            vh_count = [len(vh_list) for vh_list in vh_list_series]
            vh_score = [sum([score for _, score in vh_list]) for vh_list in vh_list_series]

            vh_dfm = pd.DataFrame(data=dict(name=snps['name'],
                                            vistaEnhancerCnt=vh_count,
                                            vistaEnhancerTotalScore=vh_score))
            return vh_dfm

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', index=False, header=True)

