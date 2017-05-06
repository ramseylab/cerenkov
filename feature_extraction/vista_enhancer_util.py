import pandas
from genome_browser_client import GenomeBrowserClient
from abstract_feature_util import AbstractFeatureUtil


# def extract_vista_enhancer(src_bed, dest):
#     snps = pandas.read_csv(src_bed, sep='\t', header=None, names=['chrom', 'chromStart', 'chromEnd', 'name'])
#
#     with GenomeBrowserClient('local_hg19') as gb_client:
#         vh_list_series = snps.apply(lambda x: gb_client.select_vista_enhancer(x['chrom'], x['chromStart']), axis=1)
#         vh_count = [len(vh_list) for vh_list in vh_list_series]
#         vh_score = [sum([score for _, score in vh_list]) for vh_list in vh_list_series]
#
#         vh_dfm = pandas.DataFrame(data=dict(name=snps['name'],
#                                             vistaEnhancerCnt=vh_count,
#                                             vistaEnhancerTotalScore=vh_score))
#         vh_dfm.to_csv(dest, sep='\t', index=False, header=True)


class VistaEnhancerUtil(AbstractFeatureUtil):
    def get_feat(self, _input):
        snps = _input

        with GenomeBrowserClient(self.db_config_key) as gb_client:
            vh_list_series = snps.apply(lambda x: gb_client.select_vista_enhancer(x['chrom'], x['chromStart']), axis=1)
            vh_count = [len(vh_list) for vh_list in vh_list_series]
            vh_score = [sum([score for _, score in vh_list]) for vh_list in vh_list_series]

            vh_dfm = pandas.DataFrame(data=dict(name=snps['name'],
                                                vistaEnhancerCnt=vh_count,
                                                vistaEnhancerTotalScore=vh_score))
            return vh_dfm

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', index=False, header=True)

# if __name__ == '__main__':
#     extract_vista_enhancer('workspace/RSNP_50kb.bed', 'workspace/RSNP_50kb_vista_enhancer.tsv')
#     extract_vista_enhancer('workspace/CSNP_50kb.bed', 'workspace/CSNP_50kb_vista_enhancer.tsv')
