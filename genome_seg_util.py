from genome_browser_client import GenomeBrowserClient
from abstract_feature_util import AbstractFeatureUtil


# def extract_genome_seg_annot(src, dest):
#     """
#         Match every SNP to GSs by rsid;
#
#         E.g. A matched group for SNP `rs1` is like below
#
#             name    seg
#             rs1        ['Repr','Quies']
#     """
#
#     rsid = pandas.read_csv(src, sep='\t').loc[:, "name"]
#
#     with GenomeBrowserClient('local_hg19') as gb_client:
#         result = gb_client.identify_genome_seg(rsid)
#
#         result = CT.remove_dup_on_chrY(result)
#
#         result.to_csv(dest, sep='\t', header=True, index=False)


class GenomeSegUtil(AbstractFeatureUtil):
    def __init__(self, reproduce_osu17=False):
        super(GenomeSegUtil, self).__init__()
        self.reproduce_osu17 = reproduce_osu17

    def get_feat(self, _input):
        """

        :param _input: the SNP data frame
        :return:
        """
        snp_dfm = _input.loc[:, ['chrom', 'name']]

        with GenomeBrowserClient(self.db_config_key) as gb_client:
            result = gb_client.identify_genome_seg(snp_dfm.loc[:, 'name'])

            # result = ct.remove_dup_on_chrY(result)

            if not self.reproduce_osu17:
                # Use clearer names for osu18
                result = result.rename(columns={'ch1Name': 'ChromhmmGm12878',
                                                'ch2Name': 'ChromhmmH1hesc',
                                                'ch3Name': 'ChromhmmHelas3',
                                                'ch4Name': 'ChromhmmHepg2',
                                                'ch5Name': 'ChromhmmHuvec',
                                                'ch6Name': 'ChromhmmK562',
                                                'sw1Name': 'SegwayGm12878',
                                                'sw2Name': 'SegwayH1hesc',
                                                'sw3Name': 'SegwayHelas3',
                                                'sw4Name': 'SegwayHepg2',
                                                'sw5Name': 'SegwayHuvec',
                                                'sw6Name': 'SegwayK562'})

            snp_dfm = snp_dfm.merge(result, how='left', on=['name', 'chrom'])

            return snp_dfm.drop(['chrom'], axis=1).fillna(0)

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', header=True, index=False)
