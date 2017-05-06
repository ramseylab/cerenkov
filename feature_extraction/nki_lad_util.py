from genome_browser_client import GenomeBrowserClient
from abstract_feature_util import AbstractFeatureUtil
import pandas


# def extract_lad_overlap(src_bed, dest):
#     snps = pandas.read_csv(src_bed, sep='\t', header=None, names=['chrom', 'chromStart', 'chromEnd', 'name'])
#
#     with GenomeBrowserClient('local_hg19') as gb_client:
#         in_lad = snps.apply(lambda x: gb_client.in_nki_lad(x['chrom'], x['chromStart']), axis=1)
#
#         overlap = pandas.DataFrame(data=dict(name=snps['name'], NkiLad=in_lad))
#         # Rearrange column order
#         overlap = overlap[['name', 'NkiLad']]
#         overlap.to_csv(dest, sep='\t', index=False, header=True)


class NkiLadUtil(AbstractFeatureUtil):
    def get_feat(self, _input):
        """

        :param _input: the SNP data frame in BED format
        :return:
        """
        snps = _input

        with GenomeBrowserClient(self.db_config_key) as gb_client:
            in_lad = snps.apply(lambda x: gb_client.in_nki_lad(x['chrom'], x['chromStart']), axis=1)

            overlap = pandas.DataFrame(data=dict(name=snps['name'], NkiLad=in_lad))
            # Rearrange column order
            overlap = overlap[['name', 'NkiLad']]
            return overlap

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', index=False, header=True)


# if __name__ == '__main__':
#     extract_lad_overlap('workspace/RSNP_50kb.bed', 'workspace/RSNP_50kb_NKI_LAD.tsv')
#     extract_lad_overlap('workspace/CSNP_50kb.bed', 'workspace/CSNP_50kb_NKI_LAD.tsv')