from abstract_feature_util import AbstractFeatureUtil
from genome_browser_client import GenomeBrowserClient
from chrom_tool import remove_dup_on_chrY, CHR_LENGTH


class CoordUtil(AbstractFeatureUtil):
    @staticmethod
    def pce_filter(snp_dfm, verbose=False):
        with GenomeBrowserClient('local_hg19') as gb_client:
            in_pce = snp_dfm.apply(lambda x: gb_client.in_protein_coding_exon(x['chrom'], x['chromStart']), axis=1)

        if verbose:
            print("[pce_filter] found {n} SNP in protein-coding exons: \r\n {dfm}".
                  format(n=sum(in_pce), dfm=snp_dfm[in_pce]))

        return snp_dfm.loc[~in_pce]

    @staticmethod
    def add_norm_coord(snp_dfm):
        # Integer division will yield float
        # See PEP 238 -- Changing the Division Operator, https://www.python.org/dev/peps/pep-0238/
        snp_dfm.loc[:, "normChromCoord"] = snp_dfm.apply(lambda row: row['chromStart'] / CHR_LENGTH[row['chrom']],
                                                         axis=1)

        return snp_dfm

    def get_feat(self, _input):
        rsid = _input

        with GenomeBrowserClient(self.db_config_key) as gb_client:
            coord_dfm = gb_client.fetch_coord(rsid)
            coord_dfm = remove_dup_on_chrY(coord_dfm)
            return coord_dfm

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', index=False, header=True)
