import os
from io import StringIO
import pandas as pd
from pybedtools import BedTool
from sys_tool import find_directory
from abstract_feature_util import AbstractFeatureUtil
from genome_browser_client import GenomeBrowserClient


class UniformDhsUtil(AbstractFeatureUtil):
    def __init__(self):
        super(UniformDhsUtil, self).__init__()
        self._dhs_bed_obj = None

    @property
    def dhs_bed_obj(self):
        # Lazy Initialization
        if self._dhs_bed_obj is None:
            dhs_bed_fn = os.path.join(self.src_data_dir, self.src_data_fn)
            dhs_bed_obj = BedTool(dhs_bed_fn)

            self._dhs_bed_obj = dhs_bed_obj

        return self._dhs_bed_obj

    def get_feat(self, _input):
        """
        :param _input: the SNP data frame in BED format
        :return:
        """
        snp_dfm = _input.loc[:, ['chrom', 'chromStart', 'chromEnd', 'name']]

        # rsnp_bed_obj = BedTool(rsnp_bed_fn)
        snp_bed_obj = BedTool(snp_dfm.to_string(index=False, header=False, index_names=False), from_string=True)

        # loj == left outer join.
        # intersect(`loj=True`) will return -1 dhsScore if no match
        # Do not use loj here; will fillna(0) later
        intx = snp_bed_obj.intersect(self.dhs_bed_obj, wb=True)
        intx_dfm = pd.read_table(StringIO(str(intx)), header=None,
                                 names=['snpChrom', 'snpChromStart', 'snpChromEnd', 'snpName',
                                        'dhsChrom', 'dhsChromStart', 'dhsChromEnd', 'dhsName', 'dhsScore',
                                        'dhsStrand'])
        dhs_score_dfm = intx_dfm.groupby(by='snpName')['dhsScore'].agg({"uniformDhsScore": sum,
                                                                        "uniformDhsCount": len})
        # dhs_score_dfm = intx_dfm.groupby(by='snpName')['dhsScore'].agg(["sum"]).rename(columns={"uniformDhsScore": "sum"})
        # dhs_score_dfm = intx_dfm.groupby(by='snpName')['dhsScore'].agg(["len"]).rename(columns={"uniformDhsCount": "len"})


        # snp_dfm.name <=> dhs_score_dfm.snpName
        snp_dhs = snp_dfm.merge(dhs_score_dfm, how='left', left_on='name', right_index=True)
        snp_dhs.fillna(0, inplace=True)

        return snp_dhs.loc[:, ["name", "uniformDhsScore", "uniformDhsCount"]]

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', header=True, index=False,
                       columns=["name", "uniformDhsScore", "uniformDhsCount"])


class MasterDhsUtil(AbstractFeatureUtil):
    def get_feat(self, _input):
        """
        :param _input: the SNP data frame
        :return:
        """
        snp_dfm = _input

        with GenomeBrowserClient(self.db_config_key) as gb_client:
            dhs_df = gb_client.fetch_master_dhs(snp_dfm.loc[:, 'name'].tolist())

        snp_dfm = snp_dfm.merge(dhs_df, how='left', on=['name', 'chrom'], copy=True)

        return snp_dfm.loc[:, ["name", "masterDhsScore", "masterDhsCount"]].fillna(0)

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', header=True, index=False,
                       columns=["name", "masterDhsScore", "masterDhsCount"])


if __name__ == '__main__':
    rsnp_bed_fn = os.path.join(find_directory('dhs'), "RSNP_50kb.bed")
    rsnp_dfm = pd.read_table(rsnp_bed_fn, header=None,
                             names=['chrom', 'chromStart', 'chromEnd', 'name'])
    csnp_bed_fn = os.path.join(find_directory('dhs'), "CSNP_50kb.bed")
    csnp_dfm = pd.read_table(csnp_bed_fn, header=None,
                             names=['chrom', 'chromStart', 'chromEnd', 'name'])

    dhs_util = MasterDhsUtil()
    dhs_util.db_config_key = 'local_hg19'

    dhs_util.temp_dest = 'FOO_CSNP.tsv'
    dhs_util.extract(csnp_dfm)

    dhs_util.temp_dest = 'FOO_RSNP.tsv'
    dhs_util.extract(rsnp_dfm)

    # tf_util = UniformDhsUtil()
    # tf_util.src_data_dir = find_directory("dhs")
    # tf_util.src_data_fn = "UniformDnaseIHS"
    #
    # tf_util.temp_dest = 'foo_rsnp.txt'
    # tf_util.extract(_input=rsnp_dfm)
    #
    # tf_util.temp_dest = 'foo_csnp.txt'
    # tf_util.extract(_input=csnp_dfm)


