import os
from io import StringIO
import pandas as pd
from pybedtools import BedTool
from abstract_feature_util import AbstractFeatureUtil
from genome_browser_client import GenomeBrowserClient
from sys_tool import find_directory
from tss_tool import min_tss_dist


class ChiapetUtil(AbstractFeatureUtil):
    def __yield_feat(self, snp_bed_fn):
        snp_bed_obj = BedTool(snp_bed_fn)
        snp_bed_dfm = pd.read_table(snp_bed_fn, header=None, names=['chrom', 'chromStart', 'chromEnd', 'name'])

        for fn in self.src_data_fn:
            complement_filename = os.path.join(self.src_data_dir, fn)

            comp_bed_obj = BedTool(complement_filename)
            intx = snp_bed_obj.intersect(comp_bed_obj, wb=True)

            intx_dfm = pd.read_table(StringIO(str(intx)), header=None,
                                     names=['snpChrom', 'snpChromStart', 'snpChromEnd', 'snpName',
                                            'blockChrom', 'blockChromStart', 'blockChromEnd',
                                            'compChrom', 'compChromStart', 'compChromEnd'])
            snp_comp = intx_dfm[['snpName', 'compChrom', 'compChromStart', 'compChromEnd']]

            with GenomeBrowserClient(self.db_config_key) as gb_client:
                tss_dist = snp_comp.apply(lambda row: gb_client.select_tss_dist(row['compChrom'],
                                                                                row['compChromStart'],
                                                                                row['compChromEnd']),
                                          axis=1, reduce=True)
                tss_dist = pd.DataFrame(tss_dist, columns=['tssDist'])
                tss_dist = pd.concat([snp_comp['snpName'], tss_dist], axis=1)
                tss_dist = tss_dist.groupby('snpName').agg(sum).reset_index()

                col_name = fn[:-4] + 'TssDist' if fn.endswith(".bed") else fn + 'TssDist'

                # fn[:-4] remove the ".bed" extension
                tss_dist.loc[:, col_name] = tss_dist.loc[:, 'tssDist'].apply(min_tss_dist)

                result_dfm = pd.merge(snp_bed_dfm, tss_dist, how='left', left_on='name', right_on='snpName'). \
                    set_index('name')

                yield result_dfm[col_name]

    def get_feat(self, _input):
        """
        :param _input: the path of the SNP bed file
        :return:
        """
        yielded_dfms = list(self.__yield_feat(_input))
        result_dfm = pd.concat(yielded_dfms, axis=1)

        return result_dfm.reset_index()

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', index=False, header=True)

# if __name__ == '__main__':
#
#     rsnp_bed_fn = os.path.join(find_directory('chia_pet'), "RSNP_50kb.bed")
#     csnp_bed_fn = os.path.join(find_directory('chia_pet'), "CSNP_50kb.bed")
#
#     cp_util = ChiapetUtil()
#     cp_util.db_config_key = 'local_hg19'
#     cp_util.src_data_dir = find_directory('chia_pet')
#     cp_util.src_data_fn = [
#         "chiaPetHct116Pol2Rep1.bed",
#         "chiaPetHelas3Pol2Rep1.bed",
#         "chiaPetK562CtcfRep1.bed",
#         "chiaPetK562Pol2Rep1.bed",
#         "chiaPetK562Pol2Rep2.bed",
#         "chiaPetMcf7CtcfRep1.bed",
#         "chiaPetMcf7CtcfRep2.bed",
#         "chiaPetMcf7EraaRep1.bed",
#         "chiaPetMcf7EraaRep2.bed",
#         "chiaPetMcf7EraaRep3.bed",
#         "chiaPetMcf7Pol2Rep1.bed",
#         "chiaPetMcf7Pol2Rep2.bed",
#         "chiaPetMcf7Pol2Rep3.bed",
#         "chiaPetMcf7Pol2Rep4.bed",
#         "chiaPetNb4Pol2Rep1.bed"
#     ]
#
#     cp_util.temp_dest = "./foo_rsnp.tsv"
#     cp_util.extract(_input=rsnp_bed_fn)
#
#     cp_util.temp_dest = "./foo_csnp.tsv"
#     cp_util.extract(_input=csnp_bed_fn)
