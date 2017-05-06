from io import StringIO
import pandas
from pybedtools import BedTool
from abstract_feature_util import AbstractFeatureUtil
from genome_browser_client import GenomeBrowserClient
from sys_tool import find_directory
from tss_tool import min_tss_dist


class ChiapetUtil(AbstractFeatureUtil):
    def __yield_feat(self, snp_bed_fn):
        snp_bed_obj = BedTool(snp_bed_fn)
        snp_bed_dfm = pandas.read_table(snp_bed_fn, header=None, names=['chrom', 'chromStart', 'chromEnd', 'name'])

        for fn in self.src_data_fn:
            complement_filename = "{dir}/{fn}".format(dir=self.src_data_dir, fn=fn)

            comp_bed_obj = BedTool(complement_filename)
            intx = snp_bed_obj.intersect(comp_bed_obj, wb=True)

            intx_dfm = pandas.read_table(StringIO(str(intx)), header=None,
                                         names=['snpChrom', 'snpChromStart', 'snpChromEnd', 'snpName',
                                                'blockChrom', 'blockChromStart', 'blockChromEnd',
                                                'compChrom', 'compChromStart', 'compChromEnd'])
            snp_comp = intx_dfm[['snpName', 'compChrom', 'compChromStart', 'compChromEnd']]

            with GenomeBrowserClient(self.db_config_key) as gb_client:
                tss_dist = snp_comp.apply(lambda row: gb_client.select_tss_dist(row['compChrom'],
                                                                                row['compChromStart'],
                                                                                row['compChromEnd']),
                                          axis=1, reduce=True)
                tss_dist = pandas.DataFrame(tss_dist, columns=['tssDist'])
                tss_dist = pandas.concat([snp_comp['snpName'], tss_dist], axis=1)
                tss_dist = tss_dist.groupby('snpName').agg(sum).reset_index()

                col_name = fn[:-4] + 'TssDist' if fn.endswith(".bed") else fn + 'TssDist'

                # fn[:-4] remove the ".bed" extension
                tss_dist.loc[:, col_name] = tss_dist.loc[:, 'tssDist'].apply(min_tss_dist)

                result_dfm = pandas.merge(snp_bed_dfm, tss_dist, how='left', left_on='name', right_on='snpName'). \
                    set_index('name')

                yield result_dfm[col_name]

    def get_feat(self, _input):
        """
        :param _input: the path of the SNP bed file
        :return:
        """
        yielded_dfms = list(self.__yield_feat(_input))
        result_dfm = pandas.concat(yielded_dfms, axis=1)

        return result_dfm.reset_index()

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', index=False, header=True)

if __name__ == '__main__':
    # rsnp_dfms = list(gen_complement_tss_dist('RSNP'))
    # rsnp_dfm = pandas.concat(rsnp_dfms, axis=1)
    # rsnp_dfm.to_csv("{}/{}.tsv".format(find_directory("chia_pet"), 'RSNP_50kb_ChiA_PET'), sep='\t')
    #
    # csnp_dfms = list(gen_complement_tss_dist('CSNP'))
    # csnp_dfm = pandas.concat(csnp_dfms, axis=1)
    # csnp_dfm.to_csv("{}/{}.tsv".format(find_directory("chia_pet"), 'CSNP_50kb_ChiA_PET'), sep='\t')

    rsnp_bed_fn = "{}/RSNP_50kb.bed".format(find_directory("chia_pet"))
    csnp_bed_fn = "{}/CSNP_50kb.bed".format(find_directory("chia_pet"))

    cp_util = ChiapetUtil()
    cp_util.db_config_key = 'local_hg19'
    cp_util.src_data_dir = find_directory('chia_pet')
    cp_util.src_data_fn = [
        "chiaPetHct116Pol2Rep1.bed",
        "chiaPetHelas3Pol2Rep1.bed",
        "chiaPetK562CtcfRep1.bed",
        "chiaPetK562Pol2Rep1.bed",
        "chiaPetK562Pol2Rep2.bed",
        "chiaPetMcf7CtcfRep1.bed",
        "chiaPetMcf7CtcfRep2.bed",
        "chiaPetMcf7EraaRep1.bed",
        "chiaPetMcf7EraaRep2.bed",
        "chiaPetMcf7EraaRep3.bed",
        "chiaPetMcf7Pol2Rep1.bed",
        "chiaPetMcf7Pol2Rep2.bed",
        "chiaPetMcf7Pol2Rep3.bed",
        "chiaPetMcf7Pol2Rep4.bed",
        "chiaPetNb4Pol2Rep1.bed"
    ]

    cp_util.temp_dest = "./foo_rsnp.tsv"
    cp_util.extract(_input=rsnp_bed_fn)

    cp_util.temp_dest = "./foo_csnp.tsv"
    cp_util.extract(_input=csnp_bed_fn)


# import math
# __interactions = dict(
#     chiaPetHct116Pol2Rep1="wgEncodeGisChiaPetHct116Pol2InteractionsRep1",
#     chiaPetHelas3Pol2Rep1="wgEncodeGisChiaPetHelas3Pol2InteractionsRep1",
#     chiaPetK562CtcfRep1="wgEncodeGisChiaPetK562CtcfInteractionsRep1",
#     chiaPetK562Pol2Rep1="wgEncodeGisChiaPetK562Pol2InteractionsRep1",
#     chiaPetK562Pol2Rep2="wgEncodeGisChiaPetK562Pol2InteractionsRep2",
#     chiaPetMcf7CtcfRep1="wgEncodeGisChiaPetMcf7CtcfInteractionsRep1",
#     chiaPetMcf7CtcfRep2="wgEncodeGisChiaPetMcf7CtcfInteractionsRep2",
#     chiaPetMcf7EraaRep1="wgEncodeGisChiaPetMcf7EraaInteractionsRep1",
#     chiaPetMcf7EraaRep2="wgEncodeGisChiaPetMcf7EraaInteractionsRep2",
#     chiaPetMcf7EraaRep3="wgEncodeGisChiaPetMcf7EraaInteractionsRep3",
#     chiaPetMcf7Pol2Rep1="wgEncodeGisChiaPetMcf7Pol2InteractionsRep1",
#     chiaPetMcf7Pol2Rep2="wgEncodeGisChiaPetMcf7Pol2InteractionsRep2",
#     chiaPetMcf7Pol2Rep3="wgEncodeGisChiaPetMcf7Pol2InteractionsRep3",
#     chiaPetMcf7Pol2Rep4="wgEncodeGisChiaPetMcf7Pol2InteractionsRep4",
#     chiaPetNb4Pol2Rep1="wgEncodeGisChiaPetNb4Pol2InteractionsRep1",
# )
#
#
# def save_chia_pet_complement():
#     with GenomeBrowserClient('local_hg19') as gb_client:
#         for key, table_name in __interactions.items():
#             cluster = gb_client.select_chia_pet_cluster(table_name)
#
#             b1_complement_chrom = cluster['b2_chrom']
#             b1_complement_chrom_start = cluster[['b2_chromStart', 'b2_chromStart']].mean(axis=1).apply(math.floor)
#             b1_complement_chrom_end = b1_complement_chrom_start + 1
#             b1_complement = pandas.concat([cluster['b1_chrom'], cluster['b1_chromStart'],
#                                            cluster['b1_chromEnd'], b1_complement_chrom,
#                                            b1_complement_chrom_start, b1_complement_chrom_end], axis=1)
#             b1_complement.columns = ['chrom', 'chromStart', 'chromEnd',
#                                      'compChrom', 'compChromStart', 'compChromEnd']
#
#             b2_complement_chrom = cluster['b1_chrom']
#             b2_complement_chrom_start = cluster[['b1_chromStart', 'b1_chromStart']].mean(axis=1).apply(math.floor)
#             b2_complement_chrom_end = b2_complement_chrom_start + 1
#             b2_complement = pandas.concat([cluster['b2_chrom'], cluster['b2_chromStart'],
#                                            cluster['b2_chromEnd'], b2_complement_chrom,
#                                            b2_complement_chrom_start, b2_complement_chrom_end], axis=1)
#             b2_complement.columns = ['chrom', 'chromStart', 'chromEnd',
#                                      'compChrom', 'compChromStart', 'compChromEnd']
#
#             print("{} cluster.shape: {}".format(key, cluster.shape))
#             print("{} b1_comp.shape: {}".format(key, b1_complement.shape))
#             print("{} b2_comp.shape: {}".format(key, b2_complement.shape))
#
#             complement = pandas.concat([b1_complement, b2_complement], axis=0, ignore_index=True).\
#                 sort_values(by=['chrom', 'chromStart']).drop_duplicates()
#             print("{} comp.shape: {}".format(key, complement.shape))
#
#             print("Saving {} bed".format(key))
#             filename = "{}/{}.bed".format(find_directory("chia_pet"), key)
#             complement.to_csv(filename, sep="\t", index=False, header=False)

# def gen_complement_tss_dist(snp_type):
#     snp_filename = "{}/{}_50kb.bed".format(find_directory("chia_pet"), snp_type)
#     snp_dfm = pandas.read_table(snp_filename, header=None, names=['chrom', 'chromStart', 'chromEnd', 'name'])
#     snp_bed = BedTool(snp_filename)
#
#     for chia_pet_key in __interactions.keys():
#         complement_filename = "{}/{}.bed".format(find_directory("chia_pet"), chia_pet_key)
#
#         comp = BedTool(complement_filename)
#         intx = snp_bed.intersect(comp, wb=True)
#
#         intx_dfm = pandas.read_table(StringIO(str(intx)), header=None,
#                                      names=['snpChrom', 'snpChromStart', 'snpChromEnd', 'snpName',
#                                             'blockChrom', 'blockChromStart', 'blockChromEnd',
#                                             'compChrom', 'compChromStart', 'compChromEnd'])
#         snp_comp = intx_dfm[['snpName', 'compChrom', 'compChromStart', 'compChromEnd']]
#
#         with GenomeBrowserClient('local_hg19') as gb_client:
#             tss_dist = snp_comp.apply(lambda row: gb_client.select_tss_dist(row['compChrom'],
#                                                                             row['compChromStart'],
#                                                                             row['compChromEnd']),
#                                       axis=1, reduce=True)
#             tss_dist = pandas.DataFrame(tss_dist, columns=['tssDist'])
#             tss_dist = pandas.concat([snp_comp['snpName'], tss_dist], axis=1)
#             tss_dist = tss_dist.groupby('snpName').agg(sum).reset_index()
#             tss_dist.loc[:, chia_pet_key + 'TssDist'] = tss_dist.loc[:, 'tssDist'].apply(min_tss_dist)
#
#             result_dfm = pandas.merge(snp_dfm, tss_dist, how='left', left_on='name', right_on='snpName').\
#                 set_index('name')
#
#             yield result_dfm[chia_pet_key + 'TssDist']
