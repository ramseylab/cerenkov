import pandas
import sys_tool
from pybedtools import BedTool
from io import StringIO
from abstract_feature_util import AbstractFeatureUtil

# Extraction From Bigwigh is DEPRECATED!!!

# def extract_from_bigwig(src_bed, dest):
#     bigwig = dict(
#         FsuBg02esRep1="wgEncodeFsuRepliChipBg02esWaveSignalRep1.bigWig",
#         FsuBg02esRep2="wgEncodeFsuRepliChipBg02esWaveSignalRep2.bigWig",
#         FsuGm06990Rep1="wgEncodeFsuRepliChipGm06990WaveSignalRep1.bigWig",
#         FsuGm06990Rep2="wgEncodeFsuRepliChipGm06990WaveSignalRep2.bigWig",
#         FsuH1hescRep1="wgEncodeFsuRepliChipH1hescWaveSignalRep1.bigWig",
#         FsuH1hescRep2="wgEncodeFsuRepliChipH1hescWaveSignalRep2.bigWig",
#         FsuH1hescRep3="wgEncodeFsuRepliChipH1hescWaveSignalRep3.bigWig",
#         FsuH7esRep1="wgEncodeFsuRepliChipH7esWaveSignalRep1.bigWig",
#         FsuH7esRep2="wgEncodeFsuRepliChipH7esWaveSignalRep2.bigWig",
#         FsuH9esRep1="wgEncodeFsuRepliChipH9esWaveSignalRep1.bigWig",
#         FsuHelas3Rep1="wgEncodeFsuRepliChipHelas3WaveSignalRep1.bigWig",
#         FsuImr90Rep1="wgEncodeFsuRepliChipImr90WaveSignalRep1.bigWig",
#         FsuIpshfib2ips4Rep1="wgEncodeFsuRepliChipIpshfib2ips4WaveSignalRep1.bigWig",
#         FsuIpshfib2ips4Rep2="wgEncodeFsuRepliChipIpshfib2ips4WaveSignalRep2.bigWig",
#         FsuIpshfib2ips5Rep1="wgEncodeFsuRepliChipIpshfib2ips5WaveSignalRep1.bigWig",
#         FsuIpshfib2ips5Rep2="wgEncodeFsuRepliChipIpshfib2ips5WaveSignalRep2.bigWig"
#     )
#
#     def frc_dfm():
#         for key in bigwig.keys():
#             src_bw = "{dir}/{file}".format(dir=sys_tool.find_directory("fsu_repli_chip"), file=bigwig[key])
#
#             # This command would output a file with undesired columns. Need pruning.
#             sys_tool.run_bigWigAverageOverBed([src_bw, src_bed, dest])
#
#             """
#             The dest has no header. Its 6 columns are:
#                name - name field from bed, which should be unique
#                size - size of bed (sum of exon sizes)
#                covered - # bases within exons covered by bigWig
#                sum - sum of values over all bases covered
#                mean0 - average over bases with non-covered bases counting as zeroes
#                mean - average over just covered bases
#             """
#
#             col_names = ['name', 'size', 'covered', 'sum', 'mean0', 'mean']
#
#             # Use `dtype=str` or `dtype=object` to preserve and not interpret dtype
#             # If pandas reads scores as floats, the dest would contain float numbers with long tails.
#             # E.g. input == 0.806; dest == 0.8059999999999999
#             _frc_dfm = pandas.read_csv(dest, sep='\t', header=None, names=col_names, usecols=['name', 'sum'], dtype=str)
#
#             # Choose column `sum` to represent repli-chip scores
#             _frc_dfm = _frc_dfm.rename(columns={'sum': key}).set_index('name')
#
#             yield _frc_dfm
#
#     dfms = list(frc_dfm())
#     dfm = pandas.concat(dfms, axis=1)
#     dfm.to_csv(dest, sep='\t', index=True, header=True)
#
# if __name__ == '__main__':
#     # extract_fsu_repli_chip('workspace/RSNP_50kb.bed', 'workspace/RSNP_50kb_FSU.tsv')
#     extract_from_bigwig('workspace/CSNP_50kb.bed', 'workspace/CSNP_50kb_FSU.tsv')


class RepliChipUtil(AbstractFeatureUtil):
    def __yield_score_dfm(self, snp_dfm):
        snp_bed_obj = BedTool(snp_dfm.to_string(index=False, header=False, index_names=False), from_string=True)

        for key, bed_fn in self.src_data_fn.items():
            rep_bed_fn = "{dir}/{file}".format(dir=self.src_data_dir, file=bed_fn)
            rep_bed_obj = BedTool(rep_bed_fn)

            # Downstream scores
            closest_iu = snp_bed_obj.closest(rep_bed_obj, D='ref', iu=True)
            closest_iu_dfm = pandas.read_table(StringIO(str(closest_iu)), header=None,
                                               names=['snpChrom', 'snpChromStart', 'snpChromEnd', 'snpName',
                                                      'repChrom', 'repChromStart', 'repChromEnd', 'repScore',
                                                      'distance'],
                                               usecols=['snpName', 'repScore'])
            closest_iu_dfm = closest_iu_dfm.rename(columns={'snpName': 'name',
                                                            'repScore': 'iu_score'}).set_index('name')

            # Upstream scores
            closest_id = snp_bed_obj.closest(rep_bed_obj, D='ref', id=True)
            closest_id_dfm = pandas.read_table(StringIO(str(closest_id)), header=None,
                                               names=['snpChrom', 'snpChromStart', 'snpChromEnd', 'snpName',
                                                      'repChrom', 'repChromStart', 'repChromEnd', 'repScore',
                                                      'distance'],
                                               usecols=['snpName', 'repScore'])
            closest_id_dfm = closest_id_dfm.rename(columns={'snpName': 'name',
                                                            'repScore': 'id_score'}).set_index('name')

            score_dfm = pandas.concat([closest_iu_dfm, closest_id_dfm], axis=1)
            score_dfm = score_dfm.assign(avg_score=(score_dfm['iu_score'] + score_dfm['id_score']) / 2.0). \
                drop(['iu_score', 'id_score'], axis=1)
            score_dfm = score_dfm.rename(columns={'avg_score': key})

            yield score_dfm

    def get_feat(self, _input):
        snp_dfm = _input

        yielded_score_dfm = list(self.__yield_score_dfm(snp_dfm))
        result = pandas.concat(yielded_score_dfm, axis=1)
        return result.reset_index()

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', index=False, header=True)

if __name__ == '__main__':

    rc_util = RepliChipUtil()

    snp_bed_fn = "{dir}/{file}".format(dir=sys_tool.find_directory("fsu_repli_chip"),
                                       file='RSNP_50kb.bed')
    snp_dfm = pandas.read_table(snp_bed_fn, header=None,
                                names=['chrom', 'chromStart', 'chromEnd', 'name'])

    rc_util.src_data_dir = sys_tool.find_directory("fsu_repli_chip")
    rc_util.src_data_fn = dict(
        FsuBg02esRep1="wgEncodeFsuRepliChipBg02esWaveSignalRep1.bed",
        FsuBg02esRep2="wgEncodeFsuRepliChipBg02esWaveSignalRep2.bed",
        FsuGm06990Rep1="wgEncodeFsuRepliChipGm06990WaveSignalRep1.bed",
        FsuGm06990Rep2="wgEncodeFsuRepliChipGm06990WaveSignalRep2.bed",
        FsuH1hescRep1="wgEncodeFsuRepliChipH1hescWaveSignalRep1.bed",
        FsuH1hescRep2="wgEncodeFsuRepliChipH1hescWaveSignalRep2.bed",
        FsuH1hescRep3="wgEncodeFsuRepliChipH1hescWaveSignalRep3.bed",
        FsuH7esRep1="wgEncodeFsuRepliChipH7esWaveSignalRep1.bed",
        FsuH7esRep2="wgEncodeFsuRepliChipH7esWaveSignalRep2.bed",
        FsuH9esRep1="wgEncodeFsuRepliChipH9esWaveSignalRep1.bed",
        FsuHelas3Rep1="wgEncodeFsuRepliChipHelas3WaveSignalRep1.bed",
        FsuImr90Rep1="wgEncodeFsuRepliChipImr90WaveSignalRep1.bed",
        FsuIpshfib2ips4Rep1="wgEncodeFsuRepliChipIpshfib2ips4WaveSignalRep1.bed",
        FsuIpshfib2ips4Rep2="wgEncodeFsuRepliChipIpshfib2ips4WaveSignalRep2.bed",
        FsuIpshfib2ips5Rep1="wgEncodeFsuRepliChipIpshfib2ips5WaveSignalRep1.bed",
        FsuIpshfib2ips5Rep2="wgEncodeFsuRepliChipIpshfib2ips5WaveSignalRep2.bed",
    )
    rc_util.temp_dest = 'RC_FSU_RSNP.txt'
    fsu_result = rc_util.extract(snp_dfm)

    rc_util.src_data_dir = sys_tool.find_directory("uw_repli_chip")
    rc_util.src_data_fn = dict(
        UwBg02esRep1="wgEncodeUwRepliSeqBg02esWaveSignalRep1.bed",
        UwBjRep1="wgEncodeUwRepliSeqBjWaveSignalRep1.bed",
        UwBjRep2="wgEncodeUwRepliSeqBjWaveSignalRep2.bed",
        UwGm06990Rep1="wgEncodeUwRepliSeqGm06990WaveSignalRep1.bed",
        UwGm12801Rep1="wgEncodeUwRepliSeqGm12801WaveSignalRep1.bed",
        UwGm12812Rep1="wgEncodeUwRepliSeqGm12812WaveSignalRep1.bed",
        UwGm12813Rep1="wgEncodeUwRepliSeqGm12813WaveSignalRep1.bed",
        UwGm12878Rep1="wgEncodeUwRepliSeqGm12878WaveSignalRep1.bed",
        UwHelas3Rep1="wgEncodeUwRepliSeqHelas3WaveSignalRep1.bed",
        UwHepg2Rep1="wgEncodeUwRepliSeqHepg2WaveSignalRep1.bed",
        UwHuvecRep1="wgEncodeUwRepliSeqHuvecWaveSignalRep1.bed",
        UwImr90Rep1="wgEncodeUwRepliSeqImr90WaveSignalRep1.bed",
        UwK562Rep1="wgEncodeUwRepliSeqK562WaveSignalRep1.bed",
        UwMcf7Rep1="wgEncodeUwRepliSeqMcf7WaveSignalRep1.bed",
        UwNhekRep1="wgEncodeUwRepliSeqNhekWaveSignalRep1.bed",
        UwSknshRep1="wgEncodeUwRepliSeqSknshWaveSignalRep1.bed",
    )
    rc_util.temp_dest = 'UW_FSU_RSNP.txt'
    fsu_result = rc_util.extract(snp_dfm)

