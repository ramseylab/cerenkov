import pandas as pd
from io import StringIO
import sys
import os
import sys_tool
from pybedtools import BedTool
from abstract_feature_util import AbstractFeatureUtil

if sys.version_info[0] == 3:
    # Using Python3
    from functools import reduce
# `reduce` is a built-in function in Python2


class RepliChipUtil(AbstractFeatureUtil):
    def __yield_score_dfm(self, snp_dfm):
        snp_bed_obj = BedTool(snp_dfm.to_string(index=False, header=False, index_names=False), from_string=True)

        for key, bed_fn in self.src_data_fn.items():
            rep_bed_fn = os.path.join(self.src_data_dir, bed_fn)
            rep_bed_obj = BedTool(rep_bed_fn)

            # Downstream scores
            closest_iu = snp_bed_obj.closest(rep_bed_obj, D='ref', iu=True)
            closest_iu_dfm = pd.read_table(StringIO(str(closest_iu)), header=None,
                                           names=['snpChrom', 'snpChromStart', 'snpChromEnd', 'snpName',
                                                  'repChrom', 'repChromStart', 'repChromEnd', 'repScore', 'distance'],
                                           usecols=['snpName', 'repScore'])
            closest_iu_dfm = closest_iu_dfm.rename(columns={'snpName': 'name',
                                                            'repScore': 'iu_score'})

            # Upstream scores
            closest_id = snp_bed_obj.closest(rep_bed_obj, D='ref', id=True)
            closest_id_dfm = pd.read_table(StringIO(str(closest_id)), header=None,
                                           names=['snpChrom', 'snpChromStart', 'snpChromEnd', 'snpName',
                                                  'repChrom', 'repChromStart', 'repChromEnd', 'repScore', 'distance'],
                                           usecols=['snpName', 'repScore'])
            closest_id_dfm = closest_id_dfm.rename(columns={'snpName': 'name',
                                                            'repScore': 'id_score'})

            # score_dfm = pd.concat([closest_iu_dfm, closest_id_dfm], axis=1)
            score_dfm = closest_iu_dfm.merge(closest_id_dfm, on='name')

            score_dfm = score_dfm.assign(avg_score=(score_dfm['iu_score'] + score_dfm['id_score']) / 2.0). \
                drop(['iu_score', 'id_score'], axis=1)
            score_dfm = score_dfm.rename(columns={'avg_score': key})

            yield score_dfm

    def get_feat(self, _input):
        snp_dfm = _input.loc[:, ['chrom', 'chromStart', 'chromEnd', 'name']]

        yielded_score_dfm = list(self.__yield_score_dfm(snp_dfm))
        # result = pd.concat(yielded_score_dfm, axis=1)
        result = reduce(lambda left, right: pd.merge(left, right, on='name'), yielded_score_dfm)
        return result

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', index=False, header=True)

if __name__ == '__main__':

    rc_util = RepliChipUtil()

    snp_bed_fn = os.path.join(sys_tool.find_directory("fsu_repli_chip"), 'RSNP_50kb.bed')
    snp_dfm = pd.read_table(snp_bed_fn, header=None, names=['chrom', 'chromStart', 'chromEnd', 'name'])

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

