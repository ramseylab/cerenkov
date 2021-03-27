import sys
import os
import pandas as pd
from pybedtools import BedTool
from io import StringIO
from sys_tool import find_directory
from abstract_feature_util import AbstractFeatureUtil

if sys.version_info[0] == 3:
    # Using Python3
    from functools import reduce
# `reduce` is a built-in function in Python2


class FitconsUtil(AbstractFeatureUtil):
    def __init__(self):
        super(FitconsUtil, self).__init__()
        self._fitcons_bed_obj_dict = None

    @property
    def fitcons_bed_obj_dict(self):
        # Lazy Initialization
        if self._fitcons_bed_obj_dict is None:
            fitcons_bed_fns = [os.path.join(self.src_data_dir, fn) for fn in self.src_data_fn.values()]
            fitcons_bed_objs = [BedTool(fn) for fn in fitcons_bed_fns]

            self._fitcons_bed_obj_dict = dict(zip(self.src_data_fn.keys(), fitcons_bed_objs))

        return self._fitcons_bed_obj_dict

    def __yield_fitcons_dfm(self, snp_bed_obj):
        for key, fitcons_bed_obj in self.fitcons_bed_obj_dict.items():
            # loj == left outer join.
            # That is, for each feature in A, if no overlaps in B are found, report a NULL.
            # However, NULL in pybedtool is a dot ('.')
            # `to_numeric` is used below to coerce '.' into NaN
            intx = snp_bed_obj.intersect(fitcons_bed_obj, wb=True, loj=True)

            intx_dfm = pd.read_table(StringIO(str(intx)), header=None,
                                     names=['snpChrom', 'snpChromStart', 'snpChromEnd', 'snpName',
                                            'fcChrom', 'fcChromStart', 'fcChromEnd', 'fcScore'])
            fitcons_dfm = intx_dfm.loc[:, ['snpName', 'fcScore']].copy()
            fitcons_dfm.loc[:, 'fcScore'] = pd.to_numeric(fitcons_dfm.loc[:, 'fcScore'], errors='coerce').fillna(0)
            fitcons_dfm.rename(columns={'snpName': 'name', 'fcScore': key}, inplace=True)

            yield fitcons_dfm

    def get_feat(self, _input):
        snp_bed_dfm = _input.loc[:, ['chrom', 'chromStart', 'chromEnd', 'name']]
        snp_bed_obj = BedTool(snp_bed_dfm.to_string(index=False, header=False, index_names=False), from_string=True)
        yielded_dfms = list(self.__yield_fitcons_dfm(snp_bed_obj))

        snp_fitcons = reduce(lambda left, right: pd.merge(left, right, on='name'), yielded_dfms)

        return snp_fitcons

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', index=False, header=True)


if __name__ == '__main__':
    fitcons_util = FitconsUtil()
    fitcons_util.src_data_dir = find_directory("fitcons")
    fitcons_util.src_data_fn = dict(
        fitConsGm="fc-gm-0.bed",
        fitConsH1="fc-h1-0.bed",
        fitConsHu="fc-hu-0.bed",
        fitConsI6="fc-i6-0.bed",
    )

    rsnp_dfm = pd.read_table(os.path.join(find_directory('fitcons'), "RSNP_50kb.bed"), header=None,
                             names=['chrom', 'chromStart', 'chromEnd', 'name'])
    fitcons_util.temp_dest = 'foo_rsnp.txt'
    fitcons_util.extract(_input=rsnp_dfm)

    csnp_dfm = pd.read_table(os.path.join(find_directory('fitcons'), "CSNP_50kb.bed"), header=None,
                             names=['chrom', 'chromStart', 'chromEnd', 'name'])
    fitcons_util.temp_dest = 'foo_csnp.txt'
    fitcons_util.extract(_input=csnp_dfm)
