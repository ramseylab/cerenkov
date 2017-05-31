"""
The following code block is loosely based on a code block in the GWAVA software version 1.0 (which is Copyright EBI)
source file "gwava_annotate.py".
"""
import os
import pandas as pd
import sys_tool
from pybedtools import BedTool
from abstract_feature_util import AbstractFeatureUtil


class GwavaUtil(AbstractFeatureUtil):
    def get_feat(self, _input):
        snp_dfm = _input.loc[:, ['chrom', 'chromStart', 'chromEnd', 'name']]
        snp_bed_obj = BedTool(snp_dfm.to_string(index=False, header=False, index_names=False), from_string=True)

        seg_bed_fn = os.path.join(self.src_data_dir, self.src_data_fn)
        seg_bed_obj = BedTool(seg_bed_fn)

        results = {}

        # The 'intersect' operation is not 'left-join' style so its result might have less entries than the SNP bed
        intersection = snp_bed_obj.intersect(seg_bed_obj, wb=True)
        if len(intersection) > 0:
            annots = intersection.groupby(g=[1, 2, 3, 4], c=8, o='collapse')
            for entry in annots:
                results[entry.name] = pd.Series(entry[4].split(',')).value_counts()

        names = {
            'CTCF': 'CTCF_REG',
            'E': 'ENH',
            'PF': 'TSS_FLANK',
            'R': 'REP',
            'T': 'TRAN',
            'TSS': 'TSS',
            'WE': 'WEAK_ENH'
        }

        gwava_dfm = pd.DataFrame(results, index=names.keys()).T.rename(columns=names)

        snp_dfm = snp_dfm.merge(gwava_dfm, how='left', left_on='name', right_index=True, copy=True)

        return snp_dfm.fillna(0).drop(['chrom', 'chromStart', 'chromEnd'], axis=1)

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', index=False, header=True)


if __name__ == '__main__':
    snp_bed_fn = os.path.join(sys_tool.find_directory("fsu_repli_chip"), 'RSNP_50kb.bed')
    snp_dfm = pd.read_table(snp_bed_fn, header=None,
                            names=['chrom', 'chromStart', 'chromEnd', 'name'])

    gwava_util = GwavaUtil()
    gwava_util.src_data_dir = sys_tool.find_directory('GWAVA')
    gwava_util.src_data_fn = 'segmentation.bed.gz'
    gwava_util.temp_dest = 'GWAVA_RSNP.txt'

    gwava_util.extract(_input=snp_dfm)
