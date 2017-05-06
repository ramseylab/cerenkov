import pandas
from pybedtools import BedTool
from io import StringIO
from sys_tool import find_directory
from abstract_feature_util import AbstractFeatureUtil


class FitconsUtil(AbstractFeatureUtil):
    def __init__(self):
        self._fitcons_bed_obj_dict = None

    @property
    def fitcons_bed_obj_dict(self):
        # Lazy Initialization
        if self._fitcons_bed_obj_dict is None:
            fitcons_bed_fns = ["{dir}/{fn}".format(dir=self.src_data_dir, fn=fn) for fn in self.src_data_fn.values()]
            fitcons_bed_objs = [BedTool(fn) for fn in fitcons_bed_fns]

            self._fitcons_bed_obj_dict = dict(zip(self.src_data_fn.keys(), fitcons_bed_objs))

            # print(self.src_data_fn.keys())
            # print(fitcons_bed_fns)

        return self._fitcons_bed_obj_dict

    def __yield_fitcons_dfm(self, snp_bed_obj):
        for key, fitcons_bed_obj in self.fitcons_bed_obj_dict.items():
            # loj == left outer join.
            # That is, for each feature in A, if no overlaps in B are found, report a NULL.
            # However, NULL in pybedtool is a dot ('.')
            # `to_numeric` is used below to coerce '.' into NaN
            intx = snp_bed_obj.intersect(fitcons_bed_obj, wb=True, loj=True)

            intx_dfm = pandas.read_table(StringIO(str(intx)), header=None,
                                         names=['snpChrom', 'snpChromStart', 'snpChromEnd', 'snpName',
                                                'fcChrom', 'fcChromStart', 'fcChromEnd', 'fcScore'])
            fitcons_dfm = intx_dfm.loc[:, ['snpName', 'fcScore']].copy()
            fitcons_dfm.loc[:, 'fcScore'] = pandas.to_numeric(fitcons_dfm.loc[:, 'fcScore'], errors='coerce')
            fitcons_dfm.rename(columns={'snpName': 'name', 'fcScore': key}, inplace=True)
            fitcons_dfm.set_index('name', inplace=True)

            yield fitcons_dfm

    def get_feat(self, _input):
        snp_bed_obj = BedTool(_input.to_string(index=False, header=False, index_names=False), from_string=True)
        yielded_dfms = list(self.__yield_fitcons_dfm(snp_bed_obj))
        snp_fitcons = pandas.concat(yielded_dfms, axis=1)
        return snp_fitcons.reset_index()

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', index=False, header=True)


# __bed_files = dict(
#     fitConsGm="fc-gm-0.bed",
#     fitConsH1="fc-h1-0.bed",
#     fitConsHu="fc-hu-0.bed",
#     fitConsI6="fc-i6-0.bed",
# )


# def gen_fitcons_score(snp_type):
#     snp_filename = "{}/{}_50kb.bed".format(find_directory("fitcons"), snp_type)
#     snp_bed = BedTool(snp_filename)
#
#     for key, file in __bed_files.items():
#         fitcons_file = "{}/{}".format(find_directory("fitcons"), file)
#
#         fitcons_bed = BedTool(fitcons_file)
#         # loj == left outer join.
#         # That is, for each feature in A, if no overlaps in B are found, report a NULL.
#         # However, NULL in pybedtool is a dot ('.')
#         # `to_numeric` is used below to coerce '.' into NaN
#         intx = snp_bed.intersect(fitcons_bed, wb=True, loj=True)
#
#         intx_dfm = pandas.read_table(StringIO(str(intx)), header=None,
#                                      names=['snpChrom', 'snpChromStart', 'snpChromEnd', 'snpName',
#                                             'fcChrom', 'fcChromStart', 'fcChromEnd', 'fcScore'])
#         fitcons_dfm = intx_dfm.loc[:, ['snpName', 'fcScore']].copy()
#         fitcons_dfm.loc[:, 'fcScore'] = pandas.to_numeric(fitcons_dfm.loc[:, 'fcScore'], errors='coerce')
#         fitcons_dfm.rename(columns={'snpName': 'name', 'fcScore': key}, inplace=True)
#         fitcons_dfm.set_index('name', inplace=True)
#
#         yield fitcons_dfm

if __name__ == '__main__':
    # rsnp_dfms = list(gen_fitcons_score('RSNP'))
    # rsnp_dfm = pandas.concat(rsnp_dfms, axis=1)
    # rsnp_dfm.to_csv("{}/{}.tsv".format(find_directory("fitcons"), 'RSNP_50kb_fitCons'), sep='\t')
    #
    # csnp_dfms = list(gen_fitcons_score('CSNP'))
    # csnp_dfm = pandas.concat(csnp_dfms, axis=1)
    # csnp_dfm.to_csv("{}/{}.tsv".format(find_directory("fitcons"), 'CSNP_50kb_fitCons'), sep='\t')
    fitcons_util = FitconsUtil()
    fitcons_util.src_data_dir = find_directory("fitcons")
    fitcons_util.src_data_fn = dict(
        fitConsGm="fc-gm-0.bed",
        fitConsH1="fc-h1-0.bed",
        fitConsHu="fc-hu-0.bed",
        fitConsI6="fc-i6-0.bed",
    )

    rsnp_dfm = pandas.read_table("{}/RSNP_50kb.bed".format(find_directory('fitcons')), header=None,
                                 names=['chrom', 'chromStart', 'chromEnd', 'name'])
    fitcons_util.temp_dest = 'foo_rsnp.txt'
    fitcons_util.extract(_input=rsnp_dfm)

    csnp_dfm = pandas.read_table("{}/CSNP_50kb.bed".format(find_directory('fitcons')), header=None,
                                 names=['chrom', 'chromStart', 'chromEnd', 'name'])
    fitcons_util.temp_dest = 'foo_csnp.txt'
    fitcons_util.extract(_input=csnp_dfm)