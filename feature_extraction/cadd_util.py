import pandas
import numpy
from sys_tool import find_directory
from abstract_feature_util import AbstractFeatureUtil


class CaddUtil(AbstractFeatureUtil):
    def __init__(self):
        self._cadd_dfm = None

    @property
    def cadd_dfm(self):
        # Lazy Initialization
        if self._cadd_dfm is None:
            # First row of `1000G_phase3.tsv` is a comment starting with "##".
            # However, `pandas.read_csv` cannot set `comment` to a string whose length > 1
            # If you set `comment="#"`, the header line starting with "#Chrom" will also be skipped
            # So just skip the first row
            _cadd_dfm = pandas.read_csv("{dir}/{fn}".format(dir=self.src_data_dir, fn=self.src_data_fn),
                                        header=0, sep='\t', skiprows=1,
                                        dtype={"#Chrom": str, "Pos": numpy.int64, "Ref": str, "Alt": str,
                                               "RawScore": numpy.float64, "PHRED": str})

            # CADD uses 1-based coordinate system, so `Pos` in CADD is `chromEnd` in UCSC Genome Browser
            _cadd_dfm.rename(columns={"#Chrom": 'chrom', "Pos": 'chromEnd'}, inplace=True)
            # Adapt UCSC Genome Browser's representation of `chrom`
            _cadd_dfm.loc[:, 'chrom'] = _cadd_dfm.loc[:, 'chrom'].apply(lambda x: "chr{}".format(x))
            _cadd_dfm.set_index(['chrom', 'chromEnd'], inplace=True)

            self._cadd_dfm = _cadd_dfm

        return self._cadd_dfm

    def get_feat(self, _input):
        """
        :param _input: SNP data frame in BED format
        :return:
        """

        snp_bed_dfm = _input
        snp_bed_dfm.set_index(['chrom', 'chromEnd'], inplace=True)

        # Join on the multi-index
        result_dfm = snp_bed_dfm.join(self.cadd_dfm, how='left').reset_index()

        # In case of duplication, fetch the one the max `RawScore`
        # `Transform` means transforming each row's `RawScore` into the max values within this row's group
        #   It does not modify the `RawScore` column in-place but returns a new Series
        # `max_raw_score` have the same index of `result_dfm`
        max_raw_score = result_dfm.groupby('name')['RawScore'].transform(max)
        # max of NaN is still NaN
        # We need to treat NaN == NaN here
        flag = (max_raw_score == result_dfm.loc[:, 'RawScore']) | \
               (max_raw_score.isnull() & result_dfm.loc[:, 'RawScore'].isnull())

        result_dfm.rename(columns={'PHRED': 'CaddPhred', 'RawScore': 'CaddRawScore'}, inplace=True)
        return result_dfm.loc[flag, ['name', 'CaddRawScore', 'CaddPhred']]

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', index=False, header=True)


# def extract_phred(_snp_dfm, _cadd_dfm):
#     # Join on the multi-index
#     result_dfm = _snp_dfm.join(_cadd_dfm, how='left').reset_index()
#
#     # In case of duplication, fetch the one the max `RawScore`
#     # `Transform` means transforming each row's `RawScore` into the max values within this row's group
#     #   It does not modify the `RawScore` column in-place but returns a new Series
#     # `max_raw_score` have the same index of `result_dfm`
#     max_raw_score = result_dfm.groupby('name')['RawScore'].transform(max)
#     # max of NaN is still NaN
#     # We need to treat NaN == NaN here
#     flag = (max_raw_score == result_dfm.loc[:, 'RawScore']) | \
#            (max_raw_score.isnull() & result_dfm.loc[:, 'RawScore'].isnull())
#
#     result_dfm.rename(columns={'PHRED': 'CaddPhred', 'RawScore': 'CaddRawScore'}, inplace=True)
#     return result_dfm.loc[flag, ['name', 'CaddRawScore', 'CaddPhred']]

if __name__ == '__main__':
    # rsnp_dfm = pandas.read_table("{}/RSNP_50kb.bed".format(find_directory('CADD')), header=None,
    #                               names=['chrom', 'chromStart', 'chromEnd', 'name'])
    # rsnp_dfm.set_index(['chrom', 'chromEnd'], inplace=True)
    #
    # csnp_dfm = pandas.read_table("{}/CSNP_50kb.bed".format(find_directory('CADD')), header=None,
    #                               names=['chrom', 'chromStart', 'chromEnd', 'name'])
    # csnp_dfm.set_index(['chrom', 'chromEnd'], inplace=True)
    #
    # # First row of `1000G_phase3.tsv` is a comment starting with "##".
    # # However, `pandas.read_csv` cannot set `comment` to a string whose length > 1
    # # If you set `comment="#"`, the header line starting with "#Chrom" will also be skipped
    # # So just skip the first row
    # cadd_dfm = pandas.read_csv("{}/1000G_phase3.tsv".format(find_directory('CADD')), header=0, sep='\t', skiprows=1,
    #                            dtype={"#Chrom": str, "Pos": numpy.int64, "Ref": str, "Alt": str,
    #                                   "RawScore": numpy.float64, "PHRED": str})
    #
    # # CADD uses 1-based coordinate system, so `Pos` in CADD is `chromEnd` in UCSC Genome Browser
    # cadd_dfm.rename(columns={"#Chrom": 'chrom', "Pos": 'chromEnd'}, inplace=True)
    # # Adapt UCSC Genome Browser's representation of `chrom`
    # cadd_dfm.loc[:, 'chrom'] = cadd_dfm.loc[:, 'chrom'].apply(lambda x: "chr{}".format(x))
    # cadd_dfm.set_index(['chrom', 'chromEnd'], inplace=True)
    #
    # rsnp_result = extract_phred(rsnp_dfm, cadd_dfm)
    # rsnp_result.to_csv("{}/RSNP_50kb_CADD.tsv".format(find_directory('CADD')), sep='\t', index=False)
    #
    # csnp_result = extract_phred(csnp_dfm, cadd_dfm)
    # csnp_result.to_csv("{}/CSNP_50kb_CADD.tsv".format(find_directory('CADD')), sep='\t', index=False)

    rsnp_dfm = pandas.read_table("{}/RSNP_50kb.bed".format(find_directory('CADD')), header=None,
                                 names=['chrom', 'chromStart', 'chromEnd', 'name'])

    csnp_dfm = pandas.read_table("{}/CSNP_50kb.bed".format(find_directory('CADD')), header=None,
                                 names=['chrom', 'chromStart', 'chromEnd', 'name'])

    cadd_util = CaddUtil()
    cadd_util.src_data_dir = find_directory('CADD')
    cadd_util.src_data_fn = "1000G_phase3.tsv"

    cadd_util.temp_dest = "./foo_rsnp.tsv"
    cadd_util.extract(_input=rsnp_dfm)

    cadd_util.temp_dest = "./foo_csnp.tsv"
    cadd_util.extract(_input=csnp_dfm)
