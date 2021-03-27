import os
import pandas as pd
import numpy
from sys_tool import find_directory
from abstract_feature_util import AbstractFeatureUtil


class CaddUtil(AbstractFeatureUtil):
    def __init__(self):
        super(CaddUtil, self).__init__()
        self._cadd_dfm = None

    @property
    def cadd_dfm(self):
        # Lazy Initialization
        if self._cadd_dfm is None:
            # First row of `1000G_phase3.tsv` is a comment starting with "##".
            # However, `pandas.read_csv` cannot set `comment` to a string whose length > 1
            # If you set `comment="#"`, the header line starting with "#Chrom" will also be skipped
            # So just skip the first row
            _cadd_dfm = pd.read_csv(os.path.join(self.src_data_dir, self.src_data_fn),
                                    header=0, sep='\t', skiprows=1,
                                    dtype={"#Chrom": str, "Pos": numpy.int64, "Ref": str, "Alt": str,
                                           "RawScore": numpy.float64, "PHRED": str})

            # CADD uses 1-based coordinate system, so `Pos` in CADD is `chromEnd` in UCSC Genome Browser
            _cadd_dfm.rename(columns={"#Chrom": 'chrom', "Pos": 'chromEnd'}, inplace=True)
            # Adapt UCSC Genome Browser's representation of `chrom`
            _cadd_dfm.loc[:, 'chrom'] = _cadd_dfm.loc[:, 'chrom'].apply(lambda x: "chr{}".format(x))

            self._cadd_dfm = _cadd_dfm

        return self._cadd_dfm

    def get_feat(self, _input):
        """
        :param _input: SNP data frame in BED format
        :return:
        """
        # Make a copy; so setting an index won't affect the the input data frame
        snp_bed_dfm = _input.loc[:, ['chrom', 'chromStart', 'chromEnd', 'name']]

        result_dfm = snp_bed_dfm.merge(self.cadd_dfm, how='left', on=['chrom', 'chromEnd'])

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
        return result_dfm.loc[flag, ['name', 'CaddRawScore', 'CaddPhred']].fillna(0)

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', index=False, header=True)


if __name__ == '__main__':
    rsnp_dfm = pd.read_table(os.path.join(find_directory('CADD'), "RSNP_50kb.bed"), header=None,
                             names=['chrom', 'chromStart', 'chromEnd', 'name'])

    csnp_dfm = pd.read_table(os.path.join(find_directory('CADD'), "CSNP_50kb.bed"), header=None,
                             names=['chrom', 'chromStart', 'chromEnd', 'name'])

    cadd_util = CaddUtil()
    cadd_util.src_data_dir = find_directory('CADD')
    cadd_util.src_data_fn = "1000G_phase3.tsv"

    cadd_util.temp_dest = "./foo_rsnp.tsv"
    cadd_util.extract(_input=rsnp_dfm)

    cadd_util.temp_dest = "./foo_csnp.tsv"
    cadd_util.extract(_input=csnp_dfm)
