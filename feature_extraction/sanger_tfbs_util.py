import os
import pandas as pd
import sys_tool
import tempfile
from abstract_feature_util import AbstractFeatureUtil


class SangerTfbsUtil(AbstractFeatureUtil):
    def get_feat(self, _input):
        """

        :param _input: the path to the SNP BED file
        :return:
        """
        src_bed = _input

        src_bw = os.path.join(self.src_data_dir, self.src_data_fn)

        with tempfile.NamedTemporaryFile(prefix='cerenkov-') as dest:
            # This command would output a file with undesired columns. Need pruning.
            sys_tool.run_bigWigAverageOverBed([src_bw, src_bed, dest.name])

            """
            The dest has no header. Its 6 columns are:
               name - name field from bed, which should be unique
               size - size of bed (sum of exon sizes
               covered - # bases within exons covered by bigWig
               sum - sum of values over all bases covered
               mean0 - average over bases with non-covered bases counting as zeroes
               mean - average over just covered bases
            """
            # Actually the last 3 columns are identical in our situation
            # Choose column `sum` to represent gerp scores
            col_names = ['name', 'size', 'covered', 'sum', 'mean0', 'mean']

            # Use `dtype=str` or `dtype=object` to preserve and not interpret dtype
            # If pandas reads scores as floats, the dest would contain float numbers with long tails.
            # E.g. input == 0.806; dest == 0.8059999999999999
            # sanger_dfm = pandas.read_csv(dest, sep='\t', header=None, names=col_names,
            #                              usecols=['name', 'sum'], dtype=str)
            sanger_dfm = pd.read_csv(dest.name, sep='\t', header=None, names=col_names, usecols=['name', 'sum'])

            # Choose column `sum` to represent tfbs summary scores
            sanger_dfm.rename(columns={'sum': 'SangerTfbsSummary'}, inplace=True)

            return sanger_dfm

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', index=False, header=True)
