import os
import pandas as pd
import sys_tool
import tempfile
from abstract_feature_util import AbstractFeatureUtil


class GerpUtil(AbstractFeatureUtil):
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

            gerp_dfm = pd.read_csv(dest.name, sep='\t', header=None, names=col_names, usecols=['name', 'sum'])
            gerp_dfm.rename(columns={'sum': 'gerp'}, inplace=True)
            return gerp_dfm

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', index=False, header=True)


class AvgGerpUtil(GerpUtil):
    def __init__(self, delta_start, delta_end):
        """
        Specify the coordinate interval on which the average GERP will be calculated.
        Technically the interval is [chromStart + delta_start, chromEnd + delta_end)

        :param delta_start: a number (<=0) added to 'chromStart' to get the left endpoint of the coordinate interval.
        :param delta_end: a number (>=0) added to 'chromEnd' to get the right endpoint of the coordinate interval
        """
        super(AvgGerpUtil, self).__init__()
        self.delta_start = delta_start
        self.delta_end = delta_end

    def get_feat(self, _input):
        """

        :param _input: SNP data frame in BED format
        :return:
        """
        snp_dfm = _input.loc[:, ['chrom', 'chromStart', 'chromEnd', 'name']]
        snp_dfm.loc[:, 'chromStart'] = snp_dfm.loc[:, 'chromStart'] + self.delta_start
        snp_dfm.loc[:, 'chromEnd'] = snp_dfm.loc[:, 'chromEnd'] + self.delta_end

        # Write to a temp BED file; feed this temp BED path to `GerpUtil.get_feat()`
        with tempfile.NamedTemporaryFile(prefix='cerenkov-') as snp_temp_bed:
            snp_dfm.to_csv(snp_temp_bed.name, sep='\t', index=False, header=False)

            src_bw = os.path.join(self.src_data_dir, self.src_data_fn)

            with tempfile.NamedTemporaryFile(prefix='cerenkov-') as dest:
                # This command would output a file with undesired columns. Need pruning.
                sys_tool.run_bigWigAverageOverBed([src_bw, snp_temp_bed.name, dest.name])

                """
                The dest has no header. Its 6 columns are:
                   name - name field from bed, which should be unique
                   size - size of bed (sum of exon sizes
                   covered - # bases within exons covered by bigWig
                   sum - sum of values over all bases covered
                   mean0 - average over bases with non-covered bases counting as zeroes
                   mean - average over just covered bases
                """

                # Choose column `mean0` to represent avg gerp scores
                col_names = ['name', 'size', 'covered', 'sum', 'mean0', 'mean']

                gerp_dfm = pd.read_csv(dest.name, sep='\t', header=None, names=col_names, usecols=['name', 'mean0'])
                gerp_dfm.rename(columns={'mean0': 'avg_gerp'}, inplace=True)
                return gerp_dfm
