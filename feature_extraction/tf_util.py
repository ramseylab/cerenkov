import os
import pandas as pd
from pybedtools import BedTool
from io import StringIO
from abstract_feature_util import AbstractFeatureUtil
from genome_browser_client import GenomeBrowserClient
from sys_tool import find_directory


class TfUtil(AbstractFeatureUtil):
    def __init__(self, reproduce_osu17=False):
        super(TfUtil, self).__init__()
        self.reproduce_osu17 = reproduce_osu17

    @staticmethod
    def sync_tf_columns(df_a, df_b, tf_col_prefix='tf_'):
        """
        Two data frames after one-hot-encoding may contains different TF columns.
        This function extends these two data frames with the union of TF columns.

        :param df_a:
        :param df_b:
        :param tf_col_prefix:
        :return:
        """
        a_col_names = set([col for col in df_a.columns.values.tolist() if col.startswith(tf_col_prefix)])
        b_col_names = set([col for col in df_b.columns.values.tolist() if col.startswith(tf_col_prefix)])

        # Extend columns of `df_a`
        a_extra_col_names = (b_col_names - a_col_names)
        if a_extra_col_names:
            a_extra_col = dict.fromkeys(list(a_extra_col_names), 0)
            df_a = df_a.assign(**a_extra_col)

        # Extend columns of `df_b`
        b_extra_col_names = (a_col_names - b_col_names)
        if b_extra_col_names:
            b_extra_col = dict.fromkeys(list(b_extra_col_names), 0)
            df_b = df_b.assign(**b_extra_col)

        return df_a, df_b

    @staticmethod
    def binary_encode(snp_dfm, cat_col_name="tfName", cat_col_sep=',', bin_col_prefix=None):
        """
        Binary-encode categorical column `cat_col_name` separated by `cat_col_sep`
        in data frame `snp_dfm` to multiple binary columns with the same prefix `bin_col_prefix`

        MySQL returns `GROUP_CONCAT(tf.name)` in `tfName` column in a comma-separated string.
            E.g. `ARID3A,ATF1,ATF2` stands for 3 TFs
        This function separates this string by commas and 3 new columns,
            `tf_ARID3A`, `tf_ATF1` and `tf_ATF2` would be 1 for this SNP

        :param snp_dfm: the data frame
        :param cat_col_name: the name of the categorical column whose values would be encoded
        :param cat_col_sep: the separator of the categorical values
        :param bin_col_prefix: the prefix of the binary columns after one-hot encoding
        :return:
        """
        ohe = snp_dfm.loc[:, cat_col_name].str.get_dummies(sep=cat_col_sep)

        if bin_col_prefix is not None:
            # Add a prefix to all column names
            ohe = ohe.add_prefix(bin_col_prefix)

        snp_dfm = pd.concat([snp_dfm, ohe], axis=1).drop(cat_col_name, axis=1)

        return snp_dfm

    def get_feat(self, _input):
        """
        :param _input: the SNP data frame
        :return:
        """
        snp_dfm = _input.loc[:, ['chrom', 'chromStart', 'chromEnd', 'name']]

        with GenomeBrowserClient(self.db_config_key) as gb_client:
            tf_df = gb_client.fetch_tf(snp_dfm.loc[:, 'name'].tolist())

        snp_dfm = snp_dfm.merge(tf_df, how='left', on=['name', 'chrom'], copy=True)

        snp_dfm.loc[:, 'tfCount'] = snp_dfm.loc[:, 'tfCount'].fillna(0)

        """
        No prefix for TF column names in osu17; Would prefer a "tf_" prefix for osu18
        An all-zero TF `POLR3R`, which `gb_client` won't return in its result, is presented in osu17 feature matrix.
        """
        if self.reproduce_osu17:
            snp_dfm = self.binary_encode(snp_dfm, cat_col_name="tfName", cat_col_sep=',', bin_col_prefix=None)
            snp_dfm = snp_dfm.assign(POLR3G=0)
        else:
            snp_dfm = self.binary_encode(snp_dfm, cat_col_name="tfName", cat_col_sep=',', bin_col_prefix="tf_")

        snp_dfm = snp_dfm.drop(['chrom', 'chromStart', 'chromEnd'], axis=1)

        return snp_dfm

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', header=True, index=False)


# TODO Binary-encode jaspar TFBS. Now a single SNP may be mapped to multiple TFBS
class JasparTfbsUtil(AbstractFeatureUtil):
    """
    annotate SNPs based on TFBS motifs using JASPAR motif match data from Ensembl BioMart
    (Ensembl Release 75, archive.ensembl.org)
    """
    def __init__(self):
        raise NotImplementedError('Binary-encoding not implemented.')

    def __read_src_file(self):
        src_file = os.path.join(self.src_data_dir, self.src_data_fn)
        jaspar_df = pd.read_table(src_file, header=0,
                                  names=['id', 'chrom', 'chromStart', 'chromEnd', 'name'],
                                  usecols=['chrom', 'chromStart', 'chromEnd', 'name'])
        jaspar_df.loc[:, 'chrom'] = 'chr' + jaspar_df.loc[:, 'chrom']

        jaspar_bed_obj = BedTool(jaspar_df.to_string(index=False, header=False, index_names=False), from_string=True)

        return jaspar_bed_obj

    def get_feat(self, _input):
        """

        :param _input: SNP data frame in BED format
        :return:
        """
        snp_dfm = _input.loc[:, ['chrom', 'chromStart', 'chromEnd', 'name']]
        snp_bed_obj = BedTool(snp_dfm.to_string(index=False, header=False, index_names=False), from_string=True)

        jasper_bed_obj = self.__read_src_file()

        intx = snp_bed_obj.intersect(jasper_bed_obj, wb=True)

        intx_dfm = pd.read_table(StringIO(str(intx)), header=None,
                                 names=['snpChrom', 'snpChromStart', 'snpChromEnd', 'snpName',
                                        'tfbsChrom', 'tfbsChromStart', 'tfbsChromEnd', 'tfbsName'])
        # TODO is 'pwm' a good feature name here?
        result = intx_dfm[['snpName', 'snpChrom', 'tfbsName']].rename(columns={'snpName': 'name',
                                                                               'snpChrom': 'chrom',
                                                                               'tfbsName': 'pwm'})

        result = snp_dfm.merge(result, how='left', on=['name', 'chrom'], copy=True)

        return result.loc[:, ['name', 'pwm']]

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', index=False, header=True)


if __name__ == '__main__':
    rsnp_bed_fn = os.path.join(find_directory('dhs'), "RSNP_50kb.bed")
    rsnp_dfm = pd.read_table(rsnp_bed_fn, header=None,
                             names=['chrom', 'chromStart', 'chromEnd', 'name'])

    jaspar_tf_util = JasparTfbsUtil()
    jaspar_tf_util.src_data_dir = find_directory('Jaspar_TFBS')
    jaspar_tf_util.src_data_fn = 'jaspar_tfbs_ensembl_75_hg19.txt'

    result = jaspar_tf_util.extract(rsnp_dfm)

    result.to_csv("FOO_RSNP.tsv", sep='\t', header=True, index=False)


# if __name__ == '__main__':
#     rsnp_bed_fn = "{}/RSNP_50kb.bed".format(find_directory('dhs'))
#     rsnp_dfm = pd.read_table(rsnp_bed_fn, header=None,
#                                  names=['chrom', 'chromStart', 'chromEnd', 'name'])
#     csnp_bed_fn = "{}/CSNP_50kb.bed".format(find_directory('dhs'))
#     csnp_dfm = pd.read_table(csnp_bed_fn, header=None,
#                                  names=['chrom', 'chromStart', 'chromEnd', 'name'])
#
#     tf_util = TFUtil()
#     tf_util.db_config_key = 'local_hg19'
#
#     # tf_util.temp_dest = 'FOO_CSNP.tsv'
#     csnp_tf = tf_util.extract(csnp_dfm)
#
#     # tf_util.temp_dest = 'FOO_RSNP.tsv'
#     rsnp_tf = tf_util.extract(rsnp_dfm)
#
#     print(rsnp_tf)
#     print(csnp_tf)
#
#     rsnp_tf, csnp_tf = TFUtil.sync_tf_columns(rsnp_tf, csnp_tf, 'tf_')
#
#     rsnp_tf.to_csv("FOO_RSNP.tsv", sep='\t', header=True, index=False)
#     csnp_tf.to_csv("FOO_CSNP.tsv", sep='\t', header=True, index=False)
