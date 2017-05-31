import os
import pandas as pd
from abstract_feature_util import AbstractFeatureUtil


class EigenUtil(AbstractFeatureUtil):
    def __init__(self):
        super(EigenUtil, self).__init__()
        self._eigen_dfm = None

    @property
    def eigen_dfm(self):
        # Lazy Initialization
        if self._eigen_dfm is None:
            eigen_fn = os.path.join(self.src_data_dir, self.src_data_fn)
            eigen_dfm = pd.read_table(eigen_fn, header=None, usecols=[0, 1, 29])  # 29th column is "raw eigen"
            eigen_dfm.columns = ['chrom', 'chromEnd', 'eigen']
            eigen_dfm.drop_duplicates(inplace=True)
            eigen_dfm.loc[:, 'chrom'] = eigen_dfm.loc[:, 'chrom'].apply(lambda x: "chr{}".format(x))

            self._eigen_dfm = eigen_dfm

        return self._eigen_dfm

    def get_feat(self, _input):
        """

        :param _input: the SNP data frame in BED format
        :return:
        """
        snp_dfm = _input.loc[:, ['chrom', 'chromStart', 'chromEnd', 'name']]

        snp_eigen = snp_dfm.merge(self.eigen_dfm, how='left', on=['chrom', 'chromEnd'])

        snp_eigen.loc[:, 'eigen'] = snp_eigen.loc[:, 'eigen'].fillna(snp_eigen.loc[:, 'eigen'].mean(skipna=True))

        return snp_eigen.loc[:, ['name', 'eigen']]

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', index=False)
