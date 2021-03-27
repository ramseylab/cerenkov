import os
import pandas as pd
from abstract_feature_util import AbstractFeatureUtil


class EigenUtil(AbstractFeatureUtil):
    def __init__(self, mode):
        """
        mode == "osu17": 
            - only use "Eigen-raw", the 30th column of the Eigen Tabix results as features
            - column renamed to 'eigen'
        mode == "osu18": 
            - use both "Eigen-raw" and "Eigen-PC-raw", the 30th and 32nd columns of the Eigen Tabix results as features
            - columns renamed to 'eigen_raw' and 'eigen_pc_raw' respectively
        mode == "genome-wide":
            - same with "osu18", use both "Eigen-raw" and "Eigen-PC-raw" columns of the Eigen Tabix results as features
            - however the Eigen Tabix results are simplified to only 4 columns, "Eigen-raw" being the 3rd column and "Eigen-PC-raw" being the 4th
            - columns renamed to 'eigen_raw' and 'eigen_pc_raw' respectively
        """
        super(EigenUtil, self).__init__()
        self._eigen_dfm = None
        self.mode = mode

    @property
    def eigen_dfm(self):
        # Lazy Initialization
        if self._eigen_dfm is None:
            eigen_fn = os.path.join(self.src_data_dir, self.src_data_fn)
            
            if self.mode == "osu17":
                eigen_dfm = pd.read_table(eigen_fn, header=None, usecols=[0, 1, 29])  # 30th column is "Eigen-raw"
                eigen_dfm.columns = ['chrom', 'chromEnd', 'eigen']
            elif self.mode == "osu18":
                eigen_dfm = pd.read_table(eigen_fn, header=None, usecols=[0, 1, 29, 31])  # 32nd column is "Eigen-PC-raw"
                eigen_dfm.columns = ['chrom', 'chromEnd', 'eigen_raw', 'eigen_pc_raw']
            elif self.mode == "genome-wide":
                eigen_dfm = pd.read_table(eigen_fn, header=None, usecols=[0, 1, 2, 3])  # 3rd column is "Eigen-raw" abnd 4th is "Eigen-PC-raw"
                eigen_dfm.columns = ['chrom', 'chromEnd', 'eigen_raw', 'eigen_pc_raw']
            else:
                raise ValueError("Invalid mode. Got mode={}".format(self.mode))

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

        if self.mode == "osu17":
            snp_eigen.loc[:, 'eigen'] = snp_eigen.loc[:, 'eigen'].fillna(snp_eigen.loc[:, 'eigen'].mean(skipna=True))
            return snp_eigen.loc[:, ['name', 'eigen']]
        elif self.mode == "osu18" or self.mode == "genome-wide":
            snp_eigen.loc[:, 'eigen_raw'] = snp_eigen.loc[:, 'eigen_raw'].fillna(snp_eigen.loc[:, 'eigen_raw'].mean(skipna=True))
            snp_eigen.loc[:, 'eigen_pc_raw'] = snp_eigen.loc[:, 'eigen_pc_raw'].fillna(snp_eigen.loc[:, 'eigen_pc_raw'].mean(skipna=True))
            return snp_eigen.loc[:, ['name', 'eigen_raw', 'eigen_pc_raw']]
        else:
            raise ValueError("Invalid mode. Got mode={}".format(self.mode))

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', index=False)
