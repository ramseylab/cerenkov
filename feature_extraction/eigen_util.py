import pandas
from sys_tool import find_directory
from abstract_feature_util import AbstractFeatureUtil


class EigenUtil(AbstractFeatureUtil):
    def __init__(self):
        self._eigen_dfm = None

    @property
    def eigen_dfm(self):
        # Lazy Initialization
        if self._eigen_dfm is None:
            eigen_fn = "{dir}/{fn}".format(dir=self.src_data_dir, fn=self.src_data_fn)
            eigen_dfm = pandas.read_table(eigen_fn, header=None,
                                          usecols=[0, 1, 29])
            eigen_dfm.columns = ['chrom', 'chromEnd', 'eigenRaw']
            eigen_dfm.drop_duplicates(inplace=True)
            eigen_dfm.loc[:, 'chrom'] = eigen_dfm.loc[:, 'chrom'].apply(lambda x: "chr{}".format(x))
            eigen_dfm.set_index(['chrom', 'chromEnd'], inplace=True)

            self._eigen_dfm = eigen_dfm

        return self._eigen_dfm

    def get_feat(self, _input):
        """

        :param _input: the SNP data frame in BED format
        :return:
        """
        snp_dfm = _input
        snp_dfm.set_index(['chrom', 'chromEnd'], inplace=True)
        snp_eigen = snp_dfm.join(self.eigen_dfm, how='left').reset_index()[['name', 'eigenRaw']]

        return snp_eigen

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', index=False)


if __name__ == '__main__':
    eigen_util = EigenUtil()
    eigen_util.src_data_dir = find_directory('eigen')
    eigen_util.src_data_fn = "eigen_scores.tsv"

    rsnp_dfm = pandas.read_table("{}/RSNP_50kb.bed".format(find_directory('eigen')), header=None,
                                 names=['chrom', 'chromStart', 'chromEnd', 'name'])
    eigen_util.temp_dest = 'foo_rsnp.txt'
    eigen_util.extract(_input=rsnp_dfm)

    csnp_dfm = pandas.read_table("{}/CSNP_50kb.bed".format(find_directory('eigen')), header=None,
                                 names=['chrom', 'chromStart', 'chromEnd', 'name'])
    eigen_util.temp_dest = 'foo_csnp.txt'
    eigen_util.extract(_input=csnp_dfm)

    # rsnp_dfm = pandas.read_table("{}/RSNP_50kb.bed".format(find_directory('eigen')), header=None,
    #                               names=['chrom', 'chromStart', 'chromEnd', 'name'])
    # rsnp_dfm.set_index(['chrom', 'chromEnd'], inplace=True)
    #
    # csnp_dfm = pandas.read_table("{}/CSNP_50kb.bed".format(find_directory('eigen')), header=None,
    #                               names=['chrom', 'chromStart', 'chromEnd', 'name'])
    # csnp_dfm.set_index(['chrom', 'chromEnd'], inplace=True)
    #
    # eigen_dfm = pandas.read_table("{}/eigen_scores.tsv".format(find_directory('eigen')), header=None, usecols=[0,1,29])
    # eigen_dfm.columns = ['chrom', 'chromEnd', 'eigenRaw']
    # eigen_dfm.drop_duplicates(inplace=True)
    # eigen_dfm.loc[:, 'chrom'] = eigen_dfm.loc[:, 'chrom'].apply(lambda x: "chr{}".format(x))
    # eigen_dfm.set_index(['chrom', 'chromEnd'], inplace=True)
    #
    # rsnp_eigen = rsnp_dfm.join(eigen_dfm, how='left').reset_index()[['name', 'eigenRaw']]
    # csnp_eigen = csnp_dfm.join(eigen_dfm, how='left').reset_index()[['name', 'eigenRaw']]
    #
    # rsnp_eigen.to_csv("{}/RSNP_50kb_eigen.tsv".format(find_directory('eigen')), sep='\t', index=False)
    # csnp_eigen.to_csv("{}/CSNP_50kb_eigen.tsv".format(find_directory('eigen')), sep='\t', index=False)
