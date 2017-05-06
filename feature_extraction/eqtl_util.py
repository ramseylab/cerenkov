from abstract_feature_util import AbstractFeatureUtil
from sys_tool import find_directory
import pandas as pd
import numpy as np
from math import log10
from glob import glob


class EqtlUtil(AbstractFeatureUtil):
    @staticmethod
    def __transform_p_value(p_value):
        if np.isnan(p_value):
            return 0
        else:
            return -log10(p_value)

    @staticmethod
    def __yield_eqtl_dfm(snp_dfm, eqtl_paths):
        for eqtl_path in eqtl_paths:
            eqtl_dfm = pd.read_table(eqtl_path, header=0)
            eqtl_pvalue_dfm = snp_dfm.merge(eqtl_dfm, how='left', left_on='name',
                                            right_on='SNP', copy=True).loc[:, ['name', 'P_Val']]
            yield eqtl_pvalue_dfm

    def get_feat(self, _input):
        snp_dfm = _input

        if self.src_data_fn is None:
            # list all files with '.portal.eqtl' extension in the specified directory
            eqtl_paths = glob(pathname="{dir}/{fn}".format(dir=self.src_data_dir, fn="*.portal.eqtl"))
        else:
            eqtl_paths = ["{dir}/{fn}".format(dir=self.src_data_dir, fn=fn) for fn in self.src_data_fn]

        yielded_eqtl_dfms = list(self.__yield_eqtl_dfm(snp_dfm, eqtl_paths))
        result = pd.concat(yielded_eqtl_dfms, axis=0)

        # drop rows with 'P_Val' being NaN
        result = result.dropna(axis=0)

        # find minimum non-NAN 'P_Val' for each SNP
        result = result.groupby('name').agg(min).reset_index()

        # some SNPs has no non-NAN 'P_Val'; mark their aggregated 'P_VAL' NAN
        result = snp_dfm.merge(result, how='left', on='name', copy=True)

        result.loc[:, 'P_Val'] = result.loc[:, 'P_Val'].apply(self.__transform_p_value)
        result = result.rename(columns={'P_Val': 'eqtlPvalue'})
        return result.loc[:, ['name', 'eqtlPvalue']]

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', header=True, index=False)


if __name__ == '__main__':
    rsnp_dfm = pd.read_table("{}/RSNP_50kb.bed".format(find_directory('CADD')), header=None,
                             names=['chrom', 'chromStart', 'chromEnd', 'name'])
    print(rsnp_dfm.shape)

    eqtl_util = EqtlUtil()
    eqtl_util.src_data_dir = find_directory('eqtl')
    # eqtl_util.src_data_fn = ['Stomach.portal.eqtl', 'Heart_Left_Ventricle.portal.eqtl']
    eqtl_util.temp_dest = 'foo_rsnp.txt'

    result = eqtl_util.extract(rsnp_dfm)
    print(result)
