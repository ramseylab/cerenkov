import tempfile
import pandas as pd
from sys_tool import run_r_script
from abstract_feature_util import AbstractFeatureUtil


class DnaShapeGcContentUtil(AbstractFeatureUtil):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.r_fn = "getFlankingSequenceFeatures.R"

    @staticmethod
    def _read_dest(path):
        result_dfm = pd.read_table(path, header=0)
        return result_dfm

    def get_feat(self, _input):
        """

        :param _input: the path of SNP BED file
        :return:
        """
        snp_bed_path = _input

        if self.temp_dest is None:
            with tempfile.NamedTemporaryFile(prefix='cerenkov-') as dest:
                run_r_script(self.r_fn, [snp_bed_path, dest.name])
                return self._read_dest(dest.name)
        else:
            run_r_script(self.r_fn, [snp_bed_path, self.temp_dest])
            return self._read_dest(self.temp_dest)

    def save_temp(self, _result):
        pass

if __name__ == '__main__':
    from sys_tool import find_directory

    rsnp_bed_path = "{}/RSNP_50kb.bed".format(find_directory('CADD'))

    dsgc_util = DnaShapeGcContentUtil()

    result = dsgc_util.extract(_input=rsnp_bed_path)

    print(result)
