import tempfile
import pandas as pd
from sys_tool import run_r_script
from abstract_feature_util import AbstractFeatureUtil


class AugmentedUtil(AbstractFeatureUtil):
    def __init__(self):
        super(AugmentedUtil, self).__init__()
        self.r_fn = "augment_osu_features_v2.R"

    def get_feat(self, _input):
        """
        :param _input: the path of feature matrix
        :return:
        """
        feat_matrix_path = _input

        if self.temp_dest is None:
            with tempfile.NamedTemporaryFile(prefix='cerenkov-') as dest:
                run_r_script(self.r_fn, [self.src_data_dir, feat_matrix_path, dest.name])

                result_dfm = pd.read_table(dest.name, header=0)
        else:
            run_r_script(self.r_fn, [self.src_data_dir, feat_matrix_path, self.temp_dest])

            result_dfm = pd.read_table(self.temp_dest, header=0)

        return result_dfm

    def save_temp(self, _result):
        pass
