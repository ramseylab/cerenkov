import tempfile
import pandas as pd
from sys_tool import run_r_script
from abstract_feature_util import AbstractFeatureUtil


def cnt_ovlp_occur(string, sub):
    """
    Count overlapping occurrence of `sub` in `string`. E.g. `cnt_ovlp_occur("CCC", "CC")` would return 2.

    Python's `str.count(word)` function does NOT count overlapping occurrence. E.g. "CCC".count("CC") only returns 1.

    :param string:
    :param sub:
    :return:
    """

    count = start = 0
    while True:
        # `str.find()` returns index if found and -1 otherwise.
        start = string.find(sub, start) + 1
        if start > 0:
            count += 1
        else:
            return count


class DnaShapeGcContentUtil(AbstractFeatureUtil):
    def __init__(self, reproduce_osu17=False):
        super(DnaShapeGcContentUtil, self).__init__()
        self.r_fn = "getFlankingSequenceFeatures.R"
        self.reproduce_osu17 = reproduce_osu17

    @staticmethod
    def __transform_col_osu17(dfm, col_name='dnaSeq', drop_col=True):
        """
        Transform column `dnaSeq` into 3 features, `localGC`, `local_purine` and `local_CpG`.
        Drop column `dnaSeq` if needed

        :param dfm:
        :param col_name:
        :param drop_col:
        :return:
        """

        def __transform_seq(dna_seq):
            dna_seq = dna_seq.upper()  # use IUPAC alphabets

            n_a = dna_seq.count('A')
            n_c = dna_seq.count('C')
            n_g = dna_seq.count('G')

            n_gc = n_g + n_c
            n_purine = n_a + n_g
            """
            This is not the right way to count CpG sites.
            Just for reproducing osu17 feature matrix.
            """
            n_cpg = dna_seq.count('CG') + dna_seq.count('GC') + \
                cnt_ovlp_occur(dna_seq, 'CC') + cnt_ovlp_occur(dna_seq, 'GG')

            return pd.Series(data={'local_GC': n_gc,
                                   'local_purine': n_purine,
                                   'local_CpG': n_cpg})

        result = dfm.loc[:, col_name].apply(__transform_seq)

        dfm = pd.concat([dfm, result], axis=1)

        if drop_col:
            dfm = dfm.drop(col_name, axis=1)

        return dfm

    @staticmethod
    def transform_col(dfm, col_name='dnaSeq', drop_col=True):
        """
        Transform column `dnaSeq` into 3 features, `localGC`, `local_purine` and `local_CpG`.
        Drop column `dnaSeq` if needed

        :param dfm:
        :param col_name:
        :param drop_col:
        :return:
        """
        def __transform_seq(dna_seq):
            dna_seq = dna_seq.upper()  # use IUPAC alphabets

            n_a = dna_seq.count('A')
            n_c = dna_seq.count('C')
            n_g = dna_seq.count('G')

            n_gc = n_g + n_c
            n_purine = n_a + n_g
            n_cpg = dna_seq.count('CG') + dna_seq.count('GC')  # consider both strands

            return pd.Series(data={'local_GC': n_gc,
                                   'local_purine': n_purine,
                                   'local_CpG': n_cpg})

        result = dfm.loc[:, col_name].apply(__transform_seq)

        dfm = pd.concat([dfm, result], axis=1)

        if drop_col:
            dfm = dfm.drop(col_name, axis=1)

        return dfm

    def get_feat(self, _input):
        """

        :param _input: the path of SNP BED file
        :return:
        """
        snp_bed_path = _input

        if self.temp_dest is None:
            with tempfile.NamedTemporaryFile(prefix='cerenkov-') as dest:
                run_r_script(self.r_fn, [snp_bed_path, dest.name])

                result_dfm = pd.read_table(dest.name, header=0)
        else:
            run_r_script(self.r_fn, [snp_bed_path, self.temp_dest])

            result_dfm = pd.read_table(self.temp_dest, header=0)

        if self.reproduce_osu17:
            result_dfm = DnaShapeGcContentUtil.__transform_col_osu17(result_dfm, drop_col=False)
        else:
            result_dfm = DnaShapeGcContentUtil.transform_col(result_dfm)

        return result_dfm

    def save_temp(self, _result):
        pass

if __name__ == '__main__':

    dsgc_util = DnaShapeGcContentUtil(reproduce_osu17=True)
    dsgc_df = dsgc_util.extract("TEST.bed")

    print(dsgc_df)
