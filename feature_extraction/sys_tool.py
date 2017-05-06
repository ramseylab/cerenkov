import os
import subprocess
from io import StringIO

# E.g. /home/ramseylab/Git-repo/recomb-rsnp-data
__BASE_DIR = os.path.dirname(os.path.realpath(__file__))
# E.g. /home/ramseylab/Git-repo/recomb-rsnp-data/R
__R_DIR = __BASE_DIR + "/R"
# E.g. /home/ramseylab/Git-repo/recomb-rsnp-data/source_data
__SRC_DATA_DIR = __BASE_DIR + "/source_data"
# E.g. /home/ramseylab/Git-repo/recomb-rsnp-data/source_data/UCSC_DHS
__DHS_DIR = __SRC_DATA_DIR + "/UCSC_DHS"
# E.g. /home/ramseylab/Git-repo/recomb-rsnp-data/source_data/GERP
__GERP_DIR = __SRC_DATA_DIR + "/GERP"
# E.g. /home/ramseylab/Git-repo/recomb-rsnp-data/source_data/GTEx_Analysis_V4_eQTLs
__EQTL_DIR = __SRC_DATA_DIR + "/GTEx_Analysis_V4_eQTLs"
# E.g. /home/ramseylab/Git-repo/recomb-rsnp-data/source_data/FSU_Repli_chip
__FSU_REPLI_CHIP_DIR = __SRC_DATA_DIR + "/FSU_Repli_chip"
# E.g. /home/ramseylab/Git-repo/recomb-rsnp-data/source_data/UW_Repli_Seq
__UW_REPLI_CHIP_DIR = __SRC_DATA_DIR + "/UW_Repli_Seq"
# E.g. /home/ramseylab/Git-repo/recomb-rsnp-data/source_data/Sanger_TFBS_Summary
__SANGER_TFBS_SUMMARY_DIR = __SRC_DATA_DIR + "/Sanger_TFBS_Summary"
# E.g. /home/ramseylab/Git-repo/recomb-rsnp-data/source_data/GIS_ChIA_PET
__CHIA_PRT_DIR = __SRC_DATA_DIR + "/GIS_ChIA_PET"
# E.g. /home/ramseylab/Git-repo/recomb-rsnp-data/source_data/CADD
__CADD_DIR = __SRC_DATA_DIR + "/CADD"
# E.g. /home/ramseylab/Git-repo/recomb-rsnp-data/source_data/fitCons
__FITCONS_DIR = __SRC_DATA_DIR + "/fitCons"
# E.g. /home/ramseylab/Git-repo/recomb-rsnp-data/source_data/Columbia_Eigen_Score
__EIGEN_DIR = __SRC_DATA_DIR + "/Columbia_Eigen_Score"


__PATHS = {
    'base': __BASE_DIR,
    'r': __R_DIR,
    'src_data': __SRC_DATA_DIR,
    'dhs': __DHS_DIR,
    'gerp': __GERP_DIR,
    'eqtl': __EQTL_DIR,
    'fsu_repli_chip': __FSU_REPLI_CHIP_DIR,
    'uw_repli_chip': __UW_REPLI_CHIP_DIR,
    'sanger_tfbs': __SANGER_TFBS_SUMMARY_DIR,
    'chia_pet': __CHIA_PRT_DIR,
    'cadd': __CADD_DIR,
    'fitcons': __FITCONS_DIR,
    'eigen': __EIGEN_DIR
}


def find_directory(name):
    return __PATHS[name.lower()]


def run_r_script(filename, args):
    script_path = find_directory("r") + "/" + filename

    if not os.path.isfile(script_path):
        raise FileNotFoundError(script_path + " not found.")

    cmd = StringIO()
    cmd.write("Rscript --vanilla " + script_path)

    arg_seq = ' '.join(args)
    cmd.write(' ' + arg_seq)

    print("SysTool: running command \r\n\t %s" % cmd.getvalue())

    subprocess.call(cmd.getvalue(), shell=True)

    cmd.close()


def run_bigWigAverageOverBed(args):
    """
    `bigWigAverageOverBed` is one of the `kentUtils` commands that are used to calculate GERP scores.

    `kentUtils` is hosted on
        https://github.com/ENCODE-DCC/kentUtils,
    or downloadable from
        http://hgdownload.cse.ucsc.edu/admin/exe/

    usage: bigWigAverageOverBed in.bw in.bed out.tab

    The output columns are:
       name - name field from bed, which should be unique
       size - size of bed (sum of exon sizes
       covered - # bases within exons covered by bigWig
       sum - sum of values over all bases covered
       mean0 - average over bases with non-covered bases counting as zeroes
       mean - average over just covered bases
    Options:
       -stats=stats.ra - Output a collection of overall statistics to stat.ra file
       -bedOut=out.bed - Make output bed that is echo of input bed but with mean column appended
       -sampleAroundCenter=N - Take sample at region N bases wide centered around bed item, rather
                         than the usual sample in the bed item.
       -minMax - include two additional columns containing the min and max observed in the area.

    :param args: a list of [in.bw, snp.bed, gerp.txt] in which
        in.bw: path to `All_hg19_RS.bw`. Downloadable (>7G) from
                http://hgdownload.cse.ucsc.edu/gbdb/hg19/bbi/All_hg19_RS.bw
    :return: None
    """

    cmd = StringIO()
    cmd.write("bigWigAverageOverBed")

    arg_seq = ' '.join(args)
    cmd.write(' ' + arg_seq)

    print("SysTool: running command \r\n\t %s" % cmd.getvalue())

    subprocess.call(cmd.getvalue(), shell=True)

    cmd.close()
