import os
import subprocess
from io import StringIO

# E.g. /home/ramseylab/Git-repo/recomb-rsnp-data
__BASE_DIR = os.path.dirname(os.path.realpath(__file__))

__R_DIR = os.path.join(__BASE_DIR, "R")

__SRC_DATA_DIR = os.path.join(__BASE_DIR, "source_data")

__DHS_DIR = os.path.join(__SRC_DATA_DIR, "UCSC_DHS")
__GERP_DIR = os.path.join(__SRC_DATA_DIR, "GERP")
__EQTL_DIR = os.path.join(__SRC_DATA_DIR, "GTEx_Analysis_V4_eQTLs")
__FSU_REPLI_CHIP_DIR = os.path.join(__SRC_DATA_DIR, "FSU_Repli_chip")
__UW_REPLI_CHIP_DIR = os.path.join(__SRC_DATA_DIR, "UW_Repli_Seq")
__SANGER_TFBS_SUMMARY_DIR = os.path.join(__SRC_DATA_DIR, "Sanger_TFBS_Summary")
__CHIA_PRT_DIR = os.path.join(__SRC_DATA_DIR, "GIS_ChIA_PET")
__CADD_DIR = os.path.join(__SRC_DATA_DIR, "CADD")
__FITCONS_DIR = os.path.join(__SRC_DATA_DIR, "fitCons")
__EIGEN_DIR = os.path.join(__SRC_DATA_DIR, "Columbia_Eigen_Score")
__GWAVA_DIR = os.path.join(__SRC_DATA_DIR, "GWAVA")
__JASPAR_TFBS_DIR = os.path.join(__SRC_DATA_DIR, "Jaspar_TFBS")
__AUGMENT_DIR = os.path.join(__SRC_DATA_DIR, "augment_osu_features_datafiles")


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
    'eigen': __EIGEN_DIR,
    'gwava': __GWAVA_DIR,
    'jaspar_tfbs': __JASPAR_TFBS_DIR,
    'augment': __AUGMENT_DIR
}


def find_directory(name):
    return __PATHS[name.lower()]


def run_r_script(filename, args):
    script_path = os.path.join(find_directory("r"), filename)

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
    cmd.write("/Users/kavisen/documents/osu/2021/winter_term/gra_steve/cerenkov/genome-wide/feature_extraction/source_data/bigWigAverageOverBed")

    arg_seq = ' '.join(args)
    cmd.write(' ' + arg_seq)

    print("SysTool: running command \r\n\t %s" % cmd.getvalue())

    subprocess.call(cmd.getvalue(), shell=True)

    cmd.close()
