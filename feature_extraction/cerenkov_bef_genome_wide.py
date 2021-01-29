import os
import sys
import argparse
import tempfile
import pandas as pd
import sys_tool
from allele_util import AlleleUtil
from gerp_util import AvgGerpUtil
from coord_util import CoordUtil
from tf_util import TfUtil
from dhs_util import MasterDhsUtil, UniformDhsUtil
from phastcons_util import PhastconsUtil
from tss_util import TssDistUtil
from eqtl_util import EqtlUtil
from gerp_util import GerpUtil
from genome_seg_util import GenomeSegUtil
from dna_shape_gc_content_util import DnaShapeGcContentUtil
from eigen_util import EigenUtil
from nki_lad_util import NkiLadUtil
from repli_chip_util import RepliChipUtil
from sanger_tfbs_util import SangerTfbsUtil
from vista_enhancer_util import VistaEnhancerUtil
from gwava_util import GwavaUtil
from augmented_feature_util import AugmentedUtil


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract BEF features.", allow_abbrev=False)

    parser.add_argument('-r', '--rsid', dest='r_src', type=str, required=True,
                        help="input file of rSNP rsid")
    parser.add_argument('-f', '--feature', dest='f_dest', type=str, required=True,
                        help="output file of the feature matrix")

    args = parser.parse_args()

    print("[cerenkov_bef_external] RSID src: {}; Feature Dest: {};".
          format(args.r_src, args.f_dest))

    # ----- 0. READ RSID -----
    rsid_fn = args.r_src
    rsid_dfm = pd.read_table(rsid_fn, header=None, names=['name'])
    rsid = rsid_dfm.loc[:, 'name'].tolist()

    # ----- 1. UTILITIES INIT -----

    db_config_key = 'local_hg19'  # 'remote_hg19'

    allele_util = AlleleUtil()
    allele_util.db_config_key = db_config_key

    coord_util = CoordUtil()
    coord_util.db_config_key = db_config_key

    tf_util = TfUtil(reproduce_osu17=True)
    tf_util.db_config_key = db_config_key

    # jaspar_tfbs_util = JasparTfbsUtil()
    # jaspar_tfbs_util.src_data_dir = sys_tool.find_directory('Jaspar_TFBS')
    # jaspar_tfbs_util.src_data_fn = 'jaspar_tfbs_ensembl_75_hg19.txt'

    mst_dhs_util = MasterDhsUtil()
    mst_dhs_util.db_config_key = db_config_key

    uni_dhs_util = UniformDhsUtil()
    uni_dhs_util.src_data_dir = sys_tool.find_directory("dhs")
    uni_dhs_util.src_data_fn = "UniformDnaseIHS"

    phastcons_util = PhastconsUtil()
    phastcons_util.db_config_key = db_config_key

    tss_dist_util = TssDistUtil()
    tss_dist_util.db_config_key = db_config_key

    eqtl_util = EqtlUtil()
    eqtl_util.src_data_dir = sys_tool.find_directory('eqtl')

    gerp_util = GerpUtil()
    gerp_util.src_data_dir = sys_tool.find_directory("gerp")
    gerp_util.src_data_fn = "All_hg19_RS.bw"

    avg_gerp_util = AvgGerpUtil(delta_start=-50, delta_end=50)
    avg_gerp_util.src_data_dir = sys_tool.find_directory("gerp")
    avg_gerp_util.src_data_fn = "All_hg19_RS.bw"
    # avg_gerp_file = os.path.join(sys_tool.find_directory("gerp"), "avg_gerp_osu17.tsv")

    g_seg_util = GenomeSegUtil(reproduce_osu17=True)
    g_seg_util.db_config_key = db_config_key

    dsgc_util = DnaShapeGcContentUtil(reproduce_osu17=True)

    # cadd_util = CaddUtil()
    # cadd_util.src_data_dir = sys_tool.find_directory('CADD')
    # cadd_util.src_data_fn = "1000G_phase3.tsv"

    eigen_util = EigenUtil(mode="genome-wide")
    eigen_util.src_data_dir = sys_tool.find_directory('eigen')
    eigen_util.src_data_fn = "mart_export_hg19_chr22_SNP.score"

    # fitcons_util = FitconsUtil()
    # fitcons_util.src_data_dir = sys_tool.find_directory("fitcons")
    # fitcons_util.src_data_fn = dict(
    #     fitConsGm="fc-gm-0.bed",
    #     fitConsH1="fc-h1-0.bed",
    #     fitConsHu="fc-hu-0.bed",
    #     fitConsI6="fc-i6-0.bed",
    # )

    nki_lad_util = NkiLadUtil()
    nki_lad_util.db_config_key = db_config_key

    fsu_rc_util = RepliChipUtil()
    fsu_rc_util.src_data_dir = sys_tool.find_directory("fsu_repli_chip")
    fsu_rc_util.src_data_fn = dict(
        FsuBg02esRep1="wgEncodeFsuRepliChipBg02esWaveSignalRep1.bed",
        FsuBg02esRep2="wgEncodeFsuRepliChipBg02esWaveSignalRep2.bed",
        FsuGm06990Rep1="wgEncodeFsuRepliChipGm06990WaveSignalRep1.bed",
        FsuGm06990Rep2="wgEncodeFsuRepliChipGm06990WaveSignalRep2.bed",
        FsuH1hescRep1="wgEncodeFsuRepliChipH1hescWaveSignalRep1.bed",
        FsuH1hescRep2="wgEncodeFsuRepliChipH1hescWaveSignalRep2.bed",
        FsuH1hescRep3="wgEncodeFsuRepliChipH1hescWaveSignalRep3.bed",
        FsuH7esRep1="wgEncodeFsuRepliChipH7esWaveSignalRep1.bed",
        FsuH7esRep2="wgEncodeFsuRepliChipH7esWaveSignalRep2.bed",
        FsuH9esRep1="wgEncodeFsuRepliChipH9esWaveSignalRep1.bed",
        FsuHelas3Rep1="wgEncodeFsuRepliChipHelas3WaveSignalRep1.bed",
        FsuImr90Rep1="wgEncodeFsuRepliChipImr90WaveSignalRep1.bed",
        FsuIpshfib2ips4Rep1="wgEncodeFsuRepliChipIpshfib2ips4WaveSignalRep1.bed",
        FsuIpshfib2ips4Rep2="wgEncodeFsuRepliChipIpshfib2ips4WaveSignalRep2.bed",
        FsuIpshfib2ips5Rep1="wgEncodeFsuRepliChipIpshfib2ips5WaveSignalRep1.bed",
        FsuIpshfib2ips5Rep2="wgEncodeFsuRepliChipIpshfib2ips5WaveSignalRep2.bed",
    )

    uw_rc_util = RepliChipUtil()
    uw_rc_util.src_data_dir = sys_tool.find_directory("uw_repli_chip")
    uw_rc_util.src_data_fn = dict(
        UwBg02esRep1="wgEncodeUwRepliSeqBg02esWaveSignalRep1.bed",
        UwBjRep1="wgEncodeUwRepliSeqBjWaveSignalRep1.bed",
        UwBjRep2="wgEncodeUwRepliSeqBjWaveSignalRep2.bed",
        UwGm06990Rep1="wgEncodeUwRepliSeqGm06990WaveSignalRep1.bed",
        UwGm12801Rep1="wgEncodeUwRepliSeqGm12801WaveSignalRep1.bed",
        UwGm12812Rep1="wgEncodeUwRepliSeqGm12812WaveSignalRep1.bed",
        UwGm12813Rep1="wgEncodeUwRepliSeqGm12813WaveSignalRep1.bed",
        UwGm12878Rep1="wgEncodeUwRepliSeqGm12878WaveSignalRep1.bed",
        UwHelas3Rep1="wgEncodeUwRepliSeqHelas3WaveSignalRep1.bed",
        UwHepg2Rep1="wgEncodeUwRepliSeqHepg2WaveSignalRep1.bed",
        UwHuvecRep1="wgEncodeUwRepliSeqHuvecWaveSignalRep1.bed",
        UwImr90Rep1="wgEncodeUwRepliSeqImr90WaveSignalRep1.bed",
        UwK562Rep1="wgEncodeUwRepliSeqK562WaveSignalRep1.bed",
        UwMcf7Rep1="wgEncodeUwRepliSeqMcf7WaveSignalRep1.bed",
        UwNhekRep1="wgEncodeUwRepliSeqNhekWaveSignalRep1.bed",
        UwSknshRep1="wgEncodeUwRepliSeqSknshWaveSignalRep1.bed",
    )

    sanger_tfbs_util = SangerTfbsUtil()
    sanger_tfbs_util.src_data_dir = sys_tool.find_directory("sanger_tfbs")
    sanger_tfbs_util.src_data_fn = "all_tfbs.bw"

    vh_util = VistaEnhancerUtil()
    vh_util.db_config_key = db_config_key

    gwava_util = GwavaUtil()
    gwava_util.src_data_dir = sys_tool.find_directory('GWAVA')
    gwava_util.src_data_fn = 'segmentation.bed'

    augmented_util = AugmentedUtil()
    augmented_util.src_data_dir = sys_tool.find_directory('augment')

    # ----- 2. ALLELE + COORD -----

    allele_df = allele_util.extract(_input=rsid)
    # allele_df = AlleleUtil.maf_filter(allele_df, maf_threshold=0.05, use_biomart=True, verbose=True)

    coord_df = coord_util.extract(_input=rsid)
    # coord_df = CoordUtil.pce_filter(coord_df, verbose=True)

    snp_df = allele_df.merge(coord_df, how='inner', on=['name', 'chrom'], copy=True)
    snp_df = AlleleUtil.identify_major_vs_minor(snp_df,
                                                maf_threshold=0.05,
                                                remove_acgt_cols=True,
                                                revise_when_freqs_equal=False,  # False for osu17; True for osu18
                                                verbose=True)
    snp_df = AlleleUtil.flip_alleles_on_rev_strand(snp_df, remove_strand_col=True, verbose=True)
    snp_df = CoordUtil.add_norm_coord(snp_df)

    # ----- 3. GENERATE BED & INPUT FILE FOR AUGMENTING FEATURES-----

    # Generate BED file from this data frame
    # Some bedtools operations require SORTED BED input
    snp_df = snp_df.sort_values(by=['chrom', 'chromStart'], ascending=True)

    bed_file = tempfile.NamedTemporaryFile(prefix='cerenkov-')
    aug_input = tempfile.NamedTemporaryFile(prefix='cerenkov-')

    snp_df.loc[:, ['chrom', 'chromStart', 'chromEnd', 'name']].to_csv(bed_file.name, sep='\t',
                                                                      header=False, index=False)
    snp_df.loc[:, ['name', 'chrom', 'chromStart', 'chromEnd', 'minorAlleleFreq']].to_csv(aug_input.name, sep='\t',
                                                                                         header=True, index=False)

    # ----- 4. THE REST FEATURES -----
    tf_df = tf_util.extract(snp_df)
    # jaspar_tfbs_df = jaspar_tfbs_util.extract(snp_df)
    mst_dhs_df = mst_dhs_util.extract(snp_df)
    uni_dhs_df = uni_dhs_util.extract(snp_df)
    phastcons_df = phastcons_util.extract(snp_df)
    tss_dist_df = tss_dist_util.extract(snp_df)
    eqtl_pval_df = eqtl_util.extract(snp_df)
    gerp_df = gerp_util.extract(bed_file.name)

    """
    1) Previous coordinate interval was set to [chromStart - 50, chromEnd -1) by mistake.
    2) `bigWigAverageOverBed`, where `avg_gerp` was calculated from, does NOT preserve the row order
        of the input BED in its output, and unfortunately we did not join the the output by rsID.
    Therefore it's very difficult to reproduce this error. We decided to simply keep previous `avg_gerp` values in
        a file for reproduction
    """
    avg_gerp_df = avg_gerp_util.extract(snp_df)

    g_seg_df = g_seg_util.extract(snp_df)
    dsgc_df = dsgc_util.extract(bed_file.name)
    # cadd_df = cadd_util.extract(snp_df)
    eigen_df = eigen_util.extract(snp_df)
    # fitcons_df = fitcons_util.extract(snp_df)
    nki_lad_df = nki_lad_util.extract(snp_df)
    fsu_rc_df = fsu_rc_util.extract(snp_df)
    uw_rc_df = uw_rc_util.extract(snp_df)
    sanger_tfbs_df = sanger_tfbs_util.extract(bed_file.name)
    vh_df = vh_util.extract(snp_df)
    gwava_df = gwava_util.extract(snp_df)
    aug_df = augmented_util.extract(aug_input.name)

    # ----- 5. MERGE -----

    df_lst = [snp_df, tf_df, mst_dhs_df, uni_dhs_df, phastcons_df, tss_dist_df, eqtl_pval_df, gerp_df, avg_gerp_df,
              g_seg_df, dsgc_df, eigen_df, nki_lad_df, fsu_rc_df, uw_rc_df, sanger_tfbs_df, vh_df, gwava_df, aug_df]

    feature_matrix = reduce(lambda left, right: pd.merge(left, right, on='name'), df_lst)

    # merge `label` column
    # feature_matrix = feature_matrix.merge(rsid_dfm, on='name')

    # drop `chromStart` and `chromEnd`
    feature_matrix = feature_matrix.drop(['chromStart', 'chromEnd'], axis=1)

    feature_matrix.to_csv(args.f_dest, sep='\t', header=True, index=False)

    # ----- 6. CLEAN UP -----

    bed_file.close()
    aug_input.close()

    # ----- 7. PROFILE -----

    util_lst = [allele_util, coord_util, tf_util, mst_dhs_util, uni_dhs_util, phastcons_util, tss_dist_util,
                eqtl_util, gerp_util, g_seg_util, dsgc_util, eigen_util, nki_lad_util, fsu_rc_util, uw_rc_util,
                sanger_tfbs_util, vh_util, gwava_util, augmented_util]

    total_seconds = sum([u.last_time_elapsed.seconds for u in util_lst])
    for util in util_lst:
        print("{:<21}: {} ({:.0f}%)".format(util.__class__.__name__,
                                            util.last_time_elapsed,
                                            100.0 * util.last_time_elapsed.seconds / total_seconds))
