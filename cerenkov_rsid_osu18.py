import sys
import pandas as pd
import argparse
from coord_util import CoordUtil
from allele_util import AlleleUtil
from hgmd_client import HgmdClient


if sys.version_info[0] == 3:
    # Using Python3
    from snap_client import SnapQuery
else:
    # Using Python2
    from snap_client_py2 import SnapQuery


def get_locus_map(r_rsid):
    query_1kg = SnapQuery().dataset('onekgpilot'). \
        population(['CEU', 'YRI', 'CHBJPT']). \
        r_squared_threshold(0.8). \
        distance_limit_in_kb(50). \
        snp(r_rsid)

    locus_1kg = query_1kg.execute(verbose=True)

    query_rel21 = SnapQuery().dataset('rel21'). \
        population(['CEU', 'YRI', 'JPT+CHB']). \
        r_squared_threshold(0.8). \
        distance_limit_in_kb(50). \
        snp(r_rsid)

    locus_rel21 = query_rel21.execute(verbose=True)
    # Use entries only in Chromosome X of HapMap release 21
    locus_rel21 = locus_rel21.loc[locus_rel21['Chromosome'] == 'chrX']

    locus_map = pd.concat([locus_1kg, locus_rel21], ignore_index=True, axis=0).\
        drop_duplicates(['SNP', 'Proxy'], keep='last')
    locus_map = locus_map[~locus_map['Proxy'].isin(r_rsid)]  # remove rSNPs in cSNPs

    return locus_map


def apply_filters(rsid):
    coord_util = CoordUtil()
    coord_util.db_config_key = 'local_hg19'
    # coord_util.temp_dest = 'Coord.tsv'
    coord = coord_util.extract(_input=rsid)
    coord = CoordUtil.pce_filter(coord, verbose=True)

    allele_util = AlleleUtil()
    allele_util.db_config_key = 'local_hg19'
    # allele_util.temp_dest = 'Allele.tsv'
    allele = allele_util.extract(_input=rsid)
    allele = AlleleUtil.maf_filter(allele, maf_threshold=0.05, use_biomart=True, verbose=True)

    df = allele.merge(coord, how='inner', on=['name', 'chrom'], copy=True)

    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Fetch rsid from HGMD/SNAP.", allow_abbrev=False)

    parser.add_argument('-r', '--rsnp', dest='r_dest', type=str, required=True,
                        help="output file of rSNP rsid")
    parser.add_argument('-c', '--csnp', dest='c_dest', type=str, required=True,
                        help="output file of cSNP rsid")

    args = parser.parse_args()

    print("[cerenkov_rsid_osu18] rSNP rsid output: {}; cSNP rsid output: {}.".
          format(args.r_dest, args.c_dest))

    with HgmdClient() as hgmd_client:
        r_rsid = hgmd_client.select_rsid()

    r_df = apply_filters(r_rsid)

    _locus_map = get_locus_map(r_df['name'].tolist())

    c_rsid = _locus_map['Proxy'].unique().tolist()

    c_df = apply_filters(c_rsid)

    r_df.to_csv(args.r_dest, header=True, index=False, sep='\t')
    c_df.to_csv(args.c_dest, header=True, index=False, sep='\t')

    print("[cerenkov_rsid_osu18] Done!")

