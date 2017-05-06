import sys
import pandas as pd
from hgmd_client import HgmdClient
from coord_util import CoordUtil
from allele_util import AlleleUtil

if sys.version_info[0] == 3:
    # Using Python3
    from snap_client import SnapQuery
else:
    # Using Python2
    from snap_client_py2 import SnapQuery


def __get_locus_map(r_rsid):
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


def __filter(rsid):
    coord_util = CoordUtil()
    coord_util.db_config_key = 'local_hg19'
    # coord_util.temp_dest = 'Coord.tsv'

    allele_util = AlleleUtil()
    allele_util.db_config_key = 'local_hg19'
    # allele_util.temp_dest = 'Allele.tsv'

    coord = coord_util.extract(_input=rsid)
    allele = allele_util.extract(_input=rsid)

    df = allele.merge(coord, how='inner', on=['name', 'chrom'], copy=True)

    df = CoordUtil.pce_filter(df, verbose=True)

    return df

if __name__ == '__main__':

    with HgmdClient() as hgmd_client:
        r_rsid = hgmd_client.select_rsid()

    r_df = __filter(r_rsid)

    locus_map = __get_locus_map(r_df['name'].tolist())

    c_rsid = locus_map['Proxy'].unique().tolist()

    c_df = __filter(c_rsid)

    r_df.to_csv('OSU18_rsnp.tsv', header=True, index=False, sep='\t')
    c_df.to_csv('OSU18_csnp.tsv', header=True, index=False, sep='\t')

    print("[cerenkov_rsid_osu18] Done!")