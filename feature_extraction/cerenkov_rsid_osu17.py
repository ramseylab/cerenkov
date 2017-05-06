import sys
import pandas as pd
import chrom_tool as ct
from genome_browser_client import GenomeBrowserClient
from allele_util import AlleleUtil
from hgmd_client import HgmdClient
from coord_util import CoordUtil

if sys.version_info[0] == 3:
    # Using Python3
    from snap_client import SnapQuery
else:
    # Using Python2
    from snap_client_py2 import SnapQuery


def __merge_coord_and_allele(snp_coord, snp_allele):
    # We filter on allele, so we apply left inner join of "allele x coord"
    return snp_allele.merge(snp_coord, how='inner', on=['name', 'chrom'], copy=True)


def __get_coord(rsid, db_config_key):
    coord_util = CoordUtil()
    coord_util.db_config_key = db_config_key

    snp_coord = coord_util.extract(_input=rsid)

    return snp_coord


def __faulty_filter_on_allele(rsid, db_config_key):
    maf_thld = 0.05

    with GenomeBrowserClient(db_config_key) as gb_client:
        snp_allele = gb_client.fetch_alleles(rsid)
        snp_allele = ct.remove_dup_on_chrY(snp_allele)

    snp_allele = AlleleUtil.transform_cols(snp_allele)

    n_allele = snp_allele.loc[:, list("ATCG")].apply(lambda x: sum(x > 0), axis=1)

    # all mono-allelic are excluded
    # mono = (n_allele == 1)

    # all bi-allelic are included
    # The problem of the previous filter: we didn't apply the "MAF >= 5%" rule to biallelic SNPs
    bi = (n_allele == 2)

    # include if min freq < 0.05 (then we can just discard this min freq);
    # This is faulty because we may include an entry with freq = (0.96, 0.02, 0.02)
    tri = (n_allele == 3)

    # include if min 2 freqs < 0.05 (then we can just discard these 2 min freqs);
    # This is faulty because we may include an entry with freq = (0.94, 0.02, 0.02, 0.02)
    quad = (n_allele == 4)

    bi_dfm = snp_allele.loc[bi, :]

    tri_dfm = snp_allele.loc[tri, :]
    min_one_under_thld = tri_dfm.loc[:, list("ATCG")].apply(lambda x: min(x[x > 0]) < maf_thld, axis=1, reduce=True)
    tri_dfm = tri_dfm.loc[min_one_under_thld, :]

    quad_dfm = snp_allele.loc[quad, :]
    min_two_under_thld = quad_dfm.loc[:, list("ATCG")].\
        apply(lambda x: (x[x > 0].sort_values().iloc[[0, 1]] < maf_thld).all(), axis=1, reduce=True)
    quad_dfm = quad_dfm.loc[min_two_under_thld, :]

    # No need to query Biomart because it's known that what Biomart would return is empty

    # And PCE exclusion is also done after `get_df_1kb`

    return pd.concat([bi_dfm, tri_dfm, quad_dfm], axis=0)


def __get_locus_map(r_rsid):
    query_1kg = SnapQuery().dataset('onekgpilot'). \
        population(['CEU', 'YRI', 'CHBJPT']). \
        r_squared_threshold(0.8). \
        distance_limit_in_kb(250). \
        snp(r_rsid)

    locus_1kg = query_1kg.execute(verbose=True)

    query_rel21 = SnapQuery().dataset('rel21'). \
        population(['CEU', 'YRI', 'JPT+CHB']). \
        r_squared_threshold(0.8). \
        distance_limit_in_kb(250). \
        snp(r_rsid)

    locus_rel21 = query_rel21.execute(verbose=True)
    # Use entries only in Chromosome X of HapMap release 21
    locus_rel21 = locus_rel21.loc[locus_rel21['Chromosome'] == 'chrX']

    locus_map = pd.concat([locus_1kg, locus_rel21], ignore_index=True, axis=0).\
        drop_duplicates(['SNP', 'Proxy'], keep='last')
    locus_map = locus_map[~locus_map['Proxy'].isin(r_rsid)]  # remove rSNPs in cSNPs

    return locus_map


def __pce_filter(snp_coord):
    return CoordUtil.pce_filter(snp_coord, verbose=True)


def __faulty_1kb_window(r_coord, c_coord, locus_map):
    """
    previously `getAllDataDF1kb()`
    """
    c_coord = c_coord.assign(keep=0)
    for _, rsnp in r_coord.iterrows():
        locus = locus_map.loc[locus_map['SNP'] == rsnp['name'], 'Proxy'].tolist()
        c_coord.loc[c_coord['name'].isin(locus) &
                    (c_coord['chromStart'] - rsnp['chromStart'] <= 1000), 'keep'] = 1

    result = c_coord.loc[c_coord['keep'] == 1, :].drop_duplicates()
    result = result.drop('keep', axis=1)
    return result


def reproduce_recomb_rsid(db_config_key, local_locus_map_file=None):
    with HgmdClient() as hgmd_client:
        r_rsid = hgmd_client.select_rsid()  # a list

    if local_locus_map_file is None:
        # locus_map is generated with the unfiltered RSNP rsid
        locus_map = __get_locus_map(r_rsid)
    else:
        locus_map = pd.read_table(local_locus_map_file, header=0)

    c_rsid = locus_map.loc[:, 'Proxy'].unique().tolist()
    c_rsid = list(set(c_rsid) - set(r_rsid))

    r_allele = __faulty_filter_on_allele(r_rsid, db_config_key)
    r_coord = __get_coord(rsid=r_rsid, db_config_key='local_hg19')
    r_dfm = __merge_coord_and_allele(r_coord, r_allele)

    c_allele = __faulty_filter_on_allele(c_rsid, db_config_key)
    c_coord = __get_coord(rsid=c_rsid, db_config_key='local_hg19')
    c_dfm = __merge_coord_and_allele(c_coord, c_allele)

    c_dfm = __faulty_1kb_window(r_dfm, c_dfm, locus_map)

    r_dfm = __pce_filter(r_dfm)
    c_dfm = __pce_filter(c_dfm)

    return r_dfm, c_dfm, locus_map


def missing_maf_filter(snp_allele):
    maf_thld = 0.05
    snp_valid = snp_allele.loc[:, list('ATCG')].apply(lambda x: sum(x >= maf_thld) == 2, axis=1, reduce=True)

    return snp_allele.loc[snp_valid]


def missing_50kb_window(r_coord, c_coord, locus_map):
    c_coord = c_coord.assign(keep=0)

    for _, rsnp in r_coord.iterrows():
        locus = locus_map.loc[locus_map['SNP'] == rsnp['name'], 'Proxy'].tolist()
        # keep the 1kb window SNPs
        c_coord.loc[c_coord['name'].isin(locus) &
                    ((c_coord['chromStart'] - rsnp['chromStart']).abs() <= 50000), 'keep'] = 1

    result = c_coord.loc[c_coord['keep'] == 1, :].drop_duplicates()
    result = result.drop('keep', axis=1)
    return result


if __name__ == '__main__':
    r_recomb_dfm, c_recomb_dfm, locus_map = reproduce_recomb_rsid(db_config_key='local_hg19',
                                                                  local_locus_map_file='Locus_Map.tsv')

    c_dfm = missing_maf_filter(c_recomb_dfm)
    r_dfm = missing_maf_filter(r_recomb_dfm)

    c_dfm = missing_50kb_window(r_dfm, c_dfm, locus_map)

    r_dfm.to_csv('OSU17_rsnp.tsv', header=True, index=False, sep='\t')
    c_dfm.to_csv('OSU17_csnp.tsv', header=True, index=False, sep='\t')

    # print(r_dfm.shape)
    # print(c_dfm.shape)

    print("[cerenkov_rsid_osu17] Done!")

# import argparse

# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(description="Fetch rsid from HGMD/SNAP.", allow_abbrev=False)
#
#     parser.add_argument('-r', '--rsnp', dest='r_dest', type=str, required=True,
#                         help="output file of rSNP rsid")
#     parser.add_argument('-c', '--csnp', dest='c_dest', type=str, required=True,
#                         help="output file of cSNP rsid")
#
#     args = parser.parse_args()
#
#     print("[cerenkov_rsid_osu17] rSNP rsid output: {}; cSNP rsid output: {}".
#           format(args.r_dest, args.c_dest))
#
#     # Original RSNP rsid from HGMD
#     with HgmdClient() as hgmd_client:
#         r_rsid = hgmd_client.select_rsid()  # a list
#
#
#
#
#     # RSNP
#     pd.Series(data=r_rsid, name='name').to_csv(args.r_dest, header=True, index=False, sep='\t')
#
#     # CSNP
#     pd.Series(data=c_rsid, name='name').to_csv(args.c_dest, header=True, index=False, sep='\t')
#
#     print("[cerenkov_rsid_osu17] Done!")



