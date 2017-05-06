import argparse
import rsid_util
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Fetch rsid from HGMD/SNAP.", allow_abbrev=False)

    parser.add_argument('-r', '--rsnp', dest='r_dest', type=str, required=True,
                        help="output file of rSNP rsid")
    parser.add_argument('-c', '--csnp', dest='c_dest', type=str, required=True,
                        help="output file of cSNP rsid")

    args = parser.parse_args()

    print("[cerenkov_rsid] rSNP rsid output: {}; cSNP rsid output: {}.".
          format(args.r_dest, args.c_dest))

    # RSNP
    r_rsid = rsid_util.extract_rsnp()
    pd.Series(data=r_rsid, name='name').to_csv(args.r_dest, header=True, index=False, sep='\t')

    # CSNP
    c_rsid = rsid_util.extract_csnp(r_rsid)
    pd.Series(data=c_rsid, name='name').to_csv(args.c_dest, header=True, index=False, sep='\t')

