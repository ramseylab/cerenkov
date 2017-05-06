from biomart import BiomartServer
from io import StringIO
import pandas as pd
import numpy as np
import chrom_tool as ct
from allele_tool import flip_allele, parse_bm_alleles_col

pd.set_option('display.width', 1000)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)


class AbstractBiomartClient(object):
    def __init__(self):
        # GRCh37 is also known as hg19

        # server = BiomartServer("http://useast.ensembl.org/biomart")
        self.server = BiomartServer("http://grch37.ensembl.org/biomart")

        # set verbose to True to get some messages
        # server.verbose = True

        # server.show_databases()

        self.database = self.server.databases["ENSEMBL_MART_SNP"]

        # db.show_datasets()

        self.dataset = self.database.datasets["hsapiens_snp"]

        # dataset.show_filters()
        # dataset.show_attributes()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        del self.dataset
        del self.database
        del self.server

    def query_snp(self, rsid_list, verbose=False):
        if verbose:
            print("query_snp parameter: %s" % rsid_list)

        response = self.dataset.search({
            'filters': {
                'snp_filter': rsid_list
            },
            'attributes': [
                'refsnp_id', 'chr_name', 'chrom_strand', 'allele', "minor_allele", "minor_allele_freq"
            ]
        }, header=1)

        lines = [line.decode('utf-8') for line in response.iter_lines()]

        tsv_str = StringIO("\r\n".join(lines))

        dfm = pd.read_csv(tsv_str, header=0, sep="\t",
                              names=["name", "chrom", "strand", "alleles", "minorAllele", "minorAlleleFreq"])

        dfm.loc[:, 'alleles'] = dfm['alleles'].apply(parse_bm_alleles_col)

        flt_dfm = self.filter(dfm)
        tx_result = self.transform(flt_dfm)

        if verbose:
            print("query_snp parameter: %s" % rsid_list)
            print("query_snp raw result: \r\n %s" % dfm)
            print("query_snp filtered result: \r\n %s" % flt_dfm)
            print("query_snp transformed result: \r\n %s" % tx_result)

        return tx_result

    def filter(self, dfm):
        raise NotImplementedError("This is an abstract method.")

    def transform(self, dfm):
        raise NotImplementedError("This is an abstract method.")


class BiomartClient(AbstractBiomartClient):
    def filter(self, dfm):
        """
        Filter out non-biallelic entries.
        """

        is_biallelic = (dfm["alleles"].str.len() == 2)
        has_minor_allele = dfm["minorAllele"].notnull()
        has_valid_maf = (dfm["minorAlleleFreq"].notnull()) & \
                        (dfm["minorAlleleFreq"] != 0.0) & \
                        (dfm["minorAlleleFreq"] != 1.0)
        has_valid_chrom = dfm.loc[:, "chrom"].astype(str).isin(ct.REGULAR_CHR_ID)

        return dfm.loc[is_biallelic & has_minor_allele & has_valid_maf & has_valid_chrom, :]

    def transform(self, dfm):
        on_rev = (dfm.loc[:, "strand"] == -1)

        if not dfm.loc[on_rev, :].empty:
            dfm.loc[on_rev, "alleles"] = dfm.loc[on_rev, "alleles"].apply(flip_allele)
            dfm.loc[on_rev, "minorAllele"] = dfm.loc[on_rev, "minorAllele"].apply(flip_allele)

        dfm = dfm.assign(majorAllele=dfm.apply(lambda x: [a for a in x["alleles"] if a != x["minorAllele"]][0],
                                               axis=1, reduce=True))
        dfm = dfm.assign(majorAlleleFreq=1-dfm["minorAlleleFreq"])

        return dfm.loc[:, ["name", "alleles", "minorAllele", "minorAlleleFreq", "majorAllele", "majorAlleleFreq"]]


class BiomartClient2(AbstractBiomartClient):
    def filter(self, dfm):
        """
        Filter out non-biallelic entries.
        """

        is_biallelic = (dfm["alleles"].str.len() == 2)
        has_minor_allele = dfm["minorAllele"].notnull()
        has_valid_maf = (dfm["minorAlleleFreq"].notnull()) & \
                        (dfm["minorAlleleFreq"] != 0.0) & \
                        (dfm["minorAlleleFreq"] != 1.0)
        has_valid_chrom = dfm.loc[:, "chrom"].astype(str).isin(ct.REGULAR_CHR_ID)

        result = dfm.loc[is_biallelic & has_minor_allele & has_valid_maf & has_valid_chrom, :]

        return result

    def transform(self, dfm):

        # change 'strand' representation
        orig_strand = dfm.loc[:, 'strand']
        dfm.loc[orig_strand == 1, 'strand'] = "+"
        dfm.loc[orig_strand == -1, 'strand'] = "-"

        # change 'chrom' representation
        dfm.loc[:, 'chrom'] = 'chr' + dfm.loc[:, 'chrom'].astype(str)

        # add 'ATCG' columns
        dfm = dfm.assign(A=np.nan, T=np.nan, C=np.nan, G=np.nan)

        def _make_allele_dict(allele_lst, minor_allele, minor_freq):
            major_allele = [a for a in allele_lst if a != minor_allele][0]
            major_freq = 1.0 - minor_freq
            return {major_allele: major_freq, minor_allele: minor_freq}

        dfm = dfm.apply(lambda x: x.fillna(_make_allele_dict(x['alleles'],
                                                             x['minorAllele'],
                                                             x['minorAlleleFreq'])),
                        axis=1)

        return dfm.drop(['alleles', "minorAllele", "minorAlleleFreq"], axis=1)

