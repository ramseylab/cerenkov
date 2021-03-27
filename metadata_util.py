import os
import pandas
import numpy
from biomart_client import BiomartClient
from genome_browser_client import GenomeBrowserClient
import chrom_tool as CT
from allele_tool import flip_allele, build_allele_freq_map


def __read_id(src):
    if not os.path.isfile(src):
        raise OSError(src + " not found.")

    rsid = pandas.read_csv(src, sep='\t').loc[:, "name"]

    return rsid


def __to_csv(dfm, dest):
    dfm.to_csv(dest, sep='\t', header=True)


def __to_bed(dfm, dest):
    dfm.to_csv(dest, sep='\t', header=False, index=False, columns=['chrom', 'chromStart', 'chromEnd', 'name'])


def __print_summary(snp_dfm):
    print("\nTotal # of SNPs: " + str(snp_dfm.shape[0]))

    print("\n# of SNPs on different chromosome types: \n")

    chrom_types = ["Regular Chromosomes", "Mito Sequence (chrM)", "Haplotype Chromosomes", "Unplaced Contigs", "Sum"]

    chrom_series = snp_dfm["chrom"]

    cnt_regular = sum(chrom_series.isin(CT.REGULAR_CHR))
    cnt_mito_seq = sum(chrom_series == CT.HOMO_SAPIENS_MITO_SEQ)
    cnt_haplotype = sum(chrom_series.isin(CT.HAPLOTYPE_CHR))
    cnt_contig = sum(chrom_series.str.endswith(CT.CONTIG_LOCALIZED_SUFFIX) |
                     chrom_series.str.startswith(CT.CONTIG_UNKNOWN_PREFIX))
    cnt_total = cnt_regular + cnt_mito_seq + cnt_haplotype + cnt_contig

    counts = [cnt_regular, cnt_mito_seq, cnt_haplotype, cnt_contig, cnt_total]

    summary = pandas.DataFrame(counts, chrom_types, ["count"])

    print(summary)

    print("\n")


def __remove_non_regular_chrom(snp_dfm, verbose=False):
    """
    There are 92 distinct values in column "chrom" of table "snp142" in database "hg19":

        - Regular                                 (24): "chr1" ~ "chr22", "chrX", "chrY"
        - The Homo sapiens mitochondrion sequence (1) : "chrM"
        - Haplotype chromosomes                   (9) : E.g. "chr6_apd_hap1"
        - Unplaced contigs                        (58 = 20 + 38):
            - If an unplaced contig is localized to a chromosome,
                the contig name is appended to the regular chromosome name, as in "chr1_gl000191_random".
            - If the chromosome is unknown,
                the contig is represented with the name "chrUn" followed by the contig identifier,
                as in "chrUn_gl000211".

    See http://hgdownload.cse.ucsc.edu/gbdb/hg19/html/description.html.

    E.g. for cSNP "rs1130552", there are 1 record on regular "chr6" and 6 on haplotype chromosomes.
        We only keep the record on "chr6".

         name          chrom
    rs1130552           chr6
    rs1130552  chr6_cox_hap2
    rs1130552  chr6_dbb_hap3
    rs1130552 chr6_mann_hap4
    rs1130552  chr6_mcf_hap5
    rs1130552  chr6_qbl_hap6
    rs1130552 chr6_ssto_hap7
    """

    is_regular = snp_dfm.loc[:, "chrom"].isin(CT.REGULAR_CHR)

    if verbose:
        print("__remove_non_regular_chrom: \r\n%s" % snp_dfm.loc[~is_regular, ["name", "chrom"]])

    return snp_dfm.loc[is_regular, :]


def __remove_non_single_class(snp_dfm, verbose=False):
    """
    We do not care about insertions, deletions, in-dels at present.
    Just keep the "single" class.
    """

    is_single = (snp_dfm["class"] == "single")

    if verbose:
        print("__remove_non_single_class: \r\n%s" % snp_dfm.loc[~is_single, ["name", "chrom"]])

    return snp_dfm.loc[is_single, :]


def __normalize_allele_strand(snp_dfm):
    """
    Keep all the alleles on FWD strand.
    If `strand` is "-", flip every base in `alleles`; otherwise do not change `alleles`.
    """

    on_rev = (snp_dfm.loc[:, "strand"] == "-")
    has_alleles = (snp_dfm.loc[:, "alleles"].str.len() > 0)
    condition = (on_rev & has_alleles)

    if not snp_dfm.loc[condition, :].empty:
        snp_dfm.loc[condition, "alleles"] = snp_dfm.loc[condition, "alleles"].apply(flip_allele)

    return snp_dfm


def __build_allele_freq_map(snp_dfm):

    snp_dfm.loc[:, "afMap"] = snp_dfm.apply(lambda row: build_allele_freq_map(row['alleles'], row['alleleFreqs']),
                                            axis=1, reduce=True)

    return snp_dfm


def __identify_major_minor_alleles(snp_dfm, verbose=False):
    """
    Policies:

    P-1. If allele frequency info are complete:
        P-1-1. Mono-allelic: Query ensembl.
        P-1-2. Bi-allelic:
            P-1-2-1. If MAF=0, query ensembl;
            P-1-2-2. otherwise calculate normally
        P-1-3. Tri-allelic: Delete a allele with minimum freq if its freq < 0.05 and then treat it
            as a bi-allelic one; discard otherwise
        P-1-4. Quad-allelic: Delete 2 alleles with minimum freqs if both their freq < 0.05 and then treat it
            as a bi-allelic one; discard otherwise
    P-2. If allele frequency info are missing: Query ensembl.

    Actions after querying ensembl:

    A-1. If ensembl claimed that it is mono-allelic, discard;
    A-2. If ensembl claimed that it is bi-allelic,
        A-2-1. If MAF=0, discard;
        A-2-2. otherwise calculate normally
    A-3. If ensembl claimed that it is tri-allelic or quad-allelic, discard
        because we cannot tell which is the major allele when ensembl only reports a minor one.
    """

    maf_thld = 0.05  # MAF threshold

    undetermined = pandas.DataFrame()  # rsID in this data frame should be double-checked thru ensembl via biomart
    excluded = pandas.DataFrame()  # rsID in this data frame should excluded. just for information display
    result = pandas.DataFrame()  # rsID in this data frame is returned

    # P-2
    has_no_afmap = (snp_dfm.loc[:, "afMap"].str.len() == 0)
    undetermined = undetermined.append(snp_dfm[has_no_afmap], ignore_index=True)

    # Identify mono-, bi-, tri-, quad-allelic
    afmap_len = snp_dfm.loc[:, "afMap"].str.len()
    mono_allelic = snp_dfm[afmap_len == 1]
    bi_allelic = snp_dfm[afmap_len == 2]
    tri_allelic = snp_dfm[afmap_len == 3]
    quad_allelic = snp_dfm[afmap_len == 4]

    if verbose:
        allele_types = ["Mono-allelic", "Bi-allelic", "Tri-allelic", "Quad-allelic", 'Total']

        cnt_mono = mono_allelic.shape[0]
        cnt_bi = bi_allelic.shape[0]
        cnt_tri = tri_allelic.shape[0]
        cnt_quad = quad_allelic.shape[0]
        cnt_total = cnt_mono + cnt_bi + cnt_tri + cnt_quad
        counts = [cnt_mono, cnt_bi, cnt_tri, cnt_quad, cnt_total]

        summary = pandas.DataFrame(counts, allele_types, ["count"])
        print(summary)

    # P-1-1
    if not mono_allelic.empty:
        undetermined = undetermined.append(mono_allelic)

    # P-1-3
    if not tri_allelic.empty:
        # An entry from "afMap" is like `[("A", 0.2), ("C", 0.3), ("G", 0.5)]`
        # x[0][1] == 0.2 in this example
        # If two smallest freqs are both greater than 5%, you cannot tell which should be the minor allele
        has_ambig_maf = numpy.array([x[0][1] >= maf_thld
                                     for x in tri_allelic.loc[:, "afMap"]])
        ambig_tri_allelic = tri_allelic.loc[has_ambig_maf, :]
        excluded = excluded.append(ambig_tri_allelic, ignore_index=True)

        if verbose:
            print("__identify_major_minor_alleles: {n} tri-allelic entries excluded \r\n{entries}".format(
                n=sum(has_ambig_maf),
                entries=ambig_tri_allelic.loc[:, ["name", "alleles", "afMap"]] if any(has_ambig_maf) else ""))

        remainder = tri_allelic.loc[~has_ambig_maf, :].copy()
        # delete the first element
        remainder.loc[:, "afMap"] = remainder.loc[:, "afMap"].apply(lambda x: x[1:])

        bi_allelic = bi_allelic.append(remainder, ignore_index=True)

    # P-1-4
    if not quad_allelic.empty:
        has_ambig_maf = numpy.array([(x[0][1] >= maf_thld) or (x[1][1] >= maf_thld)
                                     for x in quad_allelic.loc[:, "afMap"]])
        ambig_quad_allelic = quad_allelic.loc[has_ambig_maf, :]
        excluded = excluded.append(ambig_quad_allelic, ignore_index=True)

        if verbose:
            print("__identify_major_minor_alleles: {n} quad-allelic entries excluded \r\n{entries}".format(
                n=sum(has_ambig_maf),
                entries=ambig_quad_allelic.loc[:, ["name", "alleles", "afMap"]] if any(has_ambig_maf) else ""))

        remainder = quad_allelic.loc[~has_ambig_maf, :].copy()
        # delete the first 2 elements
        remainder.loc[:, "afMap"] = quad_allelic.loc[:, "afMap"].apply(lambda x: x[2:])

        bi_allelic = bi_allelic.append(remainder, ignore_index=True)

    # P-1-2
    if not bi_allelic.empty:
        # P-1-2-1
        freq_eq_zero = numpy.array([(x[0][1] == 0.0) or (x[1][1] == 1.0) for x in bi_allelic.loc[:, "afMap"]])
        undetermined = undetermined.append(bi_allelic.loc[freq_eq_zero, :])

        # P-1-2-2
        remainder = bi_allelic.loc[~freq_eq_zero, :].copy()

        remainder.loc[:, "minorAllele"] = remainder.loc[:, "afMap"].apply(lambda x: x[0][0])
        remainder.loc[:, "minorAlleleFreq"] = remainder.loc[:, "afMap"].apply(lambda x: x[0][1])
        remainder.loc[:, "majorAllele"] = remainder.loc[:, "afMap"].apply(lambda x: x[1][0])
        remainder.loc[:, "majorAlleleFreq"] = remainder.loc[:, "afMap"].apply(lambda x: x[1][1])

        result = result.append(remainder, ignore_index=True)

    if not undetermined.empty:  # list not empty
        with BiomartClient() as bm_client:
            response = bm_client.query_snp(undetermined.loc[:, "name"].tolist(), verbose)
            determined = response.loc[:, ["name", "minorAllele", "minorAlleleFreq", "majorAllele", "majorAlleleFreq"]].\
                merge(undetermined, on='name', how='left', left_index=True)

            result = result.append(determined, ignore_index=True)

    if verbose:
        maf_le_thld = result.loc[result["minorAlleleFreq"] < maf_thld]
        print("__identify_major_minor_alleles: applied 5%-MAF filter to {n} entries: \r\n{entries}".format(
            n=maf_le_thld.shape[0],
            entries=maf_le_thld.loc[:, ["name", "alleles", "afMap"]] if not maf_le_thld.empty else ""))

    return result.loc[result["minorAlleleFreq"] >= maf_thld]


def __revise_alleles_with_equal_freqs(snp_dfm):
    """
    There exist SNPs with majorAlleleFreq == minorAlleleFreq.

    Steve:

    > For those 57 SNPs, I would maybe use the following approach:
    > if the last (least significant) digit of the chromosomal coordinate is even, select the first (out of two)
        alleles listed as the minor allele.
    > If the last (least significant) digit of the chromosomal coordinate is odd, select the second (out of two)
        alleles listed as the minor allele.
    > This approach is deterministic but should more or less "randomly" and with equal probability assign the
        first or second allele to be the 'minor allele'.
    """

    has_equal_freq = (snp_dfm.loc[:, "minorAlleleFreq"] == snp_dfm.loc[:, "majorAlleleFreq"])

    print("__revise_alleles_with_equal_freqs: detected \r\n %s" %
          snp_dfm.loc[has_equal_freq, ["name", "chrom", "chromStart", "chromEnd",
                                       "majorAllele", "minorAllele", "majorAlleleFreq", "minorAlleleFreq"]])

    has_odd_chrom_start = (snp_dfm.loc[:, "chromStart"] % 2 == 1)

    idx = has_equal_freq & has_odd_chrom_start

    print("__revise_alleles_with_equal_freqs: to swap \r\n %s" %
          snp_dfm.loc[idx, ["name", "chrom", "chromStart", "chromEnd",
                            "majorAllele", "minorAllele", "majorAlleleFreq", "minorAlleleFreq"]])

    # Swap
    snp_dfm.loc[idx, ["majorAllele", "minorAllele"]] = snp_dfm.loc[idx, ["minorAllele", "majorAllele"]].values

    return snp_dfm


def __drop_redundant_col(snp_dfm):
    return snp_dfm.drop(['afMap', "alleleFreqs", "alleles"], axis=1)


def __normalize_chrom_coord(snp_dfm):
    snp_dfm.loc[:, "normChromCoord"] = snp_dfm.apply(lambda row: row['chromStart'] / CT.CHR_LENGTH[row['chrom']],
                                                     axis=1)

    return snp_dfm


def extract_metadata(src, dest_csv, dest_bed):
    """
    Extract metadata for cSNPs or rSNPs by querying UCSC database

    :param src: the rsID list
    :param dest_csv: the feature matrix
    :param dest_bed: the name of bed file to be generated
    :return: None
    """

    rsid = __read_id(src)

    with GenomeBrowserClient('local_hg19') as gb_client:
        snps = gb_client.fetch_metadata(rsid)

        __print_summary(snps)

        snps = __remove_non_regular_chrom(snps, verbose=True)
        snps = __remove_non_single_class(snps, verbose=True)
        snps = __normalize_allele_strand(snps)
        snps = __build_allele_freq_map(snps)
        snps = __identify_major_minor_alleles(snps, verbose=True)
        snps = __revise_alleles_with_equal_freqs(snps)
        snps = __drop_redundant_col(snps)
        snps = __normalize_chrom_coord(snps)
        snps = CT.remove_dup_on_chrY(snps)

        snps = snps.set_index("name")
        __to_csv(snps, dest_csv)

        snps = snps.reset_index()
        __to_bed(snps, dest_bed)
