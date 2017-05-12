"""
Created on Jul 11, 2016

@author: ramseylab
"""

REGULAR_CHR_ID = set([str(x) for x in range(1, 23)] + ["X", "Y"])  # range(1,23) = 1,2,...,22
REGULAR_CHR = {("chr" + str(x)) for x in REGULAR_CHR_ID}
HOMO_SAPIENS_MITO_SEQ = "chrM"
HAPLOTYPE_CHR = {
    "chr6_apd_hap1", 
    "chr6_cox_hap2", 
    "chr6_dbb_hap3", 
    "chr6_mann_hap4", 
    "chr6_mcf_hap5", 
    "chr6_qbl_hap6", 
    "chr6_ssto_hap7", 
    "chr4_ctg9_hap1", 
    "chr17_ctg5_hap1"
}
CONTIG_LOCALIZED_SUFFIX = "random"
CONTIG_UNKNOWN_PREFIX = "chrUn"

# Assembly Statistics for GRCh38.p7 
# Release date: March 21, 2016
# See http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/
CHR_LENGTH = {
    'chr1': 248956422,
    'chr2': 243199373,
    'chr3': 198295559,
    'chr4': 190214555,
    'chr5': 181538259,
    'chr6': 170805979,
    'chr7': 159345973,
    'chr8': 145138636,
    'chr9': 138394717,
    'chr10': 133797422,
    'chr11': 135086622,
    'chr12': 133275309,
    'chr13': 114364328,
    'chr14': 107043718,
    'chr15': 101991189,
    'chr16': 90338345,
    'chr17': 83257441,
    'chr18': 80373285,
    'chr19': 58617616,
    'chr20': 64444167,
    'chr21': 46709983,
    'chr22': 50818468,
    'chrX': 156040895,
    'chrY': 57227415
}


def remove_dup_on_chrY(snp_dfm, verbose=True):
    """
        It is known that each of "rs700442", "rs5989681" and "rs4446909" has 2 entries in
        snp146, hg19 of UCSC Genome Browser, one on chrX and the other chrY.
        
        Remove the one on chrY as suggested by Steve:
        
        > I would just use the ChrX entry, for this SNP.
        > My reasoning:  we very likely don't have haplotype data for chrY.
    """
    is_dup = snp_dfm.duplicated(subset="name", keep=False)
    on_chrY = (snp_dfm.loc[:, "chrom"] == "chrY")
    
    if verbose and is_dup.any():
        print("remove_dup_on_chrY: duplication detected: \r\n%s" % snp_dfm.loc[is_dup, :])
    
    snp_dfm = snp_dfm.drop(snp_dfm[is_dup & on_chrY].index)
    
    still_dup = snp_dfm.duplicated(subset="name", keep=False)
    assert not still_dup.any(), "remove_dup_on_chrY: still duplicated after removal: \r\n%s" % snp_dfm.loc[still_dup, :]
    
    return snp_dfm
