"""
Created on Jul 11, 2016

@author: ramseylab
"""

__ALTERNATE_NUCLEOBASE = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C"
}


def parse_bm_alleles_col(_str):
    """
    Parse 'alleles' column of biomart response
    """
    if _str == '':
        return []
    else:
        return _str.split('/')


def parse_gb_alleles_col(_bytes):
    """
    Parse 'alleles' column of table snp146 in UCSC genome browser
    """
    if _bytes == b'':
        return []
    else:
        # "A,C,".split(',') = ['A', 'C', '']
        # Remove the last empty string
        return _bytes.decode("utf-8").split(',')[:-1]


def parse_gb_allele_freqs_col(_bytes):
    """
    Parse 'alleleFreqs' column of table snp146 in UCSC genome browser
    """
    if _bytes == b'':
        return []
    else:
        # "0.2,0.8,".split(',') = ['0.2', '0.8', '']
        # Remove the last empty string
        freq_str = _bytes.decode("utf-8").split(',')[:-1]
        freq_float = [float(i) for i in freq_str]
        return freq_float


def flip_allele(allele_seq):
    """
    For alleles (e.g. ['A','G']) on REV strand, return the alternate representation (e.g. ['T','C']) on FWD strand.
    """
    complement = [__ALTERNATE_NUCLEOBASE[x] for x in allele_seq]

    if isinstance(allele_seq, str):
        return ''.join(complement)
    else:
        return complement


def build_allele_freq_map(allele_seq, freq_seq):
    """
    Construct a <allele, freq> map, sorted by <freq> ascendingly.

    0 in `allele_seq`, if exists, will be removed.

    Returned data structure is a list of tuples, e.g.
        getAlleleFreqMap(['A','G'], [0.2,0.8]) => [('A', 0.2), ('G', 0.8)]
    """
    
    if len(allele_seq) == 0 or len(freq_seq) == 0:
        # if `allele_seq` or `freq_seq` is empty list
        return [] 

    if '0' in allele_seq:
        zero_index = allele_seq.index('0')
        allele_seq.pop(zero_index)
        freq_seq.pop(zero_index)

    af_map = list(zip(allele_seq, freq_seq))
    
    # import operator
    # afMap = sorted(afMap, key=operator.itemgetter(1))
    af_map = sorted(af_map, key=lambda x: x[1])
    
    return af_map
