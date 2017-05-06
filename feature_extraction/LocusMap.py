import pandas as pd


def extract_locus_map():
    """
    Requirements:  SNAPTool results and rSNPs list (rsnps.rsID.txt)
    Output: rSNPs and LD corresponding cSNPs Mapping (export to .csv file)
    """

    # read in SNAPTool results data and rsnp.id
    ceu_1000_genome = pd.read_csv('CEU_1000Genome.txt', sep='\t')
    ceu_hapmap_21 = pd.read_csv('CEU_HapMap21.txt', sep='\t')
    yri_1000_genome = pd.read_csv('YRI_1000Genome.txt', sep='\t')
    yri_hapmap_21 = pd.read_csv('YRI_HapMap21.txt', sep='\t')
    chbjpt_1000_genome = pd.read_csv('CHBJPT_1000Genome.txt', sep='\t')
    chbjpt_hapmap_21 = pd.read_csv('CHBJPT_HapMap21.txt', sep='\t')
    rSNPs = pd.read_csv('rsnps.rsID.txt', sep='\t')

    # merge according to populations
    ceu = pd.concat([ceu_1000_genome, ceu_hapmap_21[ceu_hapmap_21['Chromosome'] == 'chrX']]).drop_duplicates(
        ['SNP', 'Proxy'], keep='last')
    ceu = ceu[~ceu['Proxy'].isin(rSNPs['name'])]  # remove rSNPs in CEU
    yri = pd.concat([yri_1000_genome, yri_hapmap_21[yri_hapmap_21['Chromosome'] == 'chrX']]).drop_duplicates(
        ['SNP', 'Proxy'], keep='last')
    yri = yri[~yri['Proxy'].isin(rSNPs['name'])]  # remove rSNPs in YRI
    chbjpt = pd.concat(
        [chbjpt_1000_genome, chbjpt_hapmap_21[chbjpt_hapmap_21['Chromosome'] == 'chrX']]).drop_duplicates(
        ['SNP', 'Proxy'], keep='last')
    chbjpt = chbjpt[~chbjpt['Proxy'].isin(rSNPs['name'])]  # remove rSNPs in CHBJPT

    locus_map = pd.concat([ceu, yri, chbjpt]).drop_duplicates(['SNP', 'Proxy'], keep='last')
    locus_map = locus_map[~locus_map['Proxy'].isin(rSNPs['name'])]  # remove rSNPs in cSNPs

    # print size of each population
    print('CEU:' % ceu.shape)
    print('YRI:' % yri.shape)
    print('CHBJPT:' % chbjpt.shape)
    print('locusMap:' % locus_map.shape)

    # export .csv files
    locus_map[['SNP', 'Proxy']].to_csv('../data/SNAP_proxy/locusMap.txt', sep='\t', Header=True)
    print("export data successfully!")
