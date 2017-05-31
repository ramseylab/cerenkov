# ----- Command Line Arguments -----
args = commandArgs(trailingOnly=TRUE)
helpMsg <- "Example Invocation: Rscript --vanilla augment_osu_features_v2.R ../source_data/augment_osu_features_datafiles SNP_input.tsv SNP_output.tsv"

if (length(args) < 3) {
  stop(helpMsg, call.=FALSE)
}

src_data_dir <- args[1]
input_fn <- args[2]  ## E.g. 'alldata_osu.txt'
output_fn <- args[3]  ## E.g. 'osudata_recomb.txt'




## ----- Source File Paths -----

gencode_tss_bed_path <- file.path(src_data_dir, "gencode_19_hg19_TSS.bed")
cpg_island_bed_path <- file.path(src_data_dir, "cpgIslandExt.bed")
jaspar_tfbs_table_path <- file.path(src_data_dir, "jaspar_tfbs_ensembl_75_hg19.txt")
het_bed_path <- file.path(src_data_dir, "het_rates_filtered.bed")
daf_bed_path <- file.path(src_data_dir, "daf_filtered.bed")
simple_repeat_bed_path <- file.path(src_data_dir, "simpleRepeat.bed")
nested_repeat_bed_path <- file.path(src_data_dir, "nestedRepeat.bed")
gene_annot_bed_path <- file.path(src_data_dir, "knownGene_filt.bed")
known_gene_exons_bed_path <- file.path(src_data_dir, "knownGene_exons.bed")




## ----- Libs -----

library(Biostrings)
library(rtracklayer)
library(pbapply)


## ----- Read OSU feature data for all SNPs; fetch BED -----
snps_data <- read.table(input_fn, header=TRUE, stringsAsFactors=TRUE)

snps_bed <- snps_data[c('chrom', 'chromStart', 'chromEnd', 'name')]
snps_granges <- GRanges(snps_bed)

N <- nrow(snps_data)



## ----- Feature 1: `gencode_tss` -----

## read GENCODE data from the hg19 BED file
gencode_TSS_granges <- import.bed(gencode_tss_bed_path, genome="hg19")  # "gencode_19_hg19_TSS.bed"

## compute distance to the nearest gencode TSS
nearest_df <- as.data.frame(distanceToNearest(snps_granges, gencode_TSS_granges))

feature_gencode_TSS_dists <-  pbapply(nearest_df, 1, function(my_row) {
    snp_index <- my_row[1]
    snp_coord <- start(snps_granges[snp_index])
    gencode_index <- my_row[2]
    gencode_unsigned_dist <- my_row[3]
    gencode_range <- gencode_TSS_granges[gencode_index]
    tss_coord <- start(gencode_range)
    gencode_strand <- ifelse(as.character(strand(gencode_range)) == "+", 1, -1)
    my_sign <- 1
    if ((gencode_strand > 0 && snp_coord < tss_coord) ||
         gencode_strand < 0 && snp_coord > tss_coord)  {
        my_sign <- -1
    }
    my_sign*gencode_unsigned_dist
})




## ----- Feature 2: `maf1kb` -----

all_chroms <- as.character(snps_bed$chrom)
unique_chroms <- unique(as.character(all_chroms))
nchrom <- length(unique_chroms)
chrom_rows_list <- lapply(unique_chroms, function(chrom_name) {
    which(all_chroms == chrom_name)
})
names(chrom_rows_list) <- unique_chroms

## compute the average minor allele frequency for SNPs within a 1 kb window
snps_maf <- snps_data[c('chrom', 'chromStart', 'minorAlleleFreq', 'name')]
feature_mean_allele_frequency_1kb <- apply(snps_maf, 1,
                                function(my_row) {
                                    my_chrom <- as.character(my_row[1])
                                    my_coord <- as.integer(my_row[2])
                                    chrom_rows <- chrom_rows_list[[my_chrom]]
                                    inds_near <- chrom_rows[which(abs(snps_maf$chromStart[chrom_rows] - my_coord) <= 500)]
                                    mean(snps_maf$minorAlleleFreq[inds_near])
                                })




## ----- Feature 3: `cpg_island` -----

feature_cpg_island <- rep(0, N)

cpg_island_granges <- import.bed(cpg_island_bed_path)  # "cpgIslandExt.bed"
overlaps_res <- findOverlaps(snps_granges, cpg_island_granges)
feature_cpg_island[queryHits(overlaps_res)] <- 1




## ----- Feature 4: `pwm` -----
## annotate SNPs based on TFBS motifs using JASPAR motif match data from Ensembl BioMart (Ensembl Release 75, archive.ensembl.org)
## A SNP may be mapped to multiple TFBS here. E.g. `mxl-1::mdl-1` in the resulted dataframe means 2 mapped TFBS
## File used here is "jaspar_tfbs_ensembl_75_hg19.txt"
jaspar_data <- read.table(jaspar_tfbs_table_path,
                          header=TRUE,
                          sep="\t",
                          stringsAsFactors=FALSE)

jaspar_data_for_granges <- data.frame(jaspar_data[, 2:5], stringsAsFactors=FALSE)
names(jaspar_data_for_granges) <- c("chrom","chromStart","chromEnd","name")
jaspar_data_for_granges$chrom <- paste("chr", jaspar_data_for_granges$chrom, sep="")

jaspar_granges <- GRanges(jaspar_data_for_granges[,1:3])
overlaps_res <- findOverlaps(snps_granges, jaspar_granges)
feature_jaspar <- rep("", N)
feature_jaspar[queryHits(overlaps_res)] <- jaspar_data_for_granges$name[subjectHits(overlaps_res)]




## ----- Feature 5: `avg_het` -----

## make a GRanges object for 1 kb blocks around each of our SNPs
snps_df_1kb_blocks <- snps_bed
snps_df_1kb_blocks$chromStart <- snps_df_1kb_blocks$chromStart - 500
snps_df_1kb_blocks$chromEnd <- snps_df_1kb_blocks$chromEnd + 500
snps_1kb_blocks_granges <- GRanges(snps_df_1kb_blocks)

## compute average heterozygosity within 1 kb region around our SNPs
het_granges <- import.bed(het_bed_path)  # "het_rates_filtered.bed"
het_overlap <- findOverlaps(snps_1kb_blocks_granges, het_granges)
het_overlap_df <- as.data.frame(het_overlap)
features_avg_het <- pbsapply(1:N, function(my_ind) {
    mean(het_granges$score[het_overlap_df$subjectHits[het_overlap_df$queryHits == my_ind]])
})
features_avg_het[is.na(features_avg_het)] <- mean(na.omit(features_avg_het))




## ----- Feature 6: `avg_daf` -----

## compute the average "derived allele frequency" within 1 kb region around our SNPs
daf_granges <- import.bed(daf_bed_path)  # "daf_filtered.bed"
daf_overlap <- findOverlaps(snps_1kb_blocks_granges, daf_granges)
daf_overlap_df <- as.data.frame(daf_overlap)
features_daf <- pbsapply(1:N, function(my_ind) {
    mean(daf_granges$score[daf_overlap_df$subjectHits[daf_overlap_df$queryHits == my_ind]])
})
features_daf[which(is.na(features_daf))] <- mean(na.omit(features_daf))




## ----- Feature 7: `simplerepeat` -----

## annotate SNPs if they lie within a "simple repeat" using data from UCSC
simple_repeat_bed <- import.bed(simple_repeat_bed_path)  # "simpleRepeat.bed"
simple_repeat_overlap <- findOverlaps(snps_granges, simple_repeat_bed)
feature_simple_repeat <- rep(0, N)
feature_simple_repeat[queryHits(simple_repeat_overlap)] <- 1




## ----- Feature 8: `nestedrepeat` -----

## annotate SNPs if they lie within a "nested repeat" using data from UCSC
nested_repeat_bed <- import.bed(nested_repeat_bed_path)  # "nestedRepeat.bed"
nested_repeat_overlap <- findOverlaps(snps_granges, nested_repeat_bed)
feature_nested_repeat <- rep(0, N)
feature_nested_repeat[queryHits(nested_repeat_overlap)] <- 1




## ----- Feature 9: `geneannot` -----

## create a SNP annotation equivalent of the "EXON", "UTR5", "UTR3", "INTRON", and "CDS" features in GWAVA, based
## on the "knownGene.txt" file from UCSC
gene_annot_granges <- import.bed(gene_annot_bed_path)  # "knownGene_filt.bed"
gene_annot_overlap <- findOverlaps(snps_granges, gene_annot_granges)
feature_gene_annot <- rep("INTERGENIC",N)
feature_gene_annot[queryHits(gene_annot_overlap)] <- gene_annot_granges$name[subjectHits(gene_annot_overlap)]




## ----- Feature 10: `ss_dist` -----

## distance to nearest splice site
known_gene_exons_df <- read.table(known_gene_exons_bed_path,
                                  header=FALSE,
                                  sep="\t",
                                  stringsAsFactors=FALSE)  # "knownGene_exons.bed"
ss1 <- known_gene_exons_df[,c(1,2)]
names(ss1) <- c("chrom","chromStart")
ss2 <- known_gene_exons_df[,c(1,3)]
names(ss2) <- c("chrom","chromStart")
splice_sites_df <- rbind(ss1, ss2)
splice_sites_df2 <- data.frame(splice_sites_df,
                               chromEnd=splice_sites_df[,2]+1,
                               strand=c(known_gene_exons_df[,4],
                                 known_gene_exons_df[,4]))
splice_sites_granges <- GRanges(splice_sites_df2)

splice_nearest_df <- as.data.frame(distanceToNearest(snps_granges, splice_sites_granges))

feature_splice_site_dists <-  pbapply(splice_nearest_df, 1, function(my_row) {
    snp_index <- my_row[1]
    snp_coord <- start(snps_granges[snp_index])
    splice_index <- my_row[2]
    splice_unsigned_dist <- my_row[3]
    splice_range <- splice_sites_granges[splice_index]
    splice_coord <- start(splice_range)
    splice_strand <- ifelse(as.character(strand(splice_range)) == "+", 1, -1)
    my_sign <- 1
    if ((splice_strand > 0 && snp_coord < splice_coord) ||
         splice_strand < 0 && snp_coord > splice_coord)  {
        my_sign <- -1
    }
    my_sign*splice_unsigned_dist
})

## ----- Merge -----

## augment the OSU feature matrix with the new features, eliminating some feature columns that have been expanded
snps_data_final <- cbind(name=snps_data['name'],
                        avg_daf=features_daf,
                        avg_het=features_avg_het,
                        cpg_island=feature_cpg_island,
                        pwm=feature_jaspar,
                        simplerepeat=feature_simple_repeat,
                        nestedrepeat=feature_nested_repeat,
                        gencode_tss=feature_gencode_TSS_dists,
                        maf1kb=feature_mean_allele_frequency_1kb,
                        geneannot=feature_gene_annot,
                        ss_dist=feature_splice_site_dists)


## put chromStart column before chromEnd column
#snps_data_final <- snps_data_final[, c(1,3,2,4:ncol(snps_data_final))]

## write the feature matrix to a file
write.table(snps_data_final,
            file=output_fn,
            sep="\t",
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE)

            
