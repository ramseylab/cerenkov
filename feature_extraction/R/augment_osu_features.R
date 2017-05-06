library(Biostrings)
library(rtracklayer)
library(pbapply)

## read file of SNPs with protein-coding exons
snps_in_exons <- read.table("pce_eliminations.txt",
                            header=TRUE,
                            stringsAsFactors=FALSE)[,1]

## read OSU feature data for all SNPs (this includes a small number of protein-coding SNPs)
osu_data_with_protein_coding_snps <- read.table("alldata_osu.txt",
                     header=TRUE,
                     row.names=1,
                     stringsAsFactors=TRUE)

## eliminate protein-coding SNPs
osu_data <- osu_data_with_protein_coding_snps[which(! rownames(osu_data_with_protein_coding_snps) %in% snps_in_exons),]

## load GWAVA SNPs so we can order OSU SNPs to be the same row order as the SNPs in GWAVA.csv
gwava_data <- read.csv("GWAVA.csv",
                       header=TRUE,
                       row.names=1,
                       stringsAsFactors=TRUE)
osu_data <- osu_data[match(rownames(gwava_data), rownames(osu_data)),  ]

## sanity check that the GWAVA feature-table and the OSU feature-table are in the correct relative row order
stopifnot(all(rownames(osu_data)==rownames(gwava_data)))
stopifnot(all(osu_data$label == gwava_data$label))

## expand TF data to one column per TF
snp_tf_lists <- lapply(as.character(osu_data$tfName), function(snp_tfs) {
    snp_tfs_vec <- NULL
    if (snp_tfs != "0") {
        snp_tfs_vec <- strsplit(snp_tfs, ",")[[1]]
    }
    snp_tfs_vec
})
tfs <- unique(do.call(c, snp_tf_lists))
names(tfs) <- tfs

tf_snp_lists <- lapply(tfs, function(query_tf) {
    sapply(snp_tf_lists,
           function(snp_tfs_vec) {
               ret_value <- 0
               if (! is.null(snp_tfs_vec)) {
                   ret_value <- ifelse(query_tf %in% snp_tfs_vec, 1, 0)
               }
               ret_value
           })
})

## obtain three feature columns for the local DNA sequence data
dna_string_to_int <- function(dna_string) {
    log2(as.integer(DNAString(dna_string)))+1
}
dna_data <- lapply(osu_data$dnaSeq,
                   function(dna_seq) {
                       dna_vec <- dna_string_to_int(dna_seq)
                       dna_vec_isgc <- (dna_vec == 2 | dna_vec == 3)
                       lv <- length(dna_vec)
                       ret_vec <- cbind(length(which(dna_vec_isgc)),
                                        length(which(dna_vec == 1 | dna_vec == 3)),
                                        length(which(dna_vec_isgc[1:(lv-1)] & dna_vec_isgc[2:lv])))
                       ret_vec
                   })
dna_data_frame <- data.frame(do.call(rbind, dna_data))
names(dna_data_frame) <- c("local_GC",
                           "local_purine",
                           "local_CpG")

## for some reason the start and end coordinates are swapped; fix this
osu_data_2 <- osu_data[, c(1,3,2,4:ncol(osu_data))]

## remove tfName and add new tf data
osu_data_3 <- osu_data_2[, setdiff(colnames(osu_data_2), "tfName")]
osu_data_3 <- cbind(osu_data_3, data.frame(tf_snp_lists))

## remove dnaSeq and add the three DNA sequence-derived fields
osu_data_4 <- osu_data_3[, setdiff(colnames(osu_data_3), "dnaSeq")]
osu_data_4 <- cbind(osu_data_4, dna_data_frame)

## read GENCODE data from the hg19 BED file
gencode_TSS_granges <- import.bed("gencode_19_hg19_TSS.bed", genome="hg19")
snps_granges <- GRanges(osu_data_4[, 1:3])

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

## make a BED file so we can calculate 100-bp-averaged GERP scores
snps_plusminus50 <- osu_data_2[, 1:3]
snps_plusminus50[,2] <- snps_plusminus50[,2]-50
snps_plusminus50[,3] <- snps_plusminus50[,2]+50

## need to make a BED file of +/- 50 bp windows around SNPs, for averaging GERP++ data
write.table(data.frame(snps_plusminus50[, 1:3], rownames(osu_data_2)),
            file="osu_snps_plusminus50.bed",
            sep="\t",
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE)

## compute average GERP++ scores within 100 bp windows around SNPs
system("bin/bigWigAverageOverBed  All_hg19_RS.bw osu_snps_plusminus50.bed gerp_100bp.txt")

## read in the 100 bp-averaged GERP++ scores
feature_avg_gerp <- read.table("gerp_100bp.txt",
                               header=FALSE,
                               sep="\t",
                               stringsAsFactors=FALSE)[,5]

all_chroms <- as.character(osu_data_4$chrom)
unique_chroms <- unique(as.character(all_chroms))
nchrom <- length(unique_chroms)
chrom_rows_list <- lapply(unique_chroms, function(chrom_name) {
    which(all_chroms == chrom_name)
})
names(chrom_rows_list) <- unique_chroms

## compute the average minor allele frequency for SNPs within a 1 kb window

feature_mean_allele_frequency_1kb <- apply(osu_data_4[, c("chrom",
                                                          "chromStart",
                                                          "minorAlleleFreq")], 1,
                                function(my_row) {
                                    my_chrom <- as.character(my_row[1])
                                    my_coord <- as.integer(my_row[2])
                                    chrom_rows <- chrom_rows_list[[my_chrom]]
                                    inds_near <- chrom_rows[which(abs(osu_data_4$chromStart[chrom_rows] - my_coord) <= 500)]
                                    mean(osu_data_4$minorAlleleFreq[inds_near])
                                })

#cpg_island_data <- read.table("cpgIslandExt.txt",
#                              header=FALSE,
#                              sep="\t",
#                              stringsAsFactors=FALSE)[,c(2,3,4)]
#names(cpg_island_data) <- c("chrom","chromStart","chromEnd")
#cpg_island_granges <- GRanges(cpg_island_data)
N <- nrow(osu_data_4)
feature_cpg_island <- rep(0,N)
cpg_island_granges <- import.bed("cpgIslandExt.bed")
overlaps_res <- findOverlaps(snps_granges, cpg_island_granges)
feature_cpg_island[queryHits(overlaps_res)] <- 1

## annotate SNPs based on TFBS motifs using JASPAR motif match data from Ensembl BioMart (Ensembl Release 75, archive.ensembl.org)
jaspar_data <- read.table("jaspar_tfbs_ensembl_75_hg19.txt",
                          header=TRUE,
                          sep="\t",
                          stringsAsFactors=FALSE)

jaspar_data_for_granges <- data.frame(jaspar_data[, 2:5],
                                      stringsAsFactors=FALSE)
names(jaspar_data_for_granges) <- c("chrom","chromStart","chromEnd","name")
jaspar_data_for_granges$chrom <- paste("chr", jaspar_data_for_granges$chrom, sep="")
                          
jaspar_granges <- GRanges(jaspar_data_for_granges[,1:3])
overlaps_res <- findOverlaps(snps_granges, jaspar_granges)
feature_jaspar <- rep("", N)
feature_jaspar[queryHits(overlaps_res)] <- jaspar_data_for_granges$name[subjectHits(overlaps_res)]

## make a GRanges object for 1 kb blocks around each of our SNPs
snps_df <- osu_data_4[,1:3]
snps_df_1kb_blocks <- snps_df
snps_df_1kb_blocks$chromStart <- snps_df_1kb_blocks$chromStart - 500
snps_df_1kb_blocks$chromEnd <- snps_df_1kb_blocks$chromEnd + 500
snps_1kb_blocks_granges <- GRanges(snps_df_1kb_blocks)

## compute average heterozygosity within 1 kb region around our SNPs
het_granges <- import.bed("het_rates_filtered.bed")
het_overlap <- findOverlaps(snps_1kb_blocks_granges, het_granges)
het_overlap_df <- as.data.frame(het_overlap)
features_avg_het <- pbsapply(1:N, function(my_ind) {
    mean(het_granges$score[het_overlap_df$subjectHits[het_overlap_df$queryHits == my_ind]])
})
features_avg_het[is.na(features_avg_het)] <- mean(na.omit(features_avg_het))

## compute the average "derived allele frequency" within 1 kb region around our SNPs
daf_granges <- import.bed("daf_filtered.bed")
daf_overlap <- findOverlaps(snps_1kb_blocks_granges, daf_granges)
daf_overlap_df <- as.data.frame(daf_overlap)
features_daf <- pbsapply(1:N, function(my_ind) {
    mean(daf_granges$score[daf_overlap_df$subjectHits[daf_overlap_df$queryHits == my_ind]])
})
features_daf[which(is.na(features_daf))] <- mean(na.omit(features_daf))

## annotate SNPs if they lie within a "simple repeat" using data from UCSC
simple_repeat_bed <- import.bed("simpleRepeat.bed")
simple_repeat_overlap <- findOverlaps(snps_granges, simple_repeat_bed)
feature_simple_repeat <- rep(0, N)
feature_simple_repeat[queryHits(simple_repeat_overlap)] <- 1

## annotate SNPs if they lie within a "nested repeat" using data from UCSC
nested_repeat_bed <- import.bed("nestedRepeat.bed")
nested_repeat_overlap <- findOverlaps(snps_granges, nested_repeat_bed)
feature_nested_repeat <- rep(0, N)
feature_nested_repeat[queryHits(nested_repeat_overlap)] <- 1

## create a SNP annotation equivalent of the "EXON", "UTR5", "UTR3", "INTRON", and "CDS" features in GWAVA, based
## on the "knownGene.txt" file from UCSC
gene_annot_granges <- import.bed("knownGene_filt.bed")
gene_annot_overlap <- findOverlaps(snps_granges, gene_annot_granges)
feature_gene_annot <- rep("INTERGENIC",N)
feature_gene_annot[queryHits(gene_annot_overlap)] <- gene_annot_granges$name[subjectHits(gene_annot_overlap)]

## distance to nearest splice site
known_gene_exons_df <- read.table("knownGene_exons.bed",
                                  header=FALSE,
                                  sep="\t",
                                  stringsAsFactors=FALSE)
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

## augment the OSU feature matrix with the new features, eliminating some feature columns that have been expanded (i.e., tfName and dnaSeq)
osu_data_final <- cbind(osu_data[, setdiff(colnames(osu_data), c("label",
                                                                 "tfName",
                                                                 "dnaSeq",
                                                                 "refUCSC"))],
                        dna_data_frame,
                        data.frame(tf_snp_lists),
                        avg_daf=features_daf,
                        avg_het=features_avg_het,
                        cpg_island=feature_cpg_island,
                        pwm=feature_jaspar,
                        simplerepeat=feature_simple_repeat,
                        nestedrepeat=feature_nested_repeat,
                        gencode_tss=feature_gencode_TSS_dists,
                        maf1kb=feature_mean_allele_frequency_1kb,
                        avg_gerp=feature_avg_gerp,
                        geneannot=feature_gene_annot,
                        ss_dist=feature_splice_site_dists,
                        gwava_data[,c("WEAK_ENH",
                                      "ENH",
                                      "REP",
                                      "TSS_FLANK",
                                      "TRAN",
                                      "TSS",
                                      "CTCF_REG")],
                        label=osu_data$label)

## Yes we are using some GWAVA columns above, but Yao has python code that will give us those columns
## independently of any GWAVA code 

## put chromStart column before chromEnd column
osu_data_final <- osu_data_final[, c(1,3,2,4:ncol(osu_data_final))]

## write the feature matrix to a file
write.table(osu_data_final,
            file="osudata_recomb.txt",
            sep="\t",
            row.names=TRUE,
            col.names=TRUE,
            quote=FALSE)

            
