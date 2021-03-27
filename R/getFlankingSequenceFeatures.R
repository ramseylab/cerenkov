# Adds DNAShape + GC Content features to a list of SNPs

# Command Line Arguments
args = commandArgs(trailingOnly=TRUE)
helpMsg <- "Example Invocation: Rscript --vanilla getFlankingSequenceFeatures.R csnps.ucsc_metadata.txt csnps.SequenceFeatures.txt"
inputFName=args[1]
outputFName=args[2]
if (length(args) < 2) {
  stop(helpMsg, call.=FALSE)
} 

# Init
set.seed(1337)
GC_WINDOW = 5 # Size of window to calculate GC_Content
DNASHAPE_WINDOW = 11 # Size of window to calculate GC_Content

# Bioconductor ref: http://www.gettinggeneticsdone.com/2011/04/using-rstats-bioconductor-to-get.html
library('BSgenome.Hsapiens.UCSC.hg19')
library(stringr)


# Expects TSV file with UCSC-format columns: c("name", "chrom", "chromStart", "chromEnd")
# snps <- read.table(inputFName, sep='\t', header=T)

# Expect a BED format
snps <- read.table(inputFName, sep='\t', header=FALSE, col.names=c("chrom", "chromStart", "chromEnd", "name"))
snps <- snps[, c("name", "chrom", "chromStart", "chromEnd")]
snps$chrom <- as.character(snps$chrom)


offset <- (DNASHAPE_WINDOW - 1)/2
snps$dnaSeq <- apply(snps, 1, function(x) {
	chrom <- x['chrom']
	chromStart <- as.integer(x['chromStart']) - offset
	chromEnd <- as.integer(x['chromEnd']) + offset - 1 # UCSC coordinate conventions
	as.character(getSeq(Hsapiens, chrom, chromStart, chromEnd))
})


getGCContent <- function(dnaSeqVec, window=GC_WINDOW) {
	flank <- (window - 1)/2
	seqLength <- nchar(dnaSeqVec[1])
	cutStart <- (seqLength - 1)/2 - flank
	cutStop <- (seqLength - 1)/2 + flank
	GCs <- sapply(dnaSeqVec, FUN=substr, start=cutStart, stop=cutStop)
	return( str_count(GCs, 'G') + str_count(GCs, 'C') )
} 

snps$GC5Content <- getGCContent(snps$dnaSeq, window=5)
snps$GC7Content <- getGCContent(snps$dnaSeq, window=7)

# devtools::install_github("ramseylab/regshape")
library(regshape)
snps$DNAShapeScore <- sapply(snps$dnaSeq, getShapeScores)

# Finish
snps <- snps[, c("name", "dnaSeq", "GC5Content", "GC7Content", "DNAShapeScore")]
write.table(snps, file=outputFName, sep='\t', row.names=FALSE, quote=FALSE)
