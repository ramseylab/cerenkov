# Command Line Arguments
args = commandArgs(trailingOnly=TRUE)
helpMsg <- "Example Invocation: Rscript --vanilla getFlankingSequenceFeatures.R ../datadir/GTEx_Analysis_V4_eQTLs csnps.ucsc_metadata.txt csnps.SequenceFeatures.txt"

if (length(args) < 3) {
	stop(helpMsg, call.=FALSE)
} 

eQTLFolder <- args[1]
inputFName <- args[2]
outputFName <- args[3]

library(plyr)

##### Utils #####

read.snp.name <- function(metadata.filename) {
	snp.name <- read.table(metadata.filename, colClasses = c("character", rep("NULL", 18)), 
						header = TRUE, stringsAsFactors = FALSE, sep = "\t")
	
	# Check duplication
	# print(sum(!duplicated(csnps.ucsc.name$name)))
	
	# subsetting one column would degrade you data frame to vector; use `drop=FALSE` to prevent this.
	snp.name <- snp.name[!duplicated(snp.name), , drop=FALSE]
	
	return(snp.name)
}

read.eqtl <- function(directory) {
	merge.multi.eqtl <- function (eqtl.dir) {
		eqtl.files <- list.files(eqtl.dir, pattern = "\\.eqtl$", full.names = TRUE)
		pvalue.table <- lapply(eqtl.files, function(x) { 
			read.table(x, colClasses = c("character", rep("NULL", 7), "numeric", rep("NULL", 3)), 
					   header = TRUE, stringsAsFactors = FALSE, sep = "\t")
		}) 
		do.call(rbind, pvalue.table)
	}
	
	pvalue.table <- merge.multi.eqtl(directory)
	colnames(pvalue.table)[1] <- "name"
	
	return(pvalue.table)
}

lookup <- function(snp.name, pvalue.table) {
	# Merge by column "name"
	snp.pvalue <- merge(snp.name, pvalue.table, by.x = "name", by.y = "name", all.x = TRUE, all.y = FALSE)
	
	return(snp.pvalue)
}

dedup <- function(snp.pvalue) {
	# Eliminate duplicated entry; Keep the smallest
	snp.pvalue.dedup <- function (np.pair, keep.smallest=TRUE) {
		if (keep.smallest == TRUE) {
			np.pair <- np.pair[order(np.pair$name, np.pair$P_Val), ] #sort by `name` and `P_Val`
		} else {
			np.pair <- np.pair[order(np.pair$name, -np.pair$P_Val), ] #sort by `name` and `-P_Val`
		}
		
		np.pair[!duplicated(np.pair$name), ]   
	}
	
	snp.pvalue <- snp.pvalue.dedup(snp.pvalue)
	
	return(snp.pvalue)
}

recalibrate <- function(snp.pvalue) {
	# P_Val Recalculation
	snp.pvalue <- adply(snp.pvalue, .margins = 1, .fun = transform, P_Val = if (is.na(P_Val)) 0 else -log10(P_Val))
	
	return(snp.pvalue)
}

write.snp.pvalue <- function(snp.pvalue, filename) {
	write.table(snp.pvalue, col.names=c("name", "eqtlPvalue"), 
				file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
}

##### High-level function #####

extract <- function(snp.name, pvalue.table) {
	snp.pvalue <- lookup(snp.name, pvalue.table)
	snp.pvalue <- dedup(snp.pvalue)
	snp.pvalue <- recalibrate(snp.pvalue)
	
	return(snp.pvalue)
}

##### main() #####

`%str+%` = function(x, y) {
	if(is.character(x) || is.character(y)) {
		return(paste(x , y, sep=""))
	} else {
		.Primitive("+")(x,y)
	}
}
	
pvalue.table <- read.eqtl(eQTLFolder)

snp.name <- read.snp.name(inputFName)

write.snp.pvalue(extract(snp.name, pvalue.table), outputFName)
