# Get ensGene from UCSC/Local: mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A hg19 -e "SELECT * FROM ensGene" > ensGene.txt
# Run in workspace/

## LOGIC (per discussion with Steve 2016-08-11)
#  Comparing SNPs from snp146 with ensGene
#  If in a gene (between txStart and txEnd), exclude
#  If NOT in an exon (i.e. not between exon start-end locations), KEEP
#  If in an exon, check if gene is protein coding (i.e. cdsStart != cdsEnd) 
#    if not, eliminate. If is protein-coding, then check if snp-coordinate 
#    is in [cdsStart, cdsEnd]. 
#      If is NOT in above [...], then KEEP, else eliminate

# Below code tries to find uneligible 


# R
rm(list = ls())

## STEP 1: Flatten out the exonStart-exonEnd pairs 
ensGene <- read.csv("ensGene.txt", sep='\t', stringsAsFactors=FALSE)
ensGene <- ensGene[, c('name', 'chrom', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds')]

totalExons <- sum(ensGene[, 'exonCount'])
exonData <- matrix('', nrow=totalExons, ncol=3)
eIdx <- 1
eRows <- nrow(ensGene)
for (n in 1:eRows) {
  if (n %% 1000 == 0) { print(sprintf("Processing row %d of %d", n, eRows)) }

  exonCount <- ensGene[n, 'exonCount']
  txName <- ensGene[n, 'name']
  exonStarts <-  unlist(strsplit(ensGene[n, 'exonStarts'],','))
  exonEnds <-  unlist(strsplit(ensGene[n, 'exonEnds'],','))

  for (k in 1:exonCount) {
    exonData[eIdx, ] <- c(txName, exonStarts[k], exonEnds[k])
    eIdx = eIdx + 1
  }
}

exonData <- as.data.frame(exonData)
names(exonData) <- c("name", "exStart", "exEnd")

ensGene <- ensGene[, c('name', 'chrom', 'cdsStart', 'cdsEnd', 'exonCount')]
ensGeneExpanded <- merge(ensGene, exonData)

# Verify 
print(totalExons)
print(nrow(ensGeneExpanded))

## STEP 2: Join with rSNPs per Steve's logic 
rsnps <- read.csv("r_all.txt", sep='\t')
csnps <- read.csv("c_all.txt", sep='\t')
rsnps$label <- 1
csnps$label <- 0
snps <- rbind(rsnps, csnps)

snps <- snps[, c('name', 'chrom', 'chromStart', 'chromEnd', 'label')]

library(sqldf)

eliminations <- sqldf("SELECT s.*, e.* FROM
  ensGeneExpanded e JOIN snps s
  ON 
    s.chrom = e.chrom
    -- is in an Exon
    AND s.chromStart BETWEEN e.exStart AND (e.exEnd-1)
    -- gene is NOT protein-coding
    AND e.cdsStart != e.cdsEnd 
    AND s.chromStart BETWEEN e.cdsStart AND (e.cdsEnd-1)
  ")


# eliminations <- sqldf("SELECT s.*, e.*, (e.cdsStart != e.cdsEnd) AS exonProCoding FROM
#   ensGeneExpanded e JOIN snps s
#   ON 
#     s.chrom = e.chrom
#     -- is in an Exon
#     AND s.chromStart BETWEEN e.exStart AND (e.exEnd-1)
#     AND (
#       -- gene is NOT protein-coding
#       (e.cdsStart = e.cdsEnd) 
#       OR
#       -- in a protein-coding gene cds
#       (
#         (e.cdsStart != e.cdsEnd) 
#          AND 
#         (s.chromStart BETWEEN e.cdsStart AND (e.cdsEnd-1))
#       )
#     )
#   ")

names(eliminations)[1] <- "snpID"
names(eliminations)[2] <- "snpChrom"
write.csv(eliminations, file="snp_noncoding_violations.csv")

## ANALYSIS
print(sqldf("SELECT label, COUNT(1) numSNPs FROM eliminations GROUP BY label"))

length(eliminations$snpID) 
length(unique(eliminations$snpID)) 

uniqEliminationsByLabel <- unique(eliminations[, c('snpID', 'label')])
nrow(uniqEliminationsByLabel) # 2486
print(sqldf("SELECT label, COUNT(1) numSNPs FROM uniqEliminationsByLabel GROUP BY label"))

dups <- sqldf("SELECT snpID, COUNT(1) repeats FROM eliminations GROUP BY 1 HAVING repeats > 1")
dupExamples <- subset(eliminations, snpID %in% dups$snpID)
dupExamples <- dupExamples[with(dupExamples, order(snpID)), ]
write.csv(dupExamples, file="duplicatesExamples.csv", row.names=F)


## PCE Feature: Output for Pipeline
pce <- data.frame(name=unique(eliminations$snpID), proteinCodingExon=1)
write.table(pce, file="pce_eliminations.txt", sep='\t', row.names=F, quote=F)
