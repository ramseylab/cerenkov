# CERENKOV feature extraction

## Prerequisites

### DNAShape + GC Content

Please install libs below before running the scripts.

```r
# Bioconductor ref: http://www.gettinggeneticsdone.com/2011/04/using-rstats-bioconductor-to-get.html
source("http://www.bioconductor.org/biocLite.R")
biocLite("BSgenome")
biocLite("BSgenome.Hsapiens.UCSC.hg19") #installs the human genome (~850 MB download).

# DNAShape lib
devtools::install_github("ramseylab/regshape")
```

### CADD

Please donwload [1000G_phase3.tsv.gz](http://krishna.gs.washington.edu/download/CADD/v1.3/1000G_phase3.tsv.gz).

### fitCons

Please download [fc-gm-0.bw](http://compgen.cshl.edu/fitCons/0downloads/tracks/V1.01/gm/scores/fc-gm-0.bw), 
[fc-h1-0.bw](http://compgen.cshl.edu/fitCons/0downloads/tracks/V1.01/h1/scores/fc-h1-0.bw), 
[fc-h1-0.bw](http://compgen.cshl.edu/fitCons/0downloads/tracks/V1.01/hu/scores/fc-hu-0.bw) and 
[fc-i6-0.bw](http://compgen.cshl.edu/fitCons/0downloads/tracks/V1.01/i6/scores/fc-i6-0.bw)

### eQTL p-values

Please install libs below before running the scripts.

```r
install.packages("plyr")
```

### FSU Repli-chip

Please download all bigwig files on [http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeFsuRepliChip/](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeFsuRepliChip/).

### UW Repli-chip

Please download all the bigwig files listed below.

- [wgEncodeUwRepliSeqBg02esWaveSignalRep1.bigWig](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBg02esWaveSignalRep1.bigWig)
- [wgEncodeUwRepliSeqBjWaveSignalRep1.bigWig](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBjWaveSignalRep1.bigWig)
- [wgEncodeUwRepliSeqBjWaveSignalRep2.bigWig](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBjWaveSignalRep2.bigWig)
- [wgEncodeUwRepliSeqGm06990WaveSignalRep1.bigWig](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm06990WaveSignalRep1.bigWig)
- [wgEncodeUwRepliSeqGm12801WaveSignalRep1.bigWig](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm12801WaveSignalRep1.bigWig)
- [wgEncodeUwRepliSeqGm12812WaveSignalRep1.bigWig](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm12812WaveSignalRep1.bigWig)
- [wgEncodeUwRepliSeqGm12813WaveSignalRep1.bigWig](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm12813WaveSignalRep1.bigWig)
- [wgEncodeUwRepliSeqGm12878WaveSignalRep1.bigWig](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm12878WaveSignalRep1.bigWig)
- [wgEncodeUwRepliSeqHelas3WaveSignalRep1.bigWig](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHelas3WaveSignalRep1.bigWig)
- [wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig)
- [wgEncodeUwRepliSeqHuvecWaveSignalRep1.bigWig](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHuvecWaveSignalRep1.bigWig)
- [wgEncodeUwRepliSeqImr90WaveSignalRep1.bigWig](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqImr90WaveSignalRep1.bigWig)
- [wgEncodeUwRepliSeqK562WaveSignalRep1.bigWig](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562WaveSignalRep1.bigWig)
- [wgEncodeUwRepliSeqMcf7WaveSignalRep1.bigWig](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqMcf7WaveSignalRep1.bigWig)
- [wgEncodeUwRepliSeqNhekWaveSignalRep1.bigWig](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqNhekWaveSignalRep1.bigWig)
- [wgEncodeUwRepliSeqSknshWaveSignalRep1.bigWig](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqSknshWaveSignalRep1.bigWig)

### GERP

Please download [All_hg19_RS.bw](http://hgdownload.cse.ucsc.edu/gbdb/hg19/bbi/All_hg19_RS.bw)

### TFBS Summary

Please download [http://ngs.sanger.ac.uk/production/ensembl/regulation/hg19/overview/all_tfbs.bw](all_tfbs.bw)