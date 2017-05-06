# recomb-rsnp-data

## Uniform DHS Scores

Please copy the `UCSC_DHS` folder at [afp://netshare@sramsey-macmini.local/netshare1/Users/yao.yao](afp://netshare@sramsey-macmini.local/netshare1/Users/yao.yao) to `datadir`.

## DNAShape + GC Content

Please install libs below before running the scripts.

```r
# Bioconductor ref: http://www.gettinggeneticsdone.com/2011/04/using-rstats-bioconductor-to-get.html
source("http://www.bioconductor.org/biocLite.R")
biocLite("BSgenome")
biocLite("BSgenome.Hsapiens.UCSC.hg19") #installs the human genome (~850 MB download).

# DNAShape lib
devtools::install_github("ramseylab/regshape")
```

## eQTL p-values

Please copy the `GTEx_Analysis_V4_eQTLs` folder at [afp://netshare@sramsey-macmini.local/netshare1/Users/yao.yao](afp://netshare@sramsey-macmini.local/netshare1/Users/yao.yao) to `datadir`.

Please install libs below before running the scripts.

```r
install.packages("plyr")
```