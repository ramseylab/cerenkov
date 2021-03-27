# CERENKOV feature extraction

## Prerequisites

### BED utilities

Please install [bedtools - 2.25.0](http://bedtools.readthedocs.io/en/latest/content/installation.html) and [pybedtools - 0.7.8](https://daler.github.io/pybedtools/main.html#installing-pybedtools).

### Bigwig utilities

Please install [kentUtils](https://github.com/ENCODE-DCC/kentUtils).

### MySQL access

Please install [PyMySQL - 0.7.6](https://github.com/PyMySQL/PyMySQL) and [SQLAlchemy - 1.0.14](https://www.sqlalchemy.org/).

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

Please copy the [CADD](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201703/CADD/) folder to **_source_data_**.

### eQTL p-values

Please install libs below before running the scripts.

```r
install.packages("plyr")
```

Please copy the [GTEx_Analysis_V4_eQTLs](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201703/GTEx_Analysis_V4_eQTLs) folder to **_source_data_**.

### fitCons

Please copy the [fitCons](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201703/fitCons/) folder to **_source_data_**.

### FSU Repli-chip

Please copy the [FSU_Repli_chip](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201703/FSU_Repli_chip/) folder to **_source_data_**.

### GERP

Please copy the [GERP](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201703/GERP/) folder to **_source_data_**.

### eQTL

Please copy the [GTEx_Analysis_V4_eQTLs](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201703/GTEx_Analysis_V4_eQTLs/) folder to **_source_data_**.

### TFBS Summary

Please copy the [Sanger_TFBS_Summary](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201703/Sanger_TFBS_Summary/) folder to **_source_data_**.

### Uniform DHS

Please copy the [UCSC_DHS](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201703/UCSC_DHS/) folder to **_source_data_**.

### UW Repli-chip

Please copy the [UW_Repli_Seq](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201703/Sanger_TFBS_Summary/) folder to **_source_data_**.

### Augmented OSU Features

Please copy the [augment_osu_features_datafiles](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201703/augment_osu_features_datafiles/) folder to **_source_data_**.