# CERENKOV
## 1. Computational Elucidation of the REgulatory NonKOding Variome

CERENKOV is a software pipeline and associated machine-learning framework for
identifying regulatory single nucleotide polymorphisms (rSNPs) in the noncoding
genome for post-analysis of genetic regions identified in genome-wide
association studies (GWAS). CERENKOV was created by Yao Yao, Zheng Liu,
Satpreet Singh, Qi Wei, and Stephen Ramsey at Oregon State University. 

The March 2017 data files for CERENKOV can be accessed on the
[Ramsey Lab file server](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201703)
(see README.md files under the subdirectories of the GitHub CERENKOV project
area, for more information about which data files are used in which parts of
CERENKOV).

## 2. Reproducing the results of the CERENKOV2 

Based on the 2017 ACM-BCB version, we further developed CERENKOV2 and submitted a methodology article 
_CERENKOV2: data-space geometric features improve machine learning-based detection of functional noncoding SNPs_
to BMC Bioinformatics in July 2018. 

To reproduce the results reported in this submission, please follow the [README file in _experiments/CERENKOV2_ folder](https://github.com/ramseylab/cerenkov/blob/master/experiments/CERENKOV2/README.md).

We revised our manuscript and code in the following weeks and re-submitted with a new title 
_CERENKOV2: improved detection of functional noncoding SNPs using data-space geometric features_ in Nov. 2018.

To reproduce the results reported in revised submission, please follow the [README file in _experiments/CERENKOV2\_revision_ folder](https://github.com/ramseylab/cerenkov/blob/master/experiments/CERENKOV2_revision/README.md).

We also provided an R scripts, [`install_packages.R`](https://github.com/ramseylab/cerenkov/blob/master/experiments/install_packages.R), 
to install all dependencies of reproduction.  

## 3. Reproducing the results of the 2017 ACM-BCB article

We presented the very first version of CERENKOV at the 2017 ACM-BCB conference in Boston in
August 2017, with an accompanying full research article 
_CERENKOV: Computational Elucidation of the Regulatory Noncoding Variome_ 
in the proceedings, describing CERENKOV and demonstrating its accuracy for discriminating rSNPs from
nonfunctional SNPs. 

The corresponding code is archived in release 
[v0.1-alpha](https://github.com/ramseylab/cerenkov/releases/tag/v0.1-alpha).

To reproduce the results of our article, please follow the instructions below.

### 3.1. Guide to source code files in CERENKOV v0.1-alpha:

- `cerenkov_ml_compare_models.R`: obtains the comparative machine-learning
performance results that were used to make Figure 3 of the article (the script
`cerenkov_analyze_ml_results_compare_models.R` actually generates the plot).

The `cerenkov_ml_compare_models.R` script will require the R packages `PRROC`, `parallel`, 
`xgboost` (version 0.6-4), `ranger` (version 0.6.0), `Matrix`, and `pbapply`. The
`cerenkov_analyze_ml_results_compare_models.R` script will require the R packages `ggplot2`
and `reshape2`.

- `cerenkov_ml_xgboost_importance.R`: obtains the feature importance scores for
the CERENKOV method that were used to make Figure 2 in the article. 

The script requires the R packages `Matrix`, `parallel`, `pbapply`, and `xgboost`.

- `cerenkov_ml_tune_xgboost.R`: obtains the grid-search tuning machine-learning
performance results that were used to make Figure 1b,c in the article (the
script `cerenkov_analyze_ml_results_tune_xgboost.R` actually generates the
plot).

The `cerenkov_ml_tune_xgboost.R` script will require the R packages `PRROC`,
`parallel`, `xgboost` (version 0.6-4), `Matrix`, and `pbapply`. The
`cerenkov_analyze_ml_results_tune_xgboost.R` script will require the R packages
`ggplot2` and `reshape2`.


### 3.2 Data files for download

The following .Rdata files (the "201703 data supplement" for CERENKOV) accompany
the article *CERENKOV: Computational elucidation of the regulatory noncoding
variome* by Yao Yao, Zheng Liu, Satpreet Singh, Qi Wei, and Stephen A. Ramsey
(submitted to the ACM-BCB conference, April 2017). The data files are available
at the following links (all are HTTP links to the file server `files.cgrb.oregonstate.edu`):

- [`features_GWAVA.Rdata`](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201703/features_GWAVA.Rdata) (604 kB)
- [`features_OSU.Rdata`](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201703/features_OSU.Rdata) (3.2 MB)
- [`features_RSVP.Rdata`](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201703/features_RSVP.Rdata) (77 MB)
- [`features_danq.Rdata`](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201703/features_danq.Rdata) (206 MB)
- [`features_deepsea.Rdata`](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201703/features_deepsea.Rdata) (176 MB)
- [`features_deltaSVM.Rdata`](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201703/features_deltaSVM.Rdata) (22 MB)
- [`scores_cadd.Rdata`](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201703/scores_cadd.Rdata) (132 kB)
- [`scores_dann.Rdata`](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201703/scores_dann.Rdata) (151 kB)
- [`scores_eigen.Rdata`](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201703/scores_eigen.Rdata) (152 kB)
- [`scores_fitcons.Rdata`](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201703/scores_fitcons.Rdata) (54 kB)
- [`snp_coordinages.Rdata`](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201703/snp_coordinates.Rdata) (104 kB)

