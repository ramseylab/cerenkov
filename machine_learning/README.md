# CERENKOV R scripts #

To reproduce the results of our article *CERENKOV:* *Computational* 
*Elucidation* *of* *the* *Regulatory* *Noncoding* *Variome* 
the following R scripts can be run in the R statistical computing
environment (version 3.2.1):

- `cerenkov_ml_compare_models.R`: obtains the comparative machine-learning
performance results that were used to make Figure 3 of the article (the script
`cerenkov_analyze_ml_results_compare_models.R` actually generates the plot).

The `cerenkov_ml_compare_models.R` script will require the R packages PRROC, parallel, 
xgboost (version 0.6-4), ranger (version 0.6.0), Matrix, and pbapply. The
`cerenkov_analyze_ml_results_compare_models.R` script will require the R packages ggplot2
and reshape2.

- `cerenkov_ml_xgboost_importance.R`: obtains the feature importance scores for
the CERENKOV method that were used to make Figure 2 in the article. 

The script requires the R packages Matrix, parallel, pbapply, and xgboost.

- `cerenkov_ml_tune_xgboost.R`: obtains the grid-search tuning machine-learning
performance results that were used to make Figure 1b,c in the article (the
script `cerenkov_analyze_ml_results_tune_xgboost.R` actually generates the
plot).

The `cerenkov_ml_tune_xgboost.R` script will require the R packages PRROC,
parallel, xgboost (version 0.6-4), Matrix, and pbapply. The
`cerenkov_analyze_ml_results_tune_xgboost.R` script will require the R packages
ggplot2 and reshape2.


# Data files for download

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
