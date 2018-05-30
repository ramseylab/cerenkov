# CERENKOV R scripts #


# Guide to current CERENKOV machine-learning R code:

- `cerenkov_func_base.R`: Defines global functions for the CERENKOV project.

- `cerenkov_func_aws.R`: Define global functions related to AWS, for the
  CERENKOV project.

- `cerenkov_incscript_setup_ml.R`: Code that is sourced by a
`cerenkov_script_XXXXXX.R` script in order to setup the machine-learning job.  I
put code in this file if it would otherwise be boilerplate in *every* script of
the form `cerenkov_script_XXXXXX.R`.
          
- `cerenkov_incscript_setup_aws.R`: Code that is sourced by a
`cerenkov_script_XXXXXX.R` script in order to setup an EC2 instance for running
the machine-learning job.  I put code in this file if it would otherwise be
boilerplate in *every* script of the form `cerenkov_script_XXXXXX.R`.

- `cerenkov_incscript_run_ml.R`: Code that is sourced by a `cerenkov_script_XXXXXX.R` script in
order to run the machine-learning job.

- `cerenkov_script_run_ml_template.R`: A template script for running a CERENKOV machine-learning
analysis.

- `cerenkov_script_analyze_ml_results_template.R`: A template script for analyzing, reducing, and
plotting the results of a CERENKOV machine-learning run.

The functions defined by the `cerenkov_func_` scripts all satisfy fhe following four properties:
1. it should not not access the global environment
2. it should be side effect-free, as much as possible
3. function names should start with "g_" (for "global")
4. quietly checks for required packages at run-time

The R scripts of the form `cerenkov_incscript_` are run by being called by `source` from
a top-level `cerenkov_script_` script.

# Reproducing the results of the 2017 ACM-BCB article:

To reproduce the results of our article *CERENKOV:* *Computational*
*Elucidation* *of* *the* *Regulatory* *Noncoding* *Variome* the following R
scripts from release **v0.1-alpha** of the project can be run in the R
statistical computing environment (version 3.2.1).

## Guide to source code files in CERENKOV v0.1-alpha:

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


## Data files for download

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
