# CERENKOV R scripts #


## 1. Guide to current CERENKOV machine-learning R code

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



## 3. Reproducing the results of the 2018 BMC-Bioinformatics submission 

To reproduce the results of our 2018 BMC-Bioinformatics submission, CERENKOV2 with SNP-radius-LLR features, the following R
scripts from release **v2.0-alpha** of the project can be run in the R statistical computing environment (version 3.4.4).

### 3.1. Guide to source code files in CERENKOV v2.0-alpha

- `cerenkov_script_run_ml_llr.R`: Replicate 5-fold xgboost training and validation 10 times, for both CERENKOV and CERENKOV2, on OSU18 dataset.
- `cerenkov_script_run_ml_llr_rf.R`: Learn a Random Forest model on all OSU18 data simply for feature importance scores.
- `cerenkov_script_analyze_ml_llr.R`: obtains the comparative machine-learning
performance results that were used to make Figure 4 and Table 2 of the article.

### 3.2 Data files for download

- [`osu18_features1.1_cerenkov2.rds`](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201806/osu18_features1.1_cerenkov2.rds) (8.5 MB)
- [`osu18_intra_locus_dist.rds`](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201806/osu18_intra_locus_dist.rds) (2.2 MB)
- [`osu18_snp_coordinates.rds`](http://files.cgrb.oregonstate.edu/Ramsey_Lab/cerenkov/datafiles_201806/osu18_snp_coordinates.rds) (167 kB)
