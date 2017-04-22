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

