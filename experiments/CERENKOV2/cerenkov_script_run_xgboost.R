## ================================================================================
## This file is part of the CERENKOV (Computational Elucidation of the
## REgulatory NonKOding Variome) software program. CERENKOV is subject to terms
## and conditions defined in the file "LICENSE.txt", which is part of this
## CERENKOV software distribution.
##
## Title       : cerenkov_script_run_xgboost.R
## 
## Description : Top-level R script for running a CERENKOV
##               machine-learning job with LLR feature reducer.
##
##               Uses:  cerenkov_func_base.R
##                      cerenkov_func_aws.R
##                      cerenkov_incscript_setup_aws.R
##                      cerenkov_incscript_setup_ml.R
##                      cerenkov_incscript_run_ml.R
##
##
## Requires    : PRROC
##
## May require : xgboost, ranger, dismo, Matrix, aws.ec2, aws.signature, pbapply, methods, Rcpp, dplyr
##
## Note        : do not use this program with Ranger version 0.6.6 (stability issues);
##               use Ranger 0.6.0 only;
## 				 AWS configuration removed; run locally only;
##
## 
## Author      : Yao Yao, Oregon State University
##               https://github.com/erikyao
## ================================================================================

g_args <- commandArgs(trailingOnly=FALSE)  ## convert to using "aargh" package, some day

source("../../machine_learning/cerenkov_func_base.R")         ## load functions used for machine-learning
source("../../machine_learning/cerenkov_func_aws.R")          ## load functions used for AWS

g_trailing_args <- g_get_args_after_hyphen_hyphen_args(g_args)

## ============================ define global parameters =============================

g_par <- c(
    list(
        num_folds_cross_validation =      5,     ## we are standardizing on 5-fold CV
        num_cv_replications =             10,            ## set to anywhere from 1--200, typically
        flag_create_fork_cluster =        FALSE,   ## TRUE or FALSE
        override_num_fork_processes =     8,   ## for EC2, set to 64; for my MBP, set to 8
        show_progress_bar =               TRUE,
        parallel_use_load_balancing =     FALSE,
        flag_locus_sampling =             TRUE,         ## set to false if you want SNP-level sampling
        # flag_xgb_importance =             FALSE,        ## if you set this to TRUE, make sure you set num_cv_replications=1
                                                        ## and num_folds_cross_validation=1
        random_number_seed =              if (is.na(g_trailing_args[1]) ||
                                              is.na(suppressWarnings(as.integer(g_trailing_args[1]))))
                                              1337 else as.integer(g_trailing_args[1]),
        nthreads_per_process =            1,
        flag_randomize_classifier_order = FALSE,
        flag_create_ec2_multi_cluster =   FALSE,  ## don't set this to true if you set "flag_create_fork_cluster" to true
        analysis_label =                  "CERENKOV2",
        output_file_base_name =           "results_xgboost",
        debug_file_parallel =             "",
        debug=                            TRUE
    )
)

## ============================== load OSU feature data; check for problems with feature data =================================
print("loading OSU data")
g_annot_feat_df <- readRDS(file="osu18_features1.1_cerenkov2.rds")
stopifnot(g_feature_matrix_is_OK(g_annot_feat_df))

## workaround for bug in Matrix (
## see Stack Overflow: https://stackoverflow.com/questions/43282720/r-error-in-validobject-object-when-running-as-script-but-not-in-console)
## note: Rscript doesn't load methods by default, but R interactive does!
library(methods)  

g_annot_feat_sm <- Matrix::sparse.model.matrix(label ~ .-1, data=g_annot_feat_df)

g_snp_names <- rownames(g_annot_feat_df)
g_label_vec <- as.integer(as.character(g_annot_feat_df$label))

## ============================== load feature data; check for problems with feature data =================================

print("loading intra-locus distance data")
g_geom_feat_df <- readRDS("osu18_intra_locus_dist.rds")
stopifnot(rownames(g_geom_feat_df) == rownames(g_annot_feat_df))


## ============================== build a list of the feature matrices that we will need =================================
g_classifier_feature_matrices_list <- list(
    annot_feat=g_annot_feat_sm,
    geom_feat=g_geom_feat_df
)

## ============================== run invariant setup code =================================
source("../../machine_learning/cerenkov_incscript_setup_ml.R")


## ============================== make closures for classifiers  =================================

g_classifier_function_xgb = g_make_classifier_function_xgboost(p_nthread=g_par$nthreads_per_process,
															   p_get_perf_results=g_get_perf_results,
															   p_feature_importance_type=NULL,
															   p_make_objective_function=function(...){"binary:logistic"},
															   p_case_group_ids=g_snp_locus_ids)

g_classifier_functions_list <- list(
    xgb=g_classifier_function_xgb
)

## ============================== assemble final list of feature matrices  =================================

## free up memory
rm(g_annot_feat_df)
rm(g_annot_feat_sm)
rm(g_geom_feat_df)


## ============================== make hyperparameter lists  =================================

## ------------------ xgboost hyperparameter lists --------------

g_xgb_hp_grid <- g_make_hyperparameter_grid_list(list(eta=c(0.1),
                                                      nrounds=c(40),
                                                      gamma=c(10),
                                                      lambda=c(1),
                                                      subsample=c(1),
                                                      # colsample_bytree=c(0.75, 0.85, 1.0),
                                                      base_score=g_class_count_frac_positive,
                                                      scale_pos_weight=c(1.0),
                                                      max_depth=c(7)))


g_classifier_list_xgb_control <- lapply(g_xgb_hp_grid,
	                                    function(hp_cell) {
	                                        list(classifier_feature_matrix_name="annot_feat",
	                                             classifier_function_name="xgb",
	                                             classifier_hyperparameter_set_type_name="g_xgb_hp_grid",
	                                             classifier_set_name="CK", # Control model: CERENKOV annotation features only
	                                             classifier_hyperparameter_list=hp_cell)
	                                    })

g_classifier_list_xgb_case <- lapply(g_xgb_hp_grid,
									 function(hp_cell) {
									 	list(classifier_feature_matrix_name="annot_feat",
									 		 classifier_function_name="xgb",
									 		 classifier_hyperparameter_set_type_name="g_xgb_hp_grid",
									 		 classifier_set_name="CK+LLR", # Case model: CERENKOV annotation features + geometric LLR features
									 		 classifier_hyperparameter_list=hp_cell, 
									 		 feature_reducer_function_name="FIT_LLR", 
									 		 feature_reducer_input_matrix_name="geom_feat", 
									 		 feature_reducer_hyperparameters_list=list())
									 })


## ============================== assemble classifier list  =================================

## TODO: make a "g_check_classifier_list" function that checks for incorrect classifier function
## names, incorrect feature matrix names, etc.

g_classifier_list <- c(
	g_classifier_list_xgb_control,
    g_classifier_list_xgb_case
)

## ============================== run invariant machine-learning code =================================

                              
source("../../machine_learning/cerenkov_incscript_run_ml.R")
