## ================================================================================
## This file is part of the CERENKOV (Computational Elucidation of the
## REgulatory NonKOding Variome) software program. CERENKOV is subject to terms
## and conditions defined in the file "LICENSE.txt", which is part of this
## CERENKOV software distribution.
##
## Title       : cerenkov_script_run_ranger.R
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
##               This code assumes "cerenkov_func_base.R" and
##               "cerenkov_func_aws.R" have already been sourced.
##
## Requires    : PRROC
##
## May require : xgboost, ranger, dismo, Matrix, aws.ec2, aws.signature, pbapply, methods, Rcpp, dplyr
##
## Note        : do not use this program with Ranger version 0.6.6 (stability issues);
##               use Ranger 0.6.0 only; Ranger 0.6.0 did not support "impurity_corrected" feature importance
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
        num_folds_cross_validation =      1,     ## we are standardizing on 5-fold CV
        num_cv_replications =             1,            ## set to anywhere from 1--200, typically
        flag_create_fork_cluster =        FALSE,   ## TRUE or FALSE
        override_num_fork_processes =     NULL,   ## for EC2, set to 64; for my MBP, set to 8
        show_progress_bar =               TRUE,
        parallel_use_load_balancing =     FALSE,
        flag_locus_sampling =             TRUE,         ## set to false if you want SNP-level sampling
		snp_coord_rds =                   "osu18_snp_coordinates(anonymous).rds",
        # flag_xgb_importance =             FALSE,        ## if you set this to TRUE, make sure you set num_cv_replications=1
                                                        ## and num_folds_cross_validation=1
        random_number_seed =              if (is.na(g_trailing_args[1]) ||
                                              is.na(suppressWarnings(as.integer(g_trailing_args[1]))))
                                              1337 else as.integer(g_trailing_args[1]),
        nthreads_per_process =            1,
        flag_randomize_classifier_order = FALSE,
        flag_create_ec2_multi_cluster =   FALSE,  ## don't set this to true if you set "flag_create_fork_cluster" to true
        analysis_label =                  "CERENKOV2",
        output_file_base_name =           "records_of_feature_importance",
        debug_file_parallel =             "",
        debug=                            TRUE
    ))


## ============================== load OSU feature data; check for problems with feature data =================================
print("loading OSU data")
g_annot_feat_df <- readRDS(file="osu18_features1.1_cerenkov2(anonymous).rds")
stopifnot(g_feature_matrix_is_OK(g_annot_feat_df))

g_snp_names <- rownames(g_annot_feat_df)
g_label_vec <- as.integer(as.character(g_annot_feat_df$label))

# Remove "label" column from dataframe
g_annot_feat_df <- g_annot_feat_df[, which(names(g_annot_feat_df) != "label")]

## ============================== load feature data; check for problems with feature data =================================

print("loading intra-locus distance data")
g_geom_feat_df <- readRDS("osu18_intra_locus_dist(anonymous).rds")
stopifnot(rownames(g_geom_feat_df) == rownames(g_annot_feat_df))

## build a list of the feature matrices that we will need
g_classifier_feature_matrices_list <- list(
	annot_feat=g_annot_feat_df,
	geom_feat=g_geom_feat_df 
)

## ============================== run invariant setup code =================================
source("../../machine_learning/cerenkov_incscript_setup_ml.R")


## ============================== make closures for classifiers  =================================

g_classifier_function_ranger_impurity <- g_make_classifier_function_ranger(p_nthread=g_par$nthreads_per_process,
																  p_get_perf_results=g_get_perf_results, 
																  p_feature_importance_type="impurity")

# g_classifier_function_ranger_impurity_corrected <- g_make_classifier_function_ranger(p_nthread=g_par$nthreads_per_process,
# 																		   p_get_perf_results=g_get_perf_results, 
# 																		   p_feature_importance_type="impurity_corrected")

g_classifier_function_ranger_permutation <- g_make_classifier_function_ranger(p_nthread=g_par$nthreads_per_process,
																  p_get_perf_results=g_get_perf_results, 
																  p_feature_importance_type="permutation")

g_classifier_functions_list <- list(
	ranger_impurity=g_classifier_function_ranger_impurity,
	# ranger_impurity_corrected=g_classifier_function_ranger_impurity_corrected,
	ranger_permutation=g_classifier_function_ranger_permutation
)

## ============================== assemble final list of feature matrices  =================================

## free up memory
rm(g_annot_feat_df)
rm(g_geom_feat_df)


## ============================== make hyperparameter lists  =================================

## ------------------ RF hyperparameter lists --------------

g_ranger_hp_grid <- g_make_hyperparameter_grid_list(list(probability=FALSE,
																mtry=14,
																num.trees=100,  
																weight_positive_class=1,
																replace=TRUE,
																sample.fraction=1))

g_classifier_list_ranger_impurity <- lapply(g_ranger_hp_grid,
											function(hp_cell) {
												list(classifier_feature_matrix_name="annot_feat",
													 classifier_function_name="ranger_impurity",
													 classifier_hyperparameter_set_type_name="g_ranger_hp_grid",
													 classifier_set_name="impurity",
													 classifier_hyperparameter_list=hp_cell, 
													 feature_reducer_function_name="FIT_LLR", 
													 feature_reducer_input_matrix_name="geom_feat", 
													 feature_reducer_hyperparameters_list=list())
											})

# g_classifier_list_ranger_impurity_corrected <- lapply(g_ranger_hp_grid,
# 														  function(hp_cell) {
# 														  	list(classifier_feature_matrix_name="annot_feat",
# 														  		 classifier_function_name="ranger_impurity_corrected",
# 														  		 classifier_hyperparameter_set_type_name="g_ranger_hp_grid",
# 														  		 classifier_set_name="impurity_corrected",
# 														  		 classifier_hyperparameter_list=hp_cell, 
# 														  		 feature_reducer_function_name="FIT_LLR", 
# 														  		 feature_reducer_input_matrix_name="geom_feat", 
# 														  		 feature_reducer_hyperparameters_list=list())
														  # })

g_classifier_list_ranger_permutation <- lapply(g_ranger_hp_grid,
											   function(hp_cell) {
											   	list(classifier_feature_matrix_name="annot_feat",
											   		 classifier_function_name="ranger_permutation",
											   		 classifier_hyperparameter_set_type_name="g_ranger_hp_grid",
											   		 classifier_set_name="permutation",
											   		 classifier_hyperparameter_list=hp_cell, 
											   		 feature_reducer_function_name="FIT_LLR", 
											   		 feature_reducer_input_matrix_name="geom_feat", 
											   		 feature_reducer_hyperparameters_list=list())
											   })




## ============================== assemble classifier list  =================================

## TODO: make a "g_check_classifier_list" function that checks for incorrect classifier function
## names, incorrect feature matrix names, etc.

g_classifier_list <- c(
	g_classifier_list_ranger_impurity,
	# g_classifier_list_ranger_impurity_corrected, 
	g_classifier_list_ranger_permutation
)

## ============================== run invariant machine-learning code =================================

                              
source("../../machine_learning/cerenkov_incscript_run_ml.R")
