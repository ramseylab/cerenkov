## This file is part of the CERENKOV (Computational Elucidation of the
## REgulatory NonKOding Variome) software program. CERENKOV is free software:
## you can redistribute it and/or modify it under the terms of Apache Software
## License version 2.0 (and incorporated into this package distribution in the
## file LICENSE).
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE.  See the Apache Software License version 2.0 for 
## details.
##
## Copyright Stephen Ramsey, Oregon State University
## 2017.03

## cerenkov_ml_tune_xgboost.R -- master script for tuning, training, and testing classifiers using R/parallel
##
## Author:  Stephen Ramsey
## 
## Packages required by this script:
##   PRROC, xgboost, Matrix, pbapply
##
## Packages conditionally required by this script, depending on parameter choices:
##   dismo
##
## R scripts called by this script:
##   cerenkov_ml_base.R
##
## Note:  do not use this program with Ranger version 0.6.6 (stability issues); use Ranger 0.6.0 only

## ============================ define global parameters =============================

g_par <- list(
    num_folds_cross_validation = 10,  
    num_cv_replications = 12,         
    flag_create_fork_cluster = TRUE,   
    override_num_fork_processes = 64,  
    notify_by_text_msg = TRUE,
    show_progress_bar = TRUE,
    flag_locus_sampling = TRUE,        ## set to false if you want SNP-level sampling
    flag_xgb_importance = FALSE,       ## if you set this to TRUE, make sure you set num_cv_replications=1 and num_folds_cross_validation=1
    random_number_seed = 1337,
    nthreads_per_process = 1,
    flag_randomize_workplan_order = FALSE
    )

source("cerenkov_ml_base.R")

## ============================== set random number seed =================================
print(sprintf("setting random number seed to: %d", g_par$random_number_seed))
set.seed(g_par$random_number_seed)

## ============================== load OSU feature data; check for problems with feature data =================================

print("loading OSU data")

load(file="features_OSU.Rdata")
g_snp_names <- rownames(g_feat_osu_df)
stopifnot(g_feature_matrix_is_OK(g_feat_osu_df))
g_label_vec <- as.integer(as.character(g_feat_osu_df$label))

## ============================== set up for global performance measures =================================

if (! require(PRROC, quietly=TRUE)) {
    stop("package PRROC is missing")
}

g_calculate_auroc <- g_make_performance_getter(roc.curve,
                                               g_interval_clipper,
                                               "auc")

g_calculate_aupvr <- g_make_performance_getter(pr.curve,
                                               g_interval_clipper,
                                               "auc.davis.goadrich")

## ============================== set up for locus sampling =================================

if (g_par$flag_locus_sampling) {
    print("loading SNP coordinates data")

    load("snp_coordinates.Rdata")  ## ==> g_snp_coords_df
    g_snp_locus_ids <- g_get_snp_locus_ids(g_snp_coords_df)
    rm(g_snp_coords_df)

    g_assign_cases_to_folds_by_locus <- g_make_assign_cases_to_folds_by_group(g_snp_locus_ids)

    unique_locus_ids <- unique(g_snp_locus_ids)
    
    g_locus_to_snp_ids_map_list <- setNames(lapply(unique_locus_ids, function(p_locus_id) {
        which(g_snp_locus_ids == p_locus_id)
    }), unique_locus_ids)

    g_calculate_avgrank <- g_make_calculate_avgrank_within_groups(g_snp_locus_ids,
                                                                  g_locus_to_snp_ids_map_list,
                                                                  g_rank_by_score_decreasing)

    g_get_perf_results <- g_make_get_perf_results(g_calculate_aupvr,
                                                  g_calculate_auroc,
                                                  g_calculate_avgrank)

} else {
    g_get_perf_results <- g_make_get_perf_results(g_calculate_aupvr,
                                                  g_calculate_auroc)
}

g_assign_cases_to_folds <- ifelse( g_par$flag_locus_sampling,
                                   c(g_assign_cases_to_folds_by_locus),
                                   c(g_assign_cases_to_folds_by_case) )[[1]]

library(Matrix)

g_feat_osu_matrix_sparse <- sparse.model.matrix(label ~ .-1, data=g_feat_osu_df)

## ============================== assemble final list of feature matrices  =================================


## build a list of the feature matrices that we will need
g_classifier_feature_matrices_list <- list(
   feat_OSU_sparsematrix=g_feat_osu_matrix_sparse
)

## free up memory
rm(g_feat_osu_df)
rm(g_feat_osu_matrix_sparse)

## ============================== make feature reducer function list  =================================

g_feature_reducer_functions_list <- list()


## ============================== precompute class balance ratio  =================================

g_label_vec_table <- table(g_label_vec)
g_class_count_ratio_negative_to_positive <- setNames(g_label_vec_table["0"]/g_label_vec_table["1"], NULL)
g_class_count_frac_positive <- setNames(g_label_vec_table["1"]/(g_label_vec_table["0"] + g_label_vec_table["1"]), NULL)

print(sprintf("Class label ratio (negative/positive): %f", g_class_count_ratio_negative_to_positive))


## ============================== make closures for classifiers  =================================

g_classifier_function_ranger <- g_make_classifier_function_ranger(p_nthread=g_par$nthreads_per_process,
                                                                   p_get_perf_results=g_get_perf_results)

g_custom_objective_for_training <- g_make_custom_xgboost_objective(g_class_count_ratio_negative_to_positive)

g_classifier_function_xgboost <- g_make_classifier_function_xgboost(p_nthread=g_par$nthreads_per_process,
                                                                    g_get_perf_results,
                                                                    p_feature_importance_type=NULL,
                                                                    p_objective_function="binary:logistic",
                                                                    p_case_group_ids=g_snp_locus_ids)

rm(g_snp_locus_ids)

g_classifier_functions_list <- list(
    XGB=g_classifier_function_xgboost
)

## ============================== make hyperparameter lists  =================================

## ------------------ xgboost hyperparameter lists --------------

g_hyperparameter_grid_list_xgb <- g_make_hyperparameter_grid_list(list(eta=c(0.01, 0.05, 0.1),
                                                                       nrounds=c(10, 20, 30),
                                                                       gamma=c(0, 1, 10, 100),
                                                                       subsample=c(0.5, 0.75, 0.85, 1),
                                                                       colsample_bytree=c(0.5, 0.75, 0.85),
                                                                       base_score=g_class_count_frac_positive,
                                                                       scale_pos_weight=c(0.125, 1, 8),
                                                                       max_depth=c(6, 8, 10)))


g_workplan_list_xgb_OSU <- lapply(g_hyperparameter_grid_list_xgb,
                                     function(p_hyp) {
                                         list(classifier_feature_matrix_name="feat_OSU_sparsematrix",
                                              classifier_function_name=ifelse(g_par$flag_xgb_importance, "XGB_importance", "XGB"),
                                              classifier_hyperparameter_set_type_name="XGB",
                                              workplan_set_name="OSU_XGB",
                                              classifier_hyperparameter_list=p_hyp)
                                     })




g_workplan_list <- c(
    g_workplan_list_xgb_OSU
)

g_workplan_list <- g_workplan_list[order(sapply(g_workplan_list, "[[", "classifier_hyperparameter_set_type_name"),
                                         sapply(g_workplan_list, "[[", "workplan_set_name"))]

names(g_workplan_list) <- 1:length(g_workplan_list)

print(sprintf("Number of workplans to process:  %d", length(g_workplan_list)))

## ============================ create the parallel cluster (fork or socket/ec2) =============================

library(parallel)

library(pbapply)
pboptions(type="txt")  ## force pblapply to display progress bar, even when using Rscript


if (g_par$flag_create_fork_cluster) {
    ## our rule of thumb is to assign one process to each logical core, to keep things simple
    g_num_cores_use <- detectCores(logical=TRUE)
    print(sprintf("Number of cores detected: %d", g_num_cores_use))
    if (! is.null(g_par$override_num_fork_processes)) {
        g_num_cores_use <- g_par$override_num_fork_processes
    }
    g_cluster <- makeForkCluster(g_num_cores_use,
                                 outfile="")
} else {
    g_cluster <- NULL
}
    
## randomize order of g_workplan_list, for load-balancing purposes

if (g_par$flag_randomize_workplan_order) {
    g_order_workplan_ids <- sample(length(g_workplan_list))
} else {
    g_order_workplan_ids <- 1:length(g_workplan_list)
}

## ============================ export to cluster =============================

if (! is.null(g_cluster)) {
    print(sprintf("Setting cluster workers to use random number streams with seed %d", g_par$random_number_seed))
    clusterSetRNGStream(g_cluster, g_par$random_number_seed)

    ## use this function for load-balancing parallelism (BUT DO NOT USE FOR DEBUGGING)
    g_func_lapply_cluster_LB <- function(p_X, p_FUNC) {
        parLapplyLB(g_cluster, p_X, p_FUNC)
    }

    ## use this for testing/debugging; exactly reproducible results
    if (g_par$show_progress_bar) {
        g_func_lapply_cluster <- function(p_X, p_FUNC) {
            pblapply(p_X, p_FUNC, cl=g_cluster)
        }
    } else {
        g_func_lapply_cluster <- function(p_X, p_FUNC) {
            parLapply(g_cluster, p_X, p_FUNC)
        }
    }

    clusterExport(cl=g_cluster, varlist=ls())

} else {
    g_func_lapply_cluster <- lapply
}

g_do_cluster_cleanup <- g_make_cluster_cleanup_function(g_cluster, NULL)

## set up notifications to go to STDOUT
g_send_message_notification <- print


## bundle up all the data structures in a simple "runner" function
g_classifier_runner_func <- function() {
    g_run_mult_classifs_mult_hyperparams_cv(g_workplan_list[g_order_workplan_ids],
                                            g_classifier_functions_list,
                                            g_classifier_feature_matrices_list,
                                            g_label_vec,
                                            g_par$num_cv_replications,
                                            g_par$num_folds_cross_validation,
                                            g_func_lapply_cluster,
                                            p_func_lapply_second_level=lapply,
                                            g_feature_reducer_functions_list,
                                            g_assign_cases_to_folds)
}


## wrap error handling around our "runner" function
g_classifier_runner_func_err_handling <- g_make_classifier_runner_func_err_handling(g_classifier_runner_func,
                                                                                    g_send_message_notification,
                                                                                    g_cluster)


## ============================ run the classifier and gather results =============================

        
print(sprintf("Starting ML at time: %s", Sys.time()))

g_ml_results <- g_classifier_runner_func_err_handling()

## save the results to a file
if (! is.null(g_ml_results)) {
    save("g_par",
         "g_workplan_list",
         "g_ml_results",
         "g_order_workplan_ids",
         file="cerenkov_ml_tune_xgboost.Rdata")
}

g_send_message_notification(sprintf("Finished ML at time: %s", Sys.time()))

g_do_cluster_cleanup()
