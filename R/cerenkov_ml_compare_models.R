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

## cerenkov_ml_compare_models.R -- compares the CERENKOV classifier to 10 other
## rSNP recognition models, on the OSU ground-truth datset of 15,331 SNPs.
##
## Author:  Stephen Ramsey
## 
## Packages required by this script:
##   PRROC, xgboost, ranger, Matrix, pbapply
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
    num_cv_replications = 200,
    flag_create_fork_cluster = TRUE,
    override_num_fork_processes = 50,  
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


## ============================== load feature data; check for problems with feature data =================================

library(Matrix)

g_feat_osu_matrix_sparse <- sparse.model.matrix(label ~ .-1, data=g_feat_osu_df)

print("loading GWAVA data")

load(file="features_GWAVA.Rdata") ## creates an R object called "g_feat_gwava_df"
stopifnot(g_feature_matrix_is_OK(g_feat_gwava_df))

stopifnot(g_snp_names == rownames(g_feat_gwava_df))
stopifnot(g_feat_osu_df$label == g_feat_gwava_df$label)

print("loading deepsea data")

load(file="features_deepsea.Rdata")

stopifnot(g_feature_matrix_is_OK(g_feat_deepsea_published_df))
stopifnot(rownames(g_feat_deepsea_published_df) == g_snp_names)
g_feat_deepsea_published_matrix <- data.matrix(g_feat_deepsea_published_df[, setdiff(names(g_feat_deepsea_published_df), "label")])

print("loading RSVP data")
      
load(file="features_RSVP.Rdata")

stopifnot(g_feature_matrix_is_OK(g_feat_rsvp_df))

stopifnot(rownames(g_feat_rsvp_df) == g_snp_names)
stopifnot(g_feat_rsvp_df$label == g_feat_osu_df$label)

g_feat_rsvp_matrix <- data.matrix(g_feat_rsvp_df[, setdiff(names(g_feat_rsvp_df), "label")])
g_feat_rsvp_imputed_df <- g_impute_missing_data_average_df(g_feat_rsvp_df)
    
rownames(g_feat_rsvp_imputed_df) <- rownames(g_feat_rsvp_df)
stopifnot(g_feat_rsvp_df$label == g_feat_rsvp_imputed_df$label)

print("loading DANQ data")
load(file="features_danq.Rdata")
stopifnot(rownames(g_feat_danq_df) == g_snp_names)
stopifnot(g_feat_danq_df$label == g_label_vec)
stopifnot(g_feature_matrix_is_OK(g_feat_danq_df))

g_feat_danq_matrix <- data.matrix(g_feat_danq_df[, setdiff(names(g_feat_danq_df), "label")])

print("loading deltaSVM data")
load(file="features_deltaSVM.Rdata")
stopifnot(rownames(g_feat_deltaSVM_df) == g_snp_names)
stopifnot(g_feat_deltaSVM_df$label == g_label_vec)
g_feat_deltaSVM_matrix_sparse <- sparse.model.matrix(label ~ .-1, data=g_feat_deltaSVM_df)
stopifnot(g_feature_matrix_is_OK(g_feat_deltaSVM_df))
    
## ===== load scores matrices for passthrough classifiers ======

print("loading CADD scores")
load("scores_cadd.Rdata")
stopifnot(names(g_scores_vec_cadd) == g_snp_names)
stopifnot(! is.na(g_scores_vec_cadd))

print("loading fitCons scores")
load("scores_fitcons.Rdata")
stopifnot(names(g_scores_vec_fitcons) == g_snp_names)
stopifnot(! is.na(g_scores_vec_fitcons))

print("loading eigen scores")
load("scores_eigen.Rdata")
stopifnot(names(g_scores_vec_eigen) == g_snp_names)
stopifnot(! is.na(g_scores_vec_eigen))

print("loading DANN scores")
load("scores_dann.Rdata")
stopifnot(names(g_scores_vec_dann) == g_snp_names)
stopifnot(! is.na(g_scores_vec_dann))

## ============================== assemble final list of feature matrices  =================================


## build a list of the feature matrices that we will need
g_classifier_feature_matrices_list <- list(
   feat_OSU_sparsematrix=g_feat_osu_matrix_sparse,
   feat_GWAVA_ranger=g_feat_gwava_df[, which(names(g_feat_gwava_df) != "label")],
   feat_RSVP_matrix=g_feat_rsvp_matrix,
   feat_RSVP_ranger=g_feat_rsvp_imputed_df[, which(names(g_feat_rsvp_imputed_df) != "label")],   
   feat_deepsea_matrix=g_feat_deepsea_published_matrix,
   feat_DANQ=g_feat_danq_matrix,
   feat_deltaSVM_matrix=g_feat_deltaSVM_matrix_sparse
)

## free up memory
rm(g_feat_osu_df)
rm(g_feat_osu_matrix_sparse)
rm(g_feat_rsvp_matrix)
rm(g_feat_rsvp_df)
rm(g_feat_rsvp_imputed_df)
rm(g_feat_deepsea_published_df)
rm(g_feat_deepsea_published_matrix)
rm(g_feat_gwava_df)
rm(g_feat_danq_df)
rm(g_feat_danq_matrix)
rm(g_feat_deltaSVM_df)
rm(g_feat_deltaSVM_matrix_sparse)

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

g_classifier_function_xgboost <- g_make_classifier_function_xgboost(p_nthread=g_par$nthreads_per_process,
                                                                    g_get_perf_results,
                                                                    p_feature_importance_type=NULL,
                                                                    p_objective_function="binary:logistic",
                                                                    p_case_group_ids=g_snp_locus_ids)

rm(g_snp_locus_ids)

g_classifier_passthrough_cadd_scores <- g_make_classifier_function_passthrough(g_scores_vec_cadd,
                                                                               p_get_perf_results=g_get_perf_results)

g_classifier_passthrough_fitCons_scores <- g_make_classifier_function_passthrough(g_scores_vec_fitcons,
                                                                                  p_get_perf_results=g_get_perf_results)

g_classifier_passthrough_eigen_scores <- g_make_classifier_function_passthrough(g_scores_vec_eigen,
                                                                                p_get_perf_results=g_get_perf_results)

g_classifier_passthrough_dann_scores <- g_make_classifier_function_passthrough(g_scores_vec_dann,
                                                                               p_get_perf_results=g_get_perf_results)

## free up memory
rm(g_scores_vec_cadd)
rm(g_scores_vec_fitcons)
rm(g_scores_vec_eigen)
rm(g_scores_vec_dann)

g_classifier_functions_list <- list(
    XGB=g_classifier_function_xgboost,
    ranger=g_classifier_function_ranger,
    CADD=g_classifier_passthrough_cadd_scores,
    fitCons=g_classifier_passthrough_fitCons_scores,
    eigen=g_classifier_passthrough_eigen_scores,
    DANN=g_classifier_passthrough_dann_scores
)

## ============================== make hyperparameter lists  =================================

## ------------------ xgboost hyperparameter lists --------------

g_hyperparameter_grid_list_xgb_best_aupvr <- g_make_hyperparameter_grid_list(list(eta=0.1,
                                                                                  nrounds=30,
                                                                                  gamma=10,
                                                                                  subsample=1,
                                                                                  colsample_bytree=0.85,
                                                                                  base_score=g_class_count_frac_positive,
                                                                                  scale_pos_weight=1,
                                                                                  max_depth=6))

g_hyperparameter_grid_list_xgb_best_avgrank <- g_make_hyperparameter_grid_list(list(eta=0.1,
                                                                                  nrounds=30,
                                                                                  gamma=100,
                                                                                  subsample=0.85,
                                                                                  colsample_bytree=0.85,
                                                                                  base_score=g_class_count_frac_positive,
                                                                                  scale_pos_weight=8,
                                                                                  max_depth=6))


g_hyperparameter_grid_list_xgb <- c(g_hyperparameter_grid_list_xgb_best_avgrank,
                                    g_hyperparameter_grid_list_xgb_best_aupvr)

g_hyperparameter_grid_list_xgb_DS <- g_make_hyperparameter_grid_list(list(eta=0.1,
                                                                          nrounds=10,
                                                                          alpha=20,
                                                                          lambda=2000,
                                                                          scale_pos_weight=8))

g_workplan_list_xgb_OSU <- lapply(g_hyperparameter_grid_list_xgb,
                                     function(p_hyp) {
                                         list(classifier_feature_matrix_name="feat_OSU_sparsematrix",
                                              classifier_function_name=ifelse(g_par$flag_xgb_importance, "XGB_importance", "XGB"),
                                              classifier_hyperparameter_set_type_name="XGB",
                                              workplan_set_name="OSU_XGB",
                                              classifier_hyperparameter_list=p_hyp)
                                     })


g_workplan_list_xgb_deepsea <- lapply(g_hyperparameter_grid_list_xgb_DS,
                                     function(p_hyp) {
                                         list(classifier_feature_matrix_name="feat_deepsea_matrix",
                                              classifier_function_name=ifelse(g_par$flag_xgb_importance, "XGB_importance", "XGB"),
                                              classifier_hyperparameter_set_type_name="XGBDS",
                                              workplan_set_name="Deepsea_XGB_published",
                                              classifier_hyperparameter_list=p_hyp)
                                     })

g_workplan_list_xgb_danq <- lapply(g_hyperparameter_grid_list_xgb_DS,
                                     function(p_hyp) {
                                         list(classifier_feature_matrix_name="feat_DANQ",
                                              classifier_function_name=ifelse(g_par$flag_xgb_importance, "XGB_importance", "XGB"),
                                              classifier_hyperparameter_set_type_name="XGBdanq",
                                              workplan_set_name="DANQ_XGB",
                                              classifier_hyperparameter_list=p_hyp)
                                     })

g_workplan_list_xgb_RSVP <- lapply(g_hyperparameter_grid_list_xgb,
                                   function(p_hyp) {
                                       list(classifier_feature_matrix_name="feat_RSVP_matrix",
                                            classifier_function_name=ifelse(g_par$flag_xgb_importance, "XGB_importance", "XGB"),
                                            classifier_hyperparameter_set_type_name="XGB",
                                            workplan_set_name="RSVP_XGB",
                                            classifier_hyperparameter_list=p_hyp)
                                   })

g_workplan_list_xgb_deltaSVM <- lapply(g_hyperparameter_grid_list_xgb,
                                   function(p_hyp) {
                                       list(classifier_feature_matrix_name="feat_deltaSVM_matrix",
                                            classifier_function_name=ifelse(g_par$flag_xgb_importance, "XGB_importance", "XGB"),
                                            classifier_hyperparameter_set_type_name="XGB",
                                            workplan_set_name="deltaSVM_XGB",
                                            classifier_hyperparameter_list=p_hyp)
                                   })

g_workplan_list_fitCons <- list(list(classifier_feature_matrix_name=NULL,
                                  classifier_function_name="fitCons",
                                  classifier_hyperparameter_set_type_name="passthrough",
                                  workplan_set_name="fitCons",
                                  classifier_hyperparameter_list=NULL))

g_workplan_list_CADD <- list(list(classifier_feature_matrix_name=NULL,
                                  classifier_function_name="CADD",
                                  classifier_hyperparameter_set_type_name="passthrough",
                                  workplan_set_name="CADD",
                                  classifier_hyperparameter_list=NULL))

g_workplan_list_eigen <- list(list(classifier_feature_matrix_name=NULL,
                                  classifier_function_name="eigen",
                                  classifier_hyperparameter_set_type_name="passthrough",
                                  workplan_set_name="eigen",
                                  classifier_hyperparameter_list=NULL))

g_workplan_list_DANN <- list(list(classifier_feature_matrix_name=NULL,
                                  classifier_function_name="DANN",
                                  classifier_hyperparameter_set_type_name="passthrough",
                                  workplan_set_name="DANN",
                                  classifier_hyperparameter_list=NULL))


## ------------------ ranger hyperparameter lists --------------

## WARNING:  do not set "probability=TRUE" for ranger; memory leak badness will result

g_workplan_list_RSVP_published <- list(list(classifier_feature_matrix_name="feat_RSVP_ranger",
                                            classifier_function_name="ranger",
                                            classifier_hyperparameter_set_type_name="ranger",
                                            workplan_set_name="RSVP_RF_published",
                                            classifier_hyperparameter_list=list(probability=FALSE,
                                                                                mtry=47,
                                                                                num.trees=100,  
                                                                                weight_positive_class=8,
                                                                                replace=TRUE,
                                                                                sample.fraction=1)))


g_workplan_list_gwava_published <- list(list(classifier_feature_matrix_name="feat_GWAVA_ranger",
                                             classifier_function_name="ranger",
                                             classifier_hyperparameter_set_type_name="ranger",
                                             workplan_set_name="GWAVA_RF_published",
                                             classifier_hyperparameter_list=list(probability=FALSE,
                                                                                 mtry=14,
                                                                                 num.trees=100,  
                                                                                 weight_positive_class=1,
                                                                                 replace=TRUE,
                                                                                 sample.fraction=1)))

## ============================== assemble workplan list  =================================

g_workplan_list <- c(
    g_workplan_list_xgb_OSU,
    g_workplan_list_xgb_deepsea,
    g_workplan_list_xgb_danq,
    g_workplan_list_xgb_RSVP,
    g_workplan_list_xgb_deltaSVM,
    g_workplan_list_gwava_published,
    g_workplan_list_RSVP_published,
    g_workplan_list_CADD,
    g_workplan_list_fitCons,
    g_workplan_list_eigen,
    g_workplan_list_DANN
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
         file="cerenkov_ml_compare_models.Rdata")
}

g_send_message_notification(sprintf("Finished ML at time: %s", Sys.time()))

g_do_cluster_cleanup()
