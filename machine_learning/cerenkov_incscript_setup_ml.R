## ================================================================================
## This file is part of the CERENKOV (Computational Elucidation of the
## REgulatory NonKOding Variome) software program. CERENKOV is subject to terms
## and conditions defined in the file "LICENSE.txt", which is part of this
## CERENKOV software distribution.
##
## Title       : cerenkov_incscript_setup_ml.R
## 
## Description : Code that is sourced by a "cerenkov_script_XXXXXX.R" script in
##               order to setup the machine-learning job.
##
##               I put code in this file if it would otherwise be boilerplate in
##               *every* script of the form "cerenkov_script_XXXXXX.R".
##
##               This code assumes "cerenkov_func_base.R" (and, if EC2 is being
##               used, "cerenkov_func_aws.R") have already been sourced.
## 
## Author      : Stephen A. Ramsey, Oregon State University
##               https://github.com/saramsey
## ================================================================================

## cannot have both progress bar and load-balancing
stopifnot( !g_par$flag_create_fork_cluster || (!g_par$show_progress_bar || !g_par$parallel_use_load_balancing ) )

## cannot do both forked cluster and EC2 distributed cluster at the same time
stopifnot( !g_par$flag_create_ec2_multi_cluster || !g_par$flag_create_fork_cluster)

## ------only if you are using an EC2 distributed cluster -----
if (g_par$flag_create_ec2_multi_cluster) {
    g_ec2_par <- g_configure_ec2_instances()
}

## ============================== set random number seed =================================
print(sprintf("setting random number seed to: %d", g_par$random_number_seed))
set.seed(g_par$random_number_seed)

## ============================== set up for locus sampling =================================

if (g_par$flag_locus_sampling) {
    print("loading SNP coordinates data")

    g_snp_coords_df <- readRDS("osu18_snp_coordinates.rds")  
    g_snp_locus_ids <- g_get_snp_locus_ids(g_snp_coords_df)
    rm(g_snp_coords_df)

    g_assign_cases_to_folds_by_locus <- g_make_assign_cases_to_folds_by_group(g_snp_locus_ids)

    g_locus_to_snp_ids_map_list <- g_make_group_to_case_ids_map_list(g_snp_locus_ids)
}

g_calculate_avgrank <- if (g_par$flag_locus_sampling) {
    g_make_calculate_avgrank_within_groups(g_snp_locus_ids,
                                           g_locus_to_snp_ids_map_list,
                                           g_rank_by_score_decreasing)
} else {
    NULL
}

## ============================== set up for global performance measures =================================

g_calculate_auroc <- g_make_performance_getter(PRROC::roc.curve,
                                               g_interval_clipper,
                                               "auc")

g_calculate_aupvr <- g_make_performance_getter(PRROC::pr.curve,
                                               g_interval_clipper,
                                               "auc.davis.goadrich")

g_get_perf_results <- g_make_get_perf_results(g_calculate_aupvr,
                                              g_calculate_auroc,
                                              g_calculate_avgrank)  

g_assign_cases_to_folds <- ifelse( g_par$flag_locus_sampling,
                                   c(g_assign_cases_to_folds_by_locus),
                                   c(g_assign_cases_to_folds_by_case) )[[1]]

## ============================== make feature reducer function list  =================================
g_feature_reducer_functions_list <- list(PLS=g_feature_reducer_pls,
										 BIN_LLR=g_feature_reducer_bin_llr,
										 FIT_LLR=g_feature_reducer_fit_llr)

## ============================== precompute class balance ratio  =================================

g_label_vec_table <- table(g_label_vec)
g_class_count_ratio_negative_to_positive <- setNames(g_label_vec_table["0"]/g_label_vec_table["1"], NULL)
g_class_count_frac_positive <- setNames(g_label_vec_table["1"]/(g_label_vec_table["0"] + g_label_vec_table["1"]), NULL)
print(sprintf("Class label ratio (negative/positive): %f", g_class_count_ratio_negative_to_positive))

## ============================== setup output file name  =================================

g_output_file_name <- if (! is.null(g_par$output_file_name)) {
                          g_par$output_file_name
                      } else {
                          sprintf("%s_%s_%d.Rdata",
                              g_par$output_file_base_name,
                              g_par$analysis_label,
                              g_par$random_number_seed)
                      }
print(sprintf("Output file name: %s", g_output_file_name))

if (file.exists(g_output_file_name)) {
    file.remove(g_output_file_name)
}

g_get_temp_dir <- if (! g_par$flag_create_ec2_multi_cluster) {
                      if ( g_par$flag_create_fork_cluster ) {
                          function() { file.path(tempdir(), Sys.getpid()) }
                      } else {
                          function() { file.path(tempdir()) }
                      }
                  } else {
                      function() { file.path(tempdir()) }
                  }

