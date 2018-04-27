## ================================================================================
## This file is part of the CERENKOV (Computational Elucidation of the
## REgulatory NonKOding Variome) software program. CERENKOV is subject to terms
## and conditions defined in the file "LICENSE.txt", which is part of this
## CERENKOV software distribution.
##
## Title       : cerenkov_incscript_run_ml.R
## 
## Description : Code that is sourced by a "cerenkov_script_XXXXXX.R" script in
##               order to run the machine-learning job.
##
##               I put code in this file if it would otherwise be boilerplate in
##               *every* script of the form "cerenkov_script_XXXXXX.R".
##
##               This code assumes the following R script have already been sourced:
##                 cerenkov_func_base.R
##                 cerenkov_func_aws.R (if EC2 is being used)
##                 cerenkov_incscript_setup_aws.R (if EC2 is being used)
##                 cerenkov_incscript_setup_ml.R
## 
## Author      : Stephen A. Ramsey, Oregon State University
##               https://github.com/saramsey
## ================================================================================

g_classifier_list <- g_classifier_list[order(sapply(g_classifier_list, "[[", "classifier_hyperparameter_set_type_name"),
                                             sapply(g_classifier_list, "[[", "classifier_set_name"))]

names(g_classifier_list) <- 1:length(g_classifier_list)

print(sprintf("Number distinct combinations of (algorithm, hyperparameter-tuple, feature-matrix) to process:  %d", length(g_classifier_list)))

## ============================ create the parallel cluster (fork or socket/ec2) =============================

library(parallel)

library(pbapply)
pboptions(type="txt")  ## force pblapply to display progress bar, even when using Rscript

g_cluster_info_list <- if (! g_par$flag_create_ec2_multi_cluster) {
    if (g_par$flag_create_fork_cluster) {
        ## our rule of thumb is to assign one process to each logical core, to keep things simple
        g_num_cores_use <- detectCores(logical=TRUE)
        print(sprintf("Number of cores detected: %d", g_num_cores_use))
        
        if (! is.null(g_par$override_num_fork_processes)) { g_num_cores_use <- g_par$override_num_fork_processes }

        func_get_temp_directory <- function() { file.path(tempdir(), Sys.getpid()) }
        
        list( cluster=makeForkCluster(nnodes=g_num_cores_use,
                                      outfile=g_par$debug_file_parallel),
              terminator_func=function() {stopCluster(g_cluster)},
              gettemp=func_get_temp_directory )
    }
    else {
        list( gettemp = g_get_temp_dir )
    }
} else {
    ## we are creating EC2 instances and setting up a cluster on them
    c(g_create_ec2_instances_and_get_ips(g_ec2_par,
                                         g_run_ec2_instance,
                                         g_get_and_configure_ip_address_for_ec2_instances,
                                         p_make_cluster=TRUE),
      list(gettemp=function() { file.path(tempdir()) } ))
}

g_cluster <- g_cluster_info_list$cluster
g_ec2_instances <- g_cluster_info_list$ec2_instances
g_ip_addresses <- g_cluster_info_list$ip_addresses
g_terminate_cluster <- g_cluster_info_list$terminator_func
##g_get_temp_directory <- g_cluster_info_list$gettemp

if (g_par$flag_create_fork_cluster) {
    clusterExport(g_cluster, varlist=list("g_make_make_temp_dir_if_doesnt_exist", "g_get_temp_dir"))
    g_make_temp_dir_if_doesnt_exist <- g_make_make_temp_dir_if_doesnt_exist(g_get_temp_dir)
    parLapply(g_cluster, 1:g_num_cores_use, g_make_temp_dir_if_doesnt_exist)
}

## randomize order of g_classifier_list, for load-balancing purposes

g_order_classifier_ids <- if (g_par$flag_randomize_classifier_order) { sample(length(g_classifier_list)) } else { 1:length(g_classifier_list) }

## ============================ set cluster RNG and export objects to cluster =============================

if (! is.null(g_cluster)) {
    ## initialize random number streams for the cluster workers
    print(sprintf("Setting cluster workers to use random number streams with seed %d", g_par$random_number_seed))
    clusterSetRNGStream(g_cluster, sample(.Machine$integer.max, size=1))
    ## export global variables to the cluster
    clusterExport(cl=g_cluster, varlist=ls())
}

## ============================ setup global lapply function for the ML work-list =============================

g_func_lapply_to_use_for_ml <-
if (! is.null(g_cluster)) {
    if (g_par$show_progress_bar) {
        function(p_X, p_FUNC) { pblapply(p_X, p_FUNC, cl=g_cluster) }
    } else {
        if (g_par$parallel_use_load_balancing) {
            function(p_X, p_FUNC) { parLapplyLB(g_cluster, p_X, p_FUNC) }
        } else {
            function(p_X, p_FUNC) { parLapply(g_cluster, p_X, p_FUNC) }
        }
    }
} else {
    if (g_par$show_progress_bar) { pblapply } else { lapply }
}


## set up notifications, via stdout or SMS
g_send_message_notification <- g_make_message_notifier_function(g_par$aws_sns_topic_arn)


## bundle up all the data structures in a simple "runner" function
g_classifier_runner_func <- function() {
    g_run_mult_classifs_mult_hyperparams_cv(g_classifier_list[g_order_classifier_ids], 
                                            g_classifier_functions_list,
                                            g_classifier_feature_matrices_list,
                                            g_label_vec,
                                            g_par$num_cv_replications,
                                            g_par$num_folds_cross_validation,
                                            g_func_lapply_to_use_for_ml,
                                            g_feature_reducer_functions_list,
                                            g_assign_cases_to_folds)
}

## wrap error handling around our "runner" function
g_classifier_runner_func_use <- if (! is.null(g_cluster) && (is.null(g_par$flag_debug_mode) || ! g_par$flag_debug_mode)) {
                                    g_make_classifier_runner_func_err_handling(g_classifier_runner_func,
                                                                               g_send_message_notification)
                                } else {
                                    print("not using try/catch")
                                    g_classifier_runner_func
                                }

## ============================ run the classifier and gather results =============================

print(sprintf("Starting ML at time: %s", Sys.time()))

g_ml_results <- g_classifier_runner_func_use()

## save the results to a file
if (! is.null(g_ml_results)) {
    save("g_par",
         "g_classifier_list",
         "g_ml_results",
         "g_args",
         "g_order_classifier_ids",
         file=g_output_file_name)
}

g_send_message_notification(sprintf("Finished ML at time: %s", Sys.time()))

if (! is.null(g_terminate_cluster)) { g_terminate_cluster() }

