## ================================================================================
## This file is part of the CERENKOV (Computational Elucidation of the
## REgulatory NonKOding Variome) software program. CERENKOV is subject to terms
## and conditions defined in the file "LICENSE.txt", which is part of this
## CERENKOV software distribution.
##
## Title       : cerenkov_func_base.R
## 
## Description : Define global functions for the CERENKOV project.
##
##               Each function in this file should satisfy the
##               following criteria:
##                 1. it should not not access the global environment
##                 2. it should be side effect-free, as much as possible
##                 3. function names should start with "g_" (for "global")
##                 4. quietly check for required packages at run-time
## 
## Author      : Stephen A. Ramsey, Oregon State University
##               https://github.com/saramsey
## ================================================================================

## cerenkov_ml_base.R:  here is where I put all the functions that do not reference any global variables

## --------------- Cerenkov-specific functions start here ---------

## p_classifier_list:  a list of "classifiers".  Each classifier is a list with elements:
##    classifier$classifier_feature_matrix_name: the integer ID of the feature matrix to use for this classifier
##    classifier$classifier_hyperparameter_list:  the list of hyperparameters to use for this classifier
##    classifier$classifier_function_name:  the function to be used for training the classifier
## p_classifier_functions_list:  list of functions, each with signature: function(p_classifier_feature_matrix,
##                                                                                p_classifier_hyperparameter_list,
##                                                                                p_label_vector,
##                                                                                p_inds_cases_test), returns list with four elements
##                                                                                (train_auroc, train_aupvr, test_auroc, test_aupvr)
## p_classifier_feature_matrices_list:  one feature matrix for each type of feature matrix data structure to use **NO CASE LABELS**
## p_case_label_vec:  numeric vector containing the feature labels (0 or 1 only)
## p_func_lapply:  function(p_X, p_FUNC) returning a list

library("plyr")  # rbind.fill

g_run_mult_classifs_mult_hyperparams_cv <- function(p_classifier_list,
                                                      p_classifier_functions_list,
                                                      p_classifier_feature_matrices_list,
                                                      p_case_label_vec,
                                                      p_num_cv_replications=1,
                                                      p_num_folds=10,
                                                      p_func_lapply=lapply,
                                                      p_feature_reducer_functions_list=NULL,
                                                      p_assign_cases_to_folds) {
    
    ## get list of unique classifier hyperparameter set type names
    classifier_hyperparameter_set_type_names_unique <- sort(unique(unlist(lapply(p_classifier_list,
                                                                             function(p_classifier) {
                                                                                 p_classifier$classifier_hyperparameter_set_type_name
                                                                             }))))

    ## check if there is at least one classifier on the classifier list
    stopifnot(length(p_classifier_list) > 0)

    ## we need to know how many cases there are, in order to assign the cases to cross-validation folds
    num_cases <- unique(sapply(p_classifier_feature_matrices_list, nrow))
    if (length(num_cases) > 1) {
        stop("all classifier feature matrices must have equal numbers of cases")
    }

    ## create a list of length equal to the number of replications; each list contains fold assignments for all SNPs
    replications_fold_assignments_list <- replicate(p_num_cv_replications,
                                                    p_assign_cases_to_folds(p_num_folds=p_num_folds,
                                                                            p_case_label_vec=p_case_label_vec),
                                                    simplify=FALSE)

    ## make a workplan list containing triples of classifier ID, replication ID, and CV fold ID
    work_list <- lapply(setNames(lapply(data.frame(t(expand.grid(1:length(p_classifier_list),
                                                                 1:p_num_cv_replications,
                                                                 1:p_num_folds))),
                                 setNames, c("classifier_id", "replication_id", "fold_id")),
                                 NULL), as.list)

    ml_global_results_list <- p_func_lapply(work_list,
                                            function(p_vec_inds_work) {
                                                classifier_id <- p_vec_inds_work$classifier_id
                                                fold_id <- p_vec_inds_work$fold_id
                                                replication_id <- p_vec_inds_work$replication_id
                                                
                                                fold_assignments <- replications_fold_assignments_list[[replication_id]]
                                                inds_cases_test <- which(fold_id == fold_assignments)
                                                if (length(inds_cases_test) == length(fold_assignments)) {
                                                    ## this means we have num_folds=1, i.e., no cross-validation; train on all cases, so set test cases to empty vector
                                                    inds_cases_test <- c()
                                                }

                                                classifier_list <- p_classifier_list[[classifier_id]]

                                                ## need to know the feature matrix name, so we can retrieve the feature matrix
                                                classifier_feature_matrix_name <- classifier_list$classifier_feature_matrix_name

                                                ## need the classifier's hyperparameter list
                                                classifier_hyperparameter_list <- classifier_list$classifier_hyperparameter_list

                                                ## need the classifier's hyperparameter set name, which dictates what top-level list slot the results go into
                                                classifier_hyperparameter_set_type_name <- classifier_list$classifier_hyperparameter_set_type_name

                                                ## need the classifier's function name so we can look up the classifier function
                                                classifier_function_name <- classifier_list$classifier_function_name

                                                ## need the classifier function so we can run the train/test cycle
                                                classifier_function <- p_classifier_functions_list[[classifier_function_name]]
                                                
                                                save_hyperparameter_list <- classifier_hyperparameter_list
                                                
                                                if (is.null(classifier_list$feature_reducer_function_name)) {
                                                    ## this is the standard case, the classifier doesn't call for using a feature matrix reducer
                                                    
                                                    if (! is.null(classifier_feature_matrix_name)) {
                                                        ## this is the standard case, we don't have a null feature matrix (means a passthrough classifier)
                                                        feature_matrix <- p_classifier_feature_matrices_list[[classifier_feature_matrix_name]]
                                                        
                                                        if (is.null(feature_matrix)) { stop(sprintf("feature matrix %s missing",
                                                                                                    classifier_feature_matrix_name)) }
                                                    } else {
                                                        ## in this case we are using a passthrough classifier (like CADD, Eigen, or fitCons)
                                                        feature_matrix <- NULL
                                                        classifier_feature_matrix_name <- "NA"
                                                    }
                                                    
                                                } else {
                                                    
                                                    ## we are using a "supervised" feature reducer function, like PLS
                                                    base_feature_matrix <- p_classifier_feature_matrices_list[[classifier_feature_matrix_name]]
                                                    if (is.null(base_feature_matrix)) { stop(sprintf("feature matrix %s missing",
                                                                                                     classifier_feature_matrix_name)) }
                                                    feature_reducer_function_name <- classifier_list$feature_reducer_function_name
                                                    feature_reducer_function <- p_feature_reducer_functions_list[[feature_reducer_function_name]]
                                                    stopifnot( ! is.null(feature_reducer_function))
                                                    feature_reducer_input_matrix_name <- classifier_list$feature_reducer_input_matrix_name
                                                    stopifnot( ! is.null(feature_reducer_input_matrix_name))
                                                    feature_reducer_input_matrix <- p_classifier_feature_matrices_list[[feature_reducer_input_matrix_name]]
                                                    if (is.null(feature_reducer_input_matrix)) { stop(sprintf("feature matrix %s missing",
                                                                                                              feature_reducer_input_matrix_name)) }
                                                    stopifnot( ! is.null(feature_reducer_input_matrix))
                                                    feature_reducer_hyperparameters_list <- classifier_list$feature_reducer_hyperparameters_list
                                                    stopifnot( ! is.null(feature_reducer_hyperparameters_list))
                                                    
                                                    ## call the feature reducer
                                                    reduced_feature_matrix <- do.call(feature_reducer_function,
                                                                                      c(list(p_input_feature_matrix=feature_reducer_input_matrix,
                                                                                             p_case_label_vec=p_case_label_vec,
                                                                                             p_inds_cases_test=inds_cases_test),
                                                                                        feature_reducer_hyperparameters_list))
                                                    
                                                    ## combine the reduced feature matrix with the base feature matrix
                                                    if ("sparseMatrix" %in% methods::is(base_feature_matrix)) {
                                                    	feature_matrix <- Matrix::cBind(base_feature_matrix, reduced_feature_matrix)
                                                    } else {
                                                        feature_matrix <- cbind(base_feature_matrix, reduced_feature_matrix)
                                                    }
                                                    
                                                    classifier_feature_matrix_name <- paste(classifier_feature_matrix_name,
                                                                                            feature_reducer_input_matrix_name, sep="_")

                                                    save_hyperparameter_list <- c(save_hyperparameter_list,
                                                                                  feature_reducer_hyperparameters_list)
                                                } 

                                                classifier_custom_objective_function_parameters_list <- classifier_list$classifier_custom_objective_function_parameters_list
                                                if (! is.null(classifier_custom_objective_function_parameters_list)) {
                                                    save_hyperparameter_list <- c(save_hyperparameter_list, classifier_custom_objective_function_parameters_list)
                                                }
                                                
                                                ## train/test the classifier for the specified hyperparameters
                                                classifier_run_time <- system.time(
                                                    classifier_ret_list <- classifier_function(p_classifier_feature_matrix=feature_matrix,
                                                                                               p_classifier_hyperparameter_list=classifier_hyperparameter_list,
                                                                                               p_label_vector=p_case_label_vec,
                                                                                               p_inds_cases_test=inds_cases_test,
                                                                                               p_custom_objective_function_parameters_list=classifier_custom_objective_function_parameters_list) )

                                                feat_import_scores <- classifier_ret_list$feat_import_scores
                                                classifier_ret_list$feat_import_scores <- NULL
                                                
                                                if (is.null(save_hyperparameter_list)) {
                                                    save_hyperparameter_list <- list(classifier_hyperparameters.="")
                                                }
                                             	
                                                ## create a list of results
                                                ## WARNING:  DO **NOT** ALTER THE ORDER OF THESE LIST ELEMENTS, IT WILL BREAK BRITTLE DOWNSTREAM ANALYSIS CODE:
                                                list(performance_results=data.frame(c(classifier_ret_list,
                                                                                      list(classifier_name=classifier_function_name,
                                                                                           classifier_feature_matrix_name=classifier_feature_matrix_name,
                                                                                           classifier_hyperparameter_set_type_name=classifier_hyperparameter_set_type_name,
                                                                                           classifier_hyperparameters=save_hyperparameter_list,
                                                                                           classifier_run_time=setNames(classifier_run_time[1], NULL),
                                                                                           classifier_set_name=classifier_list$classifier_set_name,
                                                                                           classifier_id=classifier_id,
                                                                                           replication_id=replication_id,
                                                                                           cv_fold_id=fold_id)),
                                                                                    stringsAsFactors=FALSE),
                                                     feat_import_scores=feat_import_scores)
                                            })

    ## invert res_list so that the classifier ID is the top level, and the fold ID is the second level
    ml_global_results_list_wp_top <- lapply(1:length(p_classifier_list),
                                            function(p_classifier_list_id) {
                                                ml_global_results_list[which(sapply(ml_global_results_list,
                                                                                    function(p_list) {
                                                                                        p_list$performance_results$classifier_id
                                                                                    })==p_classifier_list_id)]
                                            })
    
    ## get the hyperparameter set type name for each classifier
    classifier_hyperparameter_type_names <- sapply(ml_global_results_list_wp_top,
            function(p_res_list_for_classifier) {
                p_res_list_for_classifier[[1]]$performance_results$classifier_hyperparameter_set_type_name })
    
    ## divide classifier list into performance results data frames, organized by hyperparameter-set type name
    ml_global_results_list_hpts_top <- setNames(lapply(classifier_hyperparameter_set_type_names_unique,
            function(p_classifier_hyperparameter_type_name) {
                hptds_list_list <- ml_global_results_list_wp_top[which(classifier_hyperparameter_type_names == p_classifier_hyperparameter_type_name)]
                # `feature_reducer_hyperparameters_list`` would add extra items to `performance_result` for classifiers with reducers
                # `rbind.fill` will fill missing columns with NA
                do.call(rbind.fill, lapply(hptds_list_list, function(hptds_list) {
                    do.call(rbind, lapply(hptds_list, "[[", "performance_results"))
                }))
            }), classifier_hyperparameter_set_type_names_unique)

    ## extract feature importance scores, if they were gathered
    feat_impt_scores <- setNames(lapply(classifier_hyperparameter_set_type_names_unique,
            function(p_classifier_hyperparameter_type_name) {
                hptds_list_list <- ml_global_results_list_wp_top[which(classifier_hyperparameter_type_names == p_classifier_hyperparameter_type_name)]
                setNames(lapply(hptds_list_list, function(hptds_list) {
                    lapply(hptds_list, "[[", "feat_import_scores")
                }), sapply(hptds_list_list, function(htps_list) {
                    htps_list[[1]]$performance_results$classifier_set_name
                }))
            }), classifier_hyperparameter_set_type_names_unique)

    list(performance_results=ml_global_results_list_hpts_top,
         feature_impt_scores=feat_impt_scores)
}

g_rank_by_score_decreasing <- function(x) {
    length(x) - rank(x, ties.method="average") + 1
}

g_impute_missing_data_average_df <- function(p_df) {
    setNames(do.call(cbind, lapply(p_df,
                                   function(mycol) {
                                       mycol[is.na(mycol)] <- mean(mycol, na.rm=TRUE)
                                       data.frame(mycol)
                                   })),
             names(p_df))
}

g_make_calculate_avgrank_within_groups <- function(p_case_group_ids, p_group_to_case_ids_map_list, p_rank_by_score_decreasing) {
    stopifnot(is.list(p_group_to_case_ids_map_list))

    p_case_group_ids
    p_group_to_case_ids_map_list
    
    function(p_case_scores, p_labels, p_case_ids) {
        ## need the indices of the rSNPs in the full set of 15,331 SNPs
        inds_pos <- p_case_ids[which(p_labels==1)]

        inds_map <- match(1:length(p_case_group_ids), p_case_ids)
        
        ## for each rSNP....
        mean(sapply(inds_pos, function(p_ind_pos) {
            ## get the group ID of the rSNP
            group_id <- p_case_group_ids[p_ind_pos]

            ## get the set of all SNPs in that group
            group_case_inds <- p_group_to_case_ids_map_list[[group_id]]
            stopifnot(! is.null(group_case_inds))

            inds_cases_within_fold_set <- inds_map[group_case_inds]
            stopifnot(! is.na(inds_cases_within_fold_set))
            
            ## use group_case_inds to index into 
            scores_for_group_cases <- p_case_scores[inds_cases_within_fold_set]
            stopifnot(! is.na(scores_for_group_cases))
            
            labels_for_group_cases <- p_labels[inds_cases_within_fold_set]
            stopifnot(! is.na(labels_for_group_cases))

            group_case_inds_analyze <- c(p_ind_pos, group_case_inds[which(labels_for_group_cases==0)])
            stopifnot(! is.na(group_case_inds_analyze))
            
            scores_to_analyze <- p_case_scores[inds_map[group_case_inds_analyze]]
            stopifnot(! is.na(scores_to_analyze))

            p_rank_by_score_decreasing(scores_to_analyze)[1]
        }))
    }
}



g_interval_clipper <- function(u) {
    pmax(pmin(u, 1.0), 0.0)
}

g_make_performance_getter <- function(p_performance_auc_calculator_func, p_interval_clipper_func, p_result_slot_name) {
    function(p_prediction_scores_vec, p_labels_vec) {
        p_interval_clipper_func(p_performance_auc_calculator_func(p_prediction_scores_vec[which(p_labels_vec==1)],
                                                                  p_prediction_scores_vec[which(p_labels_vec==0)])[[p_result_slot_name]])
    }
}

g_make_get_perf_results <- function(p_calculate_aupvr,
                                    p_calculate_auroc,
                                    p_calculate_avgrank=NULL) {

    ## obtain AUROC and AUPVR performance results given a vector of scores and a numeric vector of binary class labels
    function(p_prediction_scores, p_labels, p_case_inds)  {
        stopifnot(length(p_prediction_scores) == length(p_labels))
        stopifnot(is.vector(p_prediction_scores))
        stopifnot(is.numeric(p_prediction_scores))
        stopifnot(is.vector(p_labels))
        stopifnot(is.numeric(p_labels))

        ## randomize in case of weird bias due to tied scores
        rand_order <- sample(length(p_labels))
        p_prediction_scores_rand <- p_prediction_scores[rand_order]
        p_labels_rand <- p_labels[rand_order]
        
        aupvr <- p_calculate_aupvr(p_prediction_scores_rand,
                                   p_labels_rand)

        auroc <- p_calculate_auroc(p_prediction_scores_rand,
                                   p_labels_rand)
        
        res_list <- list(aupvr=aupvr,
                         auroc=auroc)

        if (! is.null(p_calculate_avgrank)) {

            avgrank <- p_calculate_avgrank(p_prediction_scores,
                                           p_labels,
                                           p_case_inds)
                                 
            res_list <- c(res_list, list(avgrank=avgrank))
        }
        
        res_list
    }
}

## function constructor for the "ranger" classifier
## NOTE:  ranger leaks memory pretty badly when "probability=TRUE"
g_make_classifier_function_ranger <- function(p_nthread=1, p_get_perf_results, p_feature_importance_type=NULL) {

    if (! suppressWarnings( require(ranger, quietly=TRUE) ) ) {
        stop("package ranger is missing")
    }

    package_version <- function(p_package_name) {
        ip_df <- installed.packages()
        ip_df[which(ip_df[,"Package"]==p_package_name),"Version"]
    }
    
    if (package_version("ranger") == "0.6.6") {
        stop("you are running Ranger version 0.6.6 -- downgrade to version 0.6.0 (from CRAN)")
    }
    
    function(p_classifier_feature_matrix,
             p_classifier_hyperparameter_list,
             p_label_vector,
             p_inds_cases_test,
             p_custom_objective_function_parameters_list=NULL) {

        if (! require(ranger, quietly=TRUE)) {
            stop("package ranger is missing")
        }

        stopifnot( ! is.null(p_classifier_feature_matrix))
        stopifnot(length(grep("data.frame|matrix", methods::is(p_classifier_feature_matrix))) > 0)
        stopifnot(is.vector(p_label_vector))
        stopifnot(is.numeric(p_label_vector))
        stopifnot(sort(unique(p_label_vector)) == c(0,1))
        stopifnot(colnames(p_classifier_feature_matrix) != "label")

        ncases <- nrow(p_classifier_feature_matrix)
        
        inds_cases_training <- setdiff(1:ncases, p_inds_cases_test)

        train_labels <- p_label_vector[inds_cases_training]
        test_labels <- p_label_vector[p_inds_cases_test]

        rf_data_matrix <- cbind(p_classifier_feature_matrix,
                                label=factor(p_label_vector))

        if (! is.null(p_classifier_hyperparameter_list$weight_positive_class) &&
            p_classifier_hyperparameter_list$weight_positive_class != 1) {
            weights_for_classes <- c(1, p_classifier_hyperparameter_list$weight_positive_class)
            case_weights_train <- weights_for_classes[p_label_vector[inds_cases_training] + 1]
        }
        else {
            case_weights_train <- NULL
        }

        p_classifier_hyperparameter_list$weight_positive_class <- NULL
        
        hyperparameter_list_for_ranger <- c(p_classifier_hyperparameter_list,
                                            list(case.weights=case_weights_train))

        probability_tree_setting <- p_classifier_hyperparameter_list$probability
        if (is.null(probability_tree_setting)) {
            probability_tree_setting <- FALSE
        }

        param_list <- c(list(dependent.variable.name="label",
                           data=rf_data_matrix[inds_cases_training, ],
                           classification=TRUE,
                           respect.unordered.factors=TRUE,
                           verbose=FALSE,
                           num.threads=p_nthread),
                        hyperparameter_list_for_ranger)
        
        if (! is.null(p_feature_importance_type)) {
            param_list <- c(param_list, list(importance=p_feature_importance_type))
        }
        
        rf_model <- do.call(ranger, param_list)

        if (! is.null(p_feature_importance_type)) {
            feat_import_scores <- importance(rf_model)
        } else {
            feat_import_scores <- NULL
        }
        
	train_pred_res <- predict(rf_model,
                                  data=rf_data_matrix[inds_cases_training,],
		                  predict.all=!probability_tree_setting)$predictions

        ## workaround because ranger returns per-tree scores that are "1 or 2" for probability=FALSE, and between 0 and 1 for probability=TRUE

        if (probability_tree_setting) {
            train_pred_scores <- train_pred_res[,2]
	} else {
            train_pred_scores <- apply(train_pred_res - 1, 1, mean)
        }

	train_perf_results <- p_get_perf_results(train_pred_scores, train_labels, inds_cases_training)

        if (length(test_labels) > 0) {
            test_pred_res <- predict(rf_model,
                                     data=rf_data_matrix[p_inds_cases_test,],
                                     predict.all=!probability_tree_setting)$predictions

            if (probability_tree_setting) {
                test_pred_scores <- test_pred_res[,2]
            } else {
                test_pred_scores <- apply(test_pred_res - 1, 1, mean)
            }

            if (any(is.nan(test_pred_scores))) {
                print(sprintf("number of NaN values in ranger prediction scores: %d", length(which(is.nan(test_pred_scores)))))
                print("ranger hyperparameter set: ")
                print(p_classifier_hyperparameter_list)
                test_pred_scores[is.nan(test_pred_scores)] <- 0.5
            }
        
            test_perf_results <- p_get_perf_results(test_pred_scores, test_labels, p_inds_cases_test)
            
        } else {
            test_perf_results <- list(auroc=NA, aupvr=NA)
            if (! is.null(train_perf_results$avgrank)) {
                test_perf_results <- c(test_perf_results,
                                       list(avgrank=NA))
            }
        }

        ret_list <- list(train_auroc=train_perf_results$auroc,
                         train_aupvr=train_perf_results$aupvr,
                         test_auroc=test_perf_results$auroc,
                         test_aupvr=test_perf_results$aupvr)

        if (! is.null(train_perf_results$avgrank)) {
            ret_list <- c(ret_list,
                          list(train_avgrank=train_perf_results$avgrank,
                               test_avgrank=test_perf_results$avgrank))
        }
        
        if (! is.null(feat_import_scores)) {
            ret_list <- c(ret_list,
                          list(feat_import_scores=feat_import_scores))
        }

        ret_list
    }
}       


## Requires: Rcpp
## The purpose of this function is to enable running an Rcpp-defined function in
## a "worker process" via R's parallel framework, without getting the error
## "NULL value passed as symbol address". This code is based on a code example
## that Roman Francois posted on the gmane.comp.lang.r.rcpp mailing list, 26th
## Sept. 2013.
g_make_cxx_function_with_lazy_compile <- function(p_code, p_func_get_temp_directory, ...) {
    function(...) {
        temp_directory <- p_func_get_temp_directory()
        do.call(Rcpp::cppFunction, list(code=p_code, cacheDir=temp_directory))(...)
    }
}

## Workaround because Rcpp cppFunction doesn't play nice with multicore parallel;
## by default, all worker processes use the same cache directory and that leads
## to badness because they are randomly erasing one another's cached CPP files.
g_make_make_temp_dir_if_doesnt_exist <- function(p_get_temp_directory) {
    function(...) {
        temp_directory <- p_get_temp_directory()
        if (! file.exists(temp_directory)) {
            dir.create(temp_directory)
        }
    }
}

## function constructor for the "xgboost" classifier; this is hard-coded to use "NA" to denote missing data
g_make_classifier_function_xgboost <- function(p_nthread=1,
                                               p_get_perf_results,
                                               p_feature_importance_type=NULL,
                                               p_make_objective_function=function(...){"binary:logistic"},
                                               p_case_group_ids=NULL,
                                               p_verbose=0) {

    if (! require(xgboost, quietly=TRUE)) {
        stop("package xgboost is missing")
    }
    
    nthread_third_level_list <- list(nthread=p_nthread)

    avg_group_size_ncases <- mean(table(p_case_group_ids))

    function(p_classifier_feature_matrix,
             p_classifier_hyperparameter_list,
             p_label_vector,
             p_inds_cases_test,
             p_custom_objective_function_parameters_list=NULL) {

        ## load xgboost package; we will need it to assign cases to cross-validation folds
        if (! require(xgboost, quietly=TRUE)) {
            stop("package xgboost is missing")
        }
    
        stopifnot( ! is.null(p_classifier_feature_matrix))
        stopifnot(length(grep("dgCMatrix|matrix", methods::is(p_classifier_feature_matrix))) > 0)
        stopifnot(is.vector(p_label_vector))
        stopifnot(is.numeric(p_label_vector))
        stopifnot(sort(unique(p_label_vector)) == c(0,1))
        stopifnot(colnames(p_classifier_feature_matrix) != "label")

        ncases <- nrow(p_classifier_feature_matrix)
        
        inds_cases_training <- setdiff(1:ncases, p_inds_cases_test)
        
        train_labels <- p_label_vector[inds_cases_training]
        test_labels <- p_label_vector[p_inds_cases_test]

        ## when we upgrade to newer xgboost that has row-based subsetting of xgb.DMatrix objects, this code to move to cerenkov_ml.R
        model_data <- xgb.DMatrix(p_classifier_feature_matrix[inds_cases_training, ], label=train_labels, missing=NA)

        if (! is.null(p_custom_objective_function_parameters_list)) {
            custom_objective_function_parameters_list <- p_custom_objective_function_parameters_list
        } else {
            custom_objective_function_parameters_list <- list()
        }
        
        objective_function <- p_make_objective_function(inds_cases_training,
                                                        p_custom_objective_function_parameters_list=custom_objective_function_parameters_list)
        
        xgb_params <- c(p_classifier_hyperparameter_list,
                       list(objective=objective_function))

        xgb_model <- xgb.train(params=xgb_params,
                               data=model_data,
                               nrounds=p_classifier_hyperparameter_list$nrounds,
                               verbose=p_verbose,
                               nthread=p_nthread)
        
        if (! is.null(p_feature_importance_type)) {
            feat_import_scores <- as.data.frame(xgb.importance(model=xgb_model))
            feat_names_vec <- colnames(p_classifier_feature_matrix)
            feat_import_scores$Feature <- feat_names_vec[1 + as.integer(feat_import_scores$Feature)]
        } else {
            feat_import_scores <- NULL
        }
        
        train_pred_scores <- predict(xgb_model,
                                     p_classifier_feature_matrix[inds_cases_training, ],
                                     missing=NA)
        
        ## \/\/\/\/   this is a work-around for bug # 1503 in xgboost (see GitHub); might have been fixed in latest release of R xgboost package
        n_train_pred_scores <- length(train_pred_scores)
        if (length(train_pred_scores) < length(train_labels)) {
            train_labels <- train_labels[1:length(train_pred_scores)]
        }
        ## /\/\/\/\   this is a work-around for bug # 1503 in xgboost (see GitHub); might have been fixed in latest release of R xgboost package

        stopifnot(length(train_pred_scores) == length(train_labels))

        train_perf_results <- p_get_perf_results(train_pred_scores, train_labels, inds_cases_training)

        if (length(test_labels) > 0) {
            test_pred_scores <- predict(xgb_model,
                                        p_classifier_feature_matrix[p_inds_cases_test,],
                                        missing=NA)

            ## \/\/\/\/   this is a work-around for bug # 1503 in xgboost (see GitHub); might have been fixed in latest release of R xgboost package
            n_test_pred_scores <- length(test_pred_scores)
            if (length(test_pred_scores) < length(test_labels)) {
                test_labels <- test_labels[1:length(test_pred_scores)]
            }
            ## /\/\/\/\   this is a work-around for bug # 1503 in xgboost (see GitHub); might have been fixed in latest release of R xgboost package

            stopifnot(length(test_labels) == length(test_pred_scores))
            
            test_perf_results <- p_get_perf_results(test_pred_scores, test_labels, p_inds_cases_test)
        } else {
            test_perf_results <- list(auroc=NA, aupvr=NA)
            if (! is.null(train_perf_results$avgrank)) {
                test_perf_results <- c(test_perf_results,
                                       list(avgrank=NA))
            }
       }
        
        ret_list <- list(train_auroc=train_perf_results$auroc,
                         train_aupvr=train_perf_results$aupvr,
                         test_auroc=test_perf_results$auroc,
                         test_aupvr=test_perf_results$aupvr)

        if (! is.null(train_perf_results$avgrank)) {
            ret_list <- c(ret_list,
                          list(train_avgrank=train_perf_results$avgrank,
                               test_avgrank=test_perf_results$avgrank))
        }
        
        if (! is.null(feat_import_scores)) {
            ret_list <- c(ret_list, list(feat_import_scores=feat_import_scores))
        }

        ret_list
    }
}

g_make_classifier_function_passthrough <- function(p_scores_vec, p_get_perf_results) {
    stopifnot(! is.null(p_scores_vec))
    score_vec_case_names <- names(p_scores_vec)
    stopifnot(! is.null(score_vec_case_names))

    function(p_classifier_feature_matrix,
             p_classifier_hyperparameter_list,
             p_label_vector,
             p_inds_cases_test) {

        ncases <- length(p_scores_vec)
        inds_cases_training <- setdiff(1:ncases, p_inds_cases_test)

        train_labels <- p_label_vector[inds_cases_training]
        test_labels <- p_label_vector[p_inds_cases_test]

        train_pred_scores <- p_scores_vec[inds_cases_training]
        train_perf_results <- p_get_perf_results(train_pred_scores, train_labels, inds_cases_training)

        test_pred_scores <- p_scores_vec[p_inds_cases_test]
        test_perf_results <- p_get_perf_results(test_pred_scores, test_labels, p_inds_cases_test)
         
        ret_list <- list(train_auroc=train_perf_results$auroc,
                         train_aupvr=train_perf_results$aupvr,
                         test_auroc=test_perf_results$auroc,
                         test_aupvr=test_perf_results$aupvr)

        if (! is.null(train_perf_results$avgrank)) {
            ret_list <- c(ret_list,
                          list(train_avgrank=train_perf_results$avgrank,
                               test_avgrank=test_perf_results$avgrank))
        }

        ret_list
    }
}

g_feature_reducer_pls <- function(p_input_feature_matrix,
                                  p_case_label_vec,
                                  p_inds_cases_test,
                                  p_num_components) {
    stopifnot( ! is.null(p_num_components) )
    
    if (! suppressWarnings( require(pls, quietly=TRUE) ) ) { stop("package pls is missing") }
    if (! require(methods, quietly=TRUE)) { stop("package methods is missing") }
    
    p_input_feature_matrix_as_df <- data.frame(pls.x=I(p_input_feature_matrix))
    p_input_feature_matrix_as_df$label <- I(model.matrix(~y-1, data.frame(y=factor(p_case_label_vec))))
    
    inds_cases_train <- setdiff(1:length(p_case_label_vec), p_inds_cases_test)
    
    pls_model <- cppls(label ~ pls.x, data=p_input_feature_matrix_as_df[inds_cases_train,], ncomp=p_num_components)
    
    pls_pred <- predict(pls_model, p_input_feature_matrix_as_df, type="scores")
}

g_feature_reducer_bin_llr <- function(p_input_feature_matrix,
								  p_case_label_vec,
								  p_inds_cases_test,
                                  p_dist_col_names,
								  p_num_bins) {
    # For CERENKOV2 first submission, 8 types of intralocus distances were used
    # 
    # p_dist_col_names <- c(
	# 	"intra_locus_dist_avg_canberra",
	# 	"intra_locus_dist_avg_canberra_scaled",
	# 	"intra_locus_dist_avg_euclidean",
	# 	"intra_locus_dist_avg_euclidean_scaled",
	# 	"intra_locus_dist_avg_manhattan",
	# 	"intra_locus_dist_avg_manhattan_scaled",
	# 	"intra_locus_dist_avg_cosine",
	# 	"intra_locus_dist_avg_pearsons"
	# )
    
    # For CERENKOV2 revision, 10 types of intralocus distances were used
    # 
    # p_dist_col_names <- c(
	# 	"intra_locus_dist_avg_canberra",
	# 	"intra_locus_dist_avg_canberra_scaled",
	# 	"intra_locus_dist_avg_euclidean",
	# 	"intra_locus_dist_avg_euclidean_scaled",
	# 	"intra_locus_dist_avg_manhattan",
	# 	"intra_locus_dist_avg_manhattan_scaled",
	# 	"intra_locus_dist_avg_cosine",
    #   "intra_locus_dist_avg_cosine_scaled",
	# 	"intra_locus_dist_avg_pearsons",
	# 	"intra_locus_dist_avg_pearsons_scaled"
	# )
    
    stopifnot( ! is.null(p_dist_col_names) )
	stopifnot( ! is.null(p_num_bins) )
	
	if (! suppressMessages( require(dplyr, quietly=TRUE) ) ) { stop("package dplyr is missing") }
	
	# Assume the input feature matrix is actually a data frame
	p_input_feature_matrix_as_df <- p_input_feature_matrix
	p_input_feature_matrix_as_df$label <- p_case_label_vec
	
	inds_cases_train <- setdiff(1:length(p_case_label_vec), p_inds_cases_test)
	
	train_df <- p_input_feature_matrix_as_df[inds_cases_train, ]
	train_R_df <- filter(train_df, label == 1)
	train_C_df <- filter(train_df, label == 0)
	
	llr_lists <- lapply(p_dist_col_names, function(dist_col) {
		## ----- Determine break points on ALL data ----- ##
		
		dist_vec <- na.omit(train_df[, dist_col])
		breaks <- seq(min(dist_vec), max(dist_vec), length.out = p_num_bins + 1)  # n bins, n+1 break points
		breaks[1] <- -Inf
		breaks[p_num_bins + 1] <- Inf
		
		## ----- Calculate likelihood for each bin on TRAINING data ----- ##
		
		train_R_dist_vec <- na.omit(train_R_df[, dist_col])
		train_C_dist_vec <- na.omit(train_C_df[, dist_col])
		
		R_hist <- hist(train_R_dist_vec, breaks = breaks, plot = FALSE)
		C_hist <- hist(train_C_dist_vec, breaks = breaks, plot = FALSE)
		
		if (any(R_hist$counts == 0)) {
			print(R_hist$counts)
			stop(paste("RSNPs'", dist_col, "has empty bins", "( p_num_bins =", p_num_bins, ")"))
		}
		
		if (any(C_hist$counts == 0)) {
			print(C_hist$counts)
			stop(paste("CSNPs'", dist_col, "has empty bins", "( p_num_bins =", p_num_bins, ")"))
		}
		
		bin_posterior_ratio <- R_hist$counts / C_hist$counts
		prior_ratio <- length(train_R_dist_vec) / length(train_C_dist_vec)
		bin_likelihood_ratio <- bin_posterior_ratio / prior_ratio
		
		## ----- Fit ALL data to bins ----- ## 
		
		# NaN and NA elements of x are mapped to NA codes, as are values outside range of breaks.
		# Start from 1
		snp_bin_id <- .bincode(p_input_feature_matrix_as_df[, dist_col], breaks = breaks, right = TRUE, include.lowest = TRUE)
		
		## ----- Calculate Log Likelihood for ALL Data ------ ##
		
		snp_likelihood_ratio <- bin_likelihood_ratio[snp_bin_id]
		snp_log_likelihood_ratio <- log(snp_likelihood_ratio)
		snp_log_likelihood_ratio[is.na(snp_log_likelihood_ratio)] <- 0  # fillna(0)
		
		## ----- Wrap return values ----- ##
		llr_name <- gsub("avg", "llr", dist_col)
		ret <- list()
		ret[[llr_name]] <- snp_log_likelihood_ratio
		
		return(ret)
	})
	
	llr_df <- data.frame(llr_lists)
	rownames(llr_df) <- rownames(p_input_feature_matrix_as_df)
	
	return(as.matrix(llr_df))
}

g_feature_reducer_fit_llr <- function(p_input_feature_matrix,
								  p_case_label_vec,
								  p_inds_cases_test,
                                  p_distribution_configs) {
    # For CERENKOV2 first submission, 8 types of intralocus distances were used.
    # And for "cosine" and "pearsons", Weibull distributions were fitted
    # 
    # p_distribution_configs <- list(
	# 	intra_locus_dist_avg_canberra = "lnorm",
	# 	intra_locus_dist_avg_canberra_scaled = "lnorm",
	# 	intra_locus_dist_avg_euclidean = "lnorm",
	# 	intra_locus_dist_avg_euclidean_scaled = "lnorm",
	# 	intra_locus_dist_avg_manhattan = "lnorm",
	# 	intra_locus_dist_avg_manhattan_scaled = "lnorm",
	# 	intra_locus_dist_avg_cosine = "weibull",
	# 	intra_locus_dist_avg_pearsons = "weibull"
	# )

    # For CERENKOV2 revision, 10 types of intralocus distances were used.
    # And for "cosine" and "pearsons", normal distributions were fitted 
    # (because AICs were not correctly calculated previously)
    # 
    # p_distribution_configs <- list(
	# 	intra_locus_dist_avg_canberra = "lnorm",
	# 	intra_locus_dist_avg_canberra_scaled = "lnorm",
	# 	intra_locus_dist_avg_euclidean = "lnorm",
	# 	intra_locus_dist_avg_euclidean_scaled = "lnorm",
	# 	intra_locus_dist_avg_manhattan = "lnorm",
	# 	intra_locus_dist_avg_manhattan_scaled = "lnorm",
	# 	intra_locus_dist_avg_cosine = "norm",
    #   intra_locus_dist_avg_cosine_scaled = "lnorm",
	# 	intra_locus_dist_avg_pearsons = "norm",
    #   intra_locus_dist_avg_pearsons_scaled = "lnorm"
	# )

    stopifnot( ! is.null(p_distribution_configs) )

	if (! suppressMessages( require(dplyr, quietly=TRUE) ) ) { stop("package dplyr is missing") }
	if (! suppressMessages( require(fitdistrplus, quietly=TRUE) ) ) { stop("package fitdistrplus is missing") }
	
	# Assume the input feature matrix is actually a data frame
	p_input_feature_matrix_as_df <- p_input_feature_matrix
	# There is already a column for label
	# p_input_feature_matrix_as_df$label <- p_case_label_vec
	
	inds_cases_train <- setdiff(1:length(p_case_label_vec), p_inds_cases_test)
	
	train_df <- p_input_feature_matrix_as_df[inds_cases_train, ]
	train_R_df <- filter(train_df, label == 1)
	train_C_df <- filter(train_df, label == 0)
	
	llr_lists <- lapply(names(p_distribution_configs), function(dist_col) {
		## ----- Parse Config ----- ##
		
		# `dist_col` is the name of one of distance columns, 
		# such as "intra_locus_dist_avg_canberra"
		distribution_name <- p_distribution_configs[[dist_col]]
		density_func_name <- paste0("d", distribution_name)
		
		## ----- Fit Distributions on TRAINING data ----- ##
		
		train_R_dist_vec <- as.vector(na.omit(train_R_df[, dist_col]))
		train_C_dist_vec <- as.vector(na.omit(train_C_df[, dist_col]))
		
		train_R_fit <- fitdist(train_R_dist_vec, distribution_name)
		train_C_fit <- fitdist(train_C_dist_vec, distribution_name)
		
		## ----- Calculate Densities for ALL data ----- ## 
		
		# do.call can find a function by its name in string
		
		R_density_args <- as.list(train_R_fit$estimate)
		R_density_args$x <- p_input_feature_matrix_as_df[, dist_col]
		R_density_vec <- do.call(density_func_name, R_density_args)
		
		C_density_args <- as.list(train_C_fit$estimate)
		C_density_args$x <- p_input_feature_matrix_as_df[, dist_col]
		C_density_vec <- do.call(density_func_name, C_density_args)
		
		## ----- Calculate Log Likelihood for ALL Data ------ ##
		
		snp_likelihood_ratio <- R_density_vec / C_density_vec
		snp_log_likelihood_ratio <- log(snp_likelihood_ratio)
		snp_log_likelihood_ratio[is.na(snp_log_likelihood_ratio)] <- 0  # fillna(0)
		
		## ----- Wrap return values ----- ##
		llr_name <- gsub("avg", "llr", dist_col)
		ret <- list()
		ret[[llr_name]] <- snp_log_likelihood_ratio
		
		return(ret)
	})
	
	llr_df <- data.frame(llr_lists)
	rownames(llr_df) <- rownames(p_input_feature_matrix_as_df)
	
	return(as.matrix(llr_df))
}

#' Grid-expands lists of possible values for each hyperparameter, into tuples
#'
#' This function accepts a list of lists (one list for each type of
#' hyperparameter). The inner list contains possible values for the
#' hyperparameter. The function grid-expands the hyperparameter values
#' and returns a list of hyperparameter tuple lists.
#'
#' @param p_param_list_values A named list of lists; names are hyperparameter
#'     names, and values are lists of possible values for the hyperparameters.
#' @return A list of hyperparameter tuple lists, obtained by grid-expanding the
#'     list of lists of hyperparameter values passed in
#'     \code{p_param_list_values}.
#' @seealso
#' @rdname g_make_hyperparameter_grid_list
#' @export
#' @importFrom
#' @author Stephen A. Ramsey, \email{stephen.ramsey@@oregonstate.edu}
g_make_hyperparameter_grid_list <- function(p_param_list_values) {
    hyperparams_df <- expand.grid(p_param_list_values, stringsAsFactors=FALSE)
    hyperparams_list <- lapply(setNames(split(hyperparams_df, seq(nrow(hyperparams_df))), NULL), as.list)
    hyperparams_list <- lapply(hyperparams_list, function(p_list) { attr(p_list, "out.attrs") <- NULL; p_list })
}

#' Returns \code{TRUE} if the argument is a valid feature matrix, or
#' \code{FALSE} if it is not.
#'
#' This function makes sure that the data frame (or matrix) contains no NaN
#' values and no factors with more than 64 levels.  Returns \code{TRUE} if the
#' data frame argument is a valid feature matrix, or \code{FALSE} if it is
#' not. Validity depends on the data frame passing all of the following tests:
#' (1) no categorical (factor) column can contain more than 64 levels; and (2)
#' no NaN values.
#'
#' @param p_feature_matrix A data frame or matrix containing the feature
#'     data. Each row corresponds to a case, and each column corresponds to a
#'     feature.
#' @return \code{TRUE} if the feature matrix is valid; \code{FALSE} if it is
#'     not.
#' @seealso
#' @rdname g_feature_matrix_is_OK
#' @export
#' @importFrom
#' @author Stephen A. Ramsey, \email{stephen.ramsey@@oregonstate.edu}
g_feature_matrix_is_OK <- function(p_feature_matrix) {
    if ("matrix" %in% methods::is(p_feature_matrix)) {
        all(! is.nan(p_feature_matrix)) &
            all(apply(p_feature_matrix, 2,
                      function(mycol) {
                          if("integer" %in% methods::is(mycol)) {
                              length(unique(mycol)) <= 64
                          } else {
                              TRUE
                          }}))
    } else {
        all(sapply(p_feature_matrix, function(mycol) { all(! is.nan(mycol)) })) &
        all(sapply(p_feature_matrix, function(mycol) {
            col_is_ok <- TRUE
            if (is.factor(mycol)) {
                if(length(levels(mycol)) > 64) {
                    col_is_ok <- FALSE
                }
            }
            col_is_ok
        }))
    }
}

#' Return string locus IDs for each of a group of SNPs
#'
#' Given chromosomal coordinate information about a set of SNPs and a
#' user-specified inter-SNP distance cutoff in base pairs (default 50,000),
#' assigns each SNP to a locus, numbers the loci, and returns a character vector
#' (of length equal to the number of SNPs) containing the string lables of the
#' loci to which the SNPs are assigned (each SNP is assigned to one locus).
#'
#' @param p_snp_locus_coords A data frame whose number of rows is equal to the
#' number of SNPs, and whose row names are the SNP IDs, and that contains
#' two columns.  First column is named "chrom", which contains the chromosome
#' IDs as strings, and the second column is named "coord", which contains the
#' SNP coordinates as integers.
#' @param p_coord_distance_cutoff_bp A positive integer specifying the inter-SNP
#' distance cutoff for assigning SNPs to loci. For each chromosome, the nearest-neighbor
#' inter-SNP distances are calculated, and the chromosome-specific SNPs
#' are partitioned into loci wherever the inter-SNP distance exceeds the cutoff.
#' @return A character vector containing the string locus IDs for the SNPs.
#' @seealso
#' @rdname g_get_snp_locus_ids
#' @export
#' @importFrom
#' @author Stephen A. Ramsey, \email{stephen.ramsey@@oregonstate.edu}
g_get_snp_locus_ids <- function(p_snp_locus_coords, 
                                p_coord_distance_cutoff_bp=50000) {
    snp_chroms <- as.character(p_snp_locus_coords$chrom)
    snp_chroms_unique <- unique(snp_chroms)
    
    unlist(lapply(snp_chroms_unique,
           function(p_snp_chrom) {
               snp_coord_data_chrom <- subset(p_snp_locus_coords, chrom==p_snp_chrom)
               snp_coord_data_chrom <- snp_coord_data_chrom[order(snp_coord_data_chrom$coord),]
               intra_snp_distances <- diff(snp_coord_data_chrom$coord)
               intra_snp_distances_thresh <- as.integer(intra_snp_distances > p_coord_distance_cutoff_bp)
               locus_vals <- cumsum(c(1, intra_snp_distances_thresh))
               locus_names <- paste(p_snp_chrom, locus_vals, sep="-")
               names(locus_names) <- rownames(snp_coord_data_chrom)
               locus_names
           }))[rownames(p_snp_locus_coords)]
}


#' Assign cases to folds using stratified sampling
#'
#' Returns an integer vector (of length \code{N} equal to the number of cases)
#' in which each entry is the cross-validation fold assigment of the corresponding
#' case. 
#'
#' @param p_num_folds A positive integer specifying the number of
#'     cross-validation folds
#' @param p_case_label_vec An integer vector (of length \code{N}) containing the
#'     case labels (0 = negative case, 1 = positive case)
#' @return An integer vector of length \code{N} whose values are all in the
#'     range \code{1:p_num_folds}, giving the fold assignment for each case
#' @seealso 
#'  \code{\link[dismo]{kfold}}
#' @rdname g_assign_cases_to_folds_by_case
#' @export 
#' @importFrom dismo kfold
#' @author Stephen A. Ramsey, \email{stephen.ramsey@@oregonstate.edu}
g_assign_cases_to_folds_by_case <- function(p_num_folds,
                                            p_case_label_vec) {
    
    stopifnot(g_check_cross_validation_num_folds_not_too_big(p_num_folds, p_case_label_vec))
    
    num_cases <- length(p_case_label_vec)
    dismo::kfold(1:num_cases,
                 k=p_num_folds,
                 by=p_case_label_vec)
}

#' Assign cases to folds by group.
#'
#' Returns a function that will generate a vector of integer fold assignments
#' for cases to cross-validation folds, while maintaining class balance in the
#' folds and ensuring that all cases from a group are assigned to the same fold
#' together.
#'
#' @param p_case_group_ids Length \code{N} vector of group IDs (which can be
#'     character or integer).
#' @param p_slop_allowed A numeric scalar indicating the maximum amount by which
#'     the number of positive cases in a group can exceed the expected number
#'     (0.5 would mean that the actual number of positive cases in group can
#'     never exceed 1.5-fold the expected number). Default: 0.5.
#' @return A function with calling signature (\code{p_num_folds},
#'     \code{p_case_label_vec}). That function returns an integer vector of
#'     length \code{N} whose values are in the range \code{1:p_num_folds}, where
#'     \code{p_num_folds} is the number of folds for the cross-validation (i.e.,
#'     in ten-fold cross-validation, \code{p_num_folds=10}.
#' @author Stephen A. Ramsey, \email{stephen.ramsey@@oregonstate.edu}
g_make_assign_cases_to_folds_by_group <- function(p_case_group_ids, p_slop_allowed=0.5) {

    p_case_group_ids ## force R to evaluate the promise argument, so it is stored with the closure

    #' Return true/false depending on whether the number of CV folds is too large
    #'
    #' Returns true if the number of CV folds (passed as an argument) is not too
    #' large for the dataset, or false if the number of CV folds is too large for
    #' the dataset (i.e., greater than the number of positive cases, so stratified
    #' CV is not possible).  The purpose of this function is just to be used as a
    #' sanity check against user misconfiguration error (i.e., user switches the
    #' number of desired replications with the number of desired CV folds, in setting
    #' configuration parameter list for CERENKOV.
    #'
    #' @param p_num_folds A positive integer indicating the number of folds for the
    #'     CV
    #' @param p_case_label_vec A binary {0,1} integer vector (of length equal to the
    #'     number of cases) containing the class labels of the cases.  # @return
    #'     \code{TRUE} if the number of CV folds (passed as an argument) is not too
    #'     large for the dataset, or \code{FALSE} if the number of CV folds is too
    #'     large for the dataset (i.e., greater than the number of positive cases,
    #'     so stratified CV is not possible).
    #' @return \code{TRUE} if the number of CV folds (passed as an argument) is not
    #'     too large for the dataset, or \code{FALSE} if the number of CV folds is
    #'     too large for the dataset
    #' @seealso
    #' @rdname g_check_cross_validation_num_folds_not_too_big
    #' @export
    #' @importFrom
    #' @author Stephen A. Ramsey, \email{stephen.ramsey@@oregonstate.edu}
    check_cross_validation_num_folds_not_too_big <- function(p_num_folds, p_case_label_vec) {
        npos <- length(which(p_case_label_vec==1))
        nneg <- length(which(p_case_label_vec==0))
        p_num_folds < length(p_case_label_vec)/(nneg/npos + 1)
    }
    
    function(p_num_folds, p_case_label_vec) {

        stopifnot(check_cross_validation_num_folds_not_too_big(p_num_folds, p_case_label_vec))
    
        table_res <- table(p_case_group_ids, p_case_label_vec)
        group_pos_counts <- table_res[,2]
        stopifnot(group_pos_counts > 0)  ## there should be an rSNP in every group!
        ncase <- length(p_case_label_vec)
        group_counts <- apply(table_res, 1, sum)
        group_names <- rownames(table_res)
        ngroup <- length(group_names)
        group_fold_assignments <- setNames(rep(NA, ngroup), group_names)
        case_counts_for_folds <- rep(0, p_num_folds)
        pos_case_counts_for_folds <- rep(0, p_num_folds)
        fold_ids <- 1:p_num_folds
        max_fold_case_count <- ceiling((1 + p_slop_allowed) * ncase / p_num_folds)
        num_pos_cases <- length(which(p_case_label_vec==1))
        max_fold_pos_cases <- ceiling((1 + p_slop_allowed) * max_fold_case_count * num_pos_cases / ncase)
        for (i in order(group_counts, decreasing=TRUE)) {
            group_count <- group_counts[i]
            group_pos_count <- group_pos_counts[i]
            inds_allowed <- which(case_counts_for_folds + group_count <= max_fold_case_count &
                                  pos_case_counts_for_folds + group_pos_count <= max_fold_pos_cases)
            stopifnot(length(inds_allowed) > 0)
            fold_assignment <- ifelse(length(inds_allowed) > 1,
                                      sample(inds_allowed,
                                             size=1,
                                             prob=(1 - (case_counts_for_folds[inds_allowed] / max_fold_case_count)) *
                                                 (1 - (pos_case_counts_for_folds[inds_allowed] / max_fold_pos_cases))),
                                      inds_allowed)
            group_fold_assignments[i] <- fold_assignment
            case_counts_for_folds[fold_assignment] <- case_counts_for_folds[fold_assignment] + group_counts[i]
            pos_case_counts_for_folds[fold_assignment] <- pos_case_counts_for_folds[fold_assignment] + group_pos_count
        }
        ret_case_fold_assignments <- setNames(group_fold_assignments[p_case_group_ids], NULL)
        stopifnot(max(ret_case_fold_assignments)==p_num_folds)
        stopifnot(min(ret_case_fold_assignments)==1)
        stopifnot(all(! is.nan(ret_case_fold_assignments)))
        stopifnot(all(! is.na(ret_case_fold_assignments)))
        stopifnot(all(is.integer(ret_case_fold_assignments)))
        ret_case_fold_assignments
    }
}

#' Make a function that runs a classification task within an error-catching context.
#'
#' Runs a user-supplied function in a tryCatch block, and if there is a warning
#' or error, passes the warning/error text to a user-specified message
#' notification function.
#'
#' @param p_classifier_runner_func A user-supplied classification task function
#'     (that is called with no arguments)
#' @param p_send_message_notification A user-specified function for error
#'     handling (which prints the error/warning text and, depending on
#'     configuration, optionally sends the error message text in an SMS message)
#' @return The wrapper function that (when called) will run the classification
#'     task within an error-catching context.
#' @author Stephen A. Ramsey, \email{stephen.ramsey@@oregonstate.edu}
g_make_classifier_runner_func_err_handling <- function(p_classifier_runner_func,
                                                       p_send_message_notification) {
    function() {
        tryCatch( { p_classifier_runner_func() },
                 warning=function(w) { p_send_message_notification(w); NULL },
                 error=function(e) { p_send_message_notification(e); NULL })
    }
}

#' Verify that singleton groups contain only positive cases
#'
#' Checks to make sure that any singleton case groups (i.e., groups with only
#' one case belonging to each of them) are only associated with cases that are
#' positive class labels.
#'
#' @param p_labels Integer vector of length equal to the number of cases. In
#'     each vector entry, value 0 means that the case is a negative case, and
#'     value 1 means that the case is a positive case.
#' @param p_group_to_case_map_list a list of length \code{unique(p_case_groups)}
#'     (see \code{\link{g_make_gorup_to_case_ids_map_list}}) in which each
#'     element (corresponding to a single group) contains an integer vector of
#'     case IDs of the cases that belong to that group.
#' @return \code{TRUE} or \code{FALSE} indicating whether or not all singleton
#'     groups contain only positive cases.
#' @author Stephen A. Ramsey, \email{stephen.ramsey@@oregonstate.edu}
g_verify_that_singleton_groups_are_all_positive_cases <-
    function(p_labels,
             p_group_to_case_map_list) {
        num_cases_per_group <- sapply(p_group_to_case_map_list, length)
        all(p_labels[unlist(p_group_to_case_map_list[names(num_cases_per_group[num_cases_per_group==1])])]==1)
    }

#' Make the group-to-case-ids mapping list
#' 
#' From a vector containing the group IDs of a set of N cases (case IDs are
#' presumed numbered 1:N), return a list (of length equal to the number of
#' unique group IDs) mapping group ID to the case IDs of the cases that are
#' associated with the group ID.
#'
#' @param p_case_groups Length N vector of group IDs (which can be character or
#'     integer).
#' @return a list of length \code{unique(p_case_groups)} in which each element
#'     (corresponding to a single group) contains an integer vector of case IDs
#'     of the cases that belong to that group.
#' @author Stephen A. Ramsey, \email{stephen.ramsey@@oregonstate.edu}
g_make_group_to_case_ids_map_list <- function(p_case_groups) {
    unique_group_ids <- unique(p_case_groups)
    setNames(lapply(unique_group_ids, function(p_group_id) {
        which(p_case_groups == p_group_id)
    }), unique_group_ids)
}

g_get_args_after_hyphen_hyphen_args <- function(p_rscript_args) {
    ind_match <- which("--args" == p_rscript_args)
    if (length(ind_match) > 0) {
        p_rscript_args[(ind_match+1):length(p_rscript_args)]
    } else {
        NA
    }
}

g_get_script_name_from_rscript_args <- function(p_rscript_args) {
    ind_match <- grep("--file=", p_rscript_args)
    stopifnot(length(ind_match) > 0)
    gsub("--file=", "", p_rscript_args[ind_match])
}

#g_get_temp_dir <- function() {
#    file.path(tempdir())
#}
