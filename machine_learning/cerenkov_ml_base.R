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
##

## cerenkov_ml_base.R -- This file defines various functions used in the CERENKOV software.
##    This file is exclusively function definitions; no evaluations are performed. No function
##    in this file references any objects in the global environment. This file is loaded
##    by the various CERENKOV scripts (e.g., cerenkov_ml_compare_models.R) using the "source"
##    function. At some point this source file will be converted into an R package.
##
## Author:  Stephen Ramsey
## 
## Packages required by this script:
##   PRROC
##
## Packages conditionally required by this script:
##   xgboost, Matrix, dismo, ranger, methods
##
## Note:  do not use this program with Ranger version 0.6.6 (stability issues); use Ranger 0.6.0 only
##

## --------------- Cerenkov-specific functions start here ---------

## p_workplan_list:  a list of "workplans".  Each workplan is a list with elements:
##    workplan$classifier_feature_matrix_name: the integer ID of the feature matrix to use for this workplan
##    workplan$classifier_hyperparameter_list:  the list of hyperparameters to use for this classifier
##    workplan$classifier_function_name:  the function to be used for training the classifier
## p_classifier_functions_list:  list of functions, each with signature: function(p_classifier_feature_matrix,
##                                                                                p_classifier_hyperparameter_list,
##                                                                                p_label_vector,
##                                                                                p_inds_cases_test), returns list with four elements
##                                                                                (train_auroc, train_aupvr, test_auroc, test_aupvr)
## p_classifier_feature_matrices_list:  one feature matrix for each type of feature matrix data structure to use **NO CASE LABELS**
## p_case_label_vec:  numeric vector containing the feature labels (0 or 1 only)
## p_func_lapply_first_level:  function(p_X, p_FUNC) returning a list
## p_func_lapply_second_level:  function(p_X, p_FUNC) returning a list

g_run_mult_classifs_mult_hyperparams_cv <- function(p_workplan_list,
                                                    p_classifier_functions_list,
                                                    p_classifier_feature_matrices_list,
                                                    p_case_label_vec,
                                                    p_num_cv_replications=1,
                                                    p_num_folds=10,
                                                    p_func_lapply_first_level=lapply,
                                                    p_func_lapply_second_level=lapply,
                                                    p_feature_reducer_functions_list=NULL,
                                                    p_assign_cases_to_folds) {

     classifier_hyperparameter_type_names_unique <- sort(unique(unlist(lapply(p_workplan_list,
                                                                              function(p_workplan) {
                                                                                  p_workplan$classifier_hyperparameter_set_type_name
                                                                              }))))

    ## check if there is at least one workplan on the workplan list
    stopifnot(length(p_workplan_list) > 0)

    ## we need to know how many cases there are, in order to assign the cases to cross-validation folds
    num_cases <- unique(sapply(p_classifier_feature_matrices_list, nrow))
    if (length(num_cases) > 1) {
        stop("all classifier feature matrices must have equal numbers of cases")
    }

    replications_fold_assignments_list <- replicate(p_num_cv_replications,
                                                    p_assign_cases_to_folds(p_num_folds=p_num_folds,
                                                                            p_case_label_vec=p_case_label_vec),
                                                    simplify=FALSE)

    replications_folds_list_list <- lapply(data.frame(t(expand.grid(1:p_num_cv_replications,
                                                               1:p_num_folds))),
                                      function(mycol) {
                                          list(replication_id=mycol[1],
                                               fold_id=mycol[2])
                                      })

    ml_global_results_list <- p_func_lapply_first_level(replications_folds_list_list,
            function(replication_fold_list) {
                replication_id <- replication_fold_list$replication_id
                fold_id <- replication_fold_list$fold_id
                fold_assignments <- replications_fold_assignments_list[[replication_id]]
                inds_cases_test <- which(fold_id == fold_assignments)
                if (length(inds_cases_test) == length(fold_assignments)) {
                    ## this means we have num_folds=1, i.e., no cross-validation; train on all cases, so set test cases to empty vector
                    inds_cases_test <- c()
                }
# ------------- debugging code ------------
#                print(sprintf("starting workplans for replication %d, fold %d at date/time: %s", replication_id, fold_id, Sys.time()))
# ------------- debugging code ------------

                ml_single_results_list <- p_func_lapply_second_level(1:length(p_workplan_list),
                       function(p_workplan_list_index) {
                           ## get the workplan
                           p_workplan <- p_workplan_list[[p_workplan_list_index]]

                           ## get the integer workplan ID
                           workplan_id <- names(p_workplan_list)[p_workplan_list_index]

# ------------- debugging code ------------
###                         print(sprintf("Starting workplan ID %s", p_workplan$workplan_set_name))
# ------------- debugging code ------------

                           ## need to know the feature matrix name, so we can retrieve the feature matrix
                           classifier_feature_matrix_name <- p_workplan$classifier_feature_matrix_name

                           ## need the classifier's hyperparameter list
                           classifier_hyperparameter_list <- p_workplan$classifier_hyperparameter_list

                           ## need the classifier's hyperparameter set name, which dictates what top-level list slot the results go into
                           classifier_hyperparameter_set_type_name <- p_workplan$classifier_hyperparameter_set_type_name

                           ## need the classifier's function name so we can look up the classifier function
                           classifier_function_name <- p_workplan$classifier_function_name

                           ## need the classifier function so we can run the train/test cycle
                           classifier_function <- p_classifier_functions_list[[classifier_function_name]]
                           
                           if (is.null(p_workplan$feature_reducer_function_name)) {
                               ## this is the standard case, the workplan doesn't call for using a feature matrix reducer
                               
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
                                
                                save_hyperparameter_list <- classifier_hyperparameter_list
                            } else {

                                ## we are using a "supervised" feature reducer function, like PLS
                                base_feature_matrix <- p_classifier_feature_matrices_list[[classifier_feature_matrix_name]]
                                if (is.null(base_feature_matrix)) { stop(sprintf("feature matrix %s missing",
                                                                                 classifier_feature_matrix_name)) }
                                feature_reducer_function_name <- p_workplan$feature_reducer_function_name
                                feature_reducer_function <- p_feature_reducer_functions_list[[feature_reducer_function_name]]
                                stopifnot( ! is.null(feature_reducer_function))
                                feature_reducer_input_matrix_name <- p_workplan$feature_reducer_input_matrix_name
                                stopifnot( ! is.null(feature_reducer_input_matrix_name))
                                feature_reducer_input_matrix <- p_classifier_feature_matrices_list[[feature_reducer_input_matrix_name]]
                                if (is.null(feature_reducer_input_matrix)) { stop(sprintf("feature matrix %s missing",
                                                                                          feature_reducer_input_matrix_name)) }
                                stopifnot( ! is.null(feature_reducer_input_matrix))
                                feature_reducer_hyperparameters_list <- p_workplan$feature_reducer_hyperparameters_list
                                stopifnot( ! is.null(feature_reducer_hyperparameters_list))

                                ## call the feature reducer
                                reduced_feature_matrix <- do.call(feature_reducer_function,
                                                                  c(list(p_input_feature_matrix=feature_reducer_input_matrix,
                                                                         p_case_label_vec=p_case_label_vec,
                                                                         p_inds_cases_test=inds_cases_test),
                                                                    feature_reducer_hyperparameters_list))
                                
                                ## combine the reduced feature matrix with the base feature matrix
                                if ("sparseMatrix" %in% is(base_feature_matrix)) {
                                    if (! require(Matrix, quietly=TRUE)) { stop("package Matrix is missing") }
                                    feature_matrix <- cBind(base_feature_matrix, reduced_feature_matrix)
                                } else {
                                    feature_matrix <- cbind(base_feature_matrix, reduced_feature_matrix)
                                }
                                
                                classifier_feature_matrix_name <- paste(classifier_feature_matrix_name,
                                                                        feature_reducer_input_matrix_name, sep="_")

                                save_hyperparameter_list <- c(classifier_hyperparameter_list,
                                                              feature_reducer_hyperparameters_list)
                            } 

                            ## train/test the classifier for the specified hyperparameters
                            classifier_run_time <- system.time(
                                classifier_ret_list <- classifier_function(feature_matrix,
                                                                           classifier_hyperparameter_list,
                                                                           p_case_label_vec,
                                                                           inds_cases_test) )

                           feat_import_scores <- classifier_ret_list$feat_import_scores
                           classifier_ret_list$feat_import_scores <- NULL
                           
                           if (is.null(save_hyperparameter_list)) {
                               save_hyperparameter_list <- list(classifier_hyperparameters.="")
                           }

                           ## create a list of results
                           ## WARNING:  DO **NOT** ALTER THE ORDER OF THESE LIST ELEMENTS:
                            list(performance_results=data.frame(c(classifier_ret_list,
                                                 list(classifier_name=classifier_function_name,
                                                      classifier_feature_matrix_name=classifier_feature_matrix_name,
                                                      classifier_hyperparameter_set_type_name=classifier_hyperparameter_set_type_name,
                                                      classifier_hyperparameters=save_hyperparameter_list,
                                                      classifier_run_time=setNames(classifier_run_time[1], NULL),
                                                      workplan_set_name=p_workplan$workplan_set_name,
                                                      workplan_id=workplan_id,
                                                      replication_id=replication_id,
                                                      cv_fold_id=fold_id)),
                                               stringsAsFactors=FALSE),
                                 feat_import_scores=feat_import_scores)
                        })
            })

    ## invert res_list so that the workplan ID is the top level, and the fold ID is the second level
    ml_global_results_list_wp_top <- lapply(1:length(p_workplan_list),
                                            function(p_workplan_list_id) {
                                                setNames(lapply(ml_global_results_list,
                                                                "[[",
                                                                p_workplan_list_id), NULL)
                                            })

    ## get the hyperparameter set type name for each workplan
    classifier_hyperparameter_type_names <- sapply(ml_global_results_list_wp_top,
            function(p_res_list_for_workplan) {
                p_res_list_for_workplan[[1]]$performance_results$classifier_hyperparameter_set_type_name })
    
    ## divide workplan list into performance results data frames, organized by hyperparameter-set type name
    ml_global_results_list_hpts_top <- setNames(lapply(classifier_hyperparameter_type_names_unique,
            function(p_classifier_hyperparameter_type_name) {
                hptds_list_list <- ml_global_results_list_wp_top[which(classifier_hyperparameter_type_names == p_classifier_hyperparameter_type_name)]
                do.call(rbind, lapply(hptds_list_list, function(hptds_list) {
                    do.call(rbind, lapply(hptds_list, "[[", "performance_results"))
                }))
            }), classifier_hyperparameter_type_names_unique)

    ## extract feature importance scores, if they were gathered
    feat_impt_scores <- setNames(lapply(classifier_hyperparameter_type_names_unique,
            function(p_classifier_hyperparameter_type_name) {
                hptds_list_list <- ml_global_results_list_wp_top[which(classifier_hyperparameter_type_names == p_classifier_hyperparameter_type_name)]
                setNames(lapply(hptds_list_list, function(hptds_list) {
                    lapply(hptds_list, "[[", "feat_import_scores")
                }), sapply(hptds_list_list, function(htps_list) {
                    htps_list[[1]]$performance_results$workplan_set_name
                }))
            }), classifier_hyperparameter_type_names_unique)

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
        
        if (! require(PRROC, quietly=TRUE)) {
            stop("package PRROC is missing")
        }

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
## NOTE:  ranger version 0.6.0 appears to leak memory when "probability=TRUE"
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
             p_inds_cases_test) {
        if (! require(ranger, quietly=TRUE)) {
            stop("package ranger is missing")
        }
        if (! require(methods, quietly=TRUE)) {
            stop("package methods is missing")
        }

        stopifnot( ! is.null(p_classifier_feature_matrix))
        stopifnot(length(grep("data.frame|matrix", is(p_classifier_feature_matrix))) > 0)
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
        } else {
            test_perf_results <- list(auroc=NA, aupvr=NA)
            if (! is.null(train_perf_results$avgrank)) {
                test_perf_results <- c(test_perf_results,
                                       list(avgrank=NA))
            }
        }
        
        if (any(is.nan(test_pred_scores))) {
            print(sprintf("number of NaN values in ranger prediction scores: %d", length(which(is.nan(test_pred_scores)))))
            print("ranger hyperparameter set: ")
            print(p_classifier_hyperparameter_list)
            test_pred_scores[is.nan(test_pred_scores)] <- 0.5
        }
        
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
        
        if (! is.null(feat_import_scores)) {
            ret_list <- c(ret_list,
                          list(feat_import_scores=feat_import_scores))
        }

        ret_list
    }
}       

g_make_custom_xgboost_objective <- function(p_weight_pos_cases) {
    function(preds, dtrain) {
        labels <- getinfo(dtrain, "label")
        weights <- rep(1, length(preds))
        weights[labels==1] <- 0.5
        my_loss_function_gradient_values <- weights*(preds - labels)
        ##  calculate your loss function gradient here
        list(grad=my_loss_function_gradient_values, hess=rep(1, length(labels)))
    }
}

## function constructor for the "xgboost" classifier; this is hard-coded to use "NA" to denote missing data
g_make_classifier_function_xgboost <- function(p_nthread=1,
                                               p_get_perf_results,
                                               p_feature_importance_type=NULL,
                                               p_objective_function="binary:logistic",
                                               p_case_group_ids=NULL) {

    if (! require(xgboost, quietly=TRUE)) {
        stop("package xgboost is missing")
    }
    
    nthread_third_level_list <- list(nthread=p_nthread)

    avg_group_size_ncases <- mean(table(p_case_group_ids))
        
    function(p_classifier_feature_matrix,
             p_classifier_hyperparameter_list,
             p_label_vector,
             p_inds_cases_test) {

        ## load xgboost package; we will need it to assign cases to cross-validation folds
        if (! require(xgboost, quietly=TRUE)) {
            stop("package xgboost is missing")
        }
    
        stopifnot( ! is.null(p_classifier_feature_matrix))
        stopifnot(length(grep("dgCMatrix|matrix", is(p_classifier_feature_matrix))) > 0)
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
        
        xgb_params <- c(p_classifier_hyperparameter_list,
                       list(objective=p_objective_function))

        xgb_model <- xgb.train(params=xgb_params,
                               data=model_data,
                               nrounds=p_classifier_hyperparameter_list$nrounds,
                               verbose=0,
                               nthread=p_nthread)
       
        if (! is.null(p_feature_importance_type)) {
            dump_file_name <- tempfile("xgb_dump")
            xgb.dump(xgb_model, fname=dump_file_name, with.stats=TRUE)
            feat_import_scores <- as.data.frame(xgb.importance(feature_names=colnames(p_classifier_feature_matrix),
                                                               filename_dump=dump_file_name))
            file.remove(dump_file_name)
        } else {
            feat_import_scores <- NULL
        }
        
        train_pred_scores <- predict(xgb_model,
                                     p_classifier_feature_matrix[inds_cases_training, ],
                                     missing=NA)
        
        ## this is a work-around for bug # 1503 in xgboost (see GitHub); might have been fixed in latest release of R xgboost package
        n_train_pred_scores <- length(train_pred_scores)
        if (length(train_pred_scores) < length(train_labels)) {
            train_labels <- train_labels[1:length(train_pred_scores)]
        }
        ## this is a work-around for bug # 1503 in xgboost (see GitHub); might have been fixed in latest release of R xgboost package

        stopifnot(length(train_pred_scores) == length(train_labels))

        train_perf_results <- p_get_perf_results(train_pred_scores, train_labels, inds_cases_training)

        if (length(test_labels) > 0) {
            test_pred_scores <- predict(xgb_model,
                                        p_classifier_feature_matrix[p_inds_cases_test,],
                                        missing=NA)

            ## this is a work-around for bug # 1503 in xgboost (see GitHub); might have been fixed in latest release of R xgboost package
            n_test_pred_scores <- length(test_pred_scores)
            if (length(test_pred_scores) < length(test_labels)) {
                test_labels <- test_labels[1:length(test_pred_scores)]
            }
            ## this is a work-around for bug # 1503 in xgboost (see GitHub); might have been fixed in latest release of R xgboost package

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

g_make_hyperparameter_grid_list <- function(p_param_list_values) {
    hyperparams_df <- expand.grid(p_param_list_values, stringsAsFactors=FALSE)
    hyperparams_list <- lapply(setNames(split(hyperparams_df, seq(nrow(hyperparams_df))), NULL), as.list)
    hyperparams_list <- lapply(hyperparams_list, function(p_list) { attr(p_list, "out.attrs") <- NULL; p_list })
}

## This function makes sure that the data frame (or matrix) contains no NaN values and
## no factors with more than 64 levels.  Note:  NA values are allowed for XGboost but not RF.
g_feature_matrix_is_OK <- function(p_feature_matrix) {
    if (! require(methods, quietly=TRUE)) { stop("package methods is missing") }
    if ("matrix" %in% is(p_feature_matrix)) {
        all(! is.nan(p_feature_matrix)) &
            all(apply(p_feature_matrix, 2,
                      function(mycol) {
                          if("integer" %in% is(mycol)) {
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

g_assign_cases_to_folds_by_case <- function(p_num_folds,
                                            p_case_label_vec) {
    
    if (! suppressWarnings( require(dismo, quietly=TRUE) ) ) {
        stop("package dismo is missing")
    }

    num_cases <- length(p_case_label_vec)
    kfold(1:num_cases,
          k=p_num_folds,
          by=p_case_label_vec)
}

g_make_assign_cases_to_folds_by_group <- function(p_case_group_ids, p_slop_allowed=0.5) {

    p_case_group_ids ## force R to evaluate the promise argument, so it is stored with the closure

    function(p_num_folds, p_case_label_vec) {

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
        ret_case_fold_assignments
    }
}

## whatever function you pass as the first argument, conditionally makes it an
## "error handling" function depending on the third argument (this is because we
## don't want to use tryCatch when we are debugging, only when we are running in
## production)
g_make_classifier_runner_func_err_handling <- function(p_classifier_runner_func,
                                                       p_send_message_notification,
                                                       p_cluster) {
    
    if (! is.null(p_cluster)) {
        res_func <- function() {
            tryCatch( { p_classifier_runner_func() },
                     warning=function(w) { p_send_message_notification(w); NULL },
                     error=function(e) { p_send_message_notification(e); NULL })
        }
    } else {
        res_func <- p_classifier_runner_func
    }

    res_func
}

g_make_cluster_cleanup_function <- function(p_cluster, p_ec2_instances=NULL) {
    function() {
        if(! is.null(p_cluster)) {
            stopCluster(p_cluster)
        }

        if(! is.null(p_ec2_instances)) {
            terminate_instances(p_ec2_instances)
        }        
    }
}

