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

## cerenkov_analyze_ml_results_tune_xgboost.R -- makes a wide TSV file
## containing the results from the xgboost grid-search hyperparameter tuning
## carried out by the "cerenkov_ml_tune_xgboost.R" script.  Plots AVGRANK
## vs. AUPVR.
##
## Stephen A. Ramsey

load("cerenkov_ml_tune_xgboost.Rdata")

g_logit <- function(x) {
    stopifnot(x <= 1.0 & x >= 0.0)
    log(x/(1.0 - x))
}

g_inv_logit <- function(x) {
    1.0/(1.0 + exp(-x))
}

g_ml_performance_results <- g_ml_results$performance_results

hyperparameter_set_type_names <- names(g_ml_performance_results)

library(reshape2)

## for each hyperparameter set type (HPST) name, make one super-wide data frame for each workplan; put on a per-HPST list
g_data_wide <- lapply(g_ml_performance_results, function(df) {
    df$workplan_id <- as.integer(as.character(df$workplan_id))
    unique_workplan_ids <- unique(df$workplan_id)

    col_classifier_name <- which(names(df) == "classifier_name")
    
    res_df <- data.frame(do.call(rbind,
                                 lapply(unique_workplan_ids,
                                        function(workplan_id) {
                                            inds <- which(df$workplan_id == workplan_id)
                                            df[inds[1], col_classifier_name:ncol(df)]})),
                         acast(df, workplan_id ~ replication_id + cv_fold_id, value.var="test_aupvr"),
                         acast(df, workplan_id ~ replication_id + cv_fold_id, value.var="test_auroc"),
                         acast(df, workplan_id ~ replication_id + cv_fold_id, value.var="train_aupvr"),
                         acast(df, workplan_id ~ replication_id + cv_fold_id, value.var="train_auroc"),
                         stringsAsFactors=FALSE)
    
    nfolds <- max(df$cv_fold_id)
    nreps <- max(df$replication_id)
    col_offset <- ncol(df) - col_classifier_name + 1
    inds <- which(df$workplan_id == unique_workplan_ids[1])
    repfolds <- paste(df$replication_id[inds], df$cv_fold_id[inds], sep="_")
    names(res_df)[(col_offset + 1):(col_offset + nfolds*nreps)] <- paste("test_aupvr", repfolds, sep="_")
    names(res_df)[(col_offset + nfolds*nreps + 1):(col_offset + 2*nfolds*nreps)] <- paste("test_auroc", repfolds, sep="_")
    names(res_df)[(col_offset + 2*nfolds*nreps + 1):(col_offset + 3*nfolds*nreps)] <- paste("train_aupvr", repfolds, sep="_")
    names(res_df)[(col_offset + 3*nfolds*nreps + 1):(col_offset + 4*nfolds*nreps)] <- paste("train_auroc", repfolds, sep="_")

    if ("train_avgrank" %in% names(df)) {
        res_df <- cbind(res_df,
                        acast(df, workplan_id ~ replication_id + cv_fold_id, value.var="test_avgrank"),
                        acast(df, workplan_id ~ replication_id + cv_fold_id, value.var="train_avgrank"))
        names(res_df)[(col_offset + 4*nfolds*nreps + 1):(col_offset + 5*nfolds*nreps)] <- paste("test_avgrank", repfolds, sep="_")
        names(res_df)[(col_offset + 5*nfolds*nreps + 1):(col_offset + 6*nfolds*nreps)] <- paste("train_avgrank", repfolds, sep="_")

    }
    
    res_df$replication_id <- NULL
    res_df$cv_fold_id <- NULL

    res_df
})

all_workplan_ids <- unique(unlist(lapply(g_data_wide, "[[", "workplan_id")))

## compute the average performance values
perf_strings <- c("test_aupvr","test_auroc","train_aupvr","train_auroc")
if (length(grep("avgrank", names(g_data_wide[[1]]))) > 0) {
    perf_strings <- c(perf_strings, "test_avgrank", "train_avgrank")
}

g_logit_mean <- function(x) {
    g_inv_logit(mean(g_logit(x)))
}

avg_results <- t(do.call(rbind, lapply(perf_strings,
                      function(string_to_search_with) {
                          sapply(all_workplan_ids,
                                 function(workplan_id) {
                                     unlist(lapply(g_data_wide, function(df) {
                                         ind <- which(df$workplan_id == workplan_id)
                                         if (length(ind) == 0) { return(NULL) }
                                         values_to_average <- unlist(df[ind, grep(string_to_search_with, names(df))])
                                         if (length(grep("aupvr|auroc", string_to_search_with)) > 0) {
                                             reducer_func <- g_logit_mean
                                         } else {
                                             reducer_func <- mean
                                         }
                                         avg_val <- reducer_func(values_to_average)
                                     }))
                                 })
                      })))
colnames(avg_results) <- paste("avg", perf_strings, sep="_")

lower_95ci_results <- t(do.call(rbind, lapply(perf_strings,
                                            function(string_to_search_with) {
                                                sapply(all_workplan_ids,
                                                       function(workplan_id) {
                                                           unlist(lapply(g_data_wide, function(df) {
                                                               ind <- which(df$workplan_id == workplan_id)
                                                               if (length(ind) == 0) { return(NULL) }
                                                               values_to_analyze <- unlist(df[ind, grep(string_to_search_with, names(df))])
                                                               if (length(grep("aupvr|auroc", string_to_search_with)) > 0) {
                                                                   reducer_func <- g_logit_mean
                                                               } else {
                                                                   reducer_func <- mean
                                                               }
                                                               quantile(replicate(1000, { reducer_func(sample(values_to_analyze, replace=TRUE)) }), probs=0.025)
                                                           }))
                                                       })
                                            })))

colnames(lower_95ci_results) <- paste("lower_95ci_", perf_strings, sep="_")

upper_95ci_results <- t(do.call(rbind, lapply(perf_strings,
                                            function(string_to_search_with) {
                                                sapply(all_workplan_ids,
                                                       function(workplan_id) {
                                                           unlist(lapply(g_data_wide, function(df) {
                                                               ind <- which(df$workplan_id == workplan_id)
                                                               if (length(ind) == 0) { return(NULL) }
                                                               values_to_analyze <- unlist(df[ind, grep(string_to_search_with, names(df))])
                                                               if (length(grep("aupvr|auroc", string_to_search_with)) > 0) {
                                                                   reducer_func <- g_logit_mean
                                                               } else {
                                                                   reducer_func <- mean
                                                               }
                                                               quantile(replicate(1000, { reducer_func(sample(values_to_analyze, replace=TRUE)) }), probs=0.975)
                                                           }))
                                                       })
                                            })))

colnames(upper_95ci_results) <- paste("upper_95ci_", perf_strings, sep="_")

## make a master data frame, but not yet with the performance values
g_wide_df_master <- do.call(rbind, lapply(g_data_wide, function(df) {
    inds_hyperparameter_columns <- grep("classifier_hyperparameters", names(df))
    single_column_hyperparameters <- setNames(apply(df[, inds_hyperparameter_columns, drop=FALSE], 1, function(myrow) {
        mylist <- as.list(myrow)
        paste(paste(names(mylist), mylist, sep="="), collapse=", ")
    }), NULL)
    data.frame(df[,setdiff(1:(min(grep("aupvr", names(df)))-1),inds_hyperparameter_columns)],
               single_column_hyperparameters,
               df[,min(grep("aupvr", names(df))):ncol(df)],
               stringsAsFactors=FALSE)
}))


## Insert The performance values into the master data frame (making sure the relative row order is not messed up)
stopifnot(g_wide_df_master$workplan_id == sort(g_wide_df_master$workplan_id))
g_wide_df_final <- cbind(g_wide_df_master[,1:7],
                         avg_results,
                         lower_95ci_results,
                         upper_95ci_results,
                         g_wide_df_master[,8:ncol(g_wide_df_master)])
rownames(g_wide_df_final) <- NULL

## order the rows so that highest average AUPVR is the first row
g_wide_df_final <- g_wide_df_final[order(g_wide_df_final$avg_test_aupvr, decreasing=TRUE),]

## final underscore makes sure we don't pick up "avg_test_aupvr"!
inds_test_aupvr <- grep("test_aupvr_", names(g_wide_df_final))

g_wide_df_final$pvalue <- NA

if (nrow(g_wide_df_final) > 1) {
    g_wide_df_final$pvalue[2:nrow(g_wide_df_final)] <- sapply(2:nrow(g_wide_df_final),
                                                              function(rowind) {
                                                                  data_reference_row <- unlist(g_logit(g_wide_df_final[1, inds_test_aupvr]))
                                                                  data_other_row <- unlist(g_logit(g_wide_df_final[rowind, inds_test_aupvr]))
                                                                  ret_pv <- NA
                                                                  if (length(data_reference_row) > 1 && all(! is.na(data_reference_row))) {
                                                                      ret_pv <- t.test(data_reference_row,
                                                                                       data_other_row,
                                                                                       paired=TRUE)$p.value
                                                                  }
                                                                  ret_pv
                                                              })
}

g_wide_df_final <- g_wide_df_final[, c(1:(min(grep("aupvr", names(g_wide_df_final)))-1),
                                       ncol(g_wide_df_final),
                                       (min(grep("aupvr", names(g_wide_df_final)))):(ncol(g_wide_df_final)-1))]

col_ind_workplan_id <- which(names(g_wide_df_final)=="workplan_id")

## make the workplan the first column
g_wide_df_final <- g_wide_df_final[, c(col_ind_workplan_id,
                                       setdiff(1:ncol(g_wide_df_final), col_ind_workplan_id))]

## save output to a TSV file
write.table(g_wide_df_final,
            file="cerenkov_ml_tune_xgboost.txt",
            sep="\t",
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE)

plot_df <- g_wide_df_final[, c("avg_test_avgrank", "avg_test_aupvr")]

library(ggplot2)

ggplot(data=plot_df, aes(x=avg_test_aupvr, y=avg_test_avgrank)) +
    geom_point(alpha=0.2, stroke=0, size=1) +
    theme_gray(base_size=18) +
    xlab("AUPVR") +
    ylab("AVGRANK") +
    ggsave("avgrank_vs_aupvr_xgboost_tuning.pdf")


plot_df_zoom <- plot_df[which(plot_df$avg_test_avgrank <= 4.75 &
                              plot_df$avg_test_aupvr >= 0.4),]

ggplot(data=plot_df_zoom, aes(x=avg_test_aupvr, y=avg_test_avgrank)) +
    geom_point(alpha=0.2, stroke=0, size=1) +
    theme_gray(base_size=18) +
    xlab("AUPVR") +
    ylab("AVGRANK") +
    ggsave("avgrank_vs_aupvr_xgboost_tuning_zoom.pdf")
