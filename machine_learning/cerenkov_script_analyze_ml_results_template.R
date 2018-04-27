## ================================================================================
## This file is part of the CERENKOV (Computational Elucidation of the
## REgulatory NonKOding Variome) software program. CERENKOV is subject to terms
## and conditions defined in the file "LICENSE.txt", which is part of this
## CERENKOV software distribution.
##
## Title       : cerenkov_script_analyze_ml_results_template.R
## 
## Description : Script for analyzing CERENKOV machine-learning results.
##
##               Reads Rdata input file specified as argument 1 on the cmdline
##
## Usage       : cerenkov_script_analyze_ml_results.R INPUT.Rdata [OUTPUT_DIR/]
##
## Note        : For OUTPUT_DIR/, don't include a file prefix as that will be
##               obtained from the input file name (by removing the ".Rdata" file
##               suffix)
##
## Requires    : reshape2, ggplot2
##
## Author      : Stephen A. Ramsey, Oregon State University
##               https://github.com/saramsey
## ================================================================================

library(reshape2)
library(ggplot2)

g_args <- commandArgs(trailingOnly=TRUE)
g_input_file <- g_args[1]
print(sprintf("reading file: %s", g_input_file))
print(sprintf("g_args: %s", g_args))

g_output_file_base_name <- strsplit(g_input_file, ".Rdata")[[1]][1]

if (! is.na(g_args[2])) {
    g_output_file_base_name <- paste(g_args[2], "/", g_output_file_base_name, sep="")
} else {
    g_output_file_base_name <- ""
}

load(g_input_file)

g_logit <- function(x) {
    stopifnot(x <= 1.0 & x >= 0.0)
    log(x/(1.0 - x))
}

g_inv_logit <- function(x) {
    1.0/(1.0 + exp(-x))
}

g_ml_performance_results <- g_ml_results$performance_results

hyperparameter_set_type_names <- names(g_ml_performance_results)

## for each hyperparameter set type (HPST) name, make one super-wide data frame for each classifier; put on a per-HPST list
g_data_wide <- lapply(g_ml_performance_results, function(df) {
    df$classifier_id <- as.integer(as.character(df$classifier_id))
    unique_classifier_ids <- unique(df$classifier_id)

    col_classifier_name <- which(names(df) == "classifier_name")
    
    res_df <- data.frame(do.call(rbind,
                                 lapply(unique_classifier_ids,
                                        function(classifier_id) {
                                            inds <- which(df$classifier_id == classifier_id)
                                            df[inds[1], col_classifier_name:ncol(df)]})),
                         acast(df, classifier_id ~ replication_id + cv_fold_id, value.var="test_aupvr"),
                         acast(df, classifier_id ~ replication_id + cv_fold_id, value.var="test_auroc"),
                         acast(df, classifier_id ~ replication_id + cv_fold_id, value.var="train_aupvr"),
                         acast(df, classifier_id ~ replication_id + cv_fold_id, value.var="train_auroc"),
                         stringsAsFactors=FALSE)
    
    nfolds <- max(df$cv_fold_id)
    nreps <- max(df$replication_id)
    col_offset <- ncol(df) - col_classifier_name + 1
    inds <- which(df$classifier_id == unique_classifier_ids[1])
    repfolds <- paste(df$replication_id[inds], df$cv_fold_id[inds], sep="_")
    names(res_df)[(col_offset + 1):(col_offset + nfolds*nreps)] <- paste("test_aupvr", repfolds, sep="_")
    names(res_df)[(col_offset + nfolds*nreps + 1):(col_offset + 2*nfolds*nreps)] <- paste("test_auroc", repfolds, sep="_")
    names(res_df)[(col_offset + 2*nfolds*nreps + 1):(col_offset + 3*nfolds*nreps)] <- paste("train_aupvr", repfolds, sep="_")
    names(res_df)[(col_offset + 3*nfolds*nreps + 1):(col_offset + 4*nfolds*nreps)] <- paste("train_auroc", repfolds, sep="_")

    if ("train_avgrank" %in% names(df)) {
        res_df <- cbind(res_df,
                        acast(df, classifier_id ~ replication_id + cv_fold_id, value.var="test_avgrank"),
                        acast(df, classifier_id ~ replication_id + cv_fold_id, value.var="train_avgrank"))
        names(res_df)[(col_offset + 4*nfolds*nreps + 1):(col_offset + 5*nfolds*nreps)] <- paste("test_avgrank", repfolds, sep="_")
        names(res_df)[(col_offset + 5*nfolds*nreps + 1):(col_offset + 6*nfolds*nreps)] <- paste("train_avgrank", repfolds, sep="_")

    }
    
    res_df$replication_id <- NULL
    res_df$cv_fold_id <- NULL

    res_df
})

all_classifier_ids <- unique(unlist(lapply(g_data_wide, "[[", "classifier_id")))

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
                          sapply(all_classifier_ids,
                                 function(classifier_id) {
                                     unlist(lapply(g_data_wide, function(df) {
                                         ind <- which(df$classifier_id == classifier_id)
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
                                                sapply(all_classifier_ids,
                                                       function(classifier_id) {
                                                           unlist(lapply(g_data_wide, function(df) {
                                                               ind <- which(df$classifier_id == classifier_id)
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
                                                sapply(all_classifier_ids,
                                                       function(classifier_id) {
                                                           unlist(lapply(g_data_wide, function(df) {
                                                               ind <- which(df$classifier_id == classifier_id)
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
stopifnot(g_wide_df_master$classifier_id == sort(g_wide_df_master$classifier_id))
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

col_ind_classifier_id <- which(names(g_wide_df_final)=="classifier_id")

## make the classifier the first column
g_wide_df_final <- g_wide_df_final[, c(col_ind_classifier_id,
                                       setdiff(1:ncol(g_wide_df_final), col_ind_classifier_id))]

output_txt_file_name <- paste(g_output_file_base_name, ".txt", sep="")
print(sprintf("Output txt file name: %s", output_txt_file_name))

## save output to a TSV file
write.table(g_wide_df_final,
            file=output_txt_file_name,
            sep="\t",
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE)

make_results_plot <- function(p_measure_name, p_results_df) {
    cols_for_plot <- c(which("single_column_hyperparameters"==names(p_results_df)),
                       which("classifier_set_name"==names(p_results_df)),
                       setdiff(grep(paste("test", p_measure_name, sep="_"), names(p_results_df)),
                               grep(paste(p_measure_name, "_", sep=""), names(p_results_df))))

    rows_use <- which(! (p_results_df$classifier_set_name %in% c()))  # exclude here
    
    df_plot_wide <- p_results_df[rows_use, cols_for_plot]

    df_plot_melted <- data.frame(melt(df_plot_wide,
                                      id.vars=c("classifier_set_name","single_column_hyperparameters"),
                                      measure.vars=c(paste("avg_test", p_measure_name, sep="_")),
                                      value.name=toupper(p_measure_name)),
                                 lower95=melt(df_plot_wide,
                                              id.vars=c("classifier_set_name"),
                                              measure.vars=c(paste("lower_95ci__test", p_measure_name, sep="_")))$value,
                                 upper95=melt(df_plot_wide,
                                              id.vars=c("classifier_set_name"),
                                              measure.vars=c(paste("upper_95ci__test", p_measure_name, sep="_")))$value)

    row_inds_use <- which(((df_plot_melted$classifier_set_name %in% c("deltaSVM_XGB",
                                                                    "OSU_XGB",
                                                                    "RSVP_XGB")) &
                                      grepl("scale_pos_weight=1", df_plot_melted$single_column_hyperparameters)) | 
                                     ! (df_plot_melted$classifier_set_name %in% c("deltaSVM_XGB",
                                                                                "OSU_XGB",
                                                                                "RSVP_XGB")))
    df_plot_melted <- df_plot_melted[row_inds_use, ]
                                     
    classifier_set_names <- as.character(df_plot_melted$classifier_set_name)
    classifier_set_names <- gsub("_published","",classifier_set_names)

    order_decreasing <- ifelse(p_measure_name == "avgrank", TRUE, FALSE)
    
    df_plot_melted$classifier_set_name <- factor(classifier_set_names,
                                               levels=classifier_set_names[order(df_plot_melted[[toupper(p_measure_name)]],
                                                                               decreasing=order_decreasing)])

    output_file_name <- paste(g_output_file_base_name, "_", p_measure_name, "_compare_models.pdf", sep="")
    print(sprintf("Saving plot file %s", output_file_name))
                              
    ggplot(df_plot_melted, aes_string(x="classifier_set_name", y=toupper(p_measure_name))) +
        geom_point() +
        theme_gray(base_size=18) +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_text(angle=45, hjust=1)) +
        geom_errorbar(aes(ymin=lower95, ymax=upper95)) +
        ggsave(output_file_name)
}                             

if (all(! is.na(g_wide_df_final$avg_test_aupvr))) {
    make_results_plot("aupvr", g_wide_df_final)
    make_results_plot("auroc", g_wide_df_final)
    make_results_plot("avgrank", g_wide_df_final)
}

