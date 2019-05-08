g_args <- commandArgs(trailingOnly=FALSE)  ## convert to using "aargh" package, some day

source("../../machine_learning/cerenkov_func_base.R")         ## load functions used for machine-learning

g_trailing_args <- g_get_args_after_hyphen_hyphen_args(g_args)

## ============================ define global parameters =============================

g_par <- c(
    list(
        num_folds_cross_validation =      5,     ## we are standardizing on 5-fold CV
        num_cv_replications =             10,            ## set to anywhere from 1--200, typically
        flag_locus_sampling =             TRUE         ## set to false if you want SNP-level sampling
    )
)

## ============================== load OSU annotation feature data =================================
print("loading OSU data")
g_annot_feat_df <- readRDS(file="osu18_features1.1_cerenkov2.rds")
stopifnot(g_feature_matrix_is_OK(g_annot_feat_df))

g_snp_names <- rownames(g_annot_feat_df)
g_label_vec <- as.integer(as.character(g_annot_feat_df$label))

g_snp_coords_df <- readRDS("osu18_snp_coordinates.rds")
g_snp_locus_ids <- g_get_snp_locus_ids(g_snp_coords_df)

## ============================== output cv partitions =================================

g_assign_cases_to_folds_by_locus <- g_make_assign_cases_to_folds_by_group(g_snp_locus_ids)
g_locus_to_snp_ids_map_list <- g_make_group_to_case_ids_map_list(g_snp_locus_ids)

g_assign_cases_to_folds <- ifelse(g_par$flag_locus_sampling,
								   c(g_assign_cases_to_folds_by_locus),
								   c(g_assign_cases_to_folds_by_case))[[1]]

replications_fold_assignments_list <- replicate(g_par$num_cv_replications,
												g_assign_cases_to_folds(p_num_folds=g_par$num_folds_cross_validation,
																		p_case_label_vec=g_label_vec),
												simplify=FALSE)

replications_fold_assignments_matrix = matrix(unlist(replications_fold_assignments_list), 
											  ncol=g_par$num_cv_replications, 
											  byrow=FALSE)
replications_fold_assignments_df <- data.frame(replications_fold_assignments_matrix, 
				                               row.names=g_snp_names)

colnames(replications_fold_assignments_df) <- paste0("replication", 1:g_par$num_cv_replications)
replications_fold_assignments_df$name <- rownames(replications_fold_assignments_df)

write.table(replications_fold_assignments_df, file='osu18_replications_fold_assignments.tsv', quote=FALSE, sep='\t', row.names=FALSE)