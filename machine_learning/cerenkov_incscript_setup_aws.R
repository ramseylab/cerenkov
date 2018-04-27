## ================================================================================
## This file is part of the CERENKOV (Computational Elucidation of the
## REgulatory NonKOding Variome) software program. CERENKOV is subject to terms
## and conditions defined in the file "LICENSE.txt", which is part of this
## CERENKOV software distribution.
##
## Title       : cerenkov_incscript_setup_aws.R
## 
## Description : Code that is sourced by a "cerenkov_script_XXXXXX.R" script in
##               order to setup an EC2 instance for running the machine-learning
##               job.
##
##               I put code in this file if it would otherwise be boilerplate in
##               *every* EC2-using script of the form
##               "cerenkov_script_XXXXXX.R".
##
##               This code assumes "cerenkov_func_base.R" and
##               "cerenkov_func_aws.R" have already been sourced.
## 
## Author      : Stephen A. Ramsey, Oregon State University
##               https://github.com/saramsey
## ================================================================================

g_cerenkov_base_file_manifest <- c("cerenkov_func_base.R",
                                   "cerenkov_func_aws.R",
                                   "cerenkov_incscript_run_ml.R",
                                   "cerenkov_incscript_setup_ml.R",
                                   "osu18_features1.1_cerenkov2.rds",
                                   "osu18_snp_coordinates.rds")

g_output_file_name <- paste(paste(g_par$output_file_base_name,
                                  g_par$analysis_label,
                                  g_par$random_number_seed,
                                  sep="_"),".Rdata",sep="")

g_script_name <- g_get_script_name_from_rscript_args(g_args)

g_create_ec2_instances_and_get_ips <- g_make_create_ec2_instances_and_get_ips(g_create_ec2_instance)

g_setup_and_run_ml_job_in_ec2 <- g_make_setup_and_run_ml_job_in_ec2(g_configure_ec2_instance_and_connection,
                                                                    g_create_ec2_instances_and_get_ips,
                                                                    g_copy_files_to_aws_instance,
                                                                    g_run_R_script_in_aws_instance,
                                                                    g_copy_files_from_aws_instance,
                                                                    g_check_if_file_exists_in_ec2,
                                                                    g_get_ssh_options_str)

