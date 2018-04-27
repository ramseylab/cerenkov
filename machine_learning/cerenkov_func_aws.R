## ================================================================================
## This file is part of the CERENKOV (Computational Elucidation of the
## REgulatory NonKOding Variome) software program. CERENKOV is subject to terms
## and conditions defined in the file "LICENSE.txt", which is part of this
## CERENKOV software distribution.
##
## Title       : cerenkov_aws_base.R
## 
## Description : Define AWS-related functions for the CERENKOV project.
##
##               Each function in this file should satisfy the
##               following criteria:
##                 1. it should not not access the global environment
##                 2. side effects should be documented
##                 3. function names should start with "g_" (for "global")
##                 4. quietly check for required packages at run-time
## 
## Author      : Stephen A. Ramsey, Oregon State University
##               https://github.com/saramsey
## ================================================================================

## Eventually there should be an AWS-specific config file for CERENKOV
g_configure_aws_sns <- function() {
    list(aws_sns_topic_arn="SNS-TOPIC-ARN-GOES-HERE")  ## set to NULL if you don't want to be texted
}

## Eventually there should be an AWS-specific config file for CERENKOV
g_configure_ec2_instance_and_connection <- function() {
    list(
        sleep_sec_cluster_start = 120,  ## how many seconds to wait for a worker EC2 instance to start up, before we can ssh in
        sleep_sec_ip_addresses = 10,    ## how many seconds to wait for a secondary IP address on a worker node to become active, after it is associated with the network interface
        cluster_size_instances = 1,     ## this number is currently limited to <20 by default limits in EC2; to increase the limit, need
                                        ## to send a request their tech support, at this link:
                                        ## https://console.aws.amazon.com/support/home?region=us-west-2#/case/create?issueType=service-limit-increase&limitType=service-code-ec2-instances
        num_ip_addresses_per_instance = 1,
        ami_id = "ami-72c6050a",        ## your AMI needs to be an Ubuntu instance with R and all the required R packages installed
        instance_type = "m4.16xlarge",   
        security_group_name = "SPECIFY-AWS-SECURITY-GROUP-HERE",  ## the security group needs to allow inbound traffic on TCP ports 11000-11999
        username = "ubuntu",
        network_interface_device_name = "ens3",  ## NOTE:  this interface name is probably ubuntu-specific
        subnet_cidr_suffix = 20,        ## means it's a /20
        pem_file_name = "/LOCATION/OF/YOUR/PEMFILE.pem",        
        ssh_options="-o StrictHostKeyChecking=no -q",  ## the "-q" option prevents asking for a password interactively
        instance_tags=list(CostCenter="BUDGET-INDEX-GOES-HERE",
                           Customer="ONID-USERNAME-GOES-HERE",
                           Project="CERENKOV",
                           Name="ONIDUSERNAME-CERENKOV",
                           ManagerPasscode="PASSWORD-GOES-HERE")
    )
}

g_get_ssh_options_str <- function(p_ssh_options,
                                  p_pem_file_name) {
    pem_file_str <- if (! is.null(p_pem_file_name)) {
                        sprintf("-i %s", p_pem_file_name)
                    } else {
                        ""
                    }
    
    ssh_options_str <- if (! is.null(p_ssh_options)) {
                           p_ssh_options
                       } else {
                           ""
                       }

    sprintf("%s %s", p_ssh_options, pem_file_str)
}

## -------------------- EC2 specific functions ----------------------
g_create_ec2_instance <- function(p_image_id,
                                  p_instance_type,
                                  p_security_group_name,
                                  p_num_ip_addresses_per_instance,
                                  p_instance_tags,
                                  p_subnet_id=NULL) {

    if (! require(aws.ec2, quietly=TRUE)) {
        stop("package aws.ec2 is missing")
    }

    aws.signature::use_credentials()  ## it appears that use_credentials() has to be called for run_instances to work
                                          ## when we are running in Rscript (but in interactive R, for some reason we don't
    ## need to create a security group
    security_group <- aws.ec2::describe_sgroups(name=p_security_group_name)

    if (p_num_ip_addresses_per_instance > 1) {
        subnet_id <- if (is.null(p_subnet_id)) {
                         aws.ec2::describe_subnets()[[1]]$subnetId
                     } else {
                         p_subnet_id
                     }
        
        query <- list(Action="RunInstances",
                        ImageId=p_image_id,
                        InstanceType=p_instance_type,
                        MinCount=1,
                        MaxCount=1,
                        NetworkInterface.1.DeviceIndex=0,
                        NetworkInterface.1.SecurityGroupId=security_group$groupId,
                        NetworkInterface.1.SecondaryPrivateIpAddressCount=(p_num_ip_addresses_per_instance-1),
                        NetworkInterface.1.SubnetId=subnet_id)
        r <- aws.ec2::ec2HTTP(query)
        r <- lapply(r$instancesSet, `class<-`, "ec2_instance")

    } else {
        ## need to call use_credentials in order to be able to call run_instances). 
        r <- aws.ec2::run_instances(image=p_image_id,
                                    type=p_instance_type,
                                    sgroup=security_group,
                                    subnet=p_subnet_id,
                                    shutdown="terminate")
        print("done calling run_instances")
    }
    
    ## add tags to the instance
    aws.ec2::create_tags(r$item$instanceId, unlist(p_instance_tags))

    r
}

g_get_and_configure_ip_addresses_for_ec2_instances <- function(p_ec2_instances,
                                                               p_ec2_username,
                                                               p_ec2_network_interface_device_name,
                                                               p_ssh_bin="ssh",
                                                               p_ssh_options,
                                                               p_debug=FALSE) {
    setNames(unlist(lapply(p_ec2_instances, function(p_ec2_instance) {
        ## get the subnet for this instance
        subnet_obj <- aws.ec2::describe_subnets(p_ec2_instance$item$subnetId)

        ## get the CIDR suffix
        subnet_cidr_suffix <- strsplit(subnet_obj[[1]]$cidrBlock, "/")[[1]][2]
        
        instance_ip_address_set <- p_ec2_instance$networkInterfaceSet[[1]]$privateIpAddressesSet
        primary_private_ip_address <- instance_ip_address_set[[1]]$privateIpAddress
        ret_ip_address <- primary_private_ip_address

        if (p_debug) {
            print(p_ec2_instance$networkInterfaceSet[[1]])
        }
        
        if (length(instance_ip_address_set) > 1) {
            ## get a list of all secondary IP addresses for this EC2 instance
            ret_ip_address <- c(ret_ip_address,
                                setNames(
                                    unlist(
                                        lapply(instance_ip_address_set[2:length(instance_ip_address_set)],
                                               function(p_ip_address_item) {
                                                   ## need to turn on the secondary IP address; this is almost certainly Ubuntu-specific
                                                   secondary_ip_address <- p_ip_address_item$privateIpAddress[[1]]
                                                   system_cmd <- sprintf("%s %s@%s sudo ip addr add %s/%d dev %s",
                                                                         p_ssh_options,
                                                                         p_ec2_username,
                                                                         primary_private_ip_address,
                                                                         secondary_ip_address,
                                                                         subnet_cidr_suffix,
                                                                         p_ec2_network_interface_device_name)

                                                   if(p_debug) {
                                                       print(sprintf("running system command: %s %s", p_ssh_bin, system_cmd))
                                                   }
                                                   
                                                   if (0 != system2(p_ssh_bin, system_cmd)) {
                                                       stop(sprintf("Unable to execute system command: %s %s",
                                                                    p_ssh_bin,
                                                                    system_cmd))
                                                   }
                                                   
                                                   secondary_ip_address
                                               })),
                                    NULL))
        }

        ret_ip_address
    })), NULL)
}

g_make_create_ec2_instances_and_get_ips <- function(p_create_ec2_instance) {
    if (! require(aws.ec2, quietly=TRUE)) { stop("missing required package aws.ec2") }

    function(p_ec2_par,
             p_get_and_configure_ip_addresses_for_ec2_instances,
             p_make_cluster=FALSE,
             p_debug=FALSE) {
        
        num_instances <- if(p_make_cluster) { p_ec2_par$cluster_size_instances } else { 1 }
        stopifnot(! is.null(num_instances))

        ## here is where we create the EC2 instances that we will use
        ec2_instances <- do.call(c, lapply(1:num_instances, function(instance_number) {

            p_create_ec2_instance(p_ec2_par$ami_id,
                                  p_ec2_par$instance_type,
                                  p_ec2_par$security_group_name,
                                  p_ec2_par$num_ip_addresses_per_instance,
                                  p_ec2_par$instance_tags)
        }))

        if (p_debug) {
            print(sprintf("waiting for %d seconds for instance(s) to start", p_ec2_par$sleep_sec_cluster_start))
        }
        
        pbapply::pbsapply(1:p_ec2_par$sleep_sec_cluster_start, function(x) { Sys.sleep(1) })

        ip_addresses <- if(p_make_cluster) {
                            get_and_configure_ip_addresses_for_ec2_instances(ec2_instances,
                                                                             p_ec2_par$username,
                                                                             p_ec2_par$network_interface_device_name,
                                                                             p_ec2_par$ssh_options)
                        } else {
                            aws.ec2::describe_instances(ec2_instances$item$instanceId[[1]])[[1]]$instancesSet[[1]]$networkInterfaceSet$association$publicIp
                        }

        if (p_debug) {
            print("IP addresses for the cluster: ")
            print(ip_addresses)
        }

        ## if we have configured secondary private IP addresses for the instances, it seems sensible to wait a few seconds for the IP addresses to become active
        if (length(p_ec2_par$num_ip_addresses_per_instance) > 1) {
            print(sprintf("waiting for %d seconds for secondary IP addresses to become active", p_ec2_par$sleep_sec_ip_addresses))
            pbsapply(1:p_ec2_par$sleep_sec_ip_addresses, function(x) { Sys.sleep(1) })
        }

        cluster <-
            if (p_make_cluster) {
                print(sprintf("creating the PSOCK cluster"))
                parallel::makeCluster(ip_addresses,
                                      type="SOCK",
                                      outfile="/dev/null",
                                      rshcmd=sprintf("ssh %s", p_ec2_par$ssh_options),
                                      useXDR=TRUE,
                                      methods=FALSE)
            } else {
                NULL
            } 

        terminator_func <- function() {
            if (! is.null(cluster)) {parallel::stopCluster(cluster)}
            aws.ec2::terminate_instances(ec2_instances)
        }
        
        list(cluster=cluster,
             ec2_instances=ec2_instances,
             ip_addresses=ip_addresses,
             terminator_func=terminator_func)
    }    
} 
    
## :TODO: switch this function to use the AWS REST API rather than the AWS CLI
## Uses the AWS SNS service to text you when your machine-learning job is done
g_make_message_notifier_function <- function(p_aws_sns_topic_arn) {
    if (! is.null(p_aws_sns_topic_arn)) {
        function(p_message_text) {
            ## using ignore.stderr to suppress some deprecation warning from the AWS CLI on macOS
            system(paste("aws sns publish --topic-arn \"",
                         p_aws_sns_topic_arn,
                         "\" --message \"",
                         p_message_text,
                         "\"",
                         sep=""),
                   ignore.stdout=TRUE,
                   ignore.stderr=TRUE)
            print(p_message_text)
        }
    } else {
        print
    }    
}

## note:  look into using ssh.utils for this (but can you turn off StrictHostKeyChecking in ssh.utils?)
g_copy_files_to_aws_instance <- function(p_instance_ip_address,
                                         p_file_names,
                                         p_username="ubuntu",
                                         p_scp_bin="scp",
                                         p_ssh_options_str,
                                         p_debug=FALSE) {

    sapply(p_file_names, function(p_file_name) {
        my_cmd <- sprintf("%s %s %s@%s:",
                          p_ssh_options_str,
                          p_file_name,
                          p_username,
                          p_instance_ip_address)
        if (p_debug) {
            print(paste(p_scp_bin, my_cmd, sep=" "))
        }
        if (0 != system2(p_scp_bin,
                       my_cmd,
                       stdout=NULL,
                       stderr=NULL)) {
            stop(sprintf("Unable to execute scp command: scp %s", my_cmd))
        }
    })
}

g_copy_files_from_aws_instance <- function(p_instance_ip_address,
                                           p_file_names,
                                           p_username="ubuntu",
                                           p_scp_bin="scp",
                                           p_ssh_options_str,
                                           p_debug=FALSE) {


    sapply(p_file_names, function(p_file_name) {
        my_cmd <- sprintf("%s %s@%s:%s .",
                          p_ssh_options_str,
                          p_username,
                          p_instance_ip_address,
                          p_file_name)
        if (p_debug) {
            print(paste(p_scp_bin, my_cmd, sep=" "))
        }
        if (0 != system2(p_scp_bin,
                         my_cmd,
                         stdout=NULL,
                         stderr=NULL)) {
            stop(sprintf("Unable to execute scp command: scp %s", my_cmd))
        }
    })
}

g_run_R_script_in_aws_instance <- function(p_instance_ip_address,
                                           p_R_script_file_name,
                                           p_R_script_args_vec,
                                           p_output_logfile_name,
                                           p_username="ubuntu",
                                           p_ssh_bin="ssh",
                                           p_ssh_options_str,
                                           p_debug=FALSE) {
    
    my_cmd <- sprintf("%s %s@%s \"Rscript %s %s >%s 2>&1 &\"",
                      p_ssh_options_str,
                      p_username,
                      p_instance_ip_address,
                      p_R_script_file_name,
                      paste(p_R_script_args_vec, collapse=" "),
                      p_output_logfile_name)
    if (p_debug) {
        print(paste(p_ssh_bin, my_cmd, sep=" "))
    }
    stopifnot(0==system2(p_ssh_bin,
                         my_cmd,
                         stdout="",
                         stderr="",
                         wait=FALSE))
}

g_check_if_file_exists_in_ec2 <- function(p_instance_ip_address,
                                          p_remote_file_name,
                                          p_username="ubuntu",
                                          p_ssh_bin="ssh",
                                          p_ssh_options_str,
                                          p_debug=FALSE) {
    my_cmd <- sprintf("%s %s@%s \"[[ -f %s ]] && echo 1 || echo 0;\"",
                      p_ssh_options_str,
                      p_username,
                      p_instance_ip_address,
                      p_remote_file_name)
    if (p_debug) {
        print(paste(p_ssh_bin, my_cmd, sep=" "))
    }

    "1" == system2(p_ssh_bin,
                   my_cmd,
                   stdout=TRUE,
                   stderr="",
                   wait=FALSE)[1]
}

g_make_setup_and_run_ml_job_in_ec2 <- function(p_configure_ec2_interface,
                                               p_create_ec2_instances_and_get_ips,
                                               p_copy_files_to_aws_instance,
                                               p_run_script_in_aws_instance,
                                               p_copy_files_from_aws_instance,
                                               p_check_if_file_exists_in_ec2,
                                               p_get_ssh_options_str) {
    function(p_ml_runner_script_name,
             p_extra_files_to_push,
             p_random_seed,
             p_output_file_name,
             p_output_logfile_name,
             p_base_file_manifest,
             p_poll_interval_seconds=300,
             p_debug=FALSE) {

        ec2_par <- p_configure_ec2_interface()

        ssh_options_str <- p_get_ssh_options_str(ec2_par$ssh_options,
                                                 ec2_par$pem_file_name)
        
        ec2_instance_list <- p_create_ec2_instances_and_get_ips(ec2_par,
                                                                p_get_and_configure_ip_addresses_for_ec2_instances,
                                                                p_make_cluster=FALSE)

        func_to_do_in_try_block <- function() {
            ec2_instance_ip <- ec2_instance_list$ip_addresses

            stopifnot(length(ec2_instance_ip) == 1)
            
            file_manifest <- c(p_base_file_manifest,
                               p_ml_runner_script_name,
                               p_extra_files_to_push)
            
            p_copy_files_to_aws_instance(ec2_instance_ip,
                                         file_manifest,
                                         p_username=ec2_par$username,
                                         p_ssh_options_str=ssh_options_str,
                                         p_debug=p_debug)

            R_script_args_vec <- c(as.character(p_random_seed))
            
            p_run_script_in_aws_instance(ec2_instance_ip,
                                         p_ml_runner_script_name,
                                         R_script_args_vec,
                                         p_output_logfile_name,
                                         ec2_par$username,
                                         p_ssh_options_str=ssh_options_str,
                                         p_debug=p_debug)

            output_file_exists <- FALSE

            while(! output_file_exists) {
                Sys.sleep(p_poll_interval_seconds)
                if (p_debug) {
                    print(sprintf("Checking for output file at time: %s", Sys.time()))
                }
                
                output_file_exists <- p_check_if_file_exists_in_ec2(ec2_instance_ip,
                                                                    p_output_file_name,
                                                                    ec2_par$username,
                                                                    p_ssh_options_str=ssh_options_str,
                                                                    p_debug=p_debug)
            }
            
            p_copy_files_from_aws_instance(ec2_instance_ip,
                                           c(p_output_logfile_name, p_output_file_name),
                                           ec2_par$username,
                                           p_ssh_options_str=ssh_options_str,
                                           p_debug=p_debug)

        }

        tryCatch( { func_to_do_in_try_block() },
                 warning=function(w) { warning(w); NULL },
                 error=function(e) { warning(e); NULL })

        ec2_instance_list$terminator_func()

        TRUE
    }
}

