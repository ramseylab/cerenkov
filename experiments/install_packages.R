install.packages(install.packages(c("arrangements", "aws.ec2", "dplyr", "extrafont", "fitdistrplus", "ggplot2", "methods", "parallel", "pbapply", "pls", "plyr", "PRROC", "reshape2", "xgboost"))

# Need "ranger" version 0.6.0 specifically
# Rager Archive: https://cran.r-project.org/src/contrib/Archive/ranger/
ranger_package_url <- "https://cran.r-project.org/src/contrib/Archive/ranger/ranger_0.6.0.tar.gz"
install.packages(ranger_package_url, repos=NULL, type="source")