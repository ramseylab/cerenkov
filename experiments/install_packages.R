install.packages(c("plyr", "parallel", "pbapply", "reshape2", "ggplot2", "methods", "aws.ec2", "ranger", "xgboost", "pls", "dplyr", "fitdistrplus", "extrafont", "arrangements"))

# Need "ranger" version 0.6.0 specifically
# Rager Archive: https://cran.r-project.org/src/contrib/Archive/ranger/
ranger_package_url <- "https://cran.r-project.org/src/contrib/Archive/ranger/ranger_0.6.0.tar.gz"
install.packages(ranger_package_url, repos=NULL, type="source")