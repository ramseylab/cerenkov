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
library(dplyr)

g_args <- commandArgs(trailingOnly=TRUE)
g_input_file <- g_args[1]
print(sprintf("reading file: %s", g_input_file))
print(sprintf("g_args: %s", g_args))

g_output_file_base_name <- strsplit(g_input_file, ".Rdata")[[1]][1]

if (! is.na(g_args[2])) {
	if (! dir.exists(g_args[2])) {
		dir.create(g_args[2])
	}
    g_output_file_base_name <- paste(g_args[2], "/", g_output_file_base_name, sep="")
    
    g_plot_pdf_name <- paste(g_args[2], "/", "plot_feat_impt.pdf", sep="")
    g_plot_png_name <- paste(g_args[2], "/", "plot_feat_impt.png", sep="")
    g_plot_embed_pdf_name <- paste(g_args[2], "/", "plot_feat_impt_embed.pdf", sep="")
} else {
    g_output_file_base_name <- ""
    
    g_plot_pdf_name <- "plot_feat_impt.pdf"
    g_plot_png_name <- "plot_feat_impt.png"
    g_plot_embed_pdf_name <- "plot_feat_impt_embed.pdf"
}

# loading `g_ml_results` and `g_par`
load(g_input_file)


## ----- Extract impurity and permutation scores from Ranger results ----- ##

# g_ml_results[["feature_impt_scores"]] has only one element, whose key is actually "g_ranger_hp_grid"
# And "g_ranger_hp_grid" is from `classifier_hyperparameter_set_type_name` in `cerenkov_script_run_ranger.R`

# Key "impurity" and "permutation" are from `classifier_set_name` in `cerenkov_script_run_ranger.R`

# TODO: Decouple the logic of key determination

impurity_df <- data.frame(g_ml_results[["feature_impt_scores"]][[1]][["impurity"]][[1]])
impurity_df <- impurity_df[order(-impurity_df[[1]]), , drop=FALSE]

output_impurity_file_name <- paste(g_output_file_base_name, "_impurity.txt", sep="")
write.table(impurity_df, file=output_impurity_file_name, quote=FALSE, sep='\t', col.names=FALSE)

permutation_df <- data.frame(g_ml_results[["feature_impt_scores"]][[1]][["permutation"]][[1]])
permutation_df <- permutation_df[order(-permutation_df[[1]]), , drop=FALSE]

output_permutation_file_name <- paste(g_output_file_base_name, "_permutation.txt", sep="")
write.table(permutation_df, file=output_permutation_file_name, quote=FALSE, sep='\t', col.names=FALSE)


## ----- Plot Feature Importance (Figure 5) ----- ##

feat_categories <- read.table("osu18_feature_names_categories.txt",
							  header=FALSE,
							  row.names=1,
							  sep="\t",
							  stringsAsFactors=TRUE,
							  comment.char="",
							  quote="")

library(ggplot2)
df_for_plot <- merge(data.frame(permutation=permutation_df,
								impurity=impurity_df),
					 feat_categories,
					 by.x=0, by.y=0)

df_for_plot <- setNames(df_for_plot,
						c("feature", "permutation", "impurity", "categ"))

df_for_plot$categ <- factor(df_for_plot$categ,
							levels=c("LLR",
									 "repliseq",
									 "geneannot",
									 "epigenome",
									 "featdist",
									 "chrom",
									 "eigen",
									 "phylogenetic",
									 "allelism",
									 "DHS",
									 "DNAcontent",
									 "eQTL",
									 "repeats",
									 "TFBS"))

library(extrafont)
loadfonts()

fancy_scientific <- function(l) {
	# turn in to character string in scientific notation
	l <- format(l, scientific = TRUE)
	# quote the part before the exponent to keep all the digits
	l <- gsub("^(.*)e", "'\\1'e", l)
	# turn the 'e+' into plotmath format
	l <- gsub("e", "%*%10^", l)
	# return this as an expression
	parse(text=l)
}

p <- ggplot(data=df_for_plot) +
	theme_classic(base_size=20) +
	theme(legend.title = element_blank()) +
	geom_point(aes(y=permutation, x=impurity, colour=categ, size=categ), alpha=0.8, stroke=0) +
	scale_x_log10(limits=c(0.8, 200)) +
	scale_y_log10(limits=c(3e-6, 1e-2), labels=fancy_scientific) +
	scale_colour_manual(values=c("LLR"="black", 
								 "repliseq"="red",
								 "geneannot"="blue",
								 "epigenome"="orange",
								 "featdist"="green",
								 "chrom"="gold",
								 "eigen"="chocolate",
								 "phylogenetic"="brown",
								 "allelism"="cyan",
								 "DHS"="magenta",
								 "DNAcontent"="tomato",
								 "eQTL"="khaki",
								 "repeats"="darkgray",
								 "TFBS"="lightgray")) +
	scale_size_manual(values=c("LLR"=3, 
							   "repliseq"=1,
							   "geneannot"=1,
							   "epigenome"=1,
							   "featdist"=1,
							   "chrom"=1,
							   "eigen"=1,
							   "phylogenetic"=1,
							   "allelism"=1,
							   "DHS"=1,
							   "DNAcontent"=1,
							   "eQTL"=1,
							   "repeats"=1,
							   "TFBS"=1)) + 
	guides(colour = guide_legend(override.aes = list(size=5)))

ggsave(g_plot_pdf_name, p, width=6, height=4)
ggsave(g_plot_png_name, p, width=6, height=4)

embed_fonts(g_plot_pdf_name, outfile=g_plot_embed_pdf_name)