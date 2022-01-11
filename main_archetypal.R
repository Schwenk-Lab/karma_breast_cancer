## Main script for archetypal analysis pipeline of KARMA1 cohort
# Leo Dahl
# 2020-01-28


## Setup

rm(list=ls())
options(stringsAsFactors = FALSE)
starting_time <- Sys.time()
set.seed(123)

## Packages
library(rmarkdown)
library(knitr)
library(openxlsx)
library(DT)         # For the datatable function, used a lot
library(umap)
library(tsne)
library(ggplot2)
library(ggsignif)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)  # For making tableGrobs that acn be put in plotting areas
library(patchwork)  # For putting multiple ggplots in same graphing area
library(BBmisc)
library(reshape2)   # For the melt function, in multiple places
library(dplyr)
library(MDimNormn)  # For the normn_MA function, MA normalisation in normalisation.Rmd
library(archetypes) # For archetypal analysis
library(archetypal) # Another, slightly newer, package for archetypal analysis
library(gprofiler2) # For ORA using the gost function
library(table1)
library(stringr)
library(survival)   # For clogit function, condition logistic regression
library(ggrepel)    # For labels in volcano plots
library(qqman)      # For simple qqplots
library(gplots)

## Load some home-made functions
source("/home/.../misc_functions.R")
source("/home/.../plot_results.R")


### ------ Choose steps and options for analysis ------ ###

# 0 to skip and 1 to run
steps <- c(
                "get_data" = 1        # Can choose if want to use new clinical data or old
  ,        "normalisation" = 1        # Can choose different options for normalisation
  ,          "antibody_qc" = 1
  ,    "lin_normalisation" = 1        # Can choose from things to adjust for, and choose not to do logistic regression and linear mixed model
  ,        "analysis_prep" = 1 
  ,  "archetypal_analysis" = 1        # Can choose if want to make scree plot
  ,    "cluster_stability" = 0
  ,         "arch_clinvar" = 1
  ,          "new_clinvar" = 1        # More clinical variable plots with new data (Oct 2020)
  ,  "enrichment_analysis" = 1        # Must choose which cluster to compare the rest with
  , "protein_summary_table"= 1
  ,   "manuscript_figures" = 1
)

# Options for normalisation
norm_opts <- c(
                  "AbsPQN" = 1
  ,                 "MA96" = 1
)

# Options for linear regression adjustment (lin_normalisation) of proteomics data
linnorm_opts <- c(
                     "age" = 1
  ,           "entry.date" = 1
  ,                  "bmi" = 1
  ,               "logreg" = 0    # Whether to do logistic regression for evaluation of case control effect, or not
  ,               "clogit" = 1    # Whether to use conditional logistic regression to evaluate case control effect, or not
)

# Adjust densities for BMI and age too?
                 adj_dens <- 1

# Archetypal analysis options
# For package, choose either:
# 1 for the "archetypes" package, or
# 2 for the "archetypal" package

# find_opt_k = 0 will just run the analysis with the number of 
# archetypes specified by n_arch
# find_opt_k = 1 will make a scree plot if archetypes package or 
# automatically find the optimal number of k for the archetypal package
arch_opts <- c(
                 "package" = 2
  ,           "find_opt_k" = 0
)

# If no scree plot is made, specify number of archetypes
                   n_arch <- 5

# For enrichment analysis, which cluster to compare to the rest
          enrichm_cluster <- 1

### ------ End of option choices ------ ###


# get names of chosen steps
step_names <- names(which(steps == 1))
step_list <- paste0(1:length(step_names), ": ", step_names)

# Tracker for counting steps
step_count <- 1

# Directory where output results folders are
out_dir <- "/home/.../results/reports/"
script_dir <- "/home/.../scripts/" 

# List to print the options in each output HTML
all_opts <- list("steps"=steps, 
                 "norm_opts"=norm_opts, 
                 "linnorm_opts"=linnorm_opts, 
                 "adj_dens"=adj_dens, 
                 "arch_opts"=arch_opts,
                 "n_arch"=n_arch,
                 "enrichm_cluster"=enrichm_cluster)

## Loading data 

if (steps["get_data"] == 1) {
  render(input=paste0(script_dir, "get_data.Rmd"),
       output_file=paste0(out_dir, Sys.Date(), "_", step_count,  "_get_data_report.html"),
       params=list(step_vec=step_list, 
                   step_nr=step_count,
                   start_time=starting_time))
  step_count <- step_count + 1
}


## Antibody quality control, criterion 1

if (steps["antibody_qc"] == 1) {
  source("/home/.../scripts/antibody_qc_crit1.R")
}


## Normalisation using MA and AbsPQN

if (steps["normalisation"] == 1) {
  render(input=paste0(script_dir, "normalisation.Rmd"),
         output_file=paste0(out_dir, Sys.Date(), "_", step_count, "_normalisation_report.html"),
         params=list(step_vec=step_list,
                     step_nr=step_count,
                     start_time=starting_time,
                     absPQN=norm_opts["AbsPQN"],
                     ma96=norm_opts["MA96"]))
  step_count <- step_count + 1
}


## Antibody quality control, criteria 2 and 3
if (steps["antibody_qc"] == 1) {
  render(input=paste0(script_dir, "antibody_qc.Rmd"),
       output_file=paste0(out_dir, Sys.Date(), "_", step_count, "_antibody_qc_report.html"),
       params=list(step_vec=step_list,
                   step_nr=step_count,
                   start_time=starting_time))
  step_count <- step_count + 1
}


# Linear adjustment of proteomics data
if (steps["lin_normalisation"] == 1) {
  render(input=paste0(script_dir, "lin_normalisation.Rmd"), 
         output_file=paste0(out_dir, Sys.Date(), "_", step_count, "_lin_normalisation_report.html"),
         params=list(step_vec=step_list,
                     step_nr=step_count,
                     start_time=starting_time,
                     norms=norm_opts))
  step_count <- step_count + 1
}


# Final preparation of the data for anlaysis, like removing antibodies that have the same targets
if (steps["analysis_prep"] == 1) {
  render(input=paste0(script_dir, "analysis_prep.Rmd"), 
         output_file=paste0(out_dir, Sys.Date(), "_", step_count, "_analysis_prep_report.html"),
         params=list(step_vec=step_list,
                     step_nr=step_count,
                     start_time=starting_time,
                     step_arch=arch_opts["make_scree"]))
  step_count <- step_count + 1
}


# Performing archetypal analysis
if (steps["archetypal_analysis"] == 1) {
  render(input=paste0(script_dir, "archetypal_analysis.Rmd"), 
         output_file=paste0(out_dir, Sys.Date(), "_", step_count, "_archetypal_analysis_report.html"),
         params=list(step_vec=step_list,
                     step_nr=step_count,
                     start_time=starting_time,
                     step_arch=arch_opts["make_scree"],
                     norms=norm_opts))
  step_count <- step_count + 1
}


# Performing cluster stability analysis
if (steps["cluster_stability"] == 1) {
  render(input=paste0(script_dir, "cluster_stability.Rmd"), 
         output_file=paste0(out_dir, Sys.Date(), "_", step_count, "_cluster_stability_report.html"))
  step_count <- step_count + 1
}


# Plot clinical variables of clusters
if (steps["arch_clinvar"] == 1) {
  render(input=paste0(script_dir, "arch_clinvar.Rmd"), 
         output_file=paste0(out_dir, Sys.Date(), "_", step_count, "_arch_clinvar_report.html"))
  step_count <- step_count + 1
}


# Plot new clinical variables
if (steps["new_clinvar"] == 1) {
  render(input=paste0(script_dir, "new_clinvar.Rmd"), 
         output_file=paste0(out_dir, Sys.Date(), "_", step_count, "_new_clinvar_report.html"))
  step_count <- step_count + 1
}

# Perform enrichment analysis on proteins between one cluster vs the rest
if (steps["enrichment_analysis"] == 1) {
  render(input=paste0(script_dir, "enrichment_analysis.Rmd"), 
         output_file=paste0(out_dir, Sys.Date(), "_", step_count, "_enrichment_analysis_report.html"))
  step_count <- step_count + 1
}


# Make a table summarising protein information
if (steps["protein_summary_table"] == 1) {
  render(input=paste0(script_dir, "prot_summary_table.Rmd"), 
         output_file=paste0(out_dir, Sys.Date(), "_", step_count, "_prot_summary.html"))
  step_count <- step_count + 1
}


# Make manuscript figures
if (steps["manuscript_figures"] == 1) {
	render(input=paste0(script_dir, "manuscript_figures.Rmd"),
				 output_file=paste0(out_dir, Sys.Date(), "_", step_count, "_manuscript_figures.html"))
	step_count <- step_count + 1
}

