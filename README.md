
### Karma case-control and clustering


### Scripts

`main_archetypal.R`  
- Pipeline script, contains commands to run scripts in correct order, as well as options for said scripts.

`get_data.Rmd`
- Loads data needed to get started.

`normalisation.Rmd`
- Normalises the data for plate effects (MA-normalisation) and sample-to-sample effects (AbsPQN).

`antibody_qc_crit1.R` and `antibody_qc.Rmd`
- Perform quality checks of antibodies and remove those that fail.

`lin_normalisation.Rmd`
- Adjusts data for confounding variables using linear regression.

`logreg_and_clogit.Rmd`
- Performs logistic regression and/or conditional logistic regression to investigate effect of case-control status.

`analysis_prep.Rmd`
- Prepares data for performing archetypal analysis.

`archetypal_analysis.Rmd`
- Performs archetypal analysis.

`cluster_stability.Rmd`
- Computes Jaccard index of clusters to evaluate stability.

`arch_clinvar.Rmd`
- Compares several clinical variables between clusters.

`new_clinvar.Rmd`
- Compares additional clinical variables between clusters, clinical data received later and where possible merged with old data.

`enrichment_analysis.Rmd`
- Performs enrichment analysis (ORA, GSEA) of a cluster of interest. Also differential expression analysis.

`protein_summary_table.Rmd`
- Makes a summary table for the antibodies and their protein targets.

`manuscript_figures.Rmd`
- Makes figures and tables for the manuscript.

`misc_functions.R`
- Contains a collection of functions used throughout the pipeline.

`plot_results.R`
- Contains a function for comparing distributions of clinical variables between clusters.

`norm_AbsPQN.R`
- Function for performing AbsPQN. Written by Mun-Gwan Hong.
