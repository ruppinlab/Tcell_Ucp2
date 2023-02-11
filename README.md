# Title

description

## 1. Prerequisites

### Softwares

* R 3.6
* R packages
  - ImNotaGit/my.utils
  - ruppinlab/Rcplex2: required to run genome-scale metabolic modeling (GEM); otherwise may be omitted
  - ruppinlab/gembox
  - other needed packages can be obtained from CRAN or Bioconductor, please see the R scripts
* IBM ILOG CPLEX Optimization Studio 12: required to run GEM; otherwise may be omitted

### Data

Download the data files from [here](url) then decompress them into the **data** folder; or it may be possible to use the dnload.sh script in the **data** folder.

## 2. Description of folders and files

### Fraietta

Prediction of metabolic pathways important for T cell function in the context of anti-CD19 CAR-T therapy, and metabolic flux analysis using the data from Fraietta et al. 2018.

* prepare.data.R: prepare data for GEM
* run.mta.R: run the MTA algorithm to predict metabolic reactions whose knockout can result in non-responsiveness in anti-CD19 CAR-T therapy
* run.flux.analysis.R: run metabolic flux analysis comparing the responders vs non-responders of anti-CD19 CAR-T therapy

### Lu

Metabolic flux analysis of the persistent and non-persistent T cell clones in adoptive cell transfer therapy using the data from Lu et al. 2019.

* prepare.data.R: prepare data for GEM
* run.flux.analysis.R: run metabolic flux analysis comparing the persistent vs non-persistent T cell clones

### TCGA

Analysis of UCP2 gene expression in different cancer types using the TCGA dataset.

* check.ucp2.association.R: analyze the association between UCP2 expression and T cell memory/stemness genes, and patient survival

### Figures

R notebooks for generating the some of the figures in the paper.

* figure1.Rmd
* figureS1.Rmd
