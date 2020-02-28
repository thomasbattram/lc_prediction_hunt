# --------------------------------------------
# Predicting lung cancer in NSHDS with DNA methylation data
# --------------------------------------------

# 76 CpGs were identified in an EWAS of lung cancer in HUNT and the authors
# are looking to test how well the predict lung cancer in an independent
# dataset --> NSHDS.

pkgs <- c("tidyverse")
lapply(pkgs, require, character.only = TRUE)



# read in NSHDS data
data_file<-"~/tb_MB009_LungCancer/RDATA/FunctionalNormalised_Betas_PvaluePassed_Samples.RData"
