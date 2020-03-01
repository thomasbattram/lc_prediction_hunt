## ---- load_data --------------------------------

pkgs <- c("tidyverse", "knitr")
lapply(pkgs, require, character.only = T)

qc_sum <- read_tsv("~/Desktop/projects/side_projects/lc_prediction_hunt/report/report_data/qc_summary.txt")
