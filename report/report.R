## ---- load_data --------------------------------

pkgs <- c("tidyverse", "knitr")
lapply(pkgs, require, character.only = T)

dir <- "~/Desktop/projects/side_projects/lc_prediction_hunt/report/report_data/"

qc_sum <- read_tsv(paste0(dir, "qc_summary.txt"))

# roc results
load(paste0(dir, "/roc_dat.RData"))

## ---- roc_setup --------------------------------

p <- pROC::ggroc(list(ROC.plus.cpg)) +
	geom_abline(intercept = 1, slope = 1, colour = "black", alpha = 0.6) +
	annotate("text", x = 0.7, y = 0.9, label = plot_text) +
	theme_bw() +
	theme(legend.position = "none")

## ---- roc_plot --------------------------------
print(p)