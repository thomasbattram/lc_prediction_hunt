## ---- load_data --------------------------------

pkgs <- c("tidyverse", "knitr")
lapply(pkgs, require, character.only = T)

dir <- "~/Desktop/projects/side_projects/lc_prediction_hunt/report/report_data/"

qc_sum <- read_tsv(paste0(dir, "qc_summary.txt"))

# roc results
load(paste0(dir, "/roc_dat.RData"))

## ---- roc_setup --------------------------------

res <- list(roc_res$all$roc_dat, roc_res$ahrr$roc_dat)
names(res) <- c("all_cpgs", "cg05575921 (AHRR)")
## why doesn't this work!!!
p_both <- pROC::ggroc(res) +
	geom_abline(intercept = 1, slope = 1, colour = "black", alpha = 0.6) +
	annotate("text", x = 0.7, y = 0.95, label = roc_res$all$plot_text) +
	annotate("text", x = 0.7, y = 0.9, label = roc_res$ahrr$plot_text) +
	theme_bw() +
	labs(colour = NULL)
	# theme(legend.position = "none")

## ---- roc_plot --------------------------------
print(p_both)