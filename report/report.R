## ---- load_data --------------------------------

pkgs <- c("tidyverse", "knitr", "pander")
lapply(pkgs, require, character.only = T)

dir <- "~/Desktop/projects/side_projects/lc_prediction_hunt/report/report_data/"

qc_sum <- read_tsv(paste0(dir, "qc_summary.txt"))

ewas_sum <- read_tsv(paste0(dir, "ewas_dat_summary.txt"))

# roc results
load(paste0(dir, "/roc_dat.RData"))

# fix this in previous scripts!!
names(all_res)[2] <- "smoking_status_change"
ewas_sum[2,1] <- "smoking_status_change"

## ---- ewas_sum --------------------------------
pandoc.table(ewas_sum)

## ---- roc_setup --------------------------------
x=1
roc_res <- lapply(1:length(all_res), function(x) {
	res_nam <- names(all_res)[x]
	all_res[[x]]$roc_dat
})
names(roc_res) <- names(all_res)

# generate AUC tables
auc_res <- map_dfr(1:length(all_res), function(x) {
	res_nam <- names(all_res)[x]
	auc_dat <- all_res[[res_nam]]$auc %>%
		mutate(cpg_set = res_nam) %>%
		dplyr::select(cpg_set, estimate, lower, upper) %>%
		rename(auc = estimate, ci_lower = lower, ci_upper = upper)
	return(auc_dat)
})

## why doesn't this work!!!
p <- pROC::ggroc(roc_res) +
	geom_abline(intercept = 1, slope = 1, colour = "black", alpha = 0.6) +
	theme_bw() +
	labs(colour = NULL)
	# theme(legend.position = "none")

## ---- roc_plot --------------------------------
print(p)

## ---- auc_tab --------------------------------
pandoc.table(auc_res)