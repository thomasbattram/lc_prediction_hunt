## ---- load_data --------------------------------

pkgs <- c("tidyverse", "knitr", "pander")
lapply(pkgs, require, character.only = T)

dir <- "~/Desktop/projects/side_projects/lc_prediction_hunt/report/report_data/"

qc_sum <- read_tsv(paste0(dir, "qc_summary.txt"))

ewas_sum <- read_tsv(paste0(dir, "ewas_dat_summary.txt"))

# roc results
load(paste0(dir, "roc_dat.RData"))

# fix this in previous scripts!!
names(all_res)[2] <- "smoking_status_change"
ewas_sum[2,1] <- "smoking_status_change"

## ---- ewas_sum --------------------------------
pandoc.table(ewas_sum)

## ---- roc_setup --------------------------------
x=1

model_type <- c("separate_sites", "score")
roc_res <- lapply(1:length(model_type), function(x) {
	res <- map(all_res, model_type[x])
	out_res <- lapply(1:length(res), function(i) {
		res[[i]]$roc_dat
	})
	names(out_res) <- names(res)
	return(out_res)
})
names(roc_res) <- model_type

# generate AUC tables
auc_res <- map_dfr(1:length(all_res), function(x) {
	res_nam <- names(all_res)[x]
	out <- map_dfr(1:length(model_type), function(i) {
		auc_dat <- all_res[[res_nam]][[model_type[i]]]$auc %>%
			mutate(cpg_set = res_nam) %>%
			mutate(model = model_type[i]) %>%
			dplyr::select(cpg_set, model, estimate, lower, upper) %>%
			rename(auc = estimate, ci_lower = lower, ci_upper = upper)
	})
	return(out)
})

p <- lapply(1:length(model_type), function(x) {
	model <- model_type[x]
	out <- pROC::ggroc(roc_res[[model]]) +
		geom_abline(intercept = 1, slope = 1, colour = "black", alpha = 0.6) +
		theme_bw() +
		labs(colour = NULL)
})
names(p) <- model_type

## ---- roc_plot_separate_sites --------------------------------
print(p$separate_sites)

## ---- roc_plot_score --------------------------------
print(p$score)

## ---- auc_tab --------------------------------
pandoc.table(auc_res)