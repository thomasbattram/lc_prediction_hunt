# --------------------------------------------
# Predicting lung cancer in NSHDS with DNA methylation data
# --------------------------------------------

# CpGs were identified in an EWAS of smoking, smoking status change
# and lung cancer in HUNT and the authors are looking to test how well
# the CpGs predict lung cancer in an independent dataset --> NSHDS. 

pkgs <- c("tidyverse", "survival", "pROC", "readxl")
lapply(pkgs, require, character.only = TRUE)

devtools::load_all("~/repos/usefunc/")

# read in NSHDS data
data_file <- "data/cleaned_pheno_and_meth_data.RData"
load(data_file)

# read in CpGs of interest
smok_ewas <- read_xlsx("data/hunt_ewas_data.xlsx", sheet = 1)
smok_stat_ewas <- read_xlsx("data/hunt_ewas_data.xlsx", sheet = 2)
lc_ewas <- read_xlsx("data/hunt_ewas_data.xlsx", sheet = 3)

# meth_data_file <- "~/tb_MB009_LungCancer/RDATA/FunctionalNormalised_Betas_PvaluePassed_Samples.RData"
# load(meth_data_file) 


# ----------------------------------- 
# extract variables needed from nshds data
# -----------------------------------

ewas <- list(smoking = smok_ewas, 
			 smoking_status = smok_stat_ewas, 
			 lung_cancer = lc_ewas)

cpgs <- list(smoking = smok_ewas$cpg, 
			 smoking_status = smok_stat_ewas$cpg, 
			 lung_cancer = lc_ewas$cpg)

unique_cpgs <- unique(unlist(cpgs))

meth_res <- lapply(1:3, function(x) {
	cpg_list <- cpgs[[x]]
	dat <- out$meth %>%
		dplyr::filter(cpg %in% cpg_list)
	return(dat)
})
names(meth_res) <- names(cpgs)

# only 51 left...
which(!unique_cpgs %in% out$meth$cpg)

# how many are there in each cpg dataset

# ewas n_cpgs n_cpgs_in450k unique_cpgs_in450k
#
#
x=1
cpg_dat_sum <- map_dfr(1:3, function(x) {
	cpg_list <- cpgs[[x]]
	other_cpg_lists <- unique(unlist(cpgs[-x]))
	ewas_nam <- names(cpgs[x])
	ncpgs <- length(cpg_list)
	list_450k <- cpg_list[cpg_list %in% out$meth$cpg]
	list_unique <- list_450k[!list_450k %in% other_cpg_lists]
	out <- data.frame(
		ewas = ewas_nam, 
		n_cpgs = ncpgs, 
		n_cpgs_in_450k = length(list_450k), 
		n_unique_cpgs_in_450k = length(list_unique)
		)
})

overall <- data.frame(
			ewas = "combined", 
			n_cpgs = length(unlist(cpgs)), 
			n_cpgs_in_450k = sum(unique_cpgs %in% out$meth$cpg), 
			n_unique_cpgs_in_450k = sum(unique_cpgs %in% out$meth$cpg)
			)
cpg_dat_sum <- bind_rows(list(cpg_dat_sum, overall))

write.table(cpg_dat_sum, "report/report_data/ewas_dat_summary.txt", 
			row.names = F, col.names = T, quote = F, sep = "\t")

# aries check
# load("/panfs/panasas01/sscm/ms13525/aries-release-v4/data/betas/data.Robj")
# sum(rownames(beta) %in% int_cpgs[["CPG Labels"]]) # 41
#

# only 51 because they used EPIC array data!!!!

# derive DNAm score for each individual
lc_ewas <- lc_ewas %>%
	mutate(Coefficient = log(ORb))

ewas[["lung_cancer"]] <- lc_ewas

meth_res2 <- lapply(1:3, function(x) {
	print(x)
	# orient the data
	dat <- meth_res[[x]]
	rownames(dat) <- dat$cpg
	dat2 <- t(dat)[-1, ]
	dat2 <- dat2 %>%
		as.data.frame %>%
		rownames_to_column(var = "sample")

	ewas_dat <- ewas[[x]]
	fin_res <- map_dfc(1:ncol(dat2), function(x) {
		col_nam <- colnames(dat2)[x]
		if (col_nam == "sample") return(dat2[[x]])
		
		vals <- as.numeric(as.character(dat2[[x]]))
		coef <- ewas_dat %>%
			dplyr::filter(cpg == col_nam) %>%
			pull(Coefficient)
		out <- coef * vals
		return(out)
	})
	colnames(fin_res) <- colnames(dat2)
	return(fin_res)
})
names(meth_res2) <- names(cpgs)


# checking smoking variables to find pack years!
grep("SMOK", colnames(out$pheno), value = TRUE)
# looks like SMOKE_INT_CIGARETTE may be a variable that has number 
# smoked per day... but need to go back to people who sent over dataset
# to be sure

phen_res <- out$pheno %>%
	dplyr::select(sentrix, CASESET, LUNG_CANCER_CASE, DOB) # %>%
	# left_join(fin_res, by = c("sentrix" = "sample"))

# -----------------------------------
# run the analysis
# -----------------------------------

### tests
cpg_site <- grep("cg", colnames(phen_res), value = T)[1]
vars1 <- paste("strata(CASESET)", sep = " + ")
form1 <- as.formula(paste0("LUNG_CANCER_CASE ~ ", vars1))
fit1 <- clogit(form1, data = phen_res)

vars2 <- paste(cpg_site, "strata(CASESET)", sep = " + ")
form2 <- as.formula(paste0("LUNG_CANCER_CASE ~ ", vars2))
fit2 <- clogit(form2, data = phen_res)

vars2_nodob <- paste(cpg_site, "strata(CASESET)", sep = " + ")
form2_nodob <- as.formula(paste0("LUNG_CANCER_CASE ~ ", vars2))
fit2_nodob <- clogit(form2, data = phen_res)

ROC.base <- roc(phen_res$LUNG_CANCER_CASE, fit1$linear.predictors, ci = T)
ROC.plus.cpg <-roc(phen_res$LUNG_CANCER_CASE, fit2$linear.predictors, ci = T)
ROC.plus.cpg.nodob <- roc(phen_res$LUNG_CANCER_CASE, fit2_nodob$linear.predictors, ci = T)
# Compare 2 ROCs:
roc.test(ROC.base, ROC.plus.cpg)
roc.test(ROC.base, ROC.plus.cpg.nodob)
###### ------------------------------------------------------------

#  generated taking into account the smoking status and pack-years
# (pyrs): [0: never smokers; 1: former ≤10.0 pyrs; 2: former 10.1-20.0 pyrs
# ; 3: former ≥20.1 pyrs; 4: current ≤10.0 pyrs; 5: current 10.1-20.0 pyrs;
# and 6: current ≥20.1 pyrs].

cpgs <- c(cpgs, list(ahrr = "cg05575921"))
model_type <- c("separate_sites", "score")
x=1
all_res <- lapply(1:length(cpgs), function(x) {
	print(x)
	ewas_nam <- names(cpgs)[x]
	if (x == 4) {
		meth_dat <- meth_res2[[1]] %>%
			dplyr::select(sample, cg05575921)
	} else {
		meth_dat <- meth_res2[[ewas_nam]]
	}
	phen_dat <- phen_res %>%
		left_join(meth_dat, by = c("sentrix" = "sample"))
	cpg_sites <- colnames(phen_dat)[colnames(phen_dat) %in% unique_cpgs]
	phen_dat$cpg_score <- rowSums(phen_dat[, cpg_sites, drop = FALSE])
	out_res <- lapply(1:2, function(i) {
		model <- model_type[i]
		if (model == "separate_sites") {
			vars <- paste(paste(cpg_sites, collapse = " + "), "strata(CASESET)", sep = " + ")
		} else if (model == "score") {
			vars <- paste("cpg_score", "strata(CASESET)", sep = " + ")
		}
		form <- as.formula(paste0("LUNG_CANCER_CASE ~ ", vars))
		fit <- clogit(form, data = phen_dat)
		# roc curve
		roc_res <- roc(phen_dat$LUNG_CANCER_CASE, fit$linear.predictors, ci = T)
		# auc
		auc_dat <- t(as.data.frame(ci.auc(roc_res))) %>%
			as.data.frame
		rownames(auc_dat) <- NULL
		colnames(auc_dat) <- c("lower", "estimate", "upper")
		plot_text <- paste0(ewas_nam, ": ", comma(auc_dat$estimate), " (95% CI: ",
						comma(auc_dat$lower), " - ", comma(auc_dat$upper), ")")
		return(list(roc_dat = roc_res, auc = auc_dat, plot_text = plot_text))
	})
	names(out_res) <- model_type
	return(out_res)
})
names(all_res) <- names(cpgs)

save(all_res, file = "report/report_data/roc_dat.RData")

# p <- pROC::ggroc(list(roc_res$ahrr$roc_dat)) +
# 	geom_abline(intercept = 1, slope = 1, colour = "black", alpha = 0.6) +
# 	annotate("text", x = 0.7, y = 0.9, label = plot_text) +
# 	theme_bw() +
# 	theme(legend.position = "none")
# ggsave("results/roc_plot.pdf", plot = p)
