# --------------------------------------------
# Predicting lung cancer in NSHDS with DNA methylation data
# --------------------------------------------

# 76 CpGs were identified in an EWAS of lung cancer in HUNT and the authors
# are looking to test how well the predict lung cancer in an independent
# dataset --> NSHDS.

pkgs <- c("tidyverse", "survival", "pROC")
lapply(pkgs, require, character.only = TRUE)

devtools::load_all("~/repos/usefunc/")

# read in NSHDS data
data_file <- "data/cleaned_pheno_and_meth_data.RData"
load(data_file)

# read in 76 CpGs
int_cpgs <- read_csv("data/76_CpGs_HUNT_YQS.csv", skip=1)

meth_data_file <- "~/tb_MB009_LungCancer/RDATA/FunctionalNormalised_Betas_PvaluePassed_Samples.RData"
load(meth_data_file) 


# ----------------------------------- 
# extract variables needed from nshds data
# -----------------------------------

meth_res <- out$meth %>%
	dplyr::filter(cpg %in% int_cpgs[["CPG Labels"]])
# only 41 left...
which(!int_cpgs[["CPG Labels"]] %in% out$meth$cpg)
# aries check
# load("/panfs/panasas01/sscm/ms13525/aries-release-v4/data/betas/data.Robj")
# sum(rownames(beta) %in% int_cpgs[["CPG Labels"]]) # 41
#

# only 41 because they used EPIC array data!!!!

# derive DNAm score for each individual
rownames(meth_res) <- meth_res$cpg
meth_res2 <- t(meth_res)[-1, ]
meth_res2 <- meth_res2 %>%
	as.data.frame %>%
	rownames_to_column(var = "sample")

fin_res <- map_dfc(1:ncol(meth_res2), function(x) {
	col_nam <- colnames(meth_res2)[x]
	if (col_nam == "sample") return(meth_res2[[x]])
	
	vals <- as.numeric(as.character(meth_res2[[x]]))
	coef <- int_cpgs %>%
		dplyr::filter(`CPG Labels` == col_nam) %>%
		pull(Coefficient)
	out <- coef * vals
	return(out)
})
colnames(fin_res) <- colnames(meth_res2)


# checking smoking variables to find pack years!
grep("SMOK", colnames(out$pheno), value = TRUE)
# looks like SMOKE_INT_CIGARETTE may be a variable that has number 
# smoked per day... but need to go back to people who sent over dataset
# to be sure

phen_res <- out$pheno %>%
	dplyr::select(sentrix, CASESET, LUNG_CANCER_CASE, DOB) %>%
	left_join(fin_res, by = c("sentrix" = "sample"))

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

cpg_sites <- grep("cg", colnames(phen_res), value = TRUE)
vars1 <- paste(paste(cpg_sites, collapse = " + "), "strata(CASESET)", sep = " + ")
form1 <- as.formula(paste0("LUNG_CANCER_CASE ~ ", vars1))
fit1 <- clogit(form1, data = phen_res)

ahrr_cpg <- "cg05575921"
vars2 <- paste(paste(ahrr_cpg, collapse = " + "), "strata(CASESET)", sep = " + ")
form2 <- as.formula(paste0("LUNG_CANCER_CASE ~ ", vars2))
fit2 <- clogit(form2, data = phen_res)

roc_all <- roc(phen_res$LUNG_CANCER_CASE, fit1$linear.predictors, ci = T)
roc_ahrr <-roc(phen_res$LUNG_CANCER_CASE, fit2$linear.predictors, ci = T)
res <- roc.test(roc_all, roc_ahrr)

roc_list <- list(all = roc_all, ahrr = roc_ahrr)
i=1
roc_res <- lapply(1:2, function(i) {
	roc_dat <- roc_list[[i]]
	auc_dat <- t(as.data.frame(ci.auc(roc_dat))) %>%
		as.data.frame
	rownames(auc_dat) <- NULL
	colnames(auc_dat) <- c("lower", "estimate", "upper")
	plot_text <- paste0("AUC = ", comma(auc_dat$estimate), " (95% CI: ",
					comma(auc_dat$lower), " - ", comma(auc_dat$upper), ")")
	return(list(roc_dat = roc_dat, auc = auc_dat, plot_text = plot_text))
})
names(roc_res) <- names(roc_list)

save(roc_res, file = "report/report_data/roc_dat.RData")

p <- pROC::ggroc(list(roc_res$ahrr$roc_dat)) +
	geom_abline(intercept = 1, slope = 1, colour = "black", alpha = 0.6) +
	annotate("text", x = 0.7, y = 0.9, label = plot_text) +
	theme_bw() +
	theme(legend.position = "none")
ggsave("results/roc_plot.pdf", plot = p)



# models
GLM <- glm(Y ~. , data = X.base, family = binomial(logit))
# Y is LC variable (0, 1); X.base is data.frame including variables sex, age and smoking (7-levels score)
ROC.base <- roc(GLM$y, GLM$fitted.values, ci = T)
GLM <- glm(Y ~. , data = X.plus.76, family = binomial(logit))
# Y is LC variable (0, 1); X.plus.76 is data.frame including variables sex, age, smoking (7-levels score) and 76 CpGs
ROC.plus.76 <-roc(GLM$y, GLM$fitted.values, ci = T)
# Compare 2 ROCs:
roc.test(ROC.base, ROC.plus.76)

