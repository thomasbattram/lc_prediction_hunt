# --------------------------------------------
# Cleaning lung cancer and DNA methylation data in NSHDS
# --------------------------------------------

pkgs <- c("tidyverse", "minfi", "survival", "CpGassoc")
lapply(pkgs, require, character.only = TRUE)

devtools::load_all("~/repos/usefunc")

#batch variables
batch_vars<-c("MSA4Plate_id","chip.id","BCD_id")

data_dir <- "~/tb_MB009_LungCancer/"

#load the methylation data
meth_data_file <- "~/tb_MB009_LungCancer/RDATA/FunctionalNormalised_Betas_PvaluePassed_Samples.RData"
load(meth_data_file) 
meth_dat <- betaFnNorm
# detection p-value data
det_p_value_data <- "~/tb_MB009_LungCancer/RDATA/pvalues.RData"
load(det_p_value_data)
# exposure data and covariates
phenotype_data <- read_csv("~/tb_MB009_LungCancer/Other/Clinical_Information_Sheet.csv", guess_max = 1e6)
batch_data <- read.csv("~/tb_MB009_LungCancer/Other/Final_Samples_Info_File.csv", header = T, stringsAsFactors = F)
#load the cell counts
load("~/tb_MB009_LungCancer/RDATA/HousemanCellCount_AllSamples.RData")

# --------------------------------------------
# remove samples and probes based on detection p values
# --------------------------------------------

ori_dim_p_values <- dim(p.values) # 485512 probes and  485 samples
# exclude probes not in the beta file due to previous QC
p_values <- p.values[rownames(p.values) %in% rownames(meth_dat), ] 
nrow(p_values) # 485512 probes left, 0 removed
# exclude samples not in the beta file due to previous QC
p_values <- p_values[, colnames(p_values) %in% colnames(meth_dat)]
ncol(p_values) # 478 samples left, 7 removed

# exclude probes that have detection p values of >0.01 on 5% or more samples
# and exclude samples that have detection p values of >0.01 for 5% or more probes
no_failed_probes <- p_values %>%
	mutate(failed_probes = rowSums(. > 0.01)) %>%
	mutate(cpg = rownames(p_values)) %>%
	dplyr::filter(failed_probes < (0.05 * ncol(.)))

nrow(no_failed_probes) # 474237 probes left, 11275 removed

failed_samples <- colSums(p_values > 0.01)

to_rm <- names(which(failed_samples > (0.05 * nrow(p_values))))

no_failed_dat <- no_failed_probes %>% 
	dplyr::select(-to_rm, -failed_probes)
ncol(no_failed_dat) # 478 (477 samples and 1 CpG column), 1 sample removed

# exclude poor probes and samples from beta file
meth_dat <- meth_dat %>%
	rownames_to_column(var = "cpg") %>%
	dplyr::filter(cpg %in% no_failed_dat$cpg) %>%
	dplyr::select(one_of(colnames(no_failed_dat)))

# make sure p_values and beta file are in same order, rowwise and columnwise
if (!all(no_failed_dat$cpg == meth_dat$cpg)) {
	index <- match(meth_dat$cpg, no_failed_dat$cpg)
	no_failed_dat <- no_failed_dat[index, ]
}

if (!all(colnames(no_failed_dat) == colnames(meth_dat))) {
	index <- match(colnames(meth_dat), colnames(no_failed_dat))
	no_failed_dat <- no_failed_dat[, index]
}

# set failed probes to missing 
i=1
meth_res <- map_dfc(1:ncol(meth_dat), function(i) {
	samp_name <- colnames(meth_dat)[i]
	if (samp_name == "cpg") return(meth_dat$cpg)
	to_na <- which(no_failed_dat[[samp_name]] > 0.01)
	res <- meth_dat[, i]
	res[to_na] <- NA
	return(res)
})
colnames(meth_res) <- colnames(no_failed_dat)

# QC summary
qc_sum <- data.frame(qc_stage = c("original_data", "det_p_removed"), 
					 n_cpg = c(nrow(betaFnNorm), nrow(meth_res) - 1), 
					 n_sample = c(ncol(betaFnNorm), ncol(meth_res) - 1))

# ------------------------------------------------
# sort the phenotype data
# ------------------------------------------------

# sort the cell counts data out
count_dat <- data.frame(cellCount$counts) %>% 
	rownames_to_column(var = "sentrix")

# join the batch, phenotype and cell counts datasets together
phen_dat <- batch_data %>%
	mutate(LC3_ID = gsub("_[A-Z][0-9]*$", "", ID)) %>%
	left_join(phenotype_data) %>%
	dplyr::filter(!duplicated(LC3_ID)) %>%
	left_join(count_dat) %>%
	dplyr::filter(sentrix %in% colnames(meth_res))


# exlude incomplete case control pairs
pairs <- phen_dat$CASESET[duplicated(phen_dat$CASESET)]
phen_dat <- phen_dat %>%
	dplyr::filter(CASESET %in% pairs)
meth_res <- meth_res %>%
	dplyr::select(cpg, one_of(phen_dat$sentrix))

out <- list(pheno = phen_dat, meth = meth_res)
save(out, file = "~/lc_prediction_hunt/data/cleaned_pheno_and_meth_data.RData")

# write out qc summary
write.table(qc_sum, file = "data/qc_summary.txt", 
			quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(qc_sum, file = "report/report_data/qc_summary.txt", 
			quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

#create chip ID
chip.id<-unlist(strsplit(phen.all$sentrix,split="_"))
length(phen.all$sentrix)
phen.all$chip.id<-chip.id[seq(1,length(chip.id),2)]

#generate BMI and normalise 
phen.all$bmi<-phen.all$WEIGHT/((phen.all$HEIGHT/100)^2)
phen.all$lnbmi<-log(phen.all$bmi)

#create binary alcohol variable
phen.all$alc_stat<-phen.all$ALC_ETOH
phen.all$alc_stat[which(phen.all$alc_stat>0)]<-1

#table(phen.all$SMOKE_INT_CIGARETTE,phen.all$SMOKE_STATUS) #102 NAs for ex smokers
#table(phen.all$SMOKE_DUR,phen.all$SMOKE_STATUS)
#summary(phen.all$SMOKE_DUR[phen.all$SMOKE_STATUS==4])
phen.all$SMOKE_DUR[phen.all$SMOKE_STATUS==4]<-0 #set duration to 0 for never smokers
#summary(phen.all$SMOKE_DUR[phen.all$SMOKE_STATUS==5]) #4 NAs for ex smokers
#summary(phen.all$SMOKE_DUR[phen.all$SMOKE_STATUS==6]) #15 NAs for current smokers
#summary(phen.all$SMOKE_AGE_QUIT[phen.all$SMOKE_STATUS==5]) #no NAs for ex smokers
#summary(phen.all$SMOKE_AGE_QUIT[phen.all$SMOKE_STATUS==6]) #267 NAs for current smokers #4 = never smoker, 5 = ex smoker, 6 current smoker (but 3 current smokers with age quit)
#table(phen.all$SMOKE_AGE_QUIT,phen.all$SMOKE_STATUS)

#sort(names(phen.all))
phen.all$yrs.since.quit<-phen.all$DRAWAGE-phen.all$SMOKE_AGE_QUIT
#######################
#define the covariates#
#######################
#recode/factorise the covariates from analysis plan: age + sex + education + BMI + smoking_status 
#excluding subjects with missing outcomes#
phen.all$SMOKE_INT_CIGARETTE<-as.factor(phen.all$SMOKE_INT_CIGARETTE)
phen.all$EDUCATION<-as.factor(phen.all$EDUCATION)
phen.all$GENDER[phen.all$GENDER==1]<-0
phen.all$GENDER[phen.all$GENDER==2]<-1
phen.all$BCD_id<-as.factor(phen.all$BCD_id) #bisulphite modification batch
phen.all$MSA4Plate_id<-as.factor(phen.all$MSA4Plate_id) #illumina processing batch
phen.all$chip.id<-as.factor(phen.all$chip.id)  #illumina chip ID 
phen.all$SMOKE_STATUS<-as.factor(phen.all$SMOKE_STATUS)


#load the map information for the manhattan plot function
load("Data/fdata_new.RData")

table.450<-read.table("Data/humanmethylation450_15017482_v1-2.csv",sep=",",head=T,stringsAsFactors=F,skip=7,quote="",fill=T)
col.names<-c("IlmnID", "Name" , "Genome_Build", "CHR","MAPINFO","Chromosome_36","Coordinate_36", "Random_Loci","UCSC_RefGene_Name",
"UCSC_RefGene_Accession" ,  "UCSC_RefGene_Group","UCSC_CpG_Islands_Name","Relation_to_UCSC_CpG_Island","Phantom" ,"DMR","Enhancer" ,"HMM_Island","Regulatory_Feature_Name" ,
"Regulatory_Feature_Group","DHS")

#make sure number and order of CpG is same across fdata.new and beta datasets
#fdata.new is the map file for the CpG sites, containing the chromosome & bp positions
fdata.new<-fdata.new[rownames(fdata.new) %in% rownames(beta.type),]
fdata.new<-fdata.new[order(match(rownames(fdata.new), rownames(beta.type))),]
names(fdata.new)[which(names(fdata.new)=="CHR37")]<-"chromosome"
names(fdata.new)[which(names(fdata.new)=="COORDINATE_37")]<-"position"
features<-fdata.new

#logit transformation, M values
N999<-NULL
N001<-NULL
mbetas<-beta.type
for(i in 1:ncol(mbetas)){
#	print(i)
	N999[[i]]<-unlist(length(which(mbetas[,i]>0.999)))
	N001[[i]]<-unlist(length(which(mbetas[,i]<0.001)))
	mbetas[which(mbetas[,i]>0.999),i]<-0.999
	mbetas[which(mbetas[,i]<0.001),i]<-0.001
}
mbetas<-log2(mbetas/(1-mbetas))

##betas<0.001 or >0.999
#png('Results/hist_extreme_large_betas.png')
#hist(N999,main="",xlab="",freq=T)
#title("Number of probes per sample with a beta >0.999")
#legend(-0.4,200,paste("median=",median(N999),
#	" IQR=",quantile(N999,0.25),"-",quantile(N999,0.75),sep=""))
#dev.off()

#png('Results/hist_extreme_small_betas.png')
#hist(N001,main="",xlab="",freq=T)
#title("Number of probes per sample with a beta <0.001")
#legend(1,100,paste("median=",median(N001),
#	" IQR=",quantile(N001,0.25),"-",quantile(N001,0.75),sep=""))
#legend(1,90,paste("mean=",round(mean(N001),2),
#	" SD=",round(sd(N001),2),sep=""))
#dev.off()

Beta<-"beta.type"
if(logit.transform){
	Beta<-"mbetas"
}

