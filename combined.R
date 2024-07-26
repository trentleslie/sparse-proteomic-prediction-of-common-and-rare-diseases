###################################################
#### Script Name: generate_input_file
#### Author: Julia Carrasco Zanini Sanchez
#### Description: generates inpute file to run ML framework in a job array 
#### Date : 20/07/2022

rm(list = ls())

library(tidyr)
library(dplyr)

## set working directory 
setwd("/home/ukbb_inc_dz_prediction/Explore_1536_Expansion_analyses/01_inc_dz_phenos/bin/") 

## read in disease list 
dz.list <- read.delim("../data_input/atlas308_subc_inc6months_inc_dz_counts_3_10yrs.txt", sep = "\t")
head(dz.list)

dz.input <- dz.list %>% 
  select(inc_dz_name, dz_date) %>% 
  mutate(inc_yrs=rep(10,nrow(dz.list)))

# dz.input <- rbind(dz.input, dz.list %>% 
#   select(inc_dz_name, dz_date) %>% 
#   mutate(inc_yrs=rep(5,nrow(dz.list)))) %>% 
#   filter(!inc_dz_name%in%c("atlas_inc_ms","atlas_inc_sle"))

write.table(dz.input, file = "../../02_inc_dz_prediction/data_input/atlas308_file_input_params.txt", sep="\t", col.names = F, row.names = F, quote = F)

## diseases with fewer cases include consortium selected cases
dz.list.cc <- read.delim("../data_input/atlas308_consort_inc6months_inc_dz_counts_3_10yrs.txt", sep = "\t")

dz.input.cc <- dz.list.cc %>% 
  select(inc_dz_name, dz_date) %>% 
  mutate(inc_yrs=rep(10,nrow(dz.list.cc)))

dz.input.cc <- rbind(dz.input.cc, dz.list.cc %>% 
        select(inc_dz_name, dz_date) %>% 
        mutate(inc_yrs=rep(5,nrow(dz.list.cc))))
write.table(dz.input.cc, file = "../../02_inc_dz_prediction/data_input/atlas308_file_input_params_consort.txt", sep="\t", col.names = F, row.names = F, quote = F)

# list for PRS prediction models
prs.input <- dz.input %>%
  filter(inc_dz_name%in%c("atlas_inc_asthma","atlas_inc_af",
                          "atlas_inc_pri_bowel","atlas_inc_pri_breast",
                          "atlas_inc_chd_nos","atlas_inc_crohns",
                          "atlas_inc_pri_ovarian", "atlas_inc_fracture_hip",
                          "atlas_inc_fracture_wrist","atlas_inc_isch_stroke",
                          "atlas_inc_pri_skin",
                          "atlas_inc_oa","atlas_inc_parkinsons",
                          "atlas_inc_glaucoma","atlas_inc_pri_prost",
                          "atlas_inc_psoriasis","atlas_inc_rha",
                          "atlas_inc_diabetes_t2","atlas_inc_ulc_colitis",
                          "atlas_inc_vte_ex_pe"))
write.table(prs.input, file = "../../02_inc_dz_prediction/data_input/prs_input_params.txt", sep="\t", col.names = F, row.names = F, quote = F)
###################################################
#### Script Name: ML_pipeline_6m_protein_pxinfo_prediction_pred_consortcases
#### Author: Julia Carrasco Zanini Sanchez
#### Description: script to ML framework for incident diseases in atlas 308 
#### Date : 23/03/2023

rm(list = ls())

## load packages 
library(dplyr)
library(tidyr)
library(RNOmni)
# library(glue)
library(DBI)
library(dbplyr)
library(magrittr)
# library(broom)
library(survival)
library(data.table)
library(glue)
library(parallel)
library(caret)
library(glmnet)
library(ROSE)

## define disease and restrict to incident cases within X years
args <- commandArgs(trailingOnly = T)

dz <- as.character(args[1])
dz.date <- as.character(args[2])
inc.yrs <- as.numeric(args[3])

paste0("Prediction:", dz, "; ", dz.date, " at ", inc.yrs, "years")

## set up working directory 
setwd("/home/ukbb_inc_dz_prediction/Explore_1536_Expansion_analyses/02_inc_dz_prediction/bin/") 

## read in data 
atlas308 <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_atlas308_annotations.txt", data.table = F)
cov <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_covariate_annotations.txt", data.table = F)

####----  Data Cleaning and wrangling   ----####
atlas308 <- as.data.frame(atlas308 %>% 
  filter(eid_20361%in%cov$eid_20361) %>% 
  pivot_wider(names_from = variable, values_from =value))

## phenotypes missing for patient-derived info model: bmi, smoking, alcohol, physical activity, family history 
pheno <- merge(cov[,-c(11:14,17)], atlas308, by="eid_20361")

## read in px info phenotypes
px.info <- read.delim("../../01_inc_dz_phenos/data_input/Patient_model_variables.txt",sep = "\t")
pheno <- merge(pheno, px.info, by="eid_20361")

## rename to olds names to match the rest fo the script
colnames(pheno)[which(colnames(pheno)=="assess_date")] <- c("date_baseline_assessment")

## process ethnicity
pheno$ethnicity <- ifelse(is.na(pheno$ge_eur),2,pheno$ge_eur)

## get protein data and metadata
ol.npx <- fread("../../00_protein_imputation/data_input/Imputed_1536_Expansion_NPX_proteins.txt", data.table = F)
ol.npx <- filter(ol.npx, exclusion.50==0 )

## exclude participants with missing age, sex
ol.npx <- ol.npx[complete.cases(ol.npx[,c("ge_sex", "ge_age")]),]

protein_list <- read.delim("../../00_protein_imputation/data_input/protein_labels_by_panel.txt", sep = "\t")
protein_list <- protein_list %>% 
  filter(panel%in%c("Cardiometabolic","Inflammation","Neurology" ,"Oncology"))
protein_list <- protein_list$mrc_olink.id

## merge phenotypes and proteins
ukbb <- merge(ol.npx, pheno, by="eid_20361")

## select random subset or consortium selected cases only
ukbb <-
  ukbb %>% 
  filter(randomized_baseline==T|!!as.name(dz)==1)

## inverse rank normal transform npx
ukbb[,protein_list] <- sapply(ukbb[,protein_list], RankNorm)

## exclude prevalent cases and incident cases up to 6 months
## restrict to incident cases within the first X years
ukbb <- ukbb %>% 
  mutate(fol = difftime(!!as.name(dz.date),date_baseline_assessment, units = "weeks")/52.25) %>% 
  filter((fol>0.5&fol<inc.yrs)|is.na(!!as.name(dz.date)))

## censor (including censoring by death)
death.data <- read.delim("../../01_inc_dz_phenos/data_input/Death_dates_proteomic_participants.txt",sep = "\t")
death.data$date_of_death <- as.Date(death.data$date_of_death, format = "%Y-%m-%d")

ukbb <- merge(ukbb, death.data[,c("eid_20361","date_of_death")], by="eid_20361", all.x=T)

if(inc.yrs==10){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2020-12-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
} else if(inc.yrs==5){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2016-05-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
}else{stop("limit to incident casese within X years not defined")}

if(inc.yrs==10){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2020-12-31"),
                         difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                        ukbb$fol)
} else if(inc.yrs==5){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2016-05-31"),
                          difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                          ukbb$fol)
}else{stop("limit to incident casese within X years not defined")}

ukbb$fol <- ifelse(is.na(ukbb$fol_d), ukbb$fol, ukbb$fol_d)

## for sex specific outcome exclude males/females
if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_pri_breast",
           "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical",
           "atlas_inc_endometrial_hyper","atlas_inc_endometriosis","atlas_inc_female_genital_prolapse",
           "atlas_inc_leiomyoma","atlas_inc_menorrhagia",
           "atlas_inc_pid","atlas_inc_pmb")){
  # ukbb <- ukbb %>% filter(ge_sex=="female")
  ukbb <- ukbb %>% filter(ge_sex==0)
} else if(dz%in%c("atlas_inc_pri_prost","atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  # ukbb <- ukbb %>% filter(ge_sex=="male")
  ukbb <- ukbb %>% filter(ge_sex==1)
} else{ukbb <- ukbb}

## recode sex to numeric
# ukbb$ge_sex <- ifelse(ukbb$ge_sex=="male",1,
#                       ifelse(ukbb$ge_sex=="female",0,NA))

n.cases <- sum(ukbb[,dz]==1)

## add noise variables
rand.vars <- lapply(1:10, function(x){
  tmp <- runif(nrow(ukbb), min = min(ukbb$NPPB_20049), 
               max = max(ukbb$NPPB_20049))
  return(tmp)
})
rand.vars <- as.data.frame(do.call(cbind, rand.vars))
colnames(rand.vars) <- paste0("rand_var_",1:10)

ukbb <- cbind(ukbb, rand.vars)

####----    ML framework    ----####

## split into train/optimize/test or train/test depending on case number
ukbb[,dz] <- as.factor(ukbb[,dz])

if(n.cases/4>200){
  ## feature selection set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.5, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.opti <- ukbb[-trainIndex,]
  
  ## optimization and testing
  set.seed(3456) ## for reproducibility
  trainIndex2 <- createDataPartition(u.opti[,dz], p=0.5, list = F, times = 1)
  u.test <- u.opti[trainIndex2,]
  u.opti <- u.opti[-trainIndex2,]
  }else {
  ## feature selection set and test set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.7, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.test <- ukbb[-trainIndex,]
  }
  
## draw 200 random subsamples of train set
jj <- lapply(1:200, function(x) sample(nrow(u.train), round(nrow(u.train)*0.50), replace=F))

##-- feature selection  --##
las.morb <- mclapply(1:length(jj), function(i){
  print(i) ## just for test purposes
  tmp <- train(u.train[jj[[i]],c(protein_list,colnames(rand.vars))],
               u.train[jj[[i]],dz],
               metric = "Kappa",
               method="glmnet",
               family="binomial",
               tuneGrid=as.data.frame(expand_grid(
                 alpha=1,lambda=10^-seq(4,.25,-.5))),
               trControl= trainControl(method = "repeatedcv",
                                       number = 3, repeats = 5,
                                       sampling = "rose",
                                       allowParallel = T)
  )
  cf <- as.matrix(coef(tmp$finalModel,
                       s = tmp$finalModel$lambdaOpt))
  gc()
  # return(list("cf" = cf , "cv.tune" = tmp))
  return(cf)
}, mc.cores = 12, mc.allow.recursive = T)

save(las.morb, file=paste0("../data_output/feature_selection/fs_1536_",dz,"_",inc.yrs,"yrs",".RData"))

## generate feature selection ranking 
# load(paste0("../data_output/feature_selection/fs_",dz,"_",inc.yrs,"yrs",".RData"))
p.select <- matrix(,nrow = nrow(las.morb[[1]]), ncol = length(las.morb))

for(i in 1:length(las.morb)){
  tryCatch({
    tmp <- las.morb[[i]][,1]
    if(is.null(tmp)){
      return(tmp <- rep(0,nrow(p.select)))
    }else {
      NULL
    }
  }, error=function(e){
    NULL
  })
  p.select[,i] <- tmp
}
row.names(p.select) <- names(las.morb[[1]][,1])

## count matrix
p.select <-
  p.select %>% 
  as_tibble(rownames = "mrc_olink.id") %>% 
  filter(mrc_olink.id!="(Intercept)") %>% 
  mutate(select = abs(rowSums(across(where(is.numeric))))) %>% 
  arrange(desc(select)) 
  
p.select <- as.data.frame(p.select)

p.select$select.perc <- (p.select$select)/max(p.select$select)

####---- optimization and performance testing ----####

## survival objects
surv.train <- Surv(u.train$fol, u.train[,dz])
## optimization 
if(n.cases/4>200){
  surv.opt <- Surv(u.opti$fol, u.opti[,dz])
}else{surv.opt <- Surv(u.train$fol, u.train[,dz])}
## testing
surv.test <- Surv(u.test$fol, u.test[,dz])

## function to run optimization and testings
source("../../../02_inc_dz_prediction/bin/functions/coxnet_optim.R")
source("../../../02_inc_dz_prediction/bin/functions/coxnet_ridge_optim_boot.R")

## px information model basic 
clin.predictors <- c("ge_age","ge_sex","bmi","ethnicity","smoke_never",
                     "smoke_previous","smoke_current","alcohol_daily",
                     "alcohol_4_week","alcohol_2_week","alcohol_3_month","alcohol_occassion",
                     "alcohol_never")

## define additional family history for specific diseases
if(!dz%in%c("atlas_inc_diabetes_t2","atlas_inc_pri_prost", "atlas_inc_copd","atlas_inc_isch_stroke",
            "atlas_inc_af","atlas_inc_av_block_3","atlas_inc_chd_nos","atlas_inc_hf", 
            "atlas_inc_peripheral_arterial_disease","atlas_inc_parkinsons","atlas_inc_pri_lung",
            "atlas_inc_pri_bowel", "atlas_inc_pri_breast","atlas_inc_pri_ovarian","atlas_inc_pri_uterine",
            "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
            "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
            "atlas_inc_hydrocele","atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb",
            "atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu")){
  clin.predictors <- clin.predictors
} else if(dz=="atlas_inc_diabetes_t2") {
  clin.predictors <- c(clin.predictors,"father_diabetes","mother_diabetes")
} else if(dz=="atlas_inc_pri_prost") {
  clin.predictors <- c(clin.predictors[-2],"father_prostate_cancer")
} else if(dz=="atlas_inc_copd") {
  clin.predictors <- c(clin.predictors,"father_COPD","mother_COPD")
} else if(dz%in%c("atlas_inc_isch_stroke","atlas_inc_af",
                  "atlas_inc_chd_nos","atlas_inc_hf","atlas_inc_peripheral_arterial_disease")) {
  clin.predictors <- c(clin.predictors,"father_heart_disease","mother_heart_disease")
} else if(dz=="atlas_inc_parkinsons") {
  clin.predictors <- c(clin.predictors,"father_parkinsons","mother_parkinsons")
} else if(dz%in%c("atlas_inc_pri_lung","atlas_inc_sec_lung")) {
  clin.predictors <- c(clin.predictors,"father_lung_cancer","mother_lung_cancer")
} else if(dz=="atlas_inc_pri_bowel") {
  clin.predictors <- c(clin.predictors,"father_bowel_cancer","mother_bowel_cancer")
} else if(dz=="atlas_inc_pri_breast") {
  clin.predictors <- c(clin.predictors[-2],"mother_breast_cancer")
} else if(dz=="atlas_inc_hypertension") {
  clin.predictors <- c(clin.predictors[-2],"mother_high_bp","father_high_bp")
} else if(dz=="atlas_inc_depression") {
  clin.predictors <- c(clin.predictors[-2],"mother_depression","father_depression")
} else if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_benign_ovary",
                  "atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
                  "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
                  "atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb")) {
  clin.predictors <- c(clin.predictors[-2])
} else if(dz%in%c("atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  clin.predictors <- c(clin.predictors[-2])
}


boot.samples <- lapply(1:1000, function(x) sample(nrow(u.test), round(nrow(u.test)*1), replace = T))

## mode px info model optimization and testing
if(n.cases/4>200){
  res.clin <- coxnet.optim.r(Train.data = u.opti,
                           pred.vec = clin.predictors,
                           Train.surv.data = surv.opt,
                           times = 1000, 
                           boot.samples = boot.samples,
                           Test.data = u.test,
                           Test.surv.data = surv.test)
} else{
  res.clin <- coxnet.optim.r(Train.data = u.train,
                           pred.vec = clin.predictors,
                           Train.surv.data = surv.opt,
                           times = 1000, 
                           boot.samples = boot.samples,
                           Test.data = u.test,
                           Test.surv.data = surv.test)
  
}
save(res.clin, file = paste0("../data_output/px_info_models/px_info_res_1536_",dz,"_",inc.yrs,"yrs",".RData"))

## selection threshold based on noise variables
sel.t <- p.select$select.perc[grep("rand", p.select$mrc_olink.id)[1]]

prots.opti <- p.select %>% 
  filter(select.perc>sel.t) %>% 
  top_n(20, select.perc)
prots.opti <- prots.opti$mrc_olink.id
if(length(prots.opti)==0){
  prots.opti <- p.select[p.select$select.per==sel.t]
  prots.opti <- prots.opti$mrc_olink.id
}

# run px info model + proteins selected above noise
if(n.cases/4>200){
    res.clin.prots.all <- coxnet.optim.r(Train.data = u.opti,
                                   pred.vec = c(clin.predictors,prots.opti),
                                   Train.surv.data = surv.opt,
                                   times = 1000, 
                                   boot.samples = boot.samples,
                                   Test.data = u.test,
                                   Test.surv.data = surv.test)
    res.prots.all <- coxnet.optim.r(Train.data = u.opti,
                                    pred.vec = prots.opti,
                                    Train.surv.data = surv.opt,
                                    times = 1000, 
                                    boot.samples = boot.samples,
                                    Test.data = u.test,
                                    Test.surv.data = surv.test)
    } else{
      res.clin.prots.all <- coxnet.optim.r(Train.data = u.train,
                                   pred.vec = c(clin.predictors,prots.opti),
                                   Train.surv.data = surv.opt,
                                   times = 1000, 
                                   boot.samples = boot.samples,
                                   Test.data = u.test,
                                   Test.surv.data = surv.test)
      res.prots.all <- coxnet.optim.r(Train.data = u.train,
                                      pred.vec = prots.opti,
                                      Train.surv.data = surv.opt,
                                      times = 1000, 
                                      boot.samples = boot.samples,
                                      Test.data = u.test,
                                      Test.surv.data = surv.test)
    }
save(res.clin.prots.all, file = paste0("../data_output/px_info_protein_models/px_info_all_prots_above_noise_res_1536_",dz,"_",inc.yrs,"yrs",".RData"))
save(res.prots.all, file = paste0("../data_output/protein_models/prot_all_res_1536_",dz,"_",inc.yrs,"yrs",".RData"))

# run px info model + proteins selected above noise in train set
if(n.cases/4>200){
  res.clin.prots.all <- coxnet.optim.r(Train.data = u.train,
                                     pred.vec = c(clin.predictors,prots.opti),
                                     Train.surv.data = surv.train,
                                     times = 1000, 
                                     boot.samples = boot.samples,
                                     Test.data = u.test,
                                     Test.surv.data = surv.test)
} else {res.clin.prots.all <- res.clin.prots.all}

## proteins top 5 and 10 only based on weights*selection_percentage
opti.prots <- res.clin.prots.all$opt.coefficients[grep("_[0-9]{5}",row.names(res.clin.prots.all$opt.coefficients)),]
opti.prots <- opti.prots %>%
  as_tibble(rownames = "prot") %>%
  mutate(abs.coef = abs(value)) %>% 
  filter(abs.coef!=0)

## merge with selection percentage
opti.prots <- merge(opti.prots, p.select[,c("mrc_olink.id","select","select.perc")],
                    by.x="prot", by.y="mrc_olink.id")

if(nrow(opti.prots)>=10){
  opti.prots10 <- opti.prots %>%
    mutate(sel.opti.score=abs.coef*select) %>%
    top_n(10, sel.opti.score)
  prots.opti10 <- opti.prots10$prot
  
  ## run px info model + 10 top proteins based on weightsXselection.perc
  if(n.cases/4>200){
    res.clin.prots10 <- coxnet.optim.r(Train.data = u.opti,
                                     pred.vec = c(clin.predictors,prots.opti10),
                                     Train.surv.data = surv.opt,
                                     times = 1000, 
                                     boot.samples = boot.samples,
                                     Test.data = u.test,
                                     Test.surv.data = surv.test)
    res.prots10 <- coxnet.optim.r(Train.data = u.opti,
                                  pred.vec = prots.opti10,
                                  Train.surv.data = surv.opt,
                                  times = 1000, 
                                  boot.samples = boot.samples,
                                  Test.data = u.test,
                                  Test.surv.data = surv.test)
  } else{
    res.clin.prots10 <- coxnet.optim.r(Train.data = u.train,
                                     pred.vec = c(clin.predictors,prots.opti10),
                                     Train.surv.data = surv.opt,
                                     times = 1000, 
                                     boot.samples = boot.samples,
                                     Test.data = u.test,
                                     Test.surv.data = surv.test)
    res.prots10 <- coxnet.optim.r(Train.data = u.train,
                                  pred.vec = prots.opti10,
                                  Train.surv.data = surv.opt,
                                  times = 1000, 
                                  boot.samples = boot.samples,
                                  Test.data = u.test,
                                  Test.surv.data = surv.test)
    
  }
  save(res.clin.prots10, file = paste0("../data_output/px_info_protein_models/px_info_top10_prots_res_1536_",dz,"_",inc.yrs,"yrs",".RData"))
  save(res.prots10, file = paste0("../data_output/protein_models/prot_top10_res_1536_",dz,"_",inc.yrs,"yrs",".RData"))
  
}

## run px info model + 5 top proteins based on weightsXselection.perc
if(nrow(opti.prots)>=5){
  opti.prots <- opti.prots %>%
  mutate(sel.opti.score=abs.coef*select.perc) %>%
  top_n(5, sel.opti.score)
  
  prots.opti <- opti.prots$prot
} else{prots.opti <- opti.prots$prot}


if(n.cases/4>200){
  res.clin.prots <- coxnet.optim.r(Train.data = u.opti,
                                 pred.vec = c(clin.predictors,prots.opti),
                                 Train.surv.data = surv.opt,
                                 times = 1000, 
                                 boot.samples = boot.samples,
                                 Test.data = u.test,
                                 Test.surv.data = surv.test)
} else{
  res.clin.prots <- coxnet.optim.r(Train.data = u.train,
                                 pred.vec = c(clin.predictors,prots.opti),
                                 Train.surv.data = surv.opt,
                                 times = 1000, 
                                 boot.samples = boot.samples,
                                 Test.data = u.test,
                                 Test.surv.data = surv.test)
  
}
save(res.clin.prots, file = paste0("../data_output/px_info_protein_models/px_info_top5_prots_res_1536_",dz,"_",inc.yrs,"yrs",".RData"))

## optimization on top 5 proteins
if(n.cases/4>200){
  res.prot <- coxnet.optim.r(Train.data = u.opti,
                           pred.vec = prots.opti,
                           Train.surv.data = surv.opt,
                           times = 1000, 
                           boot.samples = boot.samples, 
                           Test.data = u.test,
                           Test.surv.data = surv.test)
} else{
  res.prot <- coxnet.optim.r(Train.data = u.train,
                           pred.vec = prots.opti,
                           Train.surv.data = surv.opt,
                           times = 1000, 
                           boot.samples = boot.samples,
                           Test.data = u.test,
                           Test.surv.data = surv.test)
  
}
## add feature selection to the results
res.prot$feature.selection <- p.select

## add case numbers for training and testing sets
if((n.cases/4)>200){
  case.num <- data.frame(train.cases = sum(u.train[,dz]==1),
                         opti.cases = sum(u.opti[,dz]==1),
                         test.cases = sum(u.test[,dz]==1),
                         train.controls = sum(u.train[,dz]!=1),
                         opti.controls = sum(u.opti[,dz]!=1),
                         test.controls = sum(u.test[,dz]!=1))
} else{
  case.num <- data.frame(train.cases = sum(u.train[,dz]==1),
                         opti.cases = NA,
                         test.cases = sum(u.test[,dz]==1),
                         train.controls = sum(u.train[,dz]!=1),
                         opti.controls = NA,
                         test.controls = sum(u.test[,dz]!=1))
}
res.prot$case.num <- case.num
## save results 
save(res.prot, file = paste0("../data_output/protein_models/prot_top5_res_1536_",dz,"_",inc.yrs,"yrs",".RData"))

###################################################
#### Script Name: ML_pipeline_protein_pxinfo_prediction_pred
#### Author: Julia Carrasco Zanini Sanchez
#### Description: script to ML framework for incident diseases in atlas 308 
#### Date : 23/03/2023

rm(list = ls())

## load packages 
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(RNOmni)
# library(glue)
library(DBI)
library(dbplyr)
library(magrittr)
# library(broom)
library(survival)
library(data.table)
library(glue)
library(parallel)
library(caret)
library(glmnet)
library(ROSE)
# library(odbc)

## define disease and restrict to incident cases within X years
args <- commandArgs(trailingOnly = T)

dz <- as.character(args[1])
dz.date <- as.character(args[2])
inc.yrs <- as.numeric(args[3])

paste0("Prediction:", dz, "; ", dz.date, " at ", inc.yrs, "years")

## set up working directory 
setwd("/home/ukbb_inc_dz_prediction/Explore_1536_Expansion_analyses/02_inc_dz_prediction/bin/") 

## read in data 
atlas308 <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_atlas308_annotations.txt", data.table = F)
cov <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_covariate_annotations.txt", data.table = F)

####----  Data Cleaning and wrangling   ----####
atlas308 <- as.data.frame(atlas308 %>% 
  filter(eid_20361%in%cov$eid_20361) %>% 
  pivot_wider(names_from = variable, values_from =value))

## phenotypes missing for patient-derived info model: bmi, smoking, alcohol, physical activity, family history 
pheno <- merge(cov[,-c(11:14,17)], atlas308, by="eid_20361")

## read in px info phenotypes
px.info <- read.delim("../../01_inc_dz_phenos/data_input/Patient_model_variables.txt",sep = "\t")
pheno <- merge(pheno, px.info, by="eid_20361")

## rename to olds names to match the rest fo the script
colnames(pheno)[which(colnames(pheno)=="assess_date")] <- c("date_baseline_assessment")

## process ethnicity
pheno$ethnicity <- ifelse(is.na(pheno$ge_eur),2,pheno$ge_eur)

## get protein data and metadata
ol.npx <- fread("../../00_protein_imputation/data_input/Imputed_1536_Expansion_NPX_proteins.txt", data.table = F)
ol.npx <- filter(ol.npx, exclusion.50==0 )

## exclude participants with missing age, sex
ol.npx <- ol.npx[complete.cases(ol.npx[,c("ge_sex", "ge_age")]),]

protein_list <- read.delim("../../00_protein_imputation/data_input/protein_labels_by_panel.txt", sep = "\t")
protein_list <- protein_list %>% 
  filter(panel%in%c("Cardiometabolic","Inflammation","Neurology" ,"Oncology"))
protein_list <- protein_list$mrc_olink.id

## merge phenotypes and proteins
ukbb <- merge(ol.npx, pheno, by="eid_20361")
## randomised baseline
ukbb <- ukbb %>% 
  filter(randomized_baseline==T)
## inverse rank normal transform npx
ukbb[,protein_list] <- sapply(ukbb[,protein_list], RankNorm)

## exclude prevalent cases and incident cases up to 6 months
## restrict to incident cases within the first X years
ukbb <- ukbb %>% 
  mutate(fol = difftime(!!as.name(dz.date),date_baseline_assessment, units = "weeks")/52.25) %>% 
  filter((fol>0.5&fol<inc.yrs)|is.na(!!as.name(dz.date)))

## censor (including censoring by death)
death.data <- read.delim("../../01_inc_dz_phenos/data_input/Death_dates_proteomic_participants.txt",sep = "\t")
death.data$date_of_death <- as.Date(death.data$date_of_death, format = "%Y-%m-%d")

ukbb <- merge(ukbb, death.data[,c("eid_20361","date_of_death")], by="eid_20361", all.x=T)

if(inc.yrs==10){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2020-12-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
} else if(inc.yrs==5){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2016-05-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
}else{stop("limit to incident casese within X years not defined")}

if(inc.yrs==10){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2020-12-31"),
                         difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                        ukbb$fol)
} else if(inc.yrs==5){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2016-05-31"),
                          difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                          ukbb$fol)
}else{stop("limit to incident casese within X years not defined")}

ukbb$fol <- ifelse(is.na(ukbb$fol_d), ukbb$fol, ukbb$fol_d)

## for sex specific outcome exclude males/females
if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_pri_breast",
           "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical",
           "atlas_inc_endometrial_hyper","atlas_inc_endometriosis","atlas_inc_female_genital_prolapse",
           "atlas_inc_leiomyoma","atlas_inc_menorrhagia",
           "atlas_inc_pid","atlas_inc_pmb")){
  # ukbb <- ukbb %>% filter(ge_sex=="female")
  ukbb <- ukbb %>% filter(ge_sex==0)
} else if(dz%in%c("atlas_inc_pri_prost","atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  # ukbb <- ukbb %>% filter(ge_sex=="male")
  ukbb <- ukbb %>% filter(ge_sex==1)
} else{ukbb <- ukbb}

## recode sex to numeric
# ukbb$ge_sex <- ifelse(ukbb$ge_sex=="male",1,
#                       ifelse(ukbb$ge_sex=="female",0,NA))

n.cases <- sum(ukbb[,dz]==1)

## add noise variables
rand.vars <- lapply(1:10, function(x){
  tmp <- runif(nrow(ukbb), min = min(ukbb$NPPB_20049), 
               max = max(ukbb$NPPB_20049))
  return(tmp)
})
rand.vars <- as.data.frame(do.call(cbind, rand.vars))
colnames(rand.vars) <- paste0("rand_var_",1:10)

ukbb <- cbind(ukbb, rand.vars)

####----    ML framework    ----####

## split into train/optimize/test or train/test depending on case number
ukbb[,dz] <- as.factor(ukbb[,dz])

if(n.cases/4>200){
  ## feature selection set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.5, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.opti <- ukbb[-trainIndex,]
  
  ## optimization and testing
  set.seed(3456) ## for reproducibility
  trainIndex2 <- createDataPartition(u.opti[,dz], p=0.5, list = F, times = 1)
  u.test <- u.opti[trainIndex2,]
  u.opti <- u.opti[-trainIndex2,]
  }else {
  ## feature selection set and test set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.7, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.test <- ukbb[-trainIndex,]
  }
  
## draw 200 random subsamples of train set
jj <- lapply(1:200, function(x) sample(nrow(u.train), round(nrow(u.train)*0.50), replace=F))

##-- feature selection  --##
las.morb <- mclapply(1:length(jj), function(i){
  print(i) ## just for test purposes
  tmp <- train(u.train[jj[[i]],c(protein_list,colnames(rand.vars))],
               u.train[jj[[i]],dz],
               metric = "Kappa",
               method="glmnet",
               family="binomial",
               tuneGrid=as.data.frame(expand_grid(
                 alpha=1,lambda=10^-seq(4,.25,-.5))),
               trControl= trainControl(method = "repeatedcv",
                                       number = 3, repeats = 5,
                                       sampling = "rose",
                                       allowParallel = T)
  )
  cf <- as.matrix(coef(tmp$finalModel,
                       s = tmp$finalModel$lambdaOpt))
  gc()
  # return(list("cf" = cf , "cv.tune" = tmp))
  return(cf)
}, mc.cores = 12, mc.allow.recursive = T)

save(las.morb, file=paste0("../data_output/feature_selection/fs_1536_",dz,"_",inc.yrs,"yrs",".RData"))

## generate feature selection ranking 
# load(paste0("../data_output/feature_selection/fs_",dz,"_",inc.yrs,"yrs",".RData"))
p.select <- matrix(,nrow = nrow(las.morb[[1]]), ncol = length(las.morb))

for(i in 1:length(las.morb)){
  tryCatch({
    tmp <- las.morb[[i]][,1]
    if(is.null(tmp)){
      return(tmp <- rep(0,nrow(p.select)))
    }else {
      NULL
    }
  }, error=function(e){
    NULL
  })
  p.select[,i] <- tmp
}
row.names(p.select) <- names(las.morb[[1]][,1])

## count matrix
p.select <-
  p.select %>% 
  as_tibble(rownames = "mrc_olink.id") %>% 
  filter(mrc_olink.id!="(Intercept)") %>% 
  mutate(select = abs(rowSums(across(where(is.numeric))))) %>% 
  arrange(desc(select)) 
  
p.select <- as.data.frame(p.select)

p.select$select.perc <- (p.select$select)/max(p.select$select)

####---- optimization and performance testing ----####

## survival objects
surv.train <- Surv(u.train$fol, u.train[,dz])
## optimization 
if(n.cases/4>200){
  surv.opt <- Surv(u.opti$fol, u.opti[,dz])
}else{surv.opt <- Surv(u.train$fol, u.train[,dz])}
## testing
surv.test <- Surv(u.test$fol, u.test[,dz])

## function to run optimization and testings
source("../../../02_inc_dz_prediction/bin/functions/coxnet_ridge_optim_boot.R")

## px information model basic 
clin.predictors <- c("ge_age","ge_sex","bmi","ethnicity","smoke_never",
                     "smoke_previous","smoke_current","alcohol_daily",
                     "alcohol_4_week","alcohol_2_week","alcohol_3_month","alcohol_occassion",
                     "alcohol_never")

## define additional family history for specific diseases
if(!dz%in%c("atlas_inc_diabetes_t2","atlas_inc_pri_prost", "atlas_inc_copd","atlas_inc_isch_stroke",
            "atlas_inc_af","atlas_inc_av_block_3","atlas_inc_chd_nos","atlas_inc_hf", 
            "atlas_inc_peripheral_arterial_disease","atlas_inc_parkinsons","atlas_inc_pri_lung",
            "atlas_inc_pri_bowel", "atlas_inc_pri_breast","atlas_inc_pri_ovarian","atlas_inc_pri_uterine",
            "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
            "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
            "atlas_inc_hydrocele","atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb",
            "atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu")){
  clin.predictors <- clin.predictors
} else if(dz=="atlas_inc_diabetes_t2") {
  clin.predictors <- c(clin.predictors,"father_diabetes","mother_diabetes")
} else if(dz=="atlas_inc_pri_prost") {
  clin.predictors <- c(clin.predictors[-2],"father_prostate_cancer")
} else if(dz=="atlas_inc_copd") {
  clin.predictors <- c(clin.predictors,"father_COPD","mother_COPD")
} else if(dz%in%c("atlas_inc_isch_stroke","atlas_inc_af",
                  "atlas_inc_chd_nos","atlas_inc_hf","atlas_inc_peripheral_arterial_disease")) {
  clin.predictors <- c(clin.predictors,"father_heart_disease","mother_heart_disease")
} else if(dz=="atlas_inc_parkinsons") {
  clin.predictors <- c(clin.predictors,"father_parkinsons","mother_parkinsons")
} else if(dz%in%c("atlas_inc_pri_lung","atlas_inc_sec_lung")) {
  clin.predictors <- c(clin.predictors,"father_lung_cancer","mother_lung_cancer")
} else if(dz=="atlas_inc_pri_bowel") {
  clin.predictors <- c(clin.predictors,"father_bowel_cancer","mother_bowel_cancer")
} else if(dz=="atlas_inc_pri_breast") {
  clin.predictors <- c(clin.predictors[-2],"mother_breast_cancer")
} else if(dz=="atlas_inc_hypertension") {
  clin.predictors <- c(clin.predictors[-2],"mother_high_bp","father_high_bp")
} else if(dz=="atlas_inc_depression") {
  clin.predictors <- c(clin.predictors[-2],"mother_depression","father_depression")
} else if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_benign_ovary",
                  "atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
                  "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
                  "atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb")) {
  clin.predictors <- c(clin.predictors[-2])
} else if(dz%in%c("atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  clin.predictors <- c(clin.predictors[-2])
}


boot.samples <- lapply(1:1000, function(x) sample(nrow(u.test), round(nrow(u.test)*1), replace = T))

## mode px info model optimization and testing
if(n.cases/4>200){
  res.clin <- coxnet.optim.r(Train.data = u.opti,
                           pred.vec = clin.predictors,
                           Train.surv.data = surv.opt,
                           times = 1000, 
                           boot.samples = boot.samples,
                           Test.data = u.test,
                           Test.surv.data = surv.test)
} else{
  res.clin <- coxnet.optim.r(Train.data = u.train,
                           pred.vec = clin.predictors,
                           Train.surv.data = surv.opt,
                           times = 1000, 
                           boot.samples = boot.samples,
                           Test.data = u.test,
                           Test.surv.data = surv.test)
  
}
save(res.clin, file = paste0("../data_output/px_info_models/px_info_res_1536_",dz,"_",inc.yrs,"yrs",".RData"))

## selection threshold based on noise variables
sel.t <- p.select$select.perc[grep("rand", p.select$mrc_olink.id)[1]]

prots.opti <- p.select %>% 
  filter(select.perc>sel.t) %>% 
  top_n(20, select.perc)
prots.opti <- prots.opti$mrc_olink.id
if(length(prots.opti)==0){
  prots.opti <- p.select[p.select$select.per==sel.t]
  prots.opti <- prots.opti$mrc_olink.id
}

# run px info model + proteins selected above noise
if(n.cases/4>200){
    res.clin.prots.all <- coxnet.optim.r(Train.data = u.opti,
                                   pred.vec = c(clin.predictors,prots.opti),
                                   Train.surv.data = surv.opt,
                                   times = 1000, 
                                   boot.samples = boot.samples,
                                   Test.data = u.test,
                                   Test.surv.data = surv.test)
    res.prots.all <- coxnet.optim.r(Train.data = u.opti,
                                    pred.vec = prots.opti,
                                    Train.surv.data = surv.opt,
                                    times = 1000, 
                                    boot.samples = boot.samples,
                                    Test.data = u.test,
                                    Test.surv.data = surv.test)
    } else{
      res.clin.prots.all <- coxnet.optim.r(Train.data = u.train,
                                   pred.vec = c(clin.predictors,prots.opti),
                                   Train.surv.data = surv.opt,
                                   times = 1000, 
                                   boot.samples = boot.samples,
                                   Test.data = u.test,
                                   Test.surv.data = surv.test)
      res.prots.all <- coxnet.optim.r(Train.data = u.train,
                                      pred.vec = prots.opti,
                                      Train.surv.data = surv.opt,
                                      times = 1000, 
                                      boot.samples = boot.samples,
                                      Test.data = u.test,
                                      Test.surv.data = surv.test)
    }
save(res.clin.prots.all, file = paste0("../data_output/px_info_protein_models/px_info_all_prots_above_noise_res_1536_",dz,"_",inc.yrs,"yrs",".RData"))
save(res.prots.all, file = paste0("../data_output/protein_models/prot_all_res_1536_",dz,"_",inc.yrs,"yrs",".RData"))

# run px info model + proteins selected above noise in train set
if(n.cases/4>200){
  res.clin.prots.all <- coxnet.optim.r(Train.data = u.train,
                                     pred.vec = c(clin.predictors,prots.opti),
                                     Train.surv.data = surv.train,
                                     times = 1000, 
                                     boot.samples = boot.samples,
                                     Test.data = u.test,
                                     Test.surv.data = surv.test)
} else {res.clin.prots.all <- res.clin.prots.all}

## proteins top 5 and 10 only based on weights*selection_percentage
opti.prots <- res.clin.prots.all$opt.coefficients[grep("_[0-9]{5}",row.names(res.clin.prots.all$opt.coefficients)),]
opti.prots <- opti.prots %>%
  as_tibble(rownames = "prot") %>%
  mutate(abs.coef = abs(value)) %>% 
  filter(abs.coef!=0)

## merge with selection percentage
opti.prots <- merge(opti.prots, p.select[,c("mrc_olink.id","select","select.perc")],
                    by.x="prot", by.y="mrc_olink.id")

if(nrow(opti.prots)>=10){
  opti.prots10 <- opti.prots %>%
    mutate(sel.opti.score=abs.coef*select) %>%
    top_n(10, sel.opti.score)
  prots.opti10 <- opti.prots10$prot
  
  ## run px info model + 10 top proteins based on weightsXselection.perc
  if(n.cases/4>200){
    res.clin.prots10 <- coxnet.optim.r(Train.data = u.opti,
                                     pred.vec = c(clin.predictors,prots.opti10),
                                     Train.surv.data = surv.opt,
                                     times = 1000, 
                                     boot.samples = boot.samples,
                                     Test.data = u.test,
                                     Test.surv.data = surv.test)
    res.prots10 <- coxnet.optim.r(Train.data = u.opti,
                                  pred.vec = prots.opti10,
                                  Train.surv.data = surv.opt,
                                  times = 1000, 
                                  boot.samples = boot.samples,
                                  Test.data = u.test,
                                  Test.surv.data = surv.test)
  } else{
    res.clin.prots10 <- coxnet.optim.r(Train.data = u.train,
                                     pred.vec = c(clin.predictors,prots.opti10),
                                     Train.surv.data = surv.opt,
                                     times = 1000, 
                                     boot.samples = boot.samples,
                                     Test.data = u.test,
                                     Test.surv.data = surv.test)
    res.prots10 <- coxnet.optim.r(Train.data = u.train,
                                  pred.vec = prots.opti10,
                                  Train.surv.data = surv.opt,
                                  times = 1000, 
                                  boot.samples = boot.samples,
                                  Test.data = u.test,
                                  Test.surv.data = surv.test)
    
  }
  save(res.clin.prots10, file = paste0("../data_output/px_info_protein_models/px_info_top10_prots_res_1536_",dz,"_",inc.yrs,"yrs",".RData"))
  save(res.prots10, file = paste0("../data_output/protein_models/prot_top10_res_1536_",dz,"_",inc.yrs,"yrs",".RData"))
  
}

## run px info model + 5 top proteins based on weightsXselection.perc
if(nrow(opti.prots)>=5){
  opti.prots <- opti.prots %>%
  mutate(sel.opti.score=abs.coef*select.perc) %>%
  top_n(5, sel.opti.score)
  
  prots.opti <- opti.prots$prot
} else{prots.opti <- opti.prots$prot}


if(n.cases/4>200){
  res.clin.prots <- coxnet.optim.r(Train.data = u.opti,
                                 pred.vec = c(clin.predictors,prots.opti),
                                 Train.surv.data = surv.opt,
                                 times = 1000, 
                                 boot.samples = boot.samples,
                                 Test.data = u.test,
                                 Test.surv.data = surv.test)
} else{
  res.clin.prots <- coxnet.optim.r(Train.data = u.train,
                                 pred.vec = c(clin.predictors,prots.opti),
                                 Train.surv.data = surv.opt,
                                 times = 1000, 
                                 boot.samples = boot.samples,
                                 Test.data = u.test,
                                 Test.surv.data = surv.test)
  
}
save(res.clin.prots, file = paste0("../data_output/px_info_protein_models/px_info_top5_prots_res_1536_",dz,"_",inc.yrs,"yrs",".RData"))

## optimization on top 5 proteins
if(n.cases/4>200){
  res.prot <- coxnet.optim.r(Train.data = u.opti,
                           pred.vec = prots.opti,
                           Train.surv.data = surv.opt,
                           times = 1000, 
                           boot.samples = boot.samples, 
                           Test.data = u.test,
                           Test.surv.data = surv.test)
} else{
  res.prot <- coxnet.optim.r(Train.data = u.train,
                           pred.vec = prots.opti,
                           Train.surv.data = surv.opt,
                           times = 1000, 
                           boot.samples = boot.samples,
                           Test.data = u.test,
                           Test.surv.data = surv.test)
  
}
## add feature selection to the results
res.prot$feature.selection <- p.select

## add case numbers for training and testing sets
if((n.cases/4)>200){
  case.num <- data.frame(train.cases = sum(u.train[,dz]==1),
                         opti.cases = sum(u.opti[,dz]==1),
                         test.cases = sum(u.test[,dz]==1),
                         train.controls = sum(u.train[,dz]!=1),
                         opti.controls = sum(u.opti[,dz]!=1),
                         test.controls = sum(u.test[,dz]!=1))
} else{
  case.num <- data.frame(train.cases = sum(u.train[,dz]==1),
                         opti.cases = NA,
                         test.cases = sum(u.test[,dz]==1),
                         train.controls = sum(u.train[,dz]!=1),
                         opti.controls = NA,
                         test.controls = sum(u.test[,dz]!=1))
}
res.prot$case.num <- case.num
## save results 
save(res.prot, file = paste0("../data_output/protein_models/prot_top5_res_1536_",dz,"_",inc.yrs,"yrs",".RData"))

###################################################
#### Script Name: ML_pipeline_6m_protein_pxinfo_prediction_pred_consortcases
#### Author: Julia Carrasco Zanini Sanchez
#### Description: script to ML framework for incident diseases in atlas 308 
#### Date : 23/03/2023

rm(list = ls())

## load packages 
library(dplyr)
library(tidyr)
library(RNOmni)
# library(glue)
library(DBI)
library(dbplyr)
library(magrittr)
# library(broom)
library(survival)
library(data.table)
library(glue)
library(parallel)
library(caret)
library(glmnet)
library(ROSE)

## define disease and restrict to incident cases within X years
args <- commandArgs(trailingOnly = T)

dz <- as.character(args[1])
dz.date <- as.character(args[2])
inc.yrs <- as.numeric(args[3])

paste0("Prediction:", dz, "; ", dz.date, " at ", inc.yrs, "years")

## set up working directory 
setwd("/home/ukbb_inc_dz_prediction/Explore_1536_Expansion_analyses/02_inc_dz_prediction/bin/") 

## read in data 
atlas308 <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_atlas308_annotations.txt", data.table = F)
cov <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_covariate_annotations.txt", data.table = F)

####----  Data Cleaning and wrangling   ----####
atlas308 <- as.data.frame(atlas308 %>% 
  filter(eid_20361%in%cov$eid_20361) %>% 
  pivot_wider(names_from = variable, values_from =value))

## phenotypes missing for patient-derived info model: bmi, smoking, alcohol, physical activity, family history 
pheno <- merge(cov[,-c(11:14,17)], atlas308, by="eid_20361")

## read in px info phenotypes
px.info <- read.delim("../../01_inc_dz_phenos/data_input/Patient_model_variables.txt",sep = "\t")
pheno <- merge(pheno, px.info, by="eid_20361")

## rename to olds names to match the rest fo the script
colnames(pheno)[which(colnames(pheno)=="assess_date")] <- c("date_baseline_assessment")

## process ethnicity
pheno$ethnicity <- ifelse(is.na(pheno$ge_eur),2,pheno$ge_eur)

## get protein data and metadata
ol.npx <- fread("../../00_protein_imputation/data_input/Imputed_1536_Expansion_NPX_proteins.txt", data.table = F)
ol.npx <- filter(ol.npx, exclusion.50==0 )

## exclude participants with missing age, sex
ol.npx <- ol.npx[complete.cases(ol.npx[,c("ge_sex", "ge_age")]),]

protein_list <- colnames(ol.npx)[-which(colnames(ol.npx)%in%c("eid_20361","ge_sex","ge_age",
                                                              "randomized_baseline","exclusion.50"))]

## merge phenotypes and proteins
ukbb <- merge(ol.npx, pheno, by="eid_20361")

## select random subset or consortium selected cases only
ukbb <-
  ukbb %>% 
  filter(randomized_baseline==T|!!as.name(dz)==1)

## inverse rank normal transform npx
ukbb[,protein_list] <- sapply(ukbb[,protein_list], RankNorm)

## exclude prevalent cases and incident cases up to 6 months
## restrict to incident cases within the first X years
ukbb <- ukbb %>% 
  mutate(fol = difftime(!!as.name(dz.date),date_baseline_assessment, units = "weeks")/52.25) %>% 
  filter((fol>0.5&fol<inc.yrs)|is.na(!!as.name(dz.date)))

## censor (including censoring by death)
death.data <- read.delim("../../01_inc_dz_phenos/data_input/Death_dates_proteomic_participants.txt",sep = "\t")
death.data$date_of_death <- as.Date(death.data$date_of_death, format = "%Y-%m-%d")

ukbb <- merge(ukbb, death.data[,c("eid_20361","date_of_death")], by="eid_20361", all.x=T)

if(inc.yrs==10){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2020-12-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
} else if(inc.yrs==5){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2016-05-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
}else{stop("limit to incident casese within X years not defined")}

if(inc.yrs==10){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2020-12-31"),
                         difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                        ukbb$fol)
} else if(inc.yrs==5){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2016-05-31"),
                          difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                          ukbb$fol)
}else{stop("limit to incident casese within X years not defined")}

ukbb$fol <- ifelse(is.na(ukbb$fol_d), ukbb$fol, ukbb$fol_d)

## for sex specific outcome exclude males/females
if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_pri_breast",
           "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical",
           "atlas_inc_endometrial_hyper","atlas_inc_endometriosis","atlas_inc_female_genital_prolapse",
           "atlas_inc_leiomyoma","atlas_inc_menorrhagia",
           "atlas_inc_pid","atlas_inc_pmb")){
  # ukbb <- ukbb %>% filter(ge_sex=="female")
  ukbb <- ukbb %>% filter(ge_sex==0)
} else if(dz%in%c("atlas_inc_pri_prost","atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  # ukbb <- ukbb %>% filter(ge_sex=="male")
  ukbb <- ukbb %>% filter(ge_sex==1)
} else{ukbb <- ukbb}

## recode sex to numeric
# ukbb$ge_sex <- ifelse(ukbb$ge_sex=="male",1,
#                       ifelse(ukbb$ge_sex=="female",0,NA))

n.cases <- sum(ukbb[,dz]==1)

## add noise variables
rand.vars <- lapply(1:10, function(x){
  tmp <- runif(nrow(ukbb), min = min(ukbb$NPPB_20049), 
               max = max(ukbb$NPPB_20049))
  return(tmp)
})
rand.vars <- as.data.frame(do.call(cbind, rand.vars))
colnames(rand.vars) <- paste0("rand_var_",1:10)

ukbb <- cbind(ukbb, rand.vars)

####----    ML framework    ----####

## split into train/optimize/test or train/test depending on case number
ukbb[,dz] <- as.factor(ukbb[,dz])

if(n.cases/4>200){
  ## feature selection set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.5, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.opti <- ukbb[-trainIndex,]
  
  ## optimization and testing
  set.seed(3456) ## for reproducibility
  trainIndex2 <- createDataPartition(u.opti[,dz], p=0.5, list = F, times = 1)
  u.test <- u.opti[trainIndex2,]
  u.opti <- u.opti[-trainIndex2,]
  }else {
  ## feature selection set and test set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.7, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.test <- ukbb[-trainIndex,]
  }
  
## draw 200 random subsamples of train set
jj <- lapply(1:200, function(x) sample(nrow(u.train), round(nrow(u.train)*0.50), replace=F))

##-- feature selection  --##
las.morb <- mclapply(1:length(jj), function(i){
  print(i) ## just for test purposes
  tmp <- train(u.train[jj[[i]],c(protein_list,colnames(rand.vars))],
               u.train[jj[[i]],dz],
               metric = "Kappa",
               method="glmnet",
               family="binomial",
               tuneGrid=as.data.frame(expand_grid(
                 alpha=1,lambda=10^-seq(4,.25,-.5))),
               trControl= trainControl(method = "repeatedcv",
                                       number = 3, repeats = 5,
                                       sampling = "rose",
                                       allowParallel = T)
  )
  cf <- as.matrix(coef(tmp$finalModel,
                       s = tmp$finalModel$lambdaOpt))
  gc()
  # return(list("cf" = cf , "cv.tune" = tmp))
  return(cf)
}, mc.cores = 12, mc.allow.recursive = T)

save(las.morb, file=paste0("../data_output/feature_selection/fs_",dz,"_",inc.yrs,"yrs",".RData"))

## generate feature selection ranking 
# load(paste0("../data_output/feature_selection/fs_",dz,"_",inc.yrs,"yrs",".RData"))
p.select <- matrix(,nrow = nrow(las.morb[[1]]), ncol = length(las.morb))

for(i in 1:length(las.morb)){
  tryCatch({
    tmp <- las.morb[[i]][,1]
    if(is.null(tmp)){
      return(tmp <- rep(0,nrow(p.select)))
    }else {
      NULL
    }
  }, error=function(e){
    NULL
  })
  p.select[,i] <- tmp
}
row.names(p.select) <- names(las.morb[[1]][,1])

## count matrix
p.select <-
  p.select %>% 
  as_tibble(rownames = "mrc_olink.id") %>% 
  filter(mrc_olink.id!="(Intercept)") %>% 
  mutate(select = abs(rowSums(across(where(is.numeric))))) %>% 
  arrange(desc(select)) 
  
p.select <- as.data.frame(p.select)

p.select$select.perc <- (p.select$select)/max(p.select$select)

####---- optimization and performance testing ----####

## survival objects
surv.train <- Surv(u.train$fol, u.train[,dz])
## optimization 
if(n.cases/4>200){
  surv.opt <- Surv(u.opti$fol, u.opti[,dz])
}else{surv.opt <- Surv(u.train$fol, u.train[,dz])}
## testing
surv.test <- Surv(u.test$fol, u.test[,dz])

## function to run optimization and testings
source("../../../02_inc_dz_prediction/bin/functions/coxnet_ridge_optim_boot.R")

## px information model basic 
clin.predictors <- c("ge_age","ge_sex","bmi","ethnicity","smoke_never",
                     "smoke_previous","smoke_current","alcohol_daily",
                     "alcohol_4_week","alcohol_2_week","alcohol_3_month","alcohol_occassion",
                     "alcohol_never")

## define additional family history for specific diseases
if(!dz%in%c("atlas_inc_diabetes_t2","atlas_inc_pri_prost", "atlas_inc_copd","atlas_inc_isch_stroke",
            "atlas_inc_af","atlas_inc_av_block_3","atlas_inc_chd_nos","atlas_inc_hf", 
            "atlas_inc_peripheral_arterial_disease","atlas_inc_parkinsons","atlas_inc_pri_lung",
            "atlas_inc_pri_bowel", "atlas_inc_pri_breast","atlas_inc_pri_ovarian","atlas_inc_pri_uterine",
            "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
            "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
            "atlas_inc_hydrocele","atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb",
            "atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu")){
  clin.predictors <- clin.predictors
} else if(dz=="atlas_inc_diabetes_t2") {
  clin.predictors <- c(clin.predictors,"father_diabetes","mother_diabetes")
} else if(dz=="atlas_inc_pri_prost") {
  clin.predictors <- c(clin.predictors[-2],"father_prostate_cancer")
} else if(dz=="atlas_inc_copd") {
  clin.predictors <- c(clin.predictors,"father_COPD","mother_COPD")
} else if(dz%in%c("atlas_inc_isch_stroke","atlas_inc_af",
                  "atlas_inc_chd_nos","atlas_inc_hf","atlas_inc_peripheral_arterial_disease")) {
  clin.predictors <- c(clin.predictors,"father_heart_disease","mother_heart_disease")
} else if(dz=="atlas_inc_parkinsons") {
  clin.predictors <- c(clin.predictors,"father_parkinsons","mother_parkinsons")
} else if(dz%in%c("atlas_inc_pri_lung","atlas_inc_sec_lung")) {
  clin.predictors <- c(clin.predictors,"father_lung_cancer","mother_lung_cancer")
} else if(dz=="atlas_inc_pri_bowel") {
  clin.predictors <- c(clin.predictors,"father_bowel_cancer","mother_bowel_cancer")
} else if(dz=="atlas_inc_pri_breast") {
  clin.predictors <- c(clin.predictors[-2],"mother_breast_cancer")
} else if(dz=="atlas_inc_hypertension") {
  clin.predictors <- c(clin.predictors[-2],"mother_high_bp","father_high_bp")
} else if(dz=="atlas_inc_depression") {
  clin.predictors <- c(clin.predictors[-2],"mother_depression","father_depression")
} else if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_benign_ovary",
                  "atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
                  "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
                  "atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb")) {
  clin.predictors <- c(clin.predictors[-2])
} else if(dz%in%c("atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  clin.predictors <- c(clin.predictors[-2])
}


boot.samples <- lapply(1:1000, function(x) sample(nrow(u.test), round(nrow(u.test)*1), replace = T))

## mode px info model optimization and testing
if(n.cases/4>200){
  res.clin <- coxnet.optim.r(Train.data = u.opti,
                           pred.vec = clin.predictors,
                           Train.surv.data = surv.opt,
                           times = 1000, 
                           boot.samples = boot.samples,
                           Test.data = u.test,
                           Test.surv.data = surv.test)
} else{
  res.clin <- coxnet.optim.r(Train.data = u.train,
                           pred.vec = clin.predictors,
                           Train.surv.data = surv.opt,
                           times = 1000, 
                           boot.samples = boot.samples,
                           Test.data = u.test,
                           Test.surv.data = surv.test)
  
}
save(res.clin, file = paste0("../data_output/px_info_models/px_info_res_",dz,"_",inc.yrs,"yrs",".RData"))

## selection threshold based on noise variables
sel.t <- p.select$select.perc[grep("rand", p.select$mrc_olink.id)[1]]

prots.opti <- p.select %>% 
  filter(select.perc>sel.t) %>% 
  top_n(20, select.perc)
prots.opti <- prots.opti$mrc_olink.id
if(length(prots.opti)==0){
  prots.opti <- p.select[p.select$select.per==sel.t]
  prots.opti <- prots.opti$mrc_olink.id
}

# run px info model + proteins selected above noise
if(n.cases/4>200){
    res.clin.prots.all <- coxnet.optim.r(Train.data = u.opti,
                                   pred.vec = c(clin.predictors,prots.opti),
                                   Train.surv.data = surv.opt,
                                   times = 1000, 
                                   boot.samples = boot.samples,
                                   Test.data = u.test,
                                   Test.surv.data = surv.test)
    res.prots.all <- coxnet.optim.r(Train.data = u.opti,
                                    pred.vec = prots.opti,
                                    Train.surv.data = surv.opt,
                                    times = 1000, 
                                    boot.samples = boot.samples,
                                    Test.data = u.test,
                                    Test.surv.data = surv.test)
    } else{
      res.clin.prots.all <- coxnet.optim.r(Train.data = u.train,
                                   pred.vec = c(clin.predictors,prots.opti),
                                   Train.surv.data = surv.opt,
                                   times = 1000, 
                                   boot.samples = boot.samples,
                                   Test.data = u.test,
                                   Test.surv.data = surv.test)
      res.prots.all <- coxnet.optim.r(Train.data = u.train,
                                      pred.vec = prots.opti,
                                      Train.surv.data = surv.opt,
                                      times = 1000, 
                                      boot.samples = boot.samples,
                                      Test.data = u.test,
                                      Test.surv.data = surv.test)
    }
save(res.clin.prots.all, file = paste0("../data_output/px_info_protein_models/px_info_all_prots_above_noise_res_",dz,"_",inc.yrs,"yrs",".RData"))
save(res.prots.all, file = paste0("../data_output/protein_models/prot_all_res_",dz,"_",inc.yrs,"yrs",".RData"))

# run px info model + proteins selected above noise in train set
if(n.cases/4>200){
  res.clin.prots.all <- coxnet.optim.r(Train.data = u.train,
                                     pred.vec = c(clin.predictors,prots.opti),
                                     Train.surv.data = surv.train,
                                     times = 1000, 
                                     boot.samples = boot.samples,
                                     Test.data = u.test,
                                     Test.surv.data = surv.test)
} else {res.clin.prots.all <- res.clin.prots.all}

## proteins top 5 and 10 only based on weights*selection_percentage
opti.prots <- res.clin.prots.all$opt.coefficients[grep("_[0-9]{5}",row.names(res.clin.prots.all$opt.coefficients)),]
opti.prots <- opti.prots %>%
  as_tibble(rownames = "prot") %>%
  mutate(abs.coef = abs(value)) %>% 
  filter(abs.coef!=0)

## merge with selection percentage
opti.prots <- merge(opti.prots, p.select[,c("mrc_olink.id","select","select.perc")],
                    by.x="prot", by.y="mrc_olink.id")

if(nrow(opti.prots)>=10){
  opti.prots10 <- opti.prots %>%
    mutate(sel.opti.score=abs.coef*select) %>%
    top_n(10, sel.opti.score)
  prots.opti10 <- opti.prots10$prot
  
  ## run px info model + 10 top proteins based on weightsXselection.perc
  if(n.cases/4>200){
    res.clin.prots10 <- coxnet.optim.r(Train.data = u.opti,
                                     pred.vec = c(clin.predictors,prots.opti10),
                                     Train.surv.data = surv.opt,
                                     times = 1000, 
                                     boot.samples = boot.samples,
                                     Test.data = u.test,
                                     Test.surv.data = surv.test)
    res.prots10 <- coxnet.optim.r(Train.data = u.opti,
                                  pred.vec = prots.opti10,
                                  Train.surv.data = surv.opt,
                                  times = 1000, 
                                  boot.samples = boot.samples,
                                  Test.data = u.test,
                                  Test.surv.data = surv.test)
  } else{
    res.clin.prots10 <- coxnet.optim.r(Train.data = u.train,
                                     pred.vec = c(clin.predictors,prots.opti10),
                                     Train.surv.data = surv.opt,
                                     times = 1000, 
                                     boot.samples = boot.samples,
                                     Test.data = u.test,
                                     Test.surv.data = surv.test)
    res.prots10 <- coxnet.optim.r(Train.data = u.train,
                                  pred.vec = prots.opti10,
                                  Train.surv.data = surv.opt,
                                  times = 1000, 
                                  boot.samples = boot.samples,
                                  Test.data = u.test,
                                  Test.surv.data = surv.test)
    
  }
  save(res.clin.prots10, file = paste0("../data_output/px_info_protein_models/px_info_top10_prots_res_",dz,"_",inc.yrs,"yrs",".RData"))
  save(res.prots10, file = paste0("../data_output/protein_models/prot_top10_res_",dz,"_",inc.yrs,"yrs",".RData"))
  
}

## run px info model + 5 top proteins based on weightsXselection.perc
if(nrow(opti.prots)>=5){
  opti.prots <- opti.prots %>%
  mutate(sel.opti.score=abs.coef*select.perc) %>%
  top_n(5, sel.opti.score)
  
  prots.opti <- opti.prots$prot
} else{prots.opti <- opti.prots$prot}


if(n.cases/4>200){
  res.clin.prots <- coxnet.optim.r(Train.data = u.opti,
                                 pred.vec = c(clin.predictors,prots.opti),
                                 Train.surv.data = surv.opt,
                                 times = 1000, 
                                 boot.samples = boot.samples,
                                 Test.data = u.test,
                                 Test.surv.data = surv.test)
} else{
  res.clin.prots <- coxnet.optim.r(Train.data = u.train,
                                 pred.vec = c(clin.predictors,prots.opti),
                                 Train.surv.data = surv.opt,
                                 times = 1000, 
                                 boot.samples = boot.samples,
                                 Test.data = u.test,
                                 Test.surv.data = surv.test)
  
}
save(res.clin.prots, file = paste0("../data_output/px_info_protein_models/px_info_top5_prots_res_",dz,"_",inc.yrs,"yrs",".RData"))

## optimization on top 5 proteins
if(n.cases/4>200){
  res.prot <- coxnet.optim.r(Train.data = u.opti,
                           pred.vec = prots.opti,
                           Train.surv.data = surv.opt,
                           times = 1000, 
                           boot.samples = boot.samples, 
                           Test.data = u.test,
                           Test.surv.data = surv.test)
} else{
  res.prot <- coxnet.optim.r(Train.data = u.train,
                           pred.vec = prots.opti,
                           Train.surv.data = surv.opt,
                           times = 1000, 
                           boot.samples = boot.samples,
                           Test.data = u.test,
                           Test.surv.data = surv.test)
  
}
## add feature selection to the results
res.prot$feature.selection <- p.select

## add case numbers for training and testing sets
if((n.cases/4)>200){
  case.num <- data.frame(train.cases = sum(u.train[,dz]==1),
                         opti.cases = sum(u.opti[,dz]==1),
                         test.cases = sum(u.test[,dz]==1),
                         train.controls = sum(u.train[,dz]!=1),
                         opti.controls = sum(u.opti[,dz]!=1),
                         test.controls = sum(u.test[,dz]!=1))
} else{
  case.num <- data.frame(train.cases = sum(u.train[,dz]==1),
                         opti.cases = NA,
                         test.cases = sum(u.test[,dz]==1),
                         train.controls = sum(u.train[,dz]!=1),
                         opti.controls = NA,
                         test.controls = sum(u.test[,dz]!=1))
}
res.prot$case.num <- case.num
## save results 
save(res.prot, file = paste0("../data_output/protein_models/prot_top5_res_",dz,"_",inc.yrs,"yrs",".RData"))

###################################################
#### Script Name: ML_pipeline_protein_pxinfo_prediction_pred
#### Author: Julia Carrasco Zanini Sanchez
#### Description: script to ML framework for incident diseases in atlas 308 
#### Date : 23/03/2023

rm(list = ls())

## load packages 
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(RNOmni)
# library(glue)
library(DBI)
library(dbplyr)
library(magrittr)
# library(broom)
library(survival)
library(data.table)
library(glue)
library(parallel)
library(caret)
library(glmnet)
library(ROSE)
# library(odbc)

## define disease and restrict to incident cases within X years
args <- commandArgs(trailingOnly = T)

dz <- as.character(args[1])
dz.date <- as.character(args[2])
inc.yrs <- as.numeric(args[3])

paste0("Prediction:", dz, "; ", dz.date, " at ", inc.yrs, "years")

## set up working directory 
setwd("/home/ukbb_inc_dz_prediction/Explore_1536_Expansion_analyses/02_inc_dz_prediction/bin/") 

## read in data 
atlas308 <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_atlas308_annotations.txt", data.table = F)
cov <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_covariate_annotations.txt", data.table = F)

####----  Data Cleaning and wrangling   ----####
atlas308 <- as.data.frame(atlas308 %>% 
  filter(eid_20361%in%cov$eid_20361) %>% 
  pivot_wider(names_from = variable, values_from =value))

## phenotypes missing for patient-derived info model: bmi, smoking, alcohol, physical activity, family history 
pheno <- merge(cov[,-c(11:14,17)], atlas308, by="eid_20361")

## read in px info phenotypes
px.info <- read.delim("../../01_inc_dz_phenos/data_input/Patient_model_variables.txt",sep = "\t")
pheno <- merge(pheno, px.info, by="eid_20361")

## rename to olds names to match the rest fo the script
colnames(pheno)[which(colnames(pheno)=="assess_date")] <- c("date_baseline_assessment")

## process ethnicity
pheno$ethnicity <- ifelse(is.na(pheno$ge_eur),2,pheno$ge_eur)

## get protein data and metadata
ol.npx <- fread("../../00_protein_imputation/data_input/Imputed_1536_Expansion_NPX_proteins.txt", data.table = F)
ol.npx <- filter(ol.npx, exclusion.50==0 )

## exclude participants with missing age, sex
ol.npx <- ol.npx[complete.cases(ol.npx[,c("ge_sex", "ge_age")]),]

protein_list <- colnames(ol.npx)[-which(colnames(ol.npx)%in%c("eid_20361","ge_sex","ge_age",
                                                              "randomized_baseline","exclusion.50"))]

## merge phenotypes and proteins
ukbb <- merge(ol.npx, pheno, by="eid_20361")
## randomised baseline
ukbb <- ukbb %>% 
  filter(randomized_baseline==T)
## inverse rank normal transform npx
ukbb[,protein_list] <- sapply(ukbb[,protein_list], RankNorm)

## exclude prevalent cases and incident cases up to 6 months
## restrict to incident cases within the first X years
ukbb <- ukbb %>% 
  mutate(fol = difftime(!!as.name(dz.date),date_baseline_assessment, units = "weeks")/52.25) %>% 
  filter((fol>0.5&fol<inc.yrs)|is.na(!!as.name(dz.date)))

## censor (including censoring by death)
death.data <- read.delim("../../01_inc_dz_phenos/data_input/Death_dates_proteomic_participants.txt",sep = "\t")
death.data$date_of_death <- as.Date(death.data$date_of_death, format = "%Y-%m-%d")

ukbb <- merge(ukbb, death.data[,c("eid_20361","date_of_death")], by="eid_20361", all.x=T)

if(inc.yrs==10){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2020-12-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
} else if(inc.yrs==5){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2016-05-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
}else{stop("limit to incident casese within X years not defined")}

if(inc.yrs==10){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2020-12-31"),
                         difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                        ukbb$fol)
} else if(inc.yrs==5){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2016-05-31"),
                          difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                          ukbb$fol)
}else{stop("limit to incident casese within X years not defined")}

ukbb$fol <- ifelse(is.na(ukbb$fol_d), ukbb$fol, ukbb$fol_d)

## for sex specific outcome exclude males/females
if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_pri_breast",
           "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical",
           "atlas_inc_endometrial_hyper","atlas_inc_endometriosis","atlas_inc_female_genital_prolapse",
           "atlas_inc_leiomyoma","atlas_inc_menorrhagia",
           "atlas_inc_pid","atlas_inc_pmb")){
  # ukbb <- ukbb %>% filter(ge_sex=="female")
  ukbb <- ukbb %>% filter(ge_sex==0)
} else if(dz%in%c("atlas_inc_pri_prost","atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  # ukbb <- ukbb %>% filter(ge_sex=="male")
  ukbb <- ukbb %>% filter(ge_sex==1)
} else{ukbb <- ukbb}

## recode sex to numeric
# ukbb$ge_sex <- ifelse(ukbb$ge_sex=="male",1,
#                       ifelse(ukbb$ge_sex=="female",0,NA))

n.cases <- sum(ukbb[,dz]==1)

## add noise variables
rand.vars <- lapply(1:10, function(x){
  tmp <- runif(nrow(ukbb), min = min(ukbb$NPPB_20049), 
               max = max(ukbb$NPPB_20049))
  return(tmp)
})
rand.vars <- as.data.frame(do.call(cbind, rand.vars))
colnames(rand.vars) <- paste0("rand_var_",1:10)

ukbb <- cbind(ukbb, rand.vars)

####----    ML framework    ----####

## split into train/optimize/test or train/test depending on case number
ukbb[,dz] <- as.factor(ukbb[,dz])

if(n.cases/4>200){
  ## feature selection set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.5, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.opti <- ukbb[-trainIndex,]
  
  ## optimization and testing
  set.seed(3456) ## for reproducibility
  trainIndex2 <- createDataPartition(u.opti[,dz], p=0.5, list = F, times = 1)
  u.test <- u.opti[trainIndex2,]
  u.opti <- u.opti[-trainIndex2,]
  }else {
  ## feature selection set and test set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.7, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.test <- ukbb[-trainIndex,]
  }
  
## draw 200 random subsamples of train set
jj <- lapply(1:200, function(x) sample(nrow(u.train), round(nrow(u.train)*0.50), replace=F))

##-- feature selection  --##
las.morb <- mclapply(1:length(jj), function(i){
  print(i) ## just for test purposes
  tmp <- train(u.train[jj[[i]],c(protein_list,colnames(rand.vars))],
               u.train[jj[[i]],dz],
               metric = "Kappa",
               method="glmnet",
               family="binomial",
               tuneGrid=as.data.frame(expand_grid(
                 alpha=1,lambda=10^-seq(4,.25,-.5))),
               trControl= trainControl(method = "repeatedcv",
                                       number = 3, repeats = 5,
                                       sampling = "rose",
                                       allowParallel = T)
  )
  cf <- as.matrix(coef(tmp$finalModel,
                       s = tmp$finalModel$lambdaOpt))
  gc()
  # return(list("cf" = cf , "cv.tune" = tmp))
  return(cf)
}, mc.cores = 12, mc.allow.recursive = T)

save(las.morb, file=paste0("../data_output/feature_selection/fs_",dz,"_",inc.yrs,"yrs",".RData"))

## generate feature selection ranking 
# load(paste0("../data_output/feature_selection/fs_",dz,"_",inc.yrs,"yrs",".RData"))
p.select <- matrix(,nrow = nrow(las.morb[[1]]), ncol = length(las.morb))

for(i in 1:length(las.morb)){
  tryCatch({
    tmp <- las.morb[[i]][,1]
    if(is.null(tmp)){
      return(tmp <- rep(0,nrow(p.select)))
    }else {
      NULL
    }
  }, error=function(e){
    NULL
  })
  p.select[,i] <- tmp
}
row.names(p.select) <- names(las.morb[[1]][,1])

## count matrix
p.select <-
  p.select %>% 
  as_tibble(rownames = "mrc_olink.id") %>% 
  filter(mrc_olink.id!="(Intercept)") %>% 
  mutate(select = abs(rowSums(across(where(is.numeric))))) %>% 
  arrange(desc(select)) 
  
p.select <- as.data.frame(p.select)

p.select$select.perc <- (p.select$select)/max(p.select$select)

####---- optimization and performance testing ----####

## survival objects
surv.train <- Surv(u.train$fol, u.train[,dz])
## optimization 
if(n.cases/4>200){
  surv.opt <- Surv(u.opti$fol, u.opti[,dz])
}else{surv.opt <- Surv(u.train$fol, u.train[,dz])}
## testing
surv.test <- Surv(u.test$fol, u.test[,dz])

## function to run optimization and testings
source("../../../02_inc_dz_prediction/bin/functions/coxnet_ridge_optim_boot.R")

## px information model basic 
clin.predictors <- c("ge_age","ge_sex","bmi","ethnicity","smoke_never",
                     "smoke_previous","smoke_current","alcohol_daily",
                     "alcohol_4_week","alcohol_2_week","alcohol_3_month","alcohol_occassion",
                     "alcohol_never")

## define additional family history for specific diseases
if(!dz%in%c("atlas_inc_diabetes_t2","atlas_inc_pri_prost", "atlas_inc_copd","atlas_inc_isch_stroke",
            "atlas_inc_af","atlas_inc_av_block_3","atlas_inc_chd_nos","atlas_inc_hf", 
            "atlas_inc_peripheral_arterial_disease","atlas_inc_parkinsons","atlas_inc_pri_lung",
            "atlas_inc_pri_bowel", "atlas_inc_pri_breast","atlas_inc_pri_ovarian","atlas_inc_pri_uterine",
            "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
            "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
            "atlas_inc_hydrocele","atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb",
            "atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu")){
  clin.predictors <- clin.predictors
} else if(dz=="atlas_inc_diabetes_t2") {
  clin.predictors <- c(clin.predictors,"father_diabetes","mother_diabetes")
} else if(dz=="atlas_inc_pri_prost") {
  clin.predictors <- c(clin.predictors[-2],"father_prostate_cancer")
} else if(dz=="atlas_inc_copd") {
  clin.predictors <- c(clin.predictors,"father_COPD","mother_COPD")
} else if(dz%in%c("atlas_inc_isch_stroke","atlas_inc_af",
                  "atlas_inc_chd_nos","atlas_inc_hf","atlas_inc_peripheral_arterial_disease")) {
  clin.predictors <- c(clin.predictors,"father_heart_disease","mother_heart_disease")
} else if(dz=="atlas_inc_parkinsons") {
  clin.predictors <- c(clin.predictors,"father_parkinsons","mother_parkinsons")
} else if(dz%in%c("atlas_inc_pri_lung","atlas_inc_sec_lung")) {
  clin.predictors <- c(clin.predictors,"father_lung_cancer","mother_lung_cancer")
} else if(dz=="atlas_inc_pri_bowel") {
  clin.predictors <- c(clin.predictors,"father_bowel_cancer","mother_bowel_cancer")
} else if(dz=="atlas_inc_pri_breast") {
  clin.predictors <- c(clin.predictors[-2],"mother_breast_cancer")
} else if(dz=="atlas_inc_hypertension") {
  clin.predictors <- c(clin.predictors[-2],"mother_high_bp","father_high_bp")
} else if(dz=="atlas_inc_depression") {
  clin.predictors <- c(clin.predictors[-2],"mother_depression","father_depression")
} else if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_benign_ovary",
                  "atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
                  "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
                  "atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb")) {
  clin.predictors <- c(clin.predictors[-2])
} else if(dz%in%c("atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  clin.predictors <- c(clin.predictors[-2])
}


boot.samples <- lapply(1:1000, function(x) sample(nrow(u.test), round(nrow(u.test)*1), replace = T))

## mode px info model optimization and testing
if(n.cases/4>200){
  res.clin <- coxnet.optim.r(Train.data = u.opti,
                           pred.vec = clin.predictors,
                           Train.surv.data = surv.opt,
                           times = 1000, 
                           boot.samples = boot.samples,
                           Test.data = u.test,
                           Test.surv.data = surv.test)
} else{
  res.clin <- coxnet.optim.r(Train.data = u.train,
                           pred.vec = clin.predictors,
                           Train.surv.data = surv.opt,
                           times = 1000, 
                           boot.samples = boot.samples,
                           Test.data = u.test,
                           Test.surv.data = surv.test)
  
}
save(res.clin, file = paste0("../data_output/px_info_models/px_info_res_",dz,"_",inc.yrs,"yrs",".RData"))

## selection threshold based on noise variables
sel.t <- p.select$select.perc[grep("rand", p.select$mrc_olink.id)[1]]

prots.opti <- p.select %>% 
  filter(select.perc>sel.t) %>% 
  top_n(20, select.perc)
prots.opti <- prots.opti$mrc_olink.id
if(length(prots.opti)==0){
  prots.opti <- p.select[p.select$select.per==sel.t]
  prots.opti <- prots.opti$mrc_olink.id
}

# run px info model + proteins selected above noise
if(n.cases/4>200){
    res.clin.prots.all <- coxnet.optim.r(Train.data = u.opti,
                                   pred.vec = c(clin.predictors,prots.opti),
                                   Train.surv.data = surv.opt,
                                   times = 1000, 
                                   boot.samples = boot.samples,
                                   Test.data = u.test,
                                   Test.surv.data = surv.test)
    res.prots.all <- coxnet.optim.r(Train.data = u.opti,
                                    pred.vec = prots.opti,
                                    Train.surv.data = surv.opt,
                                    times = 1000, 
                                    boot.samples = boot.samples,
                                    Test.data = u.test,
                                    Test.surv.data = surv.test)
    } else{
      res.clin.prots.all <- coxnet.optim.r(Train.data = u.train,
                                   pred.vec = c(clin.predictors,prots.opti),
                                   Train.surv.data = surv.opt,
                                   times = 1000, 
                                   boot.samples = boot.samples,
                                   Test.data = u.test,
                                   Test.surv.data = surv.test)
      res.prots.all <- coxnet.optim.r(Train.data = u.train,
                                      pred.vec = prots.opti,
                                      Train.surv.data = surv.opt,
                                      times = 1000, 
                                      boot.samples = boot.samples,
                                      Test.data = u.test,
                                      Test.surv.data = surv.test)
    }
save(res.clin.prots.all, file = paste0("../data_output/px_info_protein_models/px_info_all_prots_above_noise_res_",dz,"_",inc.yrs,"yrs",".RData"))
save(res.prots.all, file = paste0("../data_output/protein_models/prot_all_res_",dz,"_",inc.yrs,"yrs",".RData"))

# run px info model + proteins selected above noise in train set
if(n.cases/4>200){
  res.clin.prots.all <- coxnet.optim.r(Train.data = u.train,
                                     pred.vec = c(clin.predictors,prots.opti),
                                     Train.surv.data = surv.train,
                                     times = 1000, 
                                     boot.samples = boot.samples,
                                     Test.data = u.test,
                                     Test.surv.data = surv.test)
} else {res.clin.prots.all <- res.clin.prots.all}

## proteins top 5 and 10 only based on weights*selection_percentage
opti.prots <- res.clin.prots.all$opt.coefficients[grep("_[0-9]{5}",row.names(res.clin.prots.all$opt.coefficients)),]
opti.prots <- opti.prots %>%
  as_tibble(rownames = "prot") %>%
  mutate(abs.coef = abs(value)) %>% 
  filter(abs.coef!=0)

## merge with selection percentage
opti.prots <- merge(opti.prots, p.select[,c("mrc_olink.id","select","select.perc")],
                    by.x="prot", by.y="mrc_olink.id")

if(nrow(opti.prots)>=10){
  opti.prots10 <- opti.prots %>%
    mutate(sel.opti.score=abs.coef*select) %>%
    top_n(10, sel.opti.score)
  prots.opti10 <- opti.prots10$prot
  
  ## run px info model + 10 top proteins based on weightsXselection.perc
  if(n.cases/4>200){
    res.clin.prots10 <- coxnet.optim.r(Train.data = u.opti,
                                     pred.vec = c(clin.predictors,prots.opti10),
                                     Train.surv.data = surv.opt,
                                     times = 1000, 
                                     boot.samples = boot.samples,
                                     Test.data = u.test,
                                     Test.surv.data = surv.test)
    res.prots10 <- coxnet.optim.r(Train.data = u.opti,
                                  pred.vec = prots.opti10,
                                  Train.surv.data = surv.opt,
                                  times = 1000, 
                                  boot.samples = boot.samples,
                                  Test.data = u.test,
                                  Test.surv.data = surv.test)
  } else{
    res.clin.prots10 <- coxnet.optim.r(Train.data = u.train,
                                     pred.vec = c(clin.predictors,prots.opti10),
                                     Train.surv.data = surv.opt,
                                     times = 1000, 
                                     boot.samples = boot.samples,
                                     Test.data = u.test,
                                     Test.surv.data = surv.test)
    res.prots10 <- coxnet.optim.r(Train.data = u.train,
                                  pred.vec = prots.opti10,
                                  Train.surv.data = surv.opt,
                                  times = 1000, 
                                  boot.samples = boot.samples,
                                  Test.data = u.test,
                                  Test.surv.data = surv.test)
    
  }
  save(res.clin.prots10, file = paste0("../data_output/px_info_protein_models/px_info_top10_prots_res_",dz,"_",inc.yrs,"yrs",".RData"))
  save(res.prots10, file = paste0("../data_output/protein_models/prot_top10_res_",dz,"_",inc.yrs,"yrs",".RData"))
  
}

## run px info model + 5 top proteins based on weightsXselection.perc
if(nrow(opti.prots)>=5){
  opti.prots <- opti.prots %>%
  mutate(sel.opti.score=abs.coef*select.perc) %>%
  top_n(5, sel.opti.score)
  
  prots.opti <- opti.prots$prot
} else{prots.opti <- opti.prots$prot}


if(n.cases/4>200){
  res.clin.prots <- coxnet.optim.r(Train.data = u.opti,
                                 pred.vec = c(clin.predictors,prots.opti),
                                 Train.surv.data = surv.opt,
                                 times = 1000, 
                                 boot.samples = boot.samples,
                                 Test.data = u.test,
                                 Test.surv.data = surv.test)
} else{
  res.clin.prots <- coxnet.optim.r(Train.data = u.train,
                                 pred.vec = c(clin.predictors,prots.opti),
                                 Train.surv.data = surv.opt,
                                 times = 1000, 
                                 boot.samples = boot.samples,
                                 Test.data = u.test,
                                 Test.surv.data = surv.test)
  
}
save(res.clin.prots, file = paste0("../data_output/px_info_protein_models/px_info_top5_prots_res_",dz,"_",inc.yrs,"yrs",".RData"))

## optimization on top 5 proteins
if(n.cases/4>200){
  res.prot <- coxnet.optim.r(Train.data = u.opti,
                           pred.vec = prots.opti,
                           Train.surv.data = surv.opt,
                           times = 1000, 
                           boot.samples = boot.samples, 
                           Test.data = u.test,
                           Test.surv.data = surv.test)
} else{
  res.prot <- coxnet.optim.r(Train.data = u.train,
                           pred.vec = prots.opti,
                           Train.surv.data = surv.opt,
                           times = 1000, 
                           boot.samples = boot.samples,
                           Test.data = u.test,
                           Test.surv.data = surv.test)
  
}
## add feature selection to the results
res.prot$feature.selection <- p.select

## add case numbers for training and testing sets
if((n.cases/4)>200){
  case.num <- data.frame(train.cases = sum(u.train[,dz]==1),
                         opti.cases = sum(u.opti[,dz]==1),
                         test.cases = sum(u.test[,dz]==1),
                         train.controls = sum(u.train[,dz]!=1),
                         opti.controls = sum(u.opti[,dz]!=1),
                         test.controls = sum(u.test[,dz]!=1))
} else{
  case.num <- data.frame(train.cases = sum(u.train[,dz]==1),
                         opti.cases = NA,
                         test.cases = sum(u.test[,dz]==1),
                         train.controls = sum(u.train[,dz]!=1),
                         opti.controls = NA,
                         test.controls = sum(u.test[,dz]!=1))
}
res.prot$case.num <- case.num
## save results 
save(res.prot, file = paste0("../data_output/protein_models/prot_top5_res_",dz,"_",inc.yrs,"yrs",".RData"))

###################################################
#### Script Name: ML_pipeline_6m_protein_pxinfo_prediction_pred_consortcases
#### Author: Julia Carrasco Zanini Sanchez
#### Description: script to ML framework for incident diseases in atlas 308 
#### Date : 23/03/2023

rm(list = ls())

## load packages 
library(dplyr)
library(tidyr)
library(RNOmni)
# library(glue)
library(DBI)
library(dbplyr)
library(magrittr)
# library(broom)
library(survival)
library(data.table)
library(glue)
library(parallel)
library(caret)
library(glmnet)
library(ROSE)

## define disease and restrict to incident cases within X years
args <- commandArgs(trailingOnly = T)

dz <- as.character(args[1])
dz.date <- as.character(args[2])
inc.yrs <- as.numeric(args[3])

paste0("Prediction:", dz, "; ", dz.date, " at ", inc.yrs, "years")

## set up working directory 
setwd("/home/ukbb_inc_dz_prediction/Explore_1536_Expansion_analyses/02_inc_dz_prediction/bin/") 

## read in data 
atlas308 <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_atlas308_annotations.txt", data.table = F)
cov <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_covariate_annotations.txt", data.table = F)

####----  Data Cleaning and wrangling   ----####
atlas308 <- as.data.frame(atlas308 %>% 
                            filter(eid_20361%in%cov$eid_20361) %>% 
                            pivot_wider(names_from = variable, values_from =value))

## phenotypes missing for patient-derived info model: bmi, smoking, alcohol, physical activity, family history 
pheno <- merge(cov[,-c(11:14,17)], atlas308, by="eid_20361")

## read in px info phenotypes
px.info <- read.delim("../../01_inc_dz_phenos/data_input/Patient_model_variables.txt",sep = "\t")
pheno <- merge(pheno, px.info, by="eid_20361")

## rename to olds names to match the rest fo the script
colnames(pheno)[which(colnames(pheno)=="assess_date")] <- c("date_baseline_assessment")

## process ethnicity
pheno$ethnicity <- ifelse(is.na(pheno$ge_eur),2,pheno$ge_eur)

## get protein data and metadata
ol.npx <- fread("../../00_protein_imputation/data_input/Imputed_1536_Expansion_NPX_proteins.txt", data.table = F)
ol.npx <- filter(ol.npx, exclusion.50==0 )

## exclude participants with missing age, sex
ol.npx <- ol.npx[complete.cases(ol.npx[,c("ge_sex", "ge_age")]),]

protein_list <- colnames(ol.npx)[-which(colnames(ol.npx)%in%c("eid_20361","ge_sex","ge_age",
                                                              "randomized_baseline","exclusion.50"))]

## merge phenotypes and proteins
ukbb <- merge(ol.npx, pheno, by="eid_20361")

## select random subset or consortium selected cases only
ukbb <-
  ukbb %>% 
  filter(randomized_baseline==T|!!as.name(dz)==1)

## inverse rank normal transform npx
ukbb[,protein_list] <- sapply(ukbb[,protein_list], RankNorm)

## exclude prevalent cases and incident cases up to 6 months
## restrict to incident cases within the first X years
ukbb <- ukbb %>% 
  mutate(fol = difftime(!!as.name(dz.date),date_baseline_assessment, units = "weeks")/52.25) %>% 
  filter((fol>0.5&fol<inc.yrs)|is.na(!!as.name(dz.date)))

## censor (including censoring by death)
death.data <- read.delim("../../01_inc_dz_phenos/data_input/Death_dates_proteomic_participants.txt",sep = "\t")
death.data$date_of_death <- as.Date(death.data$date_of_death, format = "%Y-%m-%d")

ukbb <- merge(ukbb, death.data[,c("eid_20361","date_of_death")], by="eid_20361", all.x=T)

if(inc.yrs==10){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2020-12-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
} else if(inc.yrs==5){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2016-05-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
}else{stop("limit to incident casese within X years not defined")}

if(inc.yrs==10){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2020-12-31"),
                       difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                       ukbb$fol)
} else if(inc.yrs==5){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2016-05-31"),
                       difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                       ukbb$fol)
}else{stop("limit to incident casese within X years not defined")}

ukbb$fol <- ifelse(is.na(ukbb$fol_d), ukbb$fol, ukbb$fol_d)

## for sex specific outcome exclude males/females
if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_pri_breast",
           "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical",
           "atlas_inc_endometrial_hyper","atlas_inc_endometriosis","atlas_inc_female_genital_prolapse",
           "atlas_inc_leiomyoma","atlas_inc_menorrhagia",
           "atlas_inc_pid","atlas_inc_pmb")){
  # ukbb <- ukbb %>% filter(ge_sex=="female")
  ukbb <- ukbb %>% filter(ge_sex==0)
} else if(dz%in%c("atlas_inc_pri_prost","atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  # ukbb <- ukbb %>% filter(ge_sex=="male")
  ukbb <- ukbb %>% filter(ge_sex==1)
} else{ukbb <- ukbb}

## recode sex to numeric
# ukbb$ge_sex <- ifelse(ukbb$ge_sex=="male",1,
#                       ifelse(ukbb$ge_sex=="female",0,NA))

n.cases <- sum(ukbb[,dz]==1)

## add noise variables
rand.vars <- lapply(1:10, function(x){
  tmp <- runif(nrow(ukbb), min = min(ukbb$NPPB_20049), 
               max = max(ukbb$NPPB_20049))
  return(tmp)
})
rand.vars <- as.data.frame(do.call(cbind, rand.vars))
colnames(rand.vars) <- paste0("rand_var_",1:10)

ukbb <- cbind(ukbb, rand.vars)

####----    ML framework    ----####

## split into train/optimize/test or train/test depending on case number
ukbb[,dz] <- as.factor(ukbb[,dz])

if(n.cases/4>200){
  ## feature selection set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.5, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.opti <- ukbb[-trainIndex,]
  
  ## optimization and testing
  set.seed(3456) ## for reproducibility
  trainIndex2 <- createDataPartition(u.opti[,dz], p=0.5, list = F, times = 1)
  u.test <- u.opti[trainIndex2,]
  u.opti <- u.opti[-trainIndex2,]
}else {
  ## feature selection set and test set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.7, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.test <- ukbb[-trainIndex,]
}
  
####---- optimization and performance testing ----####

## function to run optimization and testings
source("../../../02_inc_dz_prediction/bin/functions/coxnet_ridge_optim_boot.R")

## px information model basic 
clin.predictors <- c("ge_age","ge_sex","bmi","ethnicity","smoke_never",
                     "smoke_previous","smoke_current","alcohol_daily",
                     "alcohol_4_week","alcohol_2_week","alcohol_3_month","alcohol_occassion",
                     "alcohol_never")

## define additional family history for specific diseases
if(!dz%in%c("atlas_inc_diabetes_t2","atlas_inc_pri_prost", "atlas_inc_copd","atlas_inc_isch_stroke",
            "atlas_inc_af","atlas_inc_av_block_3","atlas_inc_chd_nos","atlas_inc_hf", 
            "atlas_inc_peripheral_arterial_disease","atlas_inc_parkinsons","atlas_inc_pri_lung",
            "atlas_inc_pri_bowel", "atlas_inc_pri_breast","atlas_inc_pri_ovarian","atlas_inc_pri_uterine",
            "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
            "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
            "atlas_inc_hydrocele","atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb",
            "atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu")){
  clin.predictors <- clin.predictors
} else if(dz=="atlas_inc_diabetes_t2") {
  clin.predictors <- c(clin.predictors,"father_diabetes","mother_diabetes")
} else if(dz=="atlas_inc_pri_prost") {
  clin.predictors <- c(clin.predictors[-2],"father_prostate_cancer")
} else if(dz=="atlas_inc_copd") {
  clin.predictors <- c(clin.predictors,"father_COPD","mother_COPD")
} else if(dz%in%c("atlas_inc_isch_stroke","atlas_inc_af",
                  "atlas_inc_chd_nos","atlas_inc_hf","atlas_inc_peripheral_arterial_disease")) {
  clin.predictors <- c(clin.predictors,"father_heart_disease","mother_heart_disease")
} else if(dz=="atlas_inc_parkinsons") {
  clin.predictors <- c(clin.predictors,"father_parkinsons","mother_parkinsons")
} else if(dz%in%c("atlas_inc_pri_lung","atlas_inc_sec_lung")) {
  clin.predictors <- c(clin.predictors,"father_lung_cancer","mother_lung_cancer")
} else if(dz=="atlas_inc_pri_bowel") {
  clin.predictors <- c(clin.predictors,"father_bowel_cancer","mother_bowel_cancer")
} else if(dz=="atlas_inc_pri_breast") {
  clin.predictors <- c(clin.predictors[-2],"mother_breast_cancer")
} else if(dz=="atlas_inc_hypertension") {
  clin.predictors <- c(clin.predictors[-2],"mother_high_bp","father_high_bp")
} else if(dz=="atlas_inc_depression") {
  clin.predictors <- c(clin.predictors[-2],"mother_depression","father_depression")
} else if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_benign_ovary",
                  "atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
                  "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
                  "atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb")) {
  clin.predictors <- c(clin.predictors[-2])
} else if(dz%in%c("atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  clin.predictors <- c(clin.predictors[-2])
}

## read in PRS data 
prs <- read.delim("../../01_inc_dz_phenos/data_input/PRS_standard_ukb.txt", sep = "\t")

p.list <- list("atlas_inc_sle"="prs.sle", "atlas_inc_coeliac"="prs.celiac", "atlas_inc_crohns"="prs.crohns")

if(n.cases/4>200){
  u.train <- merge(u.train, prs[,c("eid_20361",p.list[[dz]])], by="eid_20361")
  u.opti <- merge(u.opti, prs[,c("eid_20361",p.list[[dz]])], by="eid_20361")
  u.test <- merge(u.test, prs[,c("eid_20361",p.list[[dz]])], by="eid_20361")
} else{
  u.train <- merge(u.train, prs[,c("eid_20361",p.list[[dz]])], by="eid_20361")
  u.test <- merge(u.test, prs[,c("eid_20361",p.list[[dz]])], by="eid_20361")
}

## survival objects
surv.train <- Surv(u.train$fol, u.train[,dz])
## optimization 
if(n.cases/4>200){
  surv.opt <- Surv(u.opti$fol, u.opti[,dz])
}else{surv.opt <- Surv(u.train$fol, u.train[,dz])}
## testing
surv.test <- Surv(u.test$fol, u.test[,dz])


boot.samples <- lapply(1:1000, function(x) sample(nrow(u.test), round(nrow(u.test)*1), replace = T))

## mode px info model optimization and testing
if(n.cases/4>200){
  res.clin <- coxnet.optim.r(Train.data = u.opti,
                             pred.vec = clin.predictors,
                             Train.surv.data = surv.opt,
                             times = 1000,
                             boot.samples = boot.samples,
                             Test.data = u.test,
                             Test.surv.data = surv.test)
} else{
  res.clin <- coxnet.optim.r(Train.data = u.train,
                             pred.vec = clin.predictors,
                             Train.surv.data = surv.opt,
                             times = 1000,
                             boot.samples = boot.samples,
                             Test.data = u.test,
                             Test.surv.data = surv.test)

}
save(res.clin, file = paste0("../data_output/px_info_prs_models/px_info_res_",dz,"_",inc.yrs,"yrs",".RData"))

# run px info model + PRS
if(n.cases/4>200){
  res.clin.prs <- coxnet.optim.r(Train.data = u.opti,
                               pred.vec = c(clin.predictors,p.list[[dz]]),
                               Train.surv.data = surv.opt,
                               times = 1000,
                               boot.samples = boot.samples,
                               Test.data = u.test,
                               Test.surv.data = surv.test)
} else{
  res.clin.prs <- coxnet.optim.r(Train.data = u.train,
                               pred.vec = c(clin.predictors,p.list[[dz]]),
                               Train.surv.data = surv.opt,
                               times = 1000,
                               boot.samples = boot.samples,
                               Test.data = u.test,
                               Test.surv.data = surv.test)

}
save(res.clin.prs, file = paste0("../data_output/px_info_prs_models/px_info_prs_res_",dz,"_",inc.yrs,"yrs",".RData"))

###################################################
#### Script Name: ML_pipeline_protein_pxinfo_prediction_pred
#### Author: Julia Carrasco Zanini Sanchez
#### Description: script to ML framework for incident diseases in atlas 308 
#### Date : 23/03/2023

rm(list = ls())

## load packages 
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(RNOmni)
# library(glue)
library(DBI)
library(dbplyr)
library(magrittr)
# library(broom)
library(survival)
library(data.table)
library(glue)
library(parallel)
library(caret)
library(glmnet)
library(ROSE)
# library(odbc)

## define disease and restrict to incident cases within X years
args <- commandArgs(trailingOnly = T)

dz <- as.character(args[1])
dz.date <- as.character(args[2])
inc.yrs <- as.numeric(args[3])

paste0("Prediction:", dz, "; ", dz.date, " at ", inc.yrs, "years")

## set up working directory 
setwd("/home/ukbb_inc_dz_prediction/Explore_1536_Expansion_analyses/02_inc_dz_prediction/bin/") 

## read in data 
atlas308 <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_atlas308_annotations.txt", data.table = F)
cov <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_covariate_annotations.txt", data.table = F)

####----  Data Cleaning and wrangling   ----####
atlas308 <- as.data.frame(atlas308 %>% 
                            filter(eid_20361%in%cov$eid_20361) %>% 
                            pivot_wider(names_from = variable, values_from =value))

## phenotypes missing for patient-derived info model: bmi, smoking, alcohol, physical activity, family history 
pheno <- merge(cov[,-c(11:14,17)], atlas308, by="eid_20361")

## read in px info phenotypes
px.info <- read.delim("../../01_inc_dz_phenos/data_input/Patient_model_variables.txt",sep = "\t")
pheno <- merge(pheno, px.info, by="eid_20361")

## rename to olds names to match the rest fo the script
colnames(pheno)[which(colnames(pheno)=="assess_date")] <- c("date_baseline_assessment")

## process ethnicity
pheno$ethnicity <- ifelse(is.na(pheno$ge_eur),2,pheno$ge_eur)

## get protein data and metadata
ol.npx <- fread("../../00_protein_imputation/data_input/Imputed_1536_Expansion_NPX_proteins.txt", data.table = F)
ol.npx <- filter(ol.npx, exclusion.50==0 )

## exclude participants with missing age, sex
ol.npx <- ol.npx[complete.cases(ol.npx[,c("ge_sex", "ge_age")]),]

protein_list <- colnames(ol.npx)[-which(colnames(ol.npx)%in%c("eid_20361","ge_sex","ge_age",
                                                              "randomized_baseline","exclusion.50"))]

## merge phenotypes and proteins
ukbb <- merge(ol.npx, pheno, by="eid_20361")
## randomised baseline
ukbb <- ukbb %>% 
  filter(randomized_baseline==T)
## inverse rank normal transform npx
ukbb[,protein_list] <- sapply(ukbb[,protein_list], RankNorm)

## exclude prevalent cases and incident cases up to 6 months
## restrict to incident cases within the first X years
ukbb <- ukbb %>% 
  mutate(fol = difftime(!!as.name(dz.date),date_baseline_assessment, units = "weeks")/52.25) %>% 
  filter((fol>0.5&fol<inc.yrs)|is.na(!!as.name(dz.date)))

## censor (including censoring by death)
death.data <- read.delim("../../01_inc_dz_phenos/data_input/Death_dates_proteomic_participants.txt",sep = "\t")
death.data$date_of_death <- as.Date(death.data$date_of_death, format = "%Y-%m-%d")

ukbb <- merge(ukbb, death.data[,c("eid_20361","date_of_death")], by="eid_20361", all.x=T)

if(inc.yrs==10){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2020-12-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
} else if(inc.yrs==5){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2016-05-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
}else{stop("limit to incident casese within X years not defined")}

if(inc.yrs==10){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2020-12-31"),
                       difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                       ukbb$fol)
} else if(inc.yrs==5){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2016-05-31"),
                       difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                       ukbb$fol)
}else{stop("limit to incident casese within X years not defined")}

ukbb$fol <- ifelse(is.na(ukbb$fol_d), ukbb$fol, ukbb$fol_d)

## for sex specific outcome exclude males/females
if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_pri_breast",
           "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical",
           "atlas_inc_endometrial_hyper","atlas_inc_endometriosis","atlas_inc_female_genital_prolapse",
           "atlas_inc_leiomyoma","atlas_inc_menorrhagia",
           "atlas_inc_pid","atlas_inc_pmb")){
  # ukbb <- ukbb %>% filter(ge_sex=="female")
  ukbb <- ukbb %>% filter(ge_sex==0)
} else if(dz%in%c("atlas_inc_pri_prost","atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  # ukbb <- ukbb %>% filter(ge_sex=="male")
  ukbb <- ukbb %>% filter(ge_sex==1)
} else{ukbb <- ukbb}

## recode sex to numeric
# ukbb$ge_sex <- ifelse(ukbb$ge_sex=="male",1,
#                       ifelse(ukbb$ge_sex=="female",0,NA))

n.cases <- sum(ukbb[,dz]==1)

## add noise variables
rand.vars <- lapply(1:10, function(x){
  tmp <- runif(nrow(ukbb), min = min(ukbb$NPPB_20049), 
               max = max(ukbb$NPPB_20049))
  return(tmp)
})
rand.vars <- as.data.frame(do.call(cbind, rand.vars))
colnames(rand.vars) <- paste0("rand_var_",1:10)

ukbb <- cbind(ukbb, rand.vars)

####----    ML framework    ----####

## split into train/optimize/test or train/test depending on case number
ukbb[,dz] <- as.factor(ukbb[,dz])

if(n.cases/4>200){
  ## feature selection set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.5, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.opti <- ukbb[-trainIndex,]
  
  ## optimization and testing
  set.seed(3456) ## for reproducibility
  trainIndex2 <- createDataPartition(u.opti[,dz], p=0.5, list = F, times = 1)
  u.test <- u.opti[trainIndex2,]
  u.opti <- u.opti[-trainIndex2,]
  }else {
  ## feature selection set and test set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.7, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.test <- ukbb[-trainIndex,]
  }
  
####---- optimization and performance testing ----####

## function to run optimization and testings
source("../../../02_inc_dz_prediction/bin/functions/coxnet_ridge_optim_boot.R")

## px information model basic 
clin.predictors <- c("ge_age","ge_sex","bmi","ethnicity","smoke_never",
                     "smoke_previous","smoke_current","alcohol_daily",
                     "alcohol_4_week","alcohol_2_week","alcohol_3_month","alcohol_occassion",
                     "alcohol_never")

## define additional family history for specific diseases
if(!dz%in%c("atlas_inc_diabetes_t2","atlas_inc_pri_prost", "atlas_inc_copd","atlas_inc_isch_stroke",
            "atlas_inc_af","atlas_inc_av_block_3","atlas_inc_chd_nos","atlas_inc_hf", 
            "atlas_inc_peripheral_arterial_disease","atlas_inc_parkinsons","atlas_inc_pri_lung",
            "atlas_inc_pri_bowel", "atlas_inc_pri_breast","atlas_inc_pri_ovarian","atlas_inc_pri_uterine",
            "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
            "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
            "atlas_inc_hydrocele","atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb",
            "atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu")){
  clin.predictors <- clin.predictors
} else if(dz=="atlas_inc_diabetes_t2") {
  clin.predictors <- c(clin.predictors,"father_diabetes","mother_diabetes")
} else if(dz=="atlas_inc_pri_prost") {
  clin.predictors <- c(clin.predictors[-2],"father_prostate_cancer")
} else if(dz=="atlas_inc_copd") {
  clin.predictors <- c(clin.predictors,"father_COPD","mother_COPD")
} else if(dz%in%c("atlas_inc_isch_stroke","atlas_inc_af",
                  "atlas_inc_chd_nos","atlas_inc_hf","atlas_inc_peripheral_arterial_disease")) {
  clin.predictors <- c(clin.predictors,"father_heart_disease","mother_heart_disease")
} else if(dz=="atlas_inc_parkinsons") {
  clin.predictors <- c(clin.predictors,"father_parkinsons","mother_parkinsons")
} else if(dz%in%c("atlas_inc_pri_lung","atlas_inc_sec_lung")) {
  clin.predictors <- c(clin.predictors,"father_lung_cancer","mother_lung_cancer")
} else if(dz=="atlas_inc_pri_bowel") {
  clin.predictors <- c(clin.predictors,"father_bowel_cancer","mother_bowel_cancer")
} else if(dz=="atlas_inc_pri_breast") {
  clin.predictors <- c(clin.predictors[-2],"mother_breast_cancer")
} else if(dz=="atlas_inc_hypertension") {
  clin.predictors <- c(clin.predictors[-2],"mother_high_bp","father_high_bp")
} else if(dz=="atlas_inc_depression") {
  clin.predictors <- c(clin.predictors[-2],"mother_depression","father_depression")
} else if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_benign_ovary",
                  "atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
                  "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
                  "atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb")) {
  clin.predictors <- c(clin.predictors[-2])
} else if(dz%in%c("atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  clin.predictors <- c(clin.predictors[-2])
}


## read in PRS data 
prs <- read.delim("../../01_inc_dz_phenos/data_input/PRS_standard_ukb.txt", sep = "\t")

p.list <- list("atlas_inc_rha"="prs.ra","atlas_inc_diabetes_t2"="prs.t2d","atlas_inc_ulc_colitis"="prs.uc",
               "atlas_inc_vte_ex_pe"="prs.vte","atlas_inc_asthma"="prs.asthma","atlas_inc_af"="prs.af",
               "atlas_inc_pri_bowel"="prs.bowel_ca","atlas_inc_pri_breast"="prs.breast_ca",
               "atlas_inc_chd_nos"="prs.cad","atlas_inc_pri_ovarian"="prs.ovarian_ca",
               "atlas_inc_fracture_hip"="prs.bmd","atlas_inc_fracture_wrist"="prs.bmd",
               "atlas_inc_hypertension"="prs.hypertension","atlas_inc_isch_stroke"="prs.isch_stroke",
               "atlas_inc_pri_skin"="prs.skin_ca","atlas_inc_oa"="prs.oa","atlas_inc_parkinsons"="prs.parkinsons",          
               "atlas_inc_glaucoma"="prs.glaucoma","atlas_inc_pri_prost"="prs.prostate_ca",
               "atlas_inc_psoriasis"="prs.psoriasis")

if(n.cases/4>200){
  u.train <- merge(u.train, prs[,c("eid_20361",p.list[[dz]])], by="eid_20361")
  u.opti <- merge(u.opti, prs[,c("eid_20361",p.list[[dz]])], by="eid_20361")
  u.test <- merge(u.test, prs[,c("eid_20361",p.list[[dz]])], by="eid_20361")
} else{
  u.train <- merge(u.train, prs[,c("eid_20361",p.list[[dz]])], by="eid_20361")
  u.test <- merge(u.test, prs[,c("eid_20361",p.list[[dz]])], by="eid_20361")
}

## survival objects
surv.train <- Surv(u.train$fol, u.train[,dz])
## optimization 
if(n.cases/4>200){
  surv.opt <- Surv(u.opti$fol, u.opti[,dz])
}else{surv.opt <- Surv(u.train$fol, u.train[,dz])}
## testing
surv.test <- Surv(u.test$fol, u.test[,dz])

boot.samples <- lapply(1:1000, function(x) sample(nrow(u.test), round(nrow(u.test)*1), replace = T))

## mode px info model optimization and testing
if(n.cases/4>200){
  res.clin <- coxnet.optim.r(Train.data = u.opti,
                             pred.vec = clin.predictors,
                             Train.surv.data = surv.opt,
                             times = 1000,
                             boot.samples = boot.samples,
                             Test.data = u.test,
                             Test.surv.data = surv.test)
} else{
  res.clin <- coxnet.optim.r(Train.data = u.train,
                             pred.vec = clin.predictors,
                             Train.surv.data = surv.opt,
                             times = 1000,
                             boot.samples = boot.samples,
                             Test.data = u.test,
                             Test.surv.data = surv.test)

}
save(res.clin, file = paste0("../data_output/px_info_prs_models/px_info_res_",dz,"_",inc.yrs,"yrs",".RData"))

# run px info model + PRS
if(n.cases/4>200){
  res.clin.prs <- coxnet.optim.r(Train.data = u.opti,
                               pred.vec = c(clin.predictors,p.list[[dz]]),
                               Train.surv.data = surv.opt,
                               times = 1000,
                               boot.samples = boot.samples,
                               Test.data = u.test,
                               Test.surv.data = surv.test)
} else{
  res.clin.prs <- coxnet.optim.r(Train.data = u.train,
                               pred.vec = c(clin.predictors,p.list[[dz]]),
                               Train.surv.data = surv.opt,
                               times = 1000,
                               boot.samples = boot.samples,
                               Test.data = u.test,
                               Test.surv.data = surv.test)

}

save(res.clin.prs, file = paste0("../data_output/px_info_prs_models/px_info_prs_res_",dz,"_",inc.yrs,"yrs",".RData"))
###################################################
#### Script Name: ML_pipeline_protein_pxinfo_prediction_pred
#### Author: Julia Carrasco Zanini Sanchez
#### Description: script to ML framework for incident diseases in atlas 308 
#### Date : 23/03/2023

rm(list = ls())

## load packages 
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(RNOmni)
# library(glue)
library(DBI)
library(dbplyr)
library(magrittr)
# library(broom)
library(survival)
library(data.table)
library(glue)
library(parallel)
library(caret)
library(glmnet)
library(ROSE)
# library(odbc)

## define disease and restrict to incident cases within X years
args <- commandArgs(trailingOnly = T)

dz <- as.character(args[1])
dz.date <- as.character(args[2])
inc.yrs <- as.numeric(args[3])

paste0("Prediction:", dz, "; ", dz.date, " at ", inc.yrs, "years")

## set up working directory 
setwd("/home/ukbb_inc_dz_prediction/Explore_1536_Expansion_analyses/02_inc_dz_prediction/bin/") 

## read in data 
atlas308 <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_atlas308_annotations.txt", data.table = F)
cov <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_covariate_annotations.txt", data.table = F)

####----  Data Cleaning and wrangling   ----####
atlas308 <- as.data.frame(atlas308 %>% 
  filter(eid_20361%in%cov$eid_20361) %>% 
  pivot_wider(names_from = variable, values_from =value))

## phenotypes missing for patient-derived info model: bmi, smoking, alcohol, physical activity, family history 
pheno <- merge(cov[,-c(2:14,17)], atlas308, by="eid_20361")

## read in px info phenotypes
px.info <- read.delim("../../01_inc_dz_phenos/data_input/Patient_model_variables.txt",sep = "\t")
pheno <- merge(pheno, px.info, by="eid_20361")

## rename to olds names to match the rest fo the script
colnames(pheno)[which(colnames(pheno)=="assess_date")] <- c("date_baseline_assessment")

## process ethnicity
pheno$ethnicity <- ifelse(is.na(pheno$ge_eur),2,pheno$ge_eur)

## get protein data and metadata
ol.npx <- fread("../../00_protein_imputation/data_input/Imputed_1536_Expansion_NPX_proteins.txt", data.table = F)
ol.npx <- filter(ol.npx, exclusion.50==0 )

## exclude participants with missing age, sex
ol.npx <- ol.npx[complete.cases(ol.npx[,c("ge_sex", "ge_age")]),]

protein_list <- colnames(ol.npx)[-which(colnames(ol.npx)%in%c("eid_20361","ge_sex","ge_age",
                                                              "randomized_baseline","exclusion.50"))]

## merge phenotypes and proteins
ukbb <- merge(ol.npx, pheno, by="eid_20361")
## randomised baseline
ukbb <-
  ukbb %>% 
  filter(randomized_baseline==T|!!as.name(dz)==1)

## inverse rank normal transform npx
ukbb[,protein_list] <- sapply(ukbb[,protein_list], RankNorm)

## read in biomarkers
biom <- read.delim("../../01_inc_dz_phenos/data_input/biomarkers_clean_230623_ukb.txt", sep = "\ ")
biom <- biom %>% 
  select(-ge_sex, -age_0)

biomarkers <- colnames(biom)[-ncol(biom)]
# ukbb <- merge(ukbb, biom, by="eid_20361")

## exclude prevalent cases and incident cases up to 6 months
## restrict to incident cases within the first X years
ukbb <- ukbb %>% 
  mutate(fol = difftime(!!as.name(dz.date),date_baseline_assessment, units = "weeks")/52.25) %>% 
  filter((fol>0.5&fol<inc.yrs)|is.na(!!as.name(dz.date)))

## censor (including censoring by death)
death.data <- read.delim("../../01_inc_dz_phenos/data_input/Death_dates_proteomic_participants.txt",sep = "\t")
death.data$date_of_death <- as.Date(death.data$date_of_death, format = "%Y-%m-%d")

ukbb <- merge(ukbb, death.data[,c("eid_20361","date_of_death")], by="eid_20361", all.x=T)

if(inc.yrs==10){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2020-12-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
} else if(inc.yrs==5){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2016-05-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
}else{stop("limit to incident casese within X years not defined")}

if(inc.yrs==10){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2020-12-31"),
                         difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                        ukbb$fol)
} else if(inc.yrs==5){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2016-05-31"),
                          difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                          ukbb$fol)
}else{stop("limit to incident casese within X years not defined")}

ukbb$fol <- ifelse(is.na(ukbb$fol_d), ukbb$fol, ukbb$fol_d)

## for sex specific outcome exclude males/females
if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_pri_breast",
           "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical",
           "atlas_inc_endometrial_hyper","atlas_inc_endometriosis","atlas_inc_female_genital_prolapse",
           "atlas_inc_leiomyoma","atlas_inc_menorrhagia",
           "atlas_inc_pid","atlas_inc_pmb")){
  # ukbb <- ukbb %>% filter(ge_sex=="female")
  ukbb <- ukbb %>% filter(ge_sex==0)
} else if(dz%in%c("atlas_inc_pri_prost","atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  # ukbb <- ukbb %>% filter(ge_sex=="male")
  ukbb <- ukbb %>% filter(ge_sex==1)
} else{ukbb <- ukbb}

n.cases <- sum(ukbb[,dz]==1)

## add noise variables
rand.vars <- lapply(1:10, function(x){
  tmp <- runif(nrow(ukbb), min = min(ukbb$NPPB_20049), 
               max = max(ukbb$NPPB_20049))
  return(tmp)
})
rand.vars <- as.data.frame(do.call(cbind, rand.vars))
colnames(rand.vars) <- paste0("rand_var_",1:10)

ukbb <- cbind(ukbb, rand.vars)

####----    ML framework    ----####

## split into train/optimize/test or train/test depending on case number
ukbb[,dz] <- as.factor(ukbb[,dz])

if(n.cases/4>200){
  ## feature selection set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.5, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.opti <- ukbb[-trainIndex,]
  
  ## optimization and testing
  set.seed(3456) ## for reproducibility
  trainIndex2 <- createDataPartition(u.opti[,dz], p=0.5, list = F, times = 1)
  u.test <- u.opti[trainIndex2,]
  u.opti <- u.opti[-trainIndex2,]
  ## merge biomarkers  (note have preserved the same split as for other models)
  u.train <- merge(u.train, biom, by="eid_20361")
  u.opti <- merge(u.opti, biom, by="eid_20361")
  u.test <- merge(u.test, biom, by="eid_20361")
  }else {
  ## feature selection set and test set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.7, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.test <- ukbb[-trainIndex,]
  u.train <- merge(u.train, biom, by="eid_20361")
  u.test <- merge(u.test, biom, by="eid_20361")
  }
  
## draw 200 random subsamples of train set
jj <- lapply(1:200, function(x) sample(nrow(u.train), round(nrow(u.train)*0.50), replace=F))

##-- feature selection  --##
las.morb <- mclapply(1:length(jj), function(i){
  print(i) ## just for test purposes
  tmp <- train(u.train[jj[[i]],c(biomarkers,colnames(rand.vars))],
               u.train[jj[[i]],dz],
               metric = "Kappa",
               method="glmnet",
               family="binomial",
               tuneGrid=as.data.frame(expand_grid(
                 alpha=1,lambda=10^-seq(4,.25,-.5))),
               trControl= trainControl(method = "repeatedcv",
                                       number = 3, repeats = 5,
                                       sampling = "rose",
                                       allowParallel = T)
  )
  cf <- as.matrix(coef(tmp$finalModel,
                       s = tmp$finalModel$lambdaOpt))
  gc()
  # return(list("cf" = cf , "cv.tune" = tmp))
  return(cf)
}, mc.cores = 12, mc.allow.recursive = T)

save(las.morb, file=paste0("../data_output/feature_selection/biomarker_fs_",dz,"_",inc.yrs,"yrs",".RData"))

## generate feature selection ranking 
# load(paste0("../data_output/feature_selection/biomarker_fs_",dz,"_",inc.yrs,"yrs",".RData"))
p.select <- matrix(nrow = nrow(las.morb[[1]]), ncol = length(las.morb))

for(i in 1:length(las.morb)){
  tryCatch({
    tmp <- las.morb[[i]][,1]
    if(is.null(tmp)){
      return(tmp <- rep(0,nrow(p.select)))
    }else {
      NULL
    }
  }, error=function(e){
    NULL
  })
  p.select[,i] <- tmp
}
row.names(p.select) <- names(las.morb[[1]][,1])

## count matrix
p.select <-
  p.select %>% 
  as_tibble(rownames = "biomarker.id") %>% 
  filter(biomarker.id!="(Intercept)") %>% 
  mutate(select = abs(rowSums(across(where(is.numeric))))) %>% 
  arrange(desc(select)) 
  
p.select <- as.data.frame(p.select)

p.select$select.perc <- (p.select$select)/max(p.select$select)

####---- optimization and performance testing ----####

## survival objects
surv.train <- Surv(u.train$fol, u.train[,dz])
## optimization 
if(n.cases/4>200){
  surv.opt <- Surv(u.opti$fol, u.opti[,dz])
}else{surv.opt <- Surv(u.train$fol, u.train[,dz])}
## testing
surv.test <- Surv(u.test$fol, u.test[,dz])

## function to run optimization and testings
source("../../../02_inc_dz_prediction/bin/functions/coxnet_ridge_optim_boot.R")

## px information model basic 
clin.predictors <- c("ge_age","ge_sex","bmi","ethnicity","smoke_never",
                     "smoke_previous","smoke_current","alcohol_daily",
                     "alcohol_4_week","alcohol_2_week","alcohol_3_month","alcohol_occassion",
                     "alcohol_never")

## define additional family history for specific diseases
if(!dz%in%c("atlas_inc_diabetes_t2","atlas_inc_pri_prost", "atlas_inc_copd","atlas_inc_isch_stroke",
            "atlas_inc_af","atlas_inc_av_block_3","atlas_inc_chd_nos","atlas_inc_hf", 
            "atlas_inc_peripheral_arterial_disease","atlas_inc_parkinsons","atlas_inc_pri_lung",
            "atlas_inc_pri_bowel", "atlas_inc_pri_breast","atlas_inc_pri_ovarian","atlas_inc_pri_uterine",
            "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
            "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
            "atlas_inc_hydrocele","atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb",
            "atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu")){
  clin.predictors <- clin.predictors
} else if(dz=="atlas_inc_diabetes_t2") {
  clin.predictors <- c(clin.predictors,"father_diabetes","mother_diabetes")
} else if(dz=="atlas_inc_pri_prost") {
  clin.predictors <- c(clin.predictors[-2],"father_prostate_cancer")
} else if(dz=="atlas_inc_copd") {
  clin.predictors <- c(clin.predictors,"father_COPD","mother_COPD")
} else if(dz%in%c("atlas_inc_isch_stroke","atlas_inc_af",
                  "atlas_inc_chd_nos","atlas_inc_hf","atlas_inc_peripheral_arterial_disease")) {
  clin.predictors <- c(clin.predictors,"father_heart_disease","mother_heart_disease")
} else if(dz=="atlas_inc_parkinsons") {
  clin.predictors <- c(clin.predictors,"father_parkinsons","mother_parkinsons")
} else if(dz%in%c("atlas_inc_pri_lung","atlas_inc_sec_lung")) {
  clin.predictors <- c(clin.predictors,"father_lung_cancer","mother_lung_cancer")
} else if(dz=="atlas_inc_pri_bowel") {
  clin.predictors <- c(clin.predictors,"father_bowel_cancer","mother_bowel_cancer")
} else if(dz=="atlas_inc_pri_breast") {
  clin.predictors <- c(clin.predictors[-2],"mother_breast_cancer")
} else if(dz=="atlas_inc_hypertension") {
  clin.predictors <- c(clin.predictors[-2],"mother_high_bp","father_high_bp")
} else if(dz=="atlas_inc_depression") {
  clin.predictors <- c(clin.predictors[-2],"mother_depression","father_depression")
} else if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_benign_ovary",
                  "atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
                  "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
                  "atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb")) {
  clin.predictors <- c(clin.predictors[-2])
} else if(dz%in%c("atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  clin.predictors <- c(clin.predictors[-2])
}


boot.samples <- lapply(1:1000, function(x) sample(nrow(u.test), round(nrow(u.test)*1), replace = T))

## mode px info model optimization and testing
if(n.cases/4>200){
  res.clin <- coxnet.optim.r(Train.data = u.opti,
                             pred.vec = clin.predictors,
                             Train.surv.data = surv.opt,
                             times = 1000, 
                             boot.samples = boot.samples,
                             Test.data = u.test,
                             Test.surv.data = surv.test)
} else{
  res.clin <- coxnet.optim.r(Train.data = u.train,
                             pred.vec = clin.predictors,
                             Train.surv.data = surv.opt,
                             times = 1000, 
                             boot.samples = boot.samples,
                             Test.data = u.test,
                             Test.surv.data = surv.test)
  
}
save(res.clin, file = paste0("../data_output/biomarker_models_update/px_info_res_",dz,"_",inc.yrs,"yrs",".RData"))


## selection threshold based on noise variables
sel.t <- p.select$select.perc[grep("rand", p.select$biomarker.id)[1]]

biom.opti <- p.select %>% 
  filter(select.perc>sel.t) %>% 
  top_n(20, select.perc)
biom.opti <- biom.opti$biomarker.id
if(length(biom.opti)==0){
  biom.opti <- p.select[p.select$select.per==sel.t]
  biom.opti <- biom.opti$biomarker.id
}

# run px info model + bioamrkers selected above noise
if(n.cases/4>200){
    res.clin.biom.all <- coxnet.optim.r(Train.data = u.opti,
                                   pred.vec = c(clin.predictors,biom.opti),
                                   Train.surv.data = surv.opt,
                                   times = 1000, 
                                   boot.samples = boot.samples,
                                   Test.data = u.test,
                                   Test.surv.data = surv.test)
    } else{
      res.clin.biom.all <- coxnet.optim.r(Train.data = u.train,
                                   pred.vec = c(clin.predictors,biom.opti),
                                   Train.surv.data = surv.opt,
                                   times = 1000, 
                                   boot.samples = boot.samples,
                                   Test.data = u.test,
                                   Test.surv.data = surv.test)
    }
save(res.clin.biom.all, file = paste0("../data_output/biomarker_models_update/px_info_all_biom_above_noise_res_",dz,"_",inc.yrs,"yrs",".RData"))

# run px info model + biomarkers selected above noise in train set
if(n.cases/4>200){
  res.clin.biom.all <- coxnet.optim.r(Train.data = u.train,
                                     pred.vec = c(clin.predictors,biom.opti),
                                     Train.surv.data = surv.train,
                                     times = 1000, 
                                     boot.samples = boot.samples,
                                     Test.data = u.test,
                                     Test.surv.data = surv.test)
} else {res.clin.biom.all <- res.clin.biom.all}

## biomarkers top 5 and 10 only based on weights*selection_percentage
opti.biom <- res.clin.biom.all$opt.coefficients[which(row.names(res.clin.biom.all$opt.coefficients)%in%biomarkers),]
opti.biom <- opti.biom %>%
  as_tibble(rownames = "biom") %>%
  mutate(abs.coef = abs(value)) %>% 
  filter(abs.coef!=0)

## merge with selection percentage
opti.biom <- merge(opti.biom, p.select[,c("biomarker.id","select","select.perc")],
                    by.x="biom", by.y="biomarker.id")

if(nrow(opti.biom)>=10){
  opti.biom10 <- opti.biom %>%
    mutate(sel.opti.score=abs.coef*select) %>%
    top_n(10, sel.opti.score)
  biom.opti10 <- opti.biom10$biom
  
  ## run px info model + 10 top bioamarker based on weightsXselection.perc
  if(n.cases/4>200){
    res.clin.biom10 <- coxnet.optim.r(Train.data = u.opti,
                                     pred.vec = c(clin.predictors,biom.opti10),
                                     Train.surv.data = surv.opt,
                                     times = 1000, 
                                     boot.samples = boot.samples,
                                     Test.data = u.test,
                                     Test.surv.data = surv.test)
  } else{
    res.clin.biom10 <- coxnet.optim.r(Train.data = u.train,
                                     pred.vec = c(clin.predictors,biom.opti10),
                                     Train.surv.data = surv.opt,
                                     times = 1000, 
                                     boot.samples = boot.samples,
                                     Test.data = u.test,
                                     Test.surv.data = surv.test)
    
  }
  save(res.clin.biom10, file = paste0("../data_output/biomarker_models_update/px_info_top10_biom_res_",dz,"_",inc.yrs,"yrs",".RData"))
  
} else{
  biom.opti10 <- opti.biom$biom
  ## run px info model + 10 top bioamarker based on weightsXselection.perc
  if(n.cases/4>200){
    res.clin.biom10 <- coxnet.optim.r(Train.data = u.opti,
                                      pred.vec = c(clin.predictors,biom.opti10),
                                      Train.surv.data = surv.opt,
                                      times = 1000, 
                                      boot.samples = boot.samples,
                                      Test.data = u.test,
                                      Test.surv.data = surv.test)
  } else{
    res.clin.biom10 <- coxnet.optim.r(Train.data = u.train,
                                      pred.vec = c(clin.predictors,biom.opti10),
                                      Train.surv.data = surv.opt,
                                      times = 1000, 
                                      boot.samples = boot.samples,
                                      Test.data = u.test,
                                      Test.surv.data = surv.test)
    
  }
  save(res.clin.biom10, file = paste0("../data_output/biomarker_models_update/px_info_top10_biom_res_",dz,"_",inc.yrs,"yrs",".RData"))
  
}

## run px info model + 5 top biomarker based on weightsXselection.perc
if(nrow(opti.biom)>=5){
  opti.biom <- opti.biom %>%
    mutate(sel.opti.score=abs.coef*select.perc) %>%
    top_n(5, sel.opti.score)
  biom.opti <- opti.biom$biom
  } else {
    biom.opti <- opti.biom$biom
}


if(n.cases/4>200){
  res.clin.biom <- coxnet.optim.r(Train.data = u.opti,
                                 pred.vec = c(clin.predictors,biom.opti),
                                 Train.surv.data = surv.opt,
                                 times = 1000, 
                                 boot.samples = boot.samples,
                                 Test.data = u.test,
                                 Test.surv.data = surv.test)
} else {
  res.clin.biom <- coxnet.optim.r(Train.data = u.train,
                                 pred.vec = c(clin.predictors,biom.opti),
                                 Train.surv.data = surv.opt,
                                 times = 1000, 
                                 boot.samples = boot.samples,
                                 Test.data = u.test,
                                 Test.surv.data = surv.test)
  
}
save(res.clin.biom, file = paste0("../data_output/biomarker_models_update/px_info_top5_biom_res_",dz,"_",inc.yrs,"yrs",".RData"))

## optimization on top 5 biomarkers
if(n.cases/4>200){
  res.biom <- coxnet.optim.r(Train.data = u.opti,
                           pred.vec = biom.opti,
                           Train.surv.data = surv.opt,
                           times = 1000, 
                           boot.samples = boot.samples, 
                           Test.data = u.test,
                           Test.surv.data = surv.test)
} else{
  res.biom <- coxnet.optim.r(Train.data = u.train,
                           pred.vec = biom.opti,
                           Train.surv.data = surv.opt,
                           times = 1000, 
                           boot.samples = boot.samples,
                           Test.data = u.test,
                           Test.surv.data = surv.test)
  
}
## add feature selection to the results
res.biom$feature.selection <- p.select

## add case numbers for training and testing sets
if((n.cases/4)>200){
  case.num <- data.frame(train.cases = sum(u.train[,dz]==1),
                         opti.cases = sum(u.opti[,dz]==1),
                         test.cases = sum(u.test[,dz]==1),
                         train.controls = sum(u.train[,dz]!=1),
                         opti.controls = sum(u.opti[,dz]!=1),
                         test.controls = sum(u.test[,dz]!=1))
} else{
  case.num <- data.frame(train.cases = sum(u.train[,dz]==1),
                         opti.cases = NA,
                         test.cases = sum(u.test[,dz]==1),
                         train.controls = sum(u.train[,dz]!=1),
                         opti.controls = NA,
                         test.controls = sum(u.test[,dz]!=1))
}
res.biom$case.num <- case.num
## save results 
save(res.biom, file = paste0("../data_output/biomarker_models_update/biomarker_top5_res_",dz,"_",inc.yrs,"yrs",".RData"))

###################################################
#### Script Name: ML_pipeline_protein_pxinfo_prediction_pred
#### Author: Julia Carrasco Zanini Sanchez
#### Description: script to ML framework for incident diseases in atlas 308 
#### Date : 23/03/2023

rm(list = ls())

## load packages 
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(RNOmni)
# library(glue)
library(DBI)
library(dbplyr)
library(magrittr)
# library(broom)
library(survival)
library(data.table)
library(glue)
library(parallel)
library(caret)
library(glmnet)
library(ROSE)
# library(odbc)

## define disease and restrict to incident cases within X years
args <- commandArgs(trailingOnly = T)

dz <- as.character(args[1])
dz.date <- as.character(args[2])
inc.yrs <- as.numeric(args[3])

paste0("Prediction:", dz, "; ", dz.date, " at ", inc.yrs, "years")

## set up working directory 
setwd("/home/ukbb_inc_dz_prediction/Explore_1536_Expansion_analyses/02_inc_dz_prediction/bin/") 

## read in data 
atlas308 <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_atlas308_annotations.txt", data.table = F)
cov <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_covariate_annotations.txt", data.table = F)

####----  Data Cleaning and wrangling   ----####
atlas308 <- as.data.frame(atlas308 %>% 
  filter(eid_20361%in%cov$eid_20361) %>% 
  pivot_wider(names_from = variable, values_from =value))

## phenotypes missing for patient-derived info model: bmi, smoking, alcohol, physical activity, family history 
pheno <- merge(cov[,-c(2:14,17)], atlas308, by="eid_20361")

## read in px info phenotypes
px.info <- read.delim("../../01_inc_dz_phenos/data_input/Patient_model_variables.txt",sep = "\t")
pheno <- merge(pheno, px.info, by="eid_20361")

## rename to olds names to match the rest fo the script
colnames(pheno)[which(colnames(pheno)=="assess_date")] <- c("date_baseline_assessment")

## process ethnicity
pheno$ethnicity <- ifelse(is.na(pheno$ge_eur),2,pheno$ge_eur)

## get protein data and metadata
ol.npx <- fread("../../00_protein_imputation/data_input/Imputed_1536_Expansion_NPX_proteins.txt", data.table = F)
ol.npx <- filter(ol.npx, exclusion.50==0 )

## exclude participants with missing age, sex
ol.npx <- ol.npx[complete.cases(ol.npx[,c("ge_sex", "ge_age")]),]

protein_list <- colnames(ol.npx)[-which(colnames(ol.npx)%in%c("eid_20361","ge_sex","ge_age",
                                                              "randomized_baseline","exclusion.50"))]

## merge phenotypes and proteins
ukbb <- merge(ol.npx, pheno, by="eid_20361")
## randomised baseline
ukbb <- ukbb %>% 
  filter(randomized_baseline==T)
## inverse rank normal transform npx
ukbb[,protein_list] <- sapply(ukbb[,protein_list], RankNorm)

## read in biomarkers
biom <- read.delim("../../01_inc_dz_phenos/data_input/biomarkers_clean_230623_ukb.txt", sep = "\ ")
biom <- biom %>% 
  select(-ge_sex, -age_0)

biomarkers <- colnames(biom)[-ncol(biom)]
# ukbb <- merge(ukbb, biom, by="eid_20361")

## exclude prevalent cases and incident cases up to 6 months
## restrict to incident cases within the first X years
ukbb <- ukbb %>% 
  mutate(fol = difftime(!!as.name(dz.date),date_baseline_assessment, units = "weeks")/52.25) %>% 
  filter((fol>0.5&fol<inc.yrs)|is.na(!!as.name(dz.date)))

## censor (including censoring by death)
death.data <- read.delim("../../01_inc_dz_phenos/data_input/Death_dates_proteomic_participants.txt",sep = "\t")
death.data$date_of_death <- as.Date(death.data$date_of_death, format = "%Y-%m-%d")

ukbb <- merge(ukbb, death.data[,c("eid_20361","date_of_death")], by="eid_20361", all.x=T)

if(inc.yrs==10){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2020-12-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
} else if(inc.yrs==5){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2016-05-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
}else{stop("limit to incident casese within X years not defined")}

if(inc.yrs==10){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2020-12-31"),
                         difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                        ukbb$fol)
} else if(inc.yrs==5){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2016-05-31"),
                          difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                          ukbb$fol)
}else{stop("limit to incident casese within X years not defined")}

ukbb$fol <- ifelse(is.na(ukbb$fol_d), ukbb$fol, ukbb$fol_d)

## for sex specific outcome exclude males/females
if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_pri_breast",
           "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical",
           "atlas_inc_endometrial_hyper","atlas_inc_endometriosis","atlas_inc_female_genital_prolapse",
           "atlas_inc_leiomyoma","atlas_inc_menorrhagia",
           "atlas_inc_pid","atlas_inc_pmb")){
  # ukbb <- ukbb %>% filter(ge_sex=="female")
  ukbb <- ukbb %>% filter(ge_sex==0)
} else if(dz%in%c("atlas_inc_pri_prost","atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  # ukbb <- ukbb %>% filter(ge_sex=="male")
  ukbb <- ukbb %>% filter(ge_sex==1)
} else{ukbb <- ukbb}

n.cases <- sum(ukbb[,dz]==1)

## add noise variables
rand.vars <- lapply(1:10, function(x){
  tmp <- runif(nrow(ukbb), min = min(ukbb$NPPB_20049), 
               max = max(ukbb$NPPB_20049))
  return(tmp)
})
rand.vars <- as.data.frame(do.call(cbind, rand.vars))
colnames(rand.vars) <- paste0("rand_var_",1:10)

ukbb <- cbind(ukbb, rand.vars)

####----    ML framework    ----####

## split into train/optimize/test or train/test depending on case number
ukbb[,dz] <- as.factor(ukbb[,dz])

if(n.cases/4>200){
  ## feature selection set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.5, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.opti <- ukbb[-trainIndex,]
  
  ## optimization and testing
  set.seed(3456) ## for reproducibility
  trainIndex2 <- createDataPartition(u.opti[,dz], p=0.5, list = F, times = 1)
  u.test <- u.opti[trainIndex2,]
  u.opti <- u.opti[-trainIndex2,]
  ## merge biomarkers  (note have preserved the same split as for other models)
  u.train <- merge(u.train, biom, by="eid_20361")
  u.opti <- merge(u.opti, biom, by="eid_20361")
  u.test <- merge(u.test, biom, by="eid_20361")
  }else {
  ## feature selection set and test set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.7, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.test <- ukbb[-trainIndex,]
  u.train <- merge(u.train, biom, by="eid_20361")
  u.test <- merge(u.test, biom, by="eid_20361")
  }
  
## draw 200 random subsamples of train set
jj <- lapply(1:200, function(x) sample(nrow(u.train), round(nrow(u.train)*0.50), replace=F))

##-- feature selection  --##
las.morb <- mclapply(1:length(jj), function(i){
  print(i) ## just for test purposes
  tmp <- train(u.train[jj[[i]],c(biomarkers,colnames(rand.vars))],
               u.train[jj[[i]],dz],
               metric = "Kappa",
               method="glmnet",
               family="binomial",
               tuneGrid=as.data.frame(expand_grid(
                 alpha=1,lambda=10^-seq(4,.25,-.5))),
               trControl= trainControl(method = "repeatedcv",
                                       number = 3, repeats = 5,
                                       sampling = "rose",
                                       allowParallel = T)
  )
  cf <- as.matrix(coef(tmp$finalModel,
                       s = tmp$finalModel$lambdaOpt))
  gc()
  # return(list("cf" = cf , "cv.tune" = tmp))
  return(cf)
}, mc.cores = 12, mc.allow.recursive = T)

save(las.morb, file=paste0("../data_output/feature_selection/biomarker_fs_",dz,"_",inc.yrs,"yrs",".RData"))

## generate feature selection ranking 
# load(paste0("../data_output/feature_selection/biomarker_fs_",dz,"_",inc.yrs,"yrs",".RData"))
p.select <- matrix(nrow = nrow(las.morb[[1]]), ncol = length(las.morb))

for(i in 1:length(las.morb)){
  tryCatch({
    tmp <- las.morb[[i]][,1]
    if(is.null(tmp)){
      return(tmp <- rep(0,nrow(p.select)))
    }else {
      NULL
    }
  }, error=function(e){
    NULL
  })
  p.select[,i] <- tmp
}
row.names(p.select) <- names(las.morb[[1]][,1])

## count matrix
p.select <-
  p.select %>% 
  as_tibble(rownames = "biomarker.id") %>% 
  filter(biomarker.id!="(Intercept)") %>% 
  mutate(select = abs(rowSums(across(where(is.numeric))))) %>% 
  arrange(desc(select)) 
  
p.select <- as.data.frame(p.select)

p.select$select.perc <- (p.select$select)/max(p.select$select)

####---- optimization and performance testing ----####

## survival objects
surv.train <- Surv(u.train$fol, u.train[,dz])
## optimization 
if(n.cases/4>200){
  surv.opt <- Surv(u.opti$fol, u.opti[,dz])
}else{surv.opt <- Surv(u.train$fol, u.train[,dz])}
## testing
surv.test <- Surv(u.test$fol, u.test[,dz])

## function to run optimization and testings
source("../../../02_inc_dz_prediction/bin/functions/coxnet_ridge_optim_boot.R")

## px information model basic 
clin.predictors <- c("ge_age","ge_sex","bmi","ethnicity","smoke_never",
                     "smoke_previous","smoke_current","alcohol_daily",
                     "alcohol_4_week","alcohol_2_week","alcohol_3_month","alcohol_occassion",
                     "alcohol_never")

## define additional family history for specific diseases
if(!dz%in%c("atlas_inc_diabetes_t2","atlas_inc_pri_prost", "atlas_inc_copd","atlas_inc_isch_stroke",
            "atlas_inc_af","atlas_inc_av_block_3","atlas_inc_chd_nos","atlas_inc_hf", 
            "atlas_inc_peripheral_arterial_disease","atlas_inc_parkinsons","atlas_inc_pri_lung",
            "atlas_inc_pri_bowel", "atlas_inc_pri_breast","atlas_inc_pri_ovarian","atlas_inc_pri_uterine",
            "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
            "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
            "atlas_inc_hydrocele","atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb",
            "atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu")){
  clin.predictors <- clin.predictors
} else if(dz=="atlas_inc_diabetes_t2") {
  clin.predictors <- c(clin.predictors,"father_diabetes","mother_diabetes")
} else if(dz=="atlas_inc_pri_prost") {
  clin.predictors <- c(clin.predictors[-2],"father_prostate_cancer")
} else if(dz=="atlas_inc_copd") {
  clin.predictors <- c(clin.predictors,"father_COPD","mother_COPD")
} else if(dz%in%c("atlas_inc_isch_stroke","atlas_inc_af",
                  "atlas_inc_chd_nos","atlas_inc_hf","atlas_inc_peripheral_arterial_disease")) {
  clin.predictors <- c(clin.predictors,"father_heart_disease","mother_heart_disease")
} else if(dz=="atlas_inc_parkinsons") {
  clin.predictors <- c(clin.predictors,"father_parkinsons","mother_parkinsons")
} else if(dz%in%c("atlas_inc_pri_lung","atlas_inc_sec_lung")) {
  clin.predictors <- c(clin.predictors,"father_lung_cancer","mother_lung_cancer")
} else if(dz=="atlas_inc_pri_bowel") {
  clin.predictors <- c(clin.predictors,"father_bowel_cancer","mother_bowel_cancer")
} else if(dz=="atlas_inc_pri_breast") {
  clin.predictors <- c(clin.predictors[-2],"mother_breast_cancer")
} else if(dz=="atlas_inc_hypertension") {
  clin.predictors <- c(clin.predictors[-2],"mother_high_bp","father_high_bp")
} else if(dz=="atlas_inc_depression") {
  clin.predictors <- c(clin.predictors[-2],"mother_depression","father_depression")
} else if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_benign_ovary",
                  "atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
                  "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
                  "atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb")) {
  clin.predictors <- c(clin.predictors[-2])
} else if(dz%in%c("atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  clin.predictors <- c(clin.predictors[-2])
}


boot.samples <- lapply(1:1000, function(x) sample(nrow(u.test), round(nrow(u.test)*1), replace = T))

## mode px info model optimization and testing
if(n.cases/4>200){
  res.clin <- coxnet.optim.r(Train.data = u.opti,
                             pred.vec = clin.predictors,
                             Train.surv.data = surv.opt,
                             times = 1000, 
                             boot.samples = boot.samples,
                             Test.data = u.test,
                             Test.surv.data = surv.test)
} else{
  res.clin <- coxnet.optim.r(Train.data = u.train,
                             pred.vec = clin.predictors,
                             Train.surv.data = surv.opt,
                             times = 1000, 
                             boot.samples = boot.samples,
                             Test.data = u.test,
                             Test.surv.data = surv.test)
  
}
save(res.clin, file = paste0("../data_output/biomarker_models_update/px_info_res_",dz,"_",inc.yrs,"yrs",".RData"))


## selection threshold based on noise variables
sel.t <- p.select$select.perc[grep("rand", p.select$biomarker.id)[1]]

biom.opti <- p.select %>% 
  filter(select.perc>sel.t) %>% 
  top_n(20, select.perc)
biom.opti <- biom.opti$biomarker.id
if(length(biom.opti)==0){
  biom.opti <- p.select[p.select$select.per==sel.t]
  biom.opti <- biom.opti$biomarker.id
}

# run px info model + bioamrkers selected above noise
if(n.cases/4>200){
    res.clin.biom.all <- coxnet.optim.r(Train.data = u.opti,
                                   pred.vec = c(clin.predictors,biom.opti),
                                   Train.surv.data = surv.opt,
                                   times = 1000, 
                                   boot.samples = boot.samples,
                                   Test.data = u.test,
                                   Test.surv.data = surv.test)
    } else{
      res.clin.biom.all <- coxnet.optim.r(Train.data = u.train,
                                   pred.vec = c(clin.predictors,biom.opti),
                                   Train.surv.data = surv.opt,
                                   times = 1000, 
                                   boot.samples = boot.samples,
                                   Test.data = u.test,
                                   Test.surv.data = surv.test)
    }
save(res.clin.biom.all, file = paste0("../data_output/biomarker_models_update/px_info_all_biom_above_noise_res_",dz,"_",inc.yrs,"yrs",".RData"))

# run px info model + biomarkers selected above noise in train set
if(n.cases/4>200){
  res.clin.biom.all <- coxnet.optim.r(Train.data = u.train,
                                     pred.vec = c(clin.predictors,biom.opti),
                                     Train.surv.data = surv.train,
                                     times = 1000, 
                                     boot.samples = boot.samples,
                                     Test.data = u.test,
                                     Test.surv.data = surv.test)
} else {res.clin.biom.all <- res.clin.biom.all}

## biomarkers top 5 and 10 only based on weights*selection_percentage
opti.biom <- res.clin.biom.all$opt.coefficients[which(row.names(res.clin.biom.all$opt.coefficients)%in%biomarkers),]
opti.biom <- opti.biom %>%
  as_tibble(rownames = "biom") %>%
  mutate(abs.coef = abs(value)) %>% 
  filter(abs.coef!=0)

## merge with selection percentage
opti.biom <- merge(opti.biom, p.select[,c("biomarker.id","select","select.perc")],
                    by.x="biom", by.y="biomarker.id")

if(nrow(opti.biom)>=10){
  opti.biom10 <- opti.biom %>%
    mutate(sel.opti.score=abs.coef*select) %>%
    top_n(10, sel.opti.score)
  biom.opti10 <- opti.biom10$biom
  
  ## run px info model + 10 top bioamarker based on weightsXselection.perc
  if(n.cases/4>200){
    res.clin.biom10 <- coxnet.optim.r(Train.data = u.opti,
                                     pred.vec = c(clin.predictors,biom.opti10),
                                     Train.surv.data = surv.opt,
                                     times = 1000, 
                                     boot.samples = boot.samples,
                                     Test.data = u.test,
                                     Test.surv.data = surv.test)
  } else{
    res.clin.biom10 <- coxnet.optim.r(Train.data = u.train,
                                     pred.vec = c(clin.predictors,biom.opti10),
                                     Train.surv.data = surv.opt,
                                     times = 1000, 
                                     boot.samples = boot.samples,
                                     Test.data = u.test,
                                     Test.surv.data = surv.test)
    
  }
  save(res.clin.biom10, file = paste0("../data_output/biomarker_models_update/px_info_top10_biom_res_",dz,"_",inc.yrs,"yrs",".RData"))
  
} else{
  biom.opti10 <- opti.biom$biom
  ## run px info model + 10 top bioamarker based on weightsXselection.perc
  if(n.cases/4>200){
    res.clin.biom10 <- coxnet.optim.r(Train.data = u.opti,
                                      pred.vec = c(clin.predictors,biom.opti10),
                                      Train.surv.data = surv.opt,
                                      times = 1000, 
                                      boot.samples = boot.samples,
                                      Test.data = u.test,
                                      Test.surv.data = surv.test)
  } else{
    res.clin.biom10 <- coxnet.optim.r(Train.data = u.train,
                                      pred.vec = c(clin.predictors,biom.opti10),
                                      Train.surv.data = surv.opt,
                                      times = 1000, 
                                      boot.samples = boot.samples,
                                      Test.data = u.test,
                                      Test.surv.data = surv.test)
    
  }
  save(res.clin.biom10, file = paste0("../data_output/biomarker_models_update/px_info_top10_biom_res_",dz,"_",inc.yrs,"yrs",".RData"))
  
}

## run px info model + 5 top biomarker based on weightsXselection.perc
if(nrow(opti.biom)>=5){
  opti.biom <- opti.biom %>%
    mutate(sel.opti.score=abs.coef*select.perc) %>%
    top_n(5, sel.opti.score)
  biom.opti <- opti.biom$biom
  } else {
    biom.opti <- opti.biom$biom
}


if(n.cases/4>200){
  res.clin.biom <- coxnet.optim.r(Train.data = u.opti,
                                 pred.vec = c(clin.predictors,biom.opti),
                                 Train.surv.data = surv.opt,
                                 times = 1000, 
                                 boot.samples = boot.samples,
                                 Test.data = u.test,
                                 Test.surv.data = surv.test)
} else {
  res.clin.biom <- coxnet.optim.r(Train.data = u.train,
                                 pred.vec = c(clin.predictors,biom.opti),
                                 Train.surv.data = surv.opt,
                                 times = 1000, 
                                 boot.samples = boot.samples,
                                 Test.data = u.test,
                                 Test.surv.data = surv.test)
  
}
save(res.clin.biom, file = paste0("../data_output/biomarker_models_update/px_info_top5_biom_res_",dz,"_",inc.yrs,"yrs",".RData"))

## optimization on top 5 biomarkers
if(n.cases/4>200){
  res.biom <- coxnet.optim.r(Train.data = u.opti,
                           pred.vec = biom.opti,
                           Train.surv.data = surv.opt,
                           times = 1000, 
                           boot.samples = boot.samples, 
                           Test.data = u.test,
                           Test.surv.data = surv.test)
} else{
  res.biom <- coxnet.optim.r(Train.data = u.train,
                           pred.vec = biom.opti,
                           Train.surv.data = surv.opt,
                           times = 1000, 
                           boot.samples = boot.samples,
                           Test.data = u.test,
                           Test.surv.data = surv.test)
  
}
## add feature selection to the results
res.biom$feature.selection <- p.select

## add case numbers for training and testing sets
if((n.cases/4)>200){
  case.num <- data.frame(train.cases = sum(u.train[,dz]==1),
                         opti.cases = sum(u.opti[,dz]==1),
                         test.cases = sum(u.test[,dz]==1),
                         train.controls = sum(u.train[,dz]!=1),
                         opti.controls = sum(u.opti[,dz]!=1),
                         test.controls = sum(u.test[,dz]!=1))
} else{
  case.num <- data.frame(train.cases = sum(u.train[,dz]==1),
                         opti.cases = NA,
                         test.cases = sum(u.test[,dz]==1),
                         train.controls = sum(u.train[,dz]!=1),
                         opti.controls = NA,
                         test.controls = sum(u.test[,dz]!=1))
}
res.biom$case.num <- case.num
## save results 
save(res.biom, file = paste0("../data_output/biomarker_models_update/biomarker_top5_res_",dz,"_",inc.yrs,"yrs",".RData"))

###################################################
#### Script Name: Detection_rate_curves_biomarkers
#### Author: Julia Carrasco Zanini Sanchez
#### Description: script to ML framework for incident diseases in atlas 308 
#### Date : 23/03/2023

rm(list = ls())

## load packages 
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(RNOmni)
# library(glue)
library(DBI)
library(dbplyr)
library(magrittr)
# library(broom)
library(survival)
library(data.table)
library(glue)
library(parallel)
library(caret)
library(glmnet)
library(ROSE)
# library(odbc)

## define disease and restrict to incident cases within X years
args <- commandArgs(trailingOnly = T)

dz <- as.character(args[1])
dz.date <- as.character(args[2])
inc.yrs <- as.numeric(args[3])

paste0("Prediction:", dz, "; ", dz.date, " at ", inc.yrs, "years")

## set up working directory 
setwd("/home/ukbb_inc_dz_prediction/Explore_1536_Expansion_analyses/02_inc_dz_prediction/bin/") 

## read in data 
atlas308 <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_atlas308_annotations.txt", data.table = F)
cov <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_covariate_annotations.txt", data.table = F)

####----  Data Cleaning and wrangling   ----####
atlas308 <- as.data.frame(atlas308 %>% 
                            filter(eid_20361%in%cov$eid_20361) %>% 
                            pivot_wider(names_from = variable, values_from =value))

## phenotypes missing for patient-derived info model: bmi, smoking, alcohol, physical activity, family history 
pheno <- merge(cov[,-c(2:14,17)], atlas308, by="eid_20361")

## read in px info phenotypes
px.info <- read.delim("../../01_inc_dz_phenos/data_input/Patient_model_variables.txt",sep = "\t")
pheno <- merge(pheno, px.info, by="eid_20361")

## rename to olds names to match the rest fo the script
colnames(pheno)[which(colnames(pheno)=="assess_date")] <- c("date_baseline_assessment")

## process ethnicity
pheno$ethnicity <- ifelse(is.na(pheno$ge_eur),2,pheno$ge_eur)

## get protein data and metadata
ol.npx <- fread("../../00_protein_imputation/data_input/Imputed_1536_Expansion_NPX_proteins.txt", data.table = F)
ol.npx <- filter(ol.npx, exclusion.50==0 )

## exclude participants with missing age, sex
ol.npx <- ol.npx[complete.cases(ol.npx[,c("ge_sex", "ge_age")]),]

protein_list <- colnames(ol.npx)[-which(colnames(ol.npx)%in%c("eid_20361","ge_sex","ge_age",
                                                              "randomized_baseline","exclusion.50"))]

## merge phenotypes and proteins
ukbb <- merge(ol.npx, pheno, by="eid_20361")
## randomised baseline
ukbb <-
  ukbb %>% 
  filter(randomized_baseline==T|!!as.name(dz)==1)
## inverse rank normal transform npx
ukbb[,protein_list] <- sapply(ukbb[,protein_list], RankNorm)

## read in biomarkers
biom <- read.delim("../../01_inc_dz_phenos/data_input/biomarkers_clean_230623_ukb.txt", sep = "\ ")
biom <- biom %>% 
  select(-ge_sex, -age_0)

biomarkers <- colnames(biom)[-ncol(biom)]

## exclude prevalent cases and incident cases up to 6 months
## restrict to incident cases within the first X years
ukbb <- ukbb %>% 
  mutate(fol = difftime(!!as.name(dz.date),date_baseline_assessment, units = "weeks")/52.25) %>% 
  filter((fol>0.5&fol<inc.yrs)|is.na(!!as.name(dz.date)))

## censor (including censoring by death)
death.data <- read.delim("../../01_inc_dz_phenos/data_input/Death_dates_proteomic_participants.txt",sep = "\t")
death.data$date_of_death <- as.Date(death.data$date_of_death, format = "%Y-%m-%d")

ukbb <- merge(ukbb, death.data[,c("eid_20361","date_of_death")], by="eid_20361", all.x=T)

if(inc.yrs==10){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2020-12-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
} else if(inc.yrs==5){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2016-05-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
}else{stop("limit to incident casese within X years not defined")}

if(inc.yrs==10){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2020-12-31"),
                       difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                       ukbb$fol)
} else if(inc.yrs==5){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2016-05-31"),
                       difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                       ukbb$fol)
}else{stop("limit to incident casese within X years not defined")}

ukbb$fol <- ifelse(is.na(ukbb$fol_d), ukbb$fol, ukbb$fol_d)

## for sex specific outcome exclude males/females
if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_pri_breast",
           "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical",
           "atlas_inc_endometrial_hyper","atlas_inc_endometriosis","atlas_inc_female_genital_prolapse",
           "atlas_inc_leiomyoma","atlas_inc_menorrhagia",
           "atlas_inc_pid","atlas_inc_pmb")){
  # ukbb <- ukbb %>% filter(ge_sex=="female")
  ukbb <- ukbb %>% filter(ge_sex==0)
} else if(dz%in%c("atlas_inc_pri_prost","atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  # ukbb <- ukbb %>% filter(ge_sex=="male")
  ukbb <- ukbb %>% filter(ge_sex==1)
} else{ukbb <- ukbb}

n.cases <- sum(ukbb[,dz]==1)

## add noise variables
rand.vars <- lapply(1:10, function(x){
  tmp <- runif(nrow(ukbb), min = min(ukbb$NPPB_20049), 
               max = max(ukbb$NPPB_20049))
  return(tmp)
})
rand.vars <- as.data.frame(do.call(cbind, rand.vars))
colnames(rand.vars) <- paste0("rand_var_",1:10)

ukbb <- cbind(ukbb, rand.vars)

####----    ML framework    ----####

## split into train/optimize/test or train/test depending on case number
ukbb[,dz] <- as.factor(ukbb[,dz])

if(n.cases/4>200){
  ## feature selection set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.5, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.opti <- ukbb[-trainIndex,]
  
  ## optimization and testing
  set.seed(3456) ## for reproducibility
  trainIndex2 <- createDataPartition(u.opti[,dz], p=0.5, list = F, times = 1)
  u.test <- u.opti[trainIndex2,]
  u.opti <- u.opti[-trainIndex2,]
  ## merge biomarkers  (note have preserved the same split as for other models)
  u.train <- merge(u.train, biom, by="eid_20361")
  u.opti <- merge(u.opti, biom, by="eid_20361")
  u.test <- merge(u.test, biom, by="eid_20361")
}else {
  ## feature selection set and test set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.7, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.test <- ukbb[-trainIndex,]
  u.train <- merge(u.train, biom, by="eid_20361")
  u.test <- merge(u.test, biom, by="eid_20361")
}
  
####---- optimization and performance testing ----####

## survival objects
surv.train <- Surv(u.train$fol, u.train[,dz])
## optimization 
if(n.cases/4>200){
  surv.opt <- Surv(u.opti$fol, u.opti[,dz])
}else{surv.opt <- Surv(u.train$fol, u.train[,dz])}
## testing
surv.test <- Surv(u.test$fol, u.test[,dz])

## px information model basic 
clin.predictors <- c("ge_age","ge_sex","bmi","ethnicity","smoke_never",
                     "smoke_previous","smoke_current","alcohol_daily",
                     "alcohol_4_week","alcohol_2_week","alcohol_3_month","alcohol_occassion",
                     "alcohol_never")

## define additional family history for specific diseases
if(!dz%in%c("atlas_inc_diabetes_t2","atlas_inc_pri_prost", "atlas_inc_copd","atlas_inc_isch_stroke",
            "atlas_inc_af","atlas_inc_av_block_3","atlas_inc_chd_nos","atlas_inc_hf", 
            "atlas_inc_peripheral_arterial_disease","atlas_inc_parkinsons","atlas_inc_pri_lung",
            "atlas_inc_pri_bowel", "atlas_inc_pri_breast","atlas_inc_pri_ovarian","atlas_inc_pri_uterine",
            "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
            "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
            "atlas_inc_hydrocele","atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb",
            "atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu")){
  clin.predictors <- clin.predictors
} else if(dz=="atlas_inc_diabetes_t2") {
  clin.predictors <- c(clin.predictors,"father_diabetes","mother_diabetes")
} else if(dz=="atlas_inc_pri_prost") {
  clin.predictors <- c(clin.predictors[-2],"father_prostate_cancer")
} else if(dz=="atlas_inc_copd") {
  clin.predictors <- c(clin.predictors,"father_COPD","mother_COPD")
} else if(dz%in%c("atlas_inc_isch_stroke","atlas_inc_af","atlas_inc_av_block_3",
                  "atlas_inc_chd_nos","atlas_inc_hf","atlas_inc_peripheral_arterial_disease")) {
  clin.predictors <- c(clin.predictors,"father_heart_disease","mother_heart_disease")
} else if(dz=="atlas_inc_parkinsons") {
  clin.predictors <- c(clin.predictors,"father_parkinsons","mother_parkinsons")
} else if(dz=="atlas_inc_pri_lung") {
  clin.predictors <- c(clin.predictors,"father_lung_cancer","mother_lung_cancer")
} else if(dz=="atlas_inc_pri_bowel") {
  clin.predictors <- c(clin.predictors,"father_bowel_cancer","mother_bowel_cancer")
} else if(dz=="atlas_inc_pri_breast") {
  clin.predictors <- c(clin.predictors[-2],"mother_breast_cancer")
} else if(dz=="atlas_inc_hypertension") {
  clin.predictors <- c(clin.predictors[-2],"mother_high_bp","father_high_bp")
} else if(dz=="atlas_inc_depression") {
  clin.predictors <- c(clin.predictors[-2],"mother_depression","father_depression")
} else if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_benign_ovary",
                  "atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
                  "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
                  "atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb")) {
  clin.predictors <- c(clin.predictors[-2])
} else if(dz%in%c("atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  clin.predictors <- c(clin.predictors[-2])
}


## load results for best performing biomarker model 
biom.res <- read.delim("../data_output/summary_updated_biomarker_results_10_year_models.txt", sep = "\t")
biom.res <- biom.res[which(biom.res$dz==dz),]

if (biom.res$model=="clin.biom.all"){
  load(paste0("../data_output/biomarker_models_update/px_info_all_biom_above_noise_res_",dz,"_",inc.yrs,"yrs",".RData"))
  res.clin.biom <- res.clin.biom.all
} else if (biom.res$model=="clin.biom10") {
  load(paste0("../data_output/biomarker_models_update/px_info_top10_biom_res_",dz,"_",inc.yrs,"yrs",".RData"))
  res.clin.biom <- res.clin.biom10
} else if (biom.res$model=="clin.biom") {
  load(paste0("../data_output/biomarker_models_update/px_info_top5_biom_res_",dz,"_",inc.yrs,"yrs",".RData"))
  res.clin.biom <- res.clin.biom
} else if (biom.res$model=="biom") {
  load(paste0("../data_output/biomarker_models_update/biomarker_top5_res_",dz,"_",inc.yrs,"yrs",".RData"))
  res.clin.biom <- res.biom
}

####---- DR using test set only ----####
## normalise linear predictor (0-1 values)

if(min(res.clin.biom$linear.predictor)<0){
  lp.biom <- (res.clin.biom$linear.predictor+abs(min(res.clin.biom$linear.predictor)))/
    max((res.clin.biom$linear.predictor+abs(min(res.clin.biom$linear.predictor))))
} else{
  lp.biom <- (res.clin.biom$linear.predictor)/max(res.clin.biom$linear.predictor)
}

## exclude consortium selected cases
# consort <- which(u.test$randomized_baseline==F)
# u.test <- u.test[-consort,]
# lp.biom <- lp.biom[-consort]

## find threshold for detection rate at 5% false positive rate - test set only 
thresh.grid <- lapply(seq(0.0001,1,0.0001), function(x){
  
  if(sum(lp.biom<x)>0&sum(lp.biom>x)>0){
    ## biomarkers - define predicted case status based on the selected threshold
    pred.case.biom <- ifelse(lp.biom>x,1,0)

    tp.biom <- sum(u.test[,dz]==1&pred.case.biom==1)
    fp.biom <- sum(u.test[,dz]==0&pred.case.biom==1)
    fn.biom <- sum(u.test[,dz]==1&pred.case.biom==0)
    tn.biom <- sum(u.test[,dz]==0&pred.case.biom==0)
    ## fpr = fp / tn + fp
    fpr.biom <- fp.biom / (tn.biom + fp.biom)
    ## dr = tp / tp+fn
    dr.biom <- tp.biom / (tp.biom+fn.biom)
    tmp.biom <- data.frame("threshold"=x,
                            "biom.fpr"=fpr.biom*100,
                            "biom.DR"=dr.biom)
    } else{
      tmp.biom <- data.frame("threshold"=x,
                             "biom.fpr"=NA,
                             "biom.DR"=NA)
    }
  
  ## return object 
  return(tmp.biom)
  })

thresh.grid <- do.call(rbind, thresh.grid)

## subset to the closest thresholds to 5 - 50 % fpr
thresh.biom <- lapply(c(5,10,15,20,25,30,35,40,45,50), function(x){
  tmp.biom <- thresh.grid[which.min(abs(x-thresh.grid$biom.fpr)),]
  return(tmp.biom)
})
thresh.biom <- as.data.frame(do.call(rbind, thresh.biom))
thresh.biom$FPR.true <- c(5,10,15,20,25,30,35,40,45,50)

res.dr.curve.test <- thresh.biom
colnames(res.dr.curve.test)[1] <- c("threshold.biom")
res.dr.curve.test$set <- rep("test.set", nrow(res.dr.curve.test))
  
####---- DR in the entire UKB subcohort ----####
ukbb <- merge(ukbb, biom, by="eid_20361")
# consort <- which(ukbb$randomized_baseline==F)
# ukbb <- ukbb[-consort,]

lp.biom <- as.numeric(predict(res.clin.biom$glmnet.opt,
                               type = "link",
                               newx = as.matrix(ukbb[, row.names(res.clin.biom$opt.coefficients)]),
                               s = res.clin.biom$lambda.opt))

## normalise 0- 1
if(min(lp.biom)<0){
  lp.biom <- (lp.biom+abs(min(lp.biom)))/
    max((lp.biom+abs(min(lp.biom))))
} else{
  lp.biom <- (lp.biom)/max(lp.biom)
}


## find threshold for detection rate at 5% false positive rate - test set only 
thresh.grid <- lapply(seq(0.0001,1,0.0001), function(x){
  if(sum(lp.biom<x)>0&sum(lp.biom>x)>0){
    ## biomarkers - define predicted case status based on the selected threshold
    pred.case.biom <- ifelse(lp.biom>x,1,0)

    tp.biom <- sum(ukbb[,dz]==1&pred.case.biom==1)
    fp.biom <- sum(ukbb[,dz]==0&pred.case.biom==1)
    fn.biom <- sum(ukbb[,dz]==1&pred.case.biom==0)
    tn.biom <- sum(ukbb[,dz]==0&pred.case.biom==0)
    ## fpr = fp / tn + fp
    fpr.biom <- fp.biom / (tn.biom + fp.biom)
    ## dr = tp / tp+fn
    dr.biom <- tp.biom / (tp.biom+fn.biom)
    tmp.biom <- data.frame("threshold"=x,
                           "biom.fpr"=fpr.biom*100,
                           "biom.DR"=dr.biom)
  } else{
    tmp.biom <- data.frame("threshold"=x,
                           "biom.fpr"=NA,
                           "biom.DR"=NA)
  }
  ## return object 
  return(tmp.biom)
})


thresh.grid <- do.call(rbind, thresh.grid)

## subset to the closest thresholds to 5 - 50 % fpr

thresh.biom <- lapply(c(5,10,15,20,25,30,35,40,45,50), function(x){
  tmp.biom <- thresh.grid[which.min(abs(x-thresh.grid$biom.fpr)),]
  return(tmp.biom)
})
thresh.biom <- as.data.frame(do.call(rbind, thresh.biom))
thresh.biom$FPR.true <- c(5,10,15,20,25,30,35,40,45,50)

res.dr.curve.all <- thresh.biom
colnames(res.dr.curve.all)[1] <- c("threshold.biom")
res.dr.curve.all$set <- rep("all", nrow(res.dr.curve.all))
res.dr.curve <- rbind(res.dr.curve.all, res.dr.curve.test)
write.table(res.dr.curve, file=paste0("../data_output/detection_rate_curves/updated_DR_curves_biom_" ,dz, "_", inc.yrs,"years.txt"), sep = "\t", row.names = F)
###################################################
#### Script Name: Detection_rate_curves_biomarkers
#### Author: Julia Carrasco Zanini Sanchez
#### Description: script to ML framework for incident diseases in atlas 308 
#### Date : 23/03/2023

rm(list = ls())

## load packages 
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(RNOmni)
# library(glue)
library(DBI)
library(dbplyr)
library(magrittr)
# library(broom)
library(survival)
library(data.table)
library(glue)
library(parallel)
library(caret)
library(glmnet)
library(ROSE)
# library(odbc)

## define disease and restrict to incident cases within X years
args <- commandArgs(trailingOnly = T)

dz <- as.character(args[1])
dz.date <- as.character(args[2])
inc.yrs <- as.numeric(args[3])

paste0("Prediction:", dz, "; ", dz.date, " at ", inc.yrs, "years")

## set up working directory 
setwd("/home/ukbb_inc_dz_prediction/Explore_1536_Expansion_analyses/02_inc_dz_prediction/bin/") 

## read in data 
atlas308 <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_atlas308_annotations.txt", data.table = F)
cov <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_covariate_annotations.txt", data.table = F)

####----  Data Cleaning and wrangling   ----####
atlas308 <- as.data.frame(atlas308 %>% 
                            filter(eid_20361%in%cov$eid_20361) %>% 
                            pivot_wider(names_from = variable, values_from =value))

## phenotypes missing for patient-derived info model: bmi, smoking, alcohol, physical activity, family history 
pheno <- merge(cov[,-c(2:14,17)], atlas308, by="eid_20361")

## read in px info phenotypes
px.info <- read.delim("../../01_inc_dz_phenos/data_input/Patient_model_variables.txt",sep = "\t")
pheno <- merge(pheno, px.info, by="eid_20361")

## rename to olds names to match the rest fo the script
colnames(pheno)[which(colnames(pheno)=="assess_date")] <- c("date_baseline_assessment")

## process ethnicity
pheno$ethnicity <- ifelse(is.na(pheno$ge_eur),2,pheno$ge_eur)

## get protein data and metadata
ol.npx <- fread("../../00_protein_imputation/data_input/Imputed_1536_Expansion_NPX_proteins.txt", data.table = F)
ol.npx <- filter(ol.npx, exclusion.50==0 )

## exclude participants with missing age, sex
ol.npx <- ol.npx[complete.cases(ol.npx[,c("ge_sex", "ge_age")]),]

protein_list <- colnames(ol.npx)[-which(colnames(ol.npx)%in%c("eid_20361","ge_sex","ge_age",
                                                              "randomized_baseline","exclusion.50"))]

## merge phenotypes and proteins
ukbb <- merge(ol.npx, pheno, by="eid_20361")
## randomised baseline
ukbb <- ukbb %>% 
  filter(randomized_baseline==T)
## inverse rank normal transform npx
ukbb[,protein_list] <- sapply(ukbb[,protein_list], RankNorm)

## read in biomarkers
biom <- read.delim("../../01_inc_dz_phenos/data_input/biomarkers_clean_230623_ukb.txt", sep = "\ ")
biom <- biom %>% 
  select(-ge_sex, -age_0)

biomarkers <- colnames(biom)[-ncol(biom)]

## exclude prevalent cases and incident cases up to 6 months
## restrict to incident cases within the first X years
ukbb <- ukbb %>% 
  mutate(fol = difftime(!!as.name(dz.date),date_baseline_assessment, units = "weeks")/52.25) %>% 
  filter((fol>0.5&fol<inc.yrs)|is.na(!!as.name(dz.date)))

## censor (including censoring by death)
death.data <- read.delim("../../01_inc_dz_phenos/data_input/Death_dates_proteomic_participants.txt",sep = "\t")
death.data$date_of_death <- as.Date(death.data$date_of_death, format = "%Y-%m-%d")

ukbb <- merge(ukbb, death.data[,c("eid_20361","date_of_death")], by="eid_20361", all.x=T)

if(inc.yrs==10){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2020-12-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
} else if(inc.yrs==5){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2016-05-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
}else{stop("limit to incident casese within X years not defined")}

if(inc.yrs==10){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2020-12-31"),
                       difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                       ukbb$fol)
} else if(inc.yrs==5){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2016-05-31"),
                       difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                       ukbb$fol)
}else{stop("limit to incident casese within X years not defined")}

ukbb$fol <- ifelse(is.na(ukbb$fol_d), ukbb$fol, ukbb$fol_d)

## for sex specific outcome exclude males/females
if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_pri_breast",
           "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical",
           "atlas_inc_endometrial_hyper","atlas_inc_endometriosis","atlas_inc_female_genital_prolapse",
           "atlas_inc_leiomyoma","atlas_inc_menorrhagia",
           "atlas_inc_pid","atlas_inc_pmb")){
  # ukbb <- ukbb %>% filter(ge_sex=="female")
  ukbb <- ukbb %>% filter(ge_sex==0)
} else if(dz%in%c("atlas_inc_pri_prost","atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  # ukbb <- ukbb %>% filter(ge_sex=="male")
  ukbb <- ukbb %>% filter(ge_sex==1)
} else{ukbb <- ukbb}

n.cases <- sum(ukbb[,dz]==1)

## add noise variables
rand.vars <- lapply(1:10, function(x){
  tmp <- runif(nrow(ukbb), min = min(ukbb$NPPB_20049), 
               max = max(ukbb$NPPB_20049))
  return(tmp)
})
rand.vars <- as.data.frame(do.call(cbind, rand.vars))
colnames(rand.vars) <- paste0("rand_var_",1:10)

ukbb <- cbind(ukbb, rand.vars)

####----    ML framework    ----####

## split into train/optimize/test or train/test depending on case number
ukbb[,dz] <- as.factor(ukbb[,dz])

if(n.cases/4>200){
  ## feature selection set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.5, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.opti <- ukbb[-trainIndex,]
  
  ## optimization and testing
  set.seed(3456) ## for reproducibility
  trainIndex2 <- createDataPartition(u.opti[,dz], p=0.5, list = F, times = 1)
  u.test <- u.opti[trainIndex2,]
  u.opti <- u.opti[-trainIndex2,]
  ## merge biomarkers  (note have preserved the same split as for other models)
  u.train <- merge(u.train, biom, by="eid_20361")
  u.opti <- merge(u.opti, biom, by="eid_20361")
  u.test <- merge(u.test, biom, by="eid_20361")
}else {
  ## feature selection set and test set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.7, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.test <- ukbb[-trainIndex,]
  u.train <- merge(u.train, biom, by="eid_20361")
  u.test <- merge(u.test, biom, by="eid_20361")
}
  
####---- optimization and performance testing ----####

## survival objects
surv.train <- Surv(u.train$fol, u.train[,dz])
## optimization 
if(n.cases/4>200){
  surv.opt <- Surv(u.opti$fol, u.opti[,dz])
}else{surv.opt <- Surv(u.train$fol, u.train[,dz])}
## testing
surv.test <- Surv(u.test$fol, u.test[,dz])

## px information model basic 
clin.predictors <- c("ge_age","ge_sex","bmi","ethnicity","smoke_never",
                     "smoke_previous","smoke_current","alcohol_daily",
                     "alcohol_4_week","alcohol_2_week","alcohol_3_month","alcohol_occassion",
                     "alcohol_never")

## define additional family history for specific diseases
if(!dz%in%c("atlas_inc_diabetes_t2","atlas_inc_pri_prost", "atlas_inc_copd","atlas_inc_isch_stroke",
            "atlas_inc_af","atlas_inc_av_block_3","atlas_inc_chd_nos","atlas_inc_hf", 
            "atlas_inc_peripheral_arterial_disease","atlas_inc_parkinsons","atlas_inc_pri_lung",
            "atlas_inc_pri_bowel", "atlas_inc_pri_breast","atlas_inc_pri_ovarian","atlas_inc_pri_uterine",
            "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
            "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
            "atlas_inc_hydrocele","atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb",
            "atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu")){
  clin.predictors <- clin.predictors
} else if(dz=="atlas_inc_diabetes_t2") {
  clin.predictors <- c(clin.predictors,"father_diabetes","mother_diabetes")
} else if(dz=="atlas_inc_pri_prost") {
  clin.predictors <- c(clin.predictors[-2],"father_prostate_cancer")
} else if(dz=="atlas_inc_copd") {
  clin.predictors <- c(clin.predictors,"father_COPD","mother_COPD")
} else if(dz%in%c("atlas_inc_isch_stroke","atlas_inc_af","atlas_inc_av_block_3",
                  "atlas_inc_chd_nos","atlas_inc_hf","atlas_inc_peripheral_arterial_disease")) {
  clin.predictors <- c(clin.predictors,"father_heart_disease","mother_heart_disease")
} else if(dz=="atlas_inc_parkinsons") {
  clin.predictors <- c(clin.predictors,"father_parkinsons","mother_parkinsons")
} else if(dz=="atlas_inc_pri_lung") {
  clin.predictors <- c(clin.predictors,"father_lung_cancer","mother_lung_cancer")
} else if(dz=="atlas_inc_pri_bowel") {
  clin.predictors <- c(clin.predictors,"father_bowel_cancer","mother_bowel_cancer")
} else if(dz=="atlas_inc_pri_breast") {
  clin.predictors <- c(clin.predictors[-2],"mother_breast_cancer")
} else if(dz=="atlas_inc_hypertension") {
  clin.predictors <- c(clin.predictors[-2],"mother_high_bp","father_high_bp")
} else if(dz=="atlas_inc_depression") {
  clin.predictors <- c(clin.predictors[-2],"mother_depression","father_depression")
} else if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_benign_ovary",
                  "atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
                  "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
                  "atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb")) {
  clin.predictors <- c(clin.predictors[-2])
} else if(dz%in%c("atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  clin.predictors <- c(clin.predictors[-2])
}


## load results for best performing biomarker model 
biom.res <- read.delim("../data_output/summary_updated_biomarker_results_10_year_models.txt", sep = "\t")
biom.res <- biom.res[which(biom.res$dz==dz),]

if (biom.res$model=="clin.biom.all"){
  load(paste0("../data_output/biomarker_models_update/px_info_all_biom_above_noise_res_",dz,"_",inc.yrs,"yrs",".RData"))
  res.clin.biom <- res.clin.biom.all
} else if (biom.res$model=="clin.biom10") {
  load(paste0("../data_output/biomarker_models_update/px_info_top10_biom_res_",dz,"_",inc.yrs,"yrs",".RData"))
  res.clin.biom <- res.clin.biom10
} else if (biom.res$model=="clin.biom") {
  load(paste0("../data_output/biomarker_models_update/px_info_top5_biom_res_",dz,"_",inc.yrs,"yrs",".RData"))
  res.clin.biom <- res.clin.biom
} else if (biom.res$model=="biom") {
  load(paste0("../data_output/biomarker_models_update/biomarker_top5_res_",dz,"_",inc.yrs,"yrs",".RData"))
  res.clin.biom <- res.biom
}

####---- DR using test set only ----####
## normalise linear predictor (0-1 values)

if(min(res.clin.biom$linear.predictor)<0){
  lp.biom <- (res.clin.biom$linear.predictor+abs(min(res.clin.biom$linear.predictor)))/
    max((res.clin.biom$linear.predictor+abs(min(res.clin.biom$linear.predictor))))
} else{
  lp.biom <- (res.clin.biom$linear.predictor)/max(res.clin.biom$linear.predictor)
}

## find threshold for detection rate at 5% false positive rate - test set only 
thresh.grid <- lapply(seq(0.0001,1,0.0001), function(x){
  
  if(sum(lp.biom<x)>0&sum(lp.biom>x)>0){
    ## biomarkers - define predicted case status based on the selected threshold
    pred.case.biom <- ifelse(lp.biom>x,1,0)

    tp.biom <- sum(u.test[,dz]==1&pred.case.biom==1)
    fp.biom <- sum(u.test[,dz]==0&pred.case.biom==1)
    fn.biom <- sum(u.test[,dz]==1&pred.case.biom==0)
    tn.biom <- sum(u.test[,dz]==0&pred.case.biom==0)
    ## fpr = fp / tn + fp
    fpr.biom <- fp.biom / (tn.biom + fp.biom)
    ## dr = tp / tp+fn
    dr.biom <- tp.biom / (tp.biom+fn.biom)
    tmp.biom <- data.frame("threshold"=x,
                            "biom.fpr"=fpr.biom*100,
                            "biom.DR"=dr.biom)
    } else{
      tmp.biom <- data.frame("threshold"=x,
                             "biom.fpr"=NA,
                             "biom.DR"=NA)
    }
  
  ## return object 
  return(tmp.biom)
  })

thresh.grid <- do.call(rbind, thresh.grid)

## subset to the closest thresholds to 5 - 50 % fpr
thresh.biom <- lapply(c(5,10,15,20,25,30,35,40,45,50), function(x){
  tmp.biom <- thresh.grid[which.min(abs(x-thresh.grid$biom.fpr)),]
  return(tmp.biom)
})
thresh.biom <- as.data.frame(do.call(rbind, thresh.biom))
thresh.biom$FPR.true <- c(5,10,15,20,25,30,35,40,45,50)

res.dr.curve.test <- thresh.biom
colnames(res.dr.curve.test)[1] <- c("threshold.biom")
res.dr.curve.test$set <- rep("test.set", nrow(res.dr.curve.test))
  
####---- DR in the entire UKB subcohort ----####
ukbb <- merge(ukbb, biom, by="eid_20361")
lp.biom <- as.numeric(predict(res.clin.biom$glmnet.opt,
                               type = "link",
                               newx = as.matrix(ukbb[, row.names(res.clin.biom$opt.coefficients)]),
                               s = res.clin.biom$lambda.opt))

## normalise 0- 1
if(min(lp.biom)<0){
  lp.biom <- (lp.biom+abs(min(lp.biom)))/
    max((lp.biom+abs(min(lp.biom))))
} else{
  lp.biom <- (lp.biom)/max(lp.biom)
}


## find threshold for detection rate at 5% false positive rate - test set only 
thresh.grid <- lapply(seq(0.0001,1,0.0001), function(x){
  if(sum(lp.biom<x)>0&sum(lp.biom>x)>0){
    ## biomarkers - define predicted case status based on the selected threshold
    pred.case.biom <- ifelse(lp.biom>x,1,0)

    tp.biom <- sum(ukbb[,dz]==1&pred.case.biom==1)
    fp.biom <- sum(ukbb[,dz]==0&pred.case.biom==1)
    fn.biom <- sum(ukbb[,dz]==1&pred.case.biom==0)
    tn.biom <- sum(ukbb[,dz]==0&pred.case.biom==0)
    ## fpr = fp / tn + fp
    fpr.biom <- fp.biom / (tn.biom + fp.biom)
    ## dr = tp / tp+fn
    dr.biom <- tp.biom / (tp.biom+fn.biom)
    tmp.biom <- data.frame("threshold"=x,
                           "biom.fpr"=fpr.biom*100,
                           "biom.DR"=dr.biom)
  } else{
    tmp.biom <- data.frame("threshold"=x,
                           "biom.fpr"=NA,
                           "biom.DR"=NA)
  }
  ## return object 
  return(tmp.biom)
})


thresh.grid <- do.call(rbind, thresh.grid)

## subset to the closest thresholds to 5 - 50 % fpr

thresh.biom <- lapply(c(5,10,15,20,25,30,35,40,45,50), function(x){
  tmp.biom <- thresh.grid[which.min(abs(x-thresh.grid$biom.fpr)),]
  return(tmp.biom)
})
thresh.biom <- as.data.frame(do.call(rbind, thresh.biom))
thresh.biom$FPR.true <- c(5,10,15,20,25,30,35,40,45,50)

res.dr.curve.all <- thresh.biom
colnames(res.dr.curve.all)[1] <- c("threshold.biom")
res.dr.curve.all$set <- rep("all", nrow(res.dr.curve.all))
res.dr.curve <- rbind(res.dr.curve.all, res.dr.curve.test)
write.table(res.dr.curve, file=paste0("../data_output/detection_rate_curves/updated_DR_curves_biom_" ,dz, "_", inc.yrs,"years.txt"), sep = "\t", row.names = F)
###################################################
#### Script Name: Detection_rate_curves
#### Author: Julia Carrasco Zanini Sanchez
#### Description: script to ML framework for incident diseases in atlas 308 
#### Date : 23/03/2023

rm(list = ls())

## load packages 
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(RNOmni)
# library(glue)
library(DBI)
library(dbplyr)
library(magrittr)
# library(broom)
library(survival)
library(data.table)
library(glue)
library(parallel)
library(caret)
library(glmnet)
library(ROSE)
# library(odbc)

## define disease and restrict to incident cases within X years
args <- commandArgs(trailingOnly = T)

dz <- as.character(args[1])
dz.date <- as.character(args[2])
inc.yrs <- as.numeric(args[3])

paste0("Prediction:", dz, "; ", dz.date, " at ", inc.yrs, "years")

## set up working directory 
setwd("/home/ukbb_inc_dz_prediction/Explore_1536_Expansion_analyses/02_inc_dz_prediction/bin/") 

## read in data 
atlas308 <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_atlas308_annotations.txt", data.table = F)
cov <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_covariate_annotations.txt", data.table = F)

####----  Data Cleaning and wrangling   ----####
atlas308 <- as.data.frame(atlas308 %>% 
  filter(eid_20361%in%cov$eid_20361) %>% 
  pivot_wider(names_from = variable, values_from =value))

## phenotypes missing for patient-derived info model: bmi, smoking, alcohol, physical activity, family history 
pheno <- merge(cov[,-c(11:14,17)], atlas308, by="eid_20361")

## read in px info phenotypes
px.info <- read.delim("../../01_inc_dz_phenos/data_input/Patient_model_variables.txt",sep = "\t")
pheno <- merge(pheno, px.info, by="eid_20361")

## rename to olds names to match the rest fo the script
colnames(pheno)[which(colnames(pheno)=="assess_date")] <- c("date_baseline_assessment")

## process ethnicity
pheno$ethnicity <- ifelse(is.na(pheno$ge_eur),2,pheno$ge_eur)

## get protein data and metadata
ol.npx <- fread("../../00_protein_imputation/data_input/Imputed_1536_Expansion_NPX_proteins.txt", data.table = F)
ol.npx <- filter(ol.npx, exclusion.50==0 )

## exclude participants with missing age, sex
ol.npx <- ol.npx[complete.cases(ol.npx[,c("ge_sex", "ge_age")]),]

protein_list <- colnames(ol.npx)[-which(colnames(ol.npx)%in%c("eid_20361","ge_sex","ge_age",
                                                              "randomized_baseline","exclusion.50"))]

## merge phenotypes and proteins
ukbb <- merge(ol.npx, pheno, by="eid_20361")
## randomised baseline
ukbb <-
  ukbb %>% 
  filter(randomized_baseline==T|!!as.name(dz)==1)
# inverse rank normal transform npx
ukbb[,protein_list] <- sapply(ukbb[,protein_list], RankNorm)

## exclude prevalent cases and incident cases up to 6 months
## restrict to incident cases within the first X years
ukbb <- ukbb %>% 
  mutate(fol = difftime(!!as.name(dz.date),date_baseline_assessment, units = "weeks")/52.25) %>% 
  filter((fol>0.5&fol<inc.yrs)|is.na(!!as.name(dz.date)))

## censor (including censoring by death)
death.data <- read.delim("../../01_inc_dz_phenos/data_input/Death_dates_proteomic_participants.txt",sep = "\t")
death.data$date_of_death <- as.Date(death.data$date_of_death, format = "%Y-%m-%d")

ukbb <- merge(ukbb, death.data[,c("eid_20361","date_of_death")], by="eid_20361", all.x=T)

if(inc.yrs==10){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2020-12-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
} else if(inc.yrs==5){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2016-05-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
}else{stop("limit to incident casese within X years not defined")}

if(inc.yrs==10){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2020-12-31"),
                         difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                        ukbb$fol)
} else if(inc.yrs==5){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2016-05-31"),
                          difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                          ukbb$fol)
}else{stop("limit to incident casese within X years not defined")}

ukbb$fol <- ifelse(is.na(ukbb$fol_d), ukbb$fol, ukbb$fol_d)

## for sex specific outcome exclude males/females
if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_pri_breast",
           "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical",
           "atlas_inc_endometrial_hyper","atlas_inc_endometriosis","atlas_inc_female_genital_prolapse",
           "atlas_inc_leiomyoma","atlas_inc_menorrhagia",
           "atlas_inc_pid","atlas_inc_pmb")){
  # ukbb <- ukbb %>% filter(ge_sex=="female")
  ukbb <- ukbb %>% filter(ge_sex==0)
} else if(dz%in%c("atlas_inc_pri_prost","atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  # ukbb <- ukbb %>% filter(ge_sex=="male")
  ukbb <- ukbb %>% filter(ge_sex==1)
} else{ukbb <- ukbb}

## recode sex to numeric
# ukbb$ge_sex <- ifelse(ukbb$ge_sex=="male",1,
#                       ifelse(ukbb$ge_sex=="female",0,NA))

n.cases <- sum(ukbb[,dz]==1)

## add noise variables
rand.vars <- lapply(1:10, function(x){
  tmp <- runif(nrow(ukbb), min = min(ukbb$NPPB_20049), 
               max = max(ukbb$NPPB_20049))
  return(tmp)
})
rand.vars <- as.data.frame(do.call(cbind, rand.vars))
colnames(rand.vars) <- paste0("rand_var_",1:10)

ukbb <- cbind(ukbb, rand.vars)

####----    ML framework    ----####

## split into train/optimize/test or train/test depending on case number
ukbb[,dz] <- as.factor(ukbb[,dz])

if(n.cases/4>200){
  ## feature selection set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.5, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.opti <- ukbb[-trainIndex,]
  
  ## optimization and testing
  set.seed(3456) ## for reproducibility
  trainIndex2 <- createDataPartition(u.opti[,dz], p=0.5, list = F, times = 1)
  u.test <- u.opti[trainIndex2,]
  u.opti <- u.opti[-trainIndex2,]
  }else {
  ## feature selection set and test set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.7, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.test <- ukbb[-trainIndex,]
  }
  
####---- optimization and performance testing ----####

## survival objects
surv.train <- Surv(u.train$fol, u.train[,dz])
## optimization 
if(n.cases/4>200){
  surv.opt <- Surv(u.opti$fol, u.opti[,dz])
}else{surv.opt <- Surv(u.train$fol, u.train[,dz])}
## testing
surv.test <- Surv(u.test$fol, u.test[,dz])

## px information model basic 
clin.predictors <- c("ge_age","ge_sex","bmi","ethnicity","smoke_never",
                     "smoke_previous","smoke_current","alcohol_daily",
                     "alcohol_4_week","alcohol_2_week","alcohol_3_month","alcohol_occassion",
                     "alcohol_never")

## define additional family history for specific diseases
if(!dz%in%c("atlas_inc_diabetes_t2","atlas_inc_pri_prost", "atlas_inc_copd","atlas_inc_isch_stroke",
            "atlas_inc_af","atlas_inc_av_block_3","atlas_inc_chd_nos","atlas_inc_hf", 
            "atlas_inc_peripheral_arterial_disease","atlas_inc_parkinsons","atlas_inc_pri_lung",
            "atlas_inc_pri_bowel", "atlas_inc_pri_breast","atlas_inc_pri_ovarian","atlas_inc_pri_uterine",
            "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
            "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
            "atlas_inc_hydrocele","atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb",
            "atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu")){
  clin.predictors <- clin.predictors
} else if(dz=="atlas_inc_diabetes_t2") {
  clin.predictors <- c(clin.predictors,"father_diabetes","mother_diabetes")
} else if(dz=="atlas_inc_pri_prost") {
  clin.predictors <- c(clin.predictors[-2],"father_prostate_cancer")
} else if(dz=="atlas_inc_copd") {
  clin.predictors <- c(clin.predictors,"father_COPD","mother_COPD")
} else if(dz%in%c("atlas_inc_isch_stroke","atlas_inc_af","atlas_inc_av_block_3",
                  "atlas_inc_chd_nos","atlas_inc_hf","atlas_inc_peripheral_arterial_disease")) {
  clin.predictors <- c(clin.predictors,"father_heart_disease","mother_heart_disease")
} else if(dz=="atlas_inc_parkinsons") {
  clin.predictors <- c(clin.predictors,"father_parkinsons","mother_parkinsons")
} else if(dz=="atlas_inc_pri_lung") {
  clin.predictors <- c(clin.predictors,"father_lung_cancer","mother_lung_cancer")
} else if(dz=="atlas_inc_pri_bowel") {
  clin.predictors <- c(clin.predictors,"father_bowel_cancer","mother_bowel_cancer")
} else if(dz=="atlas_inc_pri_breast") {
  clin.predictors <- c(clin.predictors[-2],"mother_breast_cancer")
} else if(dz=="atlas_inc_hypertension") {
  clin.predictors <- c(clin.predictors[-2],"mother_high_bp","father_high_bp")
} else if(dz=="atlas_inc_depression") {
  clin.predictors <- c(clin.predictors[-2],"mother_depression","father_depression")
} else if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_benign_ovary",
                  "atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
                  "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
                  "atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb")) {
  clin.predictors <- c(clin.predictors[-2])
} else if(dz%in%c("atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  clin.predictors <- c(clin.predictors[-2])
}

## load results for best performing protein model 
## load results for best performing protein model 
summary.res <- read.delim("../data_output/summary_results_10_year_models.txt", sep = "\t")
summary.res <- summary.res[which(summary.res$dz==dz),]

if (summary.res$model=="clin.prots.all"){
  load(paste0("../data_output/px_info_protein_models/px_info_all_prots_above_noise_res_",dz,"_",inc.yrs,"yrs",".RData"))
  res.clin.prots <- res.clin.prots.all
} else if (summary.res$model=="clin.prots10") {
  load(paste0("../data_output/px_info_protein_models/px_info_top10_prots_res_",dz,"_",inc.yrs,"yrs",".RData"))
  res.clin.prots <- res.clin.prots10
} else if (summary.res$model=="clin.prots") {
  load(paste0("../data_output/px_info_protein_models/px_info_top5_prots_res_",dz,"_",inc.yrs,"yrs",".RData"))
  res.clin.prots <- res.clin.prots
} else if (summary.res$model=="prot") {
  load(paste0("../data_output/protein_models/prot_top5_res_",dz,"_",inc.yrs,"yrs",".RData"))
  res.clin.prots <- res.prot
}

## load results for best performing biomarker model 
# biom.res <- read.delim("../data_output/summary_biomarker_results_10_year_models.txt", sep = "\t")
# biom.res <- biom.res[which(biom.res$dz==dz),]
# 
# if (biom.res$model=="clin.biom.all"){
#   load(paste0("../data_output/px_info_biomarker_models/px_info_all_biom_above_noise_res_",dz,"_",inc.yrs,"yrs",".RData"))
#   res.clin.biom <- res.clin.biom.all
# } else if (biom.res$model=="clin.biom10") {
#   load(paste0("../data_output/px_info_biomarker_models/px_info_top10_biom_res_",dz,"_",inc.yrs,"yrs",".RData"))
#   res.clin.biom <- res.clin.biom10
# } else if (biom.res$model=="clin.biom") {
#   load(paste0("../data_output/px_info_biomarker_models/px_info_top5_biom_res_",dz,"_",inc.yrs,"yrs",".RData"))
#   res.clin.biom <- res.clin.biom
# } else if (biom.res$model=="biom") {
#   load(paste0("../data_output/biomarker_models/biomarker_top5_res_",dz,"_",inc.yrs,"yrs",".RData"))
#   res.clin.biom <- res.biom
# }

## clinical model 
load(paste0("../data_output/px_info_models/px_info_res_",dz,"_",inc.yrs,"yrs",".RData"))

####---- DR using test set only ----####
## normalise linear predictor (0-1 values)
if(min(res.clin$linear.predictor)<0){
  lp.clin <- (res.clin$linear.predictor+abs(min(res.clin$linear.predictor)))/
    max((res.clin$linear.predictor+abs(min(res.clin$linear.predictor))))
} else{
  lp.clin <- (res.clin$linear.predictor)/max(res.clin$linear.predictor)
}

if(min(res.clin.prots$linear.predictor)<0){
  lp.prots <- (res.clin.prots$linear.predictor+abs(min(res.clin.prots$linear.predictor)))/
    max((res.clin.prots$linear.predictor+abs(min(res.clin.prots$linear.predictor))))
} else{
  lp.prots <- (res.clin.prots$linear.predictor)/max(res.clin.prots$linear.predictor)
}

## exclude consortium selected cases
# consort <- which(u.test$randomized_baseline==F)
# u.test <- u.test[-consort,]
# lp.clin <- lp.clin[-consort]
# lp.prots <- lp.prots[-consort]

# if(min(res.clin.biom$linear.predictor)<0){
#   lp.biom <- (res.clin.biom$linear.predictor+abs(min(res.clin.biom$linear.predictor)))/
#     max((res.clin.biom$linear.predictor+abs(min(res.clin.biom$linear.predictor))))
# } else{
#   lp.biom <- (res.clin.biom$linear.predictor)/max(res.clin.biom$linear.predictor)
# }

## find threshold for detection rate at 5% false positive rate - test set only 
thresh.grid <- lapply(seq(0.0001,1,0.0001), function(x){
  
  if(sum(lp.clin<x)>0&sum(lp.clin>x)>0){
    ## clinical - define predicted case status based on the selected threshold
    pred.case.clin <- ifelse(lp.clin>x,1,0)
    
    tp.clin <- sum(u.test[,dz]==1&pred.case.clin==1)
    fp.clin <- sum(u.test[,dz]==0&pred.case.clin==1)
    fn.clin <- sum(u.test[,dz]==1&pred.case.clin==0)
    tn.clin <- sum(u.test[,dz]==0&pred.case.clin==0)
    ## fpr = fp / tn + fp 
    fpr.clin <- fp.clin / (tn.clin + fp.clin)
    ## dr = tp / tp+fn
    dr.clin <- tp.clin / (tp.clin+fn.clin)
    tmp.clin <- data.frame("threshold"=x, 
                           "clin.fpr"=fpr.clin*100,
                           "clin.DR"=dr.clin)
  } else{
    tmp.clin <- data.frame("threshold"=x, 
                           "clin.fpr"=NA,
                           "clin.DR"=NA)
  }
  
  ## proteins - define predicted case status based on the selected threshold
  if(sum(lp.prots<x)>0&sum(lp.prots>x)>0){
    
    pred.case.prot <- ifelse(lp.prots>x,1,0)
    
    tp.prot <- sum(u.test[,dz]==1&pred.case.prot==1)
    fp.prot <- sum(u.test[,dz]==0&pred.case.prot==1)
    fn.prot <- sum(u.test[,dz]==1&pred.case.prot==0)
    tn.prot <- sum(u.test[,dz]==0&pred.case.prot==0)
    ## fpr = fp / tn + fp 
    fpr.prot <- fp.prot / (tn.prot + fp.prot)
    ## dr = tp / tp+fn
    dr.prot <- tp.prot / (tp.prot+fn.prot)
    
    tmp.prots <- data.frame("threshold"=x, 
                            "prot.fpr"=fpr.prot*100,
                            "prot.DR"=dr.prot)
  } else{
    tmp.prots <- data.frame("threshold"=x, 
                            "prot.fpr"=NA,
                            "prot.DR"=NA)
  }
  
  # if(sum(lp.biom<x)>0&sum(lp.biom>x)>0){
  #   ## biomarkers - define predicted case status based on the selected threshold
  #   pred.case.biom <- ifelse(lp.biom>x,1,0)
  #   
  #   tp.biom <- sum(u.test[,dz]==1&pred.case.biom==1)
  #   fp.biom <- sum(u.test[,dz]==0&pred.case.biom==1)
  #   fn.biom <- sum(u.test[,dz]==1&pred.case.biom==0)
  #   tn.biom <- sum(u.test[,dz]==0&pred.case.biom==0)
  #   ## fpr = fp / tn + fp 
  #   fpr.biom <- fp.biom / (tn.biom + fp.biom)
  #   ## dr = tp / tp+fn
  #   dr.biom <- tp.biom / (tp.biom+fn.biom)
  #   tmp.biom <- data.frame("threshold"=x, 
  #                           "biom.fpr"=fpr.biom*100,
  #                           "biom.DR"=dr.biom*100)
  #   } else{
  #     tmp.biom <- data.frame("threshold"=x,
  #                            "biom.fpr"=NA,
  #                            "biom.DR"=NA)
  #   }
  tmp <- merge(tmp.clin, tmp.prots, by="threshold")
  # tmp <- merge(tmp, tmp.biom, by="threshold")
  ## return object 
  return(tmp)
})

thresh.grid <- do.call(rbind, thresh.grid)

## subset to the closest thresholds to 5 - 50 % fpr

thresh.clin <- lapply(c(5,10,15,20,25,30,35,40,45,50), function(x){
  tmp.clin <- thresh.grid[which.min(abs(x-thresh.grid$clin.fpr)),1:3]
  return(tmp.clin)
})
thresh.clin <- as.data.frame(do.call(rbind, thresh.clin))  
thresh.clin$FPR.true <- c(5,10,15,20,25,30,35,40,45,50) 

thresh.prot <- lapply(c(5,10,15,20,25,30,35,40,45,50), function(x){
  tmp.prot <- thresh.grid[which.min(abs(x-thresh.grid$prot.fpr)),c(1,4,5)]
  return(tmp.prot)
})
thresh.prot <- as.data.frame(do.call(rbind, thresh.prot))  
thresh.prot$FPR.true <- c(5,10,15,20,25,30,35,40,45,50) 

# thresh.biom <- lapply(c(5,10,15,20,25,30,35,40,45,50), function(x){
#   tmp.biom <- thresh.grid[which.min(abs(x-thresh.grid$biom.fpr)),c(1,6,7)]
#   return(tmp.biom)
# })
# thresh.biom <- as.data.frame(do.call(rbind, thresh.biom))  
# thresh.biom$FPR.true <- c(5,10,15,20,25,30,35,40,45,50) 

res.dr.curve.test <- merge(thresh.clin, thresh.prot, by="FPR.true", suffixes=c(".clin",".prot"))
# res.dr.curve.test <- merge(res.dr.curve.test, thresh.biom, by="FPR.true")
# colnames(res.dr.curve.test)[8] <- c("threshold.biom") 
res.dr.curve.test$set <- rep("test.set", nrow(res.dr.curve.test))

####---- DR in the entire UKB subcohort ----####
## exclude consotrium selected cases
# consort <- which(ukbb$randomized_baseline==F)
# ukbb <- ukbb[-consort,]

lp.clin <- as.numeric(predict(res.clin$glmnet.opt, 
                              type = "link", 
                              newx = as.matrix(ukbb[, row.names(res.clin$opt.coefficients)]),
                              s = res.clin$lambda.opt))

lp.prots <- as.numeric(predict(res.clin.prots$glmnet.opt, 
                               type = "link", 
                               newx = as.matrix(ukbb[, row.names(res.clin.prots$opt.coefficients)]),
                               s = res.clin.prots$lambda.opt))

# lp.biom <- as.numeric(predict(res.clin.biom$glmnet.opt, 
#                                type = "link", 
#                                newx = as.matrix(ukbb[, row.names(res.clin.biom$opt.coefficients)]),
#                                s = res.clin.biom$lambda.opt))

## normalise 0- 1
if(min(lp.clin)<0){
  lp.clin <- (lp.clin+abs(min(lp.clin)))/
    max((lp.clin+abs(min(lp.clin))))
} else{
  lp.clin <- lp.clin/max(lp.clin)
}

if(min(lp.prots)<0){
  lp.prots <- (lp.prots+abs(min(lp.prots)))/
    max((lp.prots+abs(min(lp.prots))))
} else{
  lp.prots <- (lp.prots)/max(lp.prots)
}

# if(min(lp.biom)<0){
#   lp.biom <- (lp.biom+abs(min(lp.biom)))/
#     max((lp.biom+abs(min(lp.biom))))
# } else{
#   lp.biom <- (lp.biom)/max(lp.biom)
# }


## find threshold for detection rate at 5% false positive rate - test set only 
thresh.grid <- lapply(seq(0.0001,1,0.0001), function(x){
  
  if(sum(lp.clin<x)>0&sum(lp.clin>x)>0){
    ## clinical - define predicted case status based on the selected threshold
    pred.case.clin <- ifelse(lp.clin>x,1,0)
    
    tp.clin <- sum(ukbb[,dz]==1&pred.case.clin==1)
    fp.clin <- sum(ukbb[,dz]==0&pred.case.clin==1)
    fn.clin <- sum(ukbb[,dz]==1&pred.case.clin==0)
    tn.clin <- sum(ukbb[,dz]==0&pred.case.clin==0)
    ## fpr = fp / tn + fp 
    fpr.clin <- fp.clin / (tn.clin + fp.clin)
    ## dr = tp / tp+fn
    dr.clin <- tp.clin / (tp.clin+fn.clin)
    tmp.clin <- data.frame("threshold"=x, 
                           "clin.fpr"=fpr.clin*100,
                           "clin.DR"=dr.clin)
  } else{
    tmp.clin <- data.frame("threshold"=x, 
                           "clin.fpr"=NA,
                           "clin.DR"=NA)
  }
  
  ## proteins - define predicted case status based on the selected threshold
  if(sum(lp.prots<x)>0&sum(lp.prots>x)>0){
    
    pred.case.prot <- ifelse(lp.prots>x,1,0)
    
    tp.prot <- sum(ukbb[,dz]==1&pred.case.prot==1)
    fp.prot <- sum(ukbb[,dz]==0&pred.case.prot==1)
    fn.prot <- sum(ukbb[,dz]==1&pred.case.prot==0)
    tn.prot <- sum(ukbb[,dz]==0&pred.case.prot==0)
    ## fpr = fp / tn + fp 
    fpr.prot <- fp.prot / (tn.prot + fp.prot)
    ## dr = tp / tp+fn
    dr.prot <- tp.prot / (tp.prot+fn.prot)
    
    tmp.prots <- data.frame("threshold"=x, 
                            "prot.fpr"=fpr.prot*100,
                            "prot.DR"=dr.prot)
  } else{
    tmp.prots <- data.frame("threshold"=x, 
                            "prot.fpr"=NA,
                            "prot.DR"=NA)
  }
  
  # if(sum(lp.biom<x)>0&sum(lp.biom>x)>0){
  #   ## biomarkers - define predicted case status based on the selected threshold
  #   pred.case.biom <- ifelse(lp.biom>x,1,0)
  #   
  #   tp.biom <- sum(u.test[,dz]==1&pred.case.biom==1)
  #   fp.biom <- sum(u.test[,dz]==0&pred.case.biom==1)
  #   fn.biom <- sum(u.test[,dz]==1&pred.case.biom==0)
  #   tn.biom <- sum(u.test[,dz]==0&pred.case.biom==0)
  #   ## fpr = fp / tn + fp 
  #   fpr.biom <- fp.biom / (tn.biom + fp.biom)
  #   ## dr = tp / tp+fn
  #   dr.biom <- tp.biom / (tp.biom+fn.biom)
  #   tmp.biom <- data.frame("threshold"=x, 
  #                          "biom.fpr"=fpr.biom*100,
  #                          "biom.DR"=dr.biom*100)
  # } else{
  #   tmp.biom <- data.frame("threshold"=x,
  #                          "biom.fpr"=NA,
  #                          "biom.DR"=NA)
  # }
  tmp <- merge(tmp.clin, tmp.prots, by="threshold")
  # tmp <- merge(tmp, tmp.biom, by="threshold")
  ## return object 
  return(tmp)
})


thresh.grid <- do.call(rbind, thresh.grid)

## subset to the closest thresholds to 5 - 50 % fpr

thresh.clin <- lapply(c(5,10,15,20,25,30,35,40,45,50), function(x){
  tmp.clin <- thresh.grid[which.min(abs(x-thresh.grid$clin.fpr)),1:3]
  return(tmp.clin)
})
thresh.clin <- as.data.frame(do.call(rbind, thresh.clin))  
thresh.clin$FPR.true <- c(5,10,15,20,25,30,35,40,45,50) 

thresh.prot <- lapply(c(5,10,15,20,25,30,35,40,45,50), function(x){
  tmp.prot <- thresh.grid[which.min(abs(x-thresh.grid$prot.fpr)),c(1,4,5)]
  return(tmp.prot)
})
thresh.prot <- as.data.frame(do.call(rbind, thresh.prot))  
thresh.prot$FPR.true <- c(5,10,15,20,25,30,35,40,45,50) 

# thresh.biom <- lapply(c(5,10,15,20,25,30,35,40,45,50), function(x){
#   tmp.biom <- thresh.grid[which.min(abs(x-thresh.grid$biom.fpr)),c(1,6,7)]
#   return(tmp.biom)
# })
# thresh.biom <- as.data.frame(do.call(rbind, thresh.biom))  
# thresh.biom$FPR.true <- c(5,10,15,20,25,30,35,40,45,50) 

res.dr.curve.all <- merge(thresh.clin, thresh.prot, by="FPR.true", suffixes=c(".clin",".prot"))
# res.dr.curve.all <- merge(res.dr.curve.all, thresh.biom, by="FPR.true")
# colnames(res.dr.curve.all)[8] <- c("threshold.biom") 
res.dr.curve.all$set <- rep("all", nrow(res.dr.curve.all))
res.dr.curve <- rbind(res.dr.curve.all, res.dr.curve.test)
write.table(res.dr.curve, file=paste0("../data_output/detection_rate_curves/DR_curves_" ,dz, "_", inc.yrs,"years.txt"), sep = "\t", row.names = F)
###################################################
#### Script Name: Detection_rate_curves
#### Author: Julia Carrasco Zanini Sanchez
#### Description: script to ML framework for incident diseases in atlas 308 
#### Date : 23/03/2023

rm(list = ls())

## load packages 
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(RNOmni)
# library(glue)
library(DBI)
library(dbplyr)
library(magrittr)
# library(broom)
library(survival)
library(data.table)
library(glue)
library(parallel)
library(caret)
library(glmnet)
library(ROSE)
# library(odbc)

## define disease and restrict to incident cases within X years
args <- commandArgs(trailingOnly = T)

dz <- as.character(args[1])
dz.date <- as.character(args[2])
inc.yrs <- as.numeric(args[3])

paste0("Prediction:", dz, "; ", dz.date, " at ", inc.yrs, "years")

## set up working directory 
setwd("/home/ukbb_inc_dz_prediction/Explore_1536_Expansion_analyses/02_inc_dz_prediction/bin/") 

## read in data 
atlas308 <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_atlas308_annotations.txt", data.table = F)
cov <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_covariate_annotations.txt", data.table = F)

####----  Data Cleaning and wrangling   ----####
atlas308 <- as.data.frame(atlas308 %>% 
                            filter(eid_20361%in%cov$eid_20361) %>% 
                            pivot_wider(names_from = variable, values_from =value))

## phenotypes missing for patient-derived info model: bmi, smoking, alcohol, physical activity, family history 
pheno <- merge(cov[,-c(11:14,17)], atlas308, by="eid_20361")

## read in px info phenotypes
px.info <- read.delim("../../01_inc_dz_phenos/data_input/Patient_model_variables.txt",sep = "\t")
pheno <- merge(pheno, px.info, by="eid_20361")

## rename to olds names to match the rest fo the script
colnames(pheno)[which(colnames(pheno)=="assess_date")] <- c("date_baseline_assessment")

## process ethnicity
pheno$ethnicity <- ifelse(is.na(pheno$ge_eur),2,pheno$ge_eur)

## get protein data and metadata
ol.npx <- fread("../../00_protein_imputation/data_input/Imputed_1536_Expansion_NPX_proteins.txt", data.table = F)
ol.npx <- filter(ol.npx, exclusion.50==0 )

## exclude participants with missing age, sex
ol.npx <- ol.npx[complete.cases(ol.npx[,c("ge_sex", "ge_age")]),]

protein_list <- colnames(ol.npx)[-which(colnames(ol.npx)%in%c("eid_20361","ge_sex","ge_age",
                                                              "randomized_baseline","exclusion.50"))]

## merge phenotypes and proteins
ukbb <- merge(ol.npx, pheno, by="eid_20361")
## randomised baseline
ukbb <- ukbb %>% 
  filter(randomized_baseline==T)
## inverse rank normal transform npx
ukbb[,protein_list] <- sapply(ukbb[,protein_list], RankNorm)

## exclude prevalent cases and incident cases up to 6 months
## restrict to incident cases within the first X years
ukbb <- ukbb %>% 
  mutate(fol = difftime(!!as.name(dz.date),date_baseline_assessment, units = "weeks")/52.25) %>% 
  filter((fol>0.5&fol<inc.yrs)|is.na(!!as.name(dz.date)))

## censor (including censoring by death)
death.data <- read.delim("../../01_inc_dz_phenos/data_input/Death_dates_proteomic_participants.txt",sep = "\t")
death.data$date_of_death <- as.Date(death.data$date_of_death, format = "%Y-%m-%d")

ukbb <- merge(ukbb, death.data[,c("eid_20361","date_of_death")], by="eid_20361", all.x=T)

if(inc.yrs==10){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2020-12-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
} else if(inc.yrs==5){
  ukbb <- ukbb %>%
    ## sanity check - exclude cases with missing incident dates
    filter(!(is.na(!!as.name(dz.date))&!!as.name(dz)==1)) %>% 
    mutate(fol = ifelse(!!as.name(dz)==0,difftime(as.Date("2016-05-31"), date_baseline_assessment, units = "weeks")/52.25,
                        fol))
}else{stop("limit to incident casese within X years not defined")}

if(inc.yrs==10){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2020-12-31"),
                       difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                       ukbb$fol)
} else if(inc.yrs==5){
  ukbb$fol_d <- ifelse(ukbb[,dz]==0&ukbb$date_of_death<as.Date("2016-05-31"),
                       difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                       ukbb$fol)
}else{stop("limit to incident casese within X years not defined")}

ukbb$fol <- ifelse(is.na(ukbb$fol_d), ukbb$fol, ukbb$fol_d)

## for sex specific outcome exclude males/females
if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_pri_breast",
           "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical",
           "atlas_inc_endometrial_hyper","atlas_inc_endometriosis","atlas_inc_female_genital_prolapse",
           "atlas_inc_leiomyoma","atlas_inc_menorrhagia",
           "atlas_inc_pid","atlas_inc_pmb")){
  # ukbb <- ukbb %>% filter(ge_sex=="female")
  ukbb <- ukbb %>% filter(ge_sex==0)
} else if(dz%in%c("atlas_inc_pri_prost","atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  # ukbb <- ukbb %>% filter(ge_sex=="male")
  ukbb <- ukbb %>% filter(ge_sex==1)
} else{ukbb <- ukbb}

## recode sex to numeric
# ukbb$ge_sex <- ifelse(ukbb$ge_sex=="male",1,
#                       ifelse(ukbb$ge_sex=="female",0,NA))

n.cases <- sum(ukbb[,dz]==1)

## add noise variables
rand.vars <- lapply(1:10, function(x){
  tmp <- runif(nrow(ukbb), min = min(ukbb$NPPB_20049), 
               max = max(ukbb$NPPB_20049))
  return(tmp)
})
rand.vars <- as.data.frame(do.call(cbind, rand.vars))
colnames(rand.vars) <- paste0("rand_var_",1:10)

ukbb <- cbind(ukbb, rand.vars)

####----    ML framework    ----####

## split into train/optimize/test or train/test depending on case number
ukbb[,dz] <- as.factor(ukbb[,dz])

if(n.cases/4>200){
  ## feature selection set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.5, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.opti <- ukbb[-trainIndex,]
  
  ## optimization and testing
  set.seed(3456) ## for reproducibility
  trainIndex2 <- createDataPartition(u.opti[,dz], p=0.5, list = F, times = 1)
  u.test <- u.opti[trainIndex2,]
  u.opti <- u.opti[-trainIndex2,]
}else {
  ## feature selection set and test set
  set.seed(3456) ## for reproducibility
  trainIndex <- createDataPartition(ukbb[,dz], p=0.7, list = F, times = 1)
  u.train <- ukbb[trainIndex,]
  u.test <- ukbb[-trainIndex,]
}
  
####---- optimization and performance testing ----####

## survival objects
surv.train <- Surv(u.train$fol, u.train[,dz])
## optimization 
if(n.cases/4>200){
  surv.opt <- Surv(u.opti$fol, u.opti[,dz])
}else{surv.opt <- Surv(u.train$fol, u.train[,dz])}
## testing
surv.test <- Surv(u.test$fol, u.test[,dz])

## px information model basic 
clin.predictors <- c("ge_age","ge_sex","bmi","ethnicity","smoke_never",
                     "smoke_previous","smoke_current","alcohol_daily",
                     "alcohol_4_week","alcohol_2_week","alcohol_3_month","alcohol_occassion",
                     "alcohol_never")

## define additional family history for specific diseases
if(!dz%in%c("atlas_inc_diabetes_t2","atlas_inc_pri_prost", "atlas_inc_copd","atlas_inc_isch_stroke",
            "atlas_inc_af","atlas_inc_av_block_3","atlas_inc_chd_nos","atlas_inc_hf", 
            "atlas_inc_peripheral_arterial_disease","atlas_inc_parkinsons","atlas_inc_pri_lung",
            "atlas_inc_pri_bowel", "atlas_inc_pri_breast","atlas_inc_pri_ovarian","atlas_inc_pri_uterine",
            "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
            "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
            "atlas_inc_hydrocele","atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb",
            "atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu")){
  clin.predictors <- clin.predictors
} else if(dz=="atlas_inc_diabetes_t2") {
  clin.predictors <- c(clin.predictors,"father_diabetes","mother_diabetes")
} else if(dz=="atlas_inc_pri_prost") {
  clin.predictors <- c(clin.predictors[-2],"father_prostate_cancer")
} else if(dz=="atlas_inc_copd") {
  clin.predictors <- c(clin.predictors,"father_COPD","mother_COPD")
} else if(dz%in%c("atlas_inc_isch_stroke","atlas_inc_af","atlas_inc_av_block_3",
                  "atlas_inc_chd_nos","atlas_inc_hf","atlas_inc_peripheral_arterial_disease")) {
  clin.predictors <- c(clin.predictors,"father_heart_disease","mother_heart_disease")
} else if(dz=="atlas_inc_parkinsons") {
  clin.predictors <- c(clin.predictors,"father_parkinsons","mother_parkinsons")
} else if(dz=="atlas_inc_pri_lung") {
  clin.predictors <- c(clin.predictors,"father_lung_cancer","mother_lung_cancer")
} else if(dz=="atlas_inc_pri_bowel") {
  clin.predictors <- c(clin.predictors,"father_bowel_cancer","mother_bowel_cancer")
} else if(dz=="atlas_inc_pri_breast") {
  clin.predictors <- c(clin.predictors[-2],"mother_breast_cancer")
} else if(dz=="atlas_inc_hypertension") {
  clin.predictors <- c(clin.predictors[-2],"mother_high_bp","father_high_bp")
} else if(dz=="atlas_inc_depression") {
  clin.predictors <- c(clin.predictors[-2],"mother_depression","father_depression")
} else if(dz%in%c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_benign_ovary",
                  "atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
                  "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
                  "atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb")) {
  clin.predictors <- c(clin.predictors[-2])
} else if(dz%in%c("atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")){
  clin.predictors <- c(clin.predictors[-2])
}

## load results for best performing protein model 
summary.res <- read.delim("../data_output/summary_results_10_year_models.txt", sep = "\t")
summary.res <- summary.res[which(summary.res$dz==dz),]

if (summary.res$model=="clin.prots.all"){
  load(paste0("../data_output/px_info_protein_models/px_info_all_prots_above_noise_res_",dz,"_",inc.yrs,"yrs",".RData"))
  res.clin.prots <- res.clin.prots.all
  } else if (summary.res$model=="clin.prots10") {
    load(paste0("../data_output/px_info_protein_models/px_info_top10_prots_res_",dz,"_",inc.yrs,"yrs",".RData"))
    res.clin.prots <- res.clin.prots10
    } else if (summary.res$model=="clin.prots") {
      load(paste0("../data_output/px_info_protein_models/px_info_top5_prots_res_",dz,"_",inc.yrs,"yrs",".RData"))
      res.clin.prots <- res.clin.prots
      } else if (summary.res$model=="prot") {
        load(paste0("../data_output/protein_models/prot_top5_res_",dz,"_",inc.yrs,"yrs",".RData"))
        res.clin.prots <- res.prot
      }

## load results for best performing biomarker model 
# biom.res <- read.delim("../data_output/summary_biomarker_results_10_year_models.txt", sep = "\t")
# biom.res <- biom.res[which(biom.res$dz==dz),]
# 
# if (biom.res$model=="clin.biom.all"){
#   load(paste0("../data_output/px_info_biomarker_models/px_info_all_biom_above_noise_res_",dz,"_",inc.yrs,"yrs",".RData"))
#   res.clin.biom <- res.clin.biom.all
# } else if (biom.res$model=="clin.biom10") {
#   load(paste0("../data_output/px_info_biomarker_models/px_info_top10_biom_res_",dz,"_",inc.yrs,"yrs",".RData"))
#   res.clin.biom <- res.clin.biom10
# } else if (biom.res$model=="clin.biom") {
#   load(paste0("../data_output/px_info_biomarker_models/px_info_top5_biom_res_",dz,"_",inc.yrs,"yrs",".RData"))
#   res.clin.biom <- res.clin.biom
# } else if (biom.res$model=="biom") {
#   load(paste0("../data_output/biomarker_models/biomarker_top5_res_",dz,"_",inc.yrs,"yrs",".RData"))
#   res.clin.biom <- res.biom
# }

## clinical model 
load(paste0("../data_output/px_info_models/px_info_res_",dz,"_",inc.yrs,"yrs",".RData"))

####---- DR using test set only ----####
## normalise linear predictor (0-1 values)
if(min(res.clin$linear.predictor)<0){
  lp.clin <- (res.clin$linear.predictor+abs(min(res.clin$linear.predictor)))/
    max((res.clin$linear.predictor+abs(min(res.clin$linear.predictor))))
} else{
  lp.clin <- (res.clin$linear.predictor)/max(res.clin$linear.predictor)
}

if(min(res.clin.prots$linear.predictor)<0){
  lp.prots <- (res.clin.prots$linear.predictor+abs(min(res.clin.prots$linear.predictor)))/
    max((res.clin.prots$linear.predictor+abs(min(res.clin.prots$linear.predictor))))
} else{
  lp.prots <- (res.clin.prots$linear.predictor)/max(res.clin.prots$linear.predictor)
}

# if(min(res.clin.biom$linear.predictor)<0){
#   lp.biom <- (res.clin.biom$linear.predictor+abs(min(res.clin.biom$linear.predictor)))/
#     max((res.clin.biom$linear.predictor+abs(min(res.clin.biom$linear.predictor))))
# } else{
#   lp.biom <- (res.clin.biom$linear.predictor)/max(res.clin.biom$linear.predictor)
# }

## find threshold for detection rate at 5% false positive rate - test set only 
thresh.grid <- lapply(seq(0.0001,1,0.0001), function(x){
  
  if(sum(lp.clin<x)>0&sum(lp.clin>x)>0){
    ## clinical - define predicted case status based on the selected threshold
    pred.case.clin <- ifelse(lp.clin>x,1,0)
    
    tp.clin <- sum(u.test[,dz]==1&pred.case.clin==1)
    fp.clin <- sum(u.test[,dz]==0&pred.case.clin==1)
    fn.clin <- sum(u.test[,dz]==1&pred.case.clin==0)
    tn.clin <- sum(u.test[,dz]==0&pred.case.clin==0)
    ## fpr = fp / tn + fp 
    fpr.clin <- fp.clin / (tn.clin + fp.clin)
    ## dr = tp / tp+fn
    dr.clin <- tp.clin / (tp.clin+fn.clin)
    tmp.clin <- data.frame("threshold"=x, 
                           "clin.fpr"=fpr.clin*100,
                           "clin.DR"=dr.clin)
  } else{
    tmp.clin <- data.frame("threshold"=x, 
                           "clin.fpr"=NA,
                           "clin.DR"=NA)
  }
  
  ## proteins - define predicted case status based on the selected threshold
  if(sum(lp.prots<x)>0&sum(lp.prots>x)>0){
    
    pred.case.prot <- ifelse(lp.prots>x,1,0)
    
    tp.prot <- sum(u.test[,dz]==1&pred.case.prot==1)
    fp.prot <- sum(u.test[,dz]==0&pred.case.prot==1)
    fn.prot <- sum(u.test[,dz]==1&pred.case.prot==0)
    tn.prot <- sum(u.test[,dz]==0&pred.case.prot==0)
    ## fpr = fp / tn + fp 
    fpr.prot <- fp.prot / (tn.prot + fp.prot)
    ## dr = tp / tp+fn
    dr.prot <- tp.prot / (tp.prot+fn.prot)
    
    tmp.prots <- data.frame("threshold"=x, 
                           "prot.fpr"=fpr.prot*100,
                           "prot.DR"=dr.prot)
    } else{
    tmp.prots <- data.frame("threshold"=x, 
                           "prot.fpr"=NA,
                           "prot.DR"=NA)
  }
    
  # if(sum(lp.biom<x)>0&sum(lp.biom>x)>0){
  #   ## biomarkers - define predicted case status based on the selected threshold
  #   pred.case.biom <- ifelse(lp.biom>x,1,0)
  #   
  #   tp.biom <- sum(u.test[,dz]==1&pred.case.biom==1)
  #   fp.biom <- sum(u.test[,dz]==0&pred.case.biom==1)
  #   fn.biom <- sum(u.test[,dz]==1&pred.case.biom==0)
  #   tn.biom <- sum(u.test[,dz]==0&pred.case.biom==0)
  #   ## fpr = fp / tn + fp 
  #   fpr.biom <- fp.biom / (tn.biom + fp.biom)
  #   ## dr = tp / tp+fn
  #   dr.biom <- tp.biom / (tp.biom+fn.biom)
  #   tmp.biom <- data.frame("threshold"=x, 
  #                           "biom.fpr"=fpr.biom*100,
  #                           "biom.DR"=dr.biom*100)
  #   } else{
  #     tmp.biom <- data.frame("threshold"=x,
  #                            "biom.fpr"=NA,
  #                            "biom.DR"=NA)
  #   }
  tmp <- merge(tmp.clin, tmp.prots, by="threshold")
  # tmp <- merge(tmp, tmp.biom, by="threshold")
  ## return object 
  return(tmp)
  })

thresh.grid <- do.call(rbind, thresh.grid)

## subset to the closest thresholds to 5 - 50 % fpr

thresh.clin <- lapply(c(5,10,15,20,25,30,35,40,45,50), function(x){
  tmp.clin <- thresh.grid[which.min(abs(x-thresh.grid$clin.fpr)),1:3]
  return(tmp.clin)
})
thresh.clin <- as.data.frame(do.call(rbind, thresh.clin))  
thresh.clin$FPR.true <- c(5,10,15,20,25,30,35,40,45,50) 
  
thresh.prot <- lapply(c(5,10,15,20,25,30,35,40,45,50), function(x){
  tmp.prot <- thresh.grid[which.min(abs(x-thresh.grid$prot.fpr)),c(1,4,5)]
  return(tmp.prot)
})
thresh.prot <- as.data.frame(do.call(rbind, thresh.prot))  
thresh.prot$FPR.true <- c(5,10,15,20,25,30,35,40,45,50) 

# thresh.biom <- lapply(c(5,10,15,20,25,30,35,40,45,50), function(x){
#   tmp.biom <- thresh.grid[which.min(abs(x-thresh.grid$biom.fpr)),c(1,6,7)]
#   return(tmp.biom)
# })
# thresh.biom <- as.data.frame(do.call(rbind, thresh.biom))  
# thresh.biom$FPR.true <- c(5,10,15,20,25,30,35,40,45,50) 

res.dr.curve.test <- merge(thresh.clin, thresh.prot, by="FPR.true", suffixes=c(".clin",".prot"))
# res.dr.curve.test <- merge(res.dr.curve.test, thresh.biom, by="FPR.true")
# colnames(res.dr.curve.test)[8] <- c("threshold.biom") 
res.dr.curve.test$set <- rep("test.set", nrow(res.dr.curve.test))
  
####---- DR in the entire UKB subcohort ----####
lp.clin <- as.numeric(predict(res.clin$glmnet.opt, 
                              type = "link", 
                              newx = as.matrix(ukbb[, row.names(res.clin$opt.coefficients)]),
                              s = res.clin$lambda.opt))

lp.prots <- as.numeric(predict(res.clin.prots$glmnet.opt, 
                                  type = "link", 
                                  newx = as.matrix(ukbb[, row.names(res.clin.prots$opt.coefficients)]),
                                  s = res.clin.prots$lambda.opt))

# lp.biom <- as.numeric(predict(res.clin.biom$glmnet.opt, 
#                                type = "link", 
#                                newx = as.matrix(ukbb[, row.names(res.clin.biom$opt.coefficients)]),
#                                s = res.clin.biom$lambda.opt))

## normalise 0- 1
if(min(lp.clin)<0){
  lp.clin <- (lp.clin+abs(min(lp.clin)))/
    max((lp.clin+abs(min(lp.clin))))
} else{
  lp.clin <- lp.clin/max(lp.clin)
}

if(min(lp.prots)<0){
  lp.prots <- (lp.prots+abs(min(lp.prots)))/
    max((lp.prots+abs(min(lp.prots))))
} else{
  lp.prots <- (lp.prots)/max(lp.prots)
}

# if(min(lp.biom)<0){
#   lp.biom <- (lp.biom+abs(min(lp.biom)))/
#     max((lp.biom+abs(min(lp.biom))))
# } else{
#   lp.biom <- (lp.biom)/max(lp.biom)
# }


## find threshold for detection rate at 5% false positive rate - test set only 
thresh.grid <- lapply(seq(0.0001,1,0.0001), function(x){
  
  if(sum(lp.clin<x)>0&sum(lp.clin>x)>0){
    ## clinical - define predicted case status based on the selected threshold
    pred.case.clin <- ifelse(lp.clin>x,1,0)
    
    tp.clin <- sum(ukbb[,dz]==1&pred.case.clin==1)
    fp.clin <- sum(ukbb[,dz]==0&pred.case.clin==1)
    fn.clin <- sum(ukbb[,dz]==1&pred.case.clin==0)
    tn.clin <- sum(ukbb[,dz]==0&pred.case.clin==0)
    ## fpr = fp / tn + fp 
    fpr.clin <- fp.clin / (tn.clin + fp.clin)
    ## dr = tp / tp+fn
    dr.clin <- tp.clin / (tp.clin+fn.clin)
    tmp.clin <- data.frame("threshold"=x, 
                           "clin.fpr"=fpr.clin*100,
                           "clin.DR"=dr.clin)
  } else{
    tmp.clin <- data.frame("threshold"=x, 
                           "clin.fpr"=NA,
                           "clin.DR"=NA)
  }
  
  ## proteins - define predicted case status based on the selected threshold
  if(sum(lp.prots<x)>0&sum(lp.prots>x)>0){
    
    pred.case.prot <- ifelse(lp.prots>x,1,0)
    
    tp.prot <- sum(ukbb[,dz]==1&pred.case.prot==1)
    fp.prot <- sum(ukbb[,dz]==0&pred.case.prot==1)
    fn.prot <- sum(ukbb[,dz]==1&pred.case.prot==0)
    tn.prot <- sum(ukbb[,dz]==0&pred.case.prot==0)
    ## fpr = fp / tn + fp 
    fpr.prot <- fp.prot / (tn.prot + fp.prot)
    ## dr = tp / tp+fn
    dr.prot <- tp.prot / (tp.prot+fn.prot)
    
    tmp.prots <- data.frame("threshold"=x, 
                            "prot.fpr"=fpr.prot*100,
                            "prot.DR"=dr.prot)
  } else{
    tmp.prots <- data.frame("threshold"=x, 
                            "prot.fpr"=NA,
                            "prot.DR"=NA)
  }
  
  # if(sum(lp.biom<x)>0&sum(lp.biom>x)>0){
  #   ## biomarkers - define predicted case status based on the selected threshold
  #   pred.case.biom <- ifelse(lp.biom>x,1,0)
  #   
  #   tp.biom <- sum(u.test[,dz]==1&pred.case.biom==1)
  #   fp.biom <- sum(u.test[,dz]==0&pred.case.biom==1)
  #   fn.biom <- sum(u.test[,dz]==1&pred.case.biom==0)
  #   tn.biom <- sum(u.test[,dz]==0&pred.case.biom==0)
  #   ## fpr = fp / tn + fp 
  #   fpr.biom <- fp.biom / (tn.biom + fp.biom)
  #   ## dr = tp / tp+fn
  #   dr.biom <- tp.biom / (tp.biom+fn.biom)
  #   tmp.biom <- data.frame("threshold"=x, 
  #                          "biom.fpr"=fpr.biom*100,
  #                          "biom.DR"=dr.biom*100)
  # } else{
  #   tmp.biom <- data.frame("threshold"=x,
  #                          "biom.fpr"=NA,
  #                          "biom.DR"=NA)
  # }
  tmp <- merge(tmp.clin, tmp.prots, by="threshold")
  # tmp <- merge(tmp, tmp.biom, by="threshold")
  ## return object 
  return(tmp)
})


thresh.grid <- do.call(rbind, thresh.grid)

## subset to the closest thresholds to 5 - 50 % fpr

thresh.clin <- lapply(c(5,10,15,20,25,30,35,40,45,50), function(x){
  tmp.clin <- thresh.grid[which.min(abs(x-thresh.grid$clin.fpr)),1:3]
  return(tmp.clin)
})
thresh.clin <- as.data.frame(do.call(rbind, thresh.clin))  
thresh.clin$FPR.true <- c(5,10,15,20,25,30,35,40,45,50) 

thresh.prot <- lapply(c(5,10,15,20,25,30,35,40,45,50), function(x){
  tmp.prot <- thresh.grid[which.min(abs(x-thresh.grid$prot.fpr)),c(1,4,5)]
  return(tmp.prot)
})
thresh.prot <- as.data.frame(do.call(rbind, thresh.prot))  
thresh.prot$FPR.true <- c(5,10,15,20,25,30,35,40,45,50) 

# thresh.biom <- lapply(c(5,10,15,20,25,30,35,40,45,50), function(x){
#   tmp.biom <- thresh.grid[which.min(abs(x-thresh.grid$biom.fpr)),c(1,6,7)]
#   return(tmp.biom)
# })
# thresh.biom <- as.data.frame(do.call(rbind, thresh.biom))  
# thresh.biom$FPR.true <- c(5,10,15,20,25,30,35,40,45,50) 

res.dr.curve.all <- merge(thresh.clin, thresh.prot, by="FPR.true", suffixes=c(".clin",".prot"))
# res.dr.curve.all <- merge(res.dr.curve.all, thresh.biom, by="FPR.true")
# colnames(res.dr.curve.all)[8] <- c("threshold.biom") 
res.dr.curve.all$set <- rep("all", nrow(res.dr.curve.all))
res.dr.curve <- rbind(res.dr.curve.all, res.dr.curve.test)
write.table(res.dr.curve, file=paste0("../data_output/detection_rate_curves/DR_curves_" ,dz, "_", inc.yrs,"years.txt"), sep = "\t", row.names = F)
