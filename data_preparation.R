# Data Preparation Script

library(tidyr)
library(dplyr)
library(data.table)
library(RNOmni)

generate_input_file <- function() {
  setwd("/home/ukbb_inc_dz_prediction/Explore_1536_Expansion_analyses/01_inc_dz_phenos/bin/")
  
  dz.list <- read.delim("../data_input/atlas308_subc_inc6months_inc_dz_counts_3_10yrs.txt", sep = "\t")
  
  dz.input <- dz.list %>% 
    select(inc_dz_name, dz_date) %>% 
    mutate(inc_yrs=rep(10,nrow(dz.list)))
  
  write.table(dz.input, file = "../../02_inc_dz_prediction/data_input/atlas308_file_input_params.txt", sep="\t", col.names = F, row.names = F, quote = F)
  
  dz.list.cc <- read.delim("../data_input/atlas308_consort_inc6months_inc_dz_counts_3_10yrs.txt", sep = "\t")
  
  dz.input.cc <- dz.list.cc %>% 
    select(inc_dz_name, dz_date) %>% 
    mutate(inc_yrs=rep(10,nrow(dz.list.cc)))
  
  dz.input.cc <- rbind(dz.input.cc, dz.list.cc %>% 
                         select(inc_dz_name, dz_date) %>% 
                         mutate(inc_yrs=rep(5,nrow(dz.list.cc))))
  write.table(dz.input.cc, file = "../../02_inc_dz_prediction/data_input/atlas308_file_input_params_consort.txt", sep="\t", col.names = F, row.names = F, quote = F)
  
  prs.input <- dz.input %>%
    filter(inc_dz_name %in% c("atlas_inc_asthma","atlas_inc_af",
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
}

prepare_data <- function(dz, dz.date, inc.yrs) {
  setwd("/home/ukbb_inc_dz_prediction/Explore_1536_Expansion_analyses/02_inc_dz_prediction/bin/")
  
  atlas308 <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_atlas308_annotations.txt", data.table = F)
  cov <- fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_covariate_annotations.txt", data.table = F)
  
  atlas308 <- as.data.frame(atlas308 %>% 
                              filter(eid_20361 %in% cov$eid_20361) %>% 
                              pivot_wider(names_from = variable, values_from = value))
  
  pheno <- merge(cov[,-c(11:14,17)], atlas308, by="eid_20361")
  
  px.info <- read.delim("../../01_inc_dz_phenos/data_input/Patient_model_variables.txt", sep = "\t")
  pheno <- merge(pheno, px.info, by="eid_20361")
  
  colnames(pheno)[which(colnames(pheno)=="assess_date")] <- c("date_baseline_assessment")
  
  pheno$ethnicity <- ifelse(is.na(pheno$ge_eur), 2, pheno$ge_eur)
  
  ol.npx <- fread("../../00_protein_imputation/data_input/Imputed_1536_Expansion_NPX_proteins.txt", data.table = F)
  ol.npx <- filter(ol.npx, exclusion.50 == 0)
  
  ol.npx <- ol.npx[complete.cases(ol.npx[,c("ge_sex", "ge_age")]),]
  
  protein_list <- colnames(ol.npx)[-which(colnames(ol.npx) %in% c("eid_20361","ge_sex","ge_age","randomized_baseline","exclusion.50"))]
  
  ukbb <- merge(ol.npx, pheno, by="eid_20361")
  
  ukbb <- ukbb %>% 
    filter(randomized_baseline == T | !!as.name(dz) == 1)
  
  ukbb[,protein_list] <- sapply(ukbb[,protein_list], RankNorm)
  
  ukbb <- ukbb %>% 
    mutate(fol = difftime(!!as.name(dz.date), date_baseline_assessment, units = "weeks")/52.25) %>% 
    filter((fol > 0.5 & fol < inc.yrs) | is.na(!!as.name(dz.date)))
  
  death.data <- read.delim("../../01_inc_dz_phenos/data_input/Death_dates_proteomic_participants.txt", sep = "\t")
  death.data$date_of_death <- as.Date(death.data$date_of_death, format = "%Y-%m-%d")
  
  ukbb <- merge(ukbb, death.data[,c("eid_20361","date_of_death")], by="eid_20361", all.x=T)
  
  if(inc.yrs == 10) {
    ukbb <- ukbb %>%
      filter(!(is.na(!!as.name(dz.date)) & !!as.name(dz) == 1)) %>% 
      mutate(fol = ifelse(!!as.name(dz) == 0, difftime(as.Date("2020-12-31"), date_baseline_assessment, units = "weeks")/52.25, fol))
  } else if(inc.yrs == 5) {
    ukbb <- ukbb %>%
      filter(!(is.na(!!as.name(dz.date)) & !!as.name(dz) == 1)) %>% 
      mutate(fol = ifelse(!!as.name(dz) == 0, difftime(as.Date("2016-05-31"), date_baseline_assessment, units = "weeks")/52.25, fol))
  } else {
    stop("limit to incident cases within X years not defined")
  }
  
  if(inc.yrs == 10) {
    ukbb$fol_d <- ifelse(ukbb[,dz] == 0 & ukbb$date_of_death < as.Date("2020-12-31"),
                         difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                         ukbb$fol)
  } else if(inc.yrs == 5) {
    ukbb$fol_d <- ifelse(ukbb[,dz] == 0 & ukbb$date_of_death < as.Date("2016-05-31"),
                         difftime(ukbb$date_of_death, ukbb$date_baseline_assessment, units = "weeks")/52.25,
                         ukbb$fol)
  } else {
    stop("limit to incident cases within X years not defined")
  }
  
  ukbb$fol <- ifelse(is.na(ukbb$fol_d), ukbb$fol, ukbb$fol_d)
  
  if(dz %in% c("atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_pri_breast",
               "atlas_inc_benign_ovary","atlas_inc_benign_uterus","atlas_inc_cin_cervical",
               "atlas_inc_endometrial_hyper","atlas_inc_endometriosis","atlas_inc_female_genital_prolapse",
               "atlas_inc_leiomyoma","atlas_inc_menorrhagia",
               "atlas_inc_pid","atlas_inc_pmb")) {
    ukbb <- ukbb %>% filter(ge_sex == 0)
  } else if(dz %in% c("atlas_inc_pri_prost","atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")) {
    ukbb <- ukbb %>% filter(ge_sex == 1)
  }
  
  rand.vars <- lapply(1:10, function(x) {
    runif(nrow(ukbb), min = min(ukbb$NPPB_20049), max = max(ukbb$NPPB_20049))
  })
  rand.vars <- as.data.frame(do.call(cbind, rand.vars))
  colnames(rand.vars) <- paste0("rand_var_", 1:10)
  
  ukbb <- cbind(ukbb, rand.vars)
  
  return(ukbb)
}