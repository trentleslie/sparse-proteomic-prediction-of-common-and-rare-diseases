# Utility Functions Script

library(dplyr)

get_clinical_predictors <- function(dz) {
  clin.predictors <- c("ge_age","ge_sex","bmi","ethnicity","smoke_never",
                       "smoke_previous","smoke_current","alcohol_daily",
                       "alcohol_4_week","alcohol_2_week","alcohol_3_month","alcohol_occassion",
                       "alcohol_never")
  
  additional_predictors <- switch(dz,
    "atlas_inc_diabetes_t2" = c("father_diabetes","mother_diabetes"),
    "atlas_inc_pri_prost" = c("father_prostate_cancer"),
    "atlas_inc_copd" = c("father_COPD","mother_COPD"),
    "atlas_inc_isch_stroke" = ,
    "atlas_inc_af" = ,
    "atlas_inc_av_block_3" = ,
    "atlas_inc_chd_nos" = ,
    "atlas_inc_hf" = ,
    "atlas_inc_peripheral_arterial_disease" = c("father_heart_disease","mother_heart_disease"),
    "atlas_inc_parkinsons" = c("father_parkinsons","mother_parkinsons"),
    "atlas_inc_pri_lung" = c("father_lung_cancer","mother_lung_cancer"),
    "atlas_inc_pri_bowel" = c("father_bowel_cancer","mother_bowel_cancer"),
    "atlas_inc_pri_breast" = c("mother_breast_cancer"),
    "atlas_inc_hypertension" = c("mother_high_bp","father_high_bp"),
    "atlas_inc_depression" = c("mother_depression","father_depression"),
    character(0)
  )
  
  if (dz %in% c("atlas_inc_pri_prost", "atlas_inc_pri_breast", "atlas_inc_pri_ovarian","atlas_inc_pri_uterine","atlas_inc_benign_ovary",
                "atlas_inc_benign_uterus","atlas_inc_cin_cervical","atlas_inc_endometrial_hyper",
                "atlas_inc_endometriosis","atlas_inc_female_genital_prolapse","atlas_inc_leiomyoma",
                "atlas_inc_menorrhagia","atlas_inc_pid","atlas_inc_pmb",
                "atlas_inc_bph","atlas_inc_ed","atlas_inc_male_gu","atlas_inc_hydrocele")) {
    clin.predictors <- clin.predictors[-which(clin.predictors == "ge_sex")]
  }
  
  c(clin.predictors, additional_predictors)
}

split_data <- function(ukbb, dz, n.cases) {
  ukbb[,dz] <- as.factor(ukbb[,dz])
  
  if(n.cases/4 > 200) {
    set.seed(3456)
    trainIndex <- createDataPartition(ukbb[,dz], p = 0.5, list = F, times = 1)
    u.train <- ukbb[trainIndex,]
    u.opti <- ukbb[-trainIndex,]
    
    set.seed(3456)
    trainIndex2 <- createDataPartition(u.opti[,dz], p = 0.5, list = F, times = 1)
    u.test <- u.opti[trainIndex2,]
    u.opti <- u.opti[-trainIndex2,]
  } else {
    set.seed(3456)
    trainIndex <- createDataPartition(ukbb[,dz], p = 0.7, list = F, times = 1)
    u.train <- ukbb[trainIndex,]
    u.test <- ukbb[-trainIndex,]
    u.opti <- u.train  # In this case, u.opti is the same as u.train
  }
  
  list(u.train = u.train, u.opti = u.opti, u.test = u.test)
}

create_survival_objects <- function(u_train, u_opti, u_test, dz) {
  surv.train <- Surv(u_train$fol, u_train[,dz])
  surv.opt <- if(nrow(u_opti) != nrow(u_train)) Surv(u_opti$fol, u_opti[,dz]) else surv.train
  surv.test <- Surv(u_test$fol, u_test[,dz])
  
  list(surv.train = surv.train, surv.opt = surv.opt, surv.test = surv.test)
}