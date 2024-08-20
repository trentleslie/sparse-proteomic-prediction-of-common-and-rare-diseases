# Model Training Script

library(survival)
library(glmnet)

train_model <- function(train_data, pred_vec, train_surv_data, times, boot_samples, test_data, test_surv_data) {
  source("../../../02_inc_dz_prediction/bin/functions/coxnet_ridge_optim_boot.R")
  
  res <- coxnet.optim.r(Train.data = train_data,
                        pred.vec = pred_vec,
                        Train.surv.data = train_surv_data,
                        times = times, 
                        boot.samples = boot_samples,
                        Test.data = test_data,
                        Test.surv.data = test_surv_data)
  
  return(res)
}

train_clinical_model <- function(u_opti, clin_predictors, surv_opt, boot_samples, u_test, surv_test) {
  res.clin <- train_model(u_opti, clin_predictors, surv_opt, 1000, boot_samples, u_test, surv_test)
  return(res.clin)
}

train_protein_model <- function(u_opti, clin_predictors, prots_opti, surv_opt, boot_samples, u_test, surv_test) {
  res.clin.prots <- train_model(u_opti, c(clin_predictors, prots_opti), surv_opt, 1000, boot_samples, u_test, surv_test)
  return(res.clin.prots)
}

train_biomarker_model <- function(u_opti, clin_predictors, biom_opti, surv_opt, boot_samples, u_test, surv_test) {
  res.clin.biom <- train_model(u_opti, c(clin_predictors, biom_opti), surv_opt, 1000, boot_samples, u_test, surv_test)
  return(res.clin.biom)
}