# Prediction and Evaluation Script

library(dplyr)

calculate_performance_metrics <- function(predicted, actual, thresholds) {
  lapply(thresholds, function(x) {
    pred_case <- ifelse(predicted > x, 1, 0)
    
    tp <- sum(actual == 1 & pred_case == 1)
    fp <- sum(actual == 0 & pred_case == 1)
    fn <- sum(actual == 1 & pred_case == 0)
    tn <- sum(actual == 0 & pred_case == 0)
    
    fpr <- fp / (tn + fp)
    dr <- tp / (tp + fn)
    
    data.frame(threshold = x, fpr = fpr * 100, dr = dr)
  }) %>% bind_rows()
}

generate_dr_curve <- function(res_clin, res_clin_prots, u_test, dz) {
  lp.clin <- normalize_linear_predictor(res_clin$linear.predictor)
  lp.prots <- normalize_linear_predictor(res_clin_prots$linear.predictor)
  
  thresholds <- seq(0.0001, 1, 0.0001)
  
  clin_metrics <- calculate_performance_metrics(lp.clin, u_test[[dz]], thresholds)
  prot_metrics <- calculate_performance_metrics(lp.prots, u_test[[dz]], thresholds)
  
  metrics <- merge(clin_metrics, prot_metrics, by = "threshold", suffixes = c(".clin", ".prot"))
  
  fpr_levels <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
  
  res.dr.curve <- lapply(fpr_levels, function(x) {
    clin <- metrics[which.min(abs(x - metrics$fpr.clin)), c("threshold", "fpr.clin", "dr.clin")]
    prot <- metrics[which.min(abs(x - metrics$fpr.prot)), c("threshold", "fpr.prot", "dr.prot")]
    merge(clin, prot, by = "threshold")
  }) %>% bind_rows()
  
  res.dr.curve$FPR.true <- fpr_levels
  res.dr.curve$set <- "test.set"
  
  return(res.dr.curve)
}

normalize_linear_predictor <- function(lp) {
  if(min(lp) < 0) {
    (lp + abs(min(lp))) / max((lp + abs(min(lp))))
  } else {
    lp / max(lp)
  }
}