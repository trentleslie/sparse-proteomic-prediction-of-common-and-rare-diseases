# UK Biobank Analysis Module

# Load necessary libraries
library(dplyr)
library(tidyr)
library(RNOmni)
library(survival)
library(caret)
library(glmnet)
library(ROSE)
library(missForest)
library(nricens)
library(survIDINRI)
library(variancePartition)
library(gprofiler2)

# Set seed for reproducibility
set.seed(123)

# Existing functions from ukbb_analysis_module.R
preprocess_data <- function(pheno, ol.npx, protein_list, dz, dz.date, inc.yrs) {
  # ... (keep the existing implementation)
}

censor_by_death <- function(ukbb, inc.yrs) {
  # ... (keep the existing implementation)
}

feature_selection <- function(train_data, predictors, dz) {
  # ... (keep the existing implementation)
}

coxnet_optim <- function(Train.data, pred.vec, Train.surv.data, times, boot.samples, Test.data, Test.surv.data) {
  # ... (keep the existing implementation)
}

calculate_detection_rate <- function(linear_predictor, true_status, thresholds) {
  # ... (keep the existing implementation)
}

run_analysis <- function(pheno, ol.npx, protein_list, dz, dz.date, inc.yrs) {
  # ... (keep the existing implementation)
}

# New combined function incorporating work from combined.R
combined_ukbb_analysis <- function(pheno, ol.npx, protein_list, dz, dz.date, inc.yrs) {
  # Step 1: Data Preprocessing (using existing preprocess_data function)
  ukbb <- preprocess_data(pheno, ol.npx, protein_list, dz, dz.date, inc.yrs)
  
  # Step 2: Impute missing values using missForest
  ukbb_imputed <- missForest(ukbb)$ximp
  
  # Step 3: Prepare the data for modeling
  x <- model.matrix(~ . - 1, data = ukbb_imputed[, -c(1, 2)])  # Exclude time and status columns from predictors
  y <- Surv(ukbb_imputed$fol, ukbb_imputed[,dz])
  
  # Step 4: Split the data into training and testing sets
  set.seed(3456)
  trainIndex <- createDataPartition(ukbb_imputed[,dz], p = 0.7, list = FALSE)
  trainData <- ukbb_imputed[trainIndex, ]
  testData <- ukbb_imputed[-trainIndex, ]
  
  x_train <- model.matrix(~ . - 1, data = trainData[, -c(1, 2)])  # Exclude time and status columns from predictors
  y_train <- Surv(trainData$fol, trainData[,dz])
  
  x_test <- model.matrix(~ . - 1, data = testData[, -c(1, 2)])  # Exclude time and status columns from predictors
  y_test <- Surv(testData$fol, testData[,dz])
  
  # Step 5: Hyperparameter Tuning using caret
  train_control <- trainControl(method = "cv", number = 10)
  model <- train(x_train, y_train, method = "glmnet", trControl = train_control)
  
  # Step 6: Addressing Class Imbalance using ROSE
  balanced_data <- ROSE(status ~ ., data = data.frame(status = trainData[,dz], trainData[, -c(1, 2)]), seed = 1)$data
  x_balanced <- model.matrix(~ . - 1, data = balanced_data[, -1])
  y_balanced <- Surv(balanced_data$time, balanced_data$status)
  
  # Step 7: Model Training with balanced data
  model_balanced <- cv.glmnet(x_balanced, y_balanced, family = "cox")
  
  # Step 8: Model Evaluation
  predictions <- predict(model_balanced, newx = x_test, s = "lambda.min", type = "response")
  
  # Example Clinical Model Predictions for Evaluation (replace with actual clinical model predictions)
  testData$clinical_model_predictions <- runif(nrow(testData), min = 0, max = 1)
  
  # Calculate Net Reclassification Improvement (NRI)
  nricens_result <- nricens(obs = testData[,dz], pred1 = testData$clinical_model_predictions, pred2 = predictions)
  
  # Calculate Integrated Discrimination Improvement (IDI)
  idinri_result <- IDI.INF(obs = testData[,dz], pred1 = testData$clinical_model_predictions, pred2 = predictions)
  
  # Step 9: Variance Partitioning
  var_part <- fitExtractVarPartModel(~ . - 1, ukbb_imputed[, c(protein_list, dz)])
  
  # Step 10: Pathway Enrichment Analysis
  query <- list(genes = protein_list)
  enrichment_results <- gost(query, organism = "hsapiens")
  
  # Calculate detection rates (from existing ukbb_analysis_module.R)
  linear_predictor <- predict(model_balanced, newx = x_test, type = "link")
  dr_results <- calculate_detection_rate(linear_predictor, testData[,dz], seq(0.01, 0.99, 0.01))
  
  # Return results
  return(list(
    model = model_balanced,
    predictions = predictions,
    detection_rates = dr_results,
    nri = nricens_result,
    idi = idinri_result,
    variance_partition = var_part,
    pathway_enrichment = enrichment_results
  ))
}