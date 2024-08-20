# Machine Learning Pipeline for Disease Prediction

This repository contains a set of R scripts that implement a machine learning pipeline for predicting diseases using protein data and other clinical information. Below is a description of the main scripts and their functions.

## data_preparation.R

This script contains functions for preparing the input data for the machine learning pipeline.

### Functions:

1. `generate_input_file()`
   - Purpose: Generates input files for the ML pipeline, including disease lists and PRS (Polygenic Risk Score) input parameters.
   - Output: Writes several text files with input parameters for different analyses.

2. `prepare_data(dz, dz.date, inc.yrs)`
   - Purpose: Prepares the main dataset for analysis by merging various data sources and applying necessary transformations.
   - Parameters:
     - `dz`: Disease name
     - `dz.date`: Disease date
     - `inc.yrs`: Incident years (5 or 10)
   - Returns: A dataframe (`ukbb`) containing the prepared data for analysis.

## feature_selection.R

This script contains functions for performing feature selection on the protein data.

### Functions:

1. `perform_feature_selection(u.train, dz, protein_list)`
   - Purpose: Performs feature selection using Lasso regression on protein data.
   - Parameters:
     - `u.train`: Training data
     - `dz`: Disease name
     - `protein_list`: List of protein names
   - Returns: A list (`las.morb`) containing the results of the Lasso regression for each subsample.

2. `generate_feature_selection_ranking(las.morb)`
   - Purpose: Generates a ranking of features based on their selection frequency in the Lasso regression.
   - Parameters:
     - `las.morb`: Output from `perform_feature_selection()`
   - Returns: A dataframe (`p.select`) with feature rankings and selection percentages.

## model_training.R

This script contains functions for training different types of prediction models.

### Functions:

1. `train_model(train_data, pred_vec, train_surv_data, times, boot_samples, test_data, test_surv_data)`
   - Purpose: Trains a Cox regression model using the specified predictors.
   - Parameters:
     - `train_data`: Training data
     - `pred_vec`: Vector of predictor names
     - `train_surv_data`: Survival data for training
     - `times`: Number of bootstrap iterations
     - `boot_samples`: Bootstrap samples
     - `test_data`: Test data
     - `test_surv_data`: Survival data for testing
   - Returns: A list containing the trained model and performance metrics.

2. `train_clinical_model(u_opti, clin_predictors, surv_opt, boot_samples, u_test, surv_test)`
   - Purpose: Trains a model using only clinical predictors.
   - Parameters: Similar to `train_model()`, but specifically for clinical predictors.
   - Returns: A trained clinical model.

3. `train_protein_model(u_opti, clin_predictors, prots_opti, surv_opt, boot_samples, u_test, surv_test)`
   - Purpose: Trains a model using both clinical predictors and selected protein features.
   - Parameters: Similar to `train_model()`, but includes both clinical and protein predictors.
   - Returns: A trained protein model.

4. `train_biomarker_model(u_opti, clin_predictors, biom_opti, surv_opt, boot_samples, u_test, surv_test)`
   - Purpose: Trains a model using clinical predictors and selected biomarkers.
   - Parameters: Similar to `train_model()`, but includes both clinical predictors and biomarkers.
   - Returns: A trained biomarker model.

## prediction_evaluation.R

This script contains functions for evaluating the performance of the trained models.

### Functions:

1. `calculate_performance_metrics(predicted, actual, thresholds)`
   - Purpose: Calculates performance metrics (false positive rate and detection rate) for different thresholds.
   - Parameters:
     - `predicted`: Predicted probabilities
     - `actual`: Actual outcomes
     - `thresholds`: Vector of threshold values
   - Returns: A dataframe with threshold, false positive rate, and detection rate for each threshold.

2. `generate_dr_curve(res_clin, res_clin_prots, u_test, dz)`
   - Purpose: Generates detection rate curves for clinical and protein models.
   - Parameters:
     - `res_clin`: Results from clinical model
     - `res_clin_prots`: Results from protein model
     - `u_test`: Test data
     - `dz`: Disease name
   - Returns: A dataframe with detection rates for different false positive rates.

3. `normalize_linear_predictor(lp)`
   - Purpose: Normalizes the linear predictor to a 0-1 scale.
   - Parameters:
     - `lp`: Linear predictor values
   - Returns: Normalized linear predictor values.

## utility_functions.R

This script contains utility functions used across the pipeline.

### Functions:

1. `get_clinical_predictors(dz)`
   - Purpose: Returns a list of clinical predictors based on the disease.
   - Parameters:
     - `dz`: Disease name
   - Returns: A vector of clinical predictor names.

2. `split_data(ukbb, dz, n.cases)`
   - Purpose: Splits the data into training, optimization, and test sets.
   - Parameters:
     - `ukbb`: Full dataset
     - `dz`: Disease name
     - `n.cases`: Number of cases
   - Returns: A list containing train, optimization, and test datasets.

3. `create_survival_objects(u_train, u_opti, u_test, dz)`
   - Purpose: Creates survival objects for training, optimization, and test sets.
   - Parameters:
     - `u_train`: Training data
     - `u_opti`: Optimization data
     - `u_test`: Test data
     - `dz`: Disease name
   - Returns: A list of survival objects for each dataset.

Each of these scripts plays a crucial role in the overall machine learning pipeline for disease prediction. They handle data preparation, feature selection, model training, prediction evaluation, and provide utility functions to support the entire process.
