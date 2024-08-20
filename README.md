# 📚 Machine Learning Pipeline for Disease Prediction

Welcome to the documentation for the Machine Learning Pipeline for Disease Prediction. This repository contains a set of R scripts that implement a machine learning pipeline for predicting diseases using protein data and other clinical information. The main modules used in this pipeline are tidyverse, survival, glmnet, and caret.

## 📑 Index

1. [data_preparation.R](#data_preparationr)
2. [feature_selection.R](#feature_selectionr)
3. [model_training.R](#model_trainingr)
4. [prediction_evaluation.R](#prediction_evaluationr)
5. [utility_functions.R](#utility_functionsr)

---

## 1. data_preparation.R

This script contains functions for preparing the input data for the machine learning pipeline.

### 1.1 `generate_input_file()`

#### Description

Generates input files for the ML pipeline, including disease lists and PRS (Polygenic Risk Score) input parameters.

#### Parameters

This function does not take any parameters.

#### Returns

| Type          | Description                                                      |
|---------------|------------------------------------------------------------------|
| `None`        | The function writes several text files with input parameters for different analyses. |

#### Example Usage

```r
generate_input_file()
```

### 1.2 `prepare_data(dz, dz.date, inc.yrs)`

#### Description

Prepares the main dataset for analysis by merging various data sources and applying necessary transformations.

#### Parameters

| Parameter | Type    | Description                   |
|-----------|---------|-------------------------------|
| `dz`      | string  | Disease name                  |
| `dz.date` | string  | Disease date                  |
| `inc.yrs` | integer | Incident years (5 or 10)      |

#### Returns

| Type        | Description                                            |
|-------------|--------------------------------------------------------|
| `data.frame`| A dataframe (`ukbb`) containing the prepared data for analysis. |

#### Example Usage

```r
ukbb_data <- prepare_data("diabetes", "2022-01-01", 10)
```

---

## 2. feature_selection.R

This script contains functions for performing feature selection on the protein data.

### 2.1 `perform_feature_selection(u.train, dz, protein_list)`

#### Description

Performs feature selection using Lasso regression on protein data.

#### Parameters

| Parameter     | Type        | Description                   |
|---------------|-------------|-------------------------------|
| `u.train`     | data.frame  | Training data                 |
| `dz`          | string      | Disease name                  |
| `protein_list`| vector      | List of protein names         |

#### Returns

| Type  | Description                                                                    |
|-------|--------------------------------------------------------------------------------|
| `list`| A list (`las.morb`) containing the results of the Lasso regression for each subsample. |

#### Example Usage

```r
las_morb <- perform_feature_selection(training_data, "diabetes", protein_names)
```

### 2.2 `generate_feature_selection_ranking(las.morb)`

#### Description

Generates a ranking of features based on their selection frequency in the Lasso regression.

#### Parameters

| Parameter | Type  | Description                                  |
|-----------|-------|----------------------------------------------|
| `las.morb`| list  | Output from `perform_feature_selection()`    |

#### Returns

| Type        | Description                                                          |
|-------------|----------------------------------------------------------------------|
| `data.frame`| A dataframe (`p.select`) with feature rankings and selection percentages. |

#### Example Usage

```r
feature_rankings <- generate_feature_selection_ranking(las_morb)
```

---

## 3. model_training.R

This script contains functions for training different types of prediction models.

### 3.1 `train_model(train_data, pred_vec, train_surv_data, times, boot_samples, test_data, test_surv_data)`

#### Description

Trains a Cox regression model using the specified predictors.

#### Parameters

| Parameter        | Type        | Description                   |
|------------------|-------------|-------------------------------|
| `train_data`     | data.frame  | Training data                 |
| `pred_vec`       | vector      | Vector of predictor names     |
| `train_surv_data`| Surv object | Survival data for training    |
| `times`          | integer     | Number of bootstrap iterations|
| `boot_samples`   | list        | Bootstrap samples             |
| `test_data`      | data.frame  | Test data                     |
| `test_surv_data` | Surv object | Survival data for testing     |

#### Returns

| Type  | Description                                                    |
|-------|----------------------------------------------------------------|
| `list`| A list containing the trained model and performance metrics.   |

#### Example Usage

```r
model_results <- train_model(train_data, predictors, train_surv, 1000, boot_samples, test_data, test_surv)
```

### 3.2 `train_clinical_model(u_opti, clin_predictors, surv_opt, boot_samples, u_test, surv_test)`

#### Description

Trains a model using only clinical predictors.

#### Parameters

| Parameter        | Type        | Description                   |
|------------------|-------------|-------------------------------|
| `u_opti`         | data.frame  | Optimization data             |
| `clin_predictors`| vector      | Clinical predictor names      |
| `surv_opt`       | Surv object | Survival data for optimization|
| `boot_samples`   | list        | Bootstrap samples             |
| `u_test`         | data.frame  | Test data                     |
| `surv_test`      | Surv object | Survival data for testing     |

#### Returns

| Type  | Description                |
|-------|----------------------------|
| `list`| A trained clinical model.  |

#### Example Usage

```r
clinical_model <- train_clinical_model(optimization_data, clinical_predictors, surv_opt, boot_samples, test_data, surv_test)
```

### 3.3 `train_protein_model(u_opti, clin_predictors, prots_opti, surv_opt, boot_samples, u_test, surv_test)`

#### Description

Trains a model using both clinical predictors and selected protein features.

#### Parameters

| Parameter        | Type        | Description                   |
|------------------|-------------|-------------------------------|
| `u_opti`         | data.frame  | Optimization data             |
| `clin_predictors`| vector      | Clinical predictor names      |
| `prots_opti`     | vector      | Optimized protein features    |
| `surv_opt`       | Surv object | Survival data for optimization|
| `boot_samples`   | list        | Bootstrap samples             |
| `u_test`         | data.frame  | Test data                     |
| `surv_test`      | Surv object | Survival data for testing     |

#### Returns

| Type  | Description                |
|-------|----------------------------|
| `list`| A trained protein model.   |

#### Example Usage

```r
protein_model <- train_protein_model(optimization_data, clinical_predictors, protein_features, surv_opt, boot_samples, test_data, surv_test)
```

### 3.4 `train_biomarker_model(u_opti, clin_predictors, biom_opti, surv_opt, boot_samples, u_test, surv_test)`

#### Description

Trains a model using clinical predictors and selected biomarkers.

#### Parameters

| Parameter        | Type        | Description                   |
|------------------|-------------|-------------------------------|
| `u_opti`         | data.frame  | Optimization data             |
| `clin_predictors`| vector      | Clinical predictor names      |
| `biom_opti`      | vector      | Optimized biomarker features  |
| `surv_opt`       | Surv object | Survival data for optimization|
| `boot_samples`   | list        | Bootstrap samples             |
| `u_test`         | data.frame  | Test data                     |
| `surv_test`      | Surv object | Survival data for testing     |

#### Returns

| Type  | Description                |
|-------|----------------------------|
| `list`| A trained biomarker model. |

#### Example Usage

```r
biomarker_model <- train_biomarker_model(optimization_data, clinical_predictors, biomarker_features, surv_opt, boot_samples, test_data, surv_test)
```

---

## 4. prediction_evaluation.R

This script contains functions for evaluating the performance of the trained models.

### 4.1 `calculate_performance_metrics(predicted, actual, thresholds)`

#### Description

Calculates performance metrics (false positive rate and detection rate) for different thresholds.

#### Parameters

| Parameter   | Type    | Description                   |
|-------------|---------|-------------------------------|
| `predicted` | vector  | Predicted probabilities       |
| `actual`    | vector  | Actual outcomes               |
| `thresholds`| vector  | Vector of threshold values    |

#### Returns

| Type        | Description                                                                           |
|-------------|---------------------------------------------------------------------------------------|
| `data.frame`| A dataframe with threshold, false positive rate, and detection rate for each threshold.|

#### Example Usage

```r
metrics <- calculate_performance_metrics(predicted_probs, actual_outcomes, seq(0, 1, by = 0.1))
```

### 4.2 `generate_dr_curve(res_clin, res_clin_prots, u_test, dz)`

#### Description

Generates detection rate curves for clinical and protein models.

#### Parameters

| Parameter       | Type        | Description                   |
|-----------------|-------------|-------------------------------|
| `res_clin`      | list        | Results from clinical model   |
| `res_clin_prots`| list        | Results from protein model    |
| `u_test`        | data.frame  | Test data                     |
| `dz`            | string      | Disease name                  |

#### Returns

| Type        | Description                                                            |
|-------------|------------------------------------------------------------------------|
| `data.frame`| A dataframe with detection rates for different false positive rates.   |

#### Example Usage

```r
dr_curve <- generate_dr_curve(clinical_results, protein_results, test_data, "diabetes")
```

### 4.3 `normalize_linear_predictor(lp)`

#### Description

Normalizes the linear predictor to a 0-1 scale.

#### Parameters

| Parameter | Type    | Description           |
|-----------|---------|-----------------------|
| `lp`      | vector  | Linear predictor values |

#### Returns

| Type    | Description                           |
|---------|---------------------------------------|
| `vector`| Normalized linear predictor values.   |

#### Example Usage

```r
normalized_lp <- normalize_linear_predictor(linear_predictor)
```

---

## 5. utility_functions.R

This script contains utility functions used across the pipeline.

### 5.1 `get_clinical_predictors(dz)`

#### Description

Returns a list of clinical predictors based on the disease.

#### Parameters

| Parameter | Type    | Description |
|-----------|---------|-------------|
| `dz`      | string  | Disease name|

#### Returns

| Type    | Description                           |
|---------|---------------------------------------|
| `vector`| A vector of clinical predictor names. |

#### Example Usage

```r
clinical_predictors <- get_clinical_predictors("diabetes")
```

### 5.2 `split_data(ukbb, dz, n.cases)`

#### Description

Splits the data into training, optimization, and test sets.

#### Parameters

| Parameter | Type        | Description       |
|-----------|-------------|-------------------|
| `ukbb`    | data.frame  | Full dataset      |
| `dz`      | string      | Disease name      |
| `n.cases` | integer     | Number of cases   |

#### Returns

| Type  | Description                                                   |
|-------|---------------------------------------------------------------|
| `list`| A list containing train, optimization, and test datasets.     |

#### Example Usage

```r
data_splits <- split_data(ukbb_data, "diabetes", 1000)
```

### 5.3 `create_survival_objects(u_train, u_opti, u_test, dz)`

#### Description

Creates survival objects for training, optimization, and test sets.

#### Parameters

| Parameter | Type        | Description     |
|-----------|-------------|-----------------|
| `u_train` | data.frame  | Training data   |
| `u_opti`  | data.frame  | Optimization data|
| `u_test`  | data.frame  | Test data       |
| `dz`      | string      | Disease name    |

#### Returns

| Type  | Description                                           |
|-------|-------------------------------------------------------|
| `list`| A list of survival objects for each dataset.          |

#### Example Usage

```r
surv_objects <- create_survival_objects(train_data, opti_data, test_data, "diabetes")
```

---

Each of these scripts plays a crucial role in the overall machine learning pipeline for disease prediction. They handle data preparation, feature selection, model training, prediction evaluation, and provide utility functions to support the entire process. This documentation should provide a comprehensive understanding of the capabilities and usage of each function within the pipeline. Happy coding! 🎉
