# Load necessary libraries
library(missForest)
library(glmnet)
library(caret)
library(ROSE)
library(nricens)
library(survIDINRI)
library(variancePartition)
library(gprofiler2)

# Set seed for reproducibility
set.seed(123)

# Step 1: Data Preprocessing
# Load your data
data <- read.csv("your_data_file.csv")

# Impute missing values using missForest
imputed_data <- missForest(data)$ximp

# Step 2: Feature Selection and Model Development
# Prepare the data for modeling
x <- model.matrix(~ . - 1, data = imputed_data[, -c(1, 2)])  # Exclude time and status columns from predictors
y <- Surv(imputed_data$time, imputed_data$status)

# Split the data into training and testing sets
trainIndex <- createDataPartition(imputed_data$status, p = 0.7, list = FALSE)
trainData <- imputed_data[trainIndex, ]
testData <- imputed_data[-trainIndex, ]

x_train <- model.matrix(~ . - 1, data = trainData[, -c(1, 2)])  # Exclude time and status columns from predictors
y_train <- Surv(trainData$time, trainData$status)

x_test <- model.matrix(~ . - 1, data = testData[, -c(1, 2)])  # Exclude time and status columns from predictors
y_test <- Surv(testData$time, testData$status)

# Step 3: Hyperparameter Tuning using caret
train_control <- trainControl(method = "cv", number = 10)
model <- train(x_train, y_train, method = "glmnet", trControl = train_control)

# Step 4: Addressing Class Imbalance using ROSE
balanced_data <- ROSE(status ~ ., data = trainData, seed = 1)$data
x_balanced <- model.matrix(~ . - 1, data = balanced_data[, -c(1, 2)])  # Exclude time and status columns from predictors
y_balanced <- Surv(balanced_data$time, balanced_data$status)

# Step 5: Model Training with balanced data
model_balanced <- cv.glmnet(x_balanced, y_balanced, family = "cox")

# Step 6: Model Evaluation
# Predict on the test set
predictions <- predict(model_balanced, newx = x_test, s = "lambda.min", type = "response")

# Example Clinical Model Predictions for Evaluation (replace with actual clinical model predictions)
testData$clinical_model_predictions <- runif(nrow(testData), min = 0, max = 1)

# Calculate Net Reclassification Improvement (NRI)
nricens_result <- nricens(obs = testData$status, pred1 = testData$clinical_model_predictions, pred2 = predictions)

# Calculate Integrated Discrimination Improvement (IDI)
idinri_result <- IDI.INF(obs = testData$status, pred1 = testData$clinical_model_predictions, pred2 = predictions)

# Step 7: Variance Partitioning
var_part <- fitExtractVarPartModel(~ . - 1, imputed_data)

# Step 8: Pathway Enrichment Analysis
query <- list(genes = c("protein1", "protein2", "protein3", "protein4", "protein5"))  # Example proteins
enrichment_results <- gost(query, organism = "hsapiens")

# Print results
print("Net Reclassification Improvement (NRI):")
print(nricens_result)

print("Integrated Discrimination Improvement (IDI):")
print(idinri_result)

print("Variance Partitioning Results:")
print(var_part)

print("Pathway Enrichment Analysis Results:")
print(enrichment_results)
