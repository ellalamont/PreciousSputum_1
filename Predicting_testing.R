# First attempt on predicting outcome
# E. Lamont
# 8/12/25


# https://www.sthda.com/english/articles/36-classification-methods-essentials/149-penalized-logistic-regression-essentials-in-r-ridge-lasso-and-elastic-net/



source("Import_data.R") # sputumOnly_pipeSummary

# install.packages("caret")
library(caret)
# install.packages("glmnet")
library(glmnet)


###########################################################
##################### ORGANIZED DATA ######################

# Start with just the good sputum I guess, and remove the extra W0 and W2 that I did (13004 and 14005)
P_pipeSummary <- All_pipeSummary %>% filter(SampleID %in% SputumSampleList) %>% filter(!SampleID %in% c("W0_13004_S53", "W0_14005_S55"))

P_metadata <- P_pipeSummary %>% select(SampleID, Outcome)
P_metadata$Outcome <- factor(P_metadata$Outcome)


P_TPM <- Run1_tpm %>% select("X", all_of(SputumSampleList), -c("W0_13004_S53", "W0_14005_S55"))

W0sputumOnly_pipeSummary <- sputumOnly_pipeSummary %>% filter(!str_detect(SampleID, "W2"))
W0sputumOnly_TPM <- Run1_tpm %>% select(contains(("W0"))) %>% select(!contains("RT"))

# Put everything in one dataframe
P_TPM_t <- P_TPM %>% 
  column_to_rownames("X") %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID")
my_df <- inner_join(P_metadata, P_TPM_t, by = "SampleID")




###########################################################
############ CREATE TRAINING AND TESTING SETS #############
# Want at least 2 relapses in the testing set
# do a 2/3 split?

set.seed(42)
trainIndex <- createDataPartition(P_metadata$Outcome, p = 2/3, list = F)
train <- P_metadata[trainIndex,]
test <- P_metadata[-trainIndex,]
# Checked and this does actually split well, but probably want to do a more complicated re-sampling thing


###########################################################
###################### NEW WAY - LASSO ####################

# Separate the training and testing datasets
set.seed(42)
trainIndex <- createDataPartition(my_df$Outcome, p = 2/3, list = F)
train <- my_df[trainIndex,]
test <- my_df[-trainIndex,]
# Checked and there are 2 relapses in the test partition


# Create matrices for glmnet
x_train <- as.matrix(train[, !names(train) %in% c("SampleID", "Outcome")]) # matrix of predictor variables
# y_train <- train$Outcome # the response or outcome variable, which is a binary variable.
y_train <- ifelse(train$Outcome == "Relapse", 1, 0)
x_test <- as.matrix(test[, !names(test) %in% c("SampleID", "Outcome")]) 
# y_test <- test$Outcome 
y_test <- ifelse(test$Outcome == "Relapse", 1, 0)

# Stratify cross-validation inside glmnet??
# set.seed(42)
# folds <- createFolds(y_train, k = 5, list = F)


# https://www.sthda.com/english/articles/36-classification-methods-essentials/149-penalized-logistic-regression-essentials-in-r-ridge-lasso-and-elastic-net/

# Find the best lambda using cross-validation
# lamba: a numeric value defining the amount of shrinkage. Should be specify by analyst.
set.seed(42)
cv.lasso <- cv.glmnet(x_train, y_train, alpha = 1, family = "binomial")
# Warning messages:
# 1: In lognet(x, is.sparse, y, weights, offset, alpha, nobs,  ... :
#                one multinomial or binomial class has fewer than 8  observations; dangerous ground
plot(cv.lasso)
cv.lasso$lambda.min # 0.01959793
cv.lasso$lambda.1se # 0.2201483
# Not sure which of the above should be used....

# coef(cv.lasso, cv.lasso$lambda.min) # There are too many to see anything with this....


# Fit the final model on the training data
lasso.model <- glmnet(x_train, y_train, alpha = 1, family = "binomial", lambda = cv.lasso$lambda.min)

# Make prediction on test data
probabilities <- lasso.model %>% predict(newx = x_test)
predicted.classes <- ifelse(probabilities > 0.5, "relapse", "cure")

# Model accuracy rate
observed.classes <- test$Outcome
mean(predicted.classes == observed.classes) # 0, so not good!


###########################################################
#################### COMPUTE FULL MODEL ###################
# https://www.sthda.com/english/articles/36-classification-methods-essentials/149-penalized-logistic-regression-essentials-in-r-ridge-lasso-and-elastic-net/

# Fit the model
full.model <- glm(Outcome ~., data = x_train, family = binomial)
# Make predictions
probabilities <- full.model %>% predict(test.data, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, "pos", "neg")
# Model accuracy
observed.classes <- test.data$diabetes
mean(predicted.classes == observed.classes)



