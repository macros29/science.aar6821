
# Create model with default paramters
control <- trainControl(method="repeatedcv", number=10, repeats=3)
seed <- 7
metric <- "Accuracy"
set.seed(seed)
mtry <- sqrt(ncol(x))
tunegrid <- expand.grid(.mtry=mtry)
rf_default <- train(x = feature_bf_fs, y = cate_bf_fs, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_default)

library(randomForest)
library(mlbench)
library(caret)

# Load Dataset
dataset <- feature_bf_fs
x <- dataset[,1:60]
y <- dataset[,61]

# Random Search
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="random")
seed = 7
set.seed(seed)
mtry <- sqrt(ncol(x))
rf_random <- train(x = feature_bf_fs, y = cate_bf_fs, method="rf", metric=metric, tuneLength = 15, trControl=control)
print(rf_random)
plot(rf_random)
