# Kan svampen spises?
# In this example, we are aiming to predict whether a mushroom can be eaten or not
# Mushroom data is cited from UCI Machine Learning Repository. Bache and Lichman (2013).

# setup ------------------------------------------------------------------
source(here::here("R/set_up.R"))
source(here::here("R/package_load.R"))

set.seed(100)


# Loading data ------------------------------------------------------------

data(agaricus.train, package='xgboost')
data(agaricus.test, package='xgboost')
train <- agaricus.train
test <- agaricus.test

# i fremtiden skal man selv splitte data op

## structure of the data -------------------------------------------------
# dimensions
dim(train$data)
dim(test$data)

# class
class(train$data)[1]
class(train$label)


# 1.3.3.1 basic training --------------------------------------------------

# train decision tree model
bstSparse <- xgboost(data = train$data, 
                     label = train$label, 
                     max_depth = 2, # the trees won’t be deep, because our case is very simple 
                     eta = 1, 
                     nthread = 2, # the number of CPU threads we are going to use
                     nrounds = 2, #there will be two passes on the data, the second one will enhance the model by further reducing the difference between ground truth and prediction.
                     objective = "binary:logistic") # we will train a binary classification model

# More complex the relationship between your features and your label is, more passes you need.



# 1.3.3.2 Parameter variations --------------------------------------------

## 1.3.3.2.1 Dense matrix --------------------------------------------------

bstDense <- xgboost(
  data = as.matrix(train$data),
  label = train$label,
  max_depth = 2,
  eta = 1,
  nthread = 2,
  nrounds = 2,
  objective = "binary:logistic"
)


## 1.3.3.2.2 xgb.DMatrix ---------------------------------------------------

dtrain <- xgb.DMatrix(data = train$data, label = train$label, nthread = 2)
bstDMatrix <- xgboost(
  data = dtrain,
  max_depth = 2,
  eta = 1,
  nthread = 2,
  nrounds = 2,
  objective = "binary:logistic"
)


## 1.3.3.2.3 Verbose option ------------------------------------------------
# verbose = 0, no message
bst <- xgboost(data = dtrain, 
               max_depth = 2, 
               eta = 1, 
               nthread = 2, 
               nrounds = 2, 
               objective = "binary:logistic", 
               verbose = 1)

# verbose = 1, print evaluation metric
bst <- xgboost(data = dtrain, 
               max_depth = 2, 
               eta = 1, 
               nthread = 2, 
               nrounds = 2, 
               objective = "binary:logistic", 
               verbose = 1)
# verbose = 2, also print information about tree
bst <- xgboost(data = dtrain, 
               max_depth = 2, 
               eta = 1, 
               nthread = 2, 
               nrounds = 200, 
               objective = "binary:logistic", 
               verbose = 2)


# 1.4 Basic prediction using XGBoost --------------------------------------

# 1.5 Perform the prediction ----------------------------------------------

pred <- predict(bst, test$data)

# size of the prediction vector
print(length(pred))

# limit display of predictions to the first 10
print(head(pred))


# 1.6 Transform the regression in a binary classification -----------------
prediction <- as.numeric(pred > 0.5)
print(head(prediction))


# 1.7 Measuring model performance -----------------------------------------
# the average error

# den beregner procent andel af labels som den har lavet forkert ift. de rigtige labels
err <- mean(as.numeric(pred > 0.5) != test$label)
print(paste("test-error=", err))

length(which(as.numeric(pred > 0.5) != test$label))

(35/1611)*100

# 1.8 Advanced features ---------------------------------------------------

# Most of the features below have been implemented to help you to improve 
# your model by offering a better understanding of its content


# 1.8.1 Dataset preparation -----------------------------------------------

dtrain <- xgb.DMatrix(data = train$data, label = train$label, nthread = 2)
dtest <- xgb.DMatrix(data = test$data, label = test$label, nthread = 2)


# 1.8.2 Measure learning progress with xgb.train --------------------------

watchlist <- list(train=dtrain, test=dtest)

ggplot(bst$evaluation_log, aes(x = iter, y = train_logloss)) + 
  geom_line(color = "red") + geom_line(aes(y = test_logloss), color = "blue", alpha = 0.5) +
  scale_y_reverse()

bst <- xgb.train(data=dtrain, 
                 max_depth=2, 
                 eta=1, 
                 nthread = 2, 
                 nrounds=7, 
                 watchlist=watchlist, 
                 objective = "binary:logistic")

bst <- xgb.train(data=dtrain,
                 max_depth=2, 
                 eta=1, 
                 nthread = 2, 
                 nrounds=2, 
                 watchlist=watchlist, 
                 eval_metric = "error", 
                 eval_metric = "logloss", 
                 objective = "binary:logistic")


# 1.8.3 Linear boosting ---------------------------------------------------

bst <- xgb.train(data=dtrain, 
                 booster = "gblinear", 
                 max_depth=2, 
                 nthread = 2, 
                 nrounds=2, 
                 watchlist=watchlist, 
                 eval_metric = "error", 
                 eval_metric = "logloss", 
                 objective = "binary:logistic")


# 1.8.4 Manipulating xgb.DMatrix ------------------------------------------


## 1.8.4.1 Save / Load -----------------------------------------------------

xgb.DMatrix.save(dtrain, "data/xgboost/dtrain.buffer")

# to load it in, simply call xgb.DMatrix
dtrain2 <- xgb.DMatrix("data/xgboost/dtrain.buffer")

bst <- xgb.train(data=dtrain2, 
                 max_depth=2, 
                 eta=1, 
                 nthread = 2, 
                 nrounds=2, 
                 watchlist=watchlist, 
                 objective = "binary:logistic")


## 1.8.4.2 Information extraction ------------------------------------------

label = getinfo(dtest, "label")
pred <- predict(bst, dtest)
err <- as.numeric(sum(as.integer(pred > 0.5) != label))/length(label)
print(paste("test-error=", err))


# 1.8.5 View feature importance/influence from the learnt model -----------


## 1.8.5.1 View the trees from a model -------------------------------------

xgb.dump(bst, with_stats = TRUE)
xgb.plot.tree(model = bst)


## 1.8.5.2 Save and load models --------------------------------------------

# save model to binary local file
xgb.save(bst, "data/xgboost/xgboost.model")

# load binary model to R
bst2 <- xgb.load("data/xgboost/xgboost.model")
xgb.parameters(bst2) <- list(nthread = 2)
pred2 <- predict(bst2, test$data)

# And now the test
print(paste("sum(abs(pred2-pred))=", sum(abs(pred2-pred))))

# save model to R's raw vector
rawVec <- xgb.serialize(bst)

# print class
print(class(rawVec))

# load binary model to R
bst3 <- xgb.load(rawVec)
xgb.parameters(bst3) <- list(nthread = 2)
pred3 <- predict(bst3, test$data)

# pred2 should be identical to pred
print(paste("sum(abs(pred3-pred))=", sum(abs(pred2-pred))))
