
#Libraries
library(partDSA)
library(survival)
library(RSNNS)
library(Rcpp)
library(caret)
library(Matrix)
library(plyr)
library(dplyr)
library(doParallel)

# Data Preparation --------------------------------------------------------

paper2_data <- read.csv("H:/tklarkin_backup/Documents/My Research/Dissertation/Paper 2/Code/master_dataset_paper2.csv")
dim(paper2_data)
str(paper2_data)
prop.table(table(paper2_data$Near_Earth))

response <- paper2_data$Min_DST

predictors <- paper2_data[,28:ncol(paper2_data)]
colnames(predictors) <- c("MPA", "AW", "LS", "SOI", "SOF",
                          "SOR", "Acc", "Poor", "Very_Poor",
                          "RFlux", "SSN", "SSA", "NR", "XrayC", "XrayM", "XrayX")

metadata <- read.csv(file ="H:/tklarkin_backup/Documents/My Research/Journal Articles/Venice Conference Extended/Code/stage1_predictions.csv")
metadata <- metadata[,-1] #removing obs number
dim(metadata)

# Prior to Modeling -------------------------------------------------------

my.rqnc <- list(label = "Non-Convex Penalized Quantile Regression with SCAD Penalty",
                  library = "rqPen",
                  type = "Regression",
                  parameters = data.frame(parameter = c('lambda', 'tau'),
                                          class = c("numeric", "numeric"),
                                          label =  c('L1 Penalty', 'Quantile Level')),
                  grid = function(x, y, len = NULL, search = "grid"){
                    if(search == "grid") {
                      out <- expand.grid(lambda = c(1000, 1, 0.001),
                                         tau = c(0.1, 0.2, 0.3))
                    } else {
                      out <- data.frame(lambda = 10^runif(len, min = -5, 1),
                                        tau = 10^runif(len, min = -5, 1))
                    }
                    out
                  },
                  fit = function(x, y, wts, param, lev, last, classProbs, ...) {
                    model <- QICD(x = as.matrix(x), y = y, tau = param$tau, lambda = param$lambda, 
                              penalty = "SCAD", ...)
                    var.names <- colnames(x)
                    out <- list(fit = model, names = var.names)
                  },
                  predict = function(modelFit, newdata, submodels = NULL) {
                    predict(modelFit[["fit"]], newx = as.matrix(newdata))[,1]
                  },
                  predictors = function(x, ...) {
                    out <- coef(x[["fit"]])
                    out <- out[names(out) != "intercept"]
                    names(out)[out != 0]
                 },
                  varImp = function(object, ...) {
                    values <- object[["fit"]]$coef[-1] #taking intercept out
                    names(values) <- object[["names"]]
                    varImps <-abs(values) 
                    out <- data.frame(varImps)
                    colnames(out) <- "Overall"
                    if(!is.null(names(varImps))) rownames(out) <- names(varImps)
                    out   
                 },
                  tags = c("Linear Regression", "Quantile Regression", "Implicit Feature Selection", 
                           "L1 Regularization"),
                  prob = NULL,
                  sort = function(x) x[order(-x$lambda),])


#Weighted MAE
wmaeSummary <- function (data, lev = NULL, model = NULL) {
  w <- ifelse(data$obs <= -100, 10, 1) #strong storms weighted 10 times higher
  out <- (1/sum(w))*sum(w*(abs(data$pred-data$obs)))  
  names(out) <- "WMAE"
  out
}

#nested cross validation
nestedcv <- function(x, y, method, metric, trControl, tuneLength, nestedrepeats, verbose, ...) {
  results <- NULL
  # iterate over 
  for(repeatnum in 1:nestedrepeats) {
    # split into folds
    if(verbose!= 0) {
      cat("Processing repeat iteration", repeatnum, "/", nestedrepeats, "\n")
    }
    set.seed(repeatnum)
    cv = createFolds(y, k = trControl$number)
    # iterate over folds
    for(iter in 1:length(cv)) {
      # test indices
      if(verbose!= 0) {
        cat("Processing outer loop iteration", iter, "/", length(cv), "\n")
      }
      test_indices = cv[[iter]]
      train_indices = setdiff(seq(length(y)), test_indices)
      # train model with caret package
      if(method == "nnls"){ #produces error unless its in formula form (dunno why)
        dat <- cbind.data.frame(y = y[train_indices], x[train_indices,])
        model = train(y ~ ., data = dat, method=method,
                      metric=metric, trControl=trControl, tuneLength=tuneLength, ...)
      } else{
        model = train(x[train_indices,], y[train_indices], method=method,
                      metric=metric, trControl=trControl, tuneLength=tuneLength, ...)
      }  
      # do predictions
      predictions = predict(model, x[test_indices,])
      # format predictions into form for summary Function
      tmp = data.frame(predictions, y[test_indices], stringsAsFactors=FALSE)
      colnames(tmp) <- c("pred", "obs")
      tmp.100 = tmp[which(tmp$obs <= -100),]
      thisResample = data.frame(WMAE = wmaeSummary(data = tmp), row.names = NULL)
      thisResample$RMSE = postResample(pred = tmp$pred, obs = tmp$obs)["RMSE"]
      thisResample$RMSE.100 = postResample(pred = tmp.100$pred, obs = tmp.100$obs)["RMSE"]
      thisResample$Resample = iter
      thisResample$Imp.Vars = nnzero(rowSums(varImp(model, scale = FALSE)$importance))
      thisResample$ntrain = length(train_indices)
      thisResample$ntest = length(test_indices)
      thisResample$Missing = sum(is.na(thisResample))
      thisResample$Time = model$times$everything[3]
      row.names(thisResample) <- NULL
      # append this result to a list
      show(thisResample)
      results = rbind(results, thisResample)
    } # end for(iter in 1:length(cv))
  } # end for(iter in 1:nestedrepeats)
  #training model on all data
  if(method == "nnls"){
    dat <- cbind.data.frame(y = y, x)
    model = train(y ~ ., data = dat, method=method,
                  metric=metric, trControl=trControl, tuneLength=tuneLength, ...)
  } else{
    model = train(x, y, method=method,
                  metric=metric, trControl=trControl, tuneLength=tuneLength, ...)
  }
  return(list(results = results, fit = model))
}

#function based on formula from Bouckaert and Frank (2004)
#comparing RMSE
corrected.rep.kfold.cv.test <- function(a, b, metric) {
  #Assuming a and b must be trained over the same folds (which they should)
  #To account for "near-equal" folds, taking mean over all observations for train/test of last iteration
  n1 <- mean(a$results$ntrain)
  n2 <- mean(a$results$ntest)
  xij <- a$results[,metric] - b$results[,metric]
  #Number of folds times number of repeats
  kr <- length(xij)
  #Estimate for mean
  m <- mean(xij)
  #Estimate for variance with variance correction
  sigmasq <- (1/(kr-1))*(sum((xij - m)^2))
  #Calculating t-statistic
  t <- m/sqrt((1/kr+n2/n1)*sigmasq)
  #Computing p-value
  t.pval <- 2*pt(-abs(t), df = kr-1)
  data.frame("P-value for Corrected Repeated k-fold CV Test" = t.pval)
}

#Setting tune length and metric
tl <- 3
metric <- "WMAE"

#Running in parallel
cl <- makeCluster(4)
registerDoParallel(cl)

# Start here and run for reproducible results -----------------------------

set.seed(10)

#creating a list of seeds for each resample (large enough for both 1x10 and 10x10)
seeds <- vector(mode = "list", length = 101) #length is = (n_repeats*nresampling)+1
for(i in 1:100) seeds[[i]]<- sample.int(n=1000, 36) #large enough for all par combos
seeds[[101]]<-sample.int(1000, 1)#for the last model


fitControl <- trainControl(method = "cv", number = 10, summaryFunction = wmaeSummary,
                           seeds = seeds, allowParallel = TRUE, 
                           savePredictions = "final", search = "grid")

# Implementing Meta-learners --------------------------------------------------

meta.lm <- nestedcv(y = response, x = metadata, method = "lm", tuneLength = tl,
                 metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Quantile Regression with SCAD Penalty with QICD algorithm
metalearner <- nestedcv(y = response, x = metadata, method = my.rqnc, tuneLength = tl,
                        metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

# Training Base-learners Fully for Comparison ---------------------------------

#Bagged CART
base.1 <- nestedcv(x = predictors, y = response, method = "treebag", tuneLength = tl, 
                        metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Bayesian Lasso
base.2 <- nestedcv(x = predictors, y = response, method = "blasso", tuneLength = tl, 
                        metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Bayesian Lasso (Model Averaged)
base.3 <- nestedcv(x = predictors, y = response, method = "blassoAveraged", tuneLength = tl, 
                        metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Bayesian Regularized Neural Network
base.4 <- nestedcv(x = predictors, y = response, method = "brnn", tuneLength = tl, 
                        metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Bayesian Ridge Regression
base.5 <- nestedcv(x = predictors, y = response, method = "bridge", tuneLength = tl,
                        metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Boosted Linear Model
base.6 <- nestedcv(x = predictors, y = response, method = "glmboost", tuneLength = tl, 
                        metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Boosted Tree
base.7 <- nestedcv(x = predictors, y = response, method = "bstTree", tuneLength = tl, 
                        metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Regression Tree
base.8 <- nestedcv(x = predictors, y = response, method = "rpart1SE", tuneLength = tl,
                        metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Conditional Inference Random Forest
base.9 <- nestedcv(x = predictors, y = response, method = "cforest", tuneLength = tl, 
                        metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Conditional Inference Tree
base.10 <- nestedcv(x = predictors, y = response, method = "ctree", tuneLength = tl, 
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Cubist
base.11 <- nestedcv(x = predictors, y = response, method = "cubist", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#eXtreme Gradient Boosting with Linear Booster
base.12 <- nestedcv(x = predictors, y = response, method = "xgbLinear", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#eXtreme Gradient Boosting with Tree Booster
base.13 <- nestedcv(x = predictors, y = response, method = "xgbTree", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Extreme Learning Machine
base.14 <- nestedcv(x = predictors, y = response, method = "elm", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Generalized Additive Model using Splines
base.15 <- nestedcv(x = predictors, y = response, method = "gamSpline", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Independent Component Regression
base.16 <- nestedcv(x = predictors, y = response, method = "icr", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#k-Nearest Neighbors Regression
base.17 <- nestedcv(x = predictors, y = response, method = "kknn", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Ridge, Lasso, and Elastic Net
base.18 <- nestedcv(x = predictors, y = response, method = "glmnet", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Least Angle Regression
base.19 <- nestedcv(x = predictors, y = response, method = "lars2", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Linear Regression
base.20 <- nestedcv(x = predictors, y = response, method = "lm", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Linear Regression with Stepwise Selection
base.21 <- nestedcv(x = predictors, y = response, method = "leapSeq", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Generalized Additive Model using Loess
base.22 <- nestedcv(x = predictors, y = response, method = "gamLoess", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Multi-layer Perceptron
base.23 <- nestedcv(x = predictors, y = response, method = "mlp", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Multivariate Adaptive Regression Splines
base.24 <- nestedcv(x = predictors, y = response, method = "earth", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Bagged Multivariate Adaptive Regression Splines with GCV Pruning
base.25 <- nestedcv(x = predictors, y = response, method = "bagEarthGCV", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Neural Network
base.26 <- nestedcv(x = predictors, y = response, method = "nnet", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, linout = TRUE, nestedrepeats = 10, verbose = 1)

#Neural Network with Feature Extraction
base.27 <- nestedcv(x = predictors, y = response, method = "pcaNNet", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, linout = TRUE, nestedrepeats = 10, verbose = 1)

#Non-convex Penalized Quantile Regression
base.28 <- nestedcv(x = predictors, y = response, method = "rqnc", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1, tau = 0.2)

#Non-negative Least Squares
base.29 <- nestedcv(x = predictors, y = response, method = "nnls", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Partitioning Using Deletion, Substitution, and Addition Moves
base.30 <- nestedcv(x = predictors, y = response, method = "partDSA", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Partial Least Squares
base.31 <- nestedcv(x = predictors, y = response, method = "pls", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Principal Component Regression
base.32 <- nestedcv(x = predictors, y = response, method = "pcr", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Projection Pursuit Regression
base.33 <- nestedcv(x = predictors, y = response, method = "ppr", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Quantile Random Forest
base.34 <- nestedcv(x = predictors, y = response, method = "qrf", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Quantile Regression with Lasso Penalty
base.35 <- nestedcv(x = predictors, y = response, method = "rqlasso", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1, tau = 0.2)

#Random Forest
base.36 <- nestedcv(x = predictors, y = response, method = "ranger", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, importance = "permutation", nestedrepeats = 10, verbose = 1) 

#Relaxed Lasso
base.37 <- nestedcv(x = predictors, y = response, method = "relaxo", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Ridge Regression with Variable Selection
base.38 <- nestedcv(x = predictors, y = response, method = "foba", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Robust Linear Model
base.39 <- nestedcv(x = predictors, y = response, method = "rlm", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Self-Organizing Map
base.40 <- nestedcv(x = predictors, y = response, method = "xyf", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Spike and Slab Regression
base.41 <- nestedcv(x = predictors, y = response, method = "spikeslab", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Stacked AutoEncoder Deep Neural Network
base.42 <- nestedcv(x = predictors, y = response, method = "dnn", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Stochastic Gradient Boosting
base.43 <- nestedcv(x = predictors, y = response, method = "gbm", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Supervised Principal Component Analysis
base.44 <- nestedcv(x = predictors, y = response, method = "superpc", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#SVM with Linear Kernel
base.45 <- nestedcv(x = predictors, y = response, method = "svmLinear", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#SVM with Polynomial Kernel
base.46 <- nestedcv(x = predictors, y = response, method = "svmPoly", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#SVM with Radial Basis Function Kernel
base.47 <- nestedcv(x = predictors, y = response, method = "svmRadialSigma", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Regularized Random Forest
base.48 <- nestedcv(x = predictors, y = response, method = "RRFglobal", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, importance = TRUE, nestedrepeats = 10, verbose = 1)

#Multilayer Perceptron Network by Stochastic Gradient Descent
base.49 <- nestedcv(x = predictors, y = response, method = "mlpSGD", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

#Weighted k-Nearest Neighbors
base.50 <- nestedcv(x = predictors, y = response, method = "knn", tuneLength = tl,
                         metric = metric, maximize = FALSE, trControl = fitControl, nestedrepeats = 10, verbose = 1)

stopCluster(cl)

# CV Results -----------------------------------------------------------------

errors <- cbind.data.frame(colMeans(base.1$results, na.rm = TRUE), colMeans(base.2$results, na.rm = TRUE), 
                           colMeans(base.3$results, na.rm = TRUE), colMeans(base.4$results, na.rm = TRUE), 
                           colMeans(base.5$results, na.rm = TRUE),colMeans(base.6$results, na.rm = TRUE), 
                           colMeans(base.7$results, na.rm = TRUE), colMeans(base.8$results, na.rm = TRUE),
                           colMeans(base.9$results, na.rm = TRUE), colMeans(base.10$results, na.rm = TRUE), 
                           colMeans(base.11$results, na.rm = TRUE), colMeans(base.12$results, na.rm = TRUE),
                           colMeans(base.13$results, na.rm = TRUE), colMeans(base.14$results, na.rm = TRUE), 
                           colMeans(base.15$results, na.rm = TRUE), colMeans(base.16$results, na.rm = TRUE),
                           colMeans(base.17$results, na.rm = TRUE), colMeans(base.18$results, na.rm = TRUE),
                           colMeans(base.19$results, na.rm = TRUE), colMeans(base.20$results, na.rm = TRUE),
                           colMeans(base.21$results, na.rm = TRUE), colMeans(base.22$results, na.rm = TRUE),
                           colMeans(base.23$results, na.rm = TRUE), colMeans(base.24$results, na.rm = TRUE),
                           colMeans(base.25$results, na.rm = TRUE), colMeans(base.26$results, na.rm = TRUE),
                           colMeans(base.27$results, na.rm = TRUE), colMeans(base.28$results, na.rm = TRUE),
                           colMeans(base.29$results, na.rm = TRUE), colMeans(base.30$results, na.rm = TRUE),
                           colMeans(base.31$results, na.rm = TRUE), colMeans(base.32$results, na.rm = TRUE),
                           colMeans(base.33$results, na.rm = TRUE), colMeans(base.34$results, na.rm = TRUE),
                           colMeans(base.35$results, na.rm = TRUE), colMeans(base.36$results, na.rm = TRUE),
                           colMeans(base.37$results, na.rm = TRUE), colMeans(base.38$results, na.rm = TRUE),
                           colMeans(base.39$results, na.rm = TRUE), colMeans(base.40$results, na.rm = TRUE),
                           colMeans(base.41$results, na.rm = TRUE), colMeans(base.42$results, na.rm = TRUE),
                           colMeans(base.43$results, na.rm = TRUE), colMeans(base.44$results, na.rm = TRUE),
                           colMeans(base.45$results, na.rm = TRUE), colMeans(base.46$results, na.rm = TRUE),
                           colMeans(base.47$results, na.rm = TRUE), colMeans(base.48$results, na.rm = TRUE),
                           colMeans(base.49$results, na.rm = TRUE), colMeans(base.50$results, na.rm = TRUE),
                           colMeans(meta.lm$results, na.rm = TRUE), colMeans(metalearner$results, na.rm = TRUE))

colnames(errors) <- c("base1", "base2", "base3", "base4", "base5",
                      "base6", "base7", "base8", "base9", "base10", 
                      "base11", "base12", "base13", "base14", "base15", 
                      "base16", "base17", "base18", "base19", "base20", 
                      "base21", "base22", "base23", "base24", "base25", 
                      "base26", "base27", "base28", "base29", "base30", 
                      "base31", "base32", "base33", "base34", "base35", 
                      "base36", "base37", "base38", "base39", "base40", 
                      "base41", "base42", "base43", "base44", "base45", 
                      "base46", "base47", "base48", "base49", "base50",
                      "LM", "ML")


# Significance Tests ------------------------------------------------------

# Significance Testing using Corrected 10x10 CV Test

#meta vs base-learners
t.tests.wmae <- cbind.data.frame(
  corrected.rep.kfold.cv.test(a = metalearner, b = base.1, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.2, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.3, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.4, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.5, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.6, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.7, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.8, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.9, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.10, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.11, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.12, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.13, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.14, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.15, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.16, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.17, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.18, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.19, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.20, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.21, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.22, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.23, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.24, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.25, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.26, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.27, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.28, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.29, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.30, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.31, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.32, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.33, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.34, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.35, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.36, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.37, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.38, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.39, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.40, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.41, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.42, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.43, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.44, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.45, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.46, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.47, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.48, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.49, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.50, metric = "WMAE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = meta.lm, metric = "WMAE"), NaN)

colnames(t.tests.wmae) <- c("base1", "base2", "base3", "base4", "base5",
                       "base6", "base7", "base8", "base9", "base10", 
                       "base11", "base12", "base13", "base14", "base15", 
                       "base16", "base17", "base18", "base19", "base20", 
                       "base21", "base22", "base23", "base24", "base25", 
                       "base26", "base27", "base28", "base29", "base30", 
                       "base31", "base32", "base33", "base34", "base35", 
                       "base36", "base37", "base38", "base39", "base40", 
                       "base41", "base42", "base43", "base44", "base45", 
                       "base46", "base47", "base48", "base49", "base50",
                       "LM", "ML")

#meta vs base-learners
t.tests.rmse <- cbind.data.frame(
  corrected.rep.kfold.cv.test(a = metalearner, b = base.1, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.2, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.3, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.4, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.5, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.6, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.7, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.8, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.9, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.10, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.11, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.12, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.13, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.14, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.15, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.16, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.17, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.18, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.19, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.20, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.21, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.22, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.23, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.24, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.25, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.26, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.27, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.28, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.29, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.30, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.31, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.32, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.33, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.34, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.35, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.36, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.37, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.38, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.39, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.40, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.41, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.42, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.43, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.44, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.45, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.46, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.47, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.48, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.49, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.50, metric = "RMSE"),
  corrected.rep.kfold.cv.test(a = metalearner, b = meta.lm, metric = "RMSE"), NaN)

colnames(t.tests.rmse) <- c("base1", "base2", "base3", "base4", "base5",
                       "base6", "base7", "base8", "base9", "base10", 
                       "base11", "base12", "base13", "base14", "base15", 
                       "base16", "base17", "base18", "base19", "base20", 
                       "base21", "base22", "base23", "base24", "base25", 
                       "base26", "base27", "base28", "base29", "base30", 
                       "base31", "base32", "base33", "base34", "base35", 
                       "base36", "base37", "base38", "base39", "base40", 
                       "base41", "base42", "base43", "base44", "base45", 
                       "base46", "base47", "base48", "base49", "base50",
                       "LM", "ML")

t.tests.100 <- cbind.data.frame(
  corrected.rep.kfold.cv.test(a = metalearner, b = base.1, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.2, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.3, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.4, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.5, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.6, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.7, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.8, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.9, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.10, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.11, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.12, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.13, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.14, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.15, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.16, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.17, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.18, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.19, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.20, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.21, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.22, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.23, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.24, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.25, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.26, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.27, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.28, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.29, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.30, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.31, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.32, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.33, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.34, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.35, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.36, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.37, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.38, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.39, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.40, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.41, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.42, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.43, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.44, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.45, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.46, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.47, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.48, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.49, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = base.50, metric = "RMSE.100"),
  corrected.rep.kfold.cv.test(a = metalearner, b = meta.lm, metric = "RMSE.100"), NaN)

colnames(t.tests.100) <- c("base1", "base2", "base3", "base4", "base5",
                           "base6", "base7", "base8", "base9", "base10", 
                           "base11", "base12", "base13", "base14", "base15", 
                           "base16", "base17", "base18", "base19", "base20", 
                           "base21", "base22", "base23", "base24", "base25", 
                           "base26", "base27", "base28", "base29", "base30", 
                           "base31", "base32", "base33", "base34", "base35", 
                           "base36", "base37", "base38", "base39", "base40", 
                           "base41", "base42", "base43", "base44", "base45", 
                           "base46", "base47", "base48", "base49", "base50",
                           "LM", "ML")


# Aggregating and Saving Results -----------------------------------------------------

#Performance on big Halloween Storms
sub <- paper2_data[paper2_data$ID %in% c(8039,
                                         8048,
                                         8049,
                                         8050,
                                         8078,
                                         8109,
                                         8111,
                                         8125,
                                         8126,
                                         8127,
                                         8140,
                                         8149,
                                         8153,
                                         8154,
                                         8155,
                                         8156,
                                         8161,
                                         8163,
                                         8164,
                                         8166,
                                         8171,
                                         8172,
                                         8173,
                                         8180,
                                         8181,
                                         8182,
                                         8188
),]
sub.bl <- base.11$fit$pred[base.11$fit$pred$rowIndex %in% as.numeric(rownames(sub)),]
sub.bl <- arrange(sub.bl, rowIndex)
sub.bl <- sub.bl[,c("pred", "rowIndex")]
sub.ml <- metalearner$fit$pred[metalearner$fit$pred$rowIndex %in% as.numeric(rownames(sub)),]
sub.ml <- arrange(sub.ml, rowIndex)
sub.ml <- sub.ml[,c("pred", "rowIndex")]
sub.ml.lm <- meta.lm$fit$pred[meta.lm$fit$pred$rowIndex %in% as.numeric(rownames(sub)),]
sub.ml.lm <- arrange(sub.ml.lm, rowIndex)
sub.ml.lm <- sub.ml.lm[,c("pred", "rowIndex")]

sub.all <- cbind.data.frame("DT" = sub$LASCO_Date_Time, "DST" = sub$Min_DST,
                        "Best Base Prediction" = 
                          sub.bl[,"pred"],
                        "Meta-learner Prediction" = 
                          sub.ml[,"pred"],
                        "Meta-lm Prediction" = 
                          sub.ml.lm[,"pred"])

write.csv(sub.all, file = "H:/tklarkin_backup/Documents/My Research/Journal Articles/Venice Conference Extended/Code/stage1_halloween.csv")

all.results <- cbind.data.frame(t(errors), "WMAE" = t(t.tests.wmae), "RMSE" = t(t.tests.rmse), "RMSE.100" = t(t.tests.100))

write.csv(cbind.data.frame(round(all.results, digits = 4)), 
          file = "H:/tklarkin_backup/Documents/My Research/Journal Articles/Venice Conference Extended/Code/stage1_meta_learning_results.csv")

# Normalized Variable Importance Scheme for SG -------------------------------------

#Calculating variable importances from all base-learners

varimp.base.1 <- as.data.frame(varImp(base.1$fit)$importance)
varimp.base.1$RN <- rownames(varimp.base.1)
varimp.base.1 <- arrange(varimp.base.1, RN)

varimp.base.2 <- as.data.frame(varImp(base.2$fit)$importance)
varimp.base.2$RN <- rownames(varimp.base.2)
varimp.base.2 <- arrange(varimp.base.2, RN)

varimp.base.3 <- as.data.frame(varImp(base.3$fit)$importance)
varimp.base.3$RN <- rownames(varimp.base.3)
varimp.base.3 <- arrange(varimp.base.3, RN)

varimp.base.4 <- as.data.frame(varImp(base.4$fit)$importance)
varimp.base.4$RN <- rownames(varimp.base.4)
varimp.base.4 <- arrange(varimp.base.4, RN)

varimp.base.5 <- as.data.frame(varImp(base.5$fit)$importance)
varimp.base.5$RN <- rownames(varimp.base.5)
varimp.base.5 <- arrange(varimp.base.5, RN)

varimp.base.6 <- as.data.frame(varImp(base.6$fit)$importance)
varimp.base.6$RN <- rownames(varimp.base.6)
varimp.base.6 <- arrange(varimp.base.6, RN)

varimp.base.7 <- as.data.frame(varImp(base.7$fit)$importance)
varimp.base.7$RN <- rownames(varimp.base.7)
varimp.base.7 <- arrange(varimp.base.7, RN)

varimp.base.8 <- as.data.frame(varImp(base.8$fit)$importance)
varimp.base.8$RN <- rownames(varimp.base.8)
varimp.base.8 <- arrange(varimp.base.8, RN)

#Needs seed if ran by itself for reproducible results
varimp.base.9 <- as.data.frame(varImp(base.9$fit)$importance)
varimp.base.9$RN <- rownames(varimp.base.9)
varimp.base.9 <- arrange(varimp.base.9, RN)

varimp.base.10 <- as.data.frame(varImp(base.10$fit)$importance)
varimp.base.10$RN <- rownames(varimp.base.10)
varimp.base.10 <- arrange(varimp.base.10, RN)

varimp.base.11 <- as.data.frame(varImp(base.11$fit)$importance)
varimp.base.11$RN <- rownames(varimp.base.11)
varimp.base.11 <- arrange(varimp.base.11, RN)

varimp.base.12 <- as.data.frame(varImp(base.12$fit)$importance)
varimp.base.12$RN <- rownames(varimp.base.12)
varimp.base.12 <- arrange(varimp.base.12, RN)
varimp.base.12 <- merge(varimp.base.12, varimp.base.1, by = "RN", all = T)[,-3]
varimp.base.12[is.na(varimp.base.12)] <- 0
varimp.base.12 <- varimp.base.12[c(2,1)]
colnames(varimp.base.12) <- c("Overall", "RN")
varimp.base.12 <- arrange(varimp.base.12, RN)

varimp.base.13 <- as.data.frame(varImp(base.13$fit)$importance)
varimp.base.13$RN <- rownames(varimp.base.13)
varimp.base.13 <- arrange(varimp.base.13, RN)
varimp.base.13 <- merge(varimp.base.13, varimp.base.1, by = "RN", all = T)[,-3]
varimp.base.13[is.na(varimp.base.13)] <- 0
varimp.base.13 <- varimp.base.13[c(2,1)]
colnames(varimp.base.13) <- c("Overall", "RN")
varimp.base.13 <- arrange(varimp.base.13, RN)

varimp.base.14 <- as.data.frame(varImp(base.14$fit)$importance)
varimp.base.14$RN <- rownames(varimp.base.14)
varimp.base.14 <- arrange(varimp.base.14, RN)

varimp.base.15 <- as.data.frame(varImp(base.15$fit)$importance)
varimp.base.15$RN <- rownames(varimp.base.15)
varimp.base.15 <- arrange(varimp.base.15, RN)
varimp.base.15 <- merge(varimp.base.15, varimp.base.1, by = "RN", all = T)[,-3]
varimp.base.15[is.na(varimp.base.15)] <- 0
varimp.base.15 <- varimp.base.15[c(2,1)]
colnames(varimp.base.15) <- c("Overall", "RN")
varimp.base.15 <- arrange(varimp.base.15, RN)

varimp.base.16 <- as.data.frame(varImp(base.16$fit)$importance)
varimp.base.16$RN <- rownames(varimp.base.16)
varimp.base.16 <- arrange(varimp.base.16, RN)

varimp.base.17 <- as.data.frame(varImp(base.17$fit)$importance)
varimp.base.17$RN <- rownames(varimp.base.17)
varimp.base.17 <- arrange(varimp.base.17, RN)

varimp.base.18 <- as.data.frame(varImp(base.18$fit)$importance)
varimp.base.18$RN <- rownames(varimp.base.18)
varimp.base.18 <- arrange(varimp.base.18, RN)

varimp.base.19 <- as.data.frame(varImp(base.19$fit)$importance)
varimp.base.19$RN <- rownames(varimp.base.19)
varimp.base.19 <- arrange(varimp.base.19, RN)

varimp.base.20 <- as.data.frame(varImp(base.20$fit)$importance)
varimp.base.20$RN <- rownames(varimp.base.20)
varimp.base.20 <- arrange(varimp.base.20, RN)

varimp.base.21 <- as.data.frame(varImp(base.21$fit)$importance)
varimp.base.21$RN <- rownames(varimp.base.21)
varimp.base.21 <- arrange(varimp.base.21, RN)

varimp.base.22 <- as.data.frame(varImp(base.22$fit)$importance)
varimp.base.22$RN <- rownames(varimp.base.22)
varimp.base.22 <- arrange(varimp.base.22, RN)
varimp.base.22 <- merge(varimp.base.22, varimp.base.1, by = "RN", all = T)[,-3]
varimp.base.22[is.na(varimp.base.22)] <- 0
varimp.base.22 <- varimp.base.22[c(2,1)]
colnames(varimp.base.22) <- c("Overall", "RN")
varimp.base.22 <- arrange(varimp.base.22, RN)

varimp.base.23 <- as.data.frame(varImp(base.23$fit)$importance)
varimp.base.23$RN <- rownames(varimp.base.23)
varimp.base.23 <- arrange(varimp.base.23, RN)

varimp.base.24 <- as.data.frame(varImp(base.24$fit)$importance)
varimp.base.24$RN <- rownames(varimp.base.24)
varimp.base.24 <- arrange(varimp.base.24, RN)

varimp.base.25 <- as.data.frame(varImp(base.25$fit)$importance)
varimp.base.25$RN <- rownames(varimp.base.25)
varimp.base.25 <- arrange(varimp.base.25, RN)

varimp.base.26 <- as.data.frame(varImp(base.26$fit)$importance)
varimp.base.26$RN <- rownames(varimp.base.26)
varimp.base.26 <- arrange(varimp.base.26, RN)

varimp.base.27 <- as.data.frame(varImp(base.27$fit)$importance)
varimp.base.27$RN <- rownames(varimp.base.27)
varimp.base.27 <- arrange(varimp.base.27, RN)

varimp.base.28 <- as.data.frame(varImp(base.28$fit)$importance)
varimp.base.28$RN <- rownames(varimp.base.28)
varimp.base.28 <- arrange(varimp.base.28, RN)

varimp.base.29 <- as.data.frame(varImp(base.29$fit)$importance)
varimp.base.29$RN <- rownames(varimp.base.29)
varimp.base.29 <- arrange(varimp.base.29, RN)

varimp.base.30 <- as.data.frame(varImp(base.30$fit)$importance)
varimp.base.30$RN <- rownames(varimp.base.30)
varimp.base.30 <- arrange(varimp.base.30, RN)

varimp.base.31 <- as.data.frame(varImp(base.31$fit)$importance)
varimp.base.31$RN <- rownames(varimp.base.31)
varimp.base.31 <- arrange(varimp.base.31, RN)

varimp.base.32 <- as.data.frame(varImp(base.32$fit)$importance)
varimp.base.32$RN <- rownames(varimp.base.32)
varimp.base.32 <- arrange(varimp.base.32, RN)

varimp.base.33 <- as.data.frame(varImp(base.33$fit)$importance)
varimp.base.33$RN <- rownames(varimp.base.33)
varimp.base.33 <- arrange(varimp.base.33, RN)

varimp.base.34 <- as.data.frame(varImp(base.34$fit)$importance)
varimp.base.34$RN <- rownames(varimp.base.34)
varimp.base.34 <- arrange(varimp.base.34, RN)

varimp.base.35 <- as.data.frame(varImp(base.35$fit)$importance)
varimp.base.35$RN <- rownames(varimp.base.35)
varimp.base.35 <- arrange(varimp.base.35, RN)

varimp.base.36 <- as.data.frame(varImp(base.36$fit)$importance)
varimp.base.36$RN <- rownames(varimp.base.36)
varimp.base.36 <- arrange(varimp.base.36, RN)

varimp.base.37 <- as.data.frame(varImp(base.37$fit)$importance)
varimp.base.37$RN <- rownames(varimp.base.37)
varimp.base.37 <- arrange(varimp.base.37, RN)

varimp.base.38 <- as.data.frame(varImp(base.38$fit)$importance)
varimp.base.38$RN <- rownames(varimp.base.38)
varimp.base.38 <- arrange(varimp.base.38, RN)

varimp.base.39 <- as.data.frame(varImp(base.39$fit)$importance)
varimp.base.39$RN <- rownames(varimp.base.39)
varimp.base.39 <- arrange(varimp.base.39, RN)

varimp.base.40 <- as.data.frame(varImp(base.40$fit)$importance)
varimp.base.40$RN <- rownames(varimp.base.40)
varimp.base.40 <- arrange(varimp.base.40, RN)

varimp.base.41 <- as.data.frame(varImp(base.41$fit)$importance)
varimp.base.41$RN <- rownames(varimp.base.41)
varimp.base.41 <- arrange(varimp.base.41, RN)

varimp.base.42 <- as.data.frame(varImp(base.42$fit)$importance)
varimp.base.42$RN <- rownames(varimp.base.42)
varimp.base.42 <- arrange(varimp.base.42, RN)

varimp.base.43 <- as.data.frame(varImp(base.43$fit)$importance)
varimp.base.43$RN <- rownames(varimp.base.43)
varimp.base.43 <- arrange(varimp.base.43, RN)

varimp.base.44 <- as.data.frame(varImp(base.44$fit)$importance)
varimp.base.44$RN <- rownames(varimp.base.44)
varimp.base.44 <- arrange(varimp.base.44, RN)

varimp.base.45 <- as.data.frame(varImp(base.45$fit)$importance)
varimp.base.45$RN <- rownames(varimp.base.45)
varimp.base.45 <- arrange(varimp.base.45, RN)

varimp.base.46 <- as.data.frame(varImp(base.46$fit)$importance)
varimp.base.46$RN <- rownames(varimp.base.46)
varimp.base.46 <- arrange(varimp.base.46, RN)

varimp.base.47 <- as.data.frame(varImp(base.47$fit)$importance)
varimp.base.47$RN <- rownames(varimp.base.47)
varimp.base.47 <- arrange(varimp.base.47, RN)

varimp.base.48 <- as.data.frame(varImp(base.48$fit)$importance)
varimp.base.48$RN <- rownames(varimp.base.48)
varimp.base.48 <- arrange(varimp.base.48, RN)

varimp.base.49 <- as.data.frame(varImp(base.49$fit)$importance)
varimp.base.49$RN <- rownames(varimp.base.49)
varimp.base.49 <- arrange(varimp.base.49, RN)

varimp.base.50 <- as.data.frame(varImp(base.50$fit)$importance)
varimp.base.50$RN <- rownames(varimp.base.50)
varimp.base.50 <- arrange(varimp.base.50, RN)

#Combining the variable importance scores
varimp.all <- cbind.data.frame(varimp.base.1, varimp.base.2, varimp.base.3,
                               varimp.base.4, varimp.base.5, varimp.base.6,
                               varimp.base.7, varimp.base.8, varimp.base.9,
                               varimp.base.10, varimp.base.11, varimp.base.12,
                               varimp.base.13, varimp.base.14, varimp.base.15,
                               varimp.base.16, varimp.base.17, varimp.base.18,
                               varimp.base.19, varimp.base.20, varimp.base.21,
                               varimp.base.22, varimp.base.23, varimp.base.24,
                               varimp.base.25, varimp.base.26, varimp.base.27, 
                               varimp.base.28, varimp.base.29, varimp.base.30,
                               varimp.base.31, varimp.base.32, varimp.base.33,
                               varimp.base.34, varimp.base.35, varimp.base.36,
                               varimp.base.37, varimp.base.38, varimp.base.39,
                               varimp.base.40, varimp.base.41, varimp.base.42,
                               varimp.base.43, varimp.base.44, varimp.base.45,
                               varimp.base.46, varimp.base.47, varimp.base.48,
                               varimp.base.49, varimp.base.50)
dim(varimp.all)

vardata <- subset(varimp.all, select=-c(seq(4,ncol(varimp.all),2)))
vardata <- vardata[c(2,1,seq(3,length(vardata),1))]
vardata[is.na(vardata)] <- 0 #nnls produced all zero coefficients
colnames(vardata) <- c("Variable.Name", "base01", "base02", "base03", "base04", 
                       "base05", "base06", "base07", "base08", "base09", "base10", 
                       "base11", "base12", "base13", "base14", "base15", 
                       "base16", "base17", "base18", "base19", "base20", 
                       "base21", "base22", "base23", "base24", "base25", 
                       "base26", "base27", "base28", "base29", "base30", 
                       "base31", "base32", "base33", "base34", "base35", 
                       "base36", "base37", "base38", "base39", "base40", 
                       "base41", "base42", "base43", "base44", "base45", 
                       "base46", "base47", "base48", "base49", "base50")

#Obtaining Weights from meta
meta.imp <- as.data.frame(varImp(metalearner$fit, scale = FALSE)$importance)

#10 most important base-learners
meta.imp[order(meta.imp$Overall, decreasing = TRUE)[1:10], , drop = FALSE]

#No need to transpose because the weights are listed as a vector here
importance <- data.frame("Final Variable Importance" = 
                           (as.matrix(vardata[,-1]) %*% as.matrix(meta.imp$Overall)))
rownames(importance) <- vardata$Variable.Name
final.imp <- (importance[order(-importance), , drop = FALSE]-min(importance[order(-importance), , drop = FALSE]))/
  (max(importance[order(-importance), , drop = FALSE])-min(importance[order(-importance), , drop = FALSE])) * 100

final.imp

#Graphing Proposed Variable Importance Ranking for SG
graph.final.imp <- final.imp
graph.final.imp$Variable.Name <- rownames(final.imp)

pdf("H:/tklarkin_backup/Documents/My Research/Journal Articles/Venice Conference Extended/Latex File/stage1.pdf", 
    width = 11, height = 6)
ggplot(graph.final.imp, aes(y=Final.Variable.Importance, x = Variable.Name, fill=Final.Variable.Importance)) +
  geom_bar(colour="black", width=0.8, stat="identity") + coord_flip() + 
  scale_colour_gradient(limits=c(0, 100)) +
  guides(fill=FALSE) +
  xlab("Predictor Variables") + ylab("Importance Score") +
  theme(text = element_text(size=12),
        axis.title.y = element_text(angle=90, vjust = 1.5))
dev.off()

save.image("H:/tklarkin_backup/Documents/My Research/Journal Articles/Venice Conference Extended/Code/Stage 1 Evaluating Error for SG.RData")
