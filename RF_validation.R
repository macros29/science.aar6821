setwd("~/Data/Da_Mi/")

library(parallel)

k <- 10 # 10-fold cross validation

no_cores <- detectCores() - 1
feature <- read.csv("2018-01-11_RF_features_assigned_neuron.csv", header = T,row.names = 1, 
                    stringsAsFactors = F)
cate <- read.csv("2018-01-11_RF_cate_assigned_neuron.csv", header = T,row.names = 1, 
                 stringsAsFactors = F)
cate <- cate$x

feature <- feature[order(cate), ]
cate <- cate[order(cate)]

cl <- makeCluster(no_cores, outfile='LOG.TXT')
clusterExport(cl, c("k","feature", "cate"))

start.time <- Sys.time()
err_fs <- parLapply(cl, 1:1000, function(xx){
  library(randomForest)
  # subset <- unlist(sapply(unique(cate), function(xx){
  #   n <- floor(length(which(cate == xx)) / k)
  #   sample(c(1:nrow(feature))[which(cate == xx)], n, replace = F)
  # }))
  n <- floor(nrow(feature) / k)
  subset <- sample(1: nrow(feature), n, replace = F)
  train_feature <- feature[-subset, ]
  train_cate <- cate[-subset]
  test_feature <- feature[subset, ]
  test_cate <- cate[subset]
  retry <- -1
  # # Ensure each training set has > 5 samples
  if(min(table(train_cate)) < 5) {
    retry <- 100
    while(retry > 0){
      subset <- sample(1: nrow(feature), n, replace = F)
      train_feature <- feature[-subset, ]
      train_cate <- cate[-subset]
      test_feature <- feature[subset, ]
      test_cate<- cate[subset]
      if (min(table(train_cate)) < 5) {
        retry <- retry - 1
      }  else {
        retry <- -1
      }
    }
  }
  if (retry == 0) {
    stop("The test set is too small!")
  }
  if(retry == -1){
    rf <- randomForest(train_feature, as.factor(train_cate), importance = TRUE, 
                       proximity = TRUE)
    pred2_fs <- predict(rf, newdata = test_feature)
    mis <- test_cate == pred2_fs
    mis <- data.frame(sample = rownames(test_feature), value = mis)
    # 
    # cate <- as.factor(c(as.character(train_cate), as.character(pred2_fs)))
    # feature <- as.matrix(rbind(train_feature, test_feature))
    # set <- sample(1: nrow(feature), nrow(feature), replace=F)
    # cate <- cate[set]
    # feature <- feature[set, ]
  }
  return(mis)
})
stopCluster(cl)
runtime <- Sys.time() - start.time

temp <- merge(err_fs[[1]], err_fs[[2]], by = "sample", all = T)
for(i in 3:1000){
  temp <- merge(temp, err_fs[[i]], by = "sample", all = T)
}

res <- apply(temp[,-1], 1, function(xx){table(xx, useNA = "ifany")})
res <- sapply(res, function(xx){
  if(is.na(xx["FALSE"])){
    "Pass"
  } else if(is.na(xx["TRUE"])){
    "Fail"
  } else if(xx["FALSE"] / (xx["FALSE"] + xx["TRUE"]) < 0.05) {
    "Pass"
  } else { "Fail" }
})

res <- data.frame(sample = temp$sample, res = res)
res <- res[match(rownames(feature), as.character(res$sample)),]

table(cate[res$res == "Pass"])

write.csv(res, file = "2018-01-11_RF_validation.csv")
