###################################################################################################
###
### Clustering and Classification
###
### 
### Based on Rizi Ai's cluster-and-classification script
### 
### zhen.li.zl242@yale.edu  
### Yale University
###################################################################################################
### Required libraries
###################################################################################################
library(randomForest)
library(gplots)
library(matrixStats)
library(Seurat)
library(pcaReduce)
source("~/Scripts/R/heatmap.3.R")

## For the first run I upload de whole matrix
args <- commandArgs(trailingOnly =T)
setwd(paste(args[1], args[2], sep = "_"))
print(paste("Start analyze", args[1], args[2], sep = " "))

load("../raw_data/expr.Rdata")
load("../raw_data/GB.Rdata")

num=as.numeric(args[2]) + 1
name=paste('Round',num,'.',args[3],'.',sep='')

keep <- GB$region %in% args & GB$time_point %in% args & !(GB$IFC %in% c("189","180","181","182","186","190","191","192","193","194","195", "196"))
#GB.keep <- GB[keep, ]

clust_method = "ward.D2"

###
expression_type = 1 # "1" for log2 expression (eg. log2TPM, log2FPKM, log2Counts); 
                    # "2" for expression without taking log2 transformation

###################################################################################################
set.seed(1)
write(format(Sys.time(), "%Y-%m-%d %I-%p"), "01.Report.txt")
write(format(Sys.time(), "%Y-%m-%d %I-%p"), "log")
write(paste(rep("#", 100), collapse = ""), "01.Report.txt", append = T)
write(paste(rep("#", 100), collapse = ""), "log", append = T)
write("### Clustering and Classification\n###\n### Version 1.1", "01.Report.txt", append = T)
write("### Based on  Rizi Ai's cluster-and-classification.R script\n### zhen.li.zl242@yale.edu\n### Yale University", "01.Report.txt", append = T)
write("### Clustering and Classification\n###\n### Version 1.1", "log", append = T)
write("### Based on  Rizi Ai's cluster-and-classification.R script\n### zhen.li.zl242@yale.edu\n### Yale University", "log", append = T)
write(paste(rep("#", 100), collapse = ""), "01.Report.txt", append = T)
write(paste(rep("#", 100), collapse = ""), "log", append = T)

###################################################################################################
### 0. Checking input gene expression matrix
###################################################################################################
print("... Running 0.Checking input gene expression matrix")
write("... Running 0.Checking input gene expression matrix", "log", append = T)
# gene_mat <- as.matrix(Y[ ,keep ]) 
# if (expression_type == 2) {
#   gene_mat[gene_mat < 1] <- 1 #genes with expression level < 1 are considered not expressed
#   gene_mat <- log2(gene_mat)  
# }
# write(paste("Input file: ", ncol(gene_mat)," samples, ", nrow(gene_mat), " genes"), "01.Report.txt", append = T)

#Check empty entries
clean_mat <- function(x) {
  row_rm <- vector()
  col_rm <- vector()
  for (i in 1:nrow(x)) {
    if((sum(x[i, ] == 0) == ncol(x)) || (sum(is.na(x[i, ])) == ncol(x))){
      row_rm <- c(row_rm, i)
    }
  }
  for (i in 1:ncol(x)) {
    if((sum(x[ ,i] == 0) == nrow(x)) || (sum(is.na(x[ ,i])) == nrow(x))){
      col_rm <- c(col_rm, i)
    }
  }
  return (list(row_rm, col_rm))
}

# input_check <- clean_mat(gene_mat)
# if (length(input_check[[1]]) > 0) {
#  gene_mat <- gene_mat[-input_check[[1]], ]
# }
# if (length(input_check[[2]]) > 0) {
#  gene_mat <-gene_mat[ , -input_check[[2]]]
# }
 
print("... Done")
write("... Done", "log", append = T)

###################################################################################################
###
######COMPUTING OVERDISPERSION######
###
###################################################################################################
### 1.Calculating gene expression variation CV^2
###################################################################################################
print(paste(rep("#", 80), collapse = ""))
write(paste(rep("#", 80), collapse = ""), "log", append = T)
print("COMPUTING OVERDISPERSION")
write("COMPUTING OVERDISPERSION", "log", append = T)
print(paste(rep("#", 80), collapse = ""))
write(paste(rep("#", 80), collapse = ""), "log", append = T)
print("... Running 1.Calculating gene expression variation")
write("... Running 1.Calculating gene expression variation", "log", append = T)

# write.table(t(c("Mean", "Variance", "CV2")), "02.Expression-variation_stat.txt", quote = F, row.names = F, col.names = F, sep = "\t")
# no_cores <- detectCores() - 1
# cl <- makeCluster(no_cores)
# geneNames<- rownames(gene_mat)
# clusterExport(cl, "geneNames")
# parApply(cl, gene_mat, 1, function(x){
#   cv_avg <- mean(x, na.rm = T)
#   var <- var(x, na.rm=T)
#   cv2 <- var / (cv_avg * cv_avg)
#   write.table(t(c(round(cv_avg, digits = 3), 
#                   round(var, digits = 3), 
#                   round(cv2, digits = 3))), 
#               "02.Expression-variation_stat.txt", quote = F, row.names = F, 
#               col.names = F, sep = "\t", append = T)
# })
# stopCluster(cl)

print("Skip\n")
write("Skip\n", append = T)

# for(i in 1: nrow(gene_mat)){
#   cv_avg <- mean(gene_mat[i,], na.rm = T)
#   var <- var(gene_mat[i,], na.rm=T)
#   cv2 <- var / (cv_avg * cv_avg)
#   write.table(t(c(rownames(gene_mat)[i], 
#                   round(cv_avg, digits = 3), 
#                   round(var, digits = 3), 
#                   round(cv2, digits = 3))), 
#               "02.Expression-variation_stat.txt", quote = F, row.names = F, 
#               col.names = F, sep = "\t", append = T)
# }
print("... Writing expression variation stat to 02.Expression-variation_stat.txt")
write("... Writing expression variation stat to 02.Expression-variation_stat.txt", "log", append = T)
print("... Done")
write("... Done", "log", append = T)

###################################################################################################
### 2. Computing over-dispersed genes
###################################################################################################
print("... Running 2. Computing over-dispersed genes")
write("... Running 2. Computing over-dispersed genes", "log", append = T)
# cv_dat1 <- read.table("02.Expression-variation_stat.txt", header=T)
# index <- vector()
# for (i in 1: nrow(cv_dat1)){
#   scale <- 2 #Genes with mean expression (log2TPM) > 2 were used for fitting.
#   #This can be scaled to specific gene expression quantification. eg. quantile, mean, median, etc.
#   if ((!is.na(cv_dat1$Mean[i])) && (cv_dat1$Mean[i] > scale)){
#     index <- c(index,i)
#   }
# }
# cv_dat2 <- cv_dat1[index, ]
# inv_avg = 1 / cv_dat2$Mean
# fit <- lm(cv_dat2$CV2 ~ inv_avg) #fitting CV^2 and mean expression to an inverse distribution
# l <- fit$coefficients[2] / cv_dat1$Mean + fit$coefficients[1]
# se <- summary(fit)$sigma
# se_1 <- l + se
# pdf("03.Fitting_stat.pdf")
# plot(cv_dat1$Mean, cv_dat1$CV2, xlab = "Mean", ylab = "CV2", main = "Fitting gene expression variation", col = "blue", lwd = 2)
# points(cv_dat1$Mean, l, lty = 2, lwd = 1, col = "red")
# points(cv_dat1$Mean, se_1, col="red")
# dev.off()
# out_index <- vector()
# for(i in 1: nrow(gene_mat)){
#   a <- fit$coefficients[2] / cv_dat1$Mean[i] + fit$coefficients[1]
#   if((!is.na(cv_dat1$CV2[i])) && (cv_dat1$CV2[i] > (a+se))){
#     out_index <- c(out_index, i)
#   }
# }
# m <- gene_mat[out_index, ]

gene_mat <- expr[,keep]
df <- CreateSeuratObject(raw.data = as.vector(expr[,keep]), min.cells = 5, min.genes = 1000, project = "Da_Mi")
df <- NormalizeData(object = df, normalization.method = "LogNormalize", scale.factor = 10000)
df <- ScaleData(df)
df <- FindVariableGenes(object = df, mean.function = ExpMean, dispersion.function = LogVMR, 
                        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
m <- as.matrix(df@data[df@var.genes,])

# logMeans <- rowMeans(gene_mat)
# logVars <- rowVars(as.matrix(gene_mat))
# logCV2 <- logVars / logMeans^2
# 
# useForFitL <- logMeans > 1
# LogNcountsList <- data.frame(mean = logMeans, cv2 = logCV2)
# 
# fit_loglin = lm(log10(cv2)~mean, data = LogNcountsList[useForFitL,])
# ElogCV2 <- predict(fit_loglin, 
#                    data.frame(mean = LogNcountsList$mean)) - summary(fit_loglin)$sigma
# 
# out_index <- (ElogCV2 < log10(logCV2)) & useForFitL
# m <- gene_mat[out_index, ]

write.table(rownames(m), "04.Over-dispersed_gene-list.txt", quote = F, sep = "\t", row.names = F, col.names = F)
print("... Writing over-dispersed genes to 04.Over-dispersed_gene-list.txt")
write("... Writing over-dispersed genes to 04.Over-dispersed_gene-list.txt", "log", append = T)
print("... Done")
write("... Done", "log", append = T)

###################################################################################################
###
######CLUSTERING STEP######
###
###################################################################################################
### 3. Initial hierarchical clustering by over-dispersed genes 
###################################################################################################
print(paste(rep("#", 80), collapse = ""))
write(paste(rep("#", 80), collapse = ""), "log", append = T)
print("CLUSTERING STEP")
write("CLUSTERING STEP", "log", append = T)
print(paste(rep("#", 80), collapse = ""))
write(paste(rep("#", 80), collapse = ""), "log", append = T)
print("... Running 3. Initial hierarchical clustering by over-dispersed genes using pcaReduce")
write("... Running 3. Initial hierarchical clustering by over-dispersed genes using pcaReduce", "log", append = T)

myCol <- c("darkblue", "white", "red")
hc <- hclust(as.dist(1-cor(m, method = "pearson", use = "pairwise.complete.obs")), method = clust_method)
cl <- cutree(hc, k=2)

df@ident <- as.factor(cl)
 
# nbt = 1
# q = 2
# seed <- 1024
# set.seed(seed)
# Output_S <- PCAreduce(t(m), nbt = nbt, q = q, method = 'S')
# cl <- Output_S[[1]][,2]

write(paste("Initial cluster sizes | ", table(cl)[1], ", ", table(cl)[2]), "01.Report.txt", append = T)

if ((table(cl)[1] < 10) || (table(cl)[2] < 10 )) {
  stop("The sizes of one or both clusters are too small for further analysis!")
}

cl_a <- m[ , cl == 1]
cl_b <- m[ , cl == 2]

print("... Done")
write("... Done", "log", append = T)

###################################################################################################
### 4. Adjusting clusters using 2-fold over-dispersed DEGs (Optional)
###################################################################################################
# print("... Running 4. Adjusting clusters using 2-fold over-dispersed DEGs (Optional)")
# write("... Running 4. Adjusting clusters using 2-fold over-dispersed DEGs (Optional)", "log", append = T)
# get_deg <- function(m1, m2, fold){
#   p_all <- vector()
#   delta_all <- vector()
#   mean_all_1 <- vector()
#   mean_all_2 <- vector()
#   index_fc_1 <- vector()
#   index_fc_2 <- vector()
#   p_adj <- vector()
#   for (i in 1:nrow(m1)) {
#     ttest <- t.test(m1[i, ], m2[i, ])
#     p <- ttest$p.value
#     p_all <- c(p_all, p)
#     avg_1 <- mean(as.numeric(m1[i, ]), na.rm = T)
#     avg_2 <- mean(as.numeric(m2[i, ]), na.rm = T)
#     mean_all_1 <- c(mean_all_1, avg_1)
#     mean_all_2 <- c(mean_all_2, avg_2)
#     delta <- avg_1 - avg_2
#     delta_all <- c(delta_all, delta)
#     if ((!is.na(delta)) && (!is.na(p)) && (p < 0.05)) {
#       if (delta > log2(fold)) {
#         index_fc_1 <- c(index_fc_1, i)
#       }  else if(delta < -log2(fold)) {
#         index_fc_2 <- c(index_fc_2, i)
#       }
#     }
#   }
#   if (length(p_all) > 200) {
#     p_adj <- p.adjust(p_all, method = "BH")
#   }
#   return (list(index_fc_1, index_fc_2, mean_all_1, mean_all_2, delta_all, p_all, p_adj))
# }

# get_2f_deg <- get_deg(cl_a, cl_b, 2)
# m_hc <- m[c(get_2f_deg[[1]], get_2f_deg[[2]]), ]

deg <- FindAllMarkers(df, genes.use =  df@var.genes, only.pos = T)
m_hc <- m[deg$gene, ]

##### Check empty entries
clean_m_hc <- clean_mat(m_hc)
if (length(clean_m_hc[[2]]) > 0){
  m_hc <- m_hc[ , -clean_m_hc[[2]]]
}

# nbt = 1
# q = 2
# seed <- 1024
# set.seed(seed)
# Output_S <- PCAreduce(t(m_hc), nbt = nbt, q = q, method = 'S')
# cl_hc <- Output_S[[1]][,2]

hc_deg <- hclust(as.dist(1 - cor(m_hc, method = "pearson", use = "pairwise.complete.obs")), method = clust_method)
cl_hc <- cutree(hc_deg, k = 2)
write(paste("Cluster sizes at CLUSTERING step | ", table(cl_hc)[1], ", ", table(cl_hc)[2]), "01.Report.txt", append = T)
print("... Done")
write("... Done", "log", append = T)

# cl_hc <- cl
# m_hc <- m

###################################################################################################
###
######CLASSIFICATION STEP######
###
###################################################################################################
### 5. Random Forest on clusters from clustering step   
###################################################################################################
print(paste(rep("#", 80), collapse = ""))
write(paste(rep("#", 80), collapse = ""), "log", append = T)
print("CLASSIFICATION STEP")
write("CLASSIFICATION STEP", "log", append = T)
print(paste(rep("#", 80), collapse = ""))
write(paste(rep("#", 80), collapse = ""), "log", append = T)
print("... Running 5. Random Forest on clusters from clustering step")
write("... Running 5. Random Forest on clusters from clustering step", "log", append = T)
cl_deg_a <- t(m_hc[ , cl_hc == 1])
cl_deg_b <- t(m_hc[ , cl_hc == 2])
cate_bf_fs <- as.factor(c(rep("a", nrow(cl_deg_a)), rep("b", nrow(cl_deg_b))))
feature_bf_fs <- as.matrix(rbind(cl_deg_a, cl_deg_b))
set <- sample(1: nrow(feature_bf_fs), nrow(feature_bf_fs), replace=F)  ## Random reorder
cate_bf_fs <- cate_bf_fs[set]
feature_bf_fs <- feature_bf_fs[set, ]
rf_bf_fs <- randomForest(feature_bf_fs, cate_bf_fs, importance = TRUE, proximity = TRUE)
imp_bf_fs <- importance(rf_bf_fs, type = 1)
print("... Done")
write("... Done", "log", append = T)

###################################################################################################
### 6. Feature selection
###################################################################################################
print("... Running 6. Feature selection")
write("... Running 6. Feature selection", "log", append = T)
fs <- rfcv(feature_bf_fs, cate_bf_fs, cv.fold = 5, scale = "log", step = 0.9, recursive = F)
len <- length(fs$n.var[fs$error.cv == min(fs$error.cv)])
min_fs <- fs$n.var[fs$error.cv == min(fs$error.cv)][len] #get least features
ind <- order(-imp_bf_fs)[1: min_fs]
feature_fs <- feature_bf_fs[ , ind]
cate_fs <- cate_bf_fs
print("... Writing genes selected during random forest feature selection to 05.Selected-features.txt")
write("... Writing genes selected during random forest feature selection to 05.Selected-features.txt", "log", append = T)
write.table(colnames(feature_fs), "05.Selected-features.txt", quote = F, col.names = F, row.names = F, sep = "\t") 
write(paste("Number of genes selected by random forest feature selection | ", min_fs), "01.Report.txt", append = T)
print("... Done")
write("... Done", "log", append = T)

###################################################################################################
### 7. Selecting training set with random forest prediction internal vote probabilities > 0.6 
###################################################################################################
print("... Running 7. Selecting training set with random forest prediction internal vote probabilities > 0.6")
write("... Running 7. Selecting training set with random forest prediction internal vote probabilities > 0.6", "log", append = T)
rf_fs <- randomForest(feature_fs, as.factor(cate_fs), importance=TRUE, proximity=TRUE)
rf_fs
fea1_fs <- data.frame()
fea1_fs <- feature_fs[(rf_fs$predicted == 'a') & (rf_fs$votes[ , 1] > 0.6), , drop = FALSE]
cat1_fs <- rf_fs$predicted[(rf_fs$predicted =='a') & (rf_fs$votes[ , 1] > 0.6)]
fea2_fs <- data.frame()
fea2_fs <- feature_fs[(rf_fs$predicted == 'b') & (rf_fs$votes[ , 2] > 0.6), , drop = FALSE]
cat2_fs <- rf_fs$predicted[(rf_fs$predicted =='b') & (rf_fs$votes[ , 2] > 0.6)]
cate <- as.factor(c(as.character(cat1_fs), as.character(cat2_fs)))
feature <- as.matrix(rbind(fea1_fs, fea2_fs))
set <- sample(1: nrow(feature), nrow(feature), replace = F)
cate <- cate[set]
feature <- feature[set, ] 
print("... Done")
write("... Done", "log", append = T)

###################################################################################################
### 8. 100 runs of random forest 2-fold cross validation
###################################################################################################
print("... Running 8. 100 runs of random forest 10-fold cross validation")
write("... Running 8. 100 runs of random forest 10-fold cross validation", "log", append = T)
err_fs <- vector()
cate_table <- data.frame()
k <- 10 # 10-fold cross validation
for (run in 1: 100) {
  cate_table <- rbind(cate_table, table(cate))
  n <- floor(nrow(feature) / k)
  subset <- sample(1: nrow(feature), n, replace = F)
  train_feature <- feature[-subset, ]
  train_cate <- cate[-subset]
  test_feature <- feature[subset, ]
  test_cate <- cate[subset]
  retry <- -1
  # Ensure each training set has > 5 samples
  if (table(train_cate)[1] < 5 || table(train_cate)[2] < 5) {
    retry <- 100
    while(retry > 0){
      subset <- sample(1: nrow(feature), n, replace = F)
      train_feature <- feature[-subset, ]
      train_cate <- cate[-subset]
      test_feature <- feature[subset, ]
      test_cate<- cate[subset]
      if (table(train_cate)[1] < 5 || table(train_cate)[2] < 5) {
        retry <- retry - 1
      }  else {
        retry <- -1
      }
    }
  }
  if (retry == 0) {
    next
  }
  if(retry == -1){
    rf <- randomForest(train_feature, as.factor(train_cate), importance = TRUE, proximity = TRUE)
    pred2_fs <- predict(rf, newdata = test_feature)
    mis <- length(test_cate[test_cate != pred2_fs]) / length(test_cate)
    err_fs <- c(err_fs, mis)
    
    cate <- as.factor(c(as.character(train_cate), as.character(pred2_fs)))
    feature <- as.matrix(rbind(train_feature, test_feature))
    set <- sample(1: nrow(feature), nrow(feature), replace=F)
    cate <- cate[set]
    feature <- feature[set, ]
  }
}
print("... Done")
write("... Done", "log", append = T)

###################################################################################################
### 9. Finalizing clusters from classification step
###################################################################################################
print("... Running 9. Finalizing clusters from classification step")
write("... Running 9. Finalizing clusters from classification step", "log", append = T)
rf_whole <- randomForest(feature, as.factor(cate), importance = TRUE, proximity = TRUE)
rf_whole
pred_whole <- predict(rf_whole, newdata = feature_fs)
# pred_whole_prob <- predict(rf_whole, newdata = feature_fs, type = "prob")
# fea1 <- data.frame()
# fea2 <- data.frame()
# fea3 <- data.frame()
# cat1 <- vector()
# cat2 <- vector()
# cat3 <- vector()
# cate_whole <- vector()
# for (i in 1: length(pred_whole)) {
#   if ((pred_whole[i] == 'a') && (pred_whole_prob[i, 1] > 0.55)){
#     cate_whole[i] <- 'a'
#   }  else if ((pred_whole[i] == 'a') && (pred_whole_prob[i, 1] <= 0.55)) {
#     cate_whole[i]<-'c'
#   }  else if ((pred_whole[i] == 'b') && (pred_whole_prob[i, 2] > 0.55)) {
#     cate_whole[i] <- 'b'
#   }  else if ((pred_whole[i] == 'b') && (pred_whole_prob[i, 2] <= 0.55)){
#     cate_whole[i] <- 'c'
#   }
# }
cat1 <- as.factor(as.character(pred_whole[pred_whole == 'a']))
cat2 <- as.factor(as.character(pred_whole[pred_whole == 'b']))
#cat3 <- as.factor(as.character(cate_whole[cate_whole == 'c']))
fea1 <- feature_fs[pred_whole == 'a', , drop = FALSE]
fea2 <- feature_fs[pred_whole == 'b', , drop = FALSE]
#fea3 <- feature_fs[cate_whole == 'c', , drop = FALSE]

# sigval <- sigclust::sigclust(t(feature_fs),
#                              labflag = 1,
#                              label=as.numeric(as.factor(pred_whole)),
#                              icovest=2)
# 
# 
# if(sigval@pval < 0.01){
#   
# }

# Output classification stats
pdf("06.Classification_stat.pdf")
par(mfrow = c(2, 2))
with(fs, plot(n.var, error.cv, log = "x", type = "o", lwd = 2, xlab = "Number of features", ylab = "Error of cross-valicatoin", main = "Feature selections"))
plot(err_fs, col = "blue", main = "Misclassification_error_fs", xlab = "Run", ylab = "Error")
lines(err_fs, col = "blue", lwd = 2)
plot(cate_table[ , 1], col = "magenta", main = "Number of cases in\nClass 1", xlab = "Runs", ylab = "Number of samples")
lines(cate_table[ ,1], col = "magenta", lwd = 2)
plot(cate_table[ , 2], col = "cyan", main = "Number of cases in\nClass 2", xlab = "Runs", ylab = "Number of samples")
lines(cate_table[ , 2],col = "cyan", lwd = 2)
dev.off()

print("... Plotting classification stats to 06.Classification_stat.pdf")
write("... Plotting classification stats to 06.Classification_stat.pdf", "log", append = T)
write(paste("Cluster sizes at CLASSIFICATION step | ", length(cat1), ", ", length(cat2)), "01.Report.txt", append = T)
#write(paste("Number of outliers | ", length(cat3)), "01.Report.txt", append = T)

# Output expression stats and clusters
gene_name <- colnames(gene_mat)
fea1_name <- row.names(fea1)
fea2_name <- row.names(fea2)
index_fea1 <- vector()
index_fea2 <- vector()
for (i in 1: length(gene_name)) {
  for (j in 1:length(fea1_name)) {
    if (gene_name[i] == fea1_name[j]) {
      index_fea1 <- c(index_fea1, i)
    }
  }
  for (k in 1: length(fea2_name)) {
    if (gene_name[i] == fea2_name[k]) {
      index_fea2 <- c(index_fea2, i)
    }
  }
}
cl_1 <- gene_mat[ , index_fea1]
cl_2 <- gene_mat[ , index_fea2]

temp <- gene_mat[ , c(index_fea1, index_fea2)]

markers <- read.csv("../raw_data/20170323_Rubenstein_revised_gene_list.csv", header = F, stringsAsFactors = F)
markers <- unlist(markers[ ,1 ])
markers <- unique(sapply(markers, function(x){paste("\\b", x, "$", sep = "")}))
markers <- unlist(sapply(markers,function(x){grep(x, rownames(temp), value = T)}))

temp <- temp[ markers, ]

colsep <- length(index_fea1);

pairs.breaks <- seq(-4, 4, length.out=101);

pdf(file = paste(Sys.Date(), args[1], args[2], "marker_heatmap.pdf", sep = "_"),
    width = 20, height = 15, useDingbats = F);
heatmap.3(temp,
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8,
          main="Marker gene expression",
          col = bluered(100),
          symkey = F,
          cexRow=1, cexCol = 0.6,
          Rowv = F,
          Colv = F,#as.dendrogram(hca),
          #ColSideColors = as.matrix(ColSideColors),
          #ColSideColorsSize = 2,
          #RowSideColors = RowSideColors,
          dendrogram = "none",
          scale = "row",
          colsep = colsep,
          sepcolor = "black",
          labRow = substr(rownames(temp), start = 20, stop = 100), #c("Glast", "Tbr2"),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

print("... Writing clusters with associated samples to 07.Cluster1/2_sample-list.txt")
write("... Writing clusters with associated samples to 07.Cluster1/2_sample-list.txt", "log", append = T)
write.table(colnames(cl_1), "Round1.C1.sample-list.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(colnames(cl_2), "Round1.C2.sample-list.txt", quote = F, sep = "\t", row.names = F, col.names = F)
print("... Writing gene expression statistics to 08.All_genes_stat.txt")
write("... Writing gene expression statistics to 08.All_genes_stat.txt", "log", append = T)
# get_fs_deg <- get_deg(cl_1, cl_2, 2)
# if (length(get_fs_deg[[7]]) > 0) {
#   out <- cbind(rownames(cl_1), round(get_fs_deg[[3]], digits = 3), round(get_fs_deg[[4]], digits = 3), round(get_fs_deg[[5]], digits = 3), format(get_fs_deg[[6]], scientific = T, digits = 3), format(get_fs_deg[[7]], scientific = T, digits = 3))
#   write.table(out, "08.All_genes_stat.txt", quote = F, col.names = c("Gene", "Mean_Cluster1", "Mean_Cluster2", "log2FC", "Pvalue", "P_adj"), row.names = F, sep = "\t")
# } else {
#   out <- cbind(rownames(cl_1), round(get_fs_deg[[3]], digits = 3), round(get_fs_deg[[4]], digits = 3), round(get_fs_deg[[5]], digits = 3), format(get_fs_deg[[6]], scientific = T, digits = 3))
#   write.table(out, "08.All_genes_stat.txt", quote = F, col.names = c("Gene", "Mean_Cluster1", "Mean_Cluster2", "log2FC", "Pvalue"), row.names = F, sep = "\t")
# }
# write(paste("Total number DEGs with fold-change > 2 | ", (length(get_fs_deg[[1]]) + length(get_fs_deg[[2]]))), "01.Report.txt", append = T)
# write(paste("DEGs over-expressed in cluster 1 | ", length(get_fs_deg[[1]])), "01.Report.txt", append = T)
# write(paste("DEGs over-expressed in cluster 2 | ", length(get_fs_deg[[2]])), "01.Report.txt", append = T)

gene_mat0 <- round(cl_1, digits = 3)
save(gene_mat0, #DI.prior, 
     file =  paste("Round1.C1.AllGenes_NewInput.dat", sep = ""))
gene_mat0 <- round(cl_2, digits = 3)
save(gene_mat0, #DI.prior, 
     file =  paste("Round1.C2.AllGenes_NewInput.dat", sep = ""))

print("... Writing input files for the next round of Clustering-and-Classification to 09.Cluster1/21_AllGenes_NewInput.dat")
write("... Writing input files for the next round of Clustering-and-Classification to 09.Cluster1/21_AllGenes_NewInput.dat", "log", append = T)
print("... Done")
write("... Done", "log", append = T)

# Generating heatmap using selected-genes
#index_2fc <- c(get_fs_deg[[1]], get_fs_deg[[2]])
feature_names <- colnames(feature)
m_final <- as.matrix(cbind(cl_1[feature_names, ], cl_2[feature_names, ]))
hr_m_final <- hclust(as.dist(1 - cor(t(m_final), method = "pearson", use = "pairwise.complete.obs")), method = clust_method)
colsep <- ncol(cl_1)

pdf(paste(args[1], args[2], "features_Heatmap.pdf", sep = "_"))
heatmap.3(m_final, trace = "none", 
          breaks = pairs.breaks,
          Rowv = as.dendrogram(hr_m_final), 
          dendrogram = "row", 
          Colv = F, 
          colsep = ncol(cl_1), sepwidth = 0.1, sepcolor = "green", 
          col = colorRampPalette(c(myCol)), 
          #ColSideColors = c(rep("red", ncol(cl_1)), 
          #                  rep("blue", ncol(cl_2))), 
          main = paste(name, "Heatmp of clusters", sep = ""), 
          cexCol = 0.4, 
          cexRow = 0.4, 
          labRow = substr(rownames(m_final), start = 20, stop = 100))
#legend("topright", legend = c("Cluster 1", "Cluster 2"), fill = c("red", "blue"), border="white", bty="n", cex=1)
dev.off()
print("... Plotting heatmap using selected-genes to 10.Heatmap_features-selected.pdf")
write("... Plotting heatmap using selected-genes to 10.Heatmap_features-selected.pdf", "log", append = T)
print(paste(rep("#", 80), collapse = ""))
write(paste(rep("#", 80), collapse = ""), "log", append = T)
print("ALL DONE!")
write("ALL DONE!", "log", append = T)
print(paste(rep("#", 80), collapse = ""))
write(paste(rep("#", 80), collapse = ""), "log", append = T)
write(paste(rep("#", 100), collapse = ""), "01.Report.txt", append = T)
write(format(Sys.time(), "%Y-%m-%d %I-%p"), "01.Report.txt", append = T)
write(format(Sys.time(), "%Y-%m-%d %I-%p"), "log", append = T)
write("", "01.Report.txt", append = T)
write("", "log", append = T)
###END###
