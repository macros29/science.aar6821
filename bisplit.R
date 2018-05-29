setwd("~/Data/Da_Mi/")
load("working.RData")

set.seed(7)
keep <- GB$region == "dMGE" & GB$time_point == "E12.5" & !(GB$IFC %in% c("189","180","181","190","191","192","193","194","195", "196"))
df <- CreateSeuratObject(raw.data = as.vector(expr[,keep]), min.cells = 5, min.genes = 1000, project = "Da_Mi")
df <- NormalizeData(object = df, normalization.method = "LogNormalize", scale.factor = 10000)
df <- ScaleData(df)
df@meta.data <- cbind(df@meta.data, GB[df@cell.names,])
df <- FindVariableGenes(object = df, mean.function = ExpMean, dispersion.function = LogVMR, 
                        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5,do.plot = F)
df <- RunPCA(object = df, pc.genes=df@var.genes, pcs.compute = 30)
df <- RunTSNE(object = df, dims.use = 1:3, do.fast = TRUE, seed.use = 7)
df <- FindClusters(object = df,
                   reduction.type = "pca",
                   dims.use = 1:3, resolution = 0.9, 
                   print.output = 0, 
                   save.SNN = TRUE, 
                   force.recalc = T,
                   random.seed = 1)

pdf(file = paste("./E12.5_dMGE/", Sys.Date(), "_tSNE_seed_", seed.use, "_E12.5_dMGE.pdf", sep = ""), 
    width = 11, height = 10, useDingbats = F);
TSNEPlot(object = df, pt.size = 5)
dev.off() 

sigval <- sigclust::sigclust(t(df@scale.data[df@var.genes,]),
                             50,
                             labflag = 0,
                             #label=as.numeric(as.factor(df@meta.data$cc_clusters)),
                             icovest=2)

#### Random forest step ####
library(randomForest)
cate_bf_fs <- df@ident
feature_bf_fs <- as.matrix(t(df@scale.data[df@var.genes,]))
rf_bf_fs <- randomForest(feature_bf_fs, cate_bf_fs, importance = TRUE, proximity = TRUE, ntree = 2000,
                         mtry = 251)
imp_bf_fs <- importance(rf_bf_fs, type = 1)

temp <- df
temp@ident <- rf_bf_fs$predicted

pdf(file = paste("./E12.5_dMGE/", Sys.Date(), "_tSNE_seed_", seed.use, "_E12.5_dMGE.pdf", sep = ""), 
    width = 11, height = 10, useDingbats = F);
TSNEPlot(object = temp, pt.size = 5)
dev.off() 

fs <- rfcv(feature_bf_fs, cate_bf_fs, cv.fold = 10, scale = "log", step = 0.7, recursive = F)
len <- length(fs$n.var[fs$error.cv == min(fs$error.cv)])
min_fs <- fs$n.var[fs$error.cv == min(fs$error.cv)][len] #get least features
ind <- order(-imp_bf_fs)[1: min_fs]
feature_fs <- feature_bf_fs[ , ind]
cate_fs <- cate_bf_fs

rf_fs <- randomForest(feature_fs, as.factor(cate_fs), importance=TRUE, proximity=TRUE)

temp <- CreateSeuratObject(raw.data = t(feature_fs), min.cells = 0, min.genes = 0, project = "Da_Mi")
temp@scale.data <- temp@raw.data
temp@meta.data <- cbind(temp@meta.data, GB[temp@cell.names,])
temp <- FindVariableGenes(object = temp, mean.function = ExpMean, dispersion.function = LogVMR, 
                        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5,do.plot = F)
temp <- RunPCA(object = temp, pc.genes=rownames(temp@scale.data), pcs.compute = 30)
temp <- RunTSNE(object = temp, dims.use = 1:3, do.fast = TRUE, seed.use = 7)

temp <- df
temp@ident <- rf_whole$predicted

pdf(file = paste("./E12.5_dMGE/", Sys.Date(), "_tSNE_seed_", seed.use, "_E12.5_dMGE.pdf", sep = ""), 
                                 width = 11, height = 10, useDingbats = F);
TSNEPlot(object = temp, pt.size = 5)
dev.off() 

feature <- c()
cate <- c()
pred <- sort(unique(rf_fs$predicted))
for(i in 1:length(pred)){
  fea_fs <- feature_fs[(rf_fs$predicted == pred[i]) & (rf_fs$votes[ , i] > 0.5), , drop = FALSE]
  feature <- as.matrix(rbind(feature, fea_fs))
  cat_fs <- rf_fs$predicted[(rf_fs$predicted == pred[i]) & (rf_fs$votes[ , i] > 0.5)]
  cate <- as.factor(c(as.character(cate), as.character(cat_fs)))
}

rf_whole <- randomForest(feature, as.factor(cate), importance = TRUE, proximity = TRUE)
pred_whole <- predict(rf_whole, newdata = feature_fs)

temp@ident <- pred_whole

pdf(file = paste("./E12.5_dMGE/", Sys.Date(), "_tSNE_seed_", seed.use, "_E12.5_dMGE.pdf", sep = ""), 
    width = 11, height = 10, useDingbats = F);
TSNEPlot(object = temp, pt.size = 5)
dev.off() 

err_fs <- vector()
cate_table <- data.frame()
k <- 1 # 1-fold cross validation
for (run in 1:100) {
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

temp <- df
temp@ident <- fs$predicted
TSNEPlot(object = temp, pt.size = 5)

#### Work on cluster 0 in E12.5 dMGE ####
keep <- GB$New_Cell_ID %in% df@cell.names[df@ident == 0]
df0 <- CreateSeuratObject(raw.data = as.vector(expr[,keep]), min.cells = 5, min.genes = 1000, project = "Da_Mi")
df0 <- NormalizeData(object = df0, normalization.method = "LogNormalize", scale.factor = 10000)
df0 <- ScaleData(df0)
df0 <- FindVariableGenes(object = df0, mean.function = ExpMean, dispersion.function = LogVMR, 
                        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5,do.plot = F)
df0 <- RunPCA(object = df0, pc.genes=df0@var.genes, pcs.compute = 30)
df0 <- FindClusters(object = df0,
                   reduction.type = "pca",
                   dims.use = 1:3, resolution = 0.9, 
                   print.output = 0, 
                   save.SNN = TRUE, 
                   force.recalc = T)
df0.markers <- FindAllMarkers(object = df0, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
df0.top20 <- df0.markers %>% group_by(cluster) %>% top_n(20, avg_diff)

pdf(file = paste("./E12.5_dMGE/", Sys.Date(), "_E12.5_dMGE_df0.1_markers.pdf", sep = ""), 
    width = 10, height = 10, useDingbats = F);
DoHeatmap(object = df0, 
          genes.use = df0.top20$gene, 
          order.by.ident = T, 
          slim.col.label = TRUE,group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

write.csv(df0.markers, file = "2017-08-31_DE_E12.5_dMGE_df0.1.csv")

df0@ident <- as.factor(as.numeric(df0@ident) + length(unique(df@ident)) - 1)
names(df0@ident) <- df0@cell.names

res <- as.numeric(df@ident)
res[df@cell.names %in% df0@cell.names] <- as.numeric(df0@ident) + length(unique(df@ident))
res <- res - 1
df@ident <- as.factor(res)
names(df@ident) <- df@cell.names

pdf(file = paste("./E12.5_dMGE/", Sys.Date(), "_tSNE_seed_", seed.use, "_E12.5_dMGE_df0.1.pdf", sep = ""), 
    width = 11, height = 10, useDingbats = F);
TSNEPlot(object = df, pt.size = 5)
dev.off() 

df.markers <- FindAllMarkers(object = df, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

## Plot top 20 marker genes
top20 <- df.markers %>% group_by(cluster) %>% top_n(20, avg_diff)

pdf(file = paste("./E12.5_dMGE/", Sys.Date(), "_tSNE_E12.5_dMGE_df0.1_all_clusters.pdf", sep = ""), 
    width = 20, height = 40, useDingbats = F);
DoHeatmap(object = df, 
          genes.use = top20$gene, 
          order.by.ident = T, 
          slim.col.label = TRUE,group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

#### 
keep <- GB$New_Cell_ID %in% df@cell.names[df@ident == 1]

df1 <- CreateSeuratObject(raw.data = as.vector(expr[,keep]), min.cells = 5, min.genes = 1000, project = "Da_Mi")
df1 <- NormalizeData(object = df1, normalization.method = "LogNormalize", scale.factor = 10000)
df1 <- ScaleData(df1)
df1 <- FindVariableGenes(object = df1, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = F)
m <- as.matrix(df1@data[df1@var.genes,])

df1 <- RunPCA(object = df1, pc.genes=df1@var.genes, pcs.compute = 50)

#### Do tSNE ####
library(Rtsne)
seed <- 1024
set.seed(seed)
p = 20
m = 0.1
fm = 0.1
Data <- df1@scale.data[df1@var.genes,]
t_sne <- Rtsne(t(Data),perplexity = p, pca = T, initial_dims = 3, verbose = F, theta = 0, momentum = m, 
               final_momentum = fm,
               #Y_init = as.matrix(t_sne),
               max_iter = 10000)$Y
colnames(t_sne) <- c("X1", "X2")
rownames(t_sne) <- colnames(Data)
t_sne <- data.frame(t_sne)

#load(paste("2017-05-29", name, "tSNE.Rdata", sep = "_"))

pdf(file = paste(Sys.Date(), name,"seed",seed, "p", p, "m", m,"fm",fm,
                 "tSNE.pdf", sep = "_"), 
    width = 12, height = 10, useDingbats = F)
g <- ggplot(data = t_sne, aes(x = X1, y = X2)) +
  geom_point(size = 5, aes(x = X1, y = X2, color = df1@ident)) +
  #scale_color_brewer(palette = "Paired", name = 'Ages') +
  scale_color_discrete(name = "cell_type") +
  #geom_text(aes(label = colnames(Y2)), hjust = 1, vjust = -1) + 
  labs(x = "tSNE 1", y = "tSNE 2") +
  #ggtitle(paste("All cells, Perplexity =", p)) +
  theme_bw() +
  theme(text = element_text(size = 20),
        line = element_line(size = 1),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.key = element_blank(),
        legend.background = element_blank())
print(g)
dev.off(which = dev.cur())

df1 <- RunTSNE(object = df1, dims.use = 1:3, do.fast = TRUE,)

# note that you can set do.label=T to help label individual clusters

set.seed(1)
df1 <- FindClusters(object = df1,
                    reduction.type = "pca",
                    dims.use = 1:3, resolution = 0.9, 
                    print.output = 0, 
                    save.SNN = TRUE, 
                    force.recalc = T)

df1.markers <- FindAllMarkers(object = df1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

top20 <- df1.markers %>% group_by(cluster) %>% top_n(20, avg_diff)

DoHeatmap(object = df1, genes.use = top20$gene, order.by.ident = T, slim.col.label = TRUE, 
          remove.key = TRUE)

TSNEPlot(object = df1)

write.csv(df1.markers, file = "2017-08-31_DE_E12.5_dMGE_df1.csv")

df1@ident <- as.factor(as.numeric(df1@ident) + length(unique(df@ident)) - 1)
names(df1@ident) <- df1@cell.names

res <- as.numeric(df@ident)
res[df@cell.names %in% df1@cell.names] <- as.numeric(df1@ident) + length(unique(df@ident))
res <- res - 1
df@ident <- as.factor(res)
names(df@ident) <- df@cell.names

