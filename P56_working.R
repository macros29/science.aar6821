## Analyze P56 Adult data ##
## Project: Da Mi

#### 1. Load libraries and set up working directory ####
library("matrixStats")
library("reshape")
library("gplots")
library("ggplot2")
library("GGally")
source("~/Scripts/R/heatmap.3.R")
library("Rtsne")
#source("~/Scripts/R/SNN-Cliq-Gap.R")
source("~/Scripts/R/gg_color_hue.R")
library("RColorBrewer")
library("locfit")
library("MASS")
library("hydroGOF")
library("calibrate")
library("pcaReduce")
source("~/Scripts/R/scVEGs.R");
options(stringAsFactors = F);

setwd("~/Data/Da_Mi/")

#### 2. Load data ####
counts.allen <- read.csv("raw_data/GSE71585_RefSeq_counts.csv", header = T, stringsAsFactors = F)
rownames(counts.allen) <- counts.allen[,1]
counts.allen <- counts.allen[,-1]
GB.allen <- read.csv("raw_data/GSE71585_Clustering_Results.csv", header = T, stringsAsFactors = F)
GB.allen$sample_title <- gsub("\\-", "\\.", GB.allen$sample_title)
GB.allen$IN_broad_type <- substr(GB.allen$primary_type, 1, 4) 
counts.allen <- counts.allen[,GB.allen$sample_title]

keep <- GB.allen$broad_type == "GABA-ergic Neuron"
RP.genes <- grep("Rp", rownames(counts.allen))

df0.allen <- CreateSeuratObject(raw.data = as.vector(counts.allen[-c(RP.genes),keep]), 
                          min.cells = 5, min.genes = 1000, project = "Da_Mi")
df0.allen@meta.data <- cbind(df0.allen@meta.data, GB.allen[keep,])
df0.allen <- NormalizeData(object = df0.allen, normalization.method = "LogNormalize", 
                     scale.factor = 1000000)

df0.allen <- FindVariableGenes(object = df0.allen, mean.function = ExpMean, 
                         dispersion.function = LogVMR, 
                         x.low.cutoff = 0.5, 
                         x.high.cutoff = Inf, 
                         y.cutoff = 0.5, 
                         y.high.cutoff = 5, 
                         do.plot = T)

df0.allen <- ScaleData(df0.allen)

df0.allen <- RunPCA(object = df0.allen, pc.genes=df0.allen@var.genes, 
                    pcs.compute = 30, pcs.print = 1:10,
              genes.print = 10)

p1 <- PCAPlot(df0.allen, group.by = "primary_type",do.return = T)
p2 <- PCAPlot(df0.allen, group.by = "secondary_type", do.return = T)
p3 <- PCAPlot(df0.allen, group.by = "IN_broad_type",do.return = T)
pdf(paste(Sys.Date(), "P56_All_cells_PCA.pdf", sep = "_"), width = 26, height = 8, useDingbats = F);
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

#### 3. Run random forest ####
# cate_bf_fs <- as.factor(df0.allen@meta.data$primary_type)

## Find highly variable genes that were shared by P56 and Embryonic interneurons
# load("2017-11-22_df0_neu_Seurat_object.Rdata")
# P56.genes <- toupper(df0.allen@var.genes)
# Da.genes <- substr(df0.neu@var.genes, 20, 100)

## Only use shared highly variable genes to for Random Forest classification of P56 lineages
# feature_bf_fs <- as.matrix(t(df0.allen@scale.data[df0.allen@var.genes[P56.genes %in% Da.genes],]))
# 
# rf_bf_fs <- randomForest(feature_bf_fs, cate_bf_fs, importance = TRUE, proximity = TRUE)
# imp_bf_fs <- importance(rf_bf_fs, type = 1)
# 
# fs <- rfcv(feature_bf_fs, cate_bf_fs, cv.fold = 10, scale = "log", step = 0.9, recursive = F)
# len <- length(fs$n.var[fs$error.cv == min(fs$error.cv)])
# min_fs <- fs$n.var[fs$error.cv == min(fs$error.cv)][1] #get least features
# ind <- order(-imp_bf_fs)[1: min_fs]
# feature_fs <- feature_bf_fs[ , ind]
# cate_fs <- cate_bf_fs
# 
# temp <- t(feature_fs)
# 
# ord <- order(cate_fs)
# ColSideColors <- c(#gg_color_hue(3)[as.numeric(as.factor(df0.E12.5.prog@meta.data$region))],
#   gg_color_hue(23)[as.numeric(cate_fs)])#,
# ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
# ColSideColors <- as.matrix(ColSideColors[ord,])
# col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)
# 
# pairs.breaks <- seq(-1.5, 1.5, length.out=1001);
# 
# pdf(file = paste(Sys.Date(), "P56_int_heatmap_RF.pdf", sep = "_"),
#     width = 10, height = 14, useDingbats = F);
# heatmap.3(temp[,ord],
#           breaks = pairs.breaks,
#           #symbreaks = T,
#           keysize = 0.8,
#           key = F,
#           main=NA,
#           col = col,
#           symkey = F,
#           cexRow=0.5, 
#           cexCol = 0.6, 
#           Rowv = T, 
#           Colv = F, #as.dendrogram(clu), 
#           ColSideColors = ColSideColors,
#           ColSideColorsSize = 0.5,
#           dendrogram = "both",
#           #scale = "row",
#           #colsep = colsep,
#           sepcolor = "black",
#           labRow = substr(rownames(temp), 20, 100),
#           labCol = "",
#           na.rm = F);
# dev.off(dev.cur());
# 
# P56.RF.features <- colnames(feature_fs)
# write.csv(P56.RF.features, file = "P56.RF.features.csv")

P56.RF.features <- read.csv("2017-11-22_P56.RF.features.csv", header = T, row.names = 1)

feature_fs <- as.matrix(t(df0.allen@scale.data[as.character(P56.RF.features$x),]))
cate_fs <- as.factor(df0.allen@meta.data$primary_type)

rf_fs_P56 <- randomForest(feature_fs, cate_fs, importance=TRUE, proximity=TRUE)
save(rf_fs_P56, file = "P56.RF.model.Rdata")

# cate <- apply(rf_fs$votes, 1, function(xx){
#   colnames(rf_fs$votes)[xx == max(xx)]
# })
# cate <- as.factor(unlist(cate))
# cate <- cate[rowMax(rf_fs$votes) > 0.6]
# feature <- feature_fs[rowMax(rf_fs$votes) > 0.6,]
# 
# cate <- as.factor(c(as.character(cat1_fs), as.character(cat2_fs),as.character(cat3_fs)))
# feature <- as.matrix(rbind(fea1_fs, fea2_fs, fea3_fs))
# 
# set <- sample(1: nrow(feature), nrow(feature), replace = F)
# cate <- cate[set]
# feature <- feature[set, ] 
# 
# rf_whole <- randomForest(feature, as.factor(cate), importance = TRUE, proximity = TRUE)
# pred_whole <- predict(rf_whole, newdata = feature_fs)

# temp <- t(feature_fs)
# # dis <- as.dist(1 - cor(temp,use = "pairwise.complete.obs"))
# # clu <- hclust(dis, method = "complete")
# # cl <- cutree(clu, k=2)
# 
# ord <- order(pred_whole)
# ColSideColors <- c(#gg_color_hue(3)[as.numeric(as.factor(df0.E12.5.prog@meta.data$region))],
#   brewer.pal(3, name = "Paired")[as.numeric(pred_whole)])#,
# ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
# ColSideColors <- as.matrix(ColSideColors[ord,])
# col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)
# 
# pairs.breaks <- seq(-1.5, 1.5, length.out=1001);
# 
# pdf(file = paste(Sys.Date(), "cortex_int_heatmap_post-RF.pdf", sep = "_"),
#     width = 10, height = 14, useDingbats = F);
# heatmap.3(temp[,ord],
#           breaks = pairs.breaks,
#           #symbreaks = T,
#           keysize = 0.8,
#           key = F,
#           main=NA,
#           col = col,
#           symkey = F,
#           cexRow=0.5, 
#           cexCol = 0.6, 
#           Rowv = T, 
#           Colv = F, #as.dendrogram(clu), 
#           ColSideColors = ColSideColors,
#           ColSideColorsSize = 0.5,
#           dendrogram = "both",
#           #scale = "row",
#           #colsep = colsep,
#           sepcolor = "black",
#           labRow = substr(rownames(temp), 20, 100),
#           labCol = "",
#           na.rm = F);
# dev.off(dev.cur());

#### 5. Do CCA on Allen and Da Mi cells ####
genes.use <- union(substr(df0.neu@var.genes, 20, 100), 
                       toupper(df0.allen@var.genes))

#### 6.1 Between Da Mi and P56 datasets ####
dt1 <- df0.neu@raw.data[substr(rownames(df0.neu@raw.data), 20, 100) %in% genes.use, 
                                     df0.neu@cell.names]
rownames(dt1) <- make.names(substr(rownames(dt1), 20, 100), unique = T)
dt1 <- CreateSeuratObject(raw.data = dt1, 
                          min.cells = 0, min.genes = 0, project = "Da_Mi")
dt1 <- NormalizeData(object = dt1, normalization.method = "LogNormalize", 
                     scale.factor = 1000000)
lineage <- "Unknown"
dt1@meta.data$lineage <- lineage
dt1 <- ScaleData(dt1, vars.to.regress = c("nGene"))

dt2 <- df0.allen@raw.data[ toupper(rownames(df0.allen@raw.data)) %in% genes.use, ]
rownames(dt2) <- toupper(rownames(dt2))
dt2 <- CreateSeuratObject(raw.data = dt2, 
                          min.cells = 0, min.genes = 0, project = "P56")
primary_type <- paste("P56", df0.allen@meta.data$primary_type)
dt2@meta.data$lineage <- primary_type
dt2 <- NormalizeData(object = dt2, normalization.method = "LogNormalize", 
                     scale.factor = 1000000)
dt2@meta.data$orig.ident <- "P56"
dt2 <- ScaleData(dt2, vars.to.regress = "nGene")

dt2 <- SetAllIdent(dt2, "lineage")
P56.markers <- FindAllMarkers(dt2, 
                              genes.use = genes.use, 
                               thresh.use = 0.01, min.pct = 0.25, 
                               only.pos = T,
                              do.print = T)

pdf(file = paste(Sys.Date(),"P56_markers_heatmap.pdf", sep = "_"),
    width = 20, height = 40, useDingbats = F);
DoHeatmap(object = dt2,
          genes.use = P56.markers$gene,#rownames(dff.markers)[1:100],
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()
# dt2 <- SubsetData(dt2, ident.remove = c("P56 Smad3", "P56 Igtp", "P56 Sncg"))
# dt2 <- SubsetData(dt2, ident.use = c("P56 Pvalb Gpx3", "P56 Sst Cbln4", "P56 Sst Myh8",
#                                         "P56 Vip Chat", "P56 Vip Parm1",
#                                         "P56 Ndnf Cxcl14"))

# genes.use <- intersect(substr(df0.neu.lin.NotChodl@var.genes, 20, 100),
#                        toupper(df0.allen@var.genes))#c(toupper(P56.RF.features[,1]))

P56.RF.features

dt <- RunCCA(object = dt1, object2 = dt2, genes.use = toupper(P56.RF.features[,1]), 
             scale.data = T)

# p1 <- DimPlot(object = dt, reduction.use = "cca", group.by = "orig.ident", pt.size = 0.5, 
#               do.return = TRUE)
# p2 <- VlnPlot(object = dt, features.plot = "CC1", group.by = "orig.ident", do.return = TRUE)
# plot_grid(p1, p2)

pdf(file = paste(Sys.Date(), "CCA_Heatmaps.pdf", sep = "_"),
    width = 20, height = 30, useDingbats = F);
DimHeatmap(object = dt, reduction.type = "cca", dim.use = 1:20, 
           do.balanced = TRUE)
dev.off()

dt <- AlignSubspace(dt, reduction.type = "cca", grouping.var = "orig.ident", 
                    dims.align = 1:7)

p1 <- VlnPlot(object = dt, features.plot = "ACC1", group.by = "orig.ident", 
              do.return = TRUE)
p2 <- VlnPlot(object = dt, features.plot = "ACC2", group.by = "orig.ident", 
              do.return = TRUE)
plot_grid(p1, p2)

dt@meta.data$time_point <- substr(dt@meta.data$lineage, 1, 3)
dt@meta.data$time_point[1:1019] <- df0.neu@meta.data$time_point
dt@meta.data$RF_cate_whole <- dt@meta.data$lineage
dt@meta.data$RF_cate_whole[1:1019] <- paste("Emb", df0.neu@meta.data$RF_cate_whole)

#### Do KNN analysis ####
library(FNN)

MiDa.CCA <- dt@dr$cca.aligned@cell.embeddings[1:1019,1:3]
Allen.CCA <- dt@dr$cca.aligned@cell.embeddings[1020:1780,1:3]

p = 40
k = 30
dt@meta.data$lineage_knn <- "Uncertain"

dt <- RunTSNE(object = dt, 
              reduction.use = "cca.aligned", 
              dims.use = 1:3,
              perplexity = p,
              do.fast = TRUE,
              seed.use = 1)

MiDa.tsne <- dt@dr$tsne@cell.embeddings[1:1019, 1:2]
Allen.tsne <- dt@dr$tsne@cell.embeddings[1020:1780, 1:2]

temp <- knn.cv(as.matrix(x = Allen.tsne), cl = dt@meta.data$lineage[1020:1780], 
               k = k, prob = T)

res.prob <- attr(temp,"prob")
ave <- mean(attr(temp, "prob"))
sd <- sd(attr(temp, "prob"))

res.dist <- attr(temp, "nn.dist")
res.dist.ave <- mean(res.dist[res.prob > ave - sd,])
res.dist.sd <- sd(res.dist[res.prob > ave - sd,])
remove <- res.dist > (res.dist.ave + res.dist.sd)
res.dist[remove] <- NA

res.index <- attr(temp, "nn.index")
res.index[remove] <- NA

res <- apply(res.index, 1, function(xx){dt@meta.data$lineage[1020:1780][xx]})
res <- apply(res, 2, function(xx){
  if(length(which(!is.na(xx))) < 10){
    "Uncertain"
  } else {
    r <- table(xx)
    if(length(which(r == max(r))) > 1 | max(r) < 10){
      "Uncertain"
    } else {names(sort(r, decreasing = T))[1]}
  }
})

# res <- apply(res, 2, function(xx){
#   if(length(which(!is.na(xx))) < 10){
#     "Uncertain"
#   } else {
#     r <- table(xx)
#     if(length(which(r == max(r))) > 1 | max(r) < 10){
#       "Uncertain"
#     } else {names(sort(table(xx), decreasing = T)[1])}
#   }
# })


res[res.prob < (ave - sd)] <- "Uncertain"

res[res %in% names(which(table(res) <= 10))] <- "Uncertain"

dt@meta.data$lineage_knn[1020:1780] <- res

keep <- !(dt@meta.data$lineage_knn[1020:1780] %in% "Uncertain")

knn <- get.knnx(data <- as.matrix(x = Allen.tsne[keep,]),
                query <- as.matrix(x = MiDa.tsne), k = 5)

## Use distance between each cells to its neighbors as cutoff
res.dist <- knn$nn.dist

## Calculate mean and sd
res.dist.ave <- mean(res.dist)
res.dist.sd <- sd(res.dist)
remove <- res.dist > (res.dist.ave - res.dist.sd)

## Use one sd from mean as cutoff and remove any distance and asignment beyond it
res.dist[remove] <- NA
# res.dist.keep <- apply(res.dist, 1, function(xx){length(which(is.na(xx)))})

res.index <- knn$nn.index
res.index[remove] <- NA

## Find the 
res <- apply(res.index,1, function(xx){dt@meta.data$lineage_knn[1020:1780][keep][xx]})

## Find the number of major lineage assignment of each cell and determine average
# rec <- apply(temp, 2, function(xx){max(table(xx))})
# hist(rec)

## Find average distance 
# ave.d2 <- mean(rowMeans(na.omit(res.dist)))
# sd.d2 <- sd(rowMeans(na.omit(res.dist)))
# thresh.d2 <- rowMeans(na.omit(res.dist)) > ave.d2 + sd.d2

temp2 <- apply(res, 2, function(xx){
  if(length(which(!is.na(xx))) < 2){
    "Uncertain"
    } else {
      r <- table(xx)
      if(length(which(r == max(r))) > 1
         # | max(r) < 3
         ){
        "Uncertain"
      } else {gsub("P56", "Emb", names(sort(r, decreasing = T)[1]))}
    }
})

temp3 <- knn.cv(as.matrix(x = MiDa.tsne[temp2 != "Uncertain",]), 
               cl = temp2[temp2 != "Uncertain"], k = 5, prob = T)
# res.prob <- attr(temp3,"prob")
# ave <- mean(attr(temp3, "prob"))
# sd <- sd(attr(temp3, "prob"))
#  
# res.dist <- attr(temp3, "nn.dist")
# res.dist.ave <- mean(res.dist[res.prob > ave -sd])
# res.dist.sd <- sd(res.dist[res.prob > ave - sd])
# remove <- res.dist > (res.dist.ave)
# res.dist[remove] <- NA
# 
# res.index <- attr(temp3, "nn.index")
# res.index[remove] <- NA
# # 
# res <- apply(res.index, 1, function(xx){temp2[temp2 != "Uncertain"][xx]})
# 
# res <- apply(res, 2, function(r){
#   if(length(which(!is.na(r))) < 3){
#     "Uncertain"
#   } else {names(sort(table(r), decreasing = T))[1]}
# })
# 
# # temp2 <- apply(res, 2, function(xx){
# #   if(length(which(!is.na(xx))) < 9){
# #     "Uncertain"
# #   } else {
# #     r <- table(xx)
# #     
# #   }
# # })
# # 
# # res.dist.na <- apply(res.dist, 1, function(xx){length(which(is.na(xx)))})
# # res <- as.character(temp3)
# 
# res[res.prob < (ave - sd) 
#     # | rowMeans(attr(temp3, "nn.dist")) > (ave.d + sd.d)
#     ] <- "Uncertain"

temp4 <- temp2
temp4[temp2 != "Uncertain"] <- as.character(temp3)

dt@meta.data$lineage_knn[1:1019] <- temp4
# temp <- knn.cv(as.matrix(x = Allen.CCA), cl = temp2[temp2 == "Uncertain"], 
#                k = p, prob = T)
# res <- as.character(temp)
# ave <- mean(attr(temp, "prob"))
# sd <- sd(attr(temp, "prob"))
# 
# ave.d <- mean(rowMeans(attr(temp, "nn.dist"))) ## Average of mean distance of each cell to its 10 nearest neighbors
# sd.d <- sd(rowMeans(attr(temp, "nn.dist")))

p1 <- TSNEPlot(object = dt,colors.use = c("red", "green", "grey"), 
               group.by = "time_point", do.return = TRUE, pt.size = 4)
p1 <- p1 + guides(color = guide_legend(direction = "vertical")) +
  theme(legend.position = c(0.1,0.9),
        legend.background = element_rect(colour = "black"))

p2 <- TSNEPlot(object = dt, cells.use = dt@cell.names[1:1019],
               colors.use = c("grey90",sample(gg_color_hue(3)), "grey90", sample(gg_color_hue(6)[4:6])), 
               group.by = "RF_cate_whole", do.return = TRUE, pt.size = 4)
p2 <- p2 + guides(color = guide_legend(direction = "vertical")) +
  theme(legend.position = c(0.1,0.85),
        legend.background = element_rect(colour = "black"))

p3 <- TSNEPlot(object = dt, cells.use = dt@cell.names[1020:1780], 
               colors.use = c(sample(gg_color_hue(length(unique(dt@meta.data$lineage_knn[1020:1780]))-1)),
                              "grey90"),
               group.by = "lineage_knn", 
               do.return = TRUE, pt.size = 4)
p3 <- p3 + guides(color = guide_legend(direction = "horizontal", ncol = 3)) +
  theme(legend.position = c(0.15,0.9),
        legend.background = element_rect(colour = "black"))

p4 <- TSNEPlot(object = dt, cells.use = dt@cell.names[1:1019], 
               colors.use = c(sample(gg_color_hue(length(unique(dt@meta.data$lineage_knn[1:1019]))-1)),
                              "grey90"),
               group.by = "lineage_knn", 
               do.return = TRUE, pt.size = 4)
p4 <- p4 + guides(color = guide_legend(direction = "horizontal", ncol = 3)) +
  theme(legend.position = c(0.15, 0.9),
        legend.background = element_rect(colour = "black"))

pdf(file = paste(Sys.Date(),p, k, "CCA_tSNE_4.pdf", sep = "_"),
    width = 20, height = 20, useDingbats = F);
plot_grid(p1, p2, p3, p4, nrow = 2)
dev.off()

#### Do DE test for assigned clusters ####
df0.neu <- SetAllIdent(df0.neu, id = "lineage_assignment")

df0.neu.lin2 <- SubsetData(df0.neu, ident.remove = "Uncertain")
cell.type.order <- c(grep("Pvalb", unique(df0.neu.lin2@meta.data$lineage_assignment), value = T),
                     grep("Sst", unique(df0.neu.lin2@meta.data$lineage_assignment), value = T),
                     grep("Vip", unique(df0.neu.lin2@meta.data$lineage_assignment), value = T),
                     grep("Ndnf", unique(df0.neu.lin2@meta.data$lineage_assignment), value = T))
df0.neu.lin2@meta.data$lineage_assignment <- 
  factor(df0.neu.lin2@meta.data$lineage_assignment, levels = cell.type.order)

markers.KNN <- FindAllMarkers(df0.neu.lin2, thresh.use = 0.01, min.pct = 0.25, 
                               only.pos = T,do.print = T)

top40.KNN <- markers.KNN %>% group_by(cluster) %>% top_n(40, avg_diff)

pdf(file = paste(Sys.Date(),"P56_assigned_cells_DE_heatmap.pdf", sep = "_"),
    width = 40, height = 60, useDingbats = F);
DoHeatmap(object = df0.neu.lin2,
          genes.use = c(top40$gene),#rownames(dff.markers)[1:100],
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

top40 <- top40[order(match(top40$cluster,cell.type.order)),]
ord <- order(df0.neu.lin2@meta.data$lineage_assignment)
temp <- df0.neu.lin2@scale.data[top40$gene, ord]

ColSideColors <- c(#gg_color_hue(3)[as.numeric(as.factor(df0.E14.5.prog@meta.data$region))],
  brewer.pal(10, name = "Paired")[as.numeric(sort(df0.neu.lin2@meta.data$lineage_assignment))])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp2))
ColSideColors <- as.matrix(ColSideColors)

colsep <- cumsum(table(df0.neu.lin2@meta.data$lineage_assignment))
rowsep <- cumsum(table(factor(top40$cluster, levels = cell.type.order)))

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

pdf(file = paste(Sys.Date(), "knn_lineage_assigned_DE_heatmap_top40.pdf", sep = "_"),
    width = 10, height = 40, useDingbats = F);
heatmap.3(temp,
          breaks = pairs.breaks,
          symbreaks = T,
          keysize = 0.8,
          key = T,
          main=NA,
          col = col,
          symkey = F,
          cexRow=0.5, 
          cexCol = 0.6, 
          Rowv = F, 
          Colv = F, #as.dendrogram(clu), 
          # ColSideColors = ColSideColors,
          dendrogram = "both",
          colsep = colsep,
          rowsep = rowsep,
          sepwidth = c(0.1, 0.1),
          sepcolor = "white",
          scale = "row",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

#### 5.5 plot average gene expression ####
df0.neu.lin2.ave <- AverageExpression(df0.neu.lin2)

top40 <- top40[order(match(top40$cluster, cell.type.order)),]
temp <- df0.neu.lin2.ave[top40$gene, cell.type.order]

temp2 <- matrix(nrow =  nrow(temp))
colsep <- c()
for(i in 1:length(colnames(temp))){
  n <- length(which(df0.neu.lin2@meta.data$lineage_assignment == colnames(temp)[i]))
  res <- rep(temp[,i], times = n)
  colsep <- c(colsep, n)
  res <- matrix(res, nrow = nrow(temp))
  temp2 <-cbind(temp2, res)
}
colsep <- cumsum(colsep)
temp2 <- temp2[,-1]

ColSideColors <- c(#gg_color_hue(3)[as.numeric(as.factor(df0.E14.5.prog@meta.data$region))],
  brewer.pal(10, name = "Paired")[as.numeric(sort(df0.neu.lin2@meta.data$lineage_assignment))])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp2))
ColSideColors <- as.matrix(ColSideColors)

rowsep <- cunsum(table(top40$cluster))

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

pdf(file = paste(Sys.Date(), "knn_assigned_lineage_DE_ave_heatmap.pdf", sep = "_"),
    width = 10, height = 10, useDingbats = F);
heatmap.3(temp2,
          breaks = pairs.breaks,
          symbreaks = T,
          keysize = 0.8,
          key = T,
          main=NA,
          col = col,
          symkey = F,
          cexRow=0.5, 
          cexCol = 0.6, 
          Rowv = F, 
          Colv = F, #as.dendrogram(clu), 
          # ColSideColors = ColSideColors,
          dendrogram = "both",
          colsep = colsep,
          rowsep = rowsep,
          sepwidth = c(0.1,0.2),
          sepcolor = "white",
          scale = "row",
          labRow = "",#substr(rownames(temp), 20, 100),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

#### MetaNeighbor ####
data <- t(dt@dr$cca.aligned@cell.embeddings)
rec <- c(dt@meta.data$lineage_knn[1:1019], dt@meta.data$lineage[1020:1780])
data.ave <- sapply(unique(rec[dt@meta.data$lineage_knn != "Uncertain"]),
                   function(xx){rowMeans(data[, rec == xx])})

var.genes <- rownames(data.ave)[1:7]
cell.lab <- colnames(data.ave)
exp.lab <- substr(colnames(data.ave), 1, 3)
pheno <- data.frame(Sample_ID = colnames(data.ave), Study_ID = exp.lab, Celltype = cell.lab)

celltype.NV <- run_MetaNeighbor_US(var.genes, data.ave, cell.lab, pheno)
write.csv(celltype.NV, paste(Sys.Date(), "celltype_AUORC.csv", sep = "_"))
# top_hits <- get_top_hits(celltype.NV, pheno, threshold = 0.5, filename = "2017-11-25_neu_P56_MetaNeighbor.txt")

RowSideColors <- rainbow(2, v = 0.8)[as.numeric(as.factor(substr(colnames(celltype.NV), 1, 3)))]

breaks=seq(0,1,length=1001)

pdf(file = paste(Sys.Date(), "neu_MetaNeighbor3.pdf", sep = "_"),
    width = 10, height = 10, useDingbats = F);
heatmap.3(celltype.NV,
          trace="none",
          density.info="none",
          col=col,
          hclustfun = function(x){hclust(d = dist(x, method = "e"),
                                         method = "average")},
          RowSideColors = matrix(RowSideColors, nrow = 1),
          breaks=breaks,
          cexRow=1,cexCol=1,
          margins = c(10, 10))
legend("bottomleft",
       legend = c("Embryonic", "P56"),
       fill = rainbow(2, v = 0.8))
dev.off()

#### MetaNeighbor between P56 cells and knn assigned embyonic cells####
rec <- c(dt@meta.data$lineage_knn[1:1019], dt@meta.data$lineage[1020:1780])
data.ave <- sapply(unique(rec[rec != "Uncertain"]),
                   function(xx){rowMeans(data[, rec == xx])})

var.genes <- rownames(data.ave)[1:7]
cell.lab <- colnames(data.ave)
exp.lab <- substr(colnames(data.ave), 1, 3)
pheno <- data.frame(Sample_ID = colnames(data.ave), Study_ID = exp.lab, Celltype = cell.lab)

celltype.NV <- run_MetaNeighbor_US(var.genes, data.ave, cell.lab, pheno)
write.csv(celltype.NV, paste(Sys.Date(), "P56_Emb_neu_AUORC_knn.csv", sep = "_"))
# top_hits <- get_top_hits(celltype.NV, pheno, threshold = 0.5, filename = "2017-11-25_neu_P56_MetaNeighbor.txt")

RowSideColors <- rainbow(2, v = 0.8)[as.numeric(as.factor(substr(colnames(celltype.NV), 1, 3)))]

clu <- hclust(dist(t(data.ave)),method = "ward.D2")
plot(clu)

breaks=seq(0,1,length=1001)

pdf(file = paste(Sys.Date(), "neu_MetaNeighbor3.pdf", sep = "_"),
    width = 10, height = 10, useDingbats = F);
heatmap.3(celltype.NV,
          trace="none",
          density.info="none",
          col=col,
          Colv = as.dendrogram(clu),
          Rowv = as.dendrogram(clu),
          # hclustfun = function(x){hclust(d = dist(x, method = "e"),method = "complete")},
          RowSideColors = matrix(RowSideColors, nrow = 1),
          breaks=breaks,
          cexRow=1,cexCol=1,
          margins = c(10, 10))
legend("bottomleft",
       legend = c("Embryonic", "P56"),
       fill = rainbow(2, v = 0.8))
dev.off()

#### MetaNeighbor between knn assigned and RF predicted cells ####
df0.neu@meta.data$lineage_assignment <- dt@meta.data$lineage_knn[1:1019]
dt.assigned <- SetAllIdent(df0.neu, id = "lineage_assignment")
dt.assigned <- SubsetData(dt.assigned, ident.remove = c("Uncertain"))
dt.assigned <- AverageExpression(dt.assigned)

data.ave <- cbind(df0.neu.lin.NotChodl.ave, dt.assigned)

var.genes <- c(top40.NotChodl$gene, top40.KNN$gene) #df0.neu@var.genes
cell.lab <- colnames(data.ave)
exp.lab <- c(rep("RF", 6), rep("KNN", 11))
pheno <- data.frame(Sample_ID = colnames(data.ave), Study_ID = exp.lab, Celltype = cell.lab)

celltype.NV <- run_MetaNeighbor_US(var.genes, data.ave, cell.lab, pheno)
write.csv(celltype.NV, paste(Sys.Date(), "knn_RF_AUORC.csv", sep = "_"))
top_hits <- get_top_hits(celltype.NV, pheno, threshold = 0.5, filename = "2017-11-25_neu_P56_MetaNeighbor.txt")

RowSideColors <- rainbow(2, v = 0.8)[as.numeric(as.factor(exp.lab))]

breaks=seq(0,1,length=1001)

clu <- hclust(dist(t(data.ave[var.genes,]),method = "e"),method = "ward.D")
plot(clu)

pdf(file = paste(Sys.Date(), "neu_MetaNeighbor3.pdf", sep = "_"),
    width = 10, height = 10, useDingbats = F);
heatmap.3(celltype.NV,
          trace="none",
          density.info="none",
          col=col,
          # Colv = as.dendrogram(clu),
          # Rowv = as.dendrogram(clu),
          # hclustfun = function(x){hclust(d = dist(x, method = "e"),method = "complete")},
          RowSideColors = matrix(RowSideColors, nrow = 1),
          breaks=breaks,
          cexRow=1,cexCol=1,
          margins = c(10, 10))
legend("bottomleft",
       legend = c("KNN", "RF"),
       fill = rainbow(2, v = 0.8))
dev.off()

#### 7. P56 DE genes for KNN assignment ####
df0.allen@meta.data$lineage_knn <- dt@meta.data$lineage_knn[1020:1780]
df0.allen <- SetAllIdent(df0.allen, id = "primary_type")
df0.allen.knn <- SubsetData(df0.allen, 
                            ident.use = c(substr(unique(df0.allen@meta.data$lineage_knn)[-1], 5, 100)))
df0.allen.knn@ident

P56.knn.markers <- FindAllMarkers(df0.allen.knn,
                                     genes.use =  df0.allen.knn@var.genes,
                                     min.pct = 0.25, 
                                     return.thresh = 0.01,
                                     only.pos = T)

top20 <- P56.knn.markers %>% group_by(cluster) %>% top_n(20, avg_diff)

pdf(file = paste(Sys.Date(),"P56_knn_all_prog_heatmap_top20.pdf", sep = "_"),
    width = 20, height = 20, useDingbats = F);
DoHeatmap(object = df0.allen.knn,
          genes.use = top20$gene,#rownames(df.markers)[1:100],
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

P56.emb.knn.inter.DE <- intersect(toupper(P56.knn.markers$gene), substr(markers.KNN$gene, 20, 100))
markers.KNN.P56.intersected <- markers.KNN[substr(markers.KNN$gene, 20, 100) %in% P56.emb.knn.inter.DE,]

#### 7.2 Intersected P56 DE genes for RF assignment ####
df0.allen <- SetAllIdent(df0.allen, id = "primary_type")
df0.allen.RF <- SubsetData(df0.allen, 
                            ident.use = c(na.omit(as.character(unique(df0.neu.lin.NotChodl@meta.data$RF_cate_whole[1:1019])))))
P56.RF.markers <- FindAllMarkers(df0.allen.RF,
                                  genes.use =  df0.allen.RF@var.genes,
                                  min.pct = 0.25, 
                                  return.thresh = 0.01,
                                  only.pos = T)

top20 <- P56.knn.markers %>% group_by(cluster) %>% top_n(20, avg_diff)

pdf(file = paste(Sys.Date(),"P56_RF_all_prog_heatmap_top20.pdf", sep = "_"),
    width = 20, height = 20, useDingbats = F);
DoHeatmap(object = df0.allen.knn,
          genes.use = top20$gene,#rownames(df.markers)[1:100],
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

P56.emb.RF.inter.DE <- intersect(toupper(P56.RF.markers$gene), substr(markers.RF.NotChodl$gene, 20, 100))
markers.RF.P56.intersected <- markers.RF.NotChodl[substr(markers.RF.NotChodl$gene, 20, 100) %in% P56.emb.RF.inter.DE,]

write.csv(markers.RF.P56.intersected, file = "2017-11-30_RF_intersected_markers.csv")
