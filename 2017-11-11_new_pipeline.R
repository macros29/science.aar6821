## Setup 
options(stringsAsFactors = F)
library(Matrix)
library(dplyr)
library(Seurat)
library(WGCNA)
library(AnnotationDbi)
library(org.Mm.eg.db)
source("~/Scripts/R/heatmap.3.R")
source("~/Scripts/R/gg_color_hue.R")

setwd("~/project/Da_Mi/")
load("2017-11-13_working.RData")

keep <- !(GB$IFC %in% c("180","181","189","190","191","192","193","194","195", "196"))
mito.genes <- grep(pattern = "\\bMT-", x = rownames(expr))
RP.genes <- grep(pattern = "\\bRP", x = rownames(expr))

df0 <- CreateSeuratObject(raw.data = as.vector(expr[-c(mito.genes, RP.genes),keep]), 
                          min.cells = 5, min.genes = 1000, project = "Da_Mi")
df0 <- NormalizeData(object = df0, normalization.method = "LogNormalize", 
                     scale.factor = 1000000)
df0@meta.data <- cbind(df0@meta.data, GB[df0@cell.names,])

#### 1. Regress out cell cycle state ####
cc_genes <- read.csv("cell_cycle_genes.csv", header = F, stringsAsFactors = F)
res <- sapply(cc_genes[,1], function(xx){paste("\\b",xx,"$",sep = "")})
res <- sapply(res, function(xx){grep(xx, rownames(df0@data), value = T)})
cc_genes[,3] <- as.character(res)
g1s.genes <- cc_genes[,3][cc_genes$V2 == "G1/S"]
g2m.genes <- cc_genes[,3][cc_genes$V2 == "G2/M"]
df0 <- CellCycleScoring(df0, s.genes = g1s.genes, g2m.genes = g2m.genes, set.ident = T)
df0@meta.data$CC.Difference <- df0@meta.data$S.Score - df0@meta.data$G2M.Score
df0@meta.data$cc.state <- as.character(df0@ident)
df0@meta.data$cc.state[df0@meta.data$cc.state == "G1"] <- "postmitotic"
df0@meta.data$cc.state[df0@meta.data$cc.state == "S"] <- "G1S"
df0@meta.data$region <- factor(df0@meta.data$region, levels = c("dMGE", "vMGE", "CGE"))

df0 <- ScaleData(df0, vars.to.regress = c("nGene", "CC.Difference"))

df0 <- RunPCA(df0, pc.genes = cc_genes[,3], pcs.compute = 20, pcs.print = 1:12,
              gene.print = 10)

pdf(file = paste(Sys.Date(), "PCA_with_norm_all_cells_cc.pdf", sep = "_"), 
    width = 11, height = 10, useDingbats = F);
PCAPlot(df0)
dev.off()

df0 <- FindVariableGenes(object = df0, mean.function = ExpMean, 
                         dispersion.function = LogVMR, x.low.cutoff = 1, 
                         x.high.cutoff = Inf, y.cutoff = 0.5, y.high.cutoff = 5, 
                         do.plot = F)

df0 <- RunPCA(object = df0, pc.genes=df0@var.genes, pcs.compute = 20, pcs.print = 1:12,
              genes.print = 10)

pdf(paste(Sys.Date(), "All_cells_PCA_region.pdf", sep = "_"), 
    width = 8, height = 8, useDingbats = F);
p1 <- PCAPlot(df0, group.by = "region",do.return = T)
p1 <- p1 + theme(legend.position = c(0.9, 0.9),
                 legend.background = element_rect(size=0.5,
                                   linetype = "solid",
                                   color = "black"))
p1
dev.off()

p2 <- PCAPlot(df0, group.by = "time_point", do.return = T)
p3 <- PCAPlot(df0, group.by = "IFC",do.return = T)
pdf(paste(Sys.Date(), "All_cells_PCA.pdf", sep = "_"), width = 26, height = 8, useDingbats = F);
plot_grid(p1, p2, p3, ncol = 3)

df0 <-JackStraw(df0, num.replicate = 100, do.print = T, num.pc = 20)
JackStrawPlot(df0, PCs = 1:20)

#### 2. Use Rubenstein's list of genes to assign neuron or progenitor identity ####
markers <- read.csv("./Rubenstein PGs vs Neuruon markers 08-11-2017.csv", header = F, stringsAsFactors = F)
# markers <- read.csv("./MZ vs SVZ gene list 08-11-2017.csv", header = F, stringsAsFactors = F)
# markers[,1] <- toupper(unlist(markers[, 1 ]))
markers[,1] <- unique(sapply(c("TOP2A", "GAD1", "GAD2", markers[,1]), 
                             function(x){paste("\\b", x, "$", sep = "")}))
markers[,1] <- as.character(sapply(markers[,1],function(x){grep(x, rownames(df0@data), value = T)}))
markers <- markers[markers[,1] != "character(0)",]
markers <- markers[markers[,2] %in% c("PROGENITOR", "Neuron"),]
dis <- as.dist(1 - cor(df0@scale.data[markers[,1],]))
clu <- hclust(dis, method = "average")
cl <- cutree(clu, k=2)

temp <- df0@scale.data[markers[,1],]

ColSideColors <- c(gg_color_hue(3)[as.numeric(as.factor(df0@meta.data$region))],
                   rainbow(2)[as.numeric(as.factor(df0@meta.data$time_point))],
                   brewer.pal(3, name = "Paired")[cl])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))

col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

pdf(file = paste(Sys.Date(), "all_cells_neuron_progenitor_heatmap.pdf", sep = "_"),
    width = 10, height = 14, useDingbats = F);
heatmap.3(temp,
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8,
          key = F,
          main=NA,
          col = col,
          symkey = F,
          cexRow=0.5, 
          cexCol = 0.6, 
          Rowv = T, 
          Colv = as.dendrogram(clu), 
          ColSideColors = ColSideColors,
          ColSideColorsSize = 0.5,
          dendrogram = "both",
          #scale = "row",
          #colsep = colsep,
          sepcolor = "black",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

#### 3. Use random forest to reassign identity ####
library(randomForest)
cate_bf_fs <- as.factor(cl)
feature_bf_fs <- as.matrix(t(df0@scale.data[df0@var.genes,]))

rf_bf_fs <- randomForest(feature_bf_fs, cate_bf_fs, importance = TRUE, proximity = TRUE)
imp_bf_fs <- importance(rf_bf_fs, type = 1)

fs <- rfcv(feature_bf_fs, cate_bf_fs, cv.fold = 10, scale = "log", step = 0.9, recursive = F)
len <- length(fs$n.var[fs$error.cv == min(fs$error.cv)])
min_fs <- fs$n.var[fs$error.cv == min(fs$error.cv)][len] #get least features
ind <- order(-imp_bf_fs)[1: min_fs]
feature_fs <- feature_bf_fs[ , ind]
cate_fs <- cate_bf_fs

rf_fs <- randomForest(feature_fs, as.factor(cate_fs), importance=TRUE, proximity=TRUE)
fea1_fs <- data.frame()
fea1_fs <- feature_fs[(rf_fs$predicted == '1') & (rf_fs$votes[ , 1] > 0.6), , drop = FALSE]
cat1_fs <- rf_fs$predicted[(rf_fs$predicted =='1') & (rf_fs$votes[ , 1] > 0.6)]
fea2_fs <- data.frame()
fea2_fs <- feature_fs[(rf_fs$predicted == '2') & (rf_fs$votes[ , 2] > 0.6), , drop = FALSE]
cat2_fs <- rf_fs$predicted[(rf_fs$predicted =='2') & (rf_fs$votes[ , 2] > 0.6)]

cate <- as.factor(c(as.character(cat1_fs), as.character(cat2_fs)))
feature <- as.matrix(rbind(fea1_fs, fea2_fs))

set <- sample(1: nrow(feature), nrow(feature), replace = F)
cate <- cate[set]
feature <- feature[set, ] 

rf_whole <- randomForest(feature, as.factor(cate), importance = TRUE, proximity = TRUE)
pred_whole <- predict(rf_whole, newdata = feature_fs)

df0 <- AddMetaData(df0, metadata = pred_whole, col.name = "pred")
ord <- order(df0@meta.data$pred)
temp <- df0@scale.data[c(#E12.5_markers[2:53,1],
                         "ENSMUSG00000055435|MAF",
                         "ENSMUSG00000004366|SST",
                         "ENSMUSG00000074622|MAFB",
                         "ENSMUSG00000028222|CALB1",
                         "ENSMUSG00000000214|TH"), ]

# res <- df0.E12.5.neu@data[,df0.E12.5.neu@data["ENSMUSG00000004366|SST",] > 0]
# temp <- res[c("ENSMUSG00000004366|SST","ENSMUSG00000055435|MAF",
#               "ENSMUSG00000074622|MAFB",
#               colnames(feature)),]
ColSideColors <- c(gg_color_hue(3)[as.numeric(as.factor(df0@meta.data$region))],
                   rainbow(2)[as.numeric(as.factor(df0@meta.data$time_point))],
                   brewer.pal(3, name = "Paired")[df0@meta.data$pred])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- ColSideColors[ord,]

col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

pdf(file = paste(Sys.Date(), "temp_all_cells_neuron_progenitor_heatmap_predicted.pdf", sep = "_"),
    width = 10, height = 5, useDingbats = F);
heatmap.3(temp,
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8,
          key = F,
          main=NA,
          col = col,
          symkey = F,
          cexRow=0.5, 
          cexCol = 0.6, 
          Rowv = T, 
          Colv = F,#as.dendrogram(clu), 
          ColSideColors = ColSideColors,
          ColSideColorsSize = 3,
          # dendrogram = "both",
          scale = "row",
          #colsep = colsep,
          sepcolor = "black",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

#### 4. Analyze progenitor ####
df0 <- SetAllIdent(df0, id = "pred")
df0.prog <- SubsetData(df0, ident.use = 1)

df0.prog <- ScaleData(df0.prog, vars.to.regress = c("nGene", "CC.Difference"))
df0.prog <- FindVariableGenes(object = df0.prog, 
                                       mean.function = ExpMean, 
                                       dispersion.function = LogVMR, 
                                       x.low.cutoff = 1, 
                                       x.high.cutoff = Inf, 
                                       y.cutoff = 0.5, y.high.cutoff = 5, 
                                       do.plot = T)

genes.to.use <- unique(c(E12.5_markers[,1], E14.5_markers[,1]))
df0.prog <- RunPCA(object = df0.prog, 
                   pc.genes = genes.to.use, 
                   pcs.compute = 15, 
                   pcs.print = 1:10, 
                   genes.print = 10, 
                   do.print = T)

df0.prog <-JackStraw(df0.prog, num.replicate = 100, do.print = T, num.pc = 15)
JackStrawPlot(df0.prog, PCs = 1:15)

df0.prog <- RunTSNE(object = df0.prog, dims.use = 1:7, do.fast = TRUE, 
                         perplexity = 30, seed.use = 1)

df0.prog <- FindClusters(object = df0.prog, 
                              reduction.type = "tsne", 
                              dims.use = 1:2, 
                              resolution = 2, 
                              k.scale = 25, 
                              plot.SNN = F, 
                              print.output = F, 
                              save.SNN = F,
                              algorithm = 1,
                              force.recalc = TRUE, 
                              random.seed = 1,
                              k.param = 30) 

df0.prog@meta.data$clusters <- factor(df0.prog@meta.data$res.2, 
                                          levels = c(0,1,8,10,11,3,6,9,12,13,4,2,5,7),
                                      labels = 1:14)

p1 <- TSNEPlot(object = df0.prog, group.by = "clusters", pt.size = 3, do.return = T)
p2 <-TSNEPlot(object = df0.prog, group.by = "region", do.return = T,
              pt.size = 3)
p3 <-TSNEPlot(object = df0.prog, group.by = "time_point", do.return = T,
              pt.size = 3)
p4 <-TSNEPlot(object = df0.prog, group.by = "VZ.SVZ", do.return = T,
              pt.size = 3)

pdf(file = paste(Sys.Date(), "progenitor_tSNE_new.pdf", sep = "_"), width = 22, 
    height = 20, useDingbats = F);
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()
  
  
#### 4.1 Do DE analysis between progenitors of 2 time points ####
df0.prog <- SetAllIdent(df0.prog, id = "time_point")

markers.temp <- FindAllMarkers(df0.prog, min.pct = 0.1, return.thresh = 0.01,only.pos = T)

top100 <- markers.temp %>% group_by(cluster) %>% top_n(100, avg_diff)

pdf(file = paste(Sys.Date(),"E12.5_vs_E14.5_all_prog_heatmap_top100.pdf", sep = "_"),
    width = 20, height = 20, useDingbats = F);
DoHeatmap(object = df0.prog,
          genes.use = top100$gene,#rownames(df.markers)[1:100],
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

#### 4.2.1 Analyze E12.5 progentiors 3 regions ####
df0.prog <- SetAllIdent(df0.prog, id = "time_point")
df0.E12.5.prog <- SubsetData(df0.prog, ident.use = "E12.5")
df0.E12.5.prog@meta.data$region <- factor(df0.E12.5.prog@meta.data$region, 
                                          levels = c("dMGE", "vMGE", "CGE"))
df0.E12.5.prog <- SetAllIdent(df0.E12.5.prog, id = "region")

markers.temp <- FindAllMarkers(df0.E12.5.prog, min.pct = 0.1, return.thresh = 0.01,only.pos = T)

top20 <- markers.temp %>% group_by(cluster) %>% top_n(20, avg_diff)
top20 <- top20[order(match(top20$cluster,c("dMGE", "vMGE", "CGE"))),]

temp <- df0.E12.5.prog@scale.data[top20$gene,]
ord <- order(df0.E12.5.prog@meta.data$region, 
             df0.E12.5.prog@meta.data$pred, 
             df0.E12.5.prog@meta.data$VZ.SVZ,
             runif(ncol(temp), 0, 1))
ColSideColors <- c(gg_color_hue(3)[as.numeric(as.factor(df0.E12.5.prog@meta.data$region))])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- ColSideColors[ord,]

col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

colsep <- cumsum(table(df0.E12.5.prog@meta.data$region))
rowsep <- c(20, 40)

pdf(file = paste(Sys.Date(), "temp_E12.5_progenitor_heatmap_region.pdf", sep = "_"),
    width = 10, height = 8, useDingbats = F);
heatmap.3(temp[,ord],
          breaks = pairs.breaks,
          symbreaks = T,
          keysize = 1,
          key = T,
          main=NA,KeyValueName = NA,
          col = col,
          symkey = T,
          cexRow=0.7, 
          cexCol = 0.6, 
          Rowv = F, 
          Colv = F,#as.dendrogram(clu), 
          ColSideColors = as.matrix(ColSideColors),
          ColSideColorsSize = 3,
          # dendrogram = "both",
          scale = "row",
          colsep = colsep,
          rowsep = rowsep,
          sepcolor = "white",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          lwid = c(1, 5),
          lhei = c(2, 9),
          margins = c(1,6),
          na.rm = F);
dev.off(dev.cur());

#### 4.2.2 Use Random Forest to reclassify VZ and SVZ cells ####
cate_bf_fs <- as.factor(df0.E12.5.prog@ident)
feature_bf_fs <- as.matrix(t(df0.E12.5.prog@scale.data[df0.E12.5.prog@var.genes,]))

rf_bf_fs <- randomForest(feature_bf_fs, cate_bf_fs, importance = TRUE, proximity = TRUE)
imp_bf_fs <- importance(rf_bf_fs, type = 1)

fs <- rfcv(feature_bf_fs, cate_bf_fs, cv.fold = 10, scale = "log", step = 0.9, recursive = F)
len <- length(fs$n.var[fs$error.cv == min(fs$error.cv)])
min_fs <- fs$n.var[fs$error.cv == min(fs$error.cv)][len] #get least features
ind <- order(-imp_bf_fs)[1: min_fs]
feature_fs <- feature_bf_fs[ , ind]
cate_fs <- cate_bf_fs

# markers.temp <- FindAllMarkers(df0.E12.5.prog, min.pct = 0.1, return.thresh = 0.01,only.pos = T)
# top100 <- markers.temp %>% group_by(cluster) %>% top_n(100, avg_diff)

temp <- t(feature_fs)
# dis <- as.dist(1 - cor(temp,use = "pairwise.complete.obs"))
# clu <- hclust(dis, method = "complete")
# cl <- cutree(clu, k=2)

ord <- order(cate_fs)
ColSideColors <- c(#gg_color_hue(3)[as.numeric(as.factor(df0.E12.5.prog@meta.data$region))],
  brewer.pal(3, name = "Paired")[as.numeric(cate_fs)])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- as.matrix(ColSideColors[ord,])
col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

pdf(file = paste(Sys.Date(), "E12.5_progenitor_region_heatmap_RF.pdf", sep = "_"),
    width = 10, height = 14, useDingbats = F);
heatmap.3(temp[,ord],
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8,
          key = F,
          main=NA,
          col = col,
          symkey = F,
          cexRow=0.5, 
          cexCol = 0.6, 
          Rowv = T, 
          Colv = F, #as.dendrogram(clu), 
          ColSideColors = ColSideColors,
          ColSideColorsSize = 0.5,
          dendrogram = "both",
          #scale = "row",
          #colsep = colsep,
          sepcolor = "black",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

#### 4.2.3 Analyze E12.5 progenitors for VZ and SVZ ####
# df0.E12.5.prog <- SetAllIdent(df0.E12.5.prog, id = "region")
# df0.12.vzsvz <- SubsetData(df0.E12.5.prog, ident.remove = "CGE")

RGs <- read.csv("RGs markers.csv", header = F, stringsAsFactors = F)
RGs[,1] <- unique(sapply(RGs[,1], function(x){paste("\\b", x, "$", sep = "")}))
RGs[,1] <- as.character(sapply(RGs[,1],function(x){grep(x, rownames(df0.E12.5.prog@data), value = T)}))
RGs <- RGs[RGs[,1] != "character(0)",]

linnarson <- read.csv("cluster 1 and 2 gene list from Linnarson's paper.csv", header = F, stringsAsFactors = F)
linnarson[,1] <- toupper(linnarson[,1]) 
linnarson[,1] <- unique(sapply(linnarson[,1], function(x){paste("\\b", x, "$", sep = "")}))
linnarson[,1] <- as.character(sapply(linnarson[,1],function(x){grep(x, rownames(df0.E12.5.prog@data), value = T)}))
linnarson <- linnarson[linnarson[,1] != "character(0)",]

vzsvz <- read.csv("2017-11-14_VZ_vs_SVZ_genes.csv", header = F, stringsAsFactors = F)
vzsvz[,1] <- unique(sapply(vzsvz[,1], function(x){paste("\\b", x, "$", sep = "")}))
vzsvz[,1] <- as.character(sapply(vzsvz[,1],function(x){grep(x, rownames(df0.E12.5.prog@data), value = T)}))
vzsvz <- vzsvz[vzsvz[,1] != "character(0)",]
temp <- df0.E12.5.prog@scale.data[vzsvz,]

RGs.cor <- cor(t(temp),use = "pairwise.complete.obs")
linnarson.cor <- cor(t(temp),use = "pairwise.complete.obs")
vzsvz.cor <- cor(t(temp),use = "pairwise.complete.obs")

pdf(file = paste(Sys.Date(), "E12.5_progenitor_vzsvz_correlation.pdf", sep = "_"),
    width = 20, height = 20, useDingbats = F);
heatmap.3(vzsvz.cor,
          breaks = seq(-0.5,0.5, length.out = 1001),
          col = col,
          labRow = substr(rownames(temp), 20, 100),
          labCol = substr(rownames(temp), 20, 100))
dev.off()

dis <- as.dist(1 - cor(temp,use = "pairwise.complete.obs"))
clu <- hclust(dis, method = "ward.D2")
cl <- cutree(clu, k=2)

ColSideColors <- c(gg_color_hue(3)[as.numeric(as.factor(df0.E12.5.prog@meta.data$region))],
                   brewer.pal(3, name = "Paired")[as.numeric(cl)])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))

col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

pdf(file = paste(Sys.Date(), "E12.5_progenitor_MGE_VZ_SVZ_hierarchical_vzsvz.pdf", sep = "_"),
    width = 10, height = 14, useDingbats = F);
heatmap.3(temp,
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8,
          key = F,
          main=NA,
          col = col,
          symkey = F,
          cexRow=0.5, 
          cexCol = 0.6, 
          Rowv = T, 
          Colv = as.dendrogram(clu), 
          ColSideColors = ColSideColors,
          ColSideColorsSize = 0.5,
          dendrogram = "column",
          #scale = "row",
          #colsep = colsep,
          sepcolor = "black",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

#### 4.2.4 Use Random Forest to reclassify VZ and SVZ cells ####
cate_bf_fs <- as.factor(cl)
feature_bf_fs <- as.matrix(t(df0.E12.5.prog@scale.data[df0.E12.5.prog@var.genes,]))

rf_bf_fs <- randomForest(feature_bf_fs, cate_bf_fs, importance = TRUE, proximity = TRUE)
imp_bf_fs <- importance(rf_bf_fs, type = 1)

fs <- rfcv(feature_bf_fs, cate_bf_fs, cv.fold = 10, scale = "log", step = 0.9, recursive = F)
len <- length(fs$n.var[fs$error.cv == min(fs$error.cv)])
min_fs <- fs$n.var[fs$error.cv == min(fs$error.cv)][len] #get least features
ind <- order(-imp_bf_fs)[1: min_fs]
feature_fs <- feature_bf_fs[ , ind]
cate_fs <- cate_bf_fs

rf_fs <- randomForest(feature_fs, as.factor(cate_fs), importance=TRUE, proximity=TRUE)
fea1_fs <- data.frame()
fea1_fs <- feature_fs[(rf_fs$predicted == '1') & (rf_fs$votes[ , 1] > 0.6), , drop = FALSE]
cat1_fs <- rf_fs$predicted[(rf_fs$predicted =='1') & (rf_fs$votes[ , 1] > 0.6)]
fea2_fs <- data.frame()
fea2_fs <- feature_fs[(rf_fs$predicted == '2') & (rf_fs$votes[ , 2] > 0.6), , drop = FALSE]
cat2_fs <- rf_fs$predicted[(rf_fs$predicted =='2') & (rf_fs$votes[ , 2] > 0.6)]

cate <- as.factor(c(as.character(cat1_fs), as.character(cat2_fs)))
feature <- as.matrix(rbind(fea1_fs, fea2_fs))

set <- sample(1: nrow(feature), nrow(feature), replace = F)
cate <- cate[set]
feature <- feature[set, ] 

rf_whole <- randomForest(feature, as.factor(cate), importance = TRUE, proximity = TRUE)
pred_whole <- predict(rf_whole, newdata = feature_fs)

dis <- as.dist(1 - cor(df0.E12.5.prog@scale.data[RGs[,1],],use = "pairwise.complete.obs"))
clu <- hclust(dis, method = "average")
cl <- cutree(clu, k=2)

temp <- t(feature_fs)
# dis <- as.dist(1 - cor(temp,use = "pairwise.complete.obs"))
# clu <- hclust(dis, method = "complete")
# cl <- cutree(clu, k=2)

ord <- order()
ColSideColors <- c(gg_color_hue(3)[as.numeric(as.factor(df0.E12.5.prog@meta.data$region))],
  brewer.pal(3, name = "Paired")[as.numeric(pred_whole)])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- as.matrix(ColSideColors[ord,])
col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

pdf(file = paste(Sys.Date(), "E12.5_progenitor_VZ_SVZ_heatmap_post-RF.pdf", sep = "_"),
    width = 10, height = 14, useDingbats = F);
heatmap.3(temp[,ord],
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8,
          key = F,
          main=NA,
          col = col,
          symkey = F,
          cexRow=0.5, 
          cexCol = 0.6, 
          Rowv = T, 
          Colv = F, #as.dendrogram(clu), 
          ColSideColors = ColSideColors,
          ColSideColorsSize = 0.5,
          dendrogram = "both",
          #scale = "row",
          #colsep = colsep,
          sepcolor = "black",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

df0.E12.5.prog <- AddMetaData(df0.E12.5.prog, metadata = pred_whole, col.name = "VZ.SVZ")

#### 4.2.5 Get markers to do tSNE ####
# E12.5_markers <- read.csv("2017-11-15_E12.5_marker_genes2.csv", header = F, stringsAsFactors = F)
# E12.5_markers[,1] <- sapply(E12.5_markers[,1], function(x){paste("\\b", x, "$", sep = "")})
# E12.5_markers[,1] <- as.character(sapply(E12.5_markers[,1],function(x){grep(x, rownames(df0.E12.5.prog@data), value = T)}))
# E12.5_markers <- E12.5_markers[E12.5_markers[,1] != "character(0)",]

E12.5_markers <- read.csv("2017-11-15_E12.5_marker_genes.csv", header = T, stringsAsFactors = F)

df0.E12.5.prog <- RunPCA(df0.E12.5.prog, pc.genes = E12.5_markers[,1],
                         pcs.print = 1:10, pcs.compute = 10, genes.print = 10)

df0.E12.5.prog <-JackStraw(df0.E12.5.prog, num.replicate = 100, do.print = F,num.pc = 10)
JackStrawPlot(df.int, PCs = 1:10)

df0.E12.5.prog <- RunTSNE(df0.E12.5.prog, dims.use = c(1:2),#pcs[[i]], 
                         perplexity = 30,
                         seed.use = 3)

df0.E12.5.prog <- FindClusters(object = df0.E12.5.prog, 
                       reduction.type = "pca", 
                       dims.use = 1:2, 
                       resolution = 0.8, # dMGE, CGE = 0.8, E14.5 vMGE = 2
                       k.scale = 25, # dMGE, CGE = 50
                       prune.SNN = 0, #dMGE, CGE = 0.5
                       plot.SNN = F, 
                       print.output = F, 
                       save.SNN = F,
                       algorithm = 2,
                       force.recalc = TRUE, 
                       random.seed = 1,
                       k.param = 30) 

p1 <- TSNEPlot(object = df0.E12.5.prog, pt.size = 3, do.return = T)
p2 <-TSNEPlot(object = df0.E12.5.prog, group.by = "VZ.SVZ", do.return = T,
              pt.size = 3)
p3 <-TSNEPlot(object = df0.E12.5.prog, group.by = "region", do.return = T,
              pt.size = 3)

pdf(file = paste(Sys.Date(),"E12.5_prog_tSNE.pdf", sep = "_"), width = 33, 
    height = 10, useDingbats = F);
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

#### Clustering on E12.5 progenitors ####
df0.E12.5.prog@meta.data$clusters <- 
  factor(df0.E12.5.prog@meta.data$res.3, levels = c(0,11,5,2,1,7,3,10,6,9,12,8,4),
         labels = c(1:13))
df0.E12.5.prog@meta.data$clusters <- paste("E12.5", as.character(df0.E12.5.prog@meta.data$clusters))

df0.E12.5.prog <- SetAllIdent(df0.E12.5.prog, id = "clusters")

markers <- FindAllMarkers(object = df0.E12.5.prog, only.pos = TRUE, min.pct = 0.25,
                          return.thresh = 0.01)

top20 <- markers %>% group_by(cluster) %>% top_n(20, avg_diff)

pdf(file = paste(Sys.Date(), "E12.5_prog_heatmap_top20.pdf", sep = "_"), width = 20, 
    height = 40, useDingbats = F);
DoHeatmap(object = df0.E12.5.prog,
          genes.use = top20$gene,#rownames(df.markers)[1:100], #genes,
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

#### Spec analysis E12.5 ####
## Load data
library(parallel)
library(matrixStats)
source("~/Scripts/R/spec.R")

Data <- df0.E12.5.prog@scale.data[df0.E12.5.prog@var.genes,]
cell_type = df0.E12.5.prog@meta.data$clusters

n <- length(unique(cell_type))
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, outfile='LOG.TXT')
clusterExport(cl, c("cell_type", "Data"))

spec_scores <- parApply(cl = cl, Data, 1, function(x){
  source("~/Scripts/R/spec.R");
  gene <- as.numeric(x);
  opt_bin <- optimizebinsize(gene, header = cell_type);
  cat("=")
  return(unlist(specy(gene, header = cell_type, binsize = opt_bin)));
})
stopCluster(cl)

colnames(spec_scores) <- df0.E12.5.prog@var.genes
rownames(spec_scores)[1:n] <- unique(sort(cell_type));

res <- spec_scores[1:n,]
write.csv(res, file = paste(Sys.Date(), "E12.5_prog_spec_scores.csv", sep = "_"))

#save(spec_scores, file = "spec_scores.Rdata")

spec.genes <- apply(res,1, function(xx){
  names(xx[xx > 0.07])
  #names(sort(xx, decreasing = T))[1:10]
})

as.vector(spec.genes)

cat("\n", "DONE")

temp <- df0.E12.5.prog@scale.data[ unique(unlist(spec.genes)), ] 
ord <- order(df0.E12.5.prog@meta.data$clusters,
             df0.E12.5.prog@meta.data$VZ.SVZ,
             df0.E12.5.prog@meta.data$region, runif(ncol(temp), 0,1))

ColSideColors <- c(
  gg_color_hue(13)[as.numeric(df0.E12.5.prog@meta.data$clusters)],
  gg_color_hue(3)[as.numeric(df0.E12.5.prog@meta.data$region)],
  brewer.pal(3, name = "Paired")[df0.E12.5.prog@meta.data$VZ.SVZ])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- as.matrix(ColSideColors[ord,])
# colsep <- cumsum(table(df0.E12.5.prog@meta.data$clusters))
# rowsep <- cumsum(table(top20$cluster))

pairs.breaks <- seq(-1.5, 1.5, length.out=1001)

pdf(file = paste(Sys.Date(), "E12.5_prog_spec_genes.pdf", sep = "_"),
    width = 10, height = 8, useDingbats = F);
heatmap.3(temp[,ord],
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8,
          key = T,
          main=NA,
          col = col,
          symkey = F,
          cexRow=0.6, 
          cexCol = 1, 
          Rowv = F, 
          Colv = F, #as.dendrogram(clu), 
          ColSideColors = ColSideColors,
          dendrogram = "none",
          # colsep = colsep,
          # rowsep = rowsep,
          sepwidth = c(0.1, 0.1),
          sepcolor = "white",
          scale = "row",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          margin = c(1, 8),
          na.rm = F);
dev.off(dev.cur());

genes.to.plot <- c("ENSMUSG00000027430|DTD1","ENSMUSG00000027667|ZFP639",
                   "ENSMUSG00000028766|ALPL",
                   "ENSMUSG00000041702|BTBD7", "ENSMUSG00000032375|APH1B",
                   #"ENSMUSG00000026890|LHX6", "ENSMUSG00000028212|CCNE2",
                   "ENSMUSG00000026249|SERPINE2",
                   "ENSMUSG00000022437|SAMM50","ENSMUSG00000020121|SRGAP1",
                   "ENSMUSG00000052539|MAGI3", "ENSMUSG00000029819|NPY", 
                   "ENSMUSG00000042607|ASB4", "ENSMUSG00000024227|2610034M16RIK",
                   "ENSMUSG00000090401|SNORD123","ENSMUSG00000030605|MFGE8", 
                   "ENSMUSG00000005360|SLC1A3","ENSMUSG00000048402|GLI2",
                   "ENSMUSG00000024981|ACSL5","ENSMUSG00000029179|ZCCHC4",
                   "ENSMUSG00000038023|ATP6V0A2","ENSMUSG00000029189|SEL1L3")

pdf(file = paste(Sys.Date(), "E12.5_prog_spec_violin.pdf", sep = "_"),
    width = 5, height = 30, useDingbats = F);
VlnPlot(df0.E12.5.prog, features.plot = genes.to.plot,nCol = 1,
        point.size.use = 0, group.by = "clusters")
dev.off()

df0.E12.5.prog <- SetAllIdent(df0.E12.5.prog, id = "clusters")
pdf(file = paste(Sys.Date(), "E12.5_prog_vln_plot.pdf", sep = "_"),
    width = 50, height = 15, useDingbats = F);
VlnPlot(df0.E12.5.prog, features.plot = c("ENSMUSG00000001496|NKX2-1",
                                          "ENSMUSG00000004366|SST",
                                              "ENSMUSG00000030551|NR2F2",
                                              "ENSMUSG00000069171|NR2F1",
                                              "ENSMUSG00000041309|NKX6-2",
                                              "ENSMUSG00000048402|GLI2",
                                              "ENSMUSG00000004151|ETV1",
                                              "ENSMUSG00000034227|FOXJ1",
                                              "ENSMUSG00000028201|LHX8",
                                              "ENSMUSG00000027210|MEIS2",
                                              "ENSMUSG00000048562|SP8",
                                              "ENSMUSG00000010175|PROX1",
                                              "ENSMUSG00000055435|MAF"),#, 
            # cols.use = c("grey", "blue"),pt.size = 3
        )
dev.off()

pdf(file = paste(Sys.Date(), "E12.5_prog_vz_violin_plot_VZSVZ.pdf", sep = "_"),
    width = 20, height = 15, useDingbats = F);
VlnPlot(df0.E12.5.prog.vz, features.plot = c(#E14.5_markers[1:16,1])
  RG_markers[,1],
                                             "ENSMUSG00000001496|NKX2-1",
                                             "ENSMUSG00000030551|NR2F2",
                                             "ENSMUSG00000069171|NR2F1",
                                             "ENSMUSG00000041309|NKX6-2",
                                             "ENSMUSG00000048402|GLI2",
                                             "ENSMUSG00000004151|ETV1",
                                             "ENSMUSG00000034227|FOXJ1",
                                             "ENSMUSG00000028201|LHX8",
                                             "ENSMUSG00000027210|MEIS2",
                                             "ENSMUSG00000048562|SP8",
                                             "ENSMUSG00000010175|PROX1",
  "ENSMUSG00000074622|MAFB",
  "ENSMUSG00000055435|MAF",
  "ENSMUSG00000028222|CALB1",
  "ENSMUSG00000000214|TH")
        # cols.use = c("grey", "blue"),
        # pt.size = 3
)
dev.off()

df0.E12.5.prog <- SetAllIdent(df0.E12.5.prog, id = "clusters")
df0.E12.5.prog@meta.data$clusters2 <- df0.E12.5.prog@meta.data$clusters
df0.E12.5.prog@meta.data$clusters2[df0.E12.5.prog@meta.data$clusters2 %in% c("E12.5 13","E12.5 7")] <- "MAF"
df0.E12.5.prog@meta.data$clusters2[df0.E12.5.prog@meta.data$clusters2 %in% c("E12.5 4", "E12.5 5")] <- "Non-MAF"
df0.E12.5.prog <- SetAllIdent(df0.E12.5.prog, id = "clusters2")
markers.temp <- FindMarkers(df0.E12.5.prog,test.use = "t",
                            genes.use = c(#E14.5_markers[1:16,1])
                            "ENSMUSG00000001496|NKX2-1",
                            "ENSMUSG00000030551|NR2F2",
                            "ENSMUSG00000069171|NR2F1",
                            "ENSMUSG00000041309|NKX6-2",
                            "ENSMUSG00000048402|GLI2",
                            "ENSMUSG00000004151|ETV1",
                            "ENSMUSG00000034227|FOXJ1",
                            "ENSMUSG00000028201|LHX8",
                            "ENSMUSG00000027210|MEIS2",
                            "ENSMUSG00000048562|SP8",
                            "ENSMUSG00000010175|PROX1",
                            "ENSMUSG00000074622|MAFB",
                            "ENSMUSG00000055435|MAF",
                            "ENSMUSG00000028222|CALB1",
                            "ENSMUSG00000000214|TH"),
                            ident.1 = "E12.5 5", ident.2 = "E12.5 7",min.pct = 0)
markers.temp["ENSMUSG00000055435|MAF",]

#### 4.2.5.2 Analyze individual clusters ####
df0.E12.5.prog <- SetAllIdent(df0.E12.5.prog, id = "res.3")

markers.prog.5vs7 <- FindMarkers(df0.E12.5.prog,ident.1 = "5",ident.2 = "1",only.pos = F,
                                 min.pct = 0.2)
markers.prog.5vs7 <- markers.prog.5vs7[markers.prog.5vs7$p_val < 0.01,]

markers.prog.5vs7$gene <- rownames(markers.prog.5vs7)

markers.prog.5vs7 <- markers.prog.5vs7[order(markers.prog.5vs7$avg_diff,decreasing = T),]

pdf(file = paste(Sys.Date(),"P5_vz_P7_heatmap_top100.pdf", sep = "_"),
    width = 20, height = 40, useDingbats = F);
DoHeatmap(object = df0.E12.5.prog,
          cells.use = df0.E12.5.prog@cell.names[df0.E12.5.prog@ident %in% c("5", "7")],
          genes.use = rownames(markers.prog.5vs7),#rownames(dff.markers)[1:100],
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

grep("ETV", rownames(df0.E12.5.prog@raw.data), value = T)

pdf(file = paste(Sys.Date(),"P5_vz_P1_violin.pdf", sep = "_"),
    width = 5, height = 10, useDingbats = F);
VlnPlot(object = df0.E12.5.prog,ident.include = c("5","1"),nCol = 1,
        features.plot = c("ENSMUSG00000055435|MAF","ENSMUSG00000074622|MAFB",
                          "ENSMUSG00000004366|SST", 
                          "ENSMUSG00000028222|CALB1","ENSMUSG00000062209|ERBB4",
                          "ENSMUSG00000054252|FGFR3","ENSMUSG00000031278|ACSL4",
                          "ENSMUSG00000069171|NR2F1", "ENSMUSG00000004151|ETV1",
                          "ENSMUSG00000028201|LHX8","ENSMUSG00000055639|DACH1",
                          "ENSMUSG00000056121|FEZ2"),
        size.title.use = 10)
dev.off()

genes.to.plot <- c("ENSMUSG00000055435|MAF","ENSMUSG00000074622|MAFB",
  "ENSMUSG00000004366|SST", 
  "ENSMUSG00000028222|CALB1","ENSMUSG00000062209|ERBB4",
  "ENSMUSG00000054252|FGFR3","ENSMUSG00000031278|ACSL4",
  "ENSMUSG00000069171|NR2F1", "ENSMUSG00000004151|ETV1",
  "ENSMUSG00000028201|LHX8","ENSMUSG00000055639|DACH1",
  "ENSMUSG00000056121|FEZ2")

data.to.plot <- data.frame(t(df0.E12.5.prog@scale.data[genes.to.plot,]))
data.to.plot$clusters <- df0.E12.5.prog@meta.data$res.3

data.to.plot <- melt(data.to.plot, id.vars = "clusters")
data.to.plot2 <- data.to.plot[data.to.plot$clusters %in% c("2", "3"),]

data.to.plot2$variable <- substr(data.to.plot2$variable, 20, 100)
data.to.plot2$value[data.to.plot2$value <= 0] <- 0 
data.to.plot2$value[data.to.plot2$value >= 2] <- 2

pdf(file = paste(Sys.Date(),"P2_vs_P3_violin.pdf", sep = "_"),
    width = 5, height = 10, useDingbats = F);
ggplot(data.to.plot2, aes(x = clusters, y = value, fill = clusters)) +
  geom_jitter() +
  geom_violin() +
  facet_wrap(~variable)
dev.off()

markers.prog.2vs3[5516,]
#### 4.2.6 Analyze E12.5 VZ cells only ####
df0.E12.5.prog <- SetAllIdent(df0.E12.5.prog, id = "VZ.SVZ")
df0.E12.5.prog.vz <- SubsetData(df0.E12.5.prog, ident.use = 2)

df0.E12.5.prog.vz <- FindVariableGenes(object = df0.E12.5.prog.vz, 
                                       mean.function = ExpMean, 
                                       dispersion.function = LogVMR, 
                                       x.low.cutoff = 1, 
                                       x.high.cutoff = Inf, 
                                       y.cutoff = 0.5, y.high.cutoff = 5, 
                                       do.plot = T)

df0.E12.5.prog.vz <- RunPCA(object = df0.E12.5.prog.vz, 
                        pc.genes=E12.5_markers[,1], 
                        pcs.compute = 15, 
                        pcs.print = 1:10,
                        genes.print = 10)

df0.E12.5.prog.vz <-JackStraw(df0.E12.5.prog.vz, 
                           num.replicate = 100, do.print = F,num.pc = 15)
JackStrawPlot(df0.E14.5.prog, PCs = 1:15)

p1 <- PCAPlot(df0.E12.5.prog.vz, group.by = "region",do.return = T)
p2 <- PCAPlot(df0.E12.5.prog.vz, group.by = "time_point", do.return = T)
p3 <- PCAPlot(df0.E12.5.prog.vz, group.by = "IFC",do.return = T)
pdf(paste(Sys.Date(), "E12.5_prog_PCA.pdf", sep = "_"), width = 26, height = 8, useDingbats = F);
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

df0.E12.5.prog.vz <- RunTSNE(df0.E12.5.prog.vz, dims.use = c(1:2),#pcs[[i]], 
                         perplexity = 20,
                         seed.use = 3)

df0.E12.5.prog.vz <- FindClusters(object = df0.E12.5.prog.vz, 
                               reduction.type = "pca", 
                               dims.use = 1:2, 
                               resolution = 2, # dMGE, CGE = 0.8, E14.5 vMGE = 2
                               k.scale = 25, # dMGE, CGE = 50
                               prune.SNN = 0, #dMGE, CGE = 0.5
                               plot.SNN = F, 
                               print.output = F, 
                               save.SNN = F,
                               algorithm = 1,
                               force.recalc = TRUE, 
                               random.seed = 1,
                               k.param = 30) 

p1 <- TSNEPlot(object = df0.E12.5.prog.vz, pt.size = 3, do.return = T)
p2 <-TSNEPlot(object = df0.E12.5.prog.vz, group.by = "region", do.return = T,
              pt.size = 3)
# p3 <-TSNEPlot(object = df0.E12.5.prog.vz, group.by = "", do.return = T,
#               pt.size = 3)
pdf(file = paste(Sys.Date(),"E12.5_prog_vz_tSNE2.pdf", sep = "_"), width = 22, 
    height = 10, useDingbats = F);
plot_grid(p1, p2, 
          # p3, 
          ncol = 2)
dev.off()

markers.temp <- FindAllMarkers(df0.E12.5.prog.vz, min.pct = 0.1, return.thresh = 0.01,only.pos = T)

top20 <- markers.temp %>% group_by(cluster) %>% top_n(20, avg_diff)

pdf(file = paste(Sys.Date(),"E12.5_vz_prog_heatmap_top20.pdf", sep = "_"),
    width = 20, height = 40, useDingbats = F);
DoHeatmap(object = df0.E12.5.prog.vz,
          genes.use = top20$gene,#rownames(dff.markers)[1:100],
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

pdf(file = paste(Sys.Date(), "E12.5_prog_vz_violin_plot.pdf", sep = "_"),
    width = 20, height = 15, useDingbats = F);
VlnPlot(df0.E12.5.prog.vz, features.plot = c(RG_markers[,1],
          "ENSMUSG00000001496|NKX2-1",
            "ENSMUSG00000030551|NR2F2",
            "ENSMUSG00000069171|NR2F1",
            "ENSMUSG00000041309|NKX6-2",
            "ENSMUSG00000048402|GLI2",
            "ENSMUSG00000004151|ETV1",
            "ENSMUSG00000034227|FOXJ1",
            "ENSMUSG00000028201|LHX8",
            "ENSMUSG00000027210|MEIS2",
            "ENSMUSG00000048562|SP8",
            "ENSMUSG00000010175|PROX1"),
        # cols.use = c("grey", "blue"),
        # pt.size = 3
)
dev.off()

#### 4.2.7 Analyze E12.5 SVZ cells only ####
df0.E12.5.prog <- SetAllIdent(df0.E12.5.prog, id = "VZ.SVZ")
df0.E12.5.prog.svz <- SubsetData(df0.E12.5.prog, ident.use = 1)

df0.E12.5.prog.svz <- RunPCA(object = df0.E12.5.prog.svz, 
                            pc.genes=E12.5_markers[,1], 
                            pcs.compute = 15, 
                            pcs.print = 1:10,
                            genes.print = 10)

df0.E12.5.prog.svz <-JackStraw(df0.E12.5.prog.svz, 
                              num.replicate = 100, do.print = F,num.pc = 15)
JackStrawPlot(df0.E14.5.prog, PCs = 1:15)

p1 <- PCAPlot(df0.E12.5.prog.svz, group.by = "region",do.return = T)
p2 <- PCAPlot(df0.E12.5.prog.svz, group.by = "time_point", do.return = T)
p3 <- PCAPlot(df0.E12.5.prog.svz, group.by = "IFC",do.return = T)
pdf(paste(Sys.Date(), "E12.5_prog_PCA.pdf", sep = "_"), width = 26, height = 8, useDingbats = F);
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

df0.E12.5.prog.svz <- RunTSNE(df0.E12.5.prog.svz, dims.use = c(1:2),#pcs[[i]], 
                             perplexity = 20,
                             seed.use = 3)

df0.E12.5.prog.svz <- FindClusters(object = df0.E12.5.prog.svz, 
                                  reduction.type = "pca", 
                                  dims.use = 1:2, 
                                  resolution = 2, # dMGE, CGE = 0.8, E14.5 vMGE = 2
                                  k.scale = 25, # dMGE, CGE = 50
                                  prune.SNN = 0, #dMGE, CGE = 0.5
                                  plot.SNN = F, 
                                  print.output = F, 
                                  save.SNN = F,
                                  algorithm = 1,
                                  force.recalc = TRUE, 
                                  random.seed = 1,
                                  k.param = 30) 

p1 <- TSNEPlot(object = df0.E12.5.prog.svz, pt.size = 3, do.return = T)
p2 <-TSNEPlot(object = df0.E12.5.prog.svz, group.by = "region", do.return = T,
              pt.size = 3)
# p3 <-TSNEPlot(object = df0.E12.5.prog.svz, group.by = "", do.return = T,
#               pt.size = 3)
pdf(file = paste(Sys.Date(),"E12.5_prog_svz_tSNE2.pdf", sep = "_"), width = 22, 
    height = 10, useDingbats = F);
plot_grid(p1, p2, 
          # p3, 
          ncol = 2)
dev.off()

markers.temp <- FindAllMarkers(df0.E12.5.prog.svz, min.pct = 0.1, return.thresh = 0.01,only.pos = T)

top20 <- markers.temp %>% group_by(cluster) %>% top_n(20, avg_diff)

pdf(file = paste(Sys.Date(),"E12.5_svz_prog_heatmap_top20.pdf", sep = "_"),
    width = 20, height = 40, useDingbats = F);
DoHeatmap(object = df0.E12.5.prog.svz,
          genes.use = top20$gene,#rownames(dff.markers)[1:100],
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

pdf(file = paste(Sys.Date(), "E12.5_prog_svz_violin_plot_region.pdf", sep = "_"),
    width = 20, height = 15, useDingbats = F);
VlnPlot(df0.E12.5.prog.svz, features.plot = c(RG_markers[,1],
          "ENSMUSG00000001496|NKX2-1",
            "ENSMUSG00000030551|NR2F2",
            "ENSMUSG00000069171|NR2F1",
            "ENSMUSG00000041309|NKX6-2",
            "ENSMUSG00000048402|GLI2",
            "ENSMUSG00000004151|ETV1",
            "ENSMUSG00000034227|FOXJ1",
            "ENSMUSG00000028201|LHX8",
            "ENSMUSG00000027210|MEIS2",
            "ENSMUSG00000048562|SP8",
            "ENSMUSG00000010175|PROX1"),
        # cols.use = c("grey", "blue"),
        # pt.size = 3
)
dev.off()

#### 4.3.1 Analyze E14.5 progentiors 3 regions ####
df0.prog <- SetAllIdent(df0.prog, id = "time_point")
df0.E14.5.prog <- SubsetData(df0.prog, ident.use = "E14.5")
df0.E14.5.prog@meta.data$region <- factor(df0.E14.5.prog@meta.data$region, 
                                          levels = c("dMGE", "vMGE", "CGE"))
df0.E14.5.prog <- SetAllIdent(df0.E14.5.prog, id = "region")

markers.temp <- FindAllMarkers(df0.E14.5.prog, min.pct = 0.1, return.thresh = 0.01, only.pos = T)

top20 <- markers.temp %>% group_by(cluster) %>% top_n(20, avg_diff)
top20 <- top20[order(match(top20$cluster,c("dMGE", "vMGE", "CGE"))),]

temp <- df0.E14.5.prog@scale.data[top20$gene,]
ord <- order(df0.E14.5.prog@meta.data$region, 
             df0.E14.5.prog@meta.data$pred, 
             df0.E14.5.prog@meta.data$VZ.SVZ,
             runif(ncol(temp), 0, 1))
ColSideColors <- c(gg_color_hue(3)[as.numeric(as.factor(df0.E14.5.prog@meta.data$region))])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- ColSideColors[ord,]

col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

colsep <- cumsum(table(df0.E14.5.prog@meta.data$region))
rowsep <- c(20, 40)

pdf(file = paste(Sys.Date(), "temp_E14.5_progenitor_heatmap_region.pdf", sep = "_"),
    width = 10, height = 8, useDingbats = F);
heatmap.3(temp[,ord],
          breaks = pairs.breaks,
          symbreaks = T,
          keysize = 1,
          key = T,
          main=NA,KeyValueName = NA,
          col = col,
          symkey = T,
          cexRow=0.7, 
          cexCol = 0.6, 
          Rowv = F, 
          Colv = F,#as.dendrogram(clu), 
          ColSideColors = as.matrix(ColSideColors),
          ColSideColorsSize = 3,
          # dendrogram = "both",
          scale = "row",
          colsep = colsep,
          rowsep = rowsep,
          sepcolor = "white",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          lwid = c(1, 5),
          lhei = c(2, 9),
          margins = c(1,6),
          na.rm = F);
dev.off(dev.cur());

#### 4.3.2 Analyze E14.5 progenitors for VZ and SVZ ####
# df0.E14.5.prog <- SetAllIdent(df0.E14.5.prog, id = "region")
# df0.12.vzsvz <- SubsetData(df0.E14.5.prog, ident.remove = "CGE")

RGs <- read.csv("RGs markers.csv", header = F, stringsAsFactors = F)
RGs[,1] <- unique(sapply(RGs[,1], function(x){paste("\\b", x, "$", sep = "")}))
RGs[,1] <- as.character(sapply(RGs[,1],function(x){grep(x, rownames(df0.E14.5.prog@data), value = T)}))
RGs <- RGs[RGs[,1] != "character(0)",]

linnarson <- read.csv("cluster 1 and 2 gene list from Linnarson's paper.csv", header = F, stringsAsFactors = F)
linnarson[,1] <- toupper(linnarson[,1]) 
linnarson[,1] <- unique(sapply(linnarson[,1], function(x){paste("\\b", x, "$", sep = "")}))
linnarson[,1] <- as.character(sapply(linnarson[,1],function(x){grep(x, rownames(df0.E14.5.prog@data), value = T)}))
linnarson <- linnarson[linnarson[,1] != "character(0)",]

vzsvz <- read.csv("2017-11-14_E14.5_VZ_vs_SVZ_genes.csv", header = F, stringsAsFactors = F)
vzsvz[,1] <- unique(sapply(vzsvz[,1], function(x){paste("\\b", x, "$", sep = "")}))
vzsvz[,1] <- as.character(sapply(vzsvz[,1],function(x){grep(x, rownames(df0.E14.5.prog@data), value = T)}))
vzsvz <- vzsvz[vzsvz[,1] != "character(0)",]

temp <- df0.E14.5.prog@scale.data[vzsvz,]

RGs.cor <- cor(t(temp),use = "pairwise.complete.obs")
linnarson.cor <- cor(t(temp),use = "pairwise.complete.obs")
vzsvz.cor <- cor(t(temp),use = "pairwise.complete.obs")

pdf(file = paste(Sys.Date(), "E14.5_progenitor_vzsvz_correlation.pdf", sep = "_"),
    width = 20, height = 20, useDingbats = F);
heatmap.3(vzsvz.cor,
          breaks = seq(-0.5,0.5, length.out = 1001),
          col = col,
          labRow = substr(rownames(temp), 20, 100),
          labCol = substr(rownames(temp), 20, 100))
dev.off()

dis <- as.dist(1 - cor(temp,use = "pairwise.complete.obs"))
clu <- hclust(dis, method = "ward.D2")
cl <- cutree(clu, k=2)

temp <- (df0.E14.5.prog@scale.data[top100$gene,])
ColSideColors <- c(gg_color_hue(12)[as.numeric(as.factor(df0.E14.5.prog@meta.data$IFC))],
                   brewer.pal(3, name = "Paired")[as.numeric(as.factor(df0.E14.5.prog@meta.data$region))])
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ord <- order(as.numeric(as.factor(df0.E14.5.prog@meta.data$region), as.numeric(as.factor(df0.E14.5.prog@meta.data$IFC))))

col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

pdf(file = paste(Sys.Date(), "E14.5_progenitor_MGE_VZ_SVZ_hierarchical_vzsvz_temp.pdf", sep = "_"),
    width = 10, height = 14, useDingbats = F);
heatmap.3(temp[,ord],
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8,
          key = F,
          main=NA,
          col = col,
          symkey = F,
          cexRow=0.5, 
          cexCol = 0.6, 
          Rowv = F, 
          Colv = F,#as.dendrogram(clu), 
          ColSideColors = ColSideColors[ord,],
          ColSideColorsSize = 0.5,
          dendrogram = "column",
          #scale = "row",
          #colsep = colsep,
          sepcolor = "black",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

#### Spec analysis E14.5 ####
## Load data
library(parallel)
library(matrixStats)
source("~/Scripts/R/spec.R")

Data <- df0.E14.5.prog@scale.data[df0.E14.5.prog@var.genes,]
cell_type = df0.E14.5.prog@meta.data$clusters

n <- length(unique(cell_type))
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, outfile='LOG.TXT')
clusterExport(cl, c("cell_type", "Data"))

spec_scores <- parApply(cl = cl, Data, 1, function(x){
  source("~/Scripts/R/spec.R");
  gene <- as.numeric(x);
  opt_bin <- optimizebinsize(gene, header = cell_type);
  cat("=")
  return(unlist(specy(gene, header = cell_type, binsize = opt_bin)));
})
stopCluster(cl)

colnames(spec_scores) <- df0.E14.5.prog@var.genes
rownames(spec_scores)[1:n] <- unique(sort(cell_type));

res <- spec_scores[1:n,]
write.csv(res, file = paste(Sys.Date(), "E14.5_prog_spec_scores.csv", sep = "_"))

spec.genes <- apply(res,1, function(xx){
  #names(xx[xx > 0.07])
  names(sort(xx, decreasing = T))[1:100]
})

temp <- df0.E14.5.prog@scale.data[ unique(spec.genes), ] 
ord <- order(df0.E14.5.prog@meta.data$clusters,
             df0.E14.5.prog@meta.data$VZ.SVZ,
             df0.E14.5.prog@meta.data$region, runif(ncol(temp), 0,1))

ColSideColors <- c(
  gg_color_hue(13)[as.numeric(df0.E14.5.prog@meta.data$clusters)],
  gg_color_hue(3)[as.numeric(df0.E14.5.prog@meta.data$region)],
  brewer.pal(3, name = "Paired")[df0.E14.5.prog@meta.data$VZ.SVZ])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- as.matrix(ColSideColors[ord,])
# colsep <- cumsum(table(df0.E14.5.prog@meta.data$clusters))
# rowsep <- cumsum(table(top20$cluster))

pairs.breaks <- seq(-1.5, 1.5, length.out=1001)

pdf(file = paste(Sys.Date(), "E14.5_prog_spec_genes.pdf", sep = "_"),
    width = 10, height = 20, useDingbats = F);
heatmap.3(temp[,ord],
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8,
          key = T,
          main=NA,
          col = col,
          symkey = F,
          cexRow=0.6, 
          cexCol = 1, 
          Rowv = F, 
          Colv = F, #as.dendrogram(clu), 
          ColSideColors = ColSideColors,
          dendrogram = "none",
          # colsep = colsep,
          # rowsep = rowsep,
          sepwidth = c(0.1, 0.1),
          sepcolor = "white",
          scale = "row",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          margin = c(1, 8),
          na.rm = F);
dev.off(dev.cur());

genes.to.plot <- as.vector(spec.genes)

genes.to.plot <- c("ENSMUSG00000030088|ALDH1L1", "ENSMUSG00000046798|CLDN12",
                   "ENSMUSG00000022037|CLU", "ENSMUSG00000003518|DUSP3",
                   "ENSMUSG00000001768|RIN2", "ENSMUSG00000057329|BCL2",
                   "ENSMUSG00000032902|SLC16A1", "ENSMUSG00000020063|SIRT1", 
                   "ENSMUSG00000020869|LRRC59", "ENSMUSG00000030322|MBD4", 
                   "ENSMUSG00000025869|NOP16", "ENSMUSG00000022508|BCL6",
                   "ENSMUSG00000004558|NDRG2", "ENSMUSG00000026821|RALGDS",
                   "ENSMUSG00000026196|BARD1", "ENSMUSG00000007812|ZFP655",
                   "ENSMUSG00000030316|1500001M20RIK","ENSMUSG00000089857|ZFP882",
                   "ENSMUSG00000017499|CDC6", "ENSMUSG00000022160|METTL3", 
                   "ENSMUSG00000021250|FOS", "ENSMUSG00000055435|MAF",
                   "ENSMUSG00000044647|CSRNP3", "ENSMUSG00000038070|CNTLN",
                   "ENSMUSG00000024317|RNF138",
                   "ENSMUSG00000018171|VMP1", "ENSMUSG00000040022|RAB11FIP2",
                   "ENSMUSG00000002835|CHAF1A",
                   "ENSMUSG00000015488|5930434B04RIK", "ENSMUSG00000026930|GPSM1",
                   "ENSMUSG00000025348|ITGA7", "ENSMUSG00000012640|ZFP715",
                   "ENSMUSG00000052155|ACVR2A", "ENSMUSG00000028243|UBXN2B", 
                   "ENSMUSG00000068663|CLEC16A", "ENSMUSG00000049225|PDP1",
                   "ENSMUSG00000060708|CNO","ENSMUSG00000029681|BCL7B",
                   "ENSMUSG00000049225|PDP1")

pdf(file = paste(Sys.Date(), "E14.5_prog_spec_violin.pdf", sep = "_"),
    width = 10, height = 60, useDingbats = F);
VlnPlot(df0.E14.5.prog, features.plot = genes.to.plot, nCol = 1,
        point.size.use = 0, group.by = "clusters")
dev.off()

pdf(file = paste(Sys.Date(), "E14.5_prog_spec_violin_10-2.pdf", sep = "_"),
    width = 30, height = 60, useDingbats = F);
VlnPlot(df0.E14.5.prog, features.plot = genes.to.plot[1901:2000], #nCol = 1,
        point.size.use = 0, group.by = "clusters")
dev.off()

#### 4.3.3 Use Random Forest to reclassify VZ and SVZ cells ####
cate_bf_fs <- as.factor(cl)
feature_bf_fs <- as.matrix(t(df0.E14.5.prog@scale.data[df0.E14.5.prog@var.genes,]))

rf_bf_fs <- randomForest(feature_bf_fs, cate_bf_fs, importance = TRUE, proximity = TRUE)
imp_bf_fs <- importance(rf_bf_fs, type = 1)

fs <- rfcv(feature_bf_fs, cate_bf_fs, cv.fold = 10, scale = "log", step = 0.9, recursive = F)
len <- length(fs$n.var[fs$error.cv == min(fs$error.cv)])
min_fs <- fs$n.var[fs$error.cv == min(fs$error.cv)][len] #get least features
ind <- order(-imp_bf_fs)[1: min_fs]
feature_fs <- feature_bf_fs[ , ind]
cate_fs <- cate_bf_fs

rf_fs <- randomForest(feature_fs, as.factor(cate_fs), importance=TRUE, proximity=TRUE)
fea1_fs <- data.frame()
fea1_fs <- feature_fs[(rf_fs$predicted == '1') & (rf_fs$votes[ , 1] > 0.6), , drop = FALSE]
cat1_fs <- rf_fs$predicted[(rf_fs$predicted =='1') & (rf_fs$votes[ , 1] > 0.6)]
fea2_fs <- data.frame()
fea2_fs <- feature_fs[(rf_fs$predicted == '2') & (rf_fs$votes[ , 2] > 0.6), , drop = FALSE]
cat2_fs <- rf_fs$predicted[(rf_fs$predicted =='2') & (rf_fs$votes[ , 2] > 0.6)]

cate <- as.factor(c(as.character(cat1_fs), as.character(cat2_fs)))
feature <- as.matrix(rbind(fea1_fs, fea2_fs))

set <- sample(1: nrow(feature), nrow(feature), replace = F)
cate <- cate[set]
feature <- feature[set, ] 

rf_whole <- randomForest(feature, as.factor(cate), importance = TRUE, proximity = TRUE)
pred_whole <- predict(rf_whole, newdata = feature_fs)

dis <- as.dist(1 - cor(df0.E14.5.prog@scale.data[RGs[,1],],use = "pairwise.complete.obs"))
clu <- hclust(dis, method = "average")
cl <- cutree(clu, k=2)

temp <- t(feature_fs)
# dis <- as.dist(1 - cor(temp,use = "pairwise.complete.obs"))
# clu <- hclust(dis, method = "complete")
# cl <- cutree(clu, k=2)

ord <- order(pred_whole)
ColSideColors <- c(#gg_color_hue(3)[as.numeric(as.factor(df0.E14.5.prog@meta.data$region))],
  brewer.pal(3, name = "Paired")[as.numeric(pred_whole)])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- as.matrix(ColSideColors[ord,])
col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

pdf(file = paste(Sys.Date(), "E14.5_progenitor_VZ_SVZ_heatmap_post-RF.pdf", sep = "_"),
    width = 10, height = 14, useDingbats = F);
heatmap.3(temp[,ord],
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8,
          key = F,
          main=NA,
          col = col,
          symkey = F,
          cexRow=0.5, 
          cexCol = 0.6, 
          Rowv = T, 
          Colv = F, #as.dendrogram(clu), 
          ColSideColors = ColSideColors,
          ColSideColorsSize = 0.5,
          dendrogram = "both",
          #scale = "row",
          #colsep = colsep,
          sepcolor = "black",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

df0.E14.5.prog <- AddMetaData(df0.E14.5.prog, metadata = pred_whole, col.name = "VZ.SVZ")
pred_whole_E14.5 <- pred_whole
save.image("~/Data/Da_Mi/2017-11-13_working.RData")

#### 4.3.4 Use random forest ####
cate_bf_fs <- as.factor(df0.E14.5.prog@ident)
feature_bf_fs <- as.matrix(t(df0.E14.5.prog@scale.data[df0.E14.5.prog@var.genes,]))

rf_bf_fs <- randomForest(feature_bf_fs, cate_bf_fs, importance = TRUE, proximity = TRUE)
imp_bf_fs <- importance(rf_bf_fs, type = 1)

fs <- rfcv(feature_bf_fs, cate_bf_fs, cv.fold = 10, scale = "log", step = 0.9, recursive = F)
len <- length(fs$n.var[fs$error.cv == min(fs$error.cv)])
min_fs <- fs$n.var[fs$error.cv == min(fs$error.cv)][len] #get least features
ind <- order(-imp_bf_fs)[1: min_fs]
feature_fs <- feature_bf_fs[ , ind]
cate_fs <- cate_bf_fs

rf_fs <- randomForest(feature_fs, as.factor(cate_fs), importance=TRUE, proximity=TRUE)
fs$n.var

fea1_fs <- data.frame()
fea1_fs <- feature_fs[(rf_fs$predicted == '1') & (rf_fs$votes[ , 1] > 0.6), , drop = FALSE]
cat1_fs <- rf_fs$predicted[(rf_fs$predicted =='1') & (rf_fs$votes[ , 1] > 0.6)]
fea2_fs <- data.frame()
fea2_fs <- feature_fs[(rf_fs$predicted == '2') & (rf_fs$votes[ , 2] > 0.6), , drop = FALSE]
cat2_fs <- rf_fs$predicted[(rf_fs$predicted =='2') & (rf_fs$votes[ , 2] > 0.6)]

cate <- as.factor(c(as.character(cat1_fs), as.character(cat2_fs)))
feature <- as.matrix(rbind(fea1_fs, fea2_fs))

set <- sample(1: nrow(feature), nrow(feature), replace = F)
cate <- cate[set]
feature <- feature[set, ] 

rf_whole <- randomForest(feature, as.factor(cate), importance = TRUE, proximity = TRUE)
pred_whole <- predict(rf_whole, newdata = feature_fs)

dis <- as.dist(1 - cor(df0.E14.5.prog@scale.data[RGs[,1],],use = "pairwise.complete.obs"))
clu <- hclust(dis, method = "average")
cl <- cutree(clu, k=2)

#### 4.3.5 Get markers to do tSNE ####
# E14.5_markers <- read.csv("2017-11-15_E14.5_marker_genes.csv", header = F, stringsAsFactors = F)
# E14.5_markers[,1] <- sapply(E14.5_markers[,1], function(x){paste("\\b", x, "$", sep = "")})
# E14.5_markers[,1] <- as.character(sapply(E14.5_markers[,1],function(x){grep(x, rownames(df0.E14.5.prog@data), value = T)}))
# E14.5_markers <- E14.5_markers[E14.5_markers[,1] != "character(0)",]

E14.5_markers <- read.csv("2017-11-15_E14.5_marker_genes.csv", header = T, stringsAsFactors = F)

df0.E14.5.prog <- RunPCA(df0.E14.5.prog, pc.genes = E14.5_markers[,1],
                         pcs.print = 1:10, pcs.compute = 10, genes.print = 10)

df0.E14.5.prog <-JackStraw(df0.E14.5.prog, num.replicate = 100, do.print = F,num.pc = 10)
JackStrawPlot(df.int, PCs = 1:10)

df0.E14.5.prog <- RunTSNE(df0.E14.5.prog, dims.use = c(1:2),#pcs[[i]], 
                          perplexity = 30,
                          seed.use = 3)

df0.E14.5.prog <- FindClusters(object = df0.E14.5.prog, 
                               reduction.type = "pca", 
                               dims.use = 1:2, 
                               resolution = 3, # dMGE, CGE = 0.8, E14.5 vMGE = 2
                               k.scale = 25, # dMGE, CGE = 50
                               prune.SNN = 0, #dMGE, CGE = 0.5
                               plot.SNN = F, 
                               print.output = F, 
                               save.SNN = F,
                               algorithm = 1,
                               force.recalc = TRUE, 
                               random.seed = 1,
                               k.param = 30) 

p1 <- TSNEPlot(object = df0.E14.5.prog, pt.size = 3, do.return = T)
p2 <-TSNEPlot(object = df0.E14.5.prog, group.by = "VZ.SVZ", do.return = T,
              pt.size = 3)
p3 <-TSNEPlot(object = df0.E14.5.prog, group.by = "region", do.return = T,
              pt.size = 3)

pdf(file = paste(Sys.Date(),"E14.5_prog_tSNE_all_cells2.pdf", sep = "_"), width = 33, 
    height = 10, useDingbats = F);
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

df0.E14.5.prog@meta.data$clusters <- factor(df0.E14.5.prog@meta.data$res.3, 
                                            levels = c(10, 4, 1, 0, 2, 5, 9, 7, 8, 3, 6),
                                            labels = c(1:11))
df0.E14.5.prog@meta.data$clusters <- paste("E14.5", as.character(df0.E14.5.prog@meta.data$clusters))

df0.E14.5.prog <- SetAllIdent(df0.E14.5.prog, "clusters")

markers <- FindAllMarkers(object = df0.E14.5.prog, 
                          only.pos = TRUE, min.pct = 0.25, 
                          return.thresh = 0.01)
top20 <- markers %>% group_by(cluster) %>% top_n(20, avg_diff)

pdf(file = paste(Sys.Date(), "E14.5_prog_heatmap_top20.pdf", sep = "_"),
    width = 20, height = 40, useDingbats = F);
DoHeatmap(object = df0.E14.5.prog,
          genes.use = top20$gene,#rownames(df.markers)[1:100], #genes,
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

temp <- df0.E14.5.prog@scale.data[ E14.5_markers$V1[1:16], ] 
ord <- order(factor(df0.E14.5.prog@meta.data$VZ.SVZ, levels = c(2,1)), 
             df0.E14.5.prog@meta.data$region, runif(ncol(temp), 0,1))
ColSideColors <- c(gg_color_hue(3)[as.numeric(as.factor(df0.E14.5.prog@meta.data$region))],
  brewer.pal(3, name = "Paired")[df0.E14.5.prog@meta.data$VZ.SVZ])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- as.matrix(ColSideColors[ord,])
# colsep <- cumsum(table(df0.E14.5.prog@meta.data$clusters))
# rowsep <- cumsum(table(top20$cluster))

pairs.breaks <- seq(-1.5, 1.5, length.out=1001)

pdf(file = paste(Sys.Date(), "E14.5_VZ_SVZ_color_bar.pdf", sep = "_"),
    width = 10, height = 8, useDingbats = F);
heatmap.3(temp[,ord],
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8,
          key = T,
          main=NA,
          col = col,
          symkey = F,
          cexRow=1, 
          cexCol = 1, 
          Rowv = T, 
          Colv = F, #as.dendrogram(clu), 
          ColSideColors = ColSideColors,
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

RG_markers <- read.csv("RGs markers 16-11-2017.csv", header = F, stringsAsFactors = F)
RG_markers[,1] <- unique(sapply(RG_markers[,1], function(x){paste("\\b", x, "$", sep = "")}))
RG_markers[,1] <- as.character(sapply(RG_markers[,1],function(x){grep(x, rownames(df0.E14.5.prog@data), value = T)}))
RG_markers <- RG_markers[RG_markers[,1] != "character(0)",]


pdf(file = paste(Sys.Date(), "E14.5_prog_violin_plot_region.pdf", sep = "_"),
    width = 20, height = 15, useDingbats = F);
VlnPlot(df0.E14.5.prog, features.plot = #RG_markers[,1],
              c("ENSMUSG00000001496|NKX2-1",
                                              "ENSMUSG00000030551|NR2F2",
                                              "ENSMUSG00000069171|NR2F1",
                                              "ENSMUSG00000041309|NKX6-2",
                                              "ENSMUSG00000048402|GLI2",
                                              "ENSMUSG00000004151|ETV1",
                                              "ENSMUSG00000034227|FOXJ1",
                                              "ENSMUSG00000028201|LHX8",
                                              "ENSMUSG00000027210|MEIS2",
                                              "ENSMUSG00000048562|SP8",
                                              "ENSMUSG00000010175|PROX1"),
            # cols.use = c("grey", "blue"),
        # pt.size = 3
        )
dev.off()

#### 4.4 Analyze E14.5 SVZ cells only ####
df0.E14.5.prog <- SetAllIdent(df0.E14.5.prog, id = "VZ.SVZ")
df0.E14.5.prog.svz <- SubsetData(df0.E14.5.prog, ident.use = 1)

df0.E14.5.prog.svz <- RunPCA(object = df0.E14.5.prog.svz, 
                             pc.genes=E14.5_markers[,1], 
                             pcs.compute = 15, 
                             pcs.print = 1:10,
                             genes.print = 10)

# df0.E14.5.prog.svz <-JackStraw(df0.E14.5.prog.svz, 
#                                num.replicate = 100, do.print = F,num.pc = 15)
# JackStrawPlot(df0.E14.5.prog, PCs = 1:15)

# p1 <- PCAPlot(df0.E14.5.prog.svz, group.by = "region",do.return = T)
# p2 <- PCAPlot(df0.E14.5.prog.svz, group.by = "time_point", do.return = T)
# p3 <- PCAPlot(df0.E14.5.prog.svz, group.by = "IFC",do.return = T)
# pdf(paste(Sys.Date(), "E14.5_prog_PCA.pdf", sep = "_"), width = 26, height = 8, useDingbats = F);
# plot_grid(p1, p2, p3, ncol = 3)
# dev.off()

df0.E14.5.prog.svz <- RunTSNE(df0.E14.5.prog.svz, dims.use = c(1:2),#pcs[[i]], 
                              perplexity = 20,
                              seed.use = 3)

df0.E14.5.prog.svz <- FindClusters(object = df0.E14.5.prog.svz, 
                                   reduction.type = "pca", 
                                   dims.use = 1:2, 
                                   resolution = 2, # dMGE, CGE = 0.8, E14.5 vMGE = 2
                                   k.scale = 25, # dMGE, CGE = 50
                                   prune.SNN = 0, #dMGE, CGE = 0.5
                                   plot.SNN = F, 
                                   print.output = F, 
                                   save.SNN = F,
                                   algorithm = 1,
                                   force.recalc = TRUE, 
                                   random.seed = 1,
                                   k.param = 30) 

p1 <- TSNEPlot(object = df0.E14.5.prog.svz, pt.size = 3, do.return = T)
p2 <-TSNEPlot(object = df0.E14.5.prog.svz, group.by = "region", do.return = T,
              pt.size = 3)
# p3 <-TSNEPlot(object = df0.E14.5.prog.svz, group.by = "", do.return = T,
#               pt.size = 3)
pdf(file = paste(Sys.Date(),"E14.5_prog_svz_tSNE2.pdf", sep = "_"), width = 22, 
    height = 10, useDingbats = F);
plot_grid(p1, p2, 
          # p3, 
          ncol = 2)
dev.off()

markers.temp <- FindAllMarkers(df0.E14.5.prog.svz, min.pct = 0.1, return.thresh = 0.01,only.pos = T)

top20 <- markers.temp %>% group_by(cluster) %>% top_n(20, avg_diff)

pdf(file = paste(Sys.Date(),"E14.5_svz_prog_heatmap_top20.pdf", sep = "_"),
    width = 20, height = 40, useDingbats = F);
DoHeatmap(object = df0.E14.5.prog.svz,
          genes.use = top20$gene,#rownames(dff.markers)[1:100],
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

pdf(file = paste(Sys.Date(), "E14.5_prog_svz_violin_plot_region.pdf", sep = "_"),
    width = 20, height = 15, useDingbats = F);
VlnPlot(df0.E14.5.prog.svz, features.plot = c(RG_markers[,1],
                                              "ENSMUSG00000001496|NKX2-1",
                                              "ENSMUSG00000030551|NR2F2",
                                              "ENSMUSG00000069171|NR2F1",
                                              "ENSMUSG00000041309|NKX6-2",
                                              "ENSMUSG00000048402|GLI2",
                                              "ENSMUSG00000004151|ETV1",
                                              "ENSMUSG00000034227|FOXJ1",
                                              "ENSMUSG00000028201|LHX8",
                                              "ENSMUSG00000027210|MEIS2",
                                              "ENSMUSG00000048562|SP8",
                                              "ENSMUSG00000010175|PROX1"),
        # cols.use = c("grey", "blue"),
        # pt.size = 3
)
dev.off()

#### 4.5 Analyze E14.5 VZ cells only ####
df0.E14.5.prog <- SetAllIdent(df0.E14.5.prog, id = "VZ.SVZ")
df0.E14.5.prog.vz <- SubsetData(df0.E14.5.prog, ident.use = 2)

df0.E14.5.prog.vz <- RunPCA(object = df0.E14.5.prog.vz, 
                             pc.genes=E14.5_markers[,1], 
                             pcs.compute = 15, 
                             pcs.print = 1:10,
                             genes.print = 10)

# df0.E14.5.prog.vz <-JackStraw(df0.E14.5.prog.vz, 
#                                num.replicate = 100, do.print = F,num.pc = 15)
# JackStrawPlot(df0.E14.5.prog, PCs = 1:15)

# p1 <- PCAPlot(df0.E14.5.prog.vz, group.by = "region",do.return = T)
# p2 <- PCAPlot(df0.E14.5.prog.vz, group.by = "time_point", do.return = T)
# p3 <- PCAPlot(df0.E14.5.prog.vz, group.by = "IFC",do.return = T)
# pdf(paste(Sys.Date(), "E14.5_prog_PCA.pdf", sep = "_"), width = 26, height = 8, useDingbats = F);
# plot_grid(p1, p2, p3, ncol = 3)
# dev.off()

df0.E14.5.prog.vz <- RunTSNE(df0.E14.5.prog.vz, dims.use = c(1:2),#pcs[[i]], 
                              perplexity = 20,
                              seed.use = 3)

df0.E14.5.prog.vz <- FindClusters(object = df0.E14.5.prog.vz, 
                                   reduction.type = "pca", 
                                   dims.use = 1:2, 
                                   resolution = 2, # dMGE, CGE = 0.8, E14.5 vMGE = 2
                                   k.scale = 25, # dMGE, CGE = 50
                                   prune.SNN = 0, #dMGE, CGE = 0.5
                                   plot.SNN = F, 
                                   print.output = F, 
                                   save.SNN = F,
                                   algorithm = 1,
                                   force.recalc = TRUE, 
                                   random.seed = 1,
                                   k.param = 30) 

p1 <- TSNEPlot(object = df0.E14.5.prog.vz, pt.size = 3, do.return = T)
p2 <-TSNEPlot(object = df0.E14.5.prog.vz, group.by = "region", do.return = T,
              pt.size = 3)
# p3 <-TSNEPlot(object = df0.E14.5.prog.vz, group.by = "", do.return = T,
#               pt.size = 3)
pdf(file = paste(Sys.Date(),"E14.5_prog_vz_tSNE2.pdf", sep = "_"), width = 22, 
    height = 10, useDingbats = F);
plot_grid(p1, p2, 
          # p3, 
          ncol = 2)
dev.off()

markers.temp <- FindAllMarkers(df0.E14.5.prog.vz, min.pct = 0.1, return.thresh = 0.01,only.pos = T)

top20 <- markers.temp %>% group_by(cluster) %>% top_n(20, avg_diff)

pdf(file = paste(Sys.Date(),"E14.5_vz_prog_heatmap_top20.pdf", sep = "_"),
    width = 20, height = 40, useDingbats = F);
DoHeatmap(object = df0.E14.5.prog.vz,
          genes.use = top20$gene,#rownames(dff.markers)[1:100],
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

pdf(file = paste(Sys.Date(), "E14.5_prog_VZ_violin_plot_region.pdf", sep = "_"),
    width = 20, height = 15, useDingbats = F);
VlnPlot(df0.E14.5.prog.vz, features.plot = c(RG_markers[,1],
                                              "ENSMUSG00000001496|NKX2-1",
                                              "ENSMUSG00000030551|NR2F2",
                                              "ENSMUSG00000069171|NR2F1",
                                              "ENSMUSG00000041309|NKX6-2",
                                              "ENSMUSG00000048402|GLI2",
                                              "ENSMUSG00000004151|ETV1",
                                              "ENSMUSG00000034227|FOXJ1",
                                              "ENSMUSG00000028201|LHX8",
                                              "ENSMUSG00000027210|MEIS2",
                                              "ENSMUSG00000048562|SP8",
                                              "ENSMUSG00000010175|PROX1")
        # cols.use = c("grey", "blue"),
        # pt.size = 3
)
dev.off()

#### 4.6 MetaNeighbor analysis ####
# source("~/Scripts/R/2016-11-03-runMetaNeighbor.R")
source("~/Scripts/R/2017-08-28-runMN-US.R")

data <- MergeSeurat(df0.E12.5.prog, df0.E14.5.prog)
data@meta.data$clusters <- paste(data@meta.data$time_point, data@meta.data$res.3)
data <- FindVariableGenes(object = data, mean.function = ExpMean, 
                  dispersion.function = LogVMR, x.low.cutoff = 0.5, 
                  x.high.cutoff = Inf, y.cutoff = 0.5, y.high.cutoff = 5, 
                  do.plot = F)
data <- ScaleData(data, vars.to.regress = c("nGene", "CC.Difference"))
data <- SetAllIdent(data, id = "clusters")
data.ave <- AverageExpression(data, use.scale = T)

var.genes <- data@var.genes
cell.lab <- colnames(data.ave)
exp.lab <- substr(colnames(data.ave),1, 5)
pheno <- data.frame(Sample_ID = colnames(data.ave), Study_ID = exp.lab, Celltype = cell.lab)

celltype.NV <- run_MetaNeighbor_US(var.genes, data.ave, cell.lab, pheno)
write.csv(celltype.NV, file = "2017-11-26_E12.5_vs_E14.5_progenitor_clusters_MetaNeighbor.csv")

top_hits <- get_top_hits(celltype.NV, pheno, threshold = 0.5, filename = "2017-11-25_neu_P56_MetaNeighbor.txt")

RowSideColors <- rainbow(2, v = 0.8)[as.numeric(as.factor(substr(colnames(celltype.NV), 1, 3)))]

breaks=seq(0,1,length=1001)

pdf(file = paste(Sys.Date(), "Embryonic_progenitor_MetaNeighbor.pdf", sep = "_"),
    width = 10, height = 10, useDingbats = F);
heatmap.3(celltype.NV,
          trace="none",
          density.info="none",
          col=col, 
          RowSideColors = matrix(RowSideColors, nrow = 1),
          breaks=breaks,
          cexRow=0.6,cexCol=0.6,
          margins = c(8, 8))
legend("bottomleft",
       legend = c("E12.5", "E14.5"),
       fill = rainbow(2, v = 0.8),
       yjust = 5)
dev.off()

#### 4.7 MetaNeighbor analysis for E12.5 progenitor ####
data <- SetAllIdent(df0.E12.5.prog, id = "res.3")
data.ave <- AverageExpression(data, use.scale = T)

var.genes <- data@var.genes
cell.lab <- colnames(data.ave)
exp.lab <- "E12.5"
pheno <- data.frame(Sample_ID = colnames(data.ave), Study_ID = exp.lab, Celltype = cell.lab)

celltype.NV <- run_MetaNeighbor_US(var.genes, data.ave, cell.lab, pheno)
write.csv(celltype.NV, file = "2017-11-25_E12.5_progenitor_cluster_MetaNeighbor.csv")

top_hits <- get_top_hits(celltype.NV, pheno, threshold = 0.5, filename = "2017-11-25_neu_P56_MetaNeighbor.txt")

breaks=seq(0,1,length=1001)

pdf(file = paste(Sys.Date(), "E12.5_progenitor_MetaNeighbor.pdf", sep = "_"),
    width = 10, height = 10, useDingbats = F);
heatmap.3(celltype.NV,
          trace="none",
          density.info="none",
          col=col,
          breaks=breaks,
          cexRow=0.6,cexCol=0.6,
          margins = c(8, 8))
# legend("bottomleft",
#        legend = c("E12.5", "E14.5"),
#        fill = rainbow(2, v = 0.8),
#        yjust = 5)
dev.off()


#### 4.8 MetaNeighbor analysis for E14.5 progenitor ####
data <- SetAllIdent(df0.E14.5.prog, id = "res.3")
data.ave <- AverageExpression(data, use.scale = T)

var.genes <- data@var.genes
cell.lab <- colnames(data.ave)
exp.lab <- "E14.5"
pheno <- data.frame(Sample_ID = colnames(data.ave), Study_ID = exp.lab, Celltype = cell.lab)

celltype.NV <- run_MetaNeighbor_US(var.genes, data.ave, cell.lab, pheno)
write.csv(celltype.NV, file = "2017-11-25_E14.5_progenitor_cluster_MetaNeighbor.csv")

top_hits <- get_top_hits(celltype.NV, pheno, threshold = 0.5, 
                         filename = "2017-11-25_E14.5_progenitor_MetaNeighbor.txt")

breaks=seq(0,1,length=1001)

pdf(file = paste(Sys.Date(), "E14.5_progenitor_MetaNeighbor.pdf", sep = "_"),
    width = 10, height = 10, useDingbats = F);
heatmap.3(celltype.NV,
          trace="none",
          density.info="none",
          col=col,
          breaks=breaks,
          cexRow=0.6,cexCol=0.6,
          margins = c(8, 8))
# legend("bottomleft",
#        legend = c("E14.5", "E14.5"),
#        fill = rainbow(2, v = 0.8),
#        yjust = 5)
dev.off()

#### 5. Analyze neurons ####
#### 5.1 Do DE analysis between neurons of 2 time points ####
df0 <- SetAllIdent(df0, id = "pred")
df0.neu <- SubsetData(df0, ident.use = 2)
df0.neu <- SetAllIdent(df0.neu, id = "IFC")
df0.neu <- SubsetData(df0.neu, ident.remove = c("90","189", "190"))

df0.neu <- SetAllIdent(df0.neu, id = "time_point")

markers.temp <- FindAllMarkers(df0.neu, min.pct = 0.1, return.thresh = 0.01, only.pos = T)

top100 <- markers.temp %>% group_by(cluster) %>% top_n(100, avg_diff)

pdf(file = paste(Sys.Date(),"E12.5_vs_E14.5_all_neuron_heatmap_top100.pdf", sep = "_"),
    width = 20, height = 40, useDingbats = F);
DoHeatmap(object = df0.neu,
          genes.use = top100$gene,#rownames(df.markers)[1:100],
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

df0.neu <- FindVariableGenes(object = df0.neu, 
                                   mean.function = ExpMean, 
                                   dispersion.function = LogVMR,
                                   x.low.cutoff = 0.35, 
                                   x.high.cutoff = 8, 
                                   y.cutoff = 0.5,
                                   y.high.cutoff = 4,
                                   do.plot = F)
abline(v = 0.35, col = "red")


genes.to.use <- df0.neu@var.genes
set.seed(1); df0.neu <- RunPCA(object = df0.neu, 
                                     pc.genes = df0.neu@var.genes, 
                                     pcs.compute = 20, 
                                     pcs.print = 1:10, 
                                     genes.print = 10, 
                                     do.print = T)

df0.neu <-JackStraw(df0.neu, num.replicate = 100, 
                          do.print = T, num.pc = 20)
JackStrawPlot(df0.neu, PCs = 1:20)

df0.neu <- RunTSNE(object = df0.neu, dims.use = 1:15, do.fast = T, 
                         perplexity = 15, seed.use = 1)

df0.neu@meta.data$res.1 <- factor(df0.neu@meta.data$res.1,
                                  levels = c("0", "10", "6", "7","3","5","4", "9", "2", 
                                             "1", "8", "11", "12"))
df0.neu <- FindClusters(object = df0.neu, 
                              reduction.type = "tsne", 
                              dims.use = 1:2, 
                              resolution = 1, 
                              k.scale = 25,
                              prune.SNN =  1/15,
                              plot.SNN = F, 
                              print.output = F, 
                              save.SNN = F,
                              algorithm = 1,
                              force.recalc = TRUE, 
                              random.seed = 1,
                              k.param = 30) 

library(scales)
pdf(file = paste(Sys.Date(), "all_neu_lineage_distribution_histogram.pdf", sep = "_"), width = 11, 
    height = 10, useDingbats = F);
ggplot(df0.neu@meta.data, aes(x = res.1, fill = RF_cate_whole)) + 
  geom_bar(position = "fill") + 
  scale_y_continuous(labels = percent_format()) + 
  scale_fill_discrete(name = "Lineages")
dev.off()

p1 <- TSNEPlot(object = df0.neu, pt.size = 3, do.return = T)
p2 <-TSNEPlot(object = df0.neu, group.by = "region", do.return = T, pt.size = 3)
p3 <-TSNEPlot(object = df0.neu, group.by = "time_point", do.return = T, pt.size = 3)

pdf(file = paste(Sys.Date(), "all_neu_tSNE2.pdf", sep = "_"), width = 33, 
    height = 10, useDingbats = F);
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

pdf(file = paste(Sys.Date(), "all_neu_feature_plot.pdf", sep = "_"),
    width = 20, height = 20, useDingbats = F);
FeaturePlot(df0.neu, features.plot = c(#intersect(top20$gene, rownames(df0.E12.5.neu@scale.data))),
  "ENSMUSG00000004366|SST" ,
  "ENSMUSG00000001496|NKX2-1",
  "ENSMUSG00000030551|NR2F2",
  "ENSMUSG00000069171|NR2F1",
  "ENSMUSG00000041309|NKX6-2",
  "ENSMUSG00000048402|GLI2",
  "ENSMUSG00000004151|ETV1",
  "ENSMUSG00000034227|FOXJ1",
  "ENSMUSG00000028201|LHX8",
  "ENSMUSG00000027210|MEIS2",
  "ENSMUSG00000048562|SP8",
  "ENSMUSG00000010175|PROX1"),
  cols.use = c("grey", "blue"),pt.size = 3)
dev.off()

df0.neu <- SetAllIdent(df0.neu, id = "res.1")
df0.neu@meta.data$region <- factor(df0.neu@meta.data$region, levels = c("dMGE", "vMGE", "CGE"))
markers.temp <- FindAllMarkers(df0.neu, min.pct = 0.25, 
                               return.thresh = 0.01,
                               only.pos = T)

top10 <- markers.temp %>% group_by(cluster) %>% top_n(10, avg_diff)
top10 <- top10[order(match(top10$cluster, c("0", "10", "6", "7","3","5","4", "9", "2", "1", "8",
                                            "11", "12"))),]

temp <- df0.neu@scale.data[ top10$gene, ] 
ord <- order(df0.neu@meta.data$res.1, df0.neu@meta.data$time_point, df0.neu@meta.data$region,
             runif(ncol(temp), 0,1))
ColSideColors <- c(gg_color_hue(13)[as.numeric(df0.neu@meta.data$res.1)],
                   gg_color_hue(3)[as.numeric(factor(df0.neu@meta.data$region, 
                                                     levels = c("dMGE", "vMGE", "CGE")))],
                   rainbow(2, s = 0.8)[as.numeric(as.factor(df0.neu@meta.data$time_point))])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- as.matrix(ColSideColors[ord,])

colsep <- cumsum(table(df0.neu@meta.data$res.1))
rowsep <- seq(10, 120, by = 10)
pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

pdf(file = paste(Sys.Date(), "all_neurons_cluster_DE_heatmap.pdf", sep = "_"),
    width = 12, height = 20, useDingbats = F);
heatmap.3(temp[,ord],
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8,
          key = T,
          main=NA,
          col = col,
          symkey = F,
          cexRow= 1, 
          cexCol = 0.6, 
          Rowv = F, 
          Colv = F, #as.dendrogram(clu), 
          ColSideColors = ColSideColors,
          ColSideColorsSize = ,
          dendrogram = "both",
          colsep = colsep,
          rowsep = rowsep,
          sepwidth = c(0.5,0.1),
          sepcolor = "white",
          scale = "row",
          labRow = substr(rownames(temp), 20, 100),
          margins = c(1, 8),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

pdf(file = paste(Sys.Date(),"all_neuron_heatmap_top20.pdf", sep = "_"),
    width = 20, height = 40, useDingbats = F);
DoHeatmap(object = df0.neu,
          genes.use = top20$gene,#rownames(df.markers)[1:100],
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

#### Get GO term for marker genes
library("RDAVIDWebService") ##--Load RDAVIDWebService 
user <- DAVIDWebService$new(email = "zhen.li.zl242@yale.edu", 
                            url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")     ##--Log on to DAVID with email

BP <- list()
for(i in 1:length(unique(markers.temp$cluster))){
  genes <- markers.temp$gene[markers.temp$cluster == unique(markers.temp$cluster)[i]]
  david <- substr(genes, 1, 18)
  ##--Get a list of genes, pay attention to idType
  ##--Read a gene list
  ##--submit gene list to DAVID
  addList(user, inputIds = david, idType = "ENSEMBL_GENE_ID", listType = "Gene", 
          listName = as.character(unique(markers.temp$cluster)[i]))
  
  setAnnotationCategories(user,categories =  c("GOTERM_BP_FAT", "KEGG_PATHWAY"))
  term <- getFunctionalAnnotationChart(user)
  BP[[i]] <- term$Term
}

BP.temp <- sapply(BP, function(xx){xx[ 1:20 ]})
colnames(BP.temp) <- unique( markers.temp$cluster )
write.csv(BP.temp, file = paste(Sys.Date(), "all_neuron_clusters_BP_GO_terms.csv", sep = "_")) 

#### 5.2 Use P56 RF features to classify neurons ####
load("P56.RF.model.Rdata")
P56.RF.features <- read.csv("2017-11-22_P56.RF.features.csv", header = T, stringsAsFactors = F, row.names = 1)
P56.RF.features[,2] <- toupper(P56.RF.features[,1])
P56.RF.features[,2] <- unique(sapply(P56.RF.features[,2], function(x){paste("\\b", x, "$", sep = "")}))
P56.RF.features[,2] <- as.character(sapply(P56.RF.features[,2], function(x){grep(x, rownames(df0.neu@data), value = T)}))

feature_P56 <- t(df0.neu@scale.data[P56.RF.features[,2],])
colnames(feature_P56) <- P56.RF.features[,1]
# temp <- matrix(0, nrow = nrow(feature_P56), ncol = nrow(P56.RF.features) - ncol(feature_P56))
# colnames(temp) <- P56.RF.features[,1][P56.RF.features[,2] == "character(0)"]
# feature_P56 <- cbind(feature_P56, temp)

cate.temp <- predict(rf_fs_P56, newdata = feature_P56, "prob")
ave <- mean(rowMax(cate.temp))
sd <- sd(rowMax(cate.temp))
thresh <- ave

## Find highest probable identity, which has to be above average probability
cate_bf_fs <- unlist(apply(cate.temp, 1, function(xx){
  if(max(xx) >= thresh & length(which(xx == max(xx))) == 1 ){
    colnames(cate.temp)[xx == max(xx)]
  } else {"uncertain"}
}))

keep <- cate_bf_fs %in% names(which(table(cate_bf_fs) >= 5)) &
  cate_bf_fs != "uncertain"

feature_bf_fs <- as.matrix(t(df0.neu@scale.data[df0.neu@var.genes,keep]))
cate_bf_fs <- cate_bf_fs[keep]
cate_bf_fs <- as.factor(as.character(cate_bf_fs))

rf_bf_fs <- randomForest(feature_bf_fs, cate_bf_fs, importance = TRUE, proximity = TRUE)
imp_bf_fs <- importance(rf_bf_fs, type = 1)

fs <- rfcv(feature_bf_fs, cate_bf_fs, cv.fold = 10, scale = "log", step = 0.75, recursive = F)
len <- length(fs$n.var[fs$error.cv == min(fs$error.cv)])
min_fs <- fs$n.var[fs$error.cv == min(fs$error.cv)][len] #get least features
ind <- order(-imp_bf_fs)[1: min_fs]
feature_fs <- feature_bf_fs[ , ind]
cate_fs <- cate_bf_fs

#### validation ####
feature_fs <- read.csv("2018-01-11_RF_features_assigned_neuron.csv", header = T, row.names = 1)
cate_fs <- read.csv("2018-01-11_RF_cate_assigned_neuron.csv", header = T, row.names = 1)

colnames(feature_fs) <- gsub("\\.", "\\|", colnames(feature_fs))

feature_fs <- feature_fs[df0.neu.lin@cell.names,]
cate_fs <- cate_fs[df0.neu.lin@cell.names,]

rf_fs <- randomForest(feature_fs, as.factor(cate_fs), importance=TRUE, proximity=TRUE)

library(parallel)

n <- 9*nrow(feature_fs)/10

ncores <- detectCores() - 1
cl <- makeCluster(ncores, outfile='RF_validationLOG.TXT')
clusterExport(cl, c("n", "feature_fs", "rf_fs"))

cate.temp <- parLapply(cl, 1:10000, function(xx){
  library(randomForest)
  set <- sample(1:nrow(feature_fs), n, replace = F)
  cate.temp <- as.character(predict(rf_fs, feature_fs[set,]))
  cate.temp <- data.frame(sample = rownames(feature_fs)[set], assign = cate.temp, stringsAsFactors = F)
  return(cate.temp)
})
stopCluster(cl)

cate.temp2 <- merge(cate.temp[[1]], cate.temp[[2]], by = "sample", all = T)
for(i in 3:10000){
  cate.temp2 <- merge(cate.temp2, cate.temp[[i]], by = "sample", all = T)
}

cate.temp3 <- apply(cate.temp2[,-1], 1, function(xx){names(table(as.character(xx)))})
cate.temp3 <- data.frame(row.names = cate.temp2$sample, assign = cate.temp3)
cate.temp3 <- cate.temp3[df0.neu.lin@cell.names,]
cate.temp3 == df0.neu.lin@meta.data$RF_cate_whole ## All pass validation

# hist(rowMax(pred_whole))

temp <- t(feature_fs)
ord <- order(cate_fs)

ColSideColors <- c(#gg_color_hue(3)[as.numeric(as.factor(df0.E14.5.prog@meta.data$region))],
  brewer.pal(8, name = "Paired")[as.numeric(as.factor(cate.res))])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- as.matrix(ColSideColors[ord,])
pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

pdf(file = paste(Sys.Date(), "all_neurons_heatmap_post-RF.pdf", sep = "_"),
    width = 10, height = 14, useDingbats = F);
heatmap.3(temp[,ord],
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8,
          key = F,
          main=NA,
          col = col,
          symkey = F,
          cexRow=0.5, 
          cexCol = 0.6, 
          Rowv = T, 
          Colv = F, #as.dendrogram(clu), 
          ColSideColors = ColSideColors,
          dendrogram = "both",
          sepcolor = "black",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

dis <- dist(feature.res,method = "manhattan")
dm <- diffuse(dis, neigen = 2, eps.val = 5000)
dm <- data.frame(dm$X)
colnames(dm) <- c("dm_1", "dm_2")

# pdf(file = paste(Sys.Date(), "all_Da_Mi_int_DM_RF.pdf", sep = "_"),
#     width = 11, height = 10, useDingbats = F);
ggplot(dm, aes(x = dm_1, y = dm_2, 
               color = cate.res)) +
  geom_point() +
  scale_color_manual(name = "Cell Type",
                     values = c(gg_color_hue(7))) +
  theme_bw()
# dev.off()

df0.neu.lin.chodl@meta.data$RF_cate_whole2 <- cate.chodl

df0.neu.lin@meta.data$RF_cate_whole2 <- df0.neu.lin@meta.data$RF_cate_whole
df0.neu.lin@meta.data$RF_cate_whole2[match(df0.neu.lin.chodl@cell.names, df0.neu.lin@cell.names)] <- 
  df0.neu.lin.chodl@meta.data$RF_cate_whole2

#### 5.4 Do DE analysis between new lineage assignments ####
# df0.neu.lin <- SetAllIdent(df0.neu.lin, id = "RF_cate_whole")

markers.RF.NotChodl <- FindAllMarkers(df0.neu.lin.NotChodl, min.pct = 0.25,
                               only.pos = T,
                               return.thresh = 0.01,
                               min.cells = 0)

top40.NotChodl <- markers.RF.NotChodl %>% group_by(cluster) %>% top_n(40, avg_diff)
# top40 <- top40[order(match(top40$cluster, c("Pvalb Gpx3", "Sst Cbln4", "Sst Myh8",
#                                             "Sst Chodl", "Vip Chat", "Vip Parm1",
#                                             "Ndnf Cxcl14"))), ]
top40.NotChodl <- top40.NotChodl[order(match(top40.NotChodl$cluster, 
                                             c("Pvalb Gpx3", "Sst Cbln4", "Sst Myh8",
                                             "Vip Chat", "Vip Parm1", "Ndnf Cxcl14"))), ]

temp <- df0.neu.lin.NotChodl@scale.data[ top40$gene, ]
ord <- order(df0.neu.lin.NotChodl@meta.data$RF_cate_whole, runif(ncol(temp), 0,1))
ColSideColors <- c(#gg_color_hue(3)[as.numeric(as.factor(df0.E14.5.prog@meta.data$region))],
  brewer.pal(7, name = "Paired")[as.numeric(df0.neu.lin.NotChodl@meta.data$RF_cate_whole)])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- as.matrix(ColSideColors[ord,])

colsep <- cumsum(table(df0.neu.lin.NotChodl@meta.data$RF_cate_whole))
rowsep <- cumsum(c(40, 32, 33, 40, 40))
pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

pdf(file = paste(Sys.Date(), "all_neurons_heatmap_post-RF.pdf", sep = "_"),
    width = 10, height = 20, useDingbats = F);
heatmap.3(temp[,ord],
          breaks = pairs.breaks,
          #symbreaks = T,
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
          sepwidth = c(0.5,0.1),
          sepcolor = "white",
          scale = "row",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

#### 5.5 plot average gene expression ####
df0.neu.lin.NotChodl.ave <- AverageExpression(df0.neu.lin.NotChodl)
temp <- df0.neu.lin.NotChodl.ave[top40$gene,c("Pvalb Gpx3", "Sst Cbln4", "Sst Myh8",
                                     "Vip Chat", "Vip Parm1",
                                     "Ndnf Cxcl14")]
temp2 <- matrix(nrow =  nrow(temp))
colsep <- c()
for(i in 1:length(colnames(temp))){
  n <- length(which(df0.neu.lin.NotChodl@meta.data$RF_cate_whole == colnames(temp)[i]))
  res <- rep(temp[,i], times = n)
  colsep <- c(colsep, n)
  res <- matrix(res, nrow = nrow(temp))
  temp2 <-cbind(temp2, res)
}
colsep <- cumsum(colsep)
temp2 <- temp2[,-1]

ColSideColors <- c(#gg_color_hue(3)[as.numeric(as.factor(df0.E14.5.prog@meta.data$region))],
  brewer.pal(7, name = "Paired")[as.numeric(sort(df0.neu.lin.NotChodl@meta.data$RF_cate_whole))])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp2))
ColSideColors <- as.matrix(ColSideColors)

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

pdf(file = paste(Sys.Date(), "all_neurons_heatmap_ave_post-RF.pdf", sep = "_"),
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
          sepwidth = 2,
          sepcolor = "white",
          scale = "row",
          labRow = "",substr(rownames(temp), 20, 100),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

#### 6. MetaNeighbor analysis ####
# source("~/Scripts/R/2016-11-03-runMetaNeighbor.R")
source("~/Scripts/R/2017-08-28-runMN-US.R")

#### 6.1 Between Da Mi and P56 datasets ####
dt1 <- df0.neu.lin.NotChodl@raw.data[df0.neu.lin.NotChodl@var.genes, 
                                     df0.neu.lin.NotChodl@cell.names]
rownames(dt1) <- make.names(substr(rownames(dt1), 20, 100), unique = T)
dt1 <- CreateSeuratObject(raw.data = dt1, 
                          min.cells = 0, min.genes = 0, project = "Da_Mi")
lineage <- paste("Embryonic", df0.neu.lin.NotChodl@meta.data$RF_cate_whole)
dt1@meta.data$lineage <- lineage

dt1 <- NormalizeData(object = dt1, normalization.method = "LogNormalize", 
                     scale.factor = 1000000)
dt1 <- ScaleData(dt1)

dt2 <- df0.allen@raw.data[ df0.allen@var.genes, ]
rownames(dt2) <- toupper(rownames(dt2))
dt2 <- CreateSeuratObject(raw.data = dt2, 
                          min.cells = 0, min.genes = 0, project = "P56")
primary_type <- paste("P56", df0.allen@meta.data$primary_type)
dt2@meta.data$lineage <- primary_type
dt2 <- NormalizeData(object = dt2, normalization.method = "LogNormalize", 
                     scale.factor = 1000000)
dt2 <- ScaleData(dt2)
dt2@meta.data$orig.ident <- "P56"
# dt2 <- SetAllIdent(dt2, "lineage")
# dt2 <- SubsetData(dt2, ident.remove = c("P56 Smad3", "P56 Igtp", "P56 Sncg"))
# dt2 <- SubsetData(dt2, ident.use = c("P56 Pvalb Gpx3", "P56 Sst Cbln4", "P56 Sst Myh8",
#                                         "P56 Vip Chat", "P56 Vip Parm1",
#                                         "P56 Ndnf Cxcl14"))

genes.use <- intersect(substr(df0.neu.lin.NotChodl@var.genes, 20, 100),
                       toupper(df0.allen@var.genes))#c(toupper(P56.RF.features[,1]))

dt <- RunCCA(object = dt1, object2 = dt2, genes.use = genes.use)

p1 <- DimPlot(object = dt, reduction.use = "cca", group.by = "orig.ident", pt.size = 0.5, 
              do.return = TRUE)
p2 <- VlnPlot(object = dt, features.plot = "CC1", group.by = "orig.ident", do.return = TRUE)
plot_grid(p1, p2)

DimHeatmap(object = dt, reduction.type = "cca", cells.use = 500, dim.use = 1:18, 
           do.balanced = TRUE)

dt <- AlignSubspace(dt, reduction.type = "cca", grouping.var = "orig.ident", 
                    dims.align = 1:7)

p1 <- VlnPlot(object = dt, features.plot = "ACC1", group.by = "orig.ident", 
              do.return = TRUE)
p2 <- VlnPlot(object = dt, features.plot = "ACC2", group.by = "orig.ident", 
              do.return = TRUE)
plot_grid(p1, p2)

dt@meta.data$time_point <- substr(dt@meta.data$lineage, 1, 3)
dt@meta.data$time_point[1:244] <- df0.neu.lin.NotChodl@meta.data$time_point

dt <- RunTSNE(object = dt, reduction.use = "cca.aligned", dims.use = 1:7, 
              do.fast = TRUE)
# dt <- FindClusters(object = dt, reduction.type = "cca.aligned", dims.use = 1:15, 
#                    save.SNN = TRUE)
p1 <- TSNEPlot(object = dt,colors.use = c("red", "green", "grey"), group.by = "time_point", do.return = TRUE, pt.size = 2)
p2 <- TSNEPlot(object = dt, group.by = "orig.ident", do.return = TRUE, pt.size = 0.5)
TSNEPlot(object = dt, group.by = "lineage", do.return = TRUE, pt.size = 2)
p3 <- TSNEPlot(object = dt, cells.use = dt@cell.names[1:243], 
               group.by = "lineage", do.return = TRUE, pt.size = 2)
p4 <- TSNEPlot(object = dt, cells.use = dt@cell.names[244:1152], group.by = "lineage", 
               do.return = TRUE, pt.size = 2)

pdf(file = paste(Sys.Date(), "CCA_tSNE.pdf", sep = "_"),
    width = 10, height = 10, useDingbats = F);
plot_grid(p1, p2, p3, p4, nrow = 2, ncol = 2)
dev.off()

data <- t(dt@dr$cca.aligned@cell.embeddings)
data.ave <- sapply(unique(dt@meta.data$lineage),function(xx){
  rowMeans(data[, dt@meta.data$lineage == xx])
})

var.genes <- rownames(data.ave)[1:7]
cell.lab <- colnames(data.ave)
exp.lab <- c(rep("Da Mi", times = 6), rep("Allen", times = 23))
pheno <- data.frame(Sample_ID = colnames(data.ave), Study_ID = exp.lab, Celltype = cell.lab)

celltype.NV <- run_MetaNeighbor_US(var.genes, data.ave, cell.lab, pheno)
write.csv(celltype.NV, paste(Sys.Date(), "celltype_AUORC.csv", sep = "_"))
# top_hits <- get_top_hits(celltype.NV, pheno, threshold = 0.5, filename = "2017-11-25_neu_P56_MetaNeighbor.txt")

RowSideColors <- rainbow(2, v = 0.8)[as.numeric(as.factor(substr(colnames(celltype.NV), 1, 3)))]

breaks=seq(0,1,length=1001)

pdf(file = paste(Sys.Date(), "neu_MetaNeighbor.pdf", sep = "_"),
    width = 10, height = 10, useDingbats = F);
heatmap.3(celltype.NV,
          trace="none",
          density.info="none",
          col=col,
          hclustfun = function(x){hclust(d = dist(x, method = "e"),
                                                 method = "complete")},
          RowSideColors = matrix(RowSideColors, nrow = 1),
          breaks=breaks,
          cexRow=1,cexCol=1,
          margins = c(10, 10))
legend("bottomleft",
       legend = c("Embryonic", "P56"),
       fill = rainbow(2, v = 0.8))
dev.off()

#### 6.2 MetaNeighbor between Da Mi data lineages ####
data <- AverageExpression(df0.neu.lin, genes.use = unique(top40$gene))
var.genes <- rownames(data)
cell.lab <- colnames(data)
exp.lab <- rep("Da Mi", times = 7)
pheno <- data.frame(Sample_ID = colnames(data), Study_ID = exp.lab, Celltype = cell.lab)

celltype.NV <- run_MetaNeighbor_US(var.genes, data, cell.lab, pheno)

pdf(file = paste(Sys.Date(), "Da_Mi_neu_lineage_MetaNeighbor.pdf", sep = "_"),
    width = 10, height = 10, useDingbats = F);
heatmap.3(celltype.NV,
          trace="none",
          density.info="none",
          col=col, 
          breaks=breaks,
          cexRow=0.6,
          cexCol=0.6)
dev.off()

write.csv(celltype.NV, "2017-11-25_Da_Mi_lineage_AUROC.csv")

#### 6.2.2 Use P56 DE to classify lineages ####
# P56.RF.DE <- read.csv("2017-11-26_P56_4types_DE_intersect_genes.csv", 
#                       header = T, stringsAsFactors = F, row.names = 1)
# P56.RF.DE$gene_embryonic <- toupper(P56.RF.DE$gene)
# P56.RF.DE$gene_embryonic <- sapply(P56.RF.DE$gene_embryonic, function(x){paste("\\b", x, "$", sep = "")})
# P56.RF.DE$gene_embryonic <- as.character(sapply(P56.RF.DE$gene_embryonic, function(x){grep(x, rownames(df0.neu@data), value = T)}))
# 
# write.csv(P56.RF.DE, file = "2017-11-26_P56_4types_DE_intersec_genes.csv")

# DoHeatmap(object = df0.neu,
#           genes.use = c("ENSMUSG00000070889|PGAP1", "ENSMUSG00000073678|PGAP1"),#unlist(sapply(P56.RF.DE$gene, 
#                       #       function(x){grep(x, rownames(df0.neu@data), value = T)})),#rownames(df.markers)[1:100],
#           slim.col.label = TRUE,
#           group.spacing = 0.3,
#           remove.key = TRUE)

load("P56.RF.model.Rdata")

P56_fea <- read.csv("P56.RF.features.csv", header = T, 
                   stringsAsFactors = F, row.names = 1)
P56_fea[,1] <- rownames(rf_fs_P56$importance)
P56_fea[,2] <- toupper(P56_fea[,1])
P56_fea[,2] <- sapply(P56_fea[,2], function(x){paste("\\b", x, "$", sep = "")})
P56_fea[,2] <- as.character(sapply(P56_fea[,2], function(x){grep(x, rownames(df0.neu@data), value = T)}))


feature.temp <- t(df0.neu@scale.data[P56_fea[,2],])
colnames(feature.temp) <- P56_fea[,1]

## Find highest probable identity, which has to be above average probability
cate_bf_fs <- as.factor(as.character(predict(rf_fs_P56, newdata = feature.temp, 
                                             norm.votes = T)))
feature_bf_fs <- t(df0.neu@scale.data[df0.neu@var.genes,])

rf_bf_fs <- randomForest(feature_bf_fs, cate_bf_fs, importance = TRUE, proximity = F)

imp_bf_fs <- importance(rf_bf_fs, type = 1)

fs <- rfcv(feature_bf_fs, cate_bf_fs, cv.fold = 10, scale = "log", step = 0.75, recursive = F)
len <- length(fs$n.var[fs$error.cv == min(fs$error.cv)])
min_fs <- fs$n.var[fs$error.cv == min(fs$error.cv)][len] #get least features
ind <- order(-imp_bf_fs)[1: min_fs]
feature_fs <- feature_bf_fs[ , ind]
cate_fs <- cate_bf_fs

rf_fs <- randomForest(feature_fs, as.factor(cate_fs), importance=T, proximity=TRUE)

cate <- apply(rf_fs$votes, 1, function(xx){
  if(length(which(xx == max(xx))) == 1){
    colnames(rf_fs$votes)[xx == max(xx)]
  } else {"uncertain"}
})

ave <- mean(rowMax(rf_fs$votes))
sd <- sd(rowMax(rf_fs$votes))
thresh <- ave - sd

keep <- rowMax(rf_fs$votes) > thresh & cate != "uncertain"

cate <- cate[keep]
feature <- feature[keep,]

rf_whole <- randomForest(feature_fs, as.factor(cate_fs), importance=T, proximity=TRUE)

vote_whole.temp <- predict(rf_fs.temp, newdata = feature, type = "vote")
# hist(rowMax(pred_whole))

cate_whole.temp <- apply(vote_whole.temp, 1, function(xx){
  if(max(xx) > 0.6){
    names(xx)[ xx == max(xx) ]
  } else {"uncertain"}
})

#### 6.3 Find DE between lineages for DR plots ####
df0.neu <- AddMetaData(df0.neu, metadata = cate_whole, col.name = "RF_cate_whole")
df0.neu@meta.data$RF_cate_whole <- factor(as.character(df0.neu@meta.data$RF_cate_whole),
                                          levels = c("Pvalb Gpx3", "Sst Cbln4", "Sst Myh8",
                                                     "Sst Chodl", "Vip Chat", "Vip Parm1",
                                                     "Ndnf Cxcl14"))
df0.neu <- SetAllIdent(df0.neu, id = "RF_cate_whole")
df0.neu.lin <- SubsetData(df0.neu, ident.remove = "uncertain")

markers.temp <- FindAllMarkers(df0.neu.lin, min.pct = 0.25,
                               only.pos = T,
                               return.thresh = 0.01,
                               min.cells = 0)

top40 <- markers.temp %>% group_by(cluster) %>% top_n(40, avg_diff)
top40 <- top40[order(match(top40$cluster, c("Pvalb Gpx3", "Sst Cbln4", "Sst Myh8",
                       "Sst Chodl", "Vip Chat", "Vip Parm1",
                       "Ndnf Cxcl14"))), ]

temp <- df0.neu.lin@scale.data[ top40$gene, ] 
ord <- order(df0.neu.lin@meta.data$RF_cate_whole, runif(ncol(temp), 0,1))
ColSideColors <- c(#gg_color_hue(3)[as.numeric(as.factor(df0.E14.5.prog@meta.data$region))],
  brewer.pal(7, name = "Paired")[as.numeric(df0.neu.lin@meta.data$RF_cate_whole)])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- as.matrix(ColSideColors[ord,])

rowsep <- c(40, 80, 120, 160, 200, 240)
pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

pdf(file = paste(Sys.Date(), "all_neurons_heatmap_post-RF.pdf", sep = "_"),
    width = 10, height = 20, useDingbats = F);
heatmap.3(temp[,ord],
          breaks = pairs.breaks,
          #symbreaks = T,
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
          sepwidth = c(1,1),
          sepcolor = "white",
          scale = "row",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

pdf(file = paste(Sys.Date(),"post-RF_lineage_MAF_MAFB.pdf", sep = "_"),
    width = 20, height = 10, useDingbats = F);
VlnPlot(df0.neu.lin,features.plot = c("ENSMUSG00000055435|MAF", "ENSMUSG00000074622|MAFB"))
dev.off()

lin.genes <- read.csv("Violin plot gene for lineage assignment result 23-11-2017.csv", 
         header = F, stringsAsFactors = F)
lin.genes <- sapply(lin.genes[,1], 
                    function(x){paste("\\b", x, "$", sep = "")})
lin.genes <- as.character(sapply(lin.genes, 
                                 function(x){grep(x, rownames(df0.neu@data),
                                                  value = T)}))


pdf(file = paste(Sys.Date(), "neu_lineage_violin.pdf", sep = "_"),
    width = 20, height = 15, useDingbats = F);
VlnPlot(df0.neu.lin, features.plot = lin.genes)
dev.off()

dt <- data.frame(t(df0.neu.lin@scale.data[unique(lin.genes),]))
colnames(dt) <- substr(colnames(dt), 20, 100)
dt$RF_cate_whole <- df0.neu.lin@meta.data$RF_cate_whole

dt <- melt(dt, id.vars = "RF_cate_whole")

pdf(file = paste(Sys.Date(), "neu_lineage_violin_facet.pdf", sep = "_"),
    width = 10, height = 30, useDingbats = F);
ggplot(dt, aes(x = RF_cate_whole, y = value, fill = RF_cate_whole)) +
  geom_violin(color = "black", size = 0.1, scale = "width") +
  facet_wrap(~variable, 
             nrow = length(unique(dt$variable)),
             strip.position = "right",
             scales = "free_y") +
  theme(strip.text.y = element_text(angle = 0, size = 20),
        axis.text.x = element_text(angle = 45, vjust = 0.5))
dev.off()

#### 7. Neuron ####
df0.neu <- SetAllIdent(df0.neu, id = "time_point")
df0.E12.5.neu <- SubsetData(df0.neu, ident.use = "E12.5")
df0.E12.5.neu@meta.data$region <- factor(df0.E12.5.neu@meta.data$region, 
                                         levels = c("dMGE", "vMGE", "CGE"))

df0.neu <- SetAllIdent(df0.neu, id = "time_point")
df0.E14.5.neu <- SubsetData(df0.neu, ident.use = "E14.5")
df0.E14.5.neu@meta.data$region <- factor(df0.E14.5.neu@meta.data$region, 
                                         levels = c("dMGE", "vMGE", "CGE"))

#### 7.1 E12.5 neuron ####
df0.E12.5.neu <- FindVariableGenes(object = df0.E12.5.neu, 
                                   mean.function = ExpMean, 
                                   dispersion.function = LogVMR, 
                                   x.low.cutoff = 0.5, 
                                   x.high.cutoff = 8, 
                                   y.cutoff = 0.5,
                                   y.high.cutoff = 4,
                                   do.plot = T)
abline(v = 0.5, col = "red")

set.seed(1); df0.E12.5.neu <- RunPCA(object = df0.E12.5.neu, 
                                     pc.genes = df0.E12.5.neu@var.genes, 
                                     pcs.compute = 15, 
                                     pcs.print = 1:10, 
                                     genes.print = 10, 
                                     do.print = T)

df0.E12.5.neu <-JackStraw(df0.E12.5.neu, num.replicate = 100, 
                          do.print = T, num.pc = 15)
JackStrawPlot(df0.E12.5.neu, PCs = 1:15)

df0.E12.5.neu <- RunTSNE(object = df0.E12.5.neu, dims.use = 1:3, do.fast = TRUE, 
                         perplexity = 30, seed.use = 1)

df0.E12.5.neu <- FindClusters(object = df0.E12.5.neu, 
                              reduction.type = "tsne", 
                              dims.use = 1:2, 
                              resolution = 1, 
                              k.scale = 25, 
                              plot.SNN = F, 
                              print.output = F, 
                              save.SNN = F,
                              algorithm = 1,
                              force.recalc = TRUE, 
                              random.seed = 1,
                              k.param = 30) 

p1 <- TSNEPlot(object = df0.E12.5.neu, group.by = "cluster", pt.size = 3, do.return = T)
p2 <-TSNEPlot(object = df0.E12.5.neu, group.by = "region", do.return = T,
              pt.size = 3)
p3 <-TSNEPlot(object = df0.E12.5.neu, group.by = "IFC", do.return = T,
              pt.size = 3)

pdf(file = paste(Sys.Date(), "E12.5_neu_tSNE.pdf", sep = "_"), width = 33, 
    height = 10, useDingbats = F);
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

#### 7.1.2 Analyze E12.5 neurons 3 regions ####
markers.temp <- FindAllMarkers(df0.E12.5.neu, min.pct = 0.1, return.thresh = 0.01, 
                               only.pos = T)

top20 <- markers.temp %>% group_by(cluster) %>% top_n(20, avg_diff)

temp <- df0.E12.5.neu@scale.data[top20$gene,]
ord <- order(as.numeric(df0.E12.5.neu@meta.data$res.1), 
             df0.E12.5.neu@meta.data$region, 
             runif(ncol(temp), 0, 1))
ColSideColors <- c(brewer.pal(3, "Paired")[as.numeric(as.factor(df0.E12.5.neu@meta.data$region))],
                   gg_color_hue(length(unique(df0.E12.5.neu@meta.data$res.1)))[as.numeric(df0.E12.5.neu@meta.data$res.1) + 1])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- ColSideColors[ord,]

col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

colsep <- cumsum(table(as.numeric(df0.E12.5.neu@meta.data$res.1)[ord]))
rowsep <- cumsum(table(top20$cluster))

pdf(file = paste(Sys.Date(), "top20_E12.5_neuron_heatmap_region.pdf", sep = "_"),
    width = 10, height = 20, useDingbats = F);
heatmap.3(temp[,ord],
          breaks = pairs.breaks,
          symbreaks = T,
          keysize = 1,
          key = T,
          main=NA,KeyValueName = NA,
          col = col,
          symkey = T,
          cexRow=0.7, 
          cexCol = 0.6, 
          Rowv = F, 
          Colv = F,#as.dendrogram(clu), 
          ColSideColors = as.matrix(ColSideColors),
          ColSideColorsSize = 3,
          # dendrogram = "both",
          scale = "row",
          colsep = colsep,
          rowsep = rowsep,
          sepcolor = "white",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          lwid = c(1, 5),
          lhei = c(2, 9),
          margins = c(1,6),
          na.rm = F);
dev.off(dev.cur())

df0.E12.5.neu@meta.data$cluster <- df0.E12.5.neu@meta.data$res.1

#### 7.2 E14.5 neuron ####
df0.E14.5.neu <- FindVariableGenes(object = df0.E14.5.neu, 
                            mean.function = ExpMean, 
                            dispersion.function = LogVMR, 
                            x.low.cutoff = 0.5, 
                            x.high.cutoff = 8, 
                            y.cutoff = 0.5,
                            y.high.cutoff = 4,
                            do.plot = T)
abline(v = 0.5, col = "red")

set.seed(1); df0.E14.5.neu <- RunPCA(object = df0.E14.5.neu, 
                              pc.genes = df0.E14.5.neu@var.genes, 
                              pcs.compute = 15, 
                              pcs.print = 1:10, 
                              genes.print = 10, 
                              do.print = T)

df0.E14.5.neu <-JackStraw(df0.E14.5.neu, num.replicate = 100, 
                          do.print = T, num.pc = 15)
JackStrawPlot(df0.E14.5.neu, PCs = 1:15)

df0.E14.5.neu <- RunTSNE(object = df0.E14.5.neu, dims.use = 1:6, do.fast = TRUE, 
                  perplexity = 20, seed.use = 1)

df0.E14.5.neu <- FindClusters(object = df0.E14.5.neu, 
                       reduction.type = "tsne", 
                       dims.use = 1:2, 
                       resolution = 2, 
                       k.scale = 25, 
                       plot.SNN = F, 
                       print.output = F, 
                       save.SNN = F,
                       algorithm = 1,
                       force.recalc = TRUE, 
                       random.seed = 1,
                       k.param = 20) 

p1 <- TSNEPlot(object = df0.E14.5.neu, pt.size = 3, do.return = T)
p2 <-TSNEPlot(object = df0.E14.5.neu, group.by = "region", do.return = T,
              pt.size = 3)
p3 <-TSNEPlot(object = df0.E14.5.neu, group.by = "IFC", do.return = T,
              pt.size = 3)

pdf(file = paste(Sys.Date(), "E14.5_neu_tSNE.pdf", sep = "_"), width = 33, 
    height = 10, useDingbats = F);
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

#### 7.2.2 Analyze E14.5 neurons  ####
markers.temp <- FindAllMarkers(df0.E14.5.neu, min.pct = 0.1, return.thresh = 0.01, 
                               only.pos = T)

top20 <- markers.temp %>% group_by(cluster) %>% top_n(20, avg_diff)

<<<<<<< HEAD
temp <- df0.E14.5.neu@scale.data[top20$gene,]
ord <- order(df0.E14.5.neu@meta.data$res.2, 
             df0.E14.5.neu@meta.data$region, 
             runif(ncol(temp), 0, 1))
ColSideColors <- c(gg_color_hue(3)[as.numeric(as.factor(df0.E14.5.neu@meta.data$region))],
                   brewer.pal(length(unique(df0.E14.5.neu@meta.data$res.2)), "Paired")[as.numeric(df0.E14.5.neu@meta.data$res.2) + 1])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- ColSideColors[ord,]

col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

colsep <- cumsum(table(df0.E14.5.neu@meta.data$res.2[ord]))
rowsep <- cumsum(table(top20$cluster))

pdf(file = paste(Sys.Date(), "temp_E14.5_neuron_heatmap_region.pdf", sep = "_"),
    width = 10, height = 20, useDingbats = F);
heatmap.3(temp[,ord],
          breaks = pairs.breaks,
          symbreaks = T,
          keysize = 1,
          key = T,
          main=NA,KeyValueName = NA,
          col = col,
          symkey = T,
          cexRow=0.7, 
          cexCol = 0.6, 
          Rowv = F, 
          Colv = F,#as.dendrogram(clu), 
          ColSideColors = as.matrix(ColSideColors),
          ColSideColorsSize = 3,
          # dendrogram = "both",
          scale = "row",
          colsep = colsep,
          rowsep = rowsep,
          sepcolor = "white",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          lwid = c(1, 5),
          lhei = c(2, 9),
          margins = c(1,6),
          na.rm = F);
dev.off(dev.cur())

df0.E14.5.neu@meta.data$cluster <- df0.E14.5.neu@meta.data$res.2

#### 4.6 MetaNeighbor analysis ####
# source("~/Scripts/R/2016-11-03-runMetaNeighbor.R")
source("~/Scripts/R/2017-08-28-runMN-US.R")

data <- MergeSeurat(df0.E12.5.neu, df0.E14.5.neu)
data@meta.data$clusters <- paste(data@meta.data$time_point, data@meta.data$cluster)
data <- FindVariableGenes(object = data, mean.function = ExpMean, 
                          dispersion.function = LogVMR, x.low.cutoff = 0.5, 
                          x.high.cutoff = Inf, y.cutoff = 0.5, y.high.cutoff = 5, 
                          do.plot = F)
data <- ScaleData(data, vars.to.regress = c("nGene", "CC.Difference"))
data <- SetAllIdent(data, id = "clusters")
data.ave <- AverageExpression(data, use.scale = T)

var.genes <- df0.neu@var.genes#intersect(df0.E12.5.neu@var.genes, df0.E14.5.neu@var.genes)
cell.lab <- colnames(data.ave)
exp.lab <- substr(colnames(data.ave),1, 5)
pheno <- data.frame(Sample_ID = colnames(data.ave), Study_ID = exp.lab, Celltype = cell.lab)

celltype.NV <- run_MetaNeighbor_US(var.genes, data.ave, cell.lab, pheno)
write.csv(celltype.NV, file = paste(Sys.Date(), "_E12.5_vs_E14.5_neuron_clusters_MetaNeighbor.csv"))

top_hits <- get_top_hits(celltype.NV, pheno, threshold = 0.5, filename = "2017-11-25_neu_P56_MetaNeighbor.txt")

RowSideColors <- rainbow(2, v = 0.8)[as.numeric(as.factor(substr(colnames(celltype.NV), 1, 3)))]

breaks=seq(0,1,length=1001)

pdf(file = paste(Sys.Date(), "Embryonic_neuron_MetaNeighbor.pdf", sep = "_"),
    width = 10, height = 10, useDingbats = F);
heatmap.3(celltype.NV,
          trace="none",
          density.info="none",
          col=col, 
          RowSideColors = matrix(RowSideColors, nrow = 1),
          breaks=breaks,
          cexRow=0.6,cexCol=0.6,
          margins = c(8, 8))
legend("bottomleft",
       legend = c("E12.5", "E14.5"),
       fill = rainbow(2, v = 0.8),
       yjust = 5)
dev.off()

#### 5.2.2 Analyze E12.5 neurons 3 regions ####
df0.neu <- SetAllIdent(df0.neu, id = "time_point")
df0.E14.5.neu <- SubsetData(df0.neu, ident.use = "E14.5")
df0.E14.5.neu@meta.data$region <- factor(df0.E14.5.neu@meta.data$region, 
                                         levels = c("dMGE", "vMGE", "CGE"))
df0.E14.5.neu <- SetAllIdent(df0.E14.5.neu, id = "region")

markers.temp <- FindAllMarkers(df0.E14.5.neu, min.pct = 0.1, return.thresh = 0.01, only.pos = T)

top20 <- markers.temp %>% group_by(cluster) %>% top_n(20, avg_diff)
top20 <- top20[order(match(top20$cluster,c("dMGE", "vMGE", "CGE"))),]

temp <- df0.E12.5.neu@scale.data[top20$gene,]
ord <- order(df0.E12.5.neu@meta.data$region, 
             df0.E12.5.neu@meta.data$pred, 
             df0.E12.5.neu@meta.data$VZ.SVZ,
             runif(ncol(temp), 0, 1))
ColSideColors <- c(gg_color_hue(3)[as.numeric(as.factor(df0.E12.5.neu@meta.data$region))])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- ColSideColors[ord,]

col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

colsep <- cumsum(table(df0.E12.5.neu@meta.data$region))
rowsep <- c(20, 40)

pdf(file = paste(Sys.Date(), "temp_E12.5_neuron_heatmap_region.pdf", sep = "_"),
    width = 10, height = 8, useDingbats = F);
heatmap.3(temp[,ord],
          breaks = pairs.breaks,
          symbreaks = T,
          keysize = 1,
          key = T,
          main=NA,KeyValueName = NA,
          col = col,
          symkey = T,
          cexRow=0.7, 
          cexCol = 0.6, 
          Rowv = F, 
          Colv = F,#as.dendrogram(clu), 
          ColSideColors = as.matrix(ColSideColors),
          ColSideColorsSize = 1,
          # dendrogram = "both",
          scale = "row",
          colsep = colsep,
          rowsep = rowsep,
          sepcolor = "white",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          lwid = c(1, 5),
          lhei = c(2, 9),
          margins = c(1,6),
          na.rm = F);
dev.off(dev.cur());

pdf(file = paste(Sys.Date(), "E14.5_neu_feature_plot.pdf", sep = "_"),
    width = 20, height = 20, useDingbats = F);
FeaturePlot(df0.E14.5.neu, features.plot = c(#intersect(top20$gene, rownames(df0.E12.5.neu@scale.data))),
  "ENSMUSG00000004366|SST" ,
  "ENSMUSG00000001496|NKX2-1",
  "ENSMUSG00000030551|NR2F2",
  "ENSMUSG00000069171|NR2F1",
  "ENSMUSG00000041309|NKX6-2",
  "ENSMUSG00000048402|GLI2",
  "ENSMUSG00000004151|ETV1",
  "ENSMUSG00000034227|FOXJ1",
  "ENSMUSG00000028201|LHX8",
  "ENSMUSG00000027210|MEIS2",
  "ENSMUSG00000048562|SP8",
  "ENSMUSG00000010175|PROX1"),
  cols.use = c("grey", "blue"),pt.size = 3)
dev.off()

#### 8. Monocle ####
#### 8.1 Metaneighbor between progenitor clusters and neuronal lineages ####
neu.prog.int <- intersect(rownames(df0.neu@hvg.info)[1:2000], rownames(df0.prog@hvg.info)[1:2000])

df0.neu.lin2@meta.data$clusters1 <- as.character(df0.neu.lin2@meta.data$RF_cate_whole)
df0.neu.lin2@meta.data$clusters1[which(df0.neu.lin2@meta.data$clusters1 == "Sst Chodl")] <- NA
# df0.neu.lin2@meta.data$clusters1 <- substr(df0.neu.lin2@meta.data$clusters1, 1, 3)
df0.neu.lin2@meta.data$clusters1 <- df0.neu.lin2@meta.data$clusters1
df0.neu.lin2@meta.data$clusters2 <- df0.neu.lin2@meta.data$lineage_assignment
df0.neu.lin2@meta.data$clusters2[which(df0.neu.lin2@meta.data$clusters2 == "Emb Sst Chodl")] <- NA
df0.neu.lin2@meta.data$clusters2 <- substr(as.character(df0.neu.lin2@meta.data$lineage_assignment), 5, 100)

df0.neu.lin2@meta.data$clusters <- "unknown"
df0.neu.lin2@meta.data$clusters[which(df0.neu.lin2@meta.data$clusters1 == df0.neu.lin2@meta.data$clusters2)] <-
  df0.neu.lin2@meta.data$clusters1[which(df0.neu.lin2@meta.data$clusters1 == df0.neu.lin2@meta.data$clusters2)]

df0.neu.lin22 <- df0.neu.lin2
df0.neu.lin22 <- SetAllIdent(df0.neu.lin22, id = "clusters")
df0.neu.lin22 <- SubsetData(df0.neu.lin22, ident.remove = c("unknown"))
# df0.neu.lin22@meta.data$clusters <- substr(df0.neu.lin22@meta.data$clusters, 1, 7)

# df0.monocle <- MergeSeurat(df0.E12.5.prog, df0.E14.5.prog)
df0.monocle <- MergeSeurat(df0.E12.5.prog, df0.neu.lin22)
# df0.monocle <- MergeSeurat(df0.monocle, df0.neu.lin22)
df0.monocle <- SetAllIdent(df0.monocle, id = "clusters")
df0.monocle <- NormalizeData(object = df0.monocle, normalization.method = "LogNormalize", 
                             scale.factor = 1000000)
df0.monocle <- ScaleData(df0.monocle,vars.to.regress = c("nGene", "CC.Difference") )

df0.monocle <- FindVariableGenes(object = df0.monocle,
                                 mean.function = ExpMean,
                                 dispersion.function = LogVMR, x.low.cutoff = 0.5,
                                 x.high.cutoff = Inf, y.cutoff = 0.5, y.high.cutoff = 4,
                                 do.plot = T)
# markers.temp <- FindMarkers(df0.monocle,ident.1 = "Pvalb Gpx3", ident.2 = "Sst Cbln4",
#                             genes.use = df0.monocle@var.genes, only.pos = T,
#                             return.thresh = 0.01)

markers.temp <- FindAllMarkers(df0.monocle, genes.use = df0.monocle@var.genes,
                               min.pct = 0.1, only.pos = T, return.thresh = 0.01)

# abline(v = 0.5, col = "red")
 
df0.monocle.ave <- AverageExpression(df0.monocle, use.scale = T)
# monocle.cor <- cor(df0.monocle.ave[neu.prog.int,], method = "s")

var.genes <- markers.temp$gene #rownames(df0.monocle@hvg.info[1:700,]) #neu.prog.int   #neu.prog.int #markers.temp$gene
cell.lab <- colnames(df0.monocle.ave)
exp.lab <- c(rep("E12.5 Prog", times = 13), rep("Neu", times = 3))
pheno <- data.frame(Sample_ID = colnames(df0.monocle.ave), Study_ID = exp.lab, Celltype = cell.lab)

celltype.NV <- run_MetaNeighbor_US(var.genes, df0.monocle.ave, cell.lab, pheno)
write.csv(celltype.NV, paste(Sys.Date(), "prog_neuron_assignment_AUORC_E12.5_RF_KNN_subtype.csv", sep = "_"))

RowSideColors <- rainbow(3, v = 0.8)[as.numeric(as.factor(substr(colnames(celltype.NV), 1, 3)))]

breaks=seq(0,1,length=1001)

pdf(file = paste(Sys.Date(), "prog_neuron_MetaNeighbor_E12.5.pdf", sep = "_"),
    width = 10, height = 10, useDingbats = F);
heatmap.3(celltype.NV,
          trace="none",
          density.info="none",
          col=col, 
          RowSideColors = matrix(RowSideColors, nrow = 1),
          breaks=breaks,
          cexRow=0.6,cexCol=0.6,
          margins = c(8, 8))
# legend("bottomleft",
#        legend = c("E12.5 prog", "E14.5 prog", "Neuron"),
#        fill = rainbow(3, v = 0.8),
#        yjust = 5)
dev.off()

#### 8.2 monocle analysis on SST and Pvalb lineages ####
library(monocle)
library(igraph)
library(scales)
# source("/home/zl242/Scripts/R/order_cells.R")
df0.monocle.lin <- SubsetData(df0.monocle, ident.use = c("Sst Cbln4", "Pvalb Gpx3", 
                                                         "E12.5 7", "E12.5 5"))
# df0.monocle.lin@meta.data$clusters2 <- "Pva"
# df0.monocle.lin@meta.data$clusters2[df0.monocle.lin@meta.data$clusters %in% c("Emb Sst", "E12.5 3", "E12.5 7", "E14.5 10")] <- "Sst"
pd <- new("AnnotatedDataFrame", data = df0.monocle.lin@meta.data)
rownames(pd) <- df0.monocle.lin@cell.names
fd <- data.frame(gene_ID = rownames(df0.monocle.lin@raw.data), 
                 gene_short_name = substr(rownames(df0.monocle.lin@raw.data), 20, 100))
rownames(fd) <- fd$gene_ID
fd <- new("AnnotatedDataFrame", data = fd)
mm <- newCellDataSet(df0.monocle.lin@raw.data[,df0.monocle.lin@cell.names],
                     phenoData = pd, 
                     featureData = fd,
                     lowerDetectionLimit = 10,
                     expressionFamily = negbinomial.size())
                     # expressionFamily = gaussianff())

# Next, use it to estimate RNA counts
# rpc_matrix <- relative2abs(mm, method = "num_genes")
# rpc_matrix <- as.matrix(rpc_matrix)
# rpc_matrix[is.na(rpc_matrix)] <- 0

# Now, make a new CellDataSet using the RNA counts
# mm <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
#                        phenoData = pd,
#                        featureData = fd,
#                        lowerDetectionLimit = 0.5,
#                        expressionFamily = negbinomial.size())

mm <- detectGenes(mm, min_expr = 1)

expressed_genes <- row.names(subset(fData(mm), num_cells_expressed >= 10))

# # Log-transform each value in the expression matrix.
# L <- log(exprs(mm[expressed_genes,]))
# 
# # Standardize each gene, so that they are all on the same scale,
# # Then melt the data with plyr so we can plot it easily
# melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
# 
# # Plot the distribution of the standardized gene expression values.
# qplot(value, geom = "density", data = melted_dens_df) +
#   stat_function(fun = dnorm, size = 0.5, color = 'red') +
#   xlab("Standardized log(FPKM)") +
#   ylab("Density")

mm <- estimateSizeFactors(mm)
mm <- estimateDispersions(mm)

pData(mm)$Total_mRNAs <- Matrix::colSums(exprs(mm))

mm <- mm[,pData(mm)$Total_mRNAs < 1e6]

upper_bound <- 10^(mean(log10(pData(mm)$Total_mRNAs)) +
                     2*sd(log10(pData(mm)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(mm)$Total_mRNAs)) -
                     2*sd(log10(pData(mm)$Total_mRNAs)))

qplot(Total_mRNAs, data = pData(mm), geom = "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)

diff_test_res <- differentialGeneTest(mm[expressed_genes,], cores = 20, verbose = T, 
                                      fullModelFormulaStr = "~clusters",
                                      reducedModelFormulaStr = "~CC.Difference + num_genes_expressed")
# length(which(diff_test_res$qval < 0.01))

markers.temp <- FindAllMarkers(df0.monocle.lin, genes.use = df0.monocle.lin@var.genes,
                               print.bar = T, min.pct = 0.1, only.pos = T,logfc.threshold = 1,
                               return.thresh = 0.01)
top40.monocle <- markers.temp %>% group_by(cluster) %>% top_n(40, avg_logFC)

# markers.temp <- FindMarkers(df0.monocle.lin, ident.1 = "Emb Sst", ident.2 = "Emb Pva",
#                             genes.use = df0.monocle.lin@var.genes,only.pos = T)

# markers.temp2 <- markers.temp[markers.temp$p_val_adj < 0.01,]

ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
# ordering_genes <- top40.monocle$gene #P56.RF.features$V2 #markers.temp2$gene
# ordering_genes <- rownames(markers.temp)
# ordering_genes <- markers.temp$gene
# ordering_genes <- rownames(df0.monocle.lin@hvg.info[1:1000,])

mm <- setOrderingFilter(mm, ordering_genes)
mm <- reduceDimension(mm, max_components = 2, verbose = T,
                      residualModelFormulaStr = "~CC.Difference + num_genes_expressed", 
                      reduction_method = 'DDRTree')
mm <- orderCells(mm)
# mm <- orderCells(mm,root_state = 3,reverse = F)

pdf(file = paste(Sys.Date(), "projection_State_monocle.pdf", sep = "_"),
    width = 6, height = 7, useDingbats = F);
plot_cell_trajectory(mm, color_by = "State", show_backbone = T)
dev.off()

# cds <- mm["ENSMUSG00000052812|ATAD2B",]
# plot_genes_in_pseudotime(cds)

ddrtree_res <- DDRTree(exprs(mm)[ordering_genes, ],sigma = 0.001, dimensions = 2)

mm@reducedDimS <- ddrtree_res$Z
mm@reducedDimK <- ddrtree_res$Y
adjusted_K <- t(reducedDimK(mm))
dp <- as.matrix(dist(adjusted_K))
#colnames(dp) <- pData(mm)$New_Cell_ID
mm@cellPairwiseDistances <- dp
gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
dp_mst <- minimum.spanning.tree(gp)
mm@minSpanningTree <- dp_mst
mm@dim_reduce_type <- "DDRTree"

# cc_ordering <- extract_ddrtree_ordering(mm, root_cell = 200)
# pData(mm)$Pseudotime <- cc_ordering$pseudo_time

S_matrix <- reducedDimS(mm)
data_df <- data.frame(t(S_matrix[1:2,]))
colnames(data_df) <- c("data_dim_1", "data_dim_2")
rownames(data_df) <- df0.monocle.lin@cell.names#df@cell.names
data_df <- merge(data_df, pData(mm), by.x = "row.names", by.y = "row.names")

dp_mst <- minSpanningTree(mm)
edge_list <- as.data.frame(get.edgelist(dp_mst))
colnames(edge_list) <- c("source", "target")

reduced_dim_coords <- reducedDimK(mm)
ica_space_df <- data.frame(Matrix::t(reduced_dim_coords[1:2, ]))
colnames(ica_space_df) <- c("prin_graph_dim_1", "prin_graph_dim_2")
ica_space_df$sample_name <- row.names(ica_space_df)

edge_df <- merge(ica_space_df, edge_list, by.x = "sample_name", 
                 by.y = "source", all = TRUE)
edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = "source_prin_graph_dim_1", 
                                   prin_graph_dim_2 = "source_prin_graph_dim_2"))
edge_df <- merge(edge_df, ica_space_df[, c("sample_name", 
                                           "prin_graph_dim_1", "prin_graph_dim_2")], by.x = "target", 
                 by.y = "sample_name", all = TRUE)
edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = "target_prin_graph_dim_1", 
                                   prin_graph_dim_2 = "target_prin_graph_dim_2"))

g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2)) + 
  geom_point(aes(color = clusters, alpha = Pseudotime), na.rm = TRUE, size = 1)
g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                 y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                 yend = "target_prin_graph_dim_2"), size = 0.3, linetype = "solid", 
                      na.rm = TRUE, data = edge_df)
g <- g + theme_bw()
g <- g + facet_wrap(~clusters, nrow = 3)
pdf(file = paste(Sys.Date(), "projection.pdf", sep = "_"), width = 9, height = 9, useDingbats = F);
g
dev.off()

data_df_E14.5_neu <- subset(data_df, !(clusters %in% c("E12.5 12", "E12.5 13", "E12.5 7",
                       "Emb Sst Cbln4", "Emb Sst Chodl", "Emb Sst Myh8")) & pred == 2)

data_df.save <- data_df
data_df <- data_df.save
# data_df$Pseudotime <- abs(data_df$Pseudotime - 1.6)

g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2)) + 
  geom_point(aes(color = Pseudotime), na.rm = TRUE, size = 1, alpha = 0.8)
g <- g  + scale_color_gradientn(colours = terrain.colors(10), #limits = c(0, 7.5),
                                # breaks = c(0.1, 5), labels = c("Neuron", "Progenitor"),
                                oob = squish)
g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                 y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                 yend = "target_prin_graph_dim_2"), size = 0.3, linetype = "solid", 
                      na.rm = TRUE, data = edge_df)
g <- g + theme_bw()
#g1<- g1+ facet_wrap(~res.1, nrow = 3)

pdf(file = paste(Sys.Date(), "projection_Pseudotime.pdf", sep = "_"),
    width = 4, height = 3, useDingbats = F);
g
dev.off()

SST_lin <- data_df$clusters %in% c("E12.5 7", "Sst")
g1<- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
g1<- g1+ geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                 y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                 yend = "target_prin_graph_dim_2"), size = 0.3, linetype = "solid", 
                      na.rm = TRUE, data = edge_df)
g1<- g1+ geom_point(data = subset(data_df, !SST_lin),
                    aes(x = data_dim_1, y = data_dim_2), color = "grey90", size = 2)
g1<- g1+ geom_point(data = subset(data_df, SST_lin),
            aes(x = data_dim_1, y = data_dim_2, color = Pseudotime, 
                alpha = Pseudotime), size = 5) 
g1<- g1+ geom_point(data = subset(data_df, SST_lin),
                    aes(x = data_dim_1, y = data_dim_2, color = Pseudotime, 
                        alpha = Pseudotime), size = 5)
  # geom_text(data = subset(data_df,
  #                                clusters %in% c("E12.5 12", "E12.5 13", "E12.5 7",
  #                                                "Emb Sst Cbln4", "Emb Sst Chodl", "Emb Sst Myh8") &
  #                                  pred == 2),
  #           aes(x = data_dim_1, y = data_dim_2, label = Pseudotime)) +
g1<- g1+ scale_color_gradientn(colours = terrain.colors(10), #limits = c(0, 7.5), 
                               breaks = c(2.5, 12.5), labels = c("Progenitor", "Neuron"),
                               oob = squish)
g1<- g1+ scale_alpha_continuous(range = c(0.1, 1), limits = c(0, 5), oob = squish, 
                                guide = F)
g1<- g1+ theme(legend.position = c(0.1, 0.9))
  
pdf(file = paste(Sys.Date(), "projection_sst_lineage.pdf", sep = "_"),
    width = 10, height = 10, useDingbats = F);
g1
dev.off()

g2<- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
g2<- g2+ geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                 y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                 yend = "target_prin_graph_dim_2"), size = 0.3, linetype = "solid", 
                      na.rm = TRUE, data = edge_df)
g2<- g2+ geom_point(data = subset(data_df, SST_lin),
                    aes(x = data_dim_1, y = data_dim_2), color = "grey90", size = 2)
g2<- g2+ geom_point(data = subset(data_df, !SST_lin),
                    aes(x = data_dim_1, y = data_dim_2, color = Pseudotime, 
                        alpha = Pseudotime), size = 5) 
g2<- g2+ geom_point(data = subset(data_df, !SST_lin),
                    aes(x = data_dim_1, y = data_dim_2, color = Pseudotime, 
                        alpha = Pseudotime), size = 5)
# geom_text(data = subset(data_df,
#                                clusters %in% c("E12.5 12", "E12.5 13", "E12.5 7",
#                                                "Emb Sst Cbln4", "Emb Sst Chodl", "Emb Sst Myh8") &
#                                  pred == 2),
#           aes(x = data_dim_1, y = data_dim_2, label = Pseudotime)) +
g2 <- g2 + scale_color_gradientn(colours = terrain.colors(10), #limits = c(0, 7.5), 
                               breaks = c(2.5, 12.5), labels = c("Progenitor", "Neuron"),
                               oob = squish)
g2 <- g2 + scale_alpha_continuous(range = c(0.1,1), limits = c(0, 7.5), oob = squish, 
                                guide = F, breaks = c())
g2 <- g2 + scale_size_continuous(range = c(1,5), guide = F)
g2 <- g2 + theme(legend.position = c(0.1, 0.9))

pdf(file = paste(Sys.Date(), "projection_pv_lineage.pdf", sep = "_"),
    width = 10, height = 10, useDingbats = F);
g2
dev.off()

dt.to.plot <- data.frame(MAF = log10(df0.monocle@raw.data["ENSMUSG00000055435|MAF",] + 1), 
                         clusters = df0.monocle@meta.data$clusters2)

pdf(file = paste(Sys.Date(), "MAF_expression_clusters_major lineage.pdf", sep = "_"),
    width = 20, height = 10, useDingbats = F);
ggplot(data = dt.to.plot,aes(x = clusters, y = MAF)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#### Plot genes in pseudotime ####
data_df_sst <- subset(data_df, (clusters %in% c("E12.5 12", "E12.5 13", "E12.5 7",
                                            "Emb Sst Cbln4", "Emb Sst Chodl", "Emb Sst Myh8")))
df0.monocle.lin <- NormalizeData(df0.monocle.lin, scale.factor = 1000000)
df0.monocle.lin <- ScaleData(df0.monocle.lin,vars.to.regress = c("nGene", "CC.Difference"))
df0.monocle.lin <- FindVariableGenes(object = df0.monocle.lin, mean.function = ExpMean, 
                         dispersion.function = LogVMR, x.low.cutoff = 1, 
                         x.high.cutoff = Inf, y.cutoff = 0.5, y.high.cutoff = 5, 
                         do.plot = T)

expr_sst <- df0.monocle.lin@scale.data[df0.monocle.lin@var.genes,as.character(data_df_sst$New_Cell_ID)]
cor.sst <- apply(expr_sst, 1, function(xx){cor(xx, data_df_sst$Pseudotime)})
cor.sst <- sort(cor.sst, decreasing = F)

dt.to.plot <- data.frame(t(as.matrix(expr_sst[names(cor.sst[1:50]),])))
dt.to.plot$Pseudotime <- data_df_sst$Pseudotime
dt.to.plot <- melt(dt.to.plot,id.vars = "Pseudotime")
dt.to.plot$variable <- substr(dt.to.plot$variable, 20, 100)

pdf(file = paste(Sys.Date(), "pseudotime_correlate_sst.pdf", sep = "_"),
    width = 20, height = 40, useDingbats = F);
ggplot(dt.to.plot, aes(x = Pseudotime, y = value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~variable, ncol = 5)
dev.off()

plot_genes_jitter(mm["ENSMUSG00000055435|MAF",],
                  grouping = "clusters",
                  color_by = "clusters",
                  nrow= 1,
                  ncol = NULL,
                  plot_trend = TRUE)

#### Finding Genes that Change as a Function of Pseudotime ####
pData(mm)$clusters2 <- pData(mm)$clusters
pData(mm)$clusters2[pData(mm)$clusters  %in% c("E12.5 7", "Sst Cbln4")] <- "Sst"
pData(mm)$clusters2[!(pData(mm)$clusters  %in% c("E12.5 7", "Sst Cbln4"))] <- "Pva"
cds.sst <- mm[, pData(mm)$clusters2 == "Sst"]
cds.pva <- mm[, pData(mm)$clusters2 == "Pva"]


diff_test_pv <- differentialGeneTest(cds.pva, fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                      reducedModelFormulaStr = "~CC.Difference + num_genes_expressed", 
                                      cores = 20)

sig_gene_sst <- subset(diff_test_sst, qval < 0.01 & use_for_ordering == T)

pdf(file = paste(Sys.Date(), "pseudotime_sig_gene_pv_heatmap.pdf", sep = "_"),
    width = 20, height = 20, useDingbats = F);
plot_pseudotime_heatmap(cds.sst[c(rownames(sig_gene_pv)),],
                        num_clusters = 3,
                        cores = 20,
                        show_rownames = T)
dev.off()

plot_genes_branched_heatmap(lung[row.names(subset(BEAM_res, qval < 1e-4)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)

#### Heatmap ####
genes.to.plot <- c("ENSMUSG00000004151|ETV1",
"ENSMUSG00000028201|LHX8",
"ENSMUSG00000055639|DACH1",
"ENSMUSG00000000184|CCND2",
"ENSMUSG00000033740|ST18",
"ENSMUSG00000026701|PRDX6",
"ENSMUSG00000029245|EPHA5",
"ENSMUSG00000028926|CDK14",
"ENSMUSG00000055435|MAF",
"ENSMUSG00000004366|SST")

temp <- df0.monocle.ave[genes.to.plot, c("E12.5 5", "E12.5 7")]
temp2 <- matrix(nrow =  nrow(temp))

colsep <- 1
temp2 <- temp2[,-1]

# ColSideColors <- c(#gg_color_hue(3)[as.numeric(as.factor(df0.E14.5.prog@meta.data$region))],
#   brewer.pal(7, name = "Paired")[as.numeric(sort(df0.neu.lin.NotChodl@meta.data$RF_cate_whole))])#,
# ColSideColors <- matrix(ColSideColors, nrow = ncol(temp2))
# ColSideColors <- as.matrix(ColSideColors)

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

pdf(file = paste(Sys.Date(), "heatmap_ave_E12.5_5_7.pdf", sep = "_"),
    width = 10, height = 10, useDingbats = F);
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
          # dendrogram = "both",
          # colsep = colsep,
          sepwidth = 2,
          sepcolor = "white",
          # scale = "row",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          na.rm = F);
dev.off(dev.cur())
