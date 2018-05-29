setwd("~/Data/Da_Mi/")

load("2017-11-23_all_cells_Seurat.Rdata")
load("2017-11-23_all_neuron_Seurat.Rdata")
load("2017-11-19_E14.5_progenitor_Seurat.Rdata")
load("2017-11-19_E12.5_progenitor_Seurat.Rdata")

#### 1. PCA for region and time point of all cells ####
pdf(paste(Sys.Date(), "All_cells_PCA_region.pdf", sep = "_"), 
    width = 8, height = 8, useDingbats = F);
p <- PCAPlot(df0, group.by = "region",do.return = T)
p <- p1 + theme(legend.position = c(0.9, 0.9),
                 legend.background = element_rect(size=0.5,
                                                  linetype = "solid",
                                                  color = "black"))
p
dev.off()

save(df0.neu, file = "2017-11-23_all_neuron_Seurat_object.Rdata")

pdf(paste(Sys.Date(), "All_cells_PCA_time.pdf", sep = "_"), 
    width = 8, height = 8, useDingbats = F);
p <- PCAPlot(df0, group.by = "time_point",do.return = T)
p <- p1 + theme(legend.position = c(0.9, 0.9),
                 legend.background = element_rect(size=0.5,
                                                  linetype = "solid",
                                                  color = "black"))
p
dev.off()

#### 2. RF features between all neurons and all progenitors ####
prog.neu.RF.features <- c("\\bRUNX1T1$","\\bDCX$","\\bGAD2$","\\bGAD1$",
                          "\\bNRXN3$", "\\bDLX6OS1$", "\\bGM14204$", "\\bLHX6$",
                          "\\bMAPK10$", "\\bMAGED1$", "\\bMAPT$", "\\bERBB4$",
                          "\\bGNG2$","\\bCRMP1$","\\bTUBB3$","\\bSLAIN1$",
                          "\\bMLLT11$","\\bGAP43$","\\bGABRG2$","\\b2900011O08RIK$",
                          "\\bMYO10$","\\bGPR98$", "\\bEDNRB$", "\\bDDAH1$",
                          "\\bSLC1A3$","\\bFABP7$", "\\bCCNE2$","\\bHELLS$",
                          "\\bLIG1$","\\bHAT1$","\\bMCM3$","\\bMCM4$","\\bCCND2$",
                          "\\bNASP$","\\bMCM6$","\\bTACC3$","\\bKIF22$","\\bCDC6$",
                          "\\bCLSPN$","\\bMKI67$", "\\bTOP2A$", "\\bRRM2$",
                          "\\bSMC4$", "\\bSKA1$", "\\bNDE1$", "\\bASPM$","\\bCENPE$",
                          "\\bPFAS$", "\\bZFP367$")

prog.neu.RF.features <- as.character(sapply(prog.neu.RF.features,function(x){grep(x, rownames(df0@data), value = T)}))


temp <- df0@scale.data[prog.neu.RF.features,]
ord <- order(df0@meta.data$pred, df0@meta.data$time_point, df0@meta.data$region)
ColSideColors <- c(gg_color_hue(3)[as.numeric(as.factor(df0@meta.data$region))],
                   rainbow(2)[as.numeric(as.factor(df0@meta.data$time_point))],
                   brewer.pal(3, name = "Paired")[df0@meta.data$pred])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- ColSideColors[ord,]

col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

pdf(file = paste(Sys.Date(), "temp_all_cells_neuron_progenitor_heatmap_predicted.pdf", sep = "_"),
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
          Rowv = T, 
          Colv = F,#as.dendrogram(clu), 
          ColSideColors = ColSideColors,
          ColSideColorsSize = 6,
          # dendrogram = "both",
          scale = "row",
          #colsep = colsep,
          sepcolor = "black",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          lwid = c(1, 5),
          lhei = c(2, 9),
          margins = c(1,6),
          na.rm = F);
dev.off(dev.cur());

#### 3. tSNE for E12.5, E14.5 separately ####
load("2017-10-22_E12.5_Seurat_object_cc_regressed.Rdata")
load("2017-10-22_E14.5_Seurat_object_cc_regressed.Rdata")

#### 3.1 tSNE for E12.5 ####
df0 <- SetAllIdent(df0,id = "time_point")
df0.E12.5 <- SubsetData(df0, ident.use = "E12.5")
ord <- match(df.E12.5@cell.names, df0.E12.5@cell.names)
df.E12.5@meta.data$pred <- df0.E12.5@meta.data$pred[ord]
df.E12.5@meta.data$pred[is.na(df.E12.5@meta.data$pred)] <- 1
pdf(paste(Sys.Date(), "E12.5_tSNE_prog_vs_neu.pdf", sep = "_"), 
    width = 8, height = 8, useDingbats = F);
p <- TSNEPlot(df.E12.5, group.by = "pred",pt.size = 3,  do.return = T)
p <- p + theme(legend.position = c(0.9, 0.9),
                legend.background = element_rect(size=0.5,
                                                 linetype = "solid",
                                                 color = "black"))
p
dev.off()

#### 3.1.2 Feature plot for E12.5 ####
grad.genes <- read.csv("E12 gradience gene list.csv", header = F, stringsAsFactors = F)
grad.genes[,1] <- sapply(grad.genes[,1], function(x){paste("\\b", x, "$", sep = "")})
grad.genes[,1] <- as.character(sapply(grad.genes[,1],function(x){grep(x, rownames(df.E12.5@data), value = T)}))

pdf(file = paste(Sys.Date(), "E12.5_feature_plot_prog_vs_neu.pdf", sep = "_"),
    width = 20, height = 40, useDingbats = F);
FeaturePlot(df.E12.5, cols.use = c("grey90", "#003366"),
        features.plot = grad.genes[,1],nCol = 4, pt.size = 2)
dev.off()

#### Scatterplot 3D ####
library(scatterplot3d)
for( i in unique(dt$variable)){
  res <- subset(dt, variable == i)
  res2 <- rbind(res, res)
  res2$z <- rep(c(10, -10), each = nrow(res))
  res2$color <- c(gg_color_hue(length(unique(df.E12.5@ident)))[as.numeric(df.E12.5@ident) + 1],
                  rainbow(3,s = 0.3)[as.numeric(df.E12.5@meta.data$region)])
  scatterplot3d(x = res2$tSNE_1, y = res2$tSNE_2, z = res2$z,zlim = c(-15, 15), 
                pch = 16, color = res2$color, box = F,
                angle = 45)
  points3d(z = )
}
library(rgl)

dt <- data.frame(df0.E14.5.prog@dr$tsne@cell.embeddings)
dt$ident <- df0.E14.5.prog@meta.data$res.3
dt$region <- df0.E14.5.prog@meta.data$region
dt$z1 <- 10
dt$z2 <- -14

pdf(paste(Sys.Date(), "E14.5_progenitor_3D_tSNE_prog.pdf", sep = "_"), 
    width = 8, height = 8, useDingbats = F);
scatter3D(x = dt$tSNE_1, y = dt$tSNE_2, z = dt$z1, zlim = c(-15, 15), 
          colvar = as.integer(dt$ident), bty = "g",
          pch = 16, col = c(gg_color_hue(length(unique(dt$ident)))), colkey = F, 
          ticktype = "detailed", phi = 35, theta = 20)
scatter3D(x = dt$tSNE_1, y = dt$tSNE_2, z = dt$z2, zlim = c(-15, 15), 
          add = T,
          colvar = as.integer(factor(dt$region, levels = c("dMGE", "vMGE", "CGE"))),
          pch = 16, col = c(gg_color_hue(3)), colkey = F)
dev.off()

library(plot3D)
for( i in unique(dt$variable)){
  res <- subset(dt, variable == i)
  res$z1 <- rep(10,  time = nrow(res))
  res$z2 <- rep(-10, time = nrow(res))
  scatter3D(x = res$tSNE_1, y = res$tSNE_2, z = res$z1,zlim = c(-15, 15), 
            col.var = as.integer(res$ident),
            pch = 16, col = c(gg_color_hue(15)), colkey = F, phi = 0)
  points3d(z = )
}

#### 3.2 Feature plot for E12.5 Progenitor ####
grad.genes <- read.csv("E12.5 26 gene list.csv", header = F, stringsAsFactors = F)
grad.genes[,1] <- sapply(grad.genes[,1], function(x){paste("\\b", x, "$", sep = "")})
grad.genes[,1] <- as.character(sapply(grad.genes[,1],function(x){grep(x, rownames(df0.E12.5.prog@data), value = T)}))

pdf(file = paste(Sys.Date(), "E12.5_prog_26_gene.pdf", sep = "_"),
    width = 20, height = 35, useDingbats = F);
FeaturePlot(df0.E12.5.prog, cols.use = c("grey90", "#003366"),
            features.plot = grad.genes[,1],nCol = 4, pt.size = 2)
dev.off()

#### 3.3 Feature plot for E14.5 Progenitor ####
grad.genes <- read.csv("E14.5 26 gene list.csv", header = F, stringsAsFactors = F)
grad.genes[,1] <- sapply(grad.genes[,1], function(x){paste("\\b", x, "$", sep = "")})
grad.genes[,1] <- as.character(sapply(grad.genes[,1],function(x){grep(x, rownames(df0.E14.5.prog@data), value = T)}))

pdf(file = paste(Sys.Date(), "E14.5_prog_26_gene.pdf", sep = "_"),
    width = 20, height = 35, useDingbats = F);
FeaturePlot(df0.E14.5.prog, cols.use = c("grey90", "#003366"),
            features.plot = grad.genes[,1],nCol = 4, pt.size = 2)
dev.off()

#### 3.4 Feature plot for E14.5 ####
grad.genes <- read.csv("E14 gradience gene list.csv", header = F, stringsAsFactors = F)
grad.genes[,1] <- sapply(grad.genes[,1], function(x){paste("\\b", x, "$", sep = "")})
grad.genes[,1] <- as.character(sapply(grad.genes[,1],function(x){grep(x, rownames(df.E14.5@data), value = T)}))

pdf(file = paste(Sys.Date(), "E14.5_feature_plot_prog_vs_neu.pdf", sep = "_"),
    width = 20, height = 40, useDingbats = F);
FeaturePlot(df.E14.5, cols.use = c("grey90", "#003366"), features.plot = grad.genes[,1], pt.size = 2)
dev.off()

#### 4. Region difference ####
#### 4.1 E12.5 progenitor region difference ####

df0.E12.5.prog <- SetAllIdent(df0.E12.5.prog, id = "region")

markers <- FindAllMarkers(object = df0.E12.5.prog, 
                          only.pos = TRUE, min.pct = 0.25, 
                          return.thresh = 0.01)
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_diff) ## Choose top 10 genes based on fold change

temp <- df0.E12.5.prog@scale.data[top10$gene,]
ord <- order(df0.E12.5.prog@meta.data$region)
ColSideColors <- c(gg_color_hue(3)[as.numeric(as.factor(df0.E12.5.prog@meta.data$region))])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- ColSideColors[ord,]

col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

colsep <- cumsum(table(df0.E12.5.prog@meta.data$region))
rowsep <- c(10, 20)

pdf(file = paste(Sys.Date(), "E12.5_progenitor_region_DE.pdf", sep = "_"),
    width = 10, height = 5, useDingbats = F);
heatmap.3(temp[,ord],
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8,
          key = F,
          main=NA,
          col = col,
          symkey = F,
          cexRow=0.8, 
          cexCol = 0.6, 
          Rowv = F, 
          Colv = F,#as.dendrogram(clu), 
          ColSideColors = as.matrix(ColSideColors),
          ColSideColorsSize = 3,
          # dendrogram = "both",
          scale = "row",
          colsep = colsep,
          rowsep = rowsep,
          sepwidth = c(2,0.4),
          sepcolor = "white",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

#### 4.1 E12.5 progenitor region difference ####
df0.E12.5 <- SetAllIdent(df0.E12.5, id = "pred")
df0.E12.5.neu <- SubsetData(df0.E12.5, ident.use = 2)
df0.E12.5.neu <- SetAllIdent(df0.E12.5.neu, id = "region")

markers <- FindAllMarkers(object = df0.E12.5.neu, 
                          only.pos = TRUE, min.pct = 0.25, 
                          return.thresh = 0.01)
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_diff) ## Choose top 10 genes based on fold change

temp <- df0.E12.5.neu@scale.data[top10$gene,]
ord <- order(df0.E12.5.neu@meta.data$region)
ColSideColors <- c(gg_color_hue(3)[as.numeric(as.factor(df0.E12.5.neu@meta.data$region))])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- ColSideColors[ord,]

col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

colsep <- cumsum(table(df0.E12.5.neu@meta.data$region))
rowsep <- c(10, 20)

pdf(file = paste(Sys.Date(), "E12.5_neuron_region_DE.pdf", sep = "_"),
    width = 10, height = 5, useDingbats = F);
heatmap.3(temp[,ord],
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8,
          key = F,
          main=NA,
          col = col,
          symkey = F,
          cexRow=0.8, 
          cexCol = 0.6, 
          Rowv = F, 
          Colv = F,#as.dendrogram(clu), 
          ColSideColors = as.matrix(ColSideColors),
          ColSideColorsSize = 3,
          # dendrogram = "both",
          scale = "row",
          colsep = colsep,
          rowsep = rowsep,
          sepwidth = c(2,0.4),
          sepcolor = "white",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

#### 4.3 E14.5 progenitor region difference ####
df0.E14.5.prog <- SetAllIdent(df0.E14.5.prog, id = "region")

markers <- FindAllMarkers(object = df0.E14.5.prog, 
                          only.pos = TRUE, min.pct = 0.25, 
                          return.thresh = 0.01)
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_diff) ## Choose top 10 genes based on fold change

temp <- df0.E14.5.prog@scale.data[top10$gene,]
ord <- order(df0.E14.5.prog@meta.data$region, runif(n = ncol(temp), 0, 1))
ColSideColors <- c(gg_color_hue(3)[as.numeric(as.factor(df0.E14.5.prog@meta.data$region))])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- ColSideColors[ord,]

col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

colsep <- cumsum(table(df0.E14.5.prog@meta.data$region))
rowsep <- c(10, 20)

pdf(file = paste(Sys.Date(), "E14.5_progenitor_region_DE.pdf", sep = "_"),
    width = 10, height = 5, useDingbats = F);
heatmap.3(temp[,ord],
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8,
          key = F,
          main=NA,
          col = col,
          symkey = F,
          cexRow=0.8, 
          cexCol = 0.6, 
          Rowv = F, 
          Colv = F,#as.dendrogram(clu), 
          ColSideColors = as.matrix(ColSideColors),
          ColSideColorsSize = 3,
          # dendrogram = "both",
          scale = "row",
          colsep = colsep,
          rowsep = rowsep,
          sepwidth = c(2,0.4),
          sepcolor = "white",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

#### 4.4 E14.5 progenitor region difference ####
df0.E14.5 <- SetAllIdent(df0.E14.5, id = "pred")
df0.E14.5.neu <- SubsetData(df0.E14.5, ident.use = 2)
df0.E14.5.neu <- SetAllIdent(df0.E14.5.neu, id = "region")

markers <- FindAllMarkers(object = df0.E14.5.neu, 
                          only.pos = TRUE, min.pct = 0.25, 
                          return.thresh = 0.01)
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_diff) ## Choose top 10 genes based on fold change

temp <- df0.E14.5.neu@scale.data[top10$gene,]
ord <- order(df0.E14.5.neu@meta.data$region)
ColSideColors <- c(gg_color_hue(3)[as.numeric(as.factor(df0.E14.5.neu@meta.data$region))])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- ColSideColors[ord,]

col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

colsep <- cumsum(table(df0.E14.5.neu@meta.data$region))
rowsep <- c(10, 20)

pdf(file = paste(Sys.Date(), "E14.5_neuron_region_DE.pdf", sep = "_"),
    width = 10, height = 5, useDingbats = F);
heatmap.3(temp[,ord],
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8,
          key = F,
          main=NA,
          col = col,
          symkey = F,
          cexRow=0.8, 
          cexCol = 0.6, 
          Rowv = F, 
          Colv = F,#as.dendrogram(clu), 
          ColSideColors = as.matrix(ColSideColors),
          ColSideColorsSize = 3,
          # dendrogram = "both",
          scale = "row",
          colsep = colsep,
          rowsep = rowsep,
          sepwidth = c(2,0.4),
          sepcolor = "white",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

#### 4.1 E12.5 progenitor region difference ####

df0.E12.5.prog <- SetAllIdent(df0.E12.5.prog, id = "region")

markers <- FindAllMarkers(object = df0.E12.5.prog, 
                          only.pos = TRUE, min.pct = 0.25, 
                          return.thresh = 0.01)
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_diff) ## Choose top 10 genes based on fold change


df.E12.5@ident  <- factor(df.E12.5@ident, 
                          levels = c(paste("P", 1:7, sep = ""), paste("N", 1:8, sep = ""))) 

p1 <- TSNEPlot(df.E12.5, group.by = "region",do.return = T)
p1 <- p1 + theme(text = element_text(size = 20),
                 line = element_line(size = 1),
                 panel.grid = element_blank(),
                 panel.border = element_rect(fill = NA, colour = "black"),
                 legend.key = element_blank(),
                 legend.position = "",#c(0.9, 0.1),
                 legend.background = element_blank())
p2 <- TSNEPlot(df.E12.5, do.return = T)
p2 <- p2 + theme(text = element_text(size = 20),
                 line = element_line(size = 1),
                 panel.grid = element_blank(),
                 panel.border = element_rect(fill = NA, colour = "black"),
                 legend.key = element_blank(),
                 legend.position = "",#c(0.9, 0.1),
                 legend.background = element_blank())
p3 <- TSNEPlot(df.E12.5, group.by = "cc.state",do.return = T)
p3 <- p3 + theme(text = element_text(size = 20),
                 line = element_line(size = 1),
                 panel.grid = element_blank(),
                 panel.border = element_rect(fill = NA, colour = "black"),
                 legend.key = element_blank(),
                 legend.position = "",#c(0.9, 0.1),
                 legend.background = element_blank())
pdf(paste(Sys.Date(), "E12.5_tSNE.pdf", sep = "_"), width = 13, height = 4, useDingbats = F);
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

plot.genes <- c("Sst", "Maf", "Mafb", "Elfn1", "Neto1", "Chl1", "Npas1", "Npas3", "Reln", "Rnd3", "Satb1",
"Arx", "Dgkg", "Tacr1", "Unc5c", "Rarb", "Mest", "Tle4", "Cdh4", "Gse1", "Syt6", "Npy", 
"Nos1", "Calb1")
plot.genes <- toupper(plot.genes)
plot.genes <- unique(sapply(plot.genes, function(x){paste("\\b", x, "$", sep = "")}))
plot.genes <- as.character(unlist(sapply(plot.genes, function(x){grep(x, rownames(df@data), 
                                                         value = T)})))

pdf(paste(Sys.Date(), "E12.5_features_tSNE.pdf", sep = "_"), width = 20, height = 30, useDingbats = F);
FeaturePlot(object = df,features.plot = plot.genes,nCol = 4,pt.size = 3, 
                 cols.use = c("grey", "blue"))
dev.off()

pdf(paste(Sys.Date(), "genes_E14.5_tSNE.pdf", sep = "_"), width = 20, height = 30, useDingbats = F);
FeaturePlot(object = df.E14.5,features.plot = plot.genes,nCol = 4,pt.size = 3, 
            cols.use = c("grey", "blue"))
dev.off()


df.E12.5 <- SetAllIdent(df.E12.5, id = "")
df.E12.5@ident  <- factor(df.E12.5@ident, 
                          levels = c(paste("P", 1:7, sep = ""), paste("N", 1:8, sep = ""))) 

df.12.markers <- FindAllMarkers(genes.use = df.E12.5@var.genes, 
                              object = df.E12.5,
                              only.pos = TRUE, 
                              min.pct = 0.1, 
                              return.thresh = 0.01)

top10 <- df.12.markers %>% group_by(cluster) %>% top_n(10, avg_diff)

pdf(file = paste(Sys.Date(), 
                 "clusters_heatmap_top10.pdf", sep = "_"), 
    width = 20, height = 10, useDingbats = F);
DoHeatmap(object = df.E12.5, 
          genes.use = top10.12$gene, 
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

df.E14.5 <- SetAllIdent(df.E14.5, id = "clusters")
df.E14.5@ident  <- factor(df.E14.5@ident, 
                          levels = c(paste("P", 1:7, sep = ""), paste("N", 1:6, sep = ""))) 

p1 <- TSNEPlot(df.E14.5, group.by = "region",do.return = T)
p1 <- p1 + theme(text = element_text(size = 20),
                 line = element_line(size = 1),
                 panel.grid = element_blank(),
                 panel.border = element_rect(fill = NA, colour = "black"),
                 legend.key = element_blank(),
                 legend.position = "",#c(0.9, 0.1),
                 legend.background = element_blank())
p2 <- TSNEPlot(df.E14.5,do.return = T)
p2 <- p2 + theme(text = element_text(size = 20),
                 line = element_line(size = 1),
                 panel.grid = element_blank(),
                 panel.border = element_rect(fill = NA, colour = "black"),
                 legend.key = element_blank(),
                 legend.position = "",#c(0.9, 0.1),
                 legend.background = element_blank())
p3 <- TSNEPlot(df.E14.5, group.by = "cc.state",do.return = T)
p3 <- p3 + theme(text = element_text(size = 20),
                 line = element_line(size = 1),
                 panel.grid = element_blank(),
                 panel.border = element_rect(fill = NA, colour = "black"),
                 legend.key = element_blank(),
                 legend.position = "",#c(0.9, 0.1),
                 legend.background = element_blank())
pdf(paste(Sys.Date(), "E14.5_tSNE.pdf", sep = "_"), width = 13, height = 4, useDingbats = F);
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

df.14.markers <- FindAllMarkers(#genes.use = df.E14.5@var.genes, 
                                object = df.E14.5, only.pos = TRUE, 
                                min.pct = 0.1, 
                                return.thresh = 0.01)

top10.12 <- df.12.markers %>% group_by(cluster) %>% top_n(10, avg_diff)

ord <- order(df.E12.5@ident)

# gene.ord <- order(factor(top10.12$cluster, 
#                          levels = c(paste("P", 1:7, sep = ""), paste("N", 1:8, sep = ""))))
                           #c(5, 12, 10, 4, 2, 14, 7, 3, 9, 1, 6, 8, 13, 11, 0)))

gene.ord <- order(factor(top10.12$cluster, 
                         levels = c(paste("P", 1:7, sep = ""), paste("N", 1:8, sep = ""))))

# temp <- df.E12.5@scale.data[top10.12$gene[gene.ord], ]
# 
# temp <- sapply(sort(unique(df.E12.5@meta.data$clusters)), function(xx){
#   rowMeans(temp[, df.E12.5@meta.data$clusters %in% xx])
# })

df.E12.5 <- SetAllIdent(df.E12.5, id = "clusters")
df.E12.5@ident <- factor(df.E12.5@ident, levels = c(paste("P", 1:7, sep = ""), 
                                                    paste("N", 1:8, sep = "")))
names(df.E12.5@ident) <- df.E12.5@cell.names
mean.E12.5 <- AverageExpression(df.E12.5, return.seurat = T)

df.E14.5 <- SetAllIdent(df.E14.5, id = "clusters")
df.E14.5@ident <- factor(df.E14.5@ident, levels = c(paste("P", 1:7, sep = ""), 
                                                    paste("N", 1:6, sep = "")))
names(df.E14.5@ident) <- df.E14.5@cell.names
mean.E14.5 <- AverageExpression(df.E14.5, return.seurat = T)

# markers <- read.csv("./Rubenstein PGs vs Neuruon markers 08-11-2017.csv", header = F, stringsAsFactors = F)
markers <- read.csv("./MZ vs SVZ gene list 08-11-2017.csv", header = F, stringsAsFactors = F)
markers[,1] <- toupper(unlist(markers[, 1 ]))
markers[,1] <- unique(sapply(markers[,1], function(x){paste("\\b", x, "$", sep = "")}))
markers[,1] <- as.character(sapply(markers[,1],function(x){grep(x, rownames(df.E14.5@data), value = T)}))
markers <- markers[markers[,1] != "character(0)",]

markers <- FindAllMarkers(object = df.E12.5, only.pos = TRUE, 
                          min.pct = 0.1, 
                          return.thresh = 0.01)
top20.temp <- markers %>% group_by(cluster) %>% top_n(20, avg_diff)

temp <- mean.E12.5@scale.data[top20.temp$gene,]
# temp <- temp[,c(paste("P", 1:7, sep = ""), paste("N", 1:8, sep = ""))]
# colnames(temp) <- sort(unique(df.E12.5@meta.data$clusters))

ColSideColors <- c(gg_color_hue(length(unique(df.E12.5@ident))))#,
                   #rainbow(3)[as.numeric(as.factor(df.E12.5@meta.data$region))])
# RowSideColors <- rainbow(2)[as.numeric(as.factor(markers[,2]))]
# ColSideColors <- matrix(data = ColSideColors, nrow = ncol(temp));
# ColSideColors <- ColSideColors[ord, ];
# 
# colsep <- lapply(unique(df.E12.5@meta.data$clusters[ord]), 
#                  function(x){length(which(df.E12.5@meta.data$clusters == x))});
# colsep <- cumsum(unlist(colsep));  

pairs.breaks <- seq(-2.5, 2.5, length.out=1001);

col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)
# col.low = "#FF00FF", col.mid = "#000000", 
# col.high = "#FFFF00"

# temp2 <- temp
# temp2[temp2 < 0] <- NA
pdf(file = paste(Sys.Date(), "E12.5_heatmap_top_20.pdf", sep = "_"),
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
          Rowv = F, 
          Colv = F,#as.dendrogram(hca), 
          ColSideColors = as.matrix(ColSideColors),
          ColSideColorsSize = 0.5,
          # RowSideColors = t(as.matrix(RowSideColors)),
          dendrogram = "column",
          #scale = "row",
          #colsep = colsep,
          sepcolor = "black",
          labRow = "",#substr(rownames(temp), 20, 100),
          #labCol = "",
          na.rm = F)
dev.off(dev.cur())

