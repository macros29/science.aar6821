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

setwd("~/Data/Da_Mi/")
load("working_Da_Mi.RData")

temp <- read.csv("E12_mouse_metadata_with_phase_lineage.csv", header = T)
rownames(temp) <- temp$X
temp <- temp[,-1]
df.E12.5@meta.data <- temp

expr <- read.table("./raw_data/mouseBrain.Mida.scRNAseq_cleanData.gene.count.2Run.txt",row.names = 1, header = T)
colnames(expr) <- gsub("90_", "90.", colnames(expr))
colnames(expr) <- gsub("91_", "91.", colnames(expr))
colnames(expr) <- unlist(lapply(colnames(expr),function(x){strsplit(x, split = "_")[[1]][1]}))

GB$IFC <- sapply(rownames(GB), function(x){strsplit(x, split = "\\.")[[1]][2]})
GB$region <- "dMGE"
GB$region[ GB$IFC %in% c("160", "161", "162", "172", "175", "179", "182", 
                         "186", "187", "189", "190") ] <- "CGE"
GB$region[ GB$IFC %in% c("171", "174","177", "178", "180", "183", "185", 
                         "194", "196", "91", "93", "96", "98") ] <- "vMGE"
GB$time_point <- "E14.5"
GB$time_point[GB$IFC %in% c(163:187, 195:196)] <- "E12.5"
GB$mouse <- "170_171_172"
GB$mouse[GB$IFC %in% c(173, 174, 175)] <- "173_174_175"
GB$mouse[GB$IFC %in% c(176, 177, 179)] <- "176_177_179"
GB$mouse[GB$IFC %in% c(180, 181, 182)] <- "180_181_182"
GB$mouse[GB$IFC %in% c(183, 184)] <- "183_184"
GB$mouse[GB$IFC %in% c(185, 186, 187)] <- "185_186_187"
GB$mouse[GB$IFC %in% c(195, 196)] <- "195_196"
GB$mouse[GB$IFC %in% c(160, 161, 162)] <- "160_161_162"
GB$mouse[GB$IFC %in% c(189, 190)] <- "189_190"
GB$mouse[GB$IFC %in% c(192, 191)] <- "192_191"
GB$mouse[GB$IFC %in% c(193, 194)] <- "193_194"
GB$mouse[GB$IFC %in% c(90, 91)] <- "90_91"
GB$mouse[GB$IFC %in% c(92, 93)] <- "92_93"
GB$mouse[GB$IFC %in% c(96, 97)] <- "96_97"
GB$mouse[GB$IFC %in% c(98, 99)] <- "98_99"
GB$mouse[GB$IFC %in% c(101)] <- "101"

cc.genes <- read.table("regev_lab_cell_cycle_genes.txt", header = F)$V1
cc.genes <- paste("\\b", cc.genes, "$", sep = "")
cc.genes <- c(unlist(sapply(cc.genes,function(x){grep(x, rownames(df@scale.data), value = T)})))

s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

samples <- list(c("dMGE","E12.5"), c("vMGE","E12.5"), c("CGE","E12.5"),
                c("dMGE","E14.5"), c("vMGE","E14.5"), c("CGE","E14.5"))

keep <- !(GB$IFC %in% c("180","181","191","192","193","194","195", "196")) &
  !(GB$New_Cell_ID %in% c("C1.173.A8",  "C1.173.A9",  "C1.173.B12", "C1.173.C7",  "C1.184.A2",
                          "C1.184.A3",  "C1.184.A8", "C1.184.D1",  "C1.176.A8", "C1.90.H11", 
                          "C1.92.D4", "C1.92.A2", "C1.92.E6", "C1.97.A5","C1.184.A1",  "C1.184.E10",
                          "C1.101.A6", "C1.101.A7", "C1.101.A8","C1.98.C5","C1.98.C9",
                          "C1.171.A4", "C1.174.A1", "C1.177.E6","C1.177.A12", "C1.177.C12",
                          "C1.98.D4","C1.98.D5","C1.98.E11","C1.98.E7","C1.98.F10",
                          "C1.182.A3",  "C1.186.F3",  "C1.187.D10",
                          "C1.98.F7","C1.98.G11", "C1.91.D11", "C1.91.E11", "C1.91.E4",
                          "C1.91.H10", "C1.93.B12", "C1.93.C6","C1.93.C9", "C1.93.E10",
                          "C1.93.E8","C1.93.G5","C1.93.H1","C1.93.H5","C1.96.B6","C1.96.C6",
                          "C1.98.A3", "C1.98.A7", "C1.98.A8", "C1.93.A7", "C1.93.A6",
                          "C1.160.B4","C1.160.F12", "C1.161.B9", "C1.161.F12","C1.162.C6",
                          "C1.162.G6","C1.190.2.E4"))
df0 <- CreateSeuratObject(raw.data = as.vector(expr[,keep]), min.cells = 5, min.genes = 1000, project = "Da_Mi")
df0 <- NormalizeData(object = df0, normalization.method = "LogNormalize", scale.factor = 10000)
mito.genes <- grep(pattern = "\\bMT-", x = rownames(x = df0@data), value = TRUE)
percent.mito <- colSums(df0@raw.data[mito.genes, ])/colSums(df0@raw.data)
# AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
df0 <- AddMetaData(object = df0, metadata = percent.mito, col.name = "percent.mito")
df0@meta.data <- cbind(df0@meta.data, GB[df0@cell.names,])
#df0 <- ScaleData(df0)
df0 <- ScaleData(df0)

df0 <- FindVariableGenes(object = df0, mean.function = ExpMean,
                         dispersion.function = LogVMR, x.low.cutoff = 0.0125,
                         x.high.cutoff = 3, y.cutoff = 0.5, do.plot = F)

#### Check PCA for variation ####
df0 <- RunPCA(object = df0, pc.genes=df0@var.genes, pcs.compute = 20, pcs.print = 1:12,genes.print = 10)

pdf(file = paste(Sys.Date(), "PCA_pairs.pdf", sep = "_"), 
    width = 40, height = 40, useDingbats = F);
pairs(df0@dr$pca@cell.embeddings[,1:10], col = as.factor(df0@meta.data$mouse), pch = 16, cex = 1)
dev.off()

## Plot correlation of PCs and factors
rec <- apply( df0@meta.data[,4:ncol(df0@meta.data)], 2, function( xx ){apply( df0@dr$pca@cell.embeddings, 2, function( yy ){cor( as.numeric(as.factor(xx)), yy, method = "s" )} )} )
rec[abs(rec) < 0.3] <- NA
heatmap.3(rec, Colv = F, Rowv = F,na.color = "black")

pdf("All_cells_PCA_region.pdf", width = 11, height = 10, useDingbats = F);
PCAPlot(df0, group.by = "region")
dev.off()

df0@meta.data$PC.to.regress2 <- df0@dr$pca@cell.embeddings[,1]
df0 <- ScaleData(df0, vars.to.regress = c("PC.to.regress", "PC.to.regress2"))
df0 <- FindVariableGenes(object = df0, mean.function = ExpMean,
                         dispersion.function = LogVMR, x.low.cutoff = 0.0125,
                         x.high.cutoff = 3, y.cutoff = 0.5, do.plot = F)
df0 <- RunPCA(object = df0, pc.genes=df0@var.genes, pcs.compute = 20, pcs.print = 1:12,genes.print = 10)

#### Cell cycle correction ####
cc_genes <- read.csv("cell_cycle_genes.csv", header = F, stringsAsFactors = F)
res <- sapply(cc_genes[,1], function(xx){paste("\\b",xx,"$",sep = "")})
res <- sapply(res, function(xx){grep(xx, rownames(df0@scale.data), value = T)})
cc_genes[,3] <- as.character(res)
g1s.genes <- cc_genes[,3][cc_genes$V2 == "G1/S"]
g2m.genes <- cc_genes[,3][cc_genes$V2 == "G2/M"]
df0 <- CellCycleScoring(df0, s.genes = g1s.genes, g2m.genes = g2m.genes, set.ident = T)
df0@meta.data$CC.Difference <- df0@meta.data$S.Score - df0@meta.data$G2M.Score

df0 <- RunPCA(object = df0, pc.genes=cc_genes[,3],
             pcs.compute = 30, pcs.print = 1:12, genes.print = 10)
PCAPlot(df0)

#### Regress out PCs and cell cycle
df0 <- ScaleData(object = df0, 
                 vars.to.regress = c("PC.to.regress", "PC.to.regress2", "CC.Difference"),                display.progress = FALSE)
df0.E12.5@meta.data$cc.state <- as.character(df0@ident)
df0@meta.data$cc.state[df0@meta.data$cc.state == "G1"] <- "postmitotic"
df0@meta.data$cc.state[df0@meta.data$cc.state == "S"] <- "G1S"

#### Analyze data ####
df0 <- RunPCA(object = df0, pc.genes=df0@var.genes,
              pcs.compute = 30, pcs.print = 1:12, genes.print = 10)
p1 <- PCAPlot(df0, group.by = "region",do.return = T)
p2 <- PCAPlot(df0,do.return = T)
p3 <- PCAPlot(df0, group.by = "time_point",do.return = T)
pdf("All_cells_PCA.pdf", width = 22, height = 20, useDingbats = F);
plot_grid(p1, p2, p3, ncol = 2)
dev.off()

df0 <-JackStraw(df0, num.replicate = 100, do.print = F,num.pc = 20)

pdf(file = paste(Sys.Date(), "all_cells_JackStraw.pdf", sep = "_"), width = 11,
    height = 10, useDingbats = F);
JackStrawPlot(df0, PCs = 1:20)
dev.off()

df0 <- RunTSNE(df0, dims.use = 1:14,#pcs[[i]], 
               perplexity = 20,
               seed.use = 3)

# df@ident <- as.factor(df@meta.data$region)
# names(df@ident) <- df@cell.names

df0 <- FindClusters(object = df0, 
                    reduction.type = "tsne", 
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

p1 <- TSNEPlot(object = df0, pt.size = 3, do.return = T)
p2 <-TSNEPlot(object = df0, group.by = "time_point", do.return = T,
              pt.size = 3)
p3 <-TSNEPlot(object = df0, group.by = "region", do.return = T,
              pt.size = 3)

pdf(file = paste(Sys.Date(), "all_cells_tSNE.pdf", sep = "_"), width = 33, 
    height = 10, useDingbats = F);
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

df0.markers <- FindAllMarkers(genes.use = df0@var.genes, 
                              object = df0, only.pos = TRUE, 
                              min.pct = 0.25, 
                              return.thresh = 0.01)
#df.markers <- FindMarkers(df, ident.1 = "vMGE", ident.2 = "dMGE",only.pos = T,min.pct = 0.5)

top20 <- df0.markers %>% group_by(cluster) %>% top_n(20, avg_diff)

pdf(file = paste(Sys.Date(), "all_cells_clusters_heatmap_top10.pdf", sep = "_"),
    width = 20, height = 40, useDingbats = F);
DoHeatmap(object = df0,
          genes.use = top20$gene,#rownames(df.markers)[1:100],
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

## Switch to region and compare 
df0 <- SetAllIdent(df0, id = "region")

all.cell.region.dex <- list()
regions <- c("dMGE", "vMGE", "CGE")
for( i in regions ){
  temp <- list()
  for( ii in regions ){
    if( ii != i){
      temp[[length(temp) + 1]] <- rownames(FindMarkers(df0, i, ii, min.pct = 0.25,only.pos = T))
    }
  }
  all.cell.region.dex[[paste(i, "vs", ii)]] <- intersect(temp[[1]], temp[[2]])
}

all.cell.region.dex <- lapply(regions, function(i){
  lapply(regions, function(ii){
    temp <- list()
    if(ii != i){
      rownames(FindMarkers(df0, i, ii, min.pct = 0.25,only.pos = T))
    }
  })
})

ident.1 = "dMGE"
ident.2 = "vMGE"
temp <- FindMarkers(df0, ident.1, ident.2, min.pct = 0.25, only.pos = T)

## Switch to time point and subset data ##
df0 <- SetAllIdent(df0, id = "time_point")

## Analyze E12.5 
df0.E12.5 <- SubsetData(df0, ident.use = "E12.5")
df0.E12.5 <- FindVariableGenes(object = df0.E12.5, mean.function = ExpMean,
                               dispersion.function = LogVMR, x.low.cutoff = 0.0125,
                               x.high.cutoff = 3, y.cutoff = 0.5, do.plot = F)
df0.E12.5 <- RunPCA(object = df0.E12.5, pc.genes=df0.E12.5@var.genes,
              pcs.compute = 20, pcs.print = 1:12, genes.print = 10)

p1 <- PCAPlot(df0.E12.5, group.by = "region",do.return = T)
p2 <- PCAPlot(df0.E12.5,do.return = T)
p3 <- PCAPlot(df0.E12.5, group.by = "time_point",do.return = T)
pdf(paste(Sys.Date(), "E12.5_PCA.pdf", sep = "_"), width = 22, height = 20, useDingbats = F);
plot_grid(p1, p2, p3, ncol = 2)
dev.off()

df0.E12.5 <-JackStraw(df0.E12.5, num.replicate = 100, do.print = F, num.pc = 20)

pdf(file = paste(Sys.Date(), "E12.5_JackStraw.pdf", sep = "_"), width = 11,
    height = 10, useDingbats = F);
JackStrawPlot(df0.E12.5, PCs = 1:20)
dev.off()

# df0.E12.5 <- SetAllIdent(df0.E12.5, id = "res.1")
# TSNEPlot(df0.E12.5)

df0.E12.5 <- RunTSNE(df0.E12.5, dims.use = 1:8,#pcs[[i]], 
               perplexity = 20,
               seed.use = 3)

# df@ident <- as.factor(df@meta.data$region)
# names(df@ident) <- df@cell.names

df0.E12.5 <- FindClusters(object = df0.E12.5, 
                    reduction.type = "tsne", 
                    dims.use = 1:2, 
                    resolution = 0.8,
                    k.scale = 25,
                    prune.SNN = 0,
                    plot.SNN = F, 
                    print.output = F, 
                    save.SNN = F,
                    algorithm = 2,
                    force.recalc = TRUE, 
                    random.seed = 1,
                    k.param = 30) 

p1 <- TSNEPlot(object = df0.E12.5, pt.size = 3, do.return = T)
# p2 <-TSNEPlot(object = df0.E12.5, group.by = "time_point", do.return = T,
#               pt.size = 3)
p3 <-TSNEPlot(object = df0.E12.5, group.by = "region", do.return = T,
              pt.size = 3)

pdf(file = paste(Sys.Date(), "E12.5_tSNE.pdf", sep = "_"), width = 22, 
    height = 10, useDingbats = F);
plot_grid(p1, #p2, 
          p3, ncol = 2)
dev.off()

df0.E12.5.markers <- FindAllMarkers(genes.use = df0.E12.5@var.genes, 
                              object = df0.E12.5, only.pos = TRUE, 
                              min.pct = 0.25, 
                              return.thresh = 0.01)

top10 <- df0.E12.5.markers %>% group_by(cluster) %>% top_n(10, avg_diff)

pdf(file = paste(Sys.Date(), "E12.5_clusters_heatmap_top10.pdf", sep = "_"),
    width = 20, height = 20, useDingbats = F);
DoHeatmap(object = df0.E12.5,
          genes.use = top10$gene,#rownames(df.markers)[1:100],
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

rubenstein.markers <- read.csv("Rubenstein PGs vs Neuruon markers 25-10-2017.csv", header = F,
                               stringsAsFactors = F)
rubenstein.markers <- paste("\\b", rubenstein.markers[,1], "$", sep = "")
rubenstein.markers <- c(unlist(sapply(rubenstein.markers,
                                      function(x){grep(x, rownames(df0.E12.5@scale.data), 
                                                       value = T)})))

pdf(file = paste(Sys.Date(),
                 "E12.5_Rubenstein_genes.pdf", sep = "_"),
    width = 20, height = 30, useDingbats = F);
FeaturePlot(df0.E12.5, features.plot = rubenstein.markers, 
            nCol = 6,
            cols.use = c("lightgrey","blue"), pt.size = 1)
dev.off()

#### Random Forest prediction 
df0.E12.5 <- SetAllIdent(df0.E12.5, id = "cc.state")
df0.E12.5.n <- SubsetData(df0.E12.5, ident.use = "postmitotic")
df0.E12.5.p <- SubsetData(df0.E12.5, ident.remove = "postmitotic")

df0.E12.5.n <- SetAllIdent(df0.E12.5.n,id = "cell.ident")
df0.E12.5.p <- SetAllIdent(df0.E12.5.p,id = "cell.ident")

df0.E12.5.n <- SubsetData(df0.E12.5.n, ident.remove = "ZUncertain_E12")

#### Combine neuron with P56 cells
counts.allen <- read.csv("./raw_data/GSE71585_RefSeq_counts.csv", header = T, stringsAsFactors = F)
rownames(counts.allen) <- toupper(counts.allen[,1])
counts.allen <- counts.allen[,-1]

clu.allen <- read.csv("./raw_data/GSE71585_Clustering_Results.csv", header = T, stringsAsFactors = F)
keep <- clu.allen$broad_type == "GABA-ergic Neuron"

df.allen <- CreateSeuratObject(raw.data = as.vector(counts.allen[,keep]), min.cells = 5, 
                               min.genes = 1000, project = "Da_Mi")
df.allen <- NormalizeData(object = df.allen, normalization.method = "LogNormalize", 
                          scale.factor = 10000)
df.allen <- AddMetaData(object = df.allen, metadata = percent.mito, col.name = "percent.mito")
df.allen@meta.data <- cbind(df.allen@meta.data, clu.allen[keep,])
#df.allen <- ScaleData(df.allen)
df.allen <- ScaleData(df.allen)
df.allen <- FindVariableGenes(object = df.allen, do.plot = F)

percent.mito <- colSums(df.allen@raw.data[mito.genes, ])/colSums(df.allen@raw.data)
n.markers <- FindAllMarkers(df0.E12.5.n, only.pos = TRUE, min.pct = 0.25, 
                            return.thresh = 0.01)

n.top100 <- n.markers %>% group_by(cluster) %>% top_n(100, avg_diff)

DoHeatmap(object = df0.E12.5.n,
          genes.use = temp.top200$gene,#rownames(df.markers)[1:100],
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)

load("E12_tsne.Robject")
rownames(E12_tsne) <- rownames(E12_cca19)
E12_tsne <- data.frame(E12_tsne)
cate_bf_fs <- as.factor(df0.E12.5.n@meta.data$cell.ident)
feature_bf_fs <- E12_tsne[df0.E12.5.n@cell.names,]

rf_bf_fs <- randomForest(feature_bf_fs, cate_bf_fs, importance = TRUE, proximity = TRUE)
#imp_bf_fs <- importance(rf_bf_fs, type = 1)

newdata <- E12_tsne[intersect(rownames(E12_cca19), df0.E12.5.p@cell.names),]
pred_whole <- predict(rf_bf_fs, newdata = newdata)
pred_whole_prob <- predict(rf_bf_fs, newdata = newdata, type = "prob")

rec <- scale(pred_whole_prob,center = T,scale = T)
rec <- rec == apply(rec, 1, max) & rec > 1
rec <- unlist(apply(rec, 1, function(xx){colnames(rec)[xx]}))

df0.E12.5.p@meta.data$cell.ident <- "ZUncertain_E12"
ord <- match(names(rec), df0.E12.5.p@meta.data$New_Cell_ID)
df0.E12.5.p@meta.data$cell.ident[ord] <- as.character(rec)

df0.E12.5.p <- SetAllIdent(df0.E12.5.p, id = "cell.ident")

temp.markers <- FindAllMarkers(df0.E12.5.p, 
                               only.pos = TRUE, 
                               min.pct = 0, 
                               return.thresh = 1)

temp.top10 <- temp.markers %>% group_by(cluster) %>% top_n(10, avg_diff)

DoHeatmap(object = df0.E12.5,
          genes.use = c(temp.top10$gene,"ENSMUSG00000074622|MAFB"),#rownames(df.markers)[1:100],
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)

ord <- match(df0.E12.5.p@meta.data$New_Cell_ID, df0.E12.5@meta.data$New_Cell_ID)
df0.E12.5@meta.data$cell.ident[ord] <- df0.E12.5.p@meta.data$cell.ident
df0.E12.5@meta.data$cell.ident2 <- substr(df0.E12.5@meta.data$cell.ident, 1,3)

df0.E12.5 <- SetAllIdent(df0.E12.5, id = "cell.ident")

ddr.E12.5 <- SubsetData(df0.E12.5, ident.remove = "ZUncertain_E12")

temp.markers <- FindAllMarkers(ddr.E12.5, 
                               only.pos = TRUE, 
                               min.pct = 0.2, 
                               return.thresh = 0.001)

temp.top10 <- temp.markers %>% group_by(cluster) %>% top_n(10, avg_diff)

SstCbln4.markers <- FindMarkers(ddr.E12.5,ident.1 = "SstCbln4_E12",thresh.use = 0,
                                min.pct = 0.2)

DoHeatmap(object = ddr.E12.5,use.scaled = F,
          genes.use = c(temp.top10$gene,"ENSMUSG00000074622|MAFB"),#rownames(df.markers)[1:100],
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)

ddr.E12.5 <- SetAllIdent(ddr.E12.5, id = "New_Cell_ID")
keep <- intersect(rownames(E12_cca19), ddr.E12.5@cell.names)
temp.ddr.E12.5 <- SubsetData(ddr.E12.5, ident.use = keep)


ddr.E12.5 <- RunPCA(ddr.E12.5, pc.genes = temp.markers$gene, pcs.print = 1:5,
                    genes.print = 10)

p1 <- PCAPlot(ddr.E12.5, group.by = "region",do.return = T)
p2 <- PCAPlot(ddr.E12.5,do.return = T)
p3 <- PCAPlot(ddr.E12.5, group.by = "time_point",do.return = T)
pdf(paste(Sys.Date(), "E12.5_ddr_PCA.pdf", sep = "_"), width = 22, height = 20, 
    useDingbats = F)
plot_grid(p1, p2, p3, ncol = 2)
dev.off()

ddr.E12.5 <-JackStraw(ddr.E12.5, num.replicate = 100, do.print = F, num.pc = 10)

pdf(file = paste(Sys.Date(), "E12.5_ddr_JackStraw.pdf", sep = "_"), width = 11,
    height = 10, useDingbats = F);
JackStrawPlot(ddr.E12.5, PCs = 1:10)
dev.off()

ddr.E12.5 <- RunTSNE(ddr.E12.5, dims.use = 1:6,#pcs[[i]], 
                     perplexity = 20,
                     seed.use = 3)

p1 <- TSNEPlot(object = ddr.E12.5, pt.size = 3, do.return = T)
# p2 <-TSNEPlot(object = df0.E12.5, group.by = "time_point", do.return = T,
#               pt.size = 3)
p3 <-TSNEPlot(object = ddr.E12.5, group.by = "region", do.return = T,
              pt.size = 3)

pdf(file = paste(Sys.Date(), "E12.5_ddr_tSNE.pdf", sep = "_"), width = 22, 
    height = 10, useDingbats = F);
plot_grid(p1, #p2, 
          p3, ncol = 2)
dev.off()

expr <- t(E12_cca19[keep,])

temp <- read.csv("E12_assign_cell_identy_vip3groups_sst3groups_CCA19.csv", header = T, stringsAsFactors = F)
E12_tsne$x <- temp$x
E12_tsne <- merge(E12_tsne, df0.E12.5@meta.data, by = "row.names")

####
library(Rtsne)
tsne <- Rtsne(t(expr), pca = F)
tsne <- data.frame(tsne$Y)
rownames(tsne) <- colnames(expr)
tsne <- merge(tsne, pData(mm), by = "row.names")

ggplot(expr, aes(x = tSNE_1, y = tSNE_2, color = cell.ident)) +
  geom_point() +
  facet_wrap(~cell.ident)

expr <- E12_tsne[!(E12_tsne$cell.ident %in% "ZUncertain_E12"),]

ddrtree.results <- function(object, vars.to.use){
  pd <- object@meta.data
  rownames(pd) <- pd$New_Cell_ID
  pd <- new("AnnotatedDataFrame", data = pd)
  fd <- data.frame(gene_ID = rownames(object@data), 
                   gene_short_name = rownames(object@data))
  rownames(fd) <- fd$gene_ID
  fd <- new("AnnotatedDataFrame", data = fd)
  mm <- newCellDataSet(as(as.matrix(object@data), "sparseMatrix"),
                       phenoData = pd, 
                       featureData = fd,
                       lowerDetectionLimit = min(object@data), 
                       expressionFamily = gaussianff())
  
  ddrtree_res <- DDRTree(as.matrix(exprs(mm)[vars.to.use,]), dimensions = 2)
  mm@reducedDimS <- ddrtree_res$Z
  mm@reducedDimK <- ddrtree_res$Y
  colnames(mm@reducedDimS) <- colnames(mm@reducedDimK)<- pd$New_Cell_ID
  adjusted_K <- t(reducedDimK(mm))
  dp <- as.matrix(dist(adjusted_K))
  mm@cellPairwiseDistances <- dp
  gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
  dp_mst <- minimum.spanning.tree(gp)
  mm@minSpanningTree <- dp_mst
  mm@dim_reduce_type <- "DDRTree"
  
  S_matrix <- reducedDimS(mm)
  data_df <- data.frame(t(S_matrix[1:2,]))
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  rownames(data_df) <- pd$New_Cell_ID#res@cell.names#df@cell.names
  data_df <- merge(data_df, pData(mm), by.x = "row.names", by.y = "New_Cell_ID")
  
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
  
  return(list(data_df, edge_df, mm))
}

temp <- df0.E12.5@data[genes,]
data_df$MAF <- temp[2, as.character(data_df$New_Cell_ID)]
data_df$MAFB <- temp[1, as.character(data_df$New_Cell_ID)]

g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2)) + 
  geom_point(aes(color = cell.ident), na.rm = TRUE, size = 1, alpha = 0.8)
g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                 y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                 yend = "target_prin_graph_dim_2"), size = 0.3, linetype = "solid", 
                      na.rm = TRUE, data = edge_df)
g <- g + theme_bw()
# g <- g + facet_wrap(~res.1, nrow = 3)
pdf(file = paste(Sys.Date(), #samples[[i]][2], 
                 #samples[[i]][1],
                 "E12.5_id_cells.pdf", sep = "_"),
    width = 6, height = 3, useDingbats = F);
g
dev.off()

g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2)) + 
  geom_point(aes(color = cell.ident, alpha = cc.state, size = MAF), na.rm = TRUE)
g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                 y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                 yend = "target_prin_graph_dim_2"), size = 0.3, linetype = "solid", 
                      na.rm = TRUE, data = edge_df)
g <- g + facet_wrap(~cell.ident, nrow = 3)
g <- g + theme_bw()
g <- g + theme(strip.text = element_text(size = 40))
g <- g + scale_size_continuous(range = c(0, 20))
pdf(file = paste(Sys.Date(), #samples[[i]][2], 
                 #samples[[i]][1],
                 "E12.5_id_cells_facets_MAF.pdf", sep = "_"),
    width = 30, height = 30, useDingbats = F);
g
dev.off()

#### Do Pseudotime analysis for each individual lineage ####
for(i in unique(setdiff(df0.E12.5@ident, c("ZUncertain_E12")))){
  df0.temp <- SubsetData(df0.E12.5, ident.use = i)
  df0.temp <- SetAllIdent(df0.temp, id = "cc.state")
  
  markers <- FindMarkers(df0.temp, ident.1 = "postmitotic",min.cells = 0, random.seed = 1)
  
  results <- ddrtree.results(object = df0.temp, 
                             vars.to.use = rownames(markers)[markers$p_val < 0.01])
  
  data_df <- results[[1]]
  edge_df <- results[[2]]
  mm <- results[[3]]
  
  cc_ordering <- extract_ddrtree_ordering(mm, root_cell = NULL)
  pData(mm)$Pseudotime <- cc_ordering[data_df$Row.names,]$pseudo_time
  pData(mm)$State <- cc_ordering[data_df$Row.names,]$cell_state
  data_df$Pseudotime <-  cc_ordering[data_df$Row.names,]$pseudo_time
  data_df$State <- cc_ordering[data_df$Row.names,]$cell_state
  
  ## Optional: if Pseudotime is reversed, one can choose to correct it
  if(i %in% c("SstCbln4_E12","VipGpc3_E12", "SstChodl_E12", "Pvalb_E12")){
  root_cell <- select_root_cell(mm, reverse = T)
  mm@auxOrderingData[[mm@dim_reduce_type]]$root_cell <- root_cell
  
  cc_ordering_new_pseudotime <- extract_ddrtree_ordering(mm, root_cell) #re-calculate the pseudotime again
  pData(mm)$Pseudotime <- cc_ordering_new_pseudotime[data_df$Row.names,]$pseudo_time
  pData(mm)$State <- cc_ordering_new_pseudotime[data_df$Row.names,]$cell_state
  data_df$Pseudotime <-  cc_ordering_new_pseudotime[data_df$Row.names,]$pseudo_time
  data_df$State <- cc_ordering_new_pseudotime[data_df$Row.names,]$cell_state
  }
  
  # cor.res <- sort(apply(as.matrix(df0.temp@data)[,data_df$Row.names],
  #                       1, function(xx){cor(xx, data_df$Pseudotime, method = "p")}), decreasing = T)
  # 
  # genes <- c(names(cor.res)[1:98], "ENSMUSG00000074622|MAFB","ENSMUSG00000055435|MAF")
  # genes.df <- as.matrix(df0.temp@data[genes,])
  # data_df2 <- merge(data_df, t(genes.df), by.x = "Row.names", by.y = "row.names")
  # data_df2 <- melt(data_df2, measure.vars = genes)
  # 
  # g <- ggplot(data_df2, aes(x = Pseudotime, y = value)) +
  #   geom_point(size  = 5) +
  #   facet_wrap(~variable, nrow = 10, scales = "free") +
  #   geom_smooth() +
  #   theme(strip.text = element_text(size = 40))
  # 
  # pdf(file = paste(Sys.Date(), i, "E12.5_pseudotime_correlated_genes_facets.pdf", sep = "_"),
  #     width = 100, height = 100, useDingbats = F);
  # print(g)
  # dev.off()
  
  # genes <- c(names(cor.res)[(length(cor.res) - 99):length(cor.res)])
  # genes.df <- as.matrix(df0.temp@data[genes,])
  # data_df2 <- merge(data_df, t(genes.df), by.x = "Row.names", by.y = "row.names")
  # data_df2 <- melt(data_df2, measure.vars = genes)
  # 
  # g <- ggplot(data_df2, aes(x = Pseudotime, y = value)) +
  #   geom_point(size  = 5) +
  #   facet_wrap(~variable, nrow = 10, scales = "free") +
  #   geom_smooth() +
  #   theme(strip.text = element_text(size = 40))
  # 
  # pdf(file = paste(Sys.Date(), i, "E12.5_pseudotime_anti-correlated_genes_facets.pdf", sep = "_"),
  #     width = 100, height = 100, useDingbats = F);
  # print(g)
  # dev.off()
  # try({genes <- n.markers$gene[n.markers$cluster %in% i]
  # genes.df <- as.matrix(df0.temp@data[genes,])
  # data_df2 <- merge(data_df, t(genes.df), by.x = "Row.names", by.y = "row.names");
  # data_df2 <- melt(data_df2, measure.vars = genes);
  # 
  # g <- ggplot(data_df2, aes(x = Pseudotime, y = value)) +
  #   geom_point(size  = 5) +
  #   facet_wrap(~variable, nrow = 10, scales = "free") +
  #   geom_smooth() +
  #   theme(strip.text = element_text(size = 40))
  # 
  # pdf(file = paste(Sys.Date(), i, "E12.5_pseudotime_neuron_marker_genes_facets.pdf", sep = "_"),
  #     width = 100, height = 100, useDingbats = F)
  # print(g)
  # dev.off()});
  
  g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2)) + 
    geom_point(aes(color = cc.state), na.rm = TRUE, size = 1, alpha = 0.8)
  g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                   y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                   yend = "target_prin_graph_dim_2"), size = 0.3, linetype = "solid", 
                        na.rm = TRUE, data = edge_df)
  g <- g + theme_bw()
  
  g2 <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2)) + 
    geom_point(aes(color = Pseudotime), na.rm = TRUE, size = 1, alpha = 0.8)
  g2 <- g2 + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                   y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                   yend = "target_prin_graph_dim_2"), size = 0.3, linetype = "solid", 
                        na.rm = TRUE, data = edge_df)
  g2 <- g2 + theme_bw()
  
  # g <- g + facet_wrap(~res.1, nrow = 3)
  pdf(file = paste(Sys.Date(), i, "E12.5_pseudotime_trajectory.pdf", sep = "_"),
      width = 12, height = 5, useDingbats = F)
  print(plot_grid(g, g2,ncol = 2))
  dev.off();
}

g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2)) + 
  geom_point(aes(color = cc.state), na.rm = TRUE, size = 1, alpha = 0.8)
g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                 y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                 yend = "target_prin_graph_dim_2"), size = 0.3, linetype = "solid", 
                      na.rm = TRUE, data = edge_df)
g <- g + theme_bw()
# g <- g + facet_wrap(~res.1, nrow = 3)
pdf(file = paste(Sys.Date(), #samples[[i]][2], 
                 #samples[[i]][1],
                 "E12.5_id_cells.pdf", sep = "_"),
    width = 6, height = 3, useDingbats = F);
g
dev.off()

g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2)) + 
  geom_point(aes(color = cell.ident, alpha = cc.state, size = MAF), na.rm = TRUE)
g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                 y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                 yend = "target_prin_graph_dim_2"), size = 0.3, linetype = "solid", 
                      na.rm = TRUE, data = edge_df)
g <- g + facet_wrap(~cell.ident, nrow = 3)
g <- g + theme_bw()
g <- g + theme(strip.text = element_text(size = 40))
g <- g + scale_size_continuous(range = c(0, 20))
pdf(file = paste(Sys.Date(), #samples[[i]][2], 
                 #samples[[i]][1],
                 "E12.5_id_cells_facets_MAF.pdf", sep = "_"),
    width = 30, height = 30, useDingbats = F);
g
dev.off()



## Analyze E14.5 
df0.E14.5 <- SubsetData(df0, ident.use = "E14.5")
df0.E14.5 <- FindVariableGenes(object = df0.E14.5, mean.function = ExpMean,
                               dispersion.function = LogVMR, x.low.cutoff = 0.0125,
                               x.high.cutoff = 3, y.cutoff = 0.5, do.plot = F)
df0.E14.5 <- RunPCA(object = df0.E14.5, pc.genes=df0.E14.5@var.genes,
                    pcs.compute = 20, pcs.print = 1:12, genes.print = 10)
p1 <- PCAPlot(df0.E14.5, group.by = "region",do.return = T)
p2 <- PCAPlot(df0.E14.5,do.return = T)
p3 <- PCAPlot(df0.E14.5, group.by = "time_point",do.return = T)
pdf(paste(Sys.Date(), "E14.5_PCA.pdf", sep = "_"), width = 22, height = 20, useDingbats = F);
plot_grid(p1, p2, p3, ncol = 2)
dev.off()

df0.E14.5 <-JackStraw(df0.E14.5, num.replicate = 100, do.print = F,num.pc = 20)

pdf(file = paste(Sys.Date(), "E14.5_JackStraw.pdf", sep = "_"), width = 11,
    height = 10, useDingbats = F);
JackStrawPlot(df0.E14.5, PCs = 1:20)
dev.off()

# df0.E14.5 <- SetAllIdent(df0.E14.5, id = "res.1")
# TSNEPlot(df0.E14.5)

df0.E14.5 <- RunTSNE(df0.E14.5, dims.use = 1:8,#pcs[[i]], 
                     perplexity = 20,
                     seed.use = 3)

# df@ident <- as.factor(df@meta.data$region)
# names(df@ident) <- df@cell.names

df0.E14.5 <- FindClusters(object = df0.E14.5, 
                          reduction.type = "tsne", 
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

p1 <- TSNEPlot(object = df0.E14.5, pt.size = 3, do.return = T)
p2 <-TSNEPlot(object = df0.E14.5, group.by = "region", do.return = T,
              pt.size = 3)

pdf(file = paste(Sys.Date(), "E14.5_tSNE.pdf", sep = "_"), width = 22, 
    height = 10, useDingbats = F);
plot_grid(p1, p2, ncol = 3)
dev.off()

df0.E14.5.markers <- FindAllMarkers(genes.use = df0.E14.5@var.genes, 
                                    object = df0.E14.5, only.pos = TRUE, 
                                    min.pct = 0.25, 
                                    return.thresh = 0.01)

top10 <- df0.E14.5.markers %>% group_by(cluster) %>% top_n(10, avg_diff)

pdf(file = paste(Sys.Date(), "E14.5_clusters_heatmap_top10.pdf", sep = "_"),
    width = 20, height = 20, useDingbats = F);
DoHeatmap(object = df0.E14.5,
          genes.use = top10$gene,#rownames(df.markers)[1:100],
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

rubenstein.markers <- read.csv("Rubenstein PGs vs Neuruon markers 25-10-2017.csv", header = F,
                               stringsAsFactors = F)
rubenstein.markers <- paste("\\b", rubenstein.markers[,1], "$", sep = "")
rubenstein.markers <- c(unlist(sapply(rubenstein.markers,
                                      function(x){grep(x, rownames(df0.E14.5@scale.data), 
                                                       value = T)})))
pdf(file = paste(Sys.Date(),
                 "E14.5_Rubenstein_genes.pdf", sep = "_"),
    width = 20, height = 30, useDingbats = F);
FeaturePlot(df0.E14.5, features.plot = rubenstein.markers, 
            nCol = 6,
            cols.use = c("lightgrey","blue"), pt.size = 1)
dev.off()

####
var.intersect <- intersect( df.E12.5@var.genes, df.E14.5@var.genes )

temp <- sapply( unique( df.E12.5@meta.data$res.1 ), function( xx ){
  temp <- df.E12.5@scale.data[ df0@var.genes, df.E12.5@meta.data$res.1 == xx ]
  prcomp(temp)$x[ ,1 ]
})
temp <- temp[ var.intersect, ]

temp.2 <- sapply( unique( df.E14.5@meta.data$res.1 ), function( xx ){
  temp <- df.E14.5@scale.data[ df0@var.genes, df.E14.5@meta.data$res.1 == xx ]
  prcomp(temp)$x[ ,1 ]
})
temp.2 <- temp.2[ var.intersect, ]

res <- apply( temp, 2, function( xx ){apply( temp.2, 2, function( yy ){cor( xx, yy, method = "s" )} )} )
res[!(res == rowMax(res) | t(t(res) == rowMax(t(res))))] <- -1
heatmap.3( res, Rowv = T, Colv = T )

#### Monocle ####
library(monocle)
source("/home/zl242/Scripts/R/order_cells.R")

#### Analyze identified cells only ####
cell.ident <- read.csv("E12.5_cell_identy_vip3groups_sst3groups_CCA19.csv", header = F, stringsAsFactors = F)
cell.ident[,4] <- paste("C1", cell.ident[,1], cell.ident[,2], sep = ".")
colnames(cell.ident) <- c("time_point", "IFC", "cell.ident", "New_Cell_ID")
cell.ident <- cell.ident[ cell.ident$New_Cell_ID %in% df0.E12.5@cell.names, ]

df0.E12.5@meta.data$cell.ident <- "ZUncertain_E12"

df0.E12.5@meta.data$cell.ident[na.exclude(match(cell.ident$New_Cell_ID, df0.E12.5@cell.names))] <- cell.ident$cell.ident

df0.E12.5 <- SetAllIdent(df0.E12.5, id = "cell.ident")
res <- SubsetData(df0.E12.5,ident.remove = "ZUncertain_E12")

res <- FindVariableGenes(object = res, mean.function = ExpMean,
                               dispersion.function = LogVMR, x.low.cutoff = 0.0125,
                               x.high.cutoff = 3, y.cutoff = 0.5, do.plot = F)

res <- SetAllIdent(res, id="cell.ident")

res.markers <- FindAllMarkers(object = res, only.pos = TRUE,
                              min.pct = 0.25,
                              return.thresh = 0.01)

res.markers <- FindMarkers(res, ident.1 = "postmitotic",
                           min.pct = 0.25)

res.top100 <- res.markers %>% group_by(cluster) %>% top_n(100, avg_diff)

pdf(file = paste(Sys.Date(), "id_cells_heatmap_top10.pdf", sep = "_"),
    width = 20, height = 20, useDingbats = F);
DoHeatmap(object = res,
          genes.use = res.top10$gene,#rownames(df.markers)[1:100],
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

res <- RunPCA(object = res, pc.genes=res.top100$gene,
                    pcs.compute = 20, pcs.print = 1:12, genes.print = 10)
p1 <- PCAPlot(res, group.by = "cell.ident",do.return = T)
#p2 <- PCAPlot(res,do.return = T)
p3 <- PCAPlot(res, group.by = "region",do.return = T)
pdf(paste(Sys.Date(), "E12.5_id_cells_PCA.pdf", sep = "_"), width = 22, height = 20, useDingbats = F);
plot_grid(p1, #p2, 
          p3, ncol = 2)
dev.off()

res <-JackStraw(res, num.replicate = 100, do.print = F,num.pc = 20)

pdf(file = paste(Sys.Date(), "E12.5_id_cells_JackStraw.pdf", sep = "_"), width = 11,
    height = 10, useDingbats = F);
JackStrawPlot(res, PCs = 1:20)
dev.off()

load("E12_cca19.Robject")

new.names <- sapply(rownames(E12_cca19)[735:1798], function(xx){
  res <- strsplit(xx, split = "_")[[1]]
  paste("C1", res[3], res[4], sep = ".")},USE.NAMES = F)
rownames(E12_cca19)[735:1798] <- new.names
pd <- read.csv("E12_assign_cell_identy_vip3groups_sst3groups_CCA19.csv", header = T, stringsAsFactors = F)
pd[735:1798, 1] <- new.names
pd$source <- "allen"
pd$source[735:1798] <- "Da_Mi"

rec <- apply( df0@meta.data[,4:ncol(df0@meta.data)], 2, function( xx ){apply( df0@dr$pca@cell.embeddings, 2, function( yy ){cor( as.numeric(as.factor(xx)), yy, method = "s" )} )} )
rec[abs(rec) < 0.3] <- NA
heatmap.3(rec, Colv = F, Rowv = F,na.color = "black")

# ddrtree_res <- DDRTree(t(res@dr$pca@cell.embeddings[,1:5]), dimensions = 2)
rownames(pd) <- rownames(E12_cca19)
pd <- new("AnnotatedDataFrame", data = pd)
fd <- data.frame(gene_ID = colnames(E12_cca19), 
                 gene_short_name = colnames(E12_cca19))
rownames(fd) <- fd$gene_ID
fd <- new("AnnotatedDataFrame", data = fd)
mm <- newCellDataSet(as(as.matrix(t(E12_cca19)), "sparseMatrix"),
                     phenoData = pd, 
                     featureData = fd,
                     lowerDetectionLimit = min(E12_cca19), 
                     expressionFamily = gaussianff())

library(Rtsne)
seed = 3
set.seed(seed)
p = 30
cca.tSNE <- Rtsne(E12_cca19, pca = F, perplexity = p)
to.plot <- data.frame(cca.tSNE$Y)
colnames(to.plot) <- c("X1", "X2")
rownames(to.plot) <- rownames(E12_cca19)
to.plot <- cbind(to.plot, pd)
to.plot <- merge(x = to.plot, y = df0.E12.5@meta.data, by.x = "row.names", by.y = "New_Cell_ID", all.x = T)

pdf(file = paste(Sys.Date(), "_E12.5_cca19_tSNE.pdf", sep = ""), 
    width = 15, height = 10, useDingbats = F)
ggplot(data = subset(to.plot, source == "Da_Mi" & x != "ZUncertain_E12"), aes(x = X1, y = X2)) +
  #geom_point()
  #geom_point(data = subset(to.plot, source == "Da_Mi"), size = 3, aes(x = X1, y = X2, color = x))# +
  geom_point(aes(x = X1, y = X2, color = x), size = 3)  +
  scale_color_discrete(name = "Time point") +
  labs(x = "tSNE 1", y = "tSNE 2") +
  ggtitle(paste("E12.5 assigned cells, Perplexity =", p)) +
  theme_bw() +
  theme(text = element_text(size = 20),
        line = element_line(size = 1),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.key = element_blank(),
        #legend.position = "",#c(0.9, 0.1),
        legend.background = element_blank())
dev.off(dev.cur())

ddrtree_res <- DDRTree(t(cca.tSNE$Y), dimensions = 2)
mm@reducedDimS <- ddrtree_res$Z
mm@reducedDimK <- ddrtree_res$Y
adjusted_K <- t(reducedDimK(mm))
dp <- as.matrix(dist(adjusted_K))
mm@cellPairwiseDistances <- dp
gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
dp_mst <- minimum.spanning.tree(gp)
mm@minSpanningTree <- dp_mst
mm@dim_reduce_type <- "DDRTree"

cc_ordering <- extract_ddrtree_ordering(mm, root_cell = NULL)
pData(mm)$Pseudotime <- cc_ordering$pseudo_time
# pData(mm)$Pseudodev <- cc_ordering$pseudo_time

# ggplot(pData(mm), aes(x = Pseudotime, y = Pseudodev, color = res@meta.data$time_point)) + 
#   geom_point() +
#   scale_color_discrete(name = "")

# plot_cell_trajectory(mm, color_by = "res.1")

# mm <- setOrderingFilter(mm, ordering_genes = df@var.genes)

S_matrix <- reducedDimS(mm)
data_df <- data.frame(t(S_matrix[1:2,]))
colnames(data_df) <- c("data_dim_1", "data_dim_2")
rownames(data_df) <- rownames(E12_cca19)#res@cell.names#df@cell.names
data_df <- merge(data_df, pData(mm), by.x = "row.names", by.y = "row.names")
data_df <- merge(x = data_df, y = df0.E12.5@meta.data, by.x = "Row.names", by.y = "New_Cell_ID", all.x = T)

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

df <- subset(data_df, source == "Da_Mi" & cc.state != "postmitotic")
g <- ggplot(data = df, aes(x = data_dim_1, y = data_dim_2)) + 
  geom_point(aes(color = x), na.rm = TRUE, size = 1, alpha = 0.8)
g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                 y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                 yend = "target_prin_graph_dim_2"), size = 0.3, linetype = "solid", 
                      na.rm = TRUE, data = edge_df)
g <- g + theme_bw()
# g <- g + facet_wrap(~res.1, nrow = 3)
pdf(file = paste(Sys.Date(), #samples[[i]][2], 
                 #samples[[i]][1],
                 "E12.5_id_cells.pdf", sep = "_"),
    width = 6, height = 3, useDingbats = F);
g
dev.off()

g <- ggplot(data = df, aes(x = data_dim_1, y = data_dim_2)) + 
  geom_point(aes(color = x), na.rm = TRUE, size = 10)
g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                 y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                 yend = "target_prin_graph_dim_2"), size = 0.3, linetype = "solid", 
                      na.rm = TRUE, data = edge_df)
g <- g + facet_wrap(~x, nrow = 3)
g <- g + theme_bw()
g <- g + theme(strip.text = element_text(size = 40))
pdf(file = paste(Sys.Date(), #samples[[i]][2], 
                 #samples[[i]][1],
                 "E12.5_id_cells_facets.pdf", sep = "_"),
    width = 30, height = 30, useDingbats = F);
g
dev.off()

rubenstein.markers <- read.csv("Rubenstein PGs vs Neuruon markers 25-10-2017.csv", header = F,
                               stringsAsFactors = F)
rubenstein.markers <- paste("\\b", rubenstein.markers[,1], "$", sep = "")
rubenstein.markers <- c(unlist(sapply(rubenstein.markers,
                                      function(x){grep(x, rownames(res@scale.data), value = T)})))

rec <- cbind(data_df, t(res@scale.data[c("ENSMUSG00000031285|DCX", "ENSMUSG00000019874|FABP7"),]))
rec <- melt(rec,id.vars = colnames(data_df))

g <- ggplot(data = rec, aes(x = data_dim_1, y = data_dim_2)) + 
  geom_point(aes(color = value), na.rm = TRUE)
g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                 y = "source_prin_graph_dim_2", 
                                 xend = "target_prin_graph_dim_1", 
                                 yend = "target_prin_graph_dim_2"), 
                      size = 0.3, linetype = "solid", 
                      na.rm = TRUE, data = edge_df)
g <- g + theme_bw()
g <- g + facet_wrap(~variable, nrow = 1)
pdf(file = paste(Sys.Date(), #samples[[i]][2], 
                 #samples[[i]][1],
                 "id_cells_Reubenstein_genes_facets.pdf", sep = "_"),
    width = 40, height = 20, useDingbats = F);
g
dev.off()
