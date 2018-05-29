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

keep <- !(GB$IFC %in% c("180","181","191","192","193","194","195", "196")) #&
  # !(GB$New_Cell_ID %in% c("C1.173.A8",  "C1.173.A9",  "C1.173.B12", "C1.173.C7",  "C1.184.A2",
  #                         "C1.184.A3",  "C1.184.A8", "C1.184.D1",  "C1.176.A8", "C1.90.H11", 
  #                         "C1.92.D4", "C1.92.A2", "C1.92.E6", "C1.97.A5","C1.184.A1",  "C1.184.E10",
  #                         "C1.101.A6", "C1.101.A7", "C1.101.A8","C1.98.C5","C1.98.C9",
  #                         "C1.171.A4", "C1.174.A1", "C1.177.E6","C1.177.A12", "C1.177.C12",
  #                         "C1.98.D4","C1.98.D5","C1.98.E11","C1.98.E7","C1.98.F10",
  #                         "C1.182.A3",  "C1.186.F3",  "C1.187.D10",
  #                         "C1.98.F7","C1.98.G11", "C1.91.D11", "C1.91.E11", "C1.91.E4",
  #                         "C1.91.H10", "C1.93.B12", "C1.93.C6","C1.93.C9", "C1.93.E10",
  #                         "C1.93.E8","C1.93.G5","C1.93.H1","C1.93.H5","C1.96.B6","C1.96.C6",
  #                         "C1.98.A3", "C1.98.A7", "C1.98.A8", "C1.93.A7", "C1.93.A6",
  #                         "C1.160.B4","C1.160.F12", "C1.161.B9", "C1.161.F12","C1.162.C6",
  #                         "C1.162.G6","C1.190.2.E4"))
df0 <- CreateSeuratObject(raw.data = as.vector(expr[,keep]), min.cells = 5, min.genes = 1000, project = "Da_Mi")
df0 <- NormalizeData(object = df0, normalization.method = "LogNormalize", scale.factor = 1000000)
mito.genes <- grep(pattern = "\\bMT-", x = rownames(x = df0@data), value = TRUE)
percent.mito <- colSums(df0@raw.data[mito.genes, ])/colSums(df0@raw.data)
# AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
df0 <- AddMetaData(object = df0, metadata = percent.mito, col.name = "percent.mito")
df0@meta.data <- cbind(df0@meta.data, GB[df0@cell.names,])
#df0 <- ScaleData(df0)
df0 <- ScaleData(df0, vars.to.regress = "nGene")

df0 <- FindVariableGenes(object = df0, mean.function = ExpMean, 
                         dispersion.function = LogVMR, x.low.cutoff = 1, 
                         x.high.cutoff = Inf, y.cutoff = 0.5, y.high.cutoff = 5, 
                         do.plot = F)

df0 <- RunPCA(object = df0, pc.genes=df0@var.genes, pcs.compute = 20, pcs.print = 1:12,
              genes.print = 10)

pdf(file = paste(Sys.Date(), "PCA_pairs.pdf", sep = "_"), 
    width = 40, height = 40, useDingbats = F);
pairs(df0@dr$pca@cell.embeddings[,1:10], col = as.factor(df0@meta.data$mouse), pch = 16, cex = 1)
dev.off()

rec <- apply( df0@meta.data[,1:ncol(df0@meta.data)], 2, 
              function( xx ){apply( df0@dr$pca@cell.embeddings, 2, function( yy ){cor( as.numeric(as.factor(xx)), yy, method = "s" )} )} )
rec[abs(rec) < 0.2] <- NA
heatmap.3(rec, Colv = F, Rowv = F,na.color = "black")

p1 <- PCAPlot(df0, group.by = "region",do.return = T)
p2 <- PCAPlot(df0, group.by = "time_point", do.return = T)
p3 <- PCAPlot(df0, group.by = "IFC",do.return = T)
pdf(paste(Sys.Date(), "all_cells_PCA.pdf", sep = "_"), width = 26, height = 8, useDingbats = F);
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

# df0@meta.data$PC.to.regress2 <- df0@dr$pca@cell.embeddings[,1]
# df0 <- ScaleData(df0, vars.to.regress = c("PC.to.regress", "PC.to.regress2"))
# df0 <- FindVariableGenes(object = df0, mean.function = ExpMean,
#                         dispersion.function = LogVMR, x.low.cutoff = 0.0125,
#                         x.high.cutoff = 3, y.cutoff = 0.5, do.plot = F)
# df0 <- RunPCA(object = df0, pc.genes=df0@var.genes, pcs.compute = 20, pcs.print = 1:12,genes.print = 10)

## Correlation between factors and PCs
rec <- apply( df0@meta.data[,4:ncol(df0@meta.data)], 2, function( xx ){apply( df0@dr$pca@cell.embeddings, 2, function( yy ){cor( as.numeric(as.factor(xx)), yy, method = "s" )} )} )
rec[abs(rec) < 0.25] <- NA
heatmap.3(rec, Colv = F, Rowv = F,na.color = "black")
pdf("All_cells_PCA_region_norm.pdf", width = 11, height = 10, useDingbats = F);
PCAPlot(df0, group.by = "region")
dev.off()

grep("\\bTPBG$", rownames(df0@scale.data),value = T)
CALB1 <- df0@scale.data["ENSMUSG00000028222|CALB1",]
TPBG <- df@scale.data["ENSMUSG00000035274|TPBG",]
CALB1.corr <- apply(df0@scale.data,1,function(xx){cor(xx, CALB1)})

#### 1. Analyze samples seperately based on age and region ####
#### 1.1 Analyze individual datasets ####
#### 1.1.1 Setup Seurat ####
dfs <- list()
dfs.markers <- list()
for(i in 1:length(samples)){
  keep <- GB$time_point %in% samples[[i]] &
  #keep <- GB$region %in% samples[[i]] & GB$time_point %in% samples[[i]] & 
    !(GB$IFC %in% c("180","181","191","192","193","194","195", "196")) #&
    # !(GB$New_Cell_ID %in% c("C1.173.A8",  "C1.173.A9",  "C1.173.B12", "C1.173.C7",  "C1.184.A2",
    #                         "C1.184.A3",  "C1.184.A8", "C1.184.D1",  "C1.176.A8", "C1.90.H11", 
    #                         "C1.92.D4", "C1.92.A2", "C1.92.E6", "C1.97.A5","C1.184.A1",  "C1.184.E10",
    #                       "C1.101.A6", "C1.101.A7", "C1.101.A8","C1.98.C5","C1.98.C9",
    #                       "C1.171.A4", "C1.174.A1", "C1.177.E6","C1.177.A12", "C1.177.C12",
    #                       "C1.98.D4","C1.98.D5","C1.98.E11","C1.98.E7","C1.98.F10",
    #                       "C1.182.A3",  "C1.186.F3",  "C1.187.D10",
    #                       "C1.98.F7","C1.98.G11", "C1.91.D11", "C1.91.E11", "C1.91.E4",
    #                       "C1.91.H10", "C1.93.B12", "C1.93.C6","C1.93.C9", "C1.93.E10",
    #                       "C1.93.E8","C1.93.G5","C1.93.H1","C1.93.H5","C1.96.B6","C1.96.C6",
    #                       "C1.98.A3", "C1.98.A7", "C1.98.A8", "C1.93.A7", "C1.93.A6",
    #                       "C1.160.B4","C1.160.F12", "C1.161.B9", "C1.161.F12","C1.162.C6",
    #                       "C1.162.G6","C1.190.2.E4"))
  df <- CreateSeuratObject(raw.data = as.vector(expr[,keep]), min.cells = 5, min.genes = 1000, project = "Da_Mi")
  df <- NormalizeData(object = df, normalization.method = "LogNormalize", scale.factor = 100000)
  mito.genes <- grep(pattern = "\\bMT-", x = rownames(x = df@data), value = TRUE)
  percent.mito <- colSums(df@raw.data[mito.genes, ])/colSums(df@raw.data)
  # AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
  df <- AddMetaData(object = df, metadata = percent.mito, col.name = "percent.mito")
  df@meta.data <- cbind(df@meta.data, GB[df@cell.names,])
    #write.csv(df@scale.data, file = "E14.5_expression.csv")
  
  #### For E14.5 only
  # df@meta.data$PC1 <- df@dr$pca@cell.embeddings[,1]
  
  cc_genes <- read.csv("cell_cycle_genes.csv", header = F, stringsAsFactors = F)
  res <- sapply(cc_genes[,1], function(xx){paste("\\b",xx,"$",sep = "")})
  res <- sapply(res, function(xx){grep(xx, rownames(df@data), value = T)})
  cc_genes[,3] <- as.character(res)
  g1s.genes <- cc_genes[,3][cc_genes$V2 == "G1/S"]
  g2m.genes <- cc_genes[,3][cc_genes$V2 == "G2/M"]
  df <- CellCycleScoring(df, s.genes = g1s.genes, g2m.genes = g2m.genes, set.ident = T)
  df@meta.data$CC.Difference <- df@meta.data$S.Score - df@meta.data$G2M.Score
  df@meta.data$cc.state <- as.character(df@ident)
  df@meta.data$cc.state[df@meta.data$cc.state == "G1"] <- "postmitotic"
  df@meta.data$cc.state[df@meta.data$cc.state == "S"] <- "G1S"
  
  df <- ScaleData(df, vars.to.regress = c("nGene","CC.Difference"))
  
  df <- FindVariableGenes(object = df, mean.function = ExpMean,
                          dispersion.function = LogVMR, x.low.cutoff = 1,
                          x.high.cutoff = Inf, y.cutoff = 0.5,y.high.cutoff = 5, 
                          do.plot = F)
  
  set.seed(1)
  df <- RunPCA(object = df, pc.genes=df@var.genes,
               pcs.compute = 20, pcs.print = 1:12, genes.print = 10)
  
  # df <- RunPCA(object = df, pc.genes = cc_genes[,3], genes.print = 10, do.print = F)
  # PCAPlot(df)
  
  #save(df, file = "2017-10-08_E14.5_normalized_Seurat_object.Rdata")
    # VlnPlot(object = df, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
  # GenePlot(object = df, gene1 = "nUMI", gene2 = "percent.mito")
  # GenePlot(object = df, gene1 = "nUMI", gene2 = "nGene")
  # 
  # pdf(file = paste("./E12.5_dMGE/", Sys.Date(), "_nGene_E12.5_dMGE.pdf", sep = ""), 
  #     width = 11, height = 10, useDingbats = F);
  # ggplot(df@meta.data) + geom_boxplot(aes(x = IFC, y = nUMI, color = IFC))
  # dev.off()
  
  pdf(file = paste(Sys.Date(), samples[[i]][2], samples[[i]][1], "PCA_pairs.pdf", sep = "_"), 
      width = 20, height = 20, useDingbats = F);
  pairs(df@dr$pca@cell.embeddings[,1:5], col = as.factor(df@meta.data$mouse), pch = 16, cex = 3)
  dev.off()
  
  rec <- apply( df@meta.data[,1:ncol(df@meta.data)], 2, 
                function( xx ){apply( df@dr$pca@cell.embeddings, 2, function( yy ){cor( as.numeric(as.factor(xx)), yy, method = "s" )} )} )
  rec[abs(rec) < 0.2] <- NA
  heatmap.3(rec, Colv = F, Rowv = F,na.color = "black")
  
  p1 <- PCAPlot(df, group.by = "region",do.return = T)
  p2 <- PCAPlot(df, group.by = "time_point", do.return = T)
  p3 <- PCAPlot(df, group.by = "IFC",do.return = T)
  pdf(paste(Sys.Date(), samples[[i]][2],"PCA.pdf", sep = "_"), width = 26, height = 8, useDingbats = F);
  plot_grid(p1, p2, p3, ncol = 3)
  dev.off()
  
  rec <- apply( df@scale.data, 1, 
                function( xx ){apply( df@dr$pca@cell.embeddings, 2, function( yy ){cor( as.numeric(as.factor(xx)), yy, method = "s" )} )} )
  # rec[abs(rec) < 0.2] <- NA
  
  rec.genes <- apply(rec, 1, function(xx){
    names(sort(xx, decreasing = T))[1:100]
  })
  
  ## Get go terms for PC correlated genes
  library("RDAVIDWebService") ##--Load RDAVIDWebService 
  user <- DAVIDWebService$new(email = "zhen.li.zl242@yale.edu", 
                              url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")     ##--Log on to DAVID with email
  
  #KEGG <-list() 
  BP <- list()
  for(ii in 1:ncol(rec.genes)){
    david <- rec.genes[,ii]
    david <- substr(david, 1, 18)
    
    ##--Get a list of genes, pay attention to idType
    ##--Read a gene list
    ##--submit gene list to DAVID
    addList(user, inputIds = david, idType = "ENSEMBL_GENE_ID", listType = "Gene", 
            listName = colnames(rec.genes)[ii])
    
    ##--check a specific category
    # setAnnotationCategories(user, "KEGG_PATHWAY")                                  
    # 
    # ##--combine GO terms from each list
    # term <- getFunctionalAnnotationChart(user)
    # KEGG[[i]] <- term$Term
    
    setAnnotationCategories(user, "GOTERM_BP_FAT")
    term <- getFunctionalAnnotationChart(user)
    BP[[ii]] <- term$Term
  }
  
  BP.temp <- sapply(BP, function(xx){xx[1:20]})
  colnames(BP.temp) <- colnames(rec.genes)
  write.csv(BP.temp, file = paste(Sys.Date(), samples[[i]][2],"BP_GO_terms.csv", sep = "_")) 
  
  ##genes <- c("ENSMUSG00000070880|GAD1", "ENSMUSG00000020914|TOP2A", "ENSMUSG00000019874|FABP7")
  #FeaturePlot(df,features.plot = genes, reduction.use = "pca")
  # pdf(file = paste( Sys.Date(),samples[[i]][2], samples[[i]][1], 
  #                   "PCA_norm.pdf", sep = "_"), 
  #     width = 11, height = 10, useDingbats = F);
  g <- PCAPlot(df, dim.1 = 1, dim.2 = 2, group.by = "mouse",do.return = T)
  # g <- g + geom_hline(yintercept = 15)
  # g <- g + geom_vline(xintercept = 15)
  # g <- g + geom_text(aes(label = df@cell.names), size = 3)
  g
  # dev.off()
  
  biocLite("biomaRt")
  
  #names(which(df@dr$pca@cell.embeddings[,2] > 15 | df@dr$pca@cell.embeddings[,2] < -20))
  ## Use Jack Straw to determine significant PC
  df0 <-JackStraw(df0, num.replicate = 100, do.print = F,num.pc = 20)

  pdf(file = paste(Sys.Date(), samples[[i]][2], #samples[[i]][1],
                   "JackStraw.pdf", sep = "_"), width = 11,
      height = 10, useDingbats = F);
  JackStrawPlot(df0, PCs = 1:20)
  dev.off()

  # PCElbowPlot(df, num.pc = 20)
  
  pcs <- list(c(1:3), c(1:7), c(1:6), c(1:2), c(1:2), c(1:2))
  k.param <- c(20, 20, 20, 20, 20, 20)
  # for(ii in 1:11){
  #   pdf(file = paste(Sys.Date(), samples[[i]][2], "PC",ii,"Heatmap.pdf", sep = "_"), 
  #       width = 11, height = 10, useDingbats = F);
  #   PCHeatmap(object = df, pc.use = ii, do.balanced = T,#TRUE, 
  #             label.columns = FALSE, remove.key = TRUE)
  #   dev.off()
  # }
  
  df <- RunTSNE(df, dims.use = c(1, 6, 7),#pcs[[i]], 
                perplexity = 20,
                seed.use = 3)
  
  # df@ident <- as.factor(df@meta.data$region)
  # names(df@ident) <- df@cell.names
  
  df <- FindClusters(object = df, 
                     reduction.type = "pca", 
                     dims.use = c(1, 6, 7), 
                     resolution = 1, # dMGE, CGE = 0.8, E14.5 vMGE = 2
                     k.scale = 25, # dMGE, CGE = 50
                     prune.SNN = 0, #dMGE, CGE = 0.5
                     plot.SNN = F, 
                     print.output = F, 
                     save.SNN = F,
                     algorithm = 2,
                     force.recalc = TRUE, 
                     random.seed = 1,
                     k.param = k.param[i]) 
  
  p1 <- TSNEPlot(object = df, pt.size = 3, do.return = T)
  p2 <-TSNEPlot(object = df, group.by = "time_point", do.return = T,
                pt.size = 3)
  p3 <-TSNEPlot(object = df, group.by = "region", do.return = T,
                pt.size = 3)
  
  pdf(file = paste(Sys.Date(), samples[[i]][2], #samples[[i]][1],  
                   "tSNE.pdf", sep = "_"), width = 33, 
      height = 10, useDingbats = F);
  plot_grid(p1, p2, p3, ncol = 3)
  dev.off()
  
  df.markers <- FindAllMarkers(genes.use = df@var.genes, 
                               object = df, only.pos = TRUE, 
                               min.pct = 0.25, 
                               return.thresh = 0.01)
  #df.markers <- FindMarkers(df, ident.1 = "vMGE", ident.2 = "dMGE",only.pos = T,min.pct = 0.5)
  
  top10 <- df.markers %>% group_by(cluster) %>% top_n(10, avg_diff)
  
  #write.csv(df.markers, file = paste(samples[[i]][2], samples[[i]][1], "_markers_clusters.csv"))
  
  ## Plot top 20 marker genes
  
  # glc <- read.csv("./raw_data/Gene list for clusters.csv", header = F)
  # glc$V1 <- toupper(glc$V1)
  # glc$V1 <- sapply(glc$V1, function(x){paste("\\b", x, "$", sep = "")})
  # glc$V1 <- sapply(glc$V1,function(x){grep(x, rownames(df@scale.data), value = T)})
  # glc$V1 <- as.character(glc$V1)
  # glc <- glc[glc$V1 != "character(0)",]
  
  grep("PROX1", rownames(df@scale.data), value = T)
  genes <- c("ENSMUSG00000069171|NR2F1", "ENSMUSG00000001496|NKX2-1", "ENSMUSG00000041309|NKX6-2",
    "ENSMUSG00000048402|GLI2", "ENSMUSG00000004151|ETV1", "ENSMUSG00000020950|FOXG1",
    "ENSMUSG00000028201|LHX8", "ENSMUSG00000026890|LHX6", "ENSMUSG00000027210|MEIS2",
    "ENSMUSG00000030551|NR2F2", "ENSMUSG00000048562|SP8", "ENSMUSG00000010175|PROX1")
  
  pdf(file = paste(Sys.Date(),
                   #samples[[i]][2], #samples[[i]][1],
                   "clusters_heatmap_top10.pdf", sep = "_"),
      width = 20, height = 20, useDingbats = F);
  DoHeatmap(object = df0,
            genes.use = top10$gene,#rownames(df.markers)[1:100],
            slim.col.label = TRUE,
            group.spacing = 0.3,
            remove.key = TRUE)
  dev.off()
  
  dfs[[paste(samples[[i]][2], samples[[i]][1], sep="_")]] <- df
  dfs.markers[[paste(samples[[i]][2], samples[[i]][1], sep="_")]] <- df.markers
}

#

df <- dfs[[6]]
temp <- df.E14.5@dr$tsne@cell.embeddings
write.csv(temp, file = "E14.5_tSNE_cell_embeddings.csv")

df.E12.5@meta.data$clusters <- "N8"
df.E12.5@meta.data$clusters[df.E12.5@meta.data$res.1 == "1" ] <- "N3"
df.E12.5@meta.data$clusters[df.E12.5@meta.data$res.1 == "2" ] <- "P5"
df.E12.5@meta.data$clusters[df.E12.5@meta.data$res.1 == "3" ] <- "N1"
df.E12.5@meta.data$clusters[df.E12.5@meta.data$res.1 == "4" ] <- "P4"
df.E12.5@meta.data$clusters[df.E12.5@meta.data$res.1 == "5" ] <- "P1"
df.E12.5@meta.data$clusters[df.E12.5@meta.data$res.1 == "6" ] <- "N4"
df.E12.5@meta.data$clusters[df.E12.5@meta.data$res.1 == "7" ] <- "P7"
df.E12.5@meta.data$clusters[df.E12.5@meta.data$res.1 == "8" ] <- "N5"
df.E12.5@meta.data$clusters[df.E12.5@meta.data$res.1 == "9" ] <- "N2"
df.E12.5@meta.data$clusters[df.E12.5@meta.data$res.1 == "10" ] <- "P3"
df.E12.5@meta.data$clusters[df.E12.5@meta.data$res.1 == "11" ] <- "N7"
df.E12.5@meta.data$clusters[df.E12.5@meta.data$res.1 == "12" ] <- "P2"
df.E12.5@meta.data$clusters[df.E12.5@meta.data$res.1 == "13" ] <- "N6"
df.E12.5@meta.data$clusters[df.E12.5@meta.data$res.1 == "14" ] <- "P6"
df.E12.5@meta.data$clusters <- factor(df.E12.5@meta.data$clusters,levels = c(
  paste("P", 1:7, sep = ""), paste("N", 1:8, sep = "")
))

df.E14.5@meta.data$clusters <- "N1"
df.E14.5@meta.data$clusters[df.E14.5@meta.data$res.1 == "1" ] <- "N5"
df.E14.5@meta.data$clusters[df.E14.5@meta.data$res.1 == "2" ] <- "P5"
df.E14.5@meta.data$clusters[df.E14.5@meta.data$res.1 == "3" ] <- "P3"
df.E14.5@meta.data$clusters[df.E14.5@meta.data$res.1 == "4" ] <- "N4"
df.E14.5@meta.data$clusters[df.E14.5@meta.data$res.1 == "5" ] <- "P6"
df.E14.5@meta.data$clusters[df.E14.5@meta.data$res.1 == "6" ] <- "P2"
df.E14.5@meta.data$clusters[df.E14.5@meta.data$res.1 == "7" ] <- "N2"
df.E14.5@meta.data$clusters[df.E14.5@meta.data$res.1 == "8" ] <- "P4"
df.E14.5@meta.data$clusters[df.E14.5@meta.data$res.1 == "9" ] <- "N3"
df.E14.5@meta.data$clusters[df.E14.5@meta.data$res.1 == "10" ] <- "P1"
df.E14.5@meta.data$clusters[df.E14.5@meta.data$res.1 == "11" ] <- "N6"
df.E14.5@meta.data$clusters[df.E14.5@meta.data$res.1 == "12" ] <- "P7"
df.E14.5@meta.data$clusters <- factor(df.E14.5@meta.data$clusters,levels = c(
  paste("P", 1:7, sep = ""), paste("N", 1:6, sep = "")
))

df.E12.5@meta.data <- df.E12.5@meta.data[, !(colnames(df.E12.5@meta.data) %in% 
                                           grep("cell.id", colnames(df.E12.5@meta.data), 
                                                value = T))]

cell.id <- read.csv("E12_assign_cell_identy_vip3groups_sst3groups_CCA19.csv", 
                      header = T, stringsAsFactors = F)
colnames(cell.id) <- c("New_Cell_ID", "cell.id")

new.names <- sapply(cell.id[735:1798,1], function(xx){
  res <- strsplit(xx, split = "\\.")[[1]]
  paste("C1", res[2], res[3], sep = ".")}, USE.NAMES = F)
new.names <- make.names(new.names,unique = T)
cell.id[735:1798,1] <- new.names

df.E12.5@meta.data <- merge(df.E12.5@meta.data, cell.id, by = "New_Cell_ID", all.x = T)
rownames(df.E12.5@meta.data) <- df.E12.5@meta.data$New_Cell_ID
df.E12.5@meta.data$cell.id[is.na(df.E12.5@meta.data$cell.id)] <- "ZUncertain_E12"
df.E12.5@meta.data$cell.id2 <- substr(df.E12.5@meta.data$cell.id, 1, 3)

#### Change cell.type and cc.state for analysis 
cell.type <- "Vip"
cc.state <- "prog"

df.E12.5 <-SetAllIdent(df.E12.5, id = "cell.id2")
df.temp <- SubsetData(df.E12.5, ident.use = cell.type)
# df.temp <- df.E12.5

df.temp <- SetAllIdent(df.temp, id = "cc.state")
if(cc.state == "neuron"){
  res <- SubsetData(df.temp, ident.use = "postmitotic")
} else if(cc.state == "prog"){
  res <- SubsetData(df.temp, ident.remove = "postmitotic")
} else {cat("Error: cc.state must be either neuron or prog")}

res <- SetAllIdent(res, id = "cell.id")
# res <- SetAllIdent(res, id = "clusters")
# res@ident <- factor(res@ident, levels = c("Pva", "Sst", "Vip", "Ndn"))
markers.temp <- FindAllMarkers(object = res,
                               #genes.use = res@var.genes,
                               only.pos = TRUE,
                               min.pct = 0.1,
                               #min.diff.pct = ,
                               return.thresh = 0.1)

top20.temp <- markers.temp %>% group_by(cluster) %>% top_n(20, avg_diff)

# pdf(file = paste(Sys.Date(), "prog",
#                  "major_lineage_heatmap_top10.pdf", sep = "_"),
#     width = 20, height = 10, useDingbats = F);
# DoHeatmap(object = res,
#           genes.use = top10.temp$gene,
#           slim.col.label = TRUE,
#           group.spacing = 0.3,
#           remove.key = TRUE)
# dev.off()

mean.res <- AverageExpression(res, return.seurat = T)
# gene.ord <- order(factor(top20.temp$cluster, levels = c("Pva", "Sst", "Vip", "Ndn")))
gene.ord <- 1:nrow(top20.temp)

# rub.genes <- read.csv("Rubenstein PGs vs Neuruon markers 03-11-2017.csv", header = F, 
#                       stringsAsFactors = F)
temp <- mean.res@scale.data[top20.temp$gene[gene.ord],]

ColSideColors <- c(gg_color_hue(length(unique(mean.res@ident))))#,

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

pdf(file = paste(Sys.Date(), cc.state, 
                 cell.type, 
                 "E12.5_heatmap_genes.pdf", sep = "_"),
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
          dendrogram = "column",
          #scale = "row",
          #colsep = colsep,
          sepcolor = "black",
          labRow = substr(rownames(temp), 20, 100),
          #labCol = "",
          na.rm = F);
dev.off(dev.cur());


ident.1 = "SstChodl"
ident.2 = "SstOthers"
res <- FindMarkers(df.E12.5, ident.1 = ident.1, ident.2 = ident.2, min.pct = 0.5)

res.markers <- rownames(res)[c(order(res$avg_diff)[1:100], 
                     order(res$avg_diff, decreasing = T)[1:100])]

pdf(file = paste(Sys.Date(),
                 samples[[i]][2], #samples[[i]][1],
                 ident.1, ident.2, "clusters_heatmap.pdf", sep = "_"),
    width = 20, height = 30, useDingbats = F);
DoHeatmap(object = df.E12.5,
          cells.use = df.E12.5@cell.names[df.E12.5@ident %in% c(ident.1, ident.2)],
          genes.use = res.markers,#top10$gene,#rownames(df.markers)[1:100],
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

#### 
df1 <- df.E12.5
df2 <- df.E14.5
hvg.E12.5 <- rownames(x = head(x = df1@hvg.info, n = 500))
hvg.E14.5 <- rownames(x = head(x = df2@hvg.info, n = 500))
hvg.union <- union(x = hvg.E12.5, y = hvg.E14.5)

res <- RunCCA(object = df1, object2 = df2, genes.use = hvg.union)

unique(res@meta.data$time_point)

p1 <- DimPlot(object = res, reduction.use = "cca", group.by = "time_point", pt.size = 3, 
              do.return = TRUE)
p2 <- VlnPlot(object = res, features.plot = "CC1", group.by = "time_point", do.return = TRUE)
pdf(file = paste(Sys.Date(), #samples[[i]][2], 
                 #samples[[i]][1],
                 "cca.pdf", sep = "_"),
    width = 40, height = 20, useDingbats = F);
plot_grid(p1, p2)
dev.off()

# PrintDim(object = res, reduction.type = "cca", dims.print = 1, genes.print = 50)
pdf(file = paste(Sys.Date(), #samples[[i]][2], 
                 #samples[[i]][1], 
                 "cca_dims.pdf", sep = "_"),
    width = 40, height = 40, useDingbats = F);
DimHeatmap(object = res, reduction.type = "cca", cells.use = 500, dim.use = 1:20, 
           do.balanced = TRUE)
dev.off()

res <- CalcVarExpRatio(object = res, reduction.type = "pca", 
                       grouping.var = "time_point", dims.use = 1:9)

# We discard cells where the variance explained by CCA is <2-fold (ratio <
# 0.5) compared to PCA
res.all.save <- res
res <- SubsetData(object = res, subset.name = "var.ratio.pca", accept.low = 0.5)

res.discard <- SubsetData(object = res.all.save, subset.name = "var.ratio.pca", 
                           accept.high = 0.5)
#median(x = res@meta.data[, "nGene"])
#median(x = res.discard@meta.data[, "nGene"])

# VlnPlot(object = res.discard, features.plot = "PF4", group.by = "time_point")

res <- AlignSubspace(object = res, reduction.type = "cca", 
                     grouping.var = "time_point", dims.align = 1:9)

res <- RunTSNE(object = res, reduction.use = "cca.aligned", dims.use = 1:9, 
                do.fast = TRUE,seed.use = 1)
res <- FindClusters(object = res, reduction.type = "tsne", dims.use = 1:2, 
                     save.SNN = TRUE,random.seed = 1)
p1 <- TSNEPlot(object = res, group.by = "time_point", do.return = TRUE, pt.size = 5)
p2 <- TSNEPlot(object = res, do.return = TRUE, pt.size = 5)
p3 <- TSNEPlot(object = res, group.by = "region", do.return = TRUE, pt.size = 5)
pdf(file = paste(Sys.Date(), #samples[[i]][2], 
                 #samples[[i]][1],
                 "tSNE.pdf", sep = "_"),
    width = 60, height = 20, useDingbats = F);
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

#### Do monocle ####
library(monocle)
source("/home/zl242/Scripts/R/order_cells.R")

#### Analyze two time points ####
pd <- new("AnnotatedDataFrame", data = res@meta.data)
fd <- data.frame(gene_ID = rownames(res@scale.data), 
                 gene_short_name = substr(rownames(res@scale.data), start = 20, stop = 100))
rownames(fd) <- fd$gene_ID
fd <- new("AnnotatedDataFrame", data = fd)
mm <- newCellDataSet(as(as.matrix(res@scale.data), "sparseMatrix"),
                     phenoData = pd, 
                     featureData = fd,
                     lowerDetectionLimit = min(res@scale.data), 
                     expressionFamily = gaussianff())

res@ident <- as.factor(res@meta.data$time_point)
names(res@ident) <- res@cell.names

time.markers <- FindMarkers(res, ident.1 = "E12.5", ident.2 = "E14.5", min.pct = 0.5)
time.markers <- time.markers[c(order(time.markers$avg_diff)[1:100], 
                               order(time.markers$avg_diff, decreasing = T)[1:100]),]

ddrtree_res <- DDRTree(t(res@dr$cca.aligned@cell.embeddings[,1:8]),
                       dimensions = 2)
diffusionMap::diffuse()

dev.markers <- c("ENSMUSG00000066392|NRXN3","ENSMUSG00000031285|DCX","ENSMUSG00000006586|RUNX1T1","ENSMUSG00000026787|GAD2","ENSMUSG00000090063|DLX6OS1","ENSMUSG00000021087|RTN1",
"ENSMUSG00000042834|D0H4S114",      "ENSMUSG00000070880|GAD1",        "ENSMUSG00000062380|TUBB3",       "ENSMUSG00000092341|MALAT1",     
"ENSMUSG00000028078|DCLK2",       "ENSMUSG00000024501|DPYSL3",      "ENSMUSG00000027500|STMN2",       "ENSMUSG00000072235|TUBA1A",      "ENSMUSG00000055435|MAF",        
"ENSMUSG00000043004|GNG2",        "ENSMUSG00000026890|LHX6",        "ENSMUSG00000062209|ERBB4",       "ENSMUSG00000061911|MYT1L",       "ENSMUSG00000029121|CRMP1",      
"ENSMUSG00000051910|SOX6",        "ENSMUSG00000032740|CCDC88A",     "ENSMUSG00000015291|GDI1",        "ENSMUSG00000015222|MTAP2",       "ENSMUSG00000034647|ANKRD12",    
"ENSMUSG00000026950|NEB",         "ENSMUSG00000025892|GRIA4",       "ENSMUSG00000052727|MTAP1B",      "ENSMUSG00000035277|ARX",         "ENSMUSG00000044647|CSRNP3",     
"ENSMUSG00000035864|SYT1",        "ENSMUSG00000029819|NPY",         "ENSMUSG00000027799|NBEA",        "ENSMUSG00000000861|BCL11A",      "ENSMUSG00000036006|FAM65B",     
"ENSMUSG00000028926|CDK14",       "ENSMUSG00000005871|APC",         "ENSMUSG00000048027|RGMB",        "ENSMUSG00000020598|NRCAM",       "ENSMUSG00000028399|PTPRD",      
"ENSMUSG00000024269|5730494M16RIK", "ENSMUSG00000052516|ROBO2",       "ENSMUSG00000057182|SCN3A",       "ENSMUSG00000018451|6330403K07RIK", "ENSMUSG00000025237|PARP6",      
"ENSMUSG00000087490|A330076H08RIK", "ENSMUSG00000038170|PDE4DIP",     "ENSMUSG00000067786|NNAT",        "ENSMUSG00000050708|FTL1",        "ENSMUSG00000069601|ANK3",       
"ENSMUSG00000020649|RRM2",  "ENSMUSG00000026385|DBI",   "ENSMUSG00000032397|TIPIN", "ENSMUSG00000027018|HAT1",  "ENSMUSG00000028312|SMC2",  "ENSMUSG00000028691|PRDX1",
"ENSMUSG00000027239|MDK",   "ENSMUSG00000020914|TOP2A", "ENSMUSG00000005732|RANBP1","ENSMUSG00000022033|PBK",   "ENSMUSG00000028194|DDAH1", "ENSMUSG00000016319|SLC25A5",
"ENSMUSG00000026355|MCM6",  "ENSMUSG00000037474|DTL",   "ENSMUSG00000025001|HELLS", "ENSMUSG00000022673|MCM4",  "ENSMUSG00000004642|SLBP",  "ENSMUSG00000028693|NASP", 
"ENSMUSG00000031004|MKI67", "ENSMUSG00000063229|LDHA",  "ENSMUSG00000001270|CKB",   "ENSMUSG00000019874|FABP7", "ENSMUSG00000029247|PAICS", "ENSMUSG00000022122|EDNRB",
"ENSMUSG00000034349|SMC4",  "ENSMUSG00000024640|PSAT1", "ENSMUSG00000021714|CENPK", "ENSMUSG00000028873|CDCA8", "ENSMUSG00000028884|RPA2",  "ENSMUSG00000025395|PRIM1",
"ENSMUSG00000000184|CCND2", "ENSMUSG00000027469|TPX2",  "ENSMUSG00000002870|MCM2",  "ENSMUSG00000027306|NUSAP1","ENSMUSG00000029730|MCM7",  "ENSMUSG00000020415|PTTG1",
"ENSMUSG00000030677|KIF22", "ENSMUSG00000038943|PRC1",  "ENSMUSG00000075266|CENPW", "ENSMUSG00000026701|PRDX6", "ENSMUSG00000031353|RBBP7", "ENSMUSG00000022678|NDE1", 
"ENSMUSG00000042489|CLSPN", "ENSMUSG00000027715|CCNA2", "ENSMUSG00000015880|NCAPG", "ENSMUSG00000018362|KPNA2", "ENSMUSG00000030662|IPO5",  "ENSMUSG00000045328|CENPE",
"ENSMUSG00000017499|CDC6",  "ENSMUSG00000029910|MAD2L1")

ddrtree_res <- DDRTree(res@scale.data[dev.markers,], dimensions = 2)

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
pData(mm)$Pseudodev <- cc_ordering$pseudo_time

ggplot(pData(mm), aes(x = Pseudotime, y = Pseudodev, color = res@meta.data$time_point)) + 
  geom_point() +
  scale_color_discrete(name = "")

# plot_cell_trajectory(mm, color_by = "res.1")

# mm <- setOrderingFilter(mm, ordering_genes = df@var.genes)

S_matrix <- reducedDimS(mm)
data_df <- data.frame(t(S_matrix[1:2,]))
colnames(data_df) <- c("data_dim_1", "data_dim_2")
rownames(data_df) <- res@cell.names#df@cell.names
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
  geom_point(aes(color = region), na.rm = TRUE, size = 1)
g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                 y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                 yend = "target_prin_graph_dim_2"), size = 0.3, linetype = "solid", 
                      na.rm = TRUE, data = edge_df)
g <- g + theme_bw()
g <- g + facet_wrap(region~time_point, nrow = 3)
pdf(file = paste(Sys.Date(), #samples[[i]][2], 
                 #samples[[i]][1],
                 "projection_regionXtime_point.pdf", sep = "_"),
    width = 6, height = 9, useDingbats = F);
g
dev.off()

g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2)) + 
  geom_point(aes(color = pseudotime), na.rm = TRUE, size = 1, alpha = 0.8)
g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                 y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                 yend = "target_prin_graph_dim_2"), size = 0.3, linetype = "solid", 
                      na.rm = TRUE, data = edge_df)
g <- g + theme_bw()
# g <- g + facet_wrap(~res.1, nrow = 3)
pdf(file = paste(Sys.Date(), #samples[[i]][2], 
                 #samples[[i]][1],
                 "projection_region.pdf", sep = "_"),
    width = 4, height = 3, useDingbats = F);
g
dev.off()

#### Analyze E14.5 ####
p <- c(0,1,5,7,8) #E14.5 progenitors
df.E14.5@meta.data$cell.type <- "n"
df.E14.5@meta.data$cell.type[df.E14.5@meta.data$res.1 %in% p] <- "p"

i = 6
df <- dfs[[i]]
df <- df.E14.5
pd <- new("AnnotatedDataFrame", data = df@meta.data)
fd <- data.frame(gene_ID = rownames(df@scale.data), 
                 gene_short_name = substr(rownames(df@scale.data), start = 20, stop = 100))
rownames(fd) <- fd$gene_ID
fd <- new("AnnotatedDataFrame", data = fd)
mm <- newCellDataSet(as(as.matrix(df@scale.data), "sparseMatrix"),
                     phenoData = pd, 
                     featureData = fd,
                     lowerDetectionLimit = min(df@scale.data), 
                     expressionFamily = gaussianff())

#ddrtree_res <- DDRTree(t(df@dr$pca@cell.embeddings[,1:9]), dimensions = 2)
df@ident <- as.factor(df@meta.data$region)
names(df@ident) <- df@cell.names

DvV <- FindMarkers(df,genes.use = df@var.genes,ident.1 = "dMGE", min.pct = 0.5,
                   ident.2 = "vMGE",#min.diff.pct = 0.2,
                   thresh.use = 0.25)
DvC <- FindMarkers(df,genes.use = df@var.genes,ident.1 = "dMGE", min.pct = 0.5,
                   ident.2 = "CGE",#min.diff.pct = 0.2,
                   thresh.use = 0.5)
VvC <- FindMarkers(df,genes.use = df@var.genes,ident.1 = "vMGE", min.pct = 0.5,
                   ident.2 = "CGE",#min.diff.pct = 0.2,
                   thresh.use = 0.5)

region.markers <- unique(c(rownames(DvV), rownames(DvC), rownames(VvC)))

df@ident <- as.factor(df@meta.data$cell.type)
names(df@ident) <- df@cell.names
dev.markers <- FindMarkers(df, genes.use = df@var.genes,ident.1 = "p", ident.2 = "n", 
                           min.pct = 0.5)
dev.markers <- dev.markers[c(order(dev.markers$avg_diff)[1:50], 
                             order(dev.markers$avg_diff, decreasing = T)[1:50]),]


cds <- apply(df@scale.data[top10$gene,],1,function(xx){
  log2(xx - min(xx) + 1)
})

#ddrtree_res <- DDRTree(df@scale.data[region.markers,], dimensions = 2)
ddrtree_res <- DDRTree(df@scale.data[rownames(dev.markers),], dimensions = 2)
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
#pData(mm)$State <- cc_ordering$cell_state

ggplot(pData(mm), aes(x = Pseudotime.region, y = Pseudotime, color = df@meta.data$cell.type)) + 
  geom_point() +
  scale_color_discrete(name = "")


# plot_cell_trajectory(mm, color_by = "res.1")

# mm <- setOrderingFilter(mm, ordering_genes = df@var.genes)

S_matrix <- reducedDimS(mm)
data_df <- data.frame(t(S_matrix[1:2,]))
colnames(data_df) <- c("data_dim_1", "data_dim_2")
rownames(data_df) <- df@cell.names
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
  geom_point(aes(color = region), na.rm = TRUE)
g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                 y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                 yend = "target_prin_graph_dim_2"), size = 0.3, linetype = "solid", 
                      na.rm = TRUE, data = edge_df)
g <- g + theme_bw()
g <- g + facet_wrap(~region, nrow = 1)
g

g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2)) + 
  geom_point(aes(color = region), na.rm = TRUE)
g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                 y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                 yend = "target_prin_graph_dim_2"), size = 0.3, linetype = "solid", 
                      na.rm = TRUE, data = edge_df)
g <- g + theme_bw()
# g <- g + facet_wrap(~res.1, nrow = 3)
g

cds <- apply(df@scale.data[top10$gene,],1,function(xx){
  log2(xx - min(xx) + 1)
})

cds <- df@scale.data[top10$gene,]
cds.rowMedians <- rowMedians(cds)
cds[cds < cds.rowMedians] <- NA
cds <- t(cds)

colnames(cds) <- substr(colnames(cds), 20, 100)
cds <- cbind(cds, pData(mm)[,7:ncol(pData(mm))])
cds <- melt(cds,id.vars = colnames(pData(mm))[7:ncol(pData(mm))])


pdf(file = paste(Sys.Date(), "E14.5_genes_in_pseudotime.pdf", sep = ""), 
    width = 20, height = 20, useDingbats = F)
g <- ggplot(data = cds, aes(x = Pseudotime, y = value)) + 
  geom_point(aes(color = res.1), na.rm = TRUE) +
  geom_smooth() +
  facet_wrap(~variable)
g
dev.off()

V(minSpanningTree(mm))[which(degree(minSpanningTree(mm)) > 2)]$name

g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                 y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                 yend = "target_prin_graph_dim_2"), size = 0.3, linetype = "solid", 
                      na.rm = TRUE, data = edge_df)
g <- g + theme_bw()
# g <- g + facet_wrap(~res.1, nrow = 3)
g


cds <- mm[genes[1],]
pdf(file = paste("./", name, "/", Sys.Date(), "_", name, "_genes_in_pseudotime_", suffix, ".pdf", sep = ""), 
    width = 11, height = 10, useDingbats = F)
plot_genes_in_pseudotime(cds, color_by="Pseudotime", ncol = 5)
dev.off()

g
#### Spec analysis

spec.scores <- read.csv("2017-10-09_E12.5_spec_scores_all_clusters.csv", header = T, 
                        check.names = F)
spec.scores <- spec.scores[, -1]  
spec.genes <- apply(spec.scores,1,function(xx){
  names(sort(xx, decreasing = T))[1:10]
})
as.vector(spec.genes)

# grep("\\bTH$",rownames(df@scale.data), value = T)
  plot.gene  <- c("ENSMUSG00000055435|MAF", "ENSMUSG00000074622|MAFB", 
    "ENSMUSG00000028222|CALB1", "ENSMUSG00000000214|TH")
  grep("")
  plot.gene <- c("ENSMUSG00000010175|PROX1")
  pdf(file = paste(Sys.Date(), 
                   samples[[i]][2], #samples[[i]][1], 
                   "genes.pdf", sep = "_"),
      width = 20, height = 20, useDingbats = F);
  FeaturePlot(df, features.plot = plot.gene, nCol = 1,
              cols.use = c("lightgrey","blue"), pt.size = 5)
  dev.off()
  
  temp <- df@scale.data[as.vector(spec.genes),]
  ord <- order(df@ident)
  ColSideColors <- c(gg_color_hue(length(unique(df@ident)))[as.numeric(as.factor(df@ident))],
                     rainbow(3)[as.numeric(as.factor(df@meta.data$region))])
  
  ColSideColors <- matrix(data = ColSideColors, nrow = ncol(temp));
  ColSideColors <- ColSideColors[ord, ];
  
  colsep <- lapply(unique(df@ident[ord]), function(x){length(which(df@ident == x))});
  colsep <- cumsum(unlist(colsep));  
  
  pairs.breaks <- seq(-2, 2, length.out=101);
  
  pdf(file = paste(Sys.Date(), 
                   samples[[i]][2], #samples[[i]][1],
                   "heatmap_spec_genes.pdf", sep = "_"),
      width = 10, height = 10, useDingbats = F);
  heatmap.3(temp[,ord],
            breaks = pairs.breaks,
            #symbreaks = T,
            keysize = 0.8, 
            main="Marker gene expression",
            col = bluered(100),
            symkey = F,
            cexRow=0.5, 
            cexCol = 0.6, 
            Rowv = F, 
            Colv = F,#as.dendrogram(hca), 
            ColSideColors = as.matrix(ColSideColors),
            ColSideColorsSize = 2,
            dendrogram = "column",
            scale = "row",
            colsep = colsep,
            sepcolor = "black",
            labRow = substr(rownames(temp), 20, 100),
            labCol = "",
            na.rm = F);
  dev.off(dev.cur());
  
  
  
  ggplot(df@meta.data, aes(df@ident, fill = region)) +
    geom_histogram(aes(y = ..prop..), stat = "count")
  
  pdf(file = paste(Sys.Date(), samples[[i]][2], samples[[i]][1],"gene_list_for_clusters.pdf", 
                   sep = "_"), 
      width = 18, height = 90, useDingbats = F);
  VlnPlot(df, features.plot = plot.gene, y.log = T,point.size.use = 0)
  dev.off()
}


  m <- t(df@dr$pca@cell.embeddings[,1:9]) # E12.5 dMGE = 1:3; vMGE = 1:9; CGE = 1:7; E14.5 dMGE = 1:3; vMGE = 1:5; CGE = 1:6
  hc <- hclust(as.dist(1-cor(m, method = "pearson", use = "pairwise.complete.obs")), method = "ward.D2")
  cl <- cutree(hc, k=2)
  
  df@ident <- as.factor(cl)
  
  
  
  #### 1.1.4.1 Clustering using Seurat ####
  k.param <- c(1:20)
  B = 1000
  n = length(df@cell.names)
  gap <- matrix(nrow = length(k.param), ncol = 6)
  for(ii in 1:length(k.param)){
    gap[ii,1] <- k.param[ii]
    df <- FindClusters(object = df, 
                       reduction.type = "pca", 
                       dims.use = 1:3, 
                       resolution = 0.9, 
                       print.output = 0, 
                       save.SNN = TRUE, 
                       force.recalc = T, 
                       random.seed = 1,
                       k.param = k.param[ii])
    xs <- t(df@scale.data[df@var.genes,])
    W <- getW(dataset = xs, clus = df@ident)
    V.sx <- svd(xs)$v
    rng.x1 <- apply(xs %*% V.sx, 2, range)
    logWks <- c()
    #if (verbose) 
    #  cat("Bootstrapping, b = 1,2,..., B (= ", B, ")  [one \".\" per sample]:\n", 
    #      sep = "")
    for (b in 1:B) {
      z1 <- apply(rng.x1, 2, function(M, nn) runif(nn, min = M[1], 
                                                   max = M[2]), nn = n)
      z <- tcrossprod(z1, V.sx)
      logWks <- c(logWks, log(getW(dataset = z, clus = df@ident)))
      cat("=")}
    #if (verbose) 
    #  cat(".", if (b%%50 == 0) 
    #    paste(b, "\n"))
    gap[i, 2] <- log(W)
    gap[i, 3] <- mean(logWks)
    gap[i, 4] <- sqrt((1 + 1/B) * var(logWks))
    gap[i, 5] <- length(unique(df@ident))
    cat("\n\n")  
  }
  
  gap[,2] <- log(gap[,2])
  gap[,6] <- gap[,3] - gap[,2]
  gap <- data.frame(gap)
  colnames(gap) <- c("k.param", "logW", "logEW", "SE.sim", "nclust", "gap")
  gap <- gap[order(gap$gap), ]
  gap_elbow <- vector()
  for(i in 1:nrow(gap)){
    gap_elbow[i] <- (gap$gap[i] - gap$gap[i+1] + gap$SE.sim[i+1])
  }
  gap$gap_elbow <- gap_elbow
  
  write.csv(gap, file = "2017-09-12_gap_E12.5")
  
  #### 1.1.5 tSNE analysis ####
  seed.use = 3 # All 3
  df <- RunTSNE(object = df,
                dims.use = 1:3, # E12.5 dMGE = 1:3; vMGE = 1:9; CGE = 1:7; E14.5 dMGE = 1:3; vMGE = 1:5; CGE = 1:6
                do.fast = T, 
                seed.use = seed.use,
                perplexity = 20, # dMGE = 20
                max_iter = 1000)
  
  # note that you can set do.label=T to help label individual clusters
  pdf(paste(paste(samples[[i]][2],"_", samples[[i]][1],"/",Sys.Date(), sep = ""), clu, 
            samples[[i]][2], samples[[i]][1],"tSNE.pdf", sep = "_"), 
      width = 11, height = 10, useDingbats = F);
  TSNEPlot(object = df, pt.size = 3)
  dev.off() 
  
  pdf(file = paste("./", Sys.Date(), "_gap.pdf", sep = ""), 
      width = 10, height = 10, useDingbats = F)
  ggplot(data = gap, aes(x = k.param, y = gap_elbow)) + 
    geom_point() + 
    geom_line(color = "black") +
    stat_smooth() +
    ggtitle(paste('Gap elbow')) +
    theme(plot.title = element_text(size=20, face="bold", vjust=2),
          panel.background =  element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA, colour="grey50"), 
          panel.grid.major =  element_line(colour = "grey90", size = 0.2),
          panel.grid.minor =  element_line(colour = "grey98", size = 0.5),
          axis.text.x=element_text(angle=90, size=10),
          axis.title.x = element_blank())
  dev.off(which=dev.cur())
  
  temp <- data.frame(df@dr$dm@cell.embeddings)
  g <- ggplot(data = temp, aes(x = DM1, y = DM2)) +
    geom_point(size = 5, color = as.numeric(df@ident) + 1) +
    #scale_color_brewer(palette = "Paired", name = 'Ages') +
    scale_color_discrete(name = "cell_type") +
    #geom_text(aes(label = colnames(Y2)), hjust = 1, vjust = -1) + 
    labs(x = "DM 1", y = "DM 2") +
    #ggtitle(paste("All cells, Perplexity =", p)) +
    theme_bw() +
    theme(text = element_text(size = 20),
          line = element_line(size = 1),
          panel.grid = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black"),
          legend.key = element_blank(),
          legend.background = element_blank())
  plot(g)
  
  #### 1.1.3 Principal component analysis ####
  pdf(file = paste("./E12.5_dMGE/", Sys.Date(), "_PCA_E12.5_dMGE.pdf", sep = ""), 
      width = 11, height = 10, useDingbats = F);
  PCAPlot(df, group.by = "IFC")
  dev.off()
  
  pdf(file = paste("./E12.5_dMGE/", Sys.Date(), "_PC1_Heatmap_E12.5_dMGE.pdf", sep = ""), 
      width = 11, height = 10, useDingbats = F);
  PCHeatmap(object = df, pc.use = 1, do.balanced = TRUE, label.columns = FALSE, remove.key = TRUE)
  dev.off()
  
  ## Use Jack Straw to determine significant PC
  df <-JackStraw(df, num.replicate = 100, do.print = F,num.pc = 20)
  
  pdf(file = paste("./E12.5_dMGE/", Sys.Date(), "_JackStraw_E12.5_dMGE.pdf", sep = ""), 
      width = 11, height = 10, useDingbats = F);
  JackStrawPlot(df, PCs = 1:20)
  dev.off()
  
  PCElbowPlot(df, num.pc = 20)
  
  # #### 1.1.4.2 Clustering using clustering and classification ####
  # temp <- read.table("./", "/", "_outputfile.txt", stringsAsFactors = F)
  # clu <- lapply(unique(temp$V1), function(x){tail(which(temp$V1 == x), n = 1)})
  # clu <- unlist(clu)
  # clu <- temp[clu,]
  # clu <- clu[clu$V1 %in% df@cell.names,]
  # clu <- clu[na.omit(match( df@cell.names, clu$V1 )),]
  # 
  # df@meta.data$cc_clusters <- clu$V2
  # 
  #load(paste("2017-05-29", name, "tSNE.Rdata", sep = "_"))
  
  #### 1.1.6 Find MzVsSvz ####
  df@ident <- as.factor(clu$V2)
  names(df@ident) <- df@cell.names
  
  # VvD <- FindMarkers(df, ident.1 = "vMGE", ident.2 = "dMGE", only.pos = T)
  df.markers <- FindAllMarkers(genes.use = df@var.genes, object = df, only.pos = TRUE, 
                               min.pct = 0.1, thresh.use = 0.25, return.thresh = 0.01)
  write.csv(df.markers, file = "2017-09-18_E12.5_vMGE_markers_2_clusters.csv")
  
  ## Plot top 20 marker genes
  top10 <- df.markers %>% group_by(cluster) %>% top_n(10, avg_diff)
  
  pdf(file = paste(Sys.Date(), 
                   samples[[i]][2], samples[[i]][1],"3_clusters_heatmap_top10.pdf", sep = "_"), 
      width = 20, height = 20, useDingbats = F);
  DoHeatmap(object = df, 
            genes.use = top10$gene, 
            slim.col.label = TRUE,group.spacing = 0.3,
            remove.key = TRUE)
  dev.off()
  
  #### 1.1.2 Find marker genes, lineage genes and variable genes ####
  markers <- read.csv("./raw_data/20170323_Rubenstein_revised_gene_list.csv", header = F, stringsAsFactors = F)
  markers <- unlist(markers[, 1 ])
  markers <- unique(sapply(markers, function(x){paste("\\b", x, "$", sep = "")}))
  markers <- unlist(sapply(markers,function(x){grep(x, rownames(df@scale.data), value = T)}))
  
  MzVsSvz <- read.csv("./raw_data/INs lineage marker 28-08-17.csv", header = T, stringsAsFactors = F)
  MzVsSvz <- unlist(MzVsSvz[, 1 ])
  MzVsSvz <- unique(sapply(MzVsSvz, function(x){paste("\\b", x, "$", sep = "")}))
  MzVsSvz <- unlist(sapply(MzVsSvz,function(x){grep(x, rownames(df@scale.data), value = T)}))
  
  lin <- read.csv("./raw_data/INs lineage marker 21-09-17.csv", header = F, stringsAsFactors = F)
  lin[,1] <- sapply(lin[,1], function(x){paste("\\b", x, "$", sep = "")})
  lin[,1] <- as.character(sapply(lin[,1],function(x){grep(x, rownames(df@scale.data), value = T)}))
  lin <- lin[!(lin$V1 == "character(0)"),]
  
  pdf(file = paste(Sys.Date(), 
                   samples[[i]][2],  "lineage.pdf", sep = "_"),
      width = 20, height = 40, useDingbats = F);
  FeaturePlot(df, features.plot = lin$V1, nCol = 6,
              cols.use = c("lightgrey","blue"), pt.size = 2)
  dev.off()
  
  glc <- read.csv("./raw_data/gene list for comparison.csv", header = F)
  glc$V1 <- toupper(glc$V1)
  glc$V1 <- sapply(glc$V1, function(x){paste("\\b", x, "$", sep = "")})
  glc$V1 <- sapply(glc$V1,function(x){grep(x, rownames(df@scale.data), value = T)})
  glc$V1 <- as.character(glc$V1)
  glc[glc$V1 %in% df.markers$gene,]
  temp <- cbind(df.markers[na.omit(match(glc$V1, df.markers$gene)),], glc[sort(na.omit(match(df.markers$gene, glc$V1))),])
  write.csv(temp, file = "2017-09-18_E12.5_dMGE_Da_Mi_list.csv")
  
  ## Plot top 20 marker genes
  top10 <- df.markers %>% group_by(cluster) %>% top_n(10, avg_diff)
  
  pdf(file = paste(Sys.Date(), 
                   samples[[i]][2], samples[[i]][1],"3_clusters_heatmap_top10.pdf", sep = "_"), 
      width = 20, height = 20, useDingbats = F);
  DoHeatmap(object = df, 
            genes.use = top10$gene, 
            slim.col.label = TRUE,group.spacing = 0.3,
            remove.key = TRUE)
  dev.off()
  
  for(ii in 1:length(unique(top10$cluster))){
    clu = unique(top10$cluster)[ii]
    pdf(file = paste(paste(samples[[i]][2],"_", samples[[i]][1],"/",Sys.Date(), sep = ""), 
                     samples[[i]][2], samples[[i]][1], clu, "tSNE_Top10.pdf", sep = "_"),
        width = 20, height = 20, useDingbats = F);
    FeaturePlot(df, features.plot = subset(top10, cluster == clu)$gene, 
                cols.use = c("lightgrey","blue"), pt.size = 2)
    dev.off()
  }
  
  for(ii in 1:length(unique(lin$V2))){
    clu = unique(lin$V2)[ii]
    pdf(file = paste(paste(samples[[i]][2],"_", samples[[i]][1],"/",Sys.Date(), sep = ""), clu, 
                           samples[[i]][2], samples[[i]][1],"tSNE_Lineage.pdf", sep = "_"),
        width = 20, height = 20, useDingbats = F);
    FeaturePlot(df, features.plot = subset(lin, V2 == clu)$V1,  
                cols.use = c("lightgrey","blue"), pt.size = 2)
    dev.off()
  }   
  
  write.csv(df.markers, file = paste(paste(samples[[i]][2],"_", 
                                           samples[[i]][1],"/",Sys.Date(), 
                                           sep = ""),
                                     samples[[i]][2], samples[[i]][1],"cc_DE.csv",
                                     sep = "_"))
  df@meta.data$clusters <- df@ident
  dfs[[paste(samples[[i]][2], samples[[i]][1], sep="_")]] <- df
  dfs.markers[[paste(samples[[i]][2], samples[[i]][1], sep="_")]] <- df.markers
}

grep("\\bNR2F2$", rownames(df@scale.data), value = T)

# genes <- c("ENSMUSG00000004151|ETV1", "ENSMUSG00000013089|ETV5", "ENSMUSG00000048562|SP8",
#   "ENSMUSG00000001496|NKX2-1", "ENSMUSG00000030551|NR2F2")

for(i in 1:length(dfs)){
  df <- dfs[[i]]
  MzVsSvz <- read.csv("./raw_data/mz_vs_svz_fdr5_logcpm5 by Lim.csv", header = T, stringsAsFactors = F)
  MzVsSvz <- unlist(MzVsSvz[, 1 ])
  MzVsSvz <- unique(sapply(MzVsSvz, function(x){paste("\\b", x, "$", sep = "")}))
  MzVsSvz <- unlist(sapply(MzVsSvz,function(x){grep(x, rownames(df@scale.data), value = T)}))
  pdf(file = paste(paste(samples[[i]][2],"_", samples[[i]][1],"/",Sys.Date(), sep = ""),
                   samples[[i]][2], samples[[i]][1],"heatmap_3_clusters_MzVsSvz_genes.pdf", sep = "_"),
      width = 10, height = 10, useDingbats = F);
  g <- DoHeatmap(object = df, 
            genes.use = as.character(MzVsSvz),
            slim.col.label = TRUE,
            group.spacing = 0.3,
            remove.key = TRUE)
  print(g)
  dev.off()
}

for(i in 1:length(dfs)){
  temp <- dfs[[i]]
  temp@ident <- as.factor(temp@meta.data$region)
  names(temp@ident) <- temp@cell.names
  df.markers <- FindAllMarkers(genes.use = temp@var.genes, object = temp, only.pos = TRUE, 
                               min.pct = 0.25, thresh.use = 0.25, return.thresh = 0.01)
  top20 <- df.markers %>% group_by(cluster) %>% top_n(20, avg_diff)
  pdf(file = paste(paste(samples[[i]][2],"_", samples[[i]][1],"/",Sys.Date(), sep = ""),
                   samples[[i]][2], samples[[i]][1],"heatmap_dMGE_vs_vMGE.pdf", sep = "_"),
      width = 10, height = 10, useDingbats = F);
  g <- DoHeatmap(object = temp, 
                 genes.use = top20$gene,
                 slim.col.label = TRUE,
                 group.spacing = 0.3,
                 remove.key = TRUE)
  print(g)
  dev.off()
}

df.14v <- df

### Cell cycle ####
cc_genes <- read.csv("dropSEQ_cell_cycle_genes.csv", header = T, na.strings = "")
cc_genes <- apply(cc_genes, 2, toupper)
cc_genes <- apply(cc_genes, 2, trimws)
cc_genes <- apply(cc_genes, 2, function(x){paste("\\b",x, "$", sep = "")})
cc_genes <- as.list(data.frame(cc_genes))
cc_genes <- sapply(cc_genes, function(x){x[x != "\\bNA$"]})
cc_genes <- sapply(cc_genes,function(x){unlist(sapply(x, function(xx) {grep(xx, rownames(df@raw.data), value = T)}))})

for(i in 1:length(dfs)){
  df <- dfs[[i]]
  MzVsSvz <- read.csv("./raw_data/mz_vs_svz_fdr5_logcpm5 by Lim.csv", header = T, stringsAsFactors = F)
  MzVsSvz <- unlist(MzVsSvz[, 1 ])
  MzVsSvz <- unique(sapply(MzVsSvz, function(x){paste("\\b", x, "$", sep = "")}))
  MzVsSvz <- unlist(na.omit(sapply(MzVsSvz,function(x){grep(x, rownames(df@scale.data), value = T)})))
  temp <- df@scale.data[MzVsSvz,]
  
  ord <- order(df@ident)
  ColSideColors <- c(gg_color_hue(length(unique(df@ident)))[as.numeric(as.factor(df@ident))])
   
  ColSideColors <- matrix(data = ColSideColors, nrow = ncol(temp));
  ColSideColors <- ColSideColors[ord, ];
   
  colsep <- lapply(unique(df@ident[ord]), function(x){length(which(df@ident == x))});
  colsep <- cumsum(unlist(colsep));  
  
  pairs.breaks <- seq(-2, 2, length.out=101);
  
  pdf(file = paste(paste(samples[[i]][2],"_", samples[[i]][1],"/",Sys.Date(), sep = ""),
                   samples[[i]][2], samples[[i]][1],"heatmap_3_clusters_MzVsSvz_genes.pdf", sep = "_"),
      width = 10, height = 10, useDingbats = F);
  heatmap.3(temp[,ord],
            breaks = pairs.breaks,
            #symbreaks = T,
            keysize = 0.8, 
            main="Marker gene expression",
            col = bluered(100),
            symkey = F,
            cexRow=1, cexCol = 0.6, 
            Rowv = T, 
            Colv = F,#as.dendrogram(hca), 
            ColSideColors = as.matrix(ColSideColors),
            ColSideColorsSize = 2,
            dendrogram = "column",
            scale = "row",
            #colsep = colsep,
            sepcolor = "black",
            labRow = substr(rownames(temp), 20, 100),
            labCol = "",
            na.rm = F);
  dev.off(dev.cur());
}