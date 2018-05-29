

#### Combine datasets using SEURAT

library(dplyr)
library(Rtsne)
library(Seurat)
library(edgeR)

setwd('/Users/Gabriel/Desktop/YALE/BRAINSPAN/SingleCell/CLUSTERS_CLEAN/Seurat/')
load('../humanBrain.scRNAseq.sestanANDzhang.RData')

pheno<-as.data.frame(pheno)

### --> Select the genes to analyse

## Option 1. Specificity scores

specClust<-read.table('../CorrelLake/clustergroupspecificgenes/Marker_Genes_All_Clusters.txt')$V1
specGroup<-read.table('../CorrelLake/clustergroupspecificgenes/Marker_Genes_All_Groups_Trans.txt', sep='\t')$V1

markers<-unique(c(as.character(specGroup), as.character(specClust)))

markers <- unique(lapply(markers, function(x){paste("\\b", x, "$", sep = "")}))
markers <- unlist(sapply(markers,function(x){grep(x, rownames(geneCount), value = T)}))

geneCount[markers,] -> datam

### ---> Seurat Objects

brainspan<-which(pheno$batchList=="Sestan")
lake<-which(pheno$batchList=="Zhang")

BS <- CreateSeuratObject(raw.data = datam[,brainspan])
BS <- NormalizeData(object = BS)
BS <- ScaleData(object = BS)
BS <- FindVariableGenes(object = BS, do.plot = FALSE)

L <- CreateSeuratObject(raw.data = datam[,lake])
L <- NormalizeData(object = L)
L <- ScaleData(object = L)
L <- FindVariableGenes(object = L, do.plot = FALSE)

hvg.BS <- rownames(x = head(x = BS@hvg.info, n = 2000))
hvg.L <- rownames(x = head(x = L@hvg.info, n = 2000))
hvg.union <- union(x = hvg.BS, y = hvg.L)

BS@meta.data[, "protocol"] <- "Brainspan"
L@meta.data[, "protocol"] <- "Lake"

### ---> Calculate canonical correlation (CCA)

comb <- RunCCA(object = BS, object2 = L, genes.use = rownames(datam))

p1 <- DimPlot(object = comb, reduction.use = "cca", group.by = "protocol", pt.size = 0.5, 
              do.return = TRUE)
p2 <- VlnPlot(object = comb, features.plot = "CC1", group.by = "protocol", do.return = TRUE)
plot_grid(p1, p2)

PrintDim(object = comb, reduction.type = "cca", dims.print = 1:2, genes.print = 10)
DimHeatmap(object = comb, reduction.type = "cca", cells.use = 500, dim.use = 1:9, 
           do.balanced = TRUE)

comb <- CalcVarExpRatio(object = comb, reduction.type = "pca", grouping.var = "protocol", 
                        dims.use = 1:13)

# We discard cells where the variance explained by CCA is <2-fold (ratio <
# 0.5) compared to PCA
pbmc.all.save <- comb
pbmc <- SubsetData(object = comb, subset.name = "var.ratio.pca", accept.low = 0.5)

pbmc <- AlignSubspace(object = pbmc, reduction.type = "cca", grouping.var = "protocol", 
                      dims.align = 1:13)

p1 <- VlnPlot(object = pbmc, features.plot = "ACC1", group.by = "protocol", 
              do.return = TRUE)
p2 <- VlnPlot(object = pbmc, features.plot = "ACC2", group.by = "protocol", 
              do.return = TRUE)
plot_grid(p1, p2)

### New tSNE plots

pbmc <- RunTSNE(object = pbmc, reduction.use = "cca.aligned", dims.use = 1:13, 
                do.fast = TRUE)
pbmc <- FindClusters(object = pbmc, reduction.type = "cca.aligned", dims.use = 1:13, 
                     save.SNN = TRUE)
p1 <- TSNEPlot(object = pbmc, group.by = "protocol", do.return = TRUE, pt.size = 0.5)
p2 <- TSNEPlot(object = pbmc, do.return = TRUE, pt.size = 0.5)
plot_grid(p1, p2)

attributes(pbmc)$ident -> identity
pheno[which(pheno$sampList %in% names(identity)),] ->meta
meta[match(names(identity), meta$sampList),] -> meta

identity -> vec
vec[which(vec != 0)] <- 0

vec[which( meta$celltypeList == 'C9')] <- 1
vec[which( meta$celltypeList == 'C10')] <- 2
vec[which( meta$celltypeList == 'C11')] <- 3
vec[which( meta$celltypeList == 'C15')] <- 4
vec[which( meta$celltypeList == 'C16')] <- 5

vec[which( meta$celltypeList %in% paste0("Adult_Ex",c(1:8)))] <- 6
vec[which( meta$celltypeList %in% paste0("Adult_In",c(1:8)))] <- 7

attributes(pbmc)$ident <- vec

p1 <- TSNEPlot(object = pbmc, group.by = "protocol", do.return = TRUE, pt.size = 0.5)
p2 <- TSNEPlot(object = pbmc, do.return = TRUE, pt.size = 0.5, colors.use=c('white', 'darkgreen','green','blue','darkblue','cadetblue2', "darkgrey","black"))
plot_grid(p1, p2)
