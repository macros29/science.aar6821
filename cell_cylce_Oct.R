library(ggplot2)
library(reshape2)
library(devtools)
library(Matrix)
library(Seurat)
library(dplyr)
library(statmod)
library("biomaRt")
library(scater)
library(scran)

setwd("C:/2017_singleCell/08_Apr")
options(stringsAsFactors = FALSE)

E12_expr <- read.csv("E12.5_expression.csv", header = T,row.names=1)

#select E12 dMGE
chip_170<-E12_expr[,grep('C1.170.',colnames(E12_expr))]
chip_173<-E12_expr[,grep('C1.173.',colnames(E12_expr))]
chip_176<-E12_expr[,grep('C1.176.',colnames(E12_expr))]
chip_184<-E12_expr[,grep('C1.184.',colnames(E12_expr))]
ncol(chip_170)+ncol(chip_173)+ncol(chip_176)+ncol(chip_184) # 287 cells
E12dMGE_readCounts<-cbind(chip_170,chip_173,chip_176,chip_184) # 287 cells
E12_expr<-E12dMGE_readCounts

#E12 vMGE: 5 chips
chip_171<-E12_expr[,grep('C1.171.',colnames(E12_expr))]
chip_174<-E12_expr[,grep('C1.174.',colnames(E12_expr))]
chip_177<-E12_expr[,grep('C1.177.',colnames(E12_expr))]
chip_183<-E12_expr[,grep('C1.183.',colnames(E12_expr))]
chip_185<-E12_expr[,grep('C1.185.',colnames(E12_expr))]
ncol(chip_171)+ncol(chip_174)+ncol(chip_177)+ncol(chip_183)+ncol(chip_185) # 
E12vMGE_readCounts<-cbind(chip_171,chip_174,chip_177,chip_183,chip_185) 
E12_expr<-E12vMGE_readCounts

#E12 CGE: 6 chips
chip_172<-E12_expr[,grep('C1.172.',colnames(E12_expr))]
chip_175<-E12_expr[,grep('C1.175.',colnames(E12_expr))]
chip_179<-E12_expr[,grep('C1.179.',colnames(E12_expr))]
chip_182<-E12_expr[,grep('C1.182.',colnames(E12_expr))]
chip_186<-E12_expr[,grep('C1.186.',colnames(E12_expr))]
chip_187<-E12_expr[,grep('C1.187.',colnames(E12_expr))]
ncol(chip_172)+ncol(chip_175)+ncol(chip_179)+ncol(chip_182)+ncol(chip_186)+ncol(chip_187) 
E12CGE_readCounts<-cbind(chip_172,chip_175,chip_179,chip_182,chip_186,chip_187) 
E12_expr<-E12CGE_readCounts

#Also read in a list of cell cycle markers
cc_genes <- read.csv("cell_cycle_genes.csv", header = F)
genes<-rownames(E12_expr)
gene_ID <- unlist(lapply(genes, function(x){strsplit(as.character(x), split = "\\|")[[1]][1]}))
gene_Name <- unlist(lapply(genes, function(x){strsplit(as.character(x), split = "\\|")[[1]][2]}))
df<-cbind(genes,as.data.frame(gene_ID),as.data.frame(gene_Name))
tmp <- subset(df, (df$gene_Name %in% cc_genes$V1))
cc_genes <- tmp[match(cc_genes$V1, tmp$gene_Name),]
#We can segregate this list into markers of G2/M phase and markers of G1/S phase
g1s.genes <- cc_genes$genes[1:43]
g2m.genes <- cc_genes$genes[44:98]
#Create our Seurat object and complete the initalization steps
mouse_GE_E12 <- CreateSeuratObject(raw.data = E12_expr)
mouse_GE_E12@data <- mouse_GE_E12@raw.data
mouse_GE_E12 <- FindVariableGenes(object = mouse_GE_E12, do.plot = FALSE, display.progress = FALSE)
mouse_GE_E12@scale.data <- mouse_GE_E12@raw.data
mouse_GE_E12 <- RunPCA(object = mouse_GE_E12, pc.genes = mouse_GE_E12@var.genes, do.print = FALSE)
#Calculate the cell cycle scores
mouse_GE_E12 <- CellCycleScoring(object = mouse_GE_E12, s.genes = g1s.genes, g2m.genes = g2m.genes, set.ident = TRUE)
head(x = mouse_GE_E12@meta.data)
mouse_GE_E12@meta.data$Phase <- gsub("G1", "Postmitotic", mouse_GE_E12@meta.data$Phase)
mouse_GE_E12@meta.data$Phase <- gsub("S", "G1/S", mouse_GE_E12@meta.data$Phase)
mouse_GE_E12 <- RunPCA(object = mouse_GE_E12, pc.genes = c(g1s.genes, g2m.genes), do.print = FALSE)
PCAPlot(object = mouse_GE_E12,group.by = "Phase")

write.csv(mouse_GE_E12@meta.data,file="E12CGE_cell_cycle_info.csv")

E12_cell_cycle_info<-read.csv("C:/2017_singleCell/11_Oct/E12_cell_cycle_info.csv", header = T,row.names=1)
colnames(E12_cell_cycle_info)[6] <- "total_phase"
E12dMGE_cell_cycle_info<-read.csv("C:/2017_singleCell/11_Oct/E12CGE_cell_cycle_info.csv", header = T,row.names=1)
colnames(E12dMGE_cell_cycle_info)[6] <- "CGE_phase"
E12_cell_cycle_info<-E12_cell_cycle_info[rownames(E12dMGE_cell_cycle_info),]
combindCellInfo <- cbind(E12_cell_cycle_info[,"total_phase",drop=FALSE],E12dMGE_cell_cycle_info[,"CGE_phase",drop=FALSE])


write.csv(combindCellInfo,file="C:/2017_singleCell/11_Oct/E12CGE_Compare_cell_cycle_info.csv")
summary(combindCellInfo$total_phase == combindCellInfo$CGE_phase)
#E12dMGE 6/281
#E12vMGE 17/341
#E12CGE  31/433

#E14dMGE 8/253
#E14vMGE 13/312
#E14CGE 64/315

#####################################################################################################
setwd("C:/2017_singleCell/08_Apr")
options(stringsAsFactors = FALSE)

E14_expr <- read.csv("E14.5_expression.csv", header = T,row.names=1)

#E14 vMGE: 4 chips
chip_91<-E14_expr[,grep('C1.91.',colnames(E14_expr))]
chip_93<-E14_expr[,grep('C1.93.',colnames(E14_expr))]
chip_96<-E14_expr[,grep('C1.96.',colnames(E14_expr))]
chip_98<-E14_expr[,grep('C1.98.',colnames(E14_expr))]
ncol(chip_91)+ncol(chip_93)+ncol(chip_96)+ncol(chip_98) 
E14vMGE_readCounts<-cbind(chip_91,chip_93,chip_96,chip_98) 
E14_expr<-E14vMGE_readCounts

#E14 dMGE: 5 chips
chip_90<-E14_expr[,grep('C1.90.',colnames(E14_expr))]
chip_92<-E14_expr[,grep('C1.92.',colnames(E14_expr))]
chip_97<-E14_expr[,grep('C1.97.',colnames(E14_expr))]
chip_99<-E14_expr[,grep('C1.99.',colnames(E14_expr))]
chip_101<-E14_expr[,grep('C1.101.',colnames(E14_expr))]
ncol(chip_90)+ncol(chip_92)+ncol(chip_97)+ncol(chip_99)+ncol(chip_101)
E14dMGE_readCounts<-cbind(chip_90,chip_92,chip_97,chip_99,chip_101) 
E14_expr<-E14dMGE_readCounts

#E14 CGE: 3 chips
chip_160<-E14_expr[,grep('C1.160.',colnames(E14_expr))]
chip_161<-E14_expr[,grep('C1.161.',colnames(E14_expr))]
chip_162<-E14_expr[,grep('C1.162.',colnames(E14_expr))]
chip_189<-E14_expr[,grep('C1.189.',colnames(E14_expr))]
chip_190<-E14_expr[,grep('C1.190.',colnames(E14_expr))]
ncol(chip_160)+ncol(chip_161)+ncol(chip_162)+ncol(chip_189)+ncol(chip_190)
E14CGE_readCounts<-cbind(chip_160,chip_161,chip_162,chip_189,chip_190) 
E14_expr<-E14CGE_readCounts

#Also read in a list of cell cycle markers
cc_genes <- read.csv("cell_cycle_genes.csv", header = F)
genes<-rownames(E14_expr)
gene_ID <- unlist(lapply(genes, function(x){strsplit(as.character(x), split = "\\|")[[1]][1]}))
gene_Name <- unlist(lapply(genes, function(x){strsplit(as.character(x), split = "\\|")[[1]][2]}))
df<-cbind(genes,as.data.frame(gene_ID),as.data.frame(gene_Name))
tmp <- subset(df, (df$gene_Name %in% cc_genes$V1))
cc_genes <- tmp[match(cc_genes$V1, tmp$gene_Name),]
#We can segregate this list into markers of G2/M phase and markers of G1/S phase
g1s.genes <- cc_genes$genes[1:43]
g2m.genes <- cc_genes$genes[44:98]
#Create our Seurat object and complete the initalization steps
mouse_GE_E14 <- CreateSeuratObject(raw.data = E14_expr)
mouse_GE_E14@data <- mouse_GE_E14@raw.data
mouse_GE_E14 <- FindVariableGenes(object = mouse_GE_E14, do.plot = FALSE, display.progress = FALSE)
mouse_GE_E14@scale.data <- mouse_GE_E14@raw.data
mouse_GE_E14 <- RunPCA(object = mouse_GE_E14, pc.genes = mouse_GE_E14@var.genes, do.print = FALSE)
#Calculate the cell cycle scores
mouse_GE_E14 <- CellCycleScoring(object = mouse_GE_E14, s.genes = g1s.genes, g2m.genes = g2m.genes, set.ident = TRUE)
head(x = mouse_GE_E14@meta.data)
mouse_GE_E14@meta.data$Phase <- gsub("G1", "Postmitotic", mouse_GE_E14@meta.data$Phase)
mouse_GE_E14@meta.data$Phase <- gsub("S", "G1/S", mouse_GE_E14@meta.data$Phase)
mouse_GE_E14 <- RunPCA(object = mouse_GE_E14, pc.genes = c(g1s.genes, g2m.genes), do.print = FALSE)
PCAPlot(object = mouse_GE_E14,group.by = "Phase")

write.csv(mouse_GE_E14@meta.data,file="E14CGE_cell_cycle_info.csv")

E14_cell_cycle_info<-read.csv("C:/2017_singleCell/11_Oct/E14_cell_cycle_info.csv", header = T,row.names=1)
colnames(E14_cell_cycle_info)[6] <- "total_phase"
E14dMGE_cell_cycle_info<-read.csv("C:/2017_singleCell/11_Oct/E14CGE_cell_cycle_info.csv", header = T,row.names=1)
colnames(E14dMGE_cell_cycle_info)[6] <- "CGE_phase"
E14_cell_cycle_info<-E14_cell_cycle_info[rownames(E14dMGE_cell_cycle_info),]
combindCellInfo <- cbind(E14_cell_cycle_info[,"total_phase",drop=FALSE],E14dMGE_cell_cycle_info[,"CGE_phase",drop=FALSE])


write.csv(combindCellInfo,file="C:/2017_singleCell/11_Oct/E14CGE_Compare_cell_cycle_info.csv")
summary(combindCellInfo$total_phase == combindCellInfo$CGE_phase)