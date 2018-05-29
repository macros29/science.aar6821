#If Windows OS, change the max DLL in R environment
#user_renviron = path.expand(file.path("~", ".Renviron"))
#file.edit(user_renviron)
#R_MAX_NUM_DLLS=400

#Use E12dMGE to test lineage diversity
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
library(WGCNA)
library(flashClust)


#drop GM, RPL, MT and Y genes
setwd("C:/2017_singleCell/08_Apr")
options(stringAsFactors = F); 
readCounts <- read.csv("2017-08-17_finalCleaned_readCounts.csv", header = T,row.names=1) #14078 genes  2669 cells
#convert RPKM to TPM: TPM = RPKM / (sum of RPKM over all genes/transcripts) * 10^6
expr <- readCounts 
genes<-rownames(expr)
gene_ID <- unlist(lapply(genes, function(x){strsplit(as.character(x), split = "\\|")[[1]][1]}))
df<-cbind(genes,as.data.frame(gene_ID))
rownames(df)<-df$gene_ID
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
filters = listFilters(ensembl)  
test <- gene_ID
ensembl_gene <-getBM(attributes=c("ensembl_gene_id","external_gene_name","chromosome_name","transcript_biotype"), filters = "ensembl_gene_id", values=test, mart=ensembl) #29131
expr_withoutIntron <- ensembl_gene[!(ensembl_gene$chromosome_name=="Y" | ensembl_gene$chromosome_name=="MT"),] #29098
expr_withoutIntron <- expr_withoutIntron[which(expr_withoutIntron$transcript_biotype=="protein_coding"),] #12505
tmp <- df[expr_withoutIntron$ensembl_gene_id,]
expr_withoutIntron <- expr[tmp$genes,] 
GM <- grep("\\bGM", rownames(expr_withoutIntron), value = T) #183 genes
#RPL <- grep("\\bRPL", rownames(expr_withoutIntron), value = T) #53 genes
#expr_withoutIntron <- expr_withoutIntron[!(rownames(expr_withoutIntron)%in%GM) & !(rownames(expr_withoutIntron)%in%RPL), ] 
expr_withoutIntron <- expr_withoutIntron[!(rownames(expr_withoutIntron)%in%GM), ] #12322 genes
expr_matrix<-expr_withoutIntron 
readCounts=expr_matrix
colnames(readCounts)<-gsub("(.*)_(.*)_(.*)_*_(.*)_(.*)_(.*)_*","\\1",colnames(readCounts))
#select E12 dMGE
chip_170<-readCounts[,grep('C1_170_',colnames(readCounts))]
chip_173<-readCounts[,grep('C1_173_',colnames(readCounts))]
chip_176<-readCounts[,grep('C1_176_',colnames(readCounts))]
chip_184<-readCounts[,grep('C1_184_',colnames(readCounts))]
ncol(chip_170)+ncol(chip_173)+ncol(chip_176)+ncol(chip_184) # 298 cells
E12dMGE_readCounts<-cbind(chip_170,chip_173,chip_176,chip_184) # 298  cells

#################################################################################################
#build Seurat object,# genes expressed in >= 3 cells (~1% of the data). 
mouseGE <- CreateSeuratObject(raw.data = E12dMGE_readCounts, min.cells = 3, min.genes = 1000, project = "MiDa") #11645 genes across 298 samples
mouseGE <- NormalizeData(object = mouseGE,normalization.method = "LogNormalize", scale.factor = 1e4)
#transcripts_per_cell <- apply(E12dMGE_readCounts, 2, function(x) sum(x > 0))
#transcripts_per_cell<-as.data.frame(transcripts_per_cell)
#transcripts_per_cell$batch <-     c(rep("Batch_1",ncol(chip_170)), 
#			                        rep("Batch_2",ncol(chip_173)), 
#			                        rep("Batch_3",ncol(chip_176)), 
#			                        rep("Batch_5",ncol(chip_184))
#			                        )
#transcripts_per_cell$batch <- as.factor(transcripts_per_cell$batch)
#mouseGE <- AddMetaData(object = mouseGE,metadata = transcripts_per_cell[,2,drop=FALSE],col.name = "batch")
#mouseGE <- ScaleData(object = mouseGE,vars.to.regress = c("batch")) # https://github.com/satijalab/seurat/issues/49
E12dMGE_normalisedCounts <- as.data.frame(as.matrix(mouseGE@data)) 


#input for WGCNA
linageMaker <- read.csv("INs lineage marker 28-08-17.csv", header = F) 
linageMaker<-as.character(linageMaker[,1])
linageMaker<- unique(linageMaker) #68 genes
mart <- useMart(biomart = 'ensembl', dataset = 'mmusculus_gene_ensembl')
bm.query <- getBM(values=linageMaker,attributes=c("ensembl_gene_id", "external_gene_name"),filters=c("external_gene_name"),mart=mart)
bm.query<-paste(bm.query$ensembl_gene_id, toupper(bm.query$external_gene_name), sep="|")

#WGCNA
datExpr0 <- E12dMGE_normalisedCounts[bm.query,]
datExpr0 <- datExpr0[complete.cases(datExpr0), ] #38 genes
#genes.var <- rowVars(datExpr0);tmp=cbind(rownames(datExpr0),as.data.frame(genes.var))
#datExpr0 <- datExpr0[genes.var > 0.01, ] 
#sce <- newSCESet(exprsData=datExpr0)
#set.seed(29)
#var.cor <- correlatePairs(sce, subset.row=mouseGE@var.genes); head(var.cor)
#sig.cor <- var.cor$FDR <= 0.05; summary(sig.cor)
#chosen <- unique(c(var.cor$gene1[sig.cor], var.cor$gene2[sig.cor]))
#datExpr0 <- datExpr0[chosen,] 
datExpr0 <- t(datExpr0) #need tran the frame

ALLOW_WGCNA_THREADS=2
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK #true

#1. cluster samples/cells, detect outliers
sampleTree = hclust(dist(datExpr0), method = "average") #for samples
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
	 
abline(h = 5, col = "red"); 
clust = cutreeStatic(sampleTree, cutHeight = 5, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#if we do not delete any sample, usually in scRNA-seq, we do not drop cells in this step
datExpr = datExpr0[, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#######################################################################################################
#2. choose suitable sft
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "signed")
# Scale-free topology fit index as a function of the soft-thresholding power
sizeGrWindow(9,5);par(mfrow=c(1,2));cex1=0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
abline(h=0.80,col="red") #the criterion that the co-efficiency of log(k) and log(p(k)) is at least 0.8
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 6  #根据图像选择拓扑结构树的soft power阈值,  by cor* softPower, reduce noise and emphasize the stronger correlations
adj= adjacency(datExpr,type = "signed", power = softPower) #计算树之间的邻接性by TOM 
TOM=TOMsimilarityFromExpr(datExpr,networkType = "signed", TOMType = "signed", power = softPower); #计算树之间的相似性, 有向(signed)  
colnames(TOM) =rownames(TOM) =colnames(datExpr)
dissTOM=1-TOM
geneTree = flashClust(as.dist(dissTOM),method="complete")
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04)

#3. TREE CUTTING AND MODULE ASSIGNMENT
minModuleSize = 3; # the smallest number of genes in a module
x = 4 # 1 2 3 4 灵敏度， 默认是2 scRNA-seq 一般上调 3 or 4
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = x, pamRespectsDendro = TRUE, minClusterSize = minModuleSize)
table(dynamicMods) #dynamic cut
#dynamicMods<-cutreeStatic(geneTree, cutHeight = 0.975, minSize = minModuleSize)
dynamicColors = labels2colors(dynamicMods);plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

#BUILD OBJECT NETWORK TO VIEW THE GENES
MEs_mouse = moduleEigengenes(datExpr, dynamicColors)$eigengenes #eigengenes 
ME1_mouse<-MEs_mouse;row.names(ME1_mouse)<-row.names(datExpr)
write.table(ME1_mouse, "26_Sep_WGCNA_unmerged_modueEigenene.txt", sep="\t", quote=F, col.names=NA)
print("WGCNA_unmerged_modueEigenene.txt ....... saved")
#BUILD OBJECT NETWORK TO VIEW THE GENES
modules<-MEs_mouse;c_modules<-data.frame(dynamicColors);row.names(c_modules)<-colnames(datExpr)

#4. display corraltation to module eigengene for everymember of each module
colors<-data.frame(unique(dynamicColors, incomparables = FALSE))
datExpr_tmp <- as.data.frame(t(datExpr))
for(i in 1:dim(colors)[1]){assign(paste(colors[i,1],"data",sep="_"),as.matrix(datExpr_tmp[which(dynamicColors==colors[i,1]),]))}
for(i in 1:dim(colors)[1]){scores <- prcomp(t(get(paste(colors[i,1],"data",sep="_"))), center = F)$x[,1];dat<-rbind(get(paste(colors[i,1],"data",sep="_")),scores);assign(paste(colors[i,1],"data_cor",sep="_"),cor(t(dat)))}
i<-1
get(paste(colors[i,1],"data_cor",sep="_"))["scores",];print(colors[i,1]);i<-i+1
#matrix of genes and assigned module colors
c<-data.frame(dynamicColors);row.names(c)<-colnames(datExpr)

#5.GENE INFO
#signedKME for eigengene based connectivity in a single data set. corAndPvalue, bicorAndPvalue for two alternatives for calculating correlations and the corresponding p-values and Z scores.
modNames=substring(names(MEs_mouse),3)
geneModuleBicor=bicorAndPvalue(datExpr, MEs_mouse, use = "p") #bicorAndPvalue used here, but signedKME is the most common one in bulk RNA seq
geneModuleMembership = as.data.frame(geneModuleBicor$bicor)
MMPvalue = as.data.frame(geneModuleBicor$p)
names(geneModuleMembership) = paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")
geneInfo = data.frame(miRID = colnames(datExpr), moduleColor = dynamicColors, geneModuleMembership, MMPvalue)

gene_ID <- unlist(lapply(rownames(geneInfo), function(x){strsplit(as.character(x), split = "\\|")[[1]][2]}))
geneInfo<-cbind(geneInfo,as.data.frame(gene_ID))
linage_type <- read.csv("INs lineage marker 28-08-17.csv", header = F) 
geneInfo$miRID <- linage_type[match(geneInfo$gene_ID, linage_type$V1),]$V2
geneInfo <- geneInfo[,-ncol(geneInfo)]
geneInfo <- geneInfo[order(geneInfo$moduleColor),]
write.table(geneInfo, "26_Sep_WGCNA_unmerged_geneInfo.txt", sep="\t", quote=F, col.names=NA)
print("WGCNA_unmerged_geneInfo.txt ....... saved")

#6.MERGE CLOSE MODULES: Not do this in scRNA-seq. Choose different input gene list is much useful

#7.Plot tSNE gradients plot, assuming we got tNSE embedding from Seurat
tSNEembedding <- as.data.frame(mouseGE@dr$tsne@cell.embeddings)
exprMatrix<-ME1_mouse
predict_marker<-colnames(ME1_mouse)

for (i in 1:length(predict_marker)){
     geneName = predict_marker[i]
     marker_expr<-as.data.frame(exprMatrix[,predict_marker[i]])
     plotdata <- data.frame(rownames(tSNEembedding),tSNEembedding$tSNE_1,tSNEembedding$tSNE_2,marker_expr) 
     colnames(plotdata) <-c("sample","PC1","PC2","expr_level")

     p<-ggplot(plotdata, aes(PC1, PC2,colour=expr_level))+geom_point(size=6)+scale_colour_gradientn(colours=c("darkgrey","yellow","orange","red")) +
	           labs(x=paste("tSNE dim1"),y=paste("tSNE dim2"),title=colnames(marker_expr)) + theme_bw() +
               theme(text = element_text(size = 10),
               line = element_line(size = 1),
               axis.line = element_line(colour = "black"),
               legend.key = element_blank(),
               panel.grid = element_blank())	   
		   
pdf(file = paste(geneName, "_tSNEplot.pdf", sep = ""),width = 15, height = 10, useDingbats = F)
print(p)
dev.off()	
}
