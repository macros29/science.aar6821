library(Seurat)
library(WGCNA)
library(igraph)

setwd("~/Data/Da_Mi/")
load("working.RData")

i = 5               
keep <- GB$region %in% samples[[i]] & GB$time_point %in% samples[[i]] & !(GB$IFC %in% c("180","181","191","192","193","194","195", "196"))
df <- CreateSeuratObject(raw.data = as.vector(expr[,keep]), min.cells = 5, min.genes = 1000, project = "Da_Mi")
df <- NormalizeData(object = df, normalization.method = "LogNormalize", scale.factor = 10000)
mito.genes <- grep(pattern = "\\bMT-", x = rownames(x = df@data), value = TRUE)
percent.mito <- colSums(df@raw.data[mito.genes, ])/colSums(df@raw.data)
# AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
df <- AddMetaData(object = df, metadata = percent.mito, col.name = "percent.mito")
df@meta.data <- cbind(df@meta.data, GB[df@cell.names,])
#df <- ScaleData(df)
df <- ScaleData(df,vars.to.regress = c("nGene", "mouse"))

df <- FindVariableGenes(object = df, mean.function = ExpMean, dispersion.function = LogVMR, 
                        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = F)

df <- RunPCA(object = df, pc.genes=df@var.genes, pcs.compute = 30, pcs.print = 1:12,genes.print = 10)

pcs <- list(c(1:3), c(1:7), c(1:4), c(1:4), c(1:4), c(1:4))
df <- RunTSNE(df, dims.use = pcs[[i]], seed.use = 7)

#### Optional: WGCNA ####
datExpr <- t(df@scale.data)

# Choose a set of soft-thresholding powers
powers = c(seq(from = 1, to=20, by=1))

# Call the network topology analysis function
# sft <- pickSoftThreshold(datExpr, 
#                          powerVector = powers, 
#                          verbose = 0, 
#                          corFnc = "bicor", 
#                          corOptions = list(maxPOutliers =0.1), 
#                          networkType = "signed hybrid")

sft <- pickSoftThreshold(datExpr, 
                         powerVector = powers, 
                         verbose = 0, 
                         corFnc = "cor", 
                         corOptions = list(use = "p", method = "p"), 
                         networkType = "signed")

# Scale-free topology fit index as a function of the soft-thresholding power
p1 <-ggplot(data = sft$fitIndices, aes(x = Power, y = -sign(slope)*SFT.R.sq)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  geom_hline(aes(yintercept = 0.9), colour = "red")

p2 <- ggplot(data = sft$fitIndices, aes(x = Power, y = mean.k.)) +
  geom_point(size = 3) +
  geom_line(size = 0.5) +
  ggtitle(label = "Mean connectivity")

gridExtra::grid.arrange(grobs = list(p1, p2), ncol = 2)

# Mean connectivity as a function of the soft-thresholding power
# ADJ1=abs(cor(datExpr,use="p"))^2
# k=as.vector(apply(ADJ1,2,sum, na.rm=T))
# sizeGrWindow(10,5)
# par(mfrow=c(1,2))
# hist(k)
# scaleFreePlot(k, main="Check scale free topology\n")

softPowers = c(10, 8, 9, 6, 8, 10)
softPower = softPowers[i]
# cor <- bicor(datExpr, use = "pairwise.complete.obs", maxPOutliers = 0.1)
cor <- cor(datExpr, use = "p")
adj = adjacency.fromSimilarity(cor, type = "signed", power = softPower)
TOM = TOMsimilarity(adj, TOMDenom = "min", TOMType = "signed")
colnames(TOM) <- rownames(TOM) <- colnames(datExpr)

#load("2016-02-19_iPSC_TOM_sfp_5.Rdata")
dissTOM <- 1 - TOM

geneTree <- hclust(as.dist(dissTOM),method="average");

# Set the minimum module size
minModuleSize = 60;

# Module identification using dynamic tree cut

dynamicMods = cutreeDynamic(dendro = geneTree,
                            method = "hybrid",
                            deepSplit = 2,
                            #minAbsSplitHeight = 0.999,
                            minClusterSize = minModuleSize,
                            distM = dissTOM);

#Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)
#table(dynamicColors)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, 
                    hang = 0.03, 
                    addGuide = TRUE, 
                    guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")

restGenes= (dynamicColors != "grey")
cor.rest <- cor(datExpr[,restGenes], use="p")
adj2 <- adjacency.fromSimilarity(cor.rest, type = "signed", power = softPower)
#adj2 <- adjacency(datExpr[,restGenes], type = "signed hybrid", power = softPower)
diss1=1-TOMsimilarity(adjMat =adj2, 
                      TOMType = "signed",
                      TOMDenom = "min")

colnames(diss1) =rownames(diss1) = colnames(datExpr)[restGenes]
hier1=hclust(as.dist(diss1), method="average" )
plotDendroAndColors(hier1, 
                    dynamicColors[restGenes], 
                    "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, 
                    guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")

# Correlate genes with module eigengenes
MEs <- moduleEigengenes(datExpr[,restGenes], dynamicColors[restGenes])
MEs <- orderMEs(MEs$eigengenes)

# Merge similar modules
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
MEDissThres = 0.95
abline(h=MEDissThres, col = "red", lty = 2)
merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
MEs <- merge$newMEs
dynamicColors <- merge$colors

temp <- data.frame(df@dr$tsne@cell.embeddings)

gs <- lapply(1:ncol(MEs), function(i){
  name <- colnames(MEs)[i]
  ggplot(data = temp, aes(x = tSNE_1, y = tSNE_2, color = MEs[,i], alpha = MEs[,i])) +
  geom_point(size = 5) +
  labs(x = "tSNE_1", y = "tSNE_2") +
  scale_color_continuous(high = "purple", low = "grey", name = name, limits = c(quantile(MEs[,i], probs = 0.01), 
                                                                                quantile(MEs[,i], probs = 0.99))) +
  scale_alpha_continuous(name = name, limits = c(quantile(MEs[,i], probs = 0.01), 
                                                 quantile(MEs[,i], probs = 0.99))) +
  ggtitle(label = name) +
  theme_bw() +
  theme(text = element_text(size = 20),
        line = element_line(size = 1),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.key = element_blank(),
        legend.position = "",
        legend.background = element_blank())
})

pdf(file = paste(Sys.Date(), names(dfs)[i],"MEs_tSNE.pdf", sep = "_"),
    width = 40, height = 55, useDingbats = F);
ml <- gridExtra::grid.arrange(grobs = gs, ncol = 3)
dev.off()

# Plot gene network 
Alldegrees <- intramodularConnectivity(TOM, dynamicColors) # Intramodular connectivity of all genes
probes = colnames(datExpr)

for(module in unique(dynamicColors)){
  if(module == "grey"){
    next
  } else {tryCatch({
    # Select module probes
    inModule = (dynamicColors==module);
    modProbes = probes[inModule];
    # Select the corresponding Topological Overlap
    modTOM = TOM[inModule, inModule];
    dimnames(modTOM) = list(modProbes, modProbes)
    colnames(modTOM) <- rownames(modTOM) <- substr(rownames(modTOM), 20, 100)
    
    threshold = 0.01
    
    vis = exportNetworkToVisANT(modTOM,
                                weighted = TRUE,
                                threshold = threshold)
    
    network <- graph_from_data_frame(vis, directed = TRUE)
    fr <- layout.fruchterman.reingold(network)
    fr <- data.frame(fr)
    colnames(fr) <- c("x", "y")
    fr$network_name <- V(network)$name
    #fr$IMC <- Alldegrees[inModule, 2]
    
    # Convert SNN to 2-D plot coordinates
    g <- get.data.frame(network)
    g$from.x <- fr$x[match(g$from, fr$network_name)]
    g$from.y <- fr$y[match(g$from, fr$network_name)]
    g$to.x <- fr$x[match(g$to, fr$network_name)]
    g$to.y <- fr$y[match(g$to, fr$network_name)]
    
    fr$size <- unlist(lapply(fr$network_name, function(x){length(which(g$from == x)) + length(which(g$to == x))}))
    
    pdf(file = paste(Sys.Date(), names(dfs)[i],"Module", module, threshold, "threshold.pdf", sep = "_"),
        width = 20, height = 18, useDingbats = F)
    p <- ggplot(data = fr, aes(x = x, y = y))+
      geom_segment(data = g, aes(x = from.x, xend = to.x, y = from.y, yend = to.y), colour = "grey90") +
      geom_point(data = fr, aes(x = x, y = y, size = size), colour = "green") +
      scale_size_continuous(range = c(1,30)) +
      geom_text(data = fr, aes(x = x, y = y, label = network_name), size = 5) +
      theme_bw() +
      theme(rect = element_blank(),
            line = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = "none") +
      ggtitle(paste("Module", module, "threshold =", threshold, sep = " ")) +
      theme(plot.title = element_text(size=20, face="bold", vjust=2))
    print(p)
  }, error=function(e){cat("Module: ", module, " failed\n")})
    dev.off(dev.cur())
  }
}

# Relating modules to external factors such as clusters
moduleTraitCor <- list()
moduleTraitPvalue <- list()
header = as.character(paste(GB$Species,  GB$Time_point, GB$IDX))
for(i in unique(header)[order(unique(header))]){
  h <- header
  h[header != i] <- 0
  h[header == i] <- 1
  moduleTraitCor[[i]] <- cor(MEs, h, use = "p")
  moduleTraitPvalue[[i]] <- corPvalueStudent(moduleTraitCor[[i]], nrow(datExpr))
}

moduleTraitCor <- do.call(cbind, moduleTraitCor)
moduleTraitPvalue <- do.call(cbind, moduleTraitPvalue)
colnames(moduleTraitCor) <- unique(header)[order(unique(header))]
colnames(moduleTraitPvalue) <- unique(header)[order(unique(header))]
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");

t <- matrix(textMatrix, nrow = nrow(moduleTraitCor), ncol = ncol(moduleTraitCor))
heatmap.3(moduleTraitCor,
          keysize = 0.8,
          cellnote = t,
          notecol = "black",
          #breaks=pairs.breaks,
          col = greenWhiteRed(100),
          Rowv=F,
          Colv=F,
          dendrogram = "none",
          RowSideColors=t(as.matrix(substr(rownames(moduleTraitCor), 3, 100))),
          #ColSideColors=matrix(ColSideColors, ncol = 3),
          #ColSideColorsSize = 2,
          main = "Correlation between module and clusters")
legend('left',
       legend = c("Chimp", "Human", "D14", "D40", "D80"), 
       col = c(unique(ColSideColors)),
       fill = c(unique(ColSideColors)),
       cex = 1)

# Find top 5 genes with highest correlation to each module
moduleGeneCor <- apply(datExpr, 2, 
                       function(x){cor(MEs,#[,!(colnames(MEs)%in%"MEgrey")], 
                                       x, use = "p")})
moduleGenePvalue <- apply(moduleGeneCor, 2, 
                          function(x){corPvalueStudent(x, length(which(restGenes)))})
rownames(moduleGeneCor) <- rownames(moduleGenePvalue) <- colnames(MEs)#[,!(colnames(MEs)%in%"MEgrey")])
moduleGeneCor[moduleGenePvalue > 0.01] <- NA
module_colors= substr(colnames(MEs), 3, 100)
moduleHighCorGenes <- apply(moduleGeneCor, 1, function(x){x <- x[which(x > 0.3)]; x <- x[order(x, decreasing = T)]; if(length(x) > 50){names(x[1:50])} else {names(x)}})

module_colors <- substr(rownames(moduleGeneCor), 3, 100)
RowSideColors <- c()
for(i in 1:length(moduleHighCorGenes)){
  RowSideColors <- c(RowSideColors, rep(module_colors[i], times = length(moduleHighCorGenes[[i]])))
}
# ColSideColors <- c(topo.colors(2)[as.numeric(as.factor(GB$Species))],
#                    gg_color_hue(3)[as.numeric(as.numeric(GB$Time_point)-1)],
#                    brewer.pal(length(unique(GB$IDX)), "Paired")[as.numeric(GB$IDX)])
# ColSideColors <- matrix(ColSideColors, ncol = 3)
# colsep <- c()
# colsep <- paste(GB$Species,  GB$Time_point, GB$IDX)[order(paste(GB$Species, GB$Time_point, GB$IDX))]
# colsep <- lapply(unique(colsep[order(colsep)]), function(x){length(which(colsep == x))})
# colsep <- cumsum(unlist(colsep))
rowsep <- cumsum(unlist(lapply(moduleHighCorGenes, length)))
rec <- Y2[unlist(moduleHighCorGenes), ))]
pairs.breaks <- seq(-2.5, 2.5, length.out=101);
pdf(file = paste("./images/", Sys.Date(), "_", cell_source, "_WGCNA_genes.pdf", sep = ""), width = 16, height = 10, useDingbats = F)
heatmap.3(rec,
          key=T, keysize = 0.8,
          col = bluered(100),
          breaks=pairs.breaks,
          density.info = "histogram", 
          symkey=F, 
          main="Marker gene expression", 
          trace="none", cexRow=0.6, cexCol = 0.6, 
          Rowv = F, 
          Colv = F, #as.dendrogram(hca), 
          dendrogram = "none",
          # ColSideColors = ColSideColors[order(paste(GB$Species, GB$Time_point, GB$IDX)),],
          #ColSideColors = ColSideColors[order(t_sne$kmeans),],
          ColSideColorsSize = 3,
          RowSideColors = t(as.matrix(RowSideColors)),
          RowSideColorsSize = 0.5,
          scale = "row",
          colsep = colsep,
          rowsep = rowsep,
          #labRow = substr(rownames(rec), start = 17, stop = 100),
          labRow = "",
          labCol = "",
          na.rm = F
);
dev.off(dev.cur())

# Extract modules
for (i in 1:length(module_colors)){
  module <- na.omit(moduleGeneCor[i,][order(moduleGeneCor[i,], decreasing = T)])
  module <- module[module > 0.3]
  write.csv(names(module), file = paste0(Sys.Date(), "_module_", module_colors[i], ".csv", sep=""))
}

# Plot gene network 
Alldegrees <- intramodularConnectivity(TOM, dynamicColors[restGenes]) # Intramodular connectivity of all genes
probes = colnames(datExpr)
for(module in unique(dynamicColors)){
  if(module == "grey"){
    next
  } else {tryCatch({
    # Select module probes
    inModule = (dynamicColors==module);
    modProbes = probes[inModule];
    # Select the corresponding Topological Overlap
    modTOM = TOM[inModule, inModule];
    dimnames(modTOM) = list(modProbes, modProbes)
    colnames(modTOM) <- rownames(modTOM) <- substr(rownames(modTOM), 17, 100)
    
    threshold = 0.13
    
    vis = exportNetworkToVisANT(modTOM,
                                weighted = TRUE,
                                threshold = threshold)
    
    network <- graph_from_data_frame(vis, directed = TRUE)
    fr <- layout.fruchterman.reingold(network)
    fr <- data.frame(fr)
    colnames(fr) <- c("x", "y")
    fr$network_name <- V(network)$name
    #fr$IMC <- Alldegrees[inModule, 2]
    
    # Convert SNN to 2-D plot coordinates
    g <- get.data.frame(network)
    g$from.x <- fr$x[match(g$from, fr$network_name)]
    g$from.y <- fr$y[match(g$from, fr$network_name)]
    g$to.x <- fr$x[match(g$to, fr$network_name)]
    g$to.y <- fr$y[match(g$to, fr$network_name)]
    
    fr$size <- unlist(lapply(fr$network_name, function(x){length(which(g$from == x)) + length(which(g$to == x))}))
    
    pdf(file = paste(Sys.Date(), "Module", module, threshold, "threshold.pdf", sep = "_"),
        width = 20, height = 18, useDingbats = F)
    p <- ggplot(data = fr, aes(x = x, y = y))+
      geom_segment(data = g, aes(x = from.x, xend = to.x, y = from.y, yend = to.y), colour = "grey90") +
      geom_point(data = fr, aes(x = x, y = y, size = size), colour = "green", alpha = 0.75) +
      geom_text(data = fr, aes(x = x, y = y, label = network_name, size = size), vjust = -0.5) +
      theme_bw() +
      theme(rect = element_blank(),
            line = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = "none") +
      ggtitle(paste("Module", module, "threshold =", threshold, sep = " ")) +
      theme(plot.title = element_text(size=20, face="bold", vjust=2))
    print(p)
  }, error=function(e){cat("Module: ", module, " failed\n")})
    dev.off(dev.cur())
  }
}

datExpr = t(df@scale.data[df@var.genes,])
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
minModuleSize = 50; # the smallest number of genes in a module
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
linage_type <- read.csv("raw_data/INs lineage marker 21-09-17.csv", header = F) 
geneInfo$miRID <- linage_type[match(geneInfo$gene_ID, linage_type$V1),]$V2
geneInfo <- geneInfo[,-ncol(geneInfo)]
geneInfo <- geneInfo[order(geneInfo$moduleColor),]
write.table(geneInfo, "26_Sep_WGCNA_unmerged_geneInfo.txt", sep="\t", quote=F, col.names=NA)
print("WGCNA_unmerged_geneInfo.txt ....... saved")

#6.MERGE CLOSE MODULES: Not do this in scRNA-seq. Choose different input gene list is much useful

#7.Plot tSNE gradients plot, assuming we got tSNE embedding from Seurat
tSNEembedding <- as.data.frame(df@dr$tsne@cell.embeddings)
exprMatrix<-ME1_mouse
predict_marker<-colnames(ME1_mouse)

pl <- list()
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
  pl <- append(pl, p)
  pdf(file = paste(geneName, "_tSNEplot.pdf", sep = ""),width = 15, height = 10, useDingbats = F)
  print(p)
  dev.off()
}
ml <- gridExtra::grid.arrange(grobs = pl, ncol = 5)
