### Date: 2017-11-12
### Project: Da Mi's project
### Goal: Analyze cortical interneurons

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

# setwd(paste("C:\\Users\\zhenadmin\\Data\\Tarik"))
setwd("~/Data/Da_Mi/")
#### 2. load data ####

## Only load RPKMs since that's the only type of data we will be using

mCounts <- read.table("./raw_data/mouseBrain.singleCell.gene.count.txt", header = T)
rownames(mCounts) <- mCounts$Geneid
mCounts <- mCounts[, -1]

QCmetrics <- read.csv("./raw_data/brainSeq.singleCell.QCmetrics.mouse.csv", 
                      stringsAsFactors = F)
QCmetrics$sampName[1:24] <- colnames(mCounts)[1:24]
colnames(mCounts) <- gsub("90_", "90.", colnames(mCounts))
colnames(mCounts) <- gsub("91_", "91.", colnames(mCounts))
colnames(mCounts) <- gsub("zhen_", "", colnames(mCounts))
QCmetrics <- QCmetrics[ match(colnames(mCounts), QCmetrics$sampName), ]
goodBarcodes <- read.csv("./raw_data/2017-02-20_mouse_singleCell_goodBarcodes.csv", 
                         stringsAsFactors = F)
mouse <- goodBarcodes$time_point %in% c("E14.5", "E17.5", "P1")
GB <- goodBarcodes[ mouse, ]
QC <- QCmetrics[ mouse, ]

Counts<- mCounts[ , mouse ]

## select single cells with 100,000 cds mapped reads (at least 30% of all reads)
## covering at least 1000 genes and less than 10% reads mapped to mitochondria genome

SC <- GB$cell_number == 1 & GB$read_depth >= 100000 & GB$coverage >= 1000 &
  QC$cdsRatio > 30 & QC$chrM < 10 & GB$time_point %in% c("E14.5", "E17.5", "P1")
GB <- GB[ SC, ] 
QC <- QC[ SC, ]
Counts <- Counts[ ,SC ]
# RPKM <- RPKM[ ,SC ]

QC$sampName <- colnames(Counts) <- GB$new_cell_ID
GB$IFC <- sapply(as.character(GB$new_cell_ID), function(x){strsplit(x, split = "\\.")[[1]][2]})

## Remove duplicated single cells
dp <- !duplicated(GB$new_cell_ID)

GB <- GB[dp,]
QC <- QC[dp,]
Counts <- Counts[ ,dp ]
#RPKM <- RPKM[ ,dp ]

RP <- rownames(Counts) %in% grep("RP", rownames(Counts), value = T)
mito <- rownames(Counts) %in% grep("MT-", rownames(Counts), value = T)

df0 <- CreateSeuratObject(raw.data = Counts[!(RP | mito),], min.cells = 5, min.genes = 1000, project = "Tarik")
df0 <- NormalizeData(object = df0, normalization.method = "LogNormalize", scale.factor = 100000)
df0@meta.data <- cbind(df0@meta.data, GB)

ggplot(df0@meta.data, aes(x = IFC, y = percent.mito, color = IFC)) +
  geom_boxplot()

#### 3. Only analyze cortical cells ####
df.cortex <- SubsetData(df0, ident.remove = c("13", "26", "58", "59", "114", "115"))
df.cortex <- ScaleData(df.cortex)

#### 4. plot heatmap ####
panRGC <- c("\\bNES$","\\bCD133$", "\\bGLI3$","\\bPAX6$","\\bPROM1$", "\\bSOX1$", "\\bSOX9$",
            "\\bSOX2$", "\\bVIM$")
aRGC<- c("\\bFGFR3$", "\\bFOXP2$","\\bHES1$", "\\bLRIG1$", "\\bNHLH2$","\\bPALLD$","\\bS100B$", 
         "\\bTAGLN2$")
bRGC <- c("\\bFABP7$", "\\bFAM107A$","\\bHOPX$","\\bLIFR$","\\bMLC1$","\\bPTPRZ1$", "\\bSLC1A3$","\\bTNC$","\\bUNC5D$")
IPC <- c("\\bEOMES$","\\bNEUROD1$","\\bNEUROD2$","\\bNEUROD4$","\\bNEUROD6$","\\bPENK$","\\bPPP1R17$")
EENR <- c("\\CUX1$","\\bCUX2$","\\bFEZF2$", "\\bSOX4$","\\bSOX11$")
pNeuron <- c("\\bDCX$","\\bDPYSL3$","\\bENO2$","\\bRBFOX3$","\\bSYP$")
EEN <- c("\\bNOS1$","\\bNOS2$","\\bNOS3$","\\bBCL11A$" ,"\\bBCL11B$", "\\bBHLHE22$","\\bCELF4$","\\bDPYSL5","\\bFOXP1$", 
         "\\bROBO2$","\\bSOX5$","\\bSTMN2$","\\bTBR1$")
LEN <- c("\\bMEF2C$","\\bPOU3F2$","\\bPOU3F3$","\\bRORB$","\\bSATB2$")
IN <- c("\\bDLX1$", "\\bDLX2$", "\\bDLX5$", "\\bDLX6$", "\\bGAD1$", "\\bGAD2$", "\\bISRL2$", 
        "\\bNPY$", "\\bPOU3F1", "\\bPVALB$", "\\bCCK$", "\\bHTR3A$","\\bVIP$", "\\bSST$", "\\bLHX6$","\\bA930038C07RIK$","\\bMAF$", "\\bMAFB$")
Astro <- c("\\bATP13A4$", "\\bALDH1L1$","\\bAQP4$", "\\bBMPR1B$", "\\bEDNRB$", "\\bEGFR$","\\bDIO2$","\\bFGFR3$",
           "\\bGFAP$","\\bGJA1$", "\\bGRM3$", "\\bLRIG1$", 
           "\\bMLC1$","\\bSLC4A4$", "\\bSOX9$")
Olig <- c("\\bNKX2-2$", "\\bOLIG1$","\\bOLIG2$", "\\bSDC3$", "\\bPID1$")
micro <- c("\\bCCL2$","\\bCCL3$","\\bCD83$","\\bCX3CR1$","\\bIL1A$","\\bRUNX1$",  
           "\\bTMEM119$", "\\bSALL1$")

genes <- unique(c(#panRGC, aRGC, bRGC, IPC, EENR, pNeuron,EEN,LEN,
                  IN))
                  #,Astro,Olig,micro))
genes <- unlist(sapply(genes,function(x){grep(x, rownames(df.int@data), value = T)}))

#### 5. Get correlated genes from Sofia's dataset ####
sofia <- read.csv("Sofia_Normalized_expression_values_of_mRNA.csv", header = T,row.names = 1)
grep("VIP", rownames(sofia), value = T)

pv.sofia <- sofia["ENSMUSG00000005716|PVALB", ]
pv.sofia.cor <- apply(sofia, 1, function(xx){cor(xx, as.numeric(pv.sofia))})
pv.sofia.cor <- sort(pv.sofia.cor, decreasing = T)
pv.sofia.cor <- pv.sofia.cor[pv.sofia.cor > 0.7]

sst.sofia <- sofia["ENSMUSG00000004366|SST", ]
sst.sofia.cor <- apply(sofia, 1, function(xx){cor(xx, as.numeric(sst.sofia),method = "s")})
sst.sofia.cor <- sort(sst.sofia.cor, decreasing = T)
sst.sofia.cor <- sst.sofia.cor[sst.sofia.cor > 0.7]

vip.sofia <- sofia["ENSMUSG00000019772|VIP", ]
vip.sofia.cor <- apply(sofia, 1, function(xx){cor(xx, as.numeric(vip.sofia),method = "s")})
vip.sofia.cor <- sort(vip.sofia.cor, decreasing = T)
vip.sofia.cor <- vip.sofia.cor[vip.sofia.cor > 0.6]

genes <- unique(c(names(pv.sofia.cor), names(sst.sofia.cor), names(vip.sofia.cor)))
genes <- intersect(genes, rownames(df.int@data))

## Get go terms for PC correlated genes
library("RDAVIDWebService") ##--Load RDAVIDWebService 
user <- DAVIDWebService$new(email = "zhen.li.zl242@yale.edu", 
                            url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")     ##--Log on to DAVID with email

david <- substr(genes, 20, 100)

##--Get a list of genes, pay attention to idType
##--Read a gene list
##--submit gene list to DAVID
addList(user, inputIds = david, idType = "ENSEMBL_GENE_ID", listType = "Gene", 
        listName = "Da_Mi")

##--check a specific category
# setAnnotationCategories(user, "KEGG_PATHWAY")                                  
# 
# ##--combine GO terms from each list
# term <- getFunctionalAnnotationChart(user)
# KEGG[[i]] <- term$Term

setAnnotationCategories(user,categories =  c("GOTERM_BP_FAT", "KEGG_PATHWAY"))
term <- getFunctionalAnnotationChart(user)
goterm <- data.frame(matrix(unlist(term), nrow = length(term$Term)),stringsAsFactors = F)
colnames(goterm) <- names(term)
BP[[ii]] <- term$Term

goterm <- goterm[order(as.numeric(goterm$Count), decreasing = T), ]

g <- goterm[8, "Genes"]
g <- unlist(strsplit(g,split = ", "))

genes <- genes[david %in% g]
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

MZSVZ <- read.csv("MZ vs SVZ gene list 08-11-2017.csv", header = F, stringsAsFactors = F)
MZSVZ[,1] <- toupper(MZSVZ[,1])
MZSVZ[,1] <- unique(sapply(MZSVZ[,1], function(x){paste("\\b", x, "$", sep = "")}))
MZSVZ[,1] <- as.character(sapply(MZSVZ[,1],function(x){grep(x, rownames(df.int@data), value = T)}))
MZSVZ <- MZSVZ[MZSVZ[,1] != "character(0)",]

temp <- df.int@scale.data[ c(as.character(genes),MZSVZ[,1])[-116], ]
temp <- df.int@scale.data[ as.character(genes), ]

dis <- 1 - cor(as.matrix(temp), use = "pairwise.complete.obs")
hca <- hclust(as.dist(dis),method = "average")
ct <- cutree(hca, h = 0.875)

ColSideColors <- c(#brewer.pal(length(unique(ct)), "Paired")[ct])#,
                   gg_color_hue(length(unique(df.int@meta.data$time_point)))[as.numeric(as.factor(df.int@meta.data$time_point))])#,
                   # rainbow(n)[as.numeric(x)])
# 
ColSideColors <- matrix(data = ColSideColors, nrow = ncol(temp));
# ColSideColors <- ColSideColors[ord, ];

# colsep <- lapply(unique(x[ord]), function(xx){length(which(x == xx))});
# colsep <- cumsum(unlist(colsep));

pairs.breaks <- seq(-2.5, 2.5, length.out=101);

pdf(file = paste(Sys.Date(), "_int_marker_heatmap_clusters_temp.pdf", sep = ""), 
    width = 20, height = 40, useDingbats = F);
heatmap.3(temp,
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8,
          main="Marker gene expression",
          col = bluered(100),
          symkey = F,
          cexRow=1, cexCol = 0.6, 
          Rowv = F,
          Colv = T,#as.dendrogram(hca),
          ColSideColors = as.matrix(ColSideColors),
          ColSideColorsSize = 1,
          #RowSideColors = RowSideColors,
          dendrogram = "column",
          #scale = "row",
          # colsep = colsep,
          # sepcolor = "black",
          labRow = substr(rownames(temp), start = 20, stop = 100), #c("Glast", "Tbr2"),
          labCol = "",
          na.rm = F);
dev.off(dev.cur())

df.cortex <- AddMetaData(df.cortex, ct, col.name = "cutree")
df.cortex <- SetAllIdent(df.cortex, "cutree")
df.int <- SubsetData(df.cortex, ident.use = 3)

#### 6. Analyze cortical interneuron ####
df.int <- FindVariableGenes(object = df.int, 
                           mean.function = ExpMean, 
                           dispersion.function = LogVMR, 
                           x.low.cutoff = 0.35, 
                           x.high.cutoff = Inf, 
                           y.cutoff = 0.5,
                           y.high.cutoff = Inf,
                           do.plot = T)

lin <- read.csv("./raw_data/INs lineage marker 12_11-17.csv", header = F, stringsAsFactors = F)
lin[,1] <- sapply(lin[,1], function(x){paste("\\b", x, "$", sep = "")})
lin[,1] <- as.character(sapply(lin[,1],function(x){grep(x, rownames(df.int@scale.data), value = T)}))
lin <- lin[!(lin$V1 == "character(0)"),]

set.seed(1); df.int <- RunPCA(object = df.int, 
                             pc.genes = genes, 
                             pcs.compute = 15, 
                             pcs.print = 1:10, 
                             genes.print = 10, 
                             do.print = T)

df.int <-JackStraw(df.int, num.replicate = 100, do.print = F,num.pc = 15)

library(Rtsne)
tsne <- Rtsne(t(df.int@scale.data[genes,]),pca = T,
              perplexity = 5)$Y
colnames(tsne) <- c("tSNE_1", "tSNE_2")
rownames(tsne) <- df.int@cell.names
df.int <- SetDimReduction(object = df.int, reduction.type = "tsne", 
                          slot = "cell.embeddings", new.data = tsne)


pdf(file = paste(Sys.Date(), "cortical_int_JackStraw_temp.pdf", sep = "_"), width = 11,
    height = 10, useDingbats = F);
JackStrawPlot(df.int, PCs = 1:15)
dev.off()

df.int <- RunTSNE(object = df.int, dims.use = 1:5, do.fast = TRUE, 
                  perplexity = 10, seed.use = 1)
df.int <- FindClusters(object = df.int, 
                      reduction.type = "pca", 
                      dims.use = 1:5, 
                      resolution = 2, # dMGE, CGE = 0.8, E14.5 vMGE = 2
                      k.scale = 25, # dMGE, CGE = 50
                      prune.SNN = 0, #dMGE, CGE = 0.5
                      plot.SNN = F, 
                      print.output = F, 
                      save.SNN = F,
                      algorithm = 1,
                      force.recalc = TRUE, 
                      random.seed = 1,
                      k.param = 10) 

p1 <- TSNEPlot(object = df.int, pt.size = 3, do.return = T)
p2 <-TSNEPlot(object = df.int, group.by = "time_point", do.return = T,
              pt.size = 3)
p3 <-TSNEPlot(object = df.int, group.by = "IFC", do.return = T,
              pt.size = 3)

pdf(file = paste(Sys.Date(), "cortical_int_tSNE_temp.pdf", sep = "_"), width = 33, 
    height = 10, useDingbats = F);
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

markers <- FindAllMarkers(object = df.int, only.pos = TRUE, min.pct = 0.1, 
                          return.thresh = 0.01)
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_diff)

pdf(file = paste(Sys.Date(), "cortical_int_clusters_heatmap_top20.pdf", sep = "_"),
    width = 20, height = 40, useDingbats = F);
DoHeatmap(object = df.int,
          genes.use = top20$gene,#rownames(df.markers)[1:100], #genes,
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

pdf(file = paste(Sys.Date(), "cortical_int_features_temp.pdf", sep = "_"),
    width = 100, height = 70, useDingbats = F);
FeaturePlot(df.int, features.plot = genes,cols.use = c("grey", "blue"),
            nCol = 9,pt.size = 10)
dev.off()

temp <- data.frame(t(df.int@scale.data[top10$gene,]))
colnames(temp) <- substr(colnames(temp), start = 20, stop = 100)
temp$clusters <- df.int@meta.data$res.2
temp <- melt(temp, id.vars = "clusters")

pdf(file = paste(Sys.Date(), "cortical_int_features_temp.pdf", sep = "_"),
    width = 5, height = 35, useDingbats = F);
ggplot(temp, aes(x = clusters, y = value, fill = clusters, color = clusters)) +
  geom_violin(scale = "width") + 
  facet_wrap(~variable, nrow = length(top10$gene),
             strip.position = "right",scales = "free") +
  theme(axis.text = element_text(size = 4),
        axis.line.x = element_line(size = 0.1),
        axis.line.y = element_line(size = 0.1),
        axis.ticks.x = element_line(size = 0.1),
        axis.ticks.y = element_line(size = 0.1),
        strip.text.y = element_text(size = rel(0.5), angle = 0),
        strip.background = NULL)
dev.off()

#### 7. Identify Sst, PV, Vip cells ####
SST <- df.int@scale.data["ENSMUSG00000004366|SST", ]
VIP <- df.int@scale.data["ENSMUSG00000019772|VIP", ]
LHX <- df.int@scale.data["ENSMUSG00000026890|LHX6", ]
NDNF <- df.int@scale.data["ENSMUSG00000049001|A930038C07RIK", ]

df.int@meta.data$cell_type <- "uncertain" 
df.int@meta.data$cell_type[SST > 0 & VIP <= 0 & NDNF <= 0 & LHX > 0] <- "SST"
df.int@meta.data$cell_type[SST <= 0 & VIP > 0 & NDNF <= 0 & LHX <= 0 ] <- "VIP"
df.int@meta.data$cell_type[SST <= 0 & VIP <= 0 & NDNF <= 0 & LHX > 0 ] <- "PV"

df.int <- SetAllIdent(df.int, id = "cell_type")
df.int.3types <- SubsetData(df.int, ident.remove = "uncertain") 

df.int.3types <- FindVariableGenes(object = df.int.3types, 
                            mean.function = ExpMean, 
                            dispersion.function = LogVMR, 
                            x.low.cutoff = 0.35, 
                            x.high.cutoff = Inf, 
                            y.cutoff = 0.5,
                            y.high.cutoff = Inf,
                            do.plot = T)

temp <- df.int@scale.data[ as.character(genes), ]

dis <- 1 - cor(as.matrix(temp), use = "pairwise.complete.obs")
hca <- hclust(as.dist(dis),method = "average")
ct <- cutree(hca, h = 0.875)

ColSideColors <- c(brewer.pal(3, "Paired")[as.numeric(as.factor(df.int.3types@meta.data$cell_type))],
  gg_color_hue(length(unique(df.int.3types@meta.data$time_point)))[as.numeric(as.factor(df.int.3types@meta.data$time_point))])#,
# rainbow(n)[as.numeric(x)])
# 
ColSideColors <- matrix(data = ColSideColors, nrow = ncol(temp));
# ColSideColors <- ColSideColors[ord, ];

# colsep <- lapply(unique(x[ord]), function(xx){length(which(x == xx))});
# colsep <- cumsum(unlist(colsep));

pairs.breaks <- seq(-2.5, 2.5, length.out=101);

pdf(file = paste(Sys.Date(), "_int_3types_marker_heatmap_clusters_temp.pdf", sep = ""), 
    width = 20, height = 40, useDingbats = F);
heatmap.3(temp,
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8,
          main="Marker gene expression",
          col = bluered(100),
          symkey = F,
          cexRow=1, cexCol = 0.6, 
          Rowv = F,
          Colv = T,#as.dendrogram(hca),
          ColSideColors = as.matrix(ColSideColors),
          ColSideColorsSize = 1,
          #RowSideColors = RowSideColors,
          dendrogram = "column",
          #scale = "row",
          # colsep = colsep,
          # sepcolor = "black",
          labRow = substr(rownames(temp), start = 20, stop = 100), #c("Glast", "Tbr2"),
          labCol = "",
          na.rm = F);
dev.off(dev.cur())

#### 7.2 PCA and tSNE ####
df.int.3types <- FindVariableGenes(object = df.int.3types, 
                                       mean.function = ExpMean, 
                                       dispersion.function = LogVMR, 
                                       x.low.cutoff = 0.3, 
                                       x.high.cutoff = 6, 
                                       y.cutoff = 0.5, y.high.cutoff = 5, 
                                       do.plot = T)
abline(v = 0.3, col = "red")
df.int.3types <- RunPCA(df.int.3types,pc.genes = markers.temp$gene,
                        pcs.compute = 15, 
                        pcs.print = 1:10,
                        genes.print = 10)

df.int.3types <-JackStraw(df.int.3types, num.replicate = 100, 
                          do.print = T, num.pc = 15)
JackStrawPlot(df.int.3types, PCs = 1:15)

df.int.3types <- RunTSNE(df.int.3types, dims.use = c(1:2),#pcs[[i]], 
                             perplexity = 5,
                             seed.use = 3)

p1 <- TSNEPlot(object = df.int.3types, group.by = "cell_type", pt.size = 3, do.return = T)
p2 <- TSNEPlot(object = df.int.3types, group.by = "time_point", do.return = T,
              pt.size = 3)
# p3 <-TSNEPlot(object = df0.E12.5.prog.vz, group.by = "", do.return = T,
#               pt.size = 3)
pdf(file = paste(Sys.Date(),"cortex_int_3types_tSNE.pdf", sep = "_"), width = 22, 
    height = 10, useDingbats = F);
plot_grid(p1, p2, 
          # p3, 
          ncol = 2)
dev.off()

markers.3types <- FindAllMarkers(df.int.3types,
                               min.pct = 0.5,
                               return.thresh = 0.01,
                               only.pos = T)

sst.pv <- FindMarkers(df.int.3types, ident.1 = "SST", ident.2 = "PV",
            min.pct = 0.5, only.pos = F)
sst.pv$genes <- rownames(sst.pv)
sst.pv.top40 <- sst.pv %>% top_n(40, abs(avg_diff))
sst.pv.top40 <- sst.pv.top40[order(sst.pv.top40$avg_diff),]

top40 <- markers.temp %>% group_by(cluster) %>% top_n(40, avg_diff)

pdf(file = paste(Sys.Date(),"cortex_3types_DE_heatmap_sst.pv.pdf", sep = "_"),
    width = 20, height = 20, useDingbats = F);
DoHeatmap(object = df.int.3types,
          genes.use = sst.pv.top40$genes,#top40$gene,#rownames(df.markers)[1:100],
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

#### 7.3 Use Random Forest to reclassify VZ and SVZ cells ####
cate_bf_fs <- as.factor(df.int.3types@ident)
feature_bf_fs <- as.matrix(t(df.int.3types@scale.data[df.int.3types@var.genes,]))

rf_bf_fs <- randomForest(feature_bf_fs, cate_bf_fs, importance = TRUE, proximity = TRUE)
imp_bf_fs <- importance(rf_bf_fs, type = 1)

fs <- rfcv(feature_bf_fs, cate_bf_fs, cv.fold = 10, scale = "log", step = 0.9, recursive = F)
len <- length(fs$n.var[fs$error.cv == min(fs$error.cv)])
min_fs <- fs$n.var[fs$error.cv == min(fs$error.cv)][1] #get least features
ind <- order(-imp_bf_fs)[1: min_fs]
feature_fs <- feature_bf_fs[ , ind]
cate_fs <- cate_bf_fs

temp <- t(feature_fs)

ord <- order(cate_fs)
ColSideColors <- c(#gg_color_hue(3)[as.numeric(as.factor(df0.E12.5.prog@meta.data$region))],
  brewer.pal(3, name = "Paired")[as.numeric(cate_fs)])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- as.matrix(ColSideColors[ord,])
col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

pdf(file = paste(Sys.Date(), "cortex_int_heatmap_RF.pdf", sep = "_"),
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

rf_fs <- randomForest(feature_fs, as.factor(cate_fs), importance=TRUE, proximity=TRUE)
fea1_fs <- data.frame()
fea1_fs <- feature_fs[(rf_fs$predicted == 'PV') & (rf_fs$votes[ , 1] > 0.6), , drop = FALSE]
cat1_fs <- rf_fs$predicted[(rf_fs$predicted =='PV') & (rf_fs$votes[ , 1] > 0.6)]
fea2_fs <- data.frame()
fea2_fs <- feature_fs[(rf_fs$predicted == 'SST') & (rf_fs$votes[ , 2] > 0.6), , drop = FALSE]
cat2_fs <- rf_fs$predicted[(rf_fs$predicted =='SST') & (rf_fs$votes[ , 2] > 0.6)]
fea3_fs <- data.frame()
fea3_fs <- feature_fs[(rf_fs$predicted == 'VIP') & (rf_fs$votes[ , 3] > 0.6), , drop = FALSE]
cat3_fs <- rf_fs$predicted[(rf_fs$predicted =='VIP') & (rf_fs$votes[ , 3] > 0.6)]

cate <- as.factor(c(as.character(cat1_fs), as.character(cat2_fs),as.character(cat3_fs)))
feature <- as.matrix(rbind(fea1_fs, fea2_fs, fea3_fs))

set <- sample(1: nrow(feature), nrow(feature), replace = F)
cate <- cate[set]
feature <- feature[set, ] 

rf_whole <- randomForest(feature, as.factor(cate), importance = TRUE, proximity = TRUE)
pred_whole <- predict(rf_whole, newdata = feature_fs)

temp <- t(feature_fs)
# dis <- as.dist(1 - cor(temp,use = "pairwise.complete.obs"))
# clu <- hclust(dis, method = "complete")
# cl <- cutree(clu, k=2)

ord <- order(pred_whole)
ColSideColors <- c(#gg_color_hue(3)[as.numeric(as.factor(df0.E12.5.prog@meta.data$region))],
  brewer.pal(3, name = "Paired")[as.numeric(pred_whole)])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- as.matrix(ColSideColors[ord,])
col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

pdf(file = paste(Sys.Date(), "cortex_int_heatmap_post-RF.pdf", sep = "_"),
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

markers.3types.RF <- colnames(feature)

#### 7.3.2 Use Random Forest result to predict all interneuron  from cortex####
fea_int <- t(df.int@scale.data[markers.3types.RF,])
pred_whole <- predict(rf_whole, newdata = fea_int,type = "prob")
pred <- apply(pred_whole, 1, function(xx){
  if(max(xx) >= 0.5){
    c("PV", "SST", "VIP")[xx == max(xx)]
  } else {"uncertain"}
})

df.int <- AddMetaData(df.int,metadata = pred,col.name = "RF_cell_type")

#### 7.3.3 Diffusion Map and tSNE plot ####
library(diffusionMap)
library(Rtsne)
# dis <- as.dist(1 - cor(df.int@scale.data[markers.3types.RF,],
#                        use = "pairwise.complete.obs",method = "s"))
dis <- dist(t(df.int@scale.data[markers.3types.RF,]),method = "manhattan")
dm <- diffuse(dis, neigen = 2,eps.val = 320, delta = 10^-5)
dm <- data.frame(dm$X)
colnames(dm) <- c("dm_1", "dm_2")

pdf(file = paste(Sys.Date(), "cortex_int_DM_post-RF.pdf", sep = "_"),
    width = 11, height = 10, useDingbats = F);
ggplot(dm, aes(x = dm_1, y = dm_2, 
               color = df.int@meta.data$RF_cell_type)) +
  geom_point() +
  scale_color_manual(name = "Cell Type",
                     values = c(gg_color_hue(3)[1:2], "grey90", 
                                gg_color_hue(3)[3])) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8))
dev.off()

#### 7.4 Use Random Forest result to predict all interneuron from Da Mi data ####
fea_int <- t(df0.neu@scale.data[intersect(markers.3types.RF, rownames(df0.neu@scale.data)),])
fea_int <- cbind(fea_int, 0)
colnames(fea_int)[17] <- "ENSMUSG00000019772|VIP"
pred_whole <- predict(rf_whole, newdata = fea_int,type = "prob")
pred <- apply(pred_whole, 1, function(xx){
  if(max(xx) >= 0.5){
    c("PV", "SST", "VIP")[xx == max(xx)]
  } else {"uncertain"}
})

df0.neu <- AddMetaData(df0.neu,metadata = pred,col.name = "RF_cell_type")
df0.neu <- SetAllIdent(df0.neu, id = "RF_cell_type")
df0.neu.temp <- SubsetData(df0.neu, ident.remove = "uncertain")

#### 7.3 Use Random Forest to reclassify VZ and SVZ cells ####
cate_bf_fs <- as.factor(as.character(df0.neu.temp@meta.data$RF_cell_type))
feature_bf_fs <- as.matrix(t(df0.neu.temp@scale.data[df0.neu.temp@var.genes,]))

rf_bf_fs <- randomForest(feature_bf_fs, cate_bf_fs, importance = TRUE, proximity = TRUE)
imp_bf_fs <- importance(rf_bf_fs, type = 1)

fs <- rfcv(feature_bf_fs, cate_bf_fs, cv.fold = 10, scale = "log", step = 0.9, recursive = F)
len <- length(fs$n.var[fs$error.cv == min(fs$error.cv)])
min_fs <- fs$n.var[50] #get least features
ind <- order(-imp_bf_fs)[1: min_fs]
feature_fs <- feature_bf_fs[ , ind]
cate_fs <- cate_bf_fs

temp <- t(fea_int)

ord <- order(pred)
ColSideColors <- c(#gg_color_hue(3)[as.numeric(as.factor(df0.E12.5.prog@meta.data$region))],
  brewer.pal(4, name = "Paired")[as.numeric(as.factor(pred))])#,
ColSideColors <- matrix(ColSideColors, nrow = ncol(temp))
ColSideColors <- as.matrix(ColSideColors[ord,])
col <- colorRampPalette(colors = c("#1b75bb", "black","#faaf40"))(1000)

pairs.breaks <- seq(-1.5, 1.5, length.out=1001);

pdf(file = paste(Sys.Date(), "all_Da_Mi_int_heatmap_post-RF.pdf", sep = "_"),
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

rf_fs <- randomForest(feature_fs, as.factor(cate_fs), importance=TRUE, proximity=TRUE)
fea1_fs <- data.frame()
fea1_fs <- feature_fs[(rf_fs$predicted == 'PV') & (rf_fs$votes[ , 1] > 0.6), , drop = FALSE]
cat1_fs <- rf_fs$predicted[(rf_fs$predicted =='PV') & (rf_fs$votes[ , 1] > 0.6)]
fea2_fs <- data.frame()
fea2_fs <- feature_fs[(rf_fs$predicted == 'SST') & (rf_fs$votes[ , 2] > 0.6), , drop = FALSE]
cat2_fs <- rf_fs$predicted[(rf_fs$predicted =='SST') & (rf_fs$votes[ , 2] > 0.6)]
fea3_fs <- data.frame()
fea3_fs <- feature_fs[(rf_fs$predicted == 'VIP') & (rf_fs$votes[ , 3] > 0.6), , drop = FALSE]
cat3_fs <- rf_fs$predicted[(rf_fs$predicted =='VIP') & (rf_fs$votes[ , 3] > 0.6)]

cate <- as.factor(c(as.character(cat1_fs), as.character(cat2_fs),as.character(cat3_fs)))
feature <- as.matrix(rbind(fea1_fs, fea2_fs, fea3_fs))

set <- sample(1: nrow(feature), nrow(feature), replace = F)
cate <- cate[set]
feature <- feature[set, ] 

fea_int <- t(df0.neu@scale.data[c(colnames(feature),
                                  "ENSMUSG00000074622|MAFB",
                                  "ENSMUSG00000055435|MAF",
                                  "ENSMUSG00000028222|CALB1",
                                  "ENSMUSG00000008393|CARHSP1",
                                  "ENSMUSG00000004151|ETV1",
                                  "ENSMUSG00000030199|ETV6",
                                  "ENSMUSG00000038679|TRPS1"),])

rf_whole <- randomForest(feature, as.factor(cate), importance = TRUE, proximity = TRUE)
pred_whole <- predict(rf_whole, newdata = fea_int, type = "prob")

pred <- apply(pred_whole, 1, function(xx){
  if(max(xx) >= 0.9){
    c("PV", "SST", "VIP")[xx == max(xx)]
  } else {"uncertain"}
})

df0.neu <- AddMetaData(df0.neu, metadata = pred, col.name = "RF_cell_type")

dis <- dist(fea_int,method = "manhattan")
dm <- diffuse(dis, neigen = 2, eps.val = 600, delta = 10^-40)
dm <- data.frame(dm$X)
colnames(dm) <- c("dm_1", "dm_2")

pdf(file = paste(Sys.Date(), "all_int_Da_Mi_DM_post-RF.pdf", sep = "_"),
     width = 11, height = 10, useDingbats = F);
ggplot(dm, aes(x = dm_1, y = dm_2, 
               color = df0.neu@meta.data$RF_cell_type)) +
  geom_point() +
  scale_color_manual(name = "Cell Type",
                     values = c(gg_color_hue(3)[1:2], "grey90", 
                                gg_color_hue(3)[3])) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8))
dev.off()

df0.neu.temp <- SubsetData(df0.neu,ident.remove = "uncertain")
markers.temp <- FindAllMarkers(df0.neu.temp,
                               only.pos = T,
                               return.thresh = 1,
                               min.pct = 0.2,
                               genes.use = c(colnames(feature),
                                                          "ENSMUSG00000074622|MAFB",
                                                          "ENSMUSG00000055435|MAF",
                                                          "ENSMUSG00000028222|CALB1",
                                                          "ENSMUSG00000008393|CARHSP1",
                                                          "ENSMUSG00000004151|ETV1",
                                                          "ENSMUSG00000030199|ETV6",
                                                          "ENSMUSG00000038679|TRPS1"))

markers.temp <- FindAllMarkers(df0.neu.temp,
                               only.pos = T,
                               return.thresh = 1,
                               min.pct = 0.1)

top40 <- markers.temp %>% group_by(cluster) %>% top_n(40, avg_diff)

pdf(file = paste(Sys.Date(),"all__heatmap_top40.pdf", sep = "_"),
    width = 20, height = 40, useDingbats = F);
DoHeatmap(object = df0.neu.temp,
          genes.use = top40$gene,#rownames(dff.markers)[1:100],
          slim.col.label = TRUE,
          group.spacing = 0.3,
          remove.key = TRUE)
dev.off()


#### 
df.all.int <- RunCCA(object = df0.neu, object2 = df.int.3types, genes.use = markers.3types$gene)

DimPlot(df.all.int, reduction.use = "cca", group.by = "time_point")
DimHeatmap(df.all.int, reduction.type = "cca", dim.use = 1:10)

df.all.int  <- AlignSubspace(df.all.int, reduction.type = "cca", 
                             grouping.var = "orig.ident", dims.align = 1:9)
df.all.int <- RunPCA(df.all.int,pc.genes = markers.3types$gene)
PCAPlot(df.all.int,group.by = "cell_type")

dis <- as.dist(1 - cor(t(df.all.int@dr$cca.aligned@cell.embeddings[,1:9])))
dm <- diffuse(dis, neigen = 2)
dm <- data.frame(dm$X)
colnames(dm) <- c("dm_1", "dm_2")

ggplot(dm, aes(x = dm_1, y = dm_2, color = df.all.int@meta.data$cell_type)) +
  geom_point() +
  theme_bw()
