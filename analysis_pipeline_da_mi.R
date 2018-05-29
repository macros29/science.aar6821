### Date: 2016-02-29
### Project: Human Brain Development and Evolution
### Goal: Analyze single cells from human brain and iPSC together

#### 1. Load libraries ####c
library("matrixStats")
library("cluster")
library("SCnorm")
library("reshape")
library("Rtsne")
library("gplots")
library("ggplot2")
library("pcaReduce")
library("parallel")
library("locfit")
library("MASS")
library("hydroGOF")
library("calibrate")
library("GGally")
library("monocle")
library("randomForest")
source("~/Scripts/R/heatmap.3.R")
source("~/Scripts/R/SNN.R");
source("~/Scripts/R/SNN-Cliq-Gap.R")
source("~/Scripts/R/Waterfall.R")
source("~/Scripts/R/gg_color_hue.R")
source("~/Scripts/R/spec.R");
options(stringAsFactors = F);

#### 2. Setup working directory and load data ####

## 2.1 load data
#setwd("/Users/zhenadmin/Da_Mi/")
setwd("~/Data/Da_Mi/")
load("working.RData")

rawCounts <- read.table("./raw_data/mouseBrain.Mida.singleCell.gene.count.txt")
#counts <- read.table("./raw_data/mouseBrain.Mida.scRNAseq_cleanData.gene.count.txt", header = T, stringsAsFactors = F)
#cpm <- read.table("./raw_data/mouseBrain.Mida.scRNAseq_cleanData.gene.cpm.txt", header = T)
RPKM <- read.table("./raw_data/mouseBrain.Mida.scRNAseq_cleanData.gene.RPKM.txt", header = T)

df <- RPKM
#df <- cpm
rownames(df) <- df$Geneid
df <- df[, -1]

## 2.2 adjust sample names to be coherent with goodBarcodes
colnames(counts) <- gsub("90_", "90.", colnames(counts))
colnames(counts) <- gsub("91_", "91.", colnames(counts))
colnames(counts) = unlist(lapply(colnames(counts), function(x){strsplit(x, split = "_")[[1]][1]}))

colnames(df) <- gsub("90_", "90.", colnames(df))
colnames(df) <- gsub("91_", "91.", colnames(df))
colnames(df) = unlist(lapply(colnames(df), function(x){strsplit(x, split = "_")[[1]][1]}))

fil <- read.csv("raw_data/expr_matrix_filtering_readCounts.csv", header = T)

GB <- read.csv("./raw_data/new_goodBarcodes.csv", header = T)
GB$X <- gsub("\\_", "\\.", GB$X)

GB <- GB[ GB$X %in% colnames(df), ]
colnames(GB) <- c("cell_ID", "cell_number")
GB$IFC <- sapply(as.character(GB$cell_ID), function(x){strsplit(x, split = "\\.")[[1]][2]})
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

GB$year <- "2017"
GB$year[GB$IFC %in% c(90:101)] <- "2016"

QC <-  read.csv("./raw_data/oscar.scRNAseq.goodCell.QCmetrics.csv", header = T, stringsAsFactors = F)
QC$sampName <- gsub("90_", "90-", QC$sampName)
QC$sampName <- gsub("91_", "91-", QC$sampName)
rownames(QC) <- unlist(lapply(QC$sampName, function(x){strsplit(x, split = "_")[[1]][1]}))
rownames(QC) <- gsub("\\-", "\\.", rownames(QC))

GB <- GB[!(GB$IFC %in% c("189","180","181","190","191","192","193","194","195","196")),]
QC <- QC[GB$cell_ID, ]
df <- df[, GB$cell_ID]

counts <- counts[, GB$cell_ID]

#### Delete genes
expr <- df # input gene expression matrix 
genes<-rownames(expr)
gene_ID <- unlist(lapply(genes, function(x){strsplit(as.character(x), split = "\\|")[[1]][1]}))

library("biomaRt")
ensembl = useMart("ensembl", dataset="mmusculus_gene_ensembl")

test <- gene_ID
ensembl_gene <-getBM(attributes=c("ensembl_gene_id","external_gene_name","chromosome_name","transcript_biotype"),
                      filters = "ensembl_gene_id", values=test, mart=ensembl) #28590

GM <- grep("\\bGM", rownames(expr), value = T)
RP <- grep("\\bRP", rownames(expr), value = T)
keep <- ensembl_gene$transcript_biotype %in% c("protein_coding") & 
  !(ensembl_gene$chromosome_name %in% c("Y", "MT")) 
keep <- test %in% ensembl_gene$ensembl_gene_id[keep] & !(rownames(expr) %in% c(GM, RP))
expr <- expr[keep, ]

#### Try diffusion map ####
library(destiny)

dm <- DiffusionMap(t(Y2), verbose = T)

pdf(file = paste("./images/", Sys.Date(), "_dm_", p, ".pdf", 
                 sep = ""), width = 12, height = 10, useDingbats = F)
ggplot(data = dm, aes(x = DC1, y = DC3)) +
  geom_point(size = 3, aes(x = DC1, y = DC3, 
                 color = Normalized_GB$mouse)) +
  #scale_color_brewer(palette = "Paired", name = 'Ages') +
  scale_color_discrete(name = "Mouse") +
  #geom_text(aes(label = colnames(Y2)), hjust = 1, vjust = -1) + 
  labs(x = "DC1", y = "DC3") +
  ggtitle(paste("All cells, Perplexity =", p)) +
  theme_bw() +
  theme(text = element_text(size = 20),
        line = element_line(size = 1),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.key = element_blank(),
        legend.background = element_blank())
dev.off(which = dev.cur())

## Clean up
rm(expr, ensembl_gene, GM, RP, keep, test, ensembl)

#### 3. Use scVEGs to determine significant variably expressed genes ####
## Based on Chen et al., 2016. BMC Genomics. 2016; 17(Suppl 7): 508.
data <- Y
m <- dim(data)[1]
std <- apply(data, 1, sd)
avg <- apply(data, 1, mean)
cv <- std / avg

xdata <- log10(avg)
ydata <- log10(cv)
xdata <- xdata[is.na(ydata) != "TRUE"]
ydata <- ydata[is.na(ydata) != "TRUE"]

fitLoc <- locfit.robust(ydata ~ lp(xdata, nn = .2))

xSeq <- seq( min(xdata), max(xdata), 0.005 )
gapNum <- matrix( 0, length(xSeq), 1 )
for(i in 1:length(xSeq)) {
  cdx <- which((xdata >= xSeq[i] - 0.05) & (xdata < xSeq[i] + 0.05))
  gapNum[i,1] <- length(cdx)
}
cdx <- which(gapNum > m*0.005)
ySeq <- predict(fitLoc,xSeq)
ix <- length(ySeq) - 1
xSeq_all <- seq(min(xdata), max(xdata), 0.001)
xSeq <- xSeq[cdx[1]:ix[1] + 1]
ySeq <- ySeq[cdx[1]:ix[1] + 1]

b <- 1
a <- 0
df <- data.frame(x=xSeq, y = ySeq)
fit = nls(y ~ 0.5 * log10(b / 10^x + a), data = df, 
          start=list(b = b,a = a), nls.control(maxiter = 500), na.action =  'na.exclude')
newdf <- data.frame(x = xdata)
ydataFit <- predict(fit,newdata = newdf)
all_df <- data.frame(x = xSeq_all)
ydataAll <- predict(fit,newdata = all_df)

## Calculate CV difference
logX <- xdata
logXseq <- xSeq_all
cvDist <- matrix(abs(ydataFit - ydata),length(xdata), 1)

## use kernel density estimate to find the peak
dor <- density(cvDist, kernel = "gaussian")
distMid <-dor$x[which.max(dor$y)]
dist2 <- cvDist - distMid
tmpDist <- c(dist2[dist2 <= 0], abs(dist2[dist2 < 0])) + distMid
distFit <- fitdistr(tmpDist, "normal")
pRaw <- pnorm(cvDist, mean = distFit$estimate[1], sd = distFit$estimate[2], lower.tail = FALSE)
pVal <- 0.05
dx <- which(pRaw < pVal & (ydataFit - ydata) < 0)

pdf(file = paste("./images/", Sys.Date(), "_mouse_VEGs.pdf", sep = ""), 
    width = 10, height = 10, useDingbats = F);
plot(xdata, ydata, type = 'p', pch = 20, cex = 0.5)
abline(0, -0.5)
lines(xSeq, ySeq, col = 'cyan3')
lines(xSeq_all, ydataAll, col = "firebrick")
points(xdata[dx], ydata[dx], type = 'p', pch = 20, cex = 0.5, col = 'palegreen3')
textxy(max(xdata)-1, max(ydata, na.rm = TRUE) - 0.1,
       paste('alpha = ', round(coef(fit)[2], digit = 3), sep=''), cex = 1, offset = 0)
textxy(max(xdata)-1, max(ydata, na.rm = TRUE) - 0.2,
       paste('beta = ', round(coef(fit)[1], digit = 3), sep=''), cex = 1, offset = 0)
textxy(max(xdata)-1, max(ydata, na.rm = TRUE) - 0.3,
       paste('variable genes = ', length(dx), sep=''), cex = 1, offset = 0)
#textxy(xdata[2853], ydata[2853],
#       substr(rownames(Y)[2853], 17, 100), cex =0.5, offset = 0, col = "red")
dev.off()

Y2 <- Y[dx, ]

## Clean up
rm(a, b, i, m, fil, data, std, avg, cv, xdata, ydata, fitLoc, xSeq, gapNum, xSeq_all, ySeq, 
   all_df, df, fit, newdf, logX, logXseq, ydataAll, ydataFit, ix, cdx, cvDist, dor, 
   distMid, dist2, tmpDist, distFit, pRaw, pAdj, pVal)

#temp <- data.frame(t(Y[hk,]))
temp <- data.frame(t(Y[markers,]))
colnames(temp) <- substr(colnames(temp), 20, 100)
temp$IFC <- GB$IFC
temp$Time_point <- GB$Time_point
temp$batch <- GB$batch

temp1 <- melt(temp)
#CD <- "D14"
ggplot(temp1, aes(x = IFC, y = value, colour = batch)) + 
  geom_boxplot(outlier.size = 0) +
  geom_jitter() +
  #ggtitle(paste("Chimp")) +
  facet_wrap(~variable, nrow = 5) +
  theme_bw() +
  theme(plot.title = element_text(size=20, face="bold", vjust=2),
        #legend.position = c(0.95, 0.15),
        axis.text.x = element_blank())

#### OR use log-linear fit for mean-CV2 relation to identify high expression, high variable genes####
## log transform and quantile normalization ##
data <- df
#data <- NormData
RP <- grep("\\bRP", rownames(data), value = T)
MT <- grep("\\bMT", rownames(data), value = T)
Y <- log2(data[!(rownames(data) %in% c(RP, MT)), ] - min(data) + 1)
#Y <- log2(combat_edata[!(rownames(combat_edata) %in% c(RP, MT)), ] - min(combat_edata) + 1)
## Clean up
rm(data, RP, MT)

logMeans <- rowMeans(Y)
logVars <- rowVars(as.matrix(Y))
logCV2 <- logVars / logMeans^2

#hist(logMeans, breaks = 1000, xlim = range(0:8), col = "cyan3")
#abline(v = 1, lty = 2, col = "red", lwd = 2)
#quantile(logMeans, probs = 0.75)
useForFitL <- logMeans > 1
LogNcountsList <- data.frame(mean = logMeans, cv2 = logCV2)

#fit_loglin = nls(log10(cv2) ~ a*10^(-k*mean), LogNcountsList, start=c(a=1,k=1))
fit_loglin = lm(log10(cv2)~mean, data = LogNcountsList[useForFitL,])
#fit_loglin = lm(log10(cv2)~poly(mean, 4), data = LogNcountsList[useForFitL,])
ElogCV2 <- predict(fit_loglin, data.frame(mean = LogNcountsList$mean)) - summary(fit_loglin)$sigma
#ElogCV2 <- fit_loglin$coefficients[2]*logMeans + fit_loglin$coefficients[1]

is_het <- (ElogCV2 < log10(logCV2)) & useForFitL
#length(which(is_het))

LogNcountsList$ishet <- is_het 
LogNcountsList$ElogCV2 <- ElogCV2

pdf(file = paste("./images/", Sys.Date(), "_Mida_variable_genes.pdf", sep = ""), width = 11, height = 10, useDingbats = F)
ggplot(data = LogNcountsList, aes(x = mean, y = log10(cv2), colour = ishet)) +
  geom_point() +
  #geom_point(data = LogNcountsList[markers,], colour = "green", size = 2) +
  #eom_text(data = LogNcountsList[markers,], 
  #          label = substr(markers, 17, 100), 
  #          colour = "black", size = 3, vjust = -0.5) +
  geom_line(data = LogNcountsList, aes(x = logMeans, y = ElogCV2), colour = "black", 
            linetype = 2, size = 0.5) +
  geom_vline(xintercept = 1, colour = "black", size = 0.5, linetype = 2) +
  scale_color_discrete(name = "", labels=c("Stable", "Variable")) +
  theme_bw() +
  theme(text = element_text(size = 20),
        line = element_line(size = 1),
        panel.grid.minor = element_line(linetype = 3, colour = "grey60"),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.key = element_blank(),
        legend.position = c(0.85,0.85),
        legend.background = element_blank())
dev.off(which = dev.cur())

Y2 <- Y[is_het, ]

## Clean up
rm(logMeans, logVars, logCV2, useForFitL, LogNcountsList, fit_loglin, ElogCV2)

#### Find batch related genes with Spec analysis ####
batch <- GB$IFC
Data <- Y
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, outfile='LOG.TXT')
clusterExport(cl, c("batch", "Data"))

spec_scores_batch <- parApply(cl = cl, Data, 1, function(x){
  source("~/Scripts/R/spec.R");
  gene <- as.numeric(x);
  opt_bin <- optimizebinsize(gene, header = batch);
  unlist(specy(gene, header = batch, binsize = opt_bin));
})
stopCluster(cl)

rownames(spec_scores_batch)
rnames <- c(unlist(lapply(unique(sort(GB$IFC)),function(x){paste("Spec score", x)})),
            unlist(lapply(unique(sort(GB$IFC)),function(x){paste("mean/bin size", x)})))
rownames(spec_scores_batch) <- rnames

write.csv(spec_scores_batch, file = "spec_scores_batch.csv")

i <- 2
#spec_genes <- order(spec_scores[i,], decreasing = T)[1:100]
spec_genes <- unique(unlist(apply(spec_scores_batch[1:nrow(spec_scores_batch)/2],1,function(x){which(x > 0.3)})))

temp <- Data[spec_genes[1:100], ]
ord <- order(GB$IFC, GB$region, GB$time_point)
ColSideColors <- c(gg_color_hue(length(unique(GB$IFC)))[as.numeric(as.factor(GB$IFC))])

ColSideColors <- matrix(data = ColSideColors, nrow = ncol(temp));
ColSideColors <- ColSideColors[ord, ];

colsep <- lapply(unique(GB$IFC[ord]), 
                 function(x){length(which(GB$IFC == x))});
colsep <- cumsum(unlist(colsep));  

pairs.breaks <- seq(-4, 4, length.out=101);

pdf(file = paste("./images/", Sys.Date(), "_spec_heatmap_", i,"_batch.pdf", sep = ""), 
    width = 20, height = 10, useDingbats = F);
heatmap.3(temp[, ord],
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 1, 
          main="Marker gene expression",
          col = bluered(100),
          symkey = F,
          cexRow=1, cexCol = 0.6, 
          Rowv = F, 
          Colv = F,#as.dendrogram(hca), 
          ColSideColors = as.matrix(ColSideColors),
          ColSideColorsSize = 2,
          #RowSideColors = RowSideColors,
          dendrogram = "none",
          scale = "row",
          colsep = colsep,
          sepcolor = "black",
          labRow = substr(rownames(temp), start = 20, stop = 100), #c("Glast", "Tbr2"),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

Y2 <- Y[ -spec_genes, ]

#### 4. Principal component analysis ####
## Run PCA
pca0 <- prcomp(t(Y2), center = T, scale. = T)

sumData0=summary(pca0)[[1]];
ratio1=round(100*sumData0[1]^2/sum(sumData0^2),2);
ratio2=round(100*sumData0[2]^2/sum(sumData0^2),2);
ratio3=round(100*sumData0[3]^2/sum(sumData0^2),2);
#plot(pca, xlab = "Principle components" )

pdf(file = paste("./images/", Sys.Date(), "_Mida_pca_mouse.pdf", sep = ""), 
  width = 12, height = 10, useDingbats = F)
  ggplot(data = data.frame(pca0$x), aes(x = PC1, y = PC2, colour = as.factor(Normalized_GB$mouse))) +
  geom_point(size = 3) +
  scale_color_discrete(name = "Mouse") +
  labs(x=paste("PC1 ",ratio1,"%",sep=""),y=paste("PC2 ",ratio2,"%",sep="")) +
  theme_bw() +
  theme(text = element_text(size = 20),
        line = element_line(size = 1),
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        legend.key = element_blank(),
        panel.grid = element_blank())
dev.off(which = dev.cur())

df <- data.frame(pca0$x[,1:10])
df$IFC <- Normalized_GB$IFC
df$time_point <- Normalized_GB$time_point
df$region <- Normalized_GB$region

pdf(file = paste("./images/", Sys.Date(), "_pca_pairs_region.pdf", sep = ""), 
    width = 20, height = 20, useDingbats = F)
ggpairs(df, columns = 1:5, 
        mapping = ggplot2::aes(colour = region),
        lower = list(continuous = wrap("points", size = 1)),
        upper = "blank",
        diag = "blank")
dev.off(which = dev.cur())

## Clean up
rm(pca0, sumData0, ratio1, ratio2, ratio3)

#### 6. Continue analysis with variable genes only
#outliers <- c("C1.98.A6", "C1.98.A7", "C1.98.G3", "C1.98.G7", "C1.101.A6", "C1.101.A7",
#              "C1.98.G8", "C1.98.G9")
#GB2 <- GB[!(colnames(Y) %in% outliers),]
#save(Y,Y2,iPSC.GB, file = "2016-02-19_Y.Rdata")

#### 5. t-SNE ####
GB <- read.csv("20170521_clustering_info.csv", header = T)
for(i in c(1, 3, 5:6)){
  k <- keeps[[i]] & !(GB$outlier) # E12.5_dMGE #E14.5_dMGE
  GB.keep <- GB[k,]
  name <- names(keeps)[i]
  
  gene_mat <- as.matrix(Y[ , k ])

  logMeans <- rowMeans(gene_mat)
  logVars <- rowVars(as.matrix(gene_mat))
  logCV2 <- logVars / logMeans^2

  useForFitL <- logMeans > 1
  LogNcountsList <- data.frame(mean = logMeans, cv2 = logCV2)

  fit_loglin = lm(log10(cv2)~mean, data = LogNcountsList[useForFitL,])
  ElogCV2 <- predict(fit_loglin,
                     data.frame(mean = LogNcountsList$mean)) # - summary(fit_loglin)$sigma

  out_index <- (ElogCV2 < log10(logCV2)) & useForFitL
  Data <- gene_mat[out_index, ]
  
  seed <- 1024
  set.seed(seed)
  p = 20
  m = 0.1
  fm = 0.1
  t_sne <- Rtsne(t(Data),perplexity = p, pca = T, initial_dims = 5, verbose = F, theta = 0, momentum = m, 
                 final_momentum = fm,
                 #Y_init = as.matrix(t_sne),
                 max_iter = 10000)$Y
  colnames(t_sne) <- c("X1", "X2")
  rownames(t_sne) <- colnames(Data)
  t_sne <- data.frame(t_sne)
  
  #load(paste("2017-05-29", name, "tSNE.Rdata", sep = "_"))
  
  pdf(file = paste(Sys.Date(), name,"seed",seed, "p", p, "m", m,"fm",fm,
                   "tSNE.pdf", sep = "_"), 
      width = 12, height = 10, useDingbats = F)
  g <- ggplot(data = t_sne, aes(x = X1, y = X2)) +
    geom_point(size = 5, aes(x = X1, y = X2, color = GB.keep$cell_type)) +
    #scale_color_brewer(palette = "Paired", name = 'Ages') +
    scale_color_discrete(name = "cell_type") +
    #geom_text(aes(label = colnames(Y2)), hjust = 1, vjust = -1) + 
    labs(x = "tSNE 1", y = "tSNE 2") +
    #ggtitle(paste("All cells, Perplexity =", p)) +
    theme_bw() +
    theme(text = element_text(size = 20),
          line = element_line(size = 1),
          panel.grid = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black"),
          legend.key = element_blank(),
          legend.background = element_blank())
  print(g)
  dev.off(which = dev.cur())
  
  save(t_sne, file = paste(Sys.Date(), name, "tSNE.Rdata", sep = "_")) 
}

#### 7.1 Do pcaReduce ####
load("20170601_cluster_info_use_predicted.Rdata")
save(GB, file = "20170601_cluster_info_use_predicted.Rdata")

keep <- list(dMGE_14.5, vMGE_14.5, CGE_14.5, dMGE_12.5, vMGE_12.5, CGE_12.5)
names(keep) <- c("E14.5_dMGE", "E14.5_vMGE", "E14.5_CGE", "E12.5_dMGE", "E12.5_vMGE", "E12.5_CGE")

for(k in 1:6){
  GB$predicted[keep[[k]]] <- "outlier"
  name <- names(keep)[k]
  cc_output <- read.table(paste(name, "\\", name, "_outputfile.txt", sep = ""), header = F)
  clu <- lapply(unique(cc_output$V1), function(x){tail(which(cc_output$V1 == x), n = 1)})
  clu <- unlist(clu)
  clu <- cc_output[clu,]
  clu <- clu[na.omit(match(GB$cell_ID[keep[[k]]], clu$V1)),]
  GB$predicted[GB$cell_ID %in% clu$V1] <- clu$V2
}

for(k in 1:6){
  i <- keep[[k]]
  n <- names(keep)[k]
  temp <- na.omit(Y[ markers, i])
  
  ord <- order(GB$predicted[i], GB$region[i], GB$time_point[i])
  ColSideColors <- c(gg_color_hue(length(unique(GB$predicted)))[as.numeric(as.factor(GB$predicted[i]))],
                     rainbow(length(unique(GB$IFC[i])))[as.numeric(as.factor(GB$IFC[i]))])
  
  ColSideColors <- matrix(data = ColSideColors, nrow = ncol(temp));
  ColSideColors <- ColSideColors[ord, ];
  
  colsep <- lapply(unique(GB$predicted[i][ord]), 
                   function(x){length(which(GB$predicted[i] == x))});
  colsep <- cumsum(unlist(colsep));  
  
  pairs.breaks <- seq(-4, 4, length.out=101);
  
  pdf(file = paste("./", Sys.Date(), "_marker_heatmap_", 
                   length(unique(GB$predicted[i])),"_", n, "_predicted.pdf", 
                   sep = ""), 
      width = 20, height = 15, useDingbats = F);
  heatmap.3(temp[, ord],
            breaks = pairs.breaks,
            #symbreaks = T,
            keysize = 0.8, 
            main="Marker gene expression",
            col = bluered(100),
            symkey = F,
            cexRow=1, cexCol = 0.6, 
            Rowv = F, 
            Colv = F,#as.dendrogram(hca), 
            ColSideColors = as.matrix(ColSideColors),
            ColSideColorsSize = 2,
            #RowSideColors = RowSideColors,
            dendrogram = "none",
            scale = "row",
            colsep = colsep,
            sepcolor = "black",
            labRow = substr(rownames(temp), start = 20, stop = 100), #c("Glast", "Tbr2"),
            labCol = "",
            na.rm = F);
  dev.off(dev.cur());
}

for(k in 1:6){
  i <- keep[[k]]
  n <- names(keep)[k]
  temp <- na.omit(Y[ sf$V1, i])
  ord <- order(GB$predicted[i], GB$region[i], GB$time_point[i])
  ColSideColors <- c(gg_color_hue(length(unique(GB$predicted)))[as.numeric(as.factor(GB$predicted[i]))],
                     rainbow(length(unique(GB$IFC[i])))[as.numeric(as.factor(GB$IFC[i]))])
  
  ColSideColors <- matrix(data = ColSideColors, nrow = ncol(temp));
  ColSideColors <- ColSideColors[ord, ];
  
  colsep <- lapply(unique(GB$predicted[i][ord]), 
                   function(x){length(which(GB$predicted[i] == x))});
  colsep <- cumsum(unlist(colsep));  
  
  pairs.breaks <- seq(-4, 4, length.out=101);
  
  pdf(file = paste("./", Sys.Date(), "_marker_heatmap_", 
                   length(unique(GB$predicted[i])),"_", n, "_sf.pdf", 
                   sep = ""), 
      width = 20, height = 15, useDingbats = F);
  heatmap.3(temp[, ord],
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
            #RowSideColors = RowSideColors,
            dendrogram = "none",
            scale = "row",
            colsep = colsep,
            sepcolor = "black",
            labRow = substr(rownames(temp), start = 20, stop = 100), #c("Glast", "Tbr2"),
            labCol = "",
            na.rm = F);
  dev.off(dev.cur());
}


#### 7.1.2 use pcaReduce to do clustering ####
dMGE_14.5 <- GB$region == "dMGE" & GB$time_point == "E14.5" 
vMGE_14.5 <- GB$region == "vMGE" & GB$time_point == "E14.5"
CGE_14.5 <- GB$region == "CGE" & GB$time_point == "E14.5"
dMGE_12.5 <- GB$region == "dMGE" & GB$time_point == "E12.5"
vMGE_12.5 <- GB$region == "vMGE" & GB$time_point == "E12.5"
CGE_12.5 <- GB$region == "CGE" & GB$time_point == "E12.5"

GB$clusters <- NA
#load("E14.5_dMGE/2017-04-30_E14.5_dMGE_Output_S.Rdata") # genes selected by top 3 PCs
#load("E14.5_dMGE/2017-05-10_E14.5_dMGE_Output_S.Rdata") # genes selected by top 3 PCs
#load("E14.5_dMGE/2017-05-12_E14.5_dMGE_Output_S.Rdata") # Used M3Drop
#table(dm$Best.nc[1,])

load("E14.5_dMGE/2017-05-18_E14.5_dMGE_Output_S.Rdata") # genes selected by top 3 PCs
GB$clusters[GB$cell_ID %in% colnames(Data)] <- as.factor(Output_S[[1]][,7]) #15 clusters

#load("E14.5_vMGE/2017-04-30_E14.5_vMGE_Output_S.Rdata")
#load("E14.5_vMGE/2017-05-10_E14.5_vMGE_Output_S.Rdata") # genes selected by top 3 PCs
#load("E14.5_vMGE/2017-05-12_E14.5_vMGE_Output_S.Rdata") # Used M3Drop
#table(dm$Best.nc[1,])
load("E14.5_vMGE/2017-05-18_E14.5_vMGE_Output_S.Rdata") # genes selected by top 3 PCs
GB$clusters[GB$cell_ID %in% colnames(Data)] <- as.factor(Output_S[[1]][,12]) #10 clusters

#load("E14.5_CGE/2017-04-30_E14.5_CGE_Output_S.Rdata")
#load("E14.5_CGE/2017-05-10_E14.5_CGE_Output_S.Rdata") # genes selected by top 3 PCs
load("E14.5_CGE/2017-05-18_E14.5_CGE_Output_S.Rdata") # genes selected by top 3 PCs
GB$clusters[GB$cell_ID %in% colnames(Data)] <- as.factor(Output_S[[1]][,13]) #9 clusters


#load("E12.5_dMGE/2017-04-30_E12.5_dMGE_Output_S.Rdata")
#load("E12.5_dMGE/2017-05-10_E12.5_dMGE_Output_S.Rdata") # genes selected by top 3 PCs
load("E12.5_dMGE/2017-05-18_E12.5_dMGE_Output_S.Rdata") # genes selected by top 3 PCs
GB$clusters[GB$cell_ID %in% colnames(Data)] <- as.factor(Output_S[[1]][,12]) #10 clusters

#load("E12.5_vMGE/2017-04-30_E12.5_vMGE_Output_S.Rdata")
#load("E12.5_vMGE/2017-05-10_E12.5_vMGE_Output_S.Rdata")
load("E12.5_vMGE/2017-05-18_E12.5_vMGE_Output_S.Rdata") # genes selected by top 3 PCs
GB$clusters[GB$cell_ID %in% colnames(Data)] <- as.factor(Output_S[[1]][,10]) #12 clusters

#load("E12.5_CGE/2017-04-30_E12.5_CGE_Output_S.Rdata")
#load("E12.5_CGE/2017-05-10_E12.5_CGE_Output_S.Rdata")
load("E12.5_CGE/2017-05-18_E12.5_CGE_Output_S.Rdata") # genes selected by top 3 PCs
GB$clusters[GB$cell_ID %in% colnames(Data)] <- as.factor(Output_S[[1]][,10]) #12 clusters

GB$cell_type <- "O1"
GB$cell_type[ dMGE_14.5 & GB$clusters == 1 ] <- "T1"
GB$cell_type[ dMGE_14.5 & GB$clusters == 2 ] <- "N1"
GB$cell_type[ dMGE_14.5 & GB$clusters == 3 ] <- "N2"
GB$cell_type[ dMGE_14.5 & GB$clusters == 4 ] <- "P1"
GB$cell_type[ dMGE_14.5 & GB$clusters == 5 ] <- "N3"
GB$cell_type[ dMGE_14.5 & GB$clusters == 6 ] <- "N4"
GB$cell_type[ dMGE_14.5 & GB$clusters == 7 ] <- "T2"
GB$cell_type[ dMGE_14.5 & GB$clusters == 8 ] <- "T3"
GB$cell_type[ dMGE_14.5 & GB$clusters == 9 ] <- "T4"
GB$cell_type[ dMGE_14.5 & GB$clusters == 10 ] <- "P2"
GB$cell_type[ dMGE_14.5 & GB$clusters == 11 ] <- "T5"
GB$cell_type[ dMGE_14.5 & GB$clusters == 12 ] <- "N5"
GB$cell_type[ dMGE_14.5 & GB$clusters == 13 ] <- "T6"
GB$cell_type[ dMGE_14.5 & GB$clusters == 14 ] <- "T7"
GB$cell_type[ dMGE_14.5 & GB$clusters == 15 ] <- "P3"


GB$cell_type[ vMGE_14.5 & GB$clusters == 1 ] <- "N1"
GB$cell_type[ vMGE_14.5 & GB$clusters == 2 ] <- "P3"
GB$cell_type[ vMGE_14.5 & GB$clusters == 3 ] <- "P4"
GB$cell_type[ vMGE_14.5 & GB$clusters == 4 ] <- "P2"
GB$cell_type[ vMGE_14.5 & GB$clusters == 5 ] <- "N2"
GB$cell_type[ vMGE_14.5 & GB$clusters == 6 ] <- "N3"
GB$cell_type[ vMGE_14.5 & GB$clusters == 7 ] <- "N4"
GB$cell_type[ vMGE_14.5 & GB$clusters == 8 ] <- "P1"
GB$cell_type[ vMGE_14.5 & GB$clusters == 9 ] <- "T1"
GB$cell_type[ vMGE_14.5 & GB$clusters == 10 ] <- "T2"
#GB$cell_type[ vMGE_14.5 & GB$clusters == 11 ] <- "O1"
#GB$cell_type[ vMGE_14.5 & GB$clusters == 12 ] <- "T5"


GB$cell_type[ CGE_14.5 & GB$clusters == 1 ] <- "T1"
GB$cell_type[ CGE_14.5 & GB$clusters == 2 ] <- "T2"
GB$cell_type[ CGE_14.5 & GB$clusters == 3 ] <- "P1"
GB$cell_type[ CGE_14.5 & GB$clusters == 4 ] <- "T3"
GB$cell_type[ CGE_14.5 & GB$clusters == 5 ] <- "N1"
GB$cell_type[ CGE_14.5 & GB$clusters == 6 ] <- "N2"
GB$cell_type[ CGE_14.5 & GB$clusters == 7 ] <- "N3"
GB$cell_type[ CGE_14.5 & GB$clusters == 8 ] <- "T4"
GB$cell_type[ CGE_14.5 & GB$clusters == 9 ] <- "N4"


GB$cell_type[ dMGE_12.5 & GB$clusters == 1 ] <- "T2"
GB$cell_type[ dMGE_12.5 & GB$clusters == 2 ] <- "N3"
GB$cell_type[ dMGE_12.5 & GB$clusters == 3 ] <- "T1"
GB$cell_type[ dMGE_12.5 & GB$clusters == 4 ] <- "N4"
GB$cell_type[ dMGE_12.5 & GB$clusters == 5 ] <- "P1"
GB$cell_type[ dMGE_12.5 & GB$clusters == 6 ] <- "N2"
GB$cell_type[ dMGE_12.5 & GB$clusters == 7 ] <- "T3"
GB$cell_type[ dMGE_12.5 & GB$clusters == 8 ] <- "N1"
GB$cell_type[ dMGE_12.5 & GB$clusters == 9 ] <- "P2"
GB$cell_type[ dMGE_12.5 & GB$clusters == 10 ] <- "P3"
GB$cell_type[ dMGE_12.5 & GB$clusters == 11 ] <- "O1"
#GB$cell_type[ dMGE_12.5 & GB$clusters == 12 ] <- "T5"


GB$cell_type[ vMGE_12.5 & GB$clusters == 1 ] <- "T1" 
GB$cell_type[ vMGE_12.5 & GB$clusters == 2 ] <- "T2" #"P1"
GB$cell_type[ vMGE_12.5 & GB$clusters == 3 ] <- "P1" #"P2"
GB$cell_type[ vMGE_12.5 & GB$clusters == 4 ] <- "T3" #"T2"
GB$cell_type[ vMGE_12.5 & GB$clusters == 5 ] <- "N1" #"N1"
GB$cell_type[ vMGE_12.5 & GB$clusters == 6 ] <- "T4" #"P3"
GB$cell_type[ vMGE_12.5 & GB$clusters == 7 ] <- "P2" #"T3"
GB$cell_type[ vMGE_12.5 & GB$clusters == 8 ] <- "N2"
GB$cell_type[ vMGE_12.5 & GB$clusters == 9 ] <- "P3"#"P4"
GB$cell_type[ vMGE_12.5 & GB$clusters == 10 ] <- "N3" #"N3"
GB$cell_type[ vMGE_12.5 & GB$clusters == 11 ] <- "T5"
GB$cell_type[ vMGE_12.5 & GB$clusters == 12 ] <- "N4"

GB$cell_type[ CGE_12.5 & GB$clusters == 1 ] <- "T1"
GB$cell_type[ CGE_12.5 & GB$clusters == 2 ] <- "T2"
GB$cell_type[ CGE_12.5 & GB$clusters == 3 ] <- "batch1"  #"T1"
GB$cell_type[ CGE_12.5 & GB$clusters == 4 ] <- "N1"
GB$cell_type[ CGE_12.5 & GB$clusters == 5 ] <- "P1" #"T3"
GB$cell_type[ CGE_12.5 & GB$clusters == 6 ] <- "T3"
GB$cell_type[ CGE_12.5 & GB$clusters == 7 ] <- "N2"  #"N2"
GB$cell_type[ CGE_12.5 & GB$clusters == 8 ] <- "T4"  #"N3"
GB$cell_type[ CGE_12.5 & GB$clusters == 9 ] <- "batch2"  #"T4"
GB$cell_type[ CGE_12.5 & GB$clusters == 10 ] <- "T5"  #"T4"
GB$cell_type[ CGE_12.5 & GB$clusters == 11 ] <- "T6"  #"T4"
GB$cell_type[ CGE_12.5 & GB$clusters == 12 ] <- "P2"  #"T4"


GB$cell_type <- factor(GB$cell_type, 
                       levels = c("P1", "P2", "P3", "P4", "P5", "T1", "T2", "T3", 
                                  "T4", "T5", "T6", "T7","T8", "N1", "N2", "N3", "N4", 
                                  "N5", "O1", "batch1", "batch2"))

GB$outlier <- GB$cell_type %in% c("O1", "O2")

markers <- read.csv("20170323_Rubenstein_revised_gene_list.csv", header = F, stringsAsFactors = F)
markers <- unlist(markers[ ,1 ])
markers <- unique(sapply(markers, function(x){paste("\\b", x, "$", sep = "")}))
markers <- unlist(sapply(markers,function(x){grep(x, rownames(Y), value = T)}))

keep <- list(dMGE_14.5, vMGE_14.5, CGE_14.5, dMGE_12.5, vMGE_12.5, CGE_12.5)
names(keep) <- c("E14.5_dMGE", "E14.5_vMGE", "E14.5_CGE", "E12.5_dMGE", "E12.5_vMGE", "E12.5_CGE")

for(k in 1:6){
  i <- keep[[k]]
  n <- names(keeps)[k]
  temp <- Y[ markers, i ]
  GB.keep <- GB[i, ]
  load(paste(n, "/2017-05-18_",n, "_Output_S.Rdata", sep="")) # genes selected by top 3 PCs
  for(ii in 2:17){
    GB.keep$clusters[GB.keep$cell_ID %in% colnames(Data)] <- as.factor(Output_S[[1]][,ii])
    ord <- order(GB.keep$clusters, GB.keep$region, GB.keep$time_point)
    ColSideColors <- c(gg_color_hue(length(unique(GB.keep$clusters)))[as.numeric(GB.keep$clusters)],
                       rainbow(length(unique(GB.keep$IFC)))[as.numeric(as.factor(GB.keep$IFC))])
    
    ColSideColors <- matrix(data = ColSideColors, nrow = ncol(temp));
    ColSideColors <- ColSideColors[ord, ];
    
    colsep <- lapply(unique(GB.keep$clusters[ord]), 
                     function(x){length(which(GB.keep$clusters == x))});
    colsep <- cumsum(unlist(colsep));  
    
    pairs.breaks <- seq(-4, 4, length.out=101);
    
    #GB.keep$clusters <- as.factor(Output_S[[1]][,ii])
    pdf(file = paste(Sys.Date(), "marker_heatmap_", 
                     length(unique(GB.keep$clusters)), n, "clusters.pdf", 
                     sep = "_"), 
        width = 20, height = 15, useDingbats = F);
    heatmap.3(temp[, ord],
              breaks = pairs.breaks,
              #symbreaks = T,
              keysize = 0.8, 
              main="Marker gene expression",
              col = bluered(100),
              symkey = F,
              cexRow=1, cexCol = 0.6, 
              Rowv = F, 
              Colv = F,#as.dendrogram(hca), 
              ColSideColors = as.matrix(ColSideColors),
              ColSideColorsSize = 2,
              #RowSideColors = RowSideColors,
              dendrogram = "none",
              scale = "row",
              colsep = colsep,
              sepcolor = "black",
              labRow = substr(rownames(temp), start = 20, stop = 100), #c("Glast", "Tbr2"),
              labCol = "",
              na.rm = F);
    dev.off(dev.cur());
  }
}

for(k in 1:6){
  i <- keep[[k]]
  n <- names(keep)[k]
  temp <- na.omit(Y[ markers, i])
  ord <- order(GB$cell_type[i], GB$region[i], GB$time_point[i])
  ColSideColors <- c(gg_color_hue(length(unique(GB$cell_type)))[as.numeric(GB$cell_type[i])],
                     rainbow(length(unique(GB$IFC[i])))[as.numeric(as.factor(GB$IFC[i]))])
  
  ColSideColors <- matrix(data = ColSideColors, nrow = ncol(temp));
  ColSideColors <- ColSideColors[ord, ];
  
  colsep <- lapply(unique(GB$cell_type[i][ord]), 
                   function(x){length(which(GB$cell_type[i] == x))});
  colsep <- cumsum(unlist(colsep));  
  
  pairs.breaks <- seq(-4, 4, length.out=101);
  
  pdf(file = paste("./", Sys.Date(), "_marker_heatmap_", 
                   length(unique(GB$cell_type[i])),"_", n, "_cell_type.pdf", 
                   sep = ""), 
      width = 20, height = 15, useDingbats = F);
  heatmap.3(temp[, ord],
            breaks = pairs.breaks,
            #symbreaks = T,
            keysize = 0.8, 
            main="Marker gene expression",
            col = bluered(100),
            symkey = F,
            cexRow=1, cexCol = 0.6, 
            Rowv = F, 
            Colv = F,#as.dendrogram(hca), 
            ColSideColors = as.matrix(ColSideColors),
            ColSideColorsSize = 2,
            #RowSideColors = RowSideColors,
            dendrogram = "none",
            scale = "row",
            colsep = colsep,
            sepcolor = "black",
            labRow = substr(rownames(temp), start = 20, stop = 100), #c("Glast", "Tbr2"),
            labCol = "",
            na.rm = F);
  dev.off(dev.cur());
}

# set.seed(1024)
# p = 20
# m = 0.9
# t_sne <- Rtsne(t(Data), initial_dims = 12, perplexity = p, pca = T, verbose = F, theta = 0, 
#                momentum = m, final_momentum = 0.1, max_iter = 1000)$Y
# colnames(t_sne) <- c("X1", "X2")
# rownames(t_sne) <- GB.dMGE_14.5$cell_ID
# t_sne <- data.frame(t_sne)
# 
# pdf(file = paste("./E14.5_dMGE/", Sys.Date(), "_E14.5_dMGE_tSNE_perplexity_", p, "_momentum_", 
#                  m , ".pdf", sep = ""), 
#     width = 12, height = 10, useDingbats = F)
# ggplot(data = t_sne, aes(x = X1, y = X2)) +
#   geom_point(size = 5, aes(x = X1, y = X2, color = GB.dMGE_14.5$cell_type)) +
#   #scale_color_brewer(palette = "Paired", name = 'Ages') +
#   scale_color_discrete(name = "Cell Types") +
#   #geom_text(aes(label = colnames(Y2)), hjust = 1, vjust = -1) + 
#   labs(x = "tSNE 1", y = "tSNE 2") +
#   ggtitle(paste("All cells, Perplexity =", p)) +
#   theme_bw() +
#   theme(text = element_text(size = 20),
#         line = element_line(size = 1),
#         panel.grid = element_blank(),
#         panel.border = element_rect(fill = NA, colour = "black"),
#         legend.key = element_blank(),
#         legend.background = element_blank())
# dev.off(which = dev.cur())

## Clean up
rm(colsep, Output_S, seed, x, nbt, q, n, ord, i, GB.dMGE_14.5, dMGE_14.5)

#### 7.2 Random Forrest classification ####
GB$predicted <- "none"
for(k in 1:6){
  i <- keeps[[k]]
  name <- names(keeps)[k]
  load(paste(name, "/2017-04-30_", name, "_Output_S.Rdata", sep = "")) # genes selected by top 3 PCs
  feature_bf_fs <- as.matrix(t(Data))
  cate_bf_fs <- as.factor(as.character(GB$cell_type[i]))
  set <- sample(1: nrow(feature_bf_fs), nrow(feature_bf_fs), replace=F)  ## Random reorder
  cate_bf_fs <- cate_bf_fs[set]
  feature_bf_fs <- feature_bf_fs[set, ]
  rf_bf_fs <- randomForest(t(Data), as.factor(as.character(GB$cell_type[i])), importance = TRUE, proximity = TRUE)
  imp_bf_fs <- importance(rf_bf_fs, type = 1)

  fs <- rfcv(feature_bf_fs, cate_bf_fs, cv.fold = 10, scale = "log", step = 0.9, recursive = F)
  len <- length(fs$n.var[fs$error.cv == min(fs$error.cv)])
  min_fs <- fs$n.var[fs$error.cv == min(fs$error.cv)][len] #get least features
  ind <- order(-imp_bf_fs)[1: min_fs]
  feature_fs <- feature_bf_fs[ , ind]
  cate_fs <- cate_bf_fs

  rf_fs <- randomForest(feature_fs, as.factor(cate_fs), importance=TRUE, proximity=TRUE)
  #rf_fs

  GB$predicted[i] <- rf_fs$predicted[match(GB$cell_ID[i], names(rf_fs$predicted))]

  # feature_fs <- as.matrix(t(Data))
  # 
  # cate_fs <- as.factor(as.character(GB$cell_type[i]))
  # # set <- sample(1: nrow(feature_fs), nrow(feature_fs), replace=F)  ## Random reorder
  # # cate_fs <- cate_fs[set]
  # # feature_fs <- feature_fs[set, ]
  # 
  # rf_fs <- randomForest(feature_fs, cate_fs, importance = TRUE, proximity = TRUE)
  # feature_f <- feature_fs[rowMax(rf_fs$votes) > 0.5, ]
  # cate_f <- as.factor(as.character(rf_fs$predicted[rowMax(rf_fs$votes) > 0.5]))
  # 
  # rf_f <- randomForest(feature_f, cate_f, importance=TRUE, proximity=TRUE)
  # GB$predicted_cell_type[i] <- factor(predict(rf_f, newdata = feature_fs), 
  #                                        levels = levels(GB$predicted_cell_type))
}

for(k in 1:6){
  i <- keep[[k]]
  n <- names(keep)[k]
  temp <- Y2[ markers, i ]
  ord <- order(GB$predicted[i], GB$region[i], GB$time_point[i])
  ColSideColors <- c(gg_color_hue(length(unique(GB$predicted)))[as.numeric(GB$predicted[i])],
                     rainbow(length(unique(GB$IFC[i])))[as.numeric(as.factor(GB$IFC[i]))])
  
  ColSideColors <- matrix(data = ColSideColors, nrow = ncol(temp));
  ColSideColors <- ColSideColors[ord, ];
  
  colsep <- lapply(unique(GB$predicted[i][ord]), 
                   function(x){length(which(GB$predicted[i] == x))});
  colsep <- cumsum(unlist(colsep));  
  
  pairs.breaks <- seq(-4, 4, length.out=101);
  
  pdf(file = paste("./", Sys.Date(), "_marker_heatmap_", 
                   length(unique(GB$predicted[i])),"_", n, "_predicted.pdf", 
                   sep = ""), 
      width = 20, height = 15, useDingbats = F);
  heatmap.3(temp[, ord],
            breaks = pairs.breaks,
            #symbreaks = T,
            keysize = 0.8, 
            main="Marker gene expression",
            col = bluered(100),
            symkey = F,
            cexRow=1, cexCol = 0.6, 
            Rowv = F, 
            Colv = F,#as.dendrogram(hca), 
            ColSideColors = as.matrix(ColSideColors),
            ColSideColorsSize = 2,
            #RowSideColors = RowSideColors,
            dendrogram = "none",
            scale = "row",
            colsep = colsep,
            sepcolor = "black",
            labRow = substr(rownames(temp), start = 20, stop = 100), #c("Glast", "Tbr2"),
            labCol = "",
            na.rm = F);
  dev.off(dev.cur());
}


# feature <- data.frame()
# cate <- c()
# for(i in 1:length(unique(rf_fs$predicted))){
#   fea_fs <- feature_fs[(rf_fs$predicted == colnames(rf_fs$votes)[i]) & (rf_fs$votes[ , i] > 0.2), , drop = FALSE]
#   feature <- rbind(feature, fea_fs)
#   cat_fs <- rf_fs$predicted[(rf_fs$predicted ==colnames(rf_fs$votes)[i]) & (rf_fs$votes[ , i] > 0.2)]
#   cate <- as.factor(c(as.character(cate), as.character(cat_fs)))
# }
# 
# set <- sample(1: nrow(feature), nrow(feature), replace = F)
# cate <- cate[set]
# feature <- feature[set, ] 


#### Plot MAF and MAFB
genes <- unlist( sapply( c("\\bMAF$", "\\bMAFB$", "\\bSST$", "\\bCCK$", "\\bNPY$"),
                         function( x ){ grep( x, rownames(cpm), value = T )})) 
temp <- t( Y[ genes, dMGE_14.5 ] )
temp <- data.frame( temp )
colnames(temp) <- substr( colnames( temp ), 20, 100 )
#cor27 <- unlist(apply(temp, 2, function(x) {cor.test(as.numeric(res)[GB$time_point %in% c("E17.5", "P1")], x, method = "s")[3]}))
temp$Time_point <- GB.dMGE_14.5$time_point#[ res ]
temp$cell_type <- GB.dMGE_14.5$cell_type#[ res ]
temp1 <- melt(temp, id.vars = c("Time_point", "cell_type"))

pdf(file = paste("./images/", Sys.Date(), "_genes_cell_type.pdf", sep = ""), 
    width = 10*13/22, height = 4, useDingbats = F)
ggplot(data = temp1, aes(x = cell_type, y = value)) + 
  geom_violin(aes(fill="red"), scale = "width", color = "black", 
              size = 0.1) +
  labs(x = '', y = '') +
  #ylim(c(0, 13)) +
  #geom_jitter(aes(group = Time_point), temp1) +
  #ggtitle(paste("")) +
  #stat_smooth(method = "auto", aes(group = variable)) +
  facet_grid(variable ~ .) +
  theme_bw() +
  theme(text = element_text(size = 6),
        line = element_line(size = 0.1),
        axis.text.x=element_text(angle = 45, size = 10, hjust = 1),
        axis.title = element_blank(),
        panel.grid = element_blank(),#element_line(linetype = 3, colour = "grey60"),
        panel.border = element_rect(colour = "black"),
        legend.key = element_blank(),
        legend.position = "",
        legend.background = element_blank())
dev.off(which = dev.cur())
write.csv(GB.dMGE_14.5, file = paste(Sys.Date(), "dMGE_14.5_cluster_information.csv", 
                                     sep = "_"))

#### Monocle 2 analysis ####
#### Setup dataset

## load count data
#cds <- read.table("raw_data/mouseBrain.Mida.singleCell.gene.count.txt", header = T)
# rownames(cds) <- cds$Geneid
# cds <- cds[, -1]
# colnames(cds) <- gsub("90_", "90.", colnames(cds))
# colnames(cds) <- gsub("91_", "91.", colnames(cds))
# colnames(cds) = unlist(lapply(colnames(cds), function(x){strsplit(x, split = "_")[[1]][1]}))
#expr <- cds[ , colnames(Data)]

#expr <- df[rownames(Y), colnames(Data)]
GB$outlier <- GB$cell_type %in% c("O1", "O2")
keep <- list(dMGE_14.5, vMGE_14.5, CGE_14.5, dMGE_12.5, vMGE_12.5, CGE_12.5)
names(keep) <- c("E14.5_dMGE", "E14.5_vMGE", "E14.5_CGE", "E12.5_dMGE", "E12.5_vMGE", "E12.5_CGE")

set.seed(1024)
mms <- list()
ordering_genes <- list()
for(k in 1:6){
  df  <- dfs[k]
  pd <- new("AnnotatedDataFrame", data = df@meta.data)
  rownames(pd) <- df@meta.data$cell_ID
  fd <- data.frame(gene_ID = rownames(expr), 
                   gene_short_name = substr(rownames(expr), start = 20, stop = 100))
  rownames(fd) <- fd$gene_ID
  fd <- new("AnnotatedDataFrame", data = fd)
  mm <- newCellDataSet(as(as.matrix(expr), "sparseMatrix"),
                       phenoData = pd, 
                       featureData = fd,
                       lowerDetectionLimit = 0.1, 
                       expressionFamily = tobit(Lower = 0.1))
  
  
  ## convert RPKM to RPC
  #rpc_matrix <- relative2abs(mm)
  
  # mm <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
  #                      phenoData = pd, featureData = fd,
  #                      lowerDetectionLimit = 0.5,
  #                      expressionFamily = negbinomial.size())
  
  ## clean up
  rm(expr, rpc_matrix, pd, fd)
  
  ## calculate size factor and dispersions
  #mm <- estimateSizeFactors(mm)
  #mm <- estimateDispersions(mm)
  
  ## filtering low-quality cells
  mm <- detectGenes(mm, min_expr = 0.1)
  #head(fData(mm))
  
  expr_genes <- row.names(subset(fData(mm), num_cells_expressed >= 10))
  
  ## take a look at the distribution of mRNA totals across the cells
  pData(mm)$Total_mRNAs <- colSums(exprs(mm))
  mm <- mm[,pData(mm)$Total_mRNAs < 1e6]
  upper_bound <- 10^(mean(log10(pData(mm)$Total_mRNAs)) + 2*sd(log10(pData(mm)$Total_mRNAs)))
  lower_bound <- 10^(mean(log10(pData(mm)$Total_mRNAs)) - 2*sd(log10(pData(mm)$Total_mRNAs)))
  
  pdf(file = paste("./", name, "/", Sys.Date(), "_", name, "_mRNA_totals.pdf", sep = ""), 
      width = 10, height = 10, useDingbats = F)
  qplot(Total_mRNAs, data=pData(mm), color=time_point, geom="density") +
    geom_vline(xintercept=lower_bound) +
    geom_vline(xintercept=upper_bound)
  dev.off()
  
  mm <- mm[,pData(mm)$Total_mRNAs > lower_bound &
             pData(mm)$Total_mRNAs < upper_bound]
  mm <- detectGenes(mm, min_expr = 0.1)
  
  #mm <- reduceDimension(mm, max_components = 2, num_dim = 11, reduction_method = 'tSNE',
  #                      verbose = T)
  #pData(mm)$cell_type <- GB.keep$cell_type[GB.keep$cell_ID %in% pData(mm)$cell_ID]
  #
  #mm <- clusterCells(mm, num_clusters = 12)
  
  #pdf(file = paste("./", name, "/", Sys.Date(), "_", name, "_cell_type.pdf", sep = ""), 
  #    width = 10, height = 10, useDingbats = F)
  #plot_cell_clusters(mm, 1, 2, color = "cell_type", cell_size = 3)
  #dev.off(which = dev.cur())
  
  # data<- log2(exprs(mm) + 1)
  # logMeans <- rowMeans(data)
  # logVars <- rowVars(as.matrix(data))
  # logCV2 <- logVars / logMeans^2
  # 
  # useForFitL <- logMeans > 1
  # LogNcountsList <- data.frame(mean = logMeans, cv2 = logCV2)
  # 
  # fit_loglin = lm(log10(cv2)~mean, data = LogNcountsList[useForFitL,])
  # ElogCV2 <- predict(fit_loglin, data.frame(mean = LogNcountsList$mean)) - summary(fit_loglin)$sigma
  # 
  # is_het <- (ElogCV2 < log10(logCV2)) & useForFitL
  # 
  # LogNcountsList$ishet <- is_het 
  # LogNcountsList$ElogCV2 <- ElogCV2
  
  # pdf(file = paste("./images/", Sys.Date(), "_Mida_variable_genes.pdf", sep = ""), width = 11, height = 10, useDingbats = F)
  # ggplot(data = LogNcountsList, aes(x = mean, y = log10(cv2), colour = ishet)) +
  #   geom_point() +
  #   #geom_point(data = LogNcountsList[markers,], colour = "green", size = 2) +
  #   #eom_text(data = LogNcountsList[markers,], 
  #   #          label = substr(markers, 17, 100), 
  #   #          colour = "black", size = 3, vjust = -0.5) +
  #   geom_line(data = LogNcountsList, aes(x = logMeans, y = ElogCV2), colour = "black", 
  #             linetype = 2, size = 0.5) +
  #   geom_vline(xintercept = 1, colour = "black", size = 0.5, linetype = 2) +
  #   scale_color_discrete(name = "", labels=c("Stable", "Variable")) +
  #   theme_bw() +
  #   theme(text = element_text(size = 20),
  #         line = element_line(size = 1),
  #         panel.grid.minor = element_line(linetype = 3, colour = "grey60"),
  #         panel.border = element_rect(fill = NA, colour = "black"),
  #         legend.key = element_blank(),
  #         legend.position = c(0.85,0.85),
  #         legend.background = element_blank())
  # dev.off(which = dev.cur())
  ordering_genes[[k]] <- differentialGeneTest(mm, fullModelFormulaStr = "~cell_type")
  mms[[k]] <- mm
}  

# diff_test_res <- list()
# for(k in 1:6){
#   mm <- mms[[k]]
#   dtr <- differentialGeneTest(cds, fullModelFormulaStr = "~predicted_cell_type")
#   diff_test_res[[k]] <- dtr
# }
save(mms, file = "2017-05-18_mms.Rdata")

for(k in 1:6){
  name <- names(keep)[k]
  mm <- mms[[k]]
  mm <- reduceDimension(mm, max_components = 2, verbose = F, norm_method = "none")
  mm <- setOrderingFilter(mm, ordering_genes = ordering_genes[[k]])
  mm <- orderCells(mm)
  mms[[k]] <- mm 
  
  ## plot trajectory colored by cell type
  pdf(file = paste("C:/Users/zhenadmin/Google Drive/scRNA_analysis/18_May/monocle/", name, "/", Sys.Date(), "_", name, "_trajectory_cell_type.pdf", sep = ""), 
      width = 11, height = 10, useDingbats = F)
  print(plot_cell_trajectory(mm, color_by="cell_type", cell_size = 3))
  dev.off()
  
  ## plot trajectory colored by pseudotime
  pdf(file = paste("C:/Users/zhenadmin/Google Drive/scRNA_analysis/18_May/monocle/", name, "/", Sys.Date(), "_", name, "_trajectory_pseudotime.pdf", sep = ""), 
      width = 11, height = 10, useDingbats = F)
  print(plot_cell_trajectory(mm, color_by="Pseudotime", cell_size = 3))
  dev.off()
  
  ## plot trajectory colored by individual cell type
  pdf(file = paste("C:/Users/zhenadmin/Google Drive/scRNA_analysis/18_May/monocle/", name, "/", Sys.Date(), "_", name, "_trajectory_cell_type_wrap.pdf", sep = ""), 
      width = 11, height = 10, useDingbats = F)
  print(plot_cell_trajectory(mm, color_by="cell_type") + facet_wrap(~cell_type, nrow=3))
  dev.off()
  
  ## plot trajectory colored by predicted cell type
  # pdf(file = paste("C:/Users/zhenadmin/Google Drive/scRNA_analysis/18_May/monocle/", name, "/", Sys.Date(), "_", name, "_trajectory_predicted_wrap.pdf", sep = ""), 
  #     width = 11, height = 10, useDingbats = F)
  # print(plot_cell_trajectory(mm, color_by="predicted") + facet_wrap(~predicted, nrow=3))
  # dev.off()
  
  ## plot trajectory colored by state
  pdf(file = paste("C:/Users/zhenadmin/Google Drive/scRNA_analysis/18_May/monocle/", name, "/", Sys.Date(), "_", name, "_trajectory_state.pdf", sep = ""), 
      width = 11, height = 10, useDingbats = F)
  print(plot_cell_trajectory(mm, color_by="State") + facet_wrap(~State, nrow=3))
  dev.off()
  
  ## plot trajectory colored by predicted cell type
  # pdf(file = paste("C:/Users/zhenadmin/Google Drive/scRNA_analysis/18_May/monocle/", name, "/", Sys.Date(), "_", name, "_trajectory_predicted_cell_type.pdf", sep = ""), 
  #     width = 11, height = 10, useDingbats = F)
  # print(plot_cell_trajectory(mm, color_by="predicted_cell_type") + 
  #         facet_wrap(~predicted_cell_type, nrow=3))
  # dev.off()
}

corrs <- list()
for(k in 1:6){
  name <- names(keep)[k]
  mm <- mms[[k]]
  
  dt <- exprs(mm)
  logMeans <- rowMeans(dt)
  logVars <- rowVars(as.matrix(dt))
  logCV2 <- logVars / logMeans^2
  
  useForFitL <- logMeans > 1
  LogNcountsList <- data.frame(mean = logMeans, cv2 = logCV2)
  
  fit_loglin = lm(log10(cv2)~mean, data = LogNcountsList[useForFitL,])
  ElogCV2 <- predict(fit_loglin, data.frame(mean = LogNcountsList$mean)) - summary(fit_loglin)$sigma
  
  fData(mm)$is_het <- (ElogCV2 < log10(logCV2)) & useForFitL
  
  pt <- pData(mm)$Pseudotime
  corr <- unlist(apply(exprs(mm[fData(mm)$is_het,]) ,1 , function(x){cor(x, pt)}))
  corr <- sort(corr, decreasing = T)
  corrs[[k]] <- corr
}

for(k in 1:6){
  name <- names(keep)[k]
  mm <- mms[[k]]
  markers <- names(corrs[[k]])[1:100]
  
  cds <- mm[markers,]
  pdf(file = paste("C:/Users/zhenadmin/Google Drive/scRNA_analysis/18_May/monocle/", name, "/", Sys.Date(), "_", name, "_correlted_genes_in_pseudotime.pdf", sep = ""), 
      width = 20, height = 20, useDingbats = F)
  g <- plot_genes_in_pseudotime(cds, color_by="cell_type", ncol = 10)
  print(g)
  dev.off()
}

for(k in 1:6){
  name <- names(keep)[k]
  mm <- mms[[k]]
  markers <- names(tail(corrs[[k]], n = 100))
  
  cds <- mm[markers,]
  pdf(file = paste("C:/Users/zhenadmin/Google Drive/scRNA_analysis/18_May/monocle/", name, "/", Sys.Date(), "_", name, "_anti-correlted_genes_in_pseudotime.pdf", sep = ""), 
      width = 20, height = 20, useDingbats = F)
  g <- plot_genes_in_pseudotime(cds, color_by="cell_type", ncol = 10)
  print(g)
  dev.off()
}

for(k in 1:6){
  name <- names(keep)[k]
  write.csv(corrs[[k]], file = paste("C:/Users/zhenadmin/Google Drive/scRNA_analysis/18_May/monocle/", name, "_gene_correlation.csv", sep = ""))
}



#### predicted cell type ####
dMGE_14.5 <- GB$region == "dMGE" & GB$time_point == "E14.5" 
vMGE_14.5 <- GB$region == "vMGE" & GB$time_point == "E14.5"
CGE_14.5 <- GB$region == "CGE" & GB$time_point == "E14.5"
dMGE_12.5 <- GB$region == "dMGE" & GB$time_point == "E12.5"
vMGE_12.5 <- GB$region == "vMGE" & GB$time_point == "E12.5"
CGE_12.5 <- GB$region == "CGE" & GB$time_point == "E12.5"

GB$predicted_cell_type <- "none"
GB$predicted_cell_type[ dMGE_14.5 & GB$predicted == 1 ] <- "N1"
GB$predicted_cell_type[ dMGE_14.5 & GB$predicted == 2 ] <- "N2"
GB$predicted_cell_type[ dMGE_14.5 & GB$predicted == 3 ] <- "T3"
GB$predicted_cell_type[ dMGE_14.5 & GB$predicted == 4 ] <- "N3"
GB$predicted_cell_type[ dMGE_14.5 & GB$predicted == 5 ] <- "T1"
GB$predicted_cell_type[ dMGE_14.5 & GB$predicted == 6 ] <- "T2"
GB$predicted_cell_type[ dMGE_14.5 & GB$predicted == 7 ] <- "N4"
GB$predicted_cell_type[ dMGE_14.5 & GB$predicted == 8 ] <- "N5"
GB$predicted_cell_type[ dMGE_14.5 & GB$predicted == 9 ] <- "P1"
GB$predicted_cell_type[ dMGE_14.5 & GB$predicted == 10 ] <- "P2"
GB$predicted_cell_type[ dMGE_14.5 & GB$predicted == 11 ] <- "P3"
GB$predicted_cell_type[ dMGE_14.5 & GB$predicted == 12 ] <- "T4"

GB$predicted_cell_type[ vMGE_14.5 & GB$predicted == 1 ] <- "N1"
GB$predicted_cell_type[ vMGE_14.5 & GB$predicted == 2 ] <- "T1"
GB$predicted_cell_type[ vMGE_14.5 & GB$predicted == 3 ] <- "N2"
GB$predicted_cell_type[ vMGE_14.5 & GB$predicted == 4 ] <- "N3"
GB$predicted_cell_type[ vMGE_14.5 & GB$predicted == 5 ] <- "N4"
GB$predicted_cell_type[ vMGE_14.5 & GB$predicted == 6 ] <- "P1"
GB$predicted_cell_type[ vMGE_14.5 & GB$predicted == 7 ] <- "P2"
GB$predicted_cell_type[ vMGE_14.5 & GB$predicted == 8 ] <- "P3"
GB$predicted_cell_type[ vMGE_14.5 & GB$predicted == 9 ] <- "T2"
GB$predicted_cell_type[ vMGE_14.5 & GB$predicted == 10 ] <- "P4"
GB$predicted_cell_type[ vMGE_14.5 & GB$predicted == 11 ] <- "N5"
GB$predicted_cell_type[ vMGE_14.5 & GB$predicted == 12 ] <- "T3"
GB$predicted_cell_type[ vMGE_14.5 & GB$predicted == 13 ] <- "T4"

GB$predicted_cell_type[ CGE_14.5 & GB$predicted == 1 ] <- "N1"
GB$predicted_cell_type[ CGE_14.5 & GB$predicted == 2 ] <- "N2"
GB$predicted_cell_type[ CGE_14.5 & GB$predicted == 3 ] <- "N3"
GB$predicted_cell_type[ CGE_14.5 & GB$predicted == 4 ] <- "P1"
GB$predicted_cell_type[ CGE_14.5 & GB$predicted == 5 ] <- "T2"
GB$predicted_cell_type[ CGE_14.5 & GB$predicted == 6 ] <- "T1"
GB$predicted_cell_type[ CGE_14.5 & GB$predicted == 7 ] <- "P2"

GB$predicted_cell_type[ dMGE_12.5 & GB$predicted == 1 ] <- "N1"
GB$predicted_cell_type[ dMGE_12.5 & GB$predicted == 2 ] <- "N2"
GB$predicted_cell_type[ dMGE_12.5 & GB$predicted == 3 ] <- "N3"
GB$predicted_cell_type[ dMGE_12.5 & GB$predicted == 4 ] <- "N4"
GB$predicted_cell_type[ dMGE_12.5 & GB$predicted == 5 ] <- "P1"
GB$predicted_cell_type[ dMGE_12.5 & GB$predicted == 6 ] <- "P2"
GB$predicted_cell_type[ dMGE_12.5 & GB$predicted == 7 ] <- "T1"
GB$predicted_cell_type[ dMGE_12.5 & GB$predicted == 8 ] <- "T2"
GB$predicted_cell_type[ dMGE_12.5 & GB$predicted == 9 ] <- "T3"
GB$predicted_cell_type[ dMGE_12.5 & GB$predicted == 10 ] <- "N5"
GB$predicted_cell_type[ dMGE_12.5 & GB$predicted == 11 ] <- "T4"

GB$predicted_cell_type[ vMGE_12.5 & GB$predicted == 1 ] <- "N1" 
GB$predicted_cell_type[ vMGE_12.5 & GB$predicted == 2 ] <- "N2" #"P1"
GB$predicted_cell_type[ vMGE_12.5 & GB$predicted == 3 ] <- "N3" #"P2"
GB$predicted_cell_type[ vMGE_12.5 & GB$predicted == 4 ] <- "N4" #"T2"
GB$predicted_cell_type[ vMGE_12.5 & GB$predicted == 5 ] <- "N5" #"N1"
GB$predicted_cell_type[ vMGE_12.5 & GB$predicted == 6 ] <- "P1" #"P3"
GB$predicted_cell_type[ vMGE_12.5 & GB$predicted == 7 ] <- "P2" #"T3"
GB$predicted_cell_type[ vMGE_12.5 & GB$predicted == 8 ] <- "T1"
GB$predicted_cell_type[ vMGE_12.5 & GB$predicted == 9 ] <- "T2"#"P4"
GB$predicted_cell_type[ vMGE_12.5 & GB$predicted == 10 ] <- "T3" #"N3"
GB$predicted_cell_type[ vMGE_12.5 & GB$predicted == 11 ] <- "T4"
GB$predicted_cell_type[ vMGE_12.5 & GB$predicted == 12 ] <- "T5"

GB$predicted_cell_type[ CGE_12.5 & GB$predicted == 1 ] <- "N1"
GB$predicted_cell_type[ CGE_12.5 & GB$predicted == 2 ] <- "N2"
GB$predicted_cell_type[ CGE_12.5 & GB$predicted == 3 ] <- "P1"  #"T1"
GB$predicted_cell_type[ CGE_12.5 & GB$predicted == 4 ] <- "P2"
GB$predicted_cell_type[ CGE_12.5 & GB$predicted == 5 ] <- "P3" #"T3"
GB$predicted_cell_type[ CGE_12.5 & GB$predicted == 6 ] <- "P4"
GB$predicted_cell_type[ CGE_12.5 & GB$predicted == 7 ] <- "T1"  #"N2"
GB$predicted_cell_type[ CGE_12.5 & GB$predicted == 8 ] <- "T2"  #"N3"
GB$predicted_cell_type[ CGE_12.5 & GB$predicted == 9 ] <- "T3"  #"T4"

GB$predicted_cell_type <- factor(GB$predicted_cell_type, 
                       levels = c("P1", "P2", "P3", "P4", "T1", "T2", "T3", 
                                  "T4", "T5", "N1", "N2", "N3", "N4", "N5"))

for(k in 1:6){
  i <- keep[[k]]
  n <- names(keep)[k]
  temp <- Y[ markers, i ]
  ord <- order(GB$predicted_cell_type[i], GB$region[i], GB$time_point[i])
  ColSideColors <- c(gg_color_hue(length(unique(GB$predicted_cell_type)))[as.numeric(GB$predicted_cell_type[i])],
                     rainbow(length(unique(GB$IFC[i])))[as.numeric(as.factor(GB$IFC[i]))])
  
  ColSideColors <- matrix(data = ColSideColors, nrow = ncol(temp));
  ColSideColors <- ColSideColors[ord, ];
  
  colsep <- lapply(unique(GB$predicted_cell_type[i][ord]), 
                   function(x){length(which(GB$predicted_cell_type[i] == x))});
  colsep <- cumsum(unlist(colsep));  
  
  pairs.breaks <- seq(-4, 4, length.out=101);
  
  pdf(file = paste("./", n, "/", Sys.Date(), "_marker_heatmap_", 
                   length(unique(GB$predicted_cell_type[i])),"_", n, "_predicted_cell_type.pdf", 
                   sep = ""), 
      width = 20, height = 15, useDingbats = F);
  heatmap.3(temp[, ord],
            breaks = pairs.breaks,
            #symbreaks = T,
            keysize = 0.8, 
            main="Marker gene expression",
            col = bluered(100),
            symkey = F,
            cexRow=1, cexCol = 0.6, 
            Rowv = F, 
            Colv = F,#as.dendrogram(hca), 
            ColSideColors = as.matrix(ColSideColors),
            ColSideColorsSize = 2,
            #RowSideColors = RowSideColors,
            dendrogram = "none",
            scale = "row",
            colsep = colsep,
            sepcolor = "black",
            labRow = substr(rownames(temp), start = 20, stop = 100), #c("Glast", "Tbr2"),
            labCol = "",
            na.rm = F);
  dev.off(dev.cur());
}

#### 7.2.3 Do spec analysis on E14.5 dMGE ####
spec_scores <- read.csv("2017-06-23_E12.5_spec_scores_all_clusters.csv", header = T)

spec_scores <- read.csv("2017-06-26_E14.5_spec_scores_all_clusters.csv", header = T)

rownames(spec_scores) <- spec_scores$X
spec_scores <- spec_scores[,-1]

spec_genes <- c(order(spec_scores[1,], decreasing = T)[1:50], 
             order(spec_scores[2,], decreasing = T)[1:50],
             order(spec_scores[3,], decreasing = T)[1:50])
             
keep <- GB$time_point %in% "E14.5"
GB.keep <- GB[keep,]
temp <- Y[unique(spec_genes), keep]

ord <- order(GB.keep$region)

ColSideColors <- c(gg_color_hue(length(unique(GB.keep$IFC)))[as.numeric(as.factor(GB.keep$IFC))],
                   rainbow(length(unique(GB.keep$region)))[as.numeric(as.factor(GB.keep$region))])

ColSideColors <- matrix(data = ColSideColors, nrow = ncol(temp));
ColSideColors <- ColSideColors[ord, ];

colsep <- lapply(unique(GB.keep$region[ord]), 
                 function(x){length(which(GB.keep$region == x))});
colsep <- cumsum(unlist(colsep));  

pairs.breaks <- seq(-4, 4, length.out=101);

pdf(file = paste("./images/", Sys.Date(), "_spec_heatmap_E14.5.pdf", sep = ""), 
    width = 10, height = 15, useDingbats = F);
heatmap.3(temp[, ord],
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.5, 
          main="Marker gene expression",
          col = bluered(100),
          symkey = F,
          cexRow=1, cexCol = 1, 
          Rowv = F, 
          Colv = F,#as.dendrogram(hca), 
          ColSideColors = as.matrix(ColSideColors),
          ColSideColorsSize = 2,
          #RowSideColors = RowSideColors,
          dendrogram = "none",
          scale = "row",
          colsep = colsep,
          sepcolor = "black",
          labRow = substr(rownames(temp), start = 20, stop = 100), #c("Glast", "Tbr2"),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

#### 7.3.1 Do pcaReduce on E14.5 vMGE cells ####
vMGE_14.5 <- GB$region == "vMGE" & GB$time_point == "E14.5" #& GB$year == "2016"
GB.vMGE_14.5 <- GB[vMGE_14.5, ]

#seed <- 3
#set.seed(seed)

#nbt = 1
#q = 15
#Output_S <- PCAreduce(t(Data[genes, ]), nbt = nbt, q = q, method = 'S')
#load("2017-04-04_PCA_genes_E14.5_vMGE_Output_S.Rdata")
load("E14.5_vMGE/2017-04-30_E14.5_vMGE_Output_S.Rdata") # genes selected by top 3 PCs

GB.vMGE_14.5$clusters <- as.factor(Output_S[[1]][,9]) #13 clusters

set.seed(1024)
p = 20
m = 0.9
t_sne <- Rtsne(t(Data), initial_dims = 12, perplexity = p, pca = T, verbose = F, theta = 0, 
               momentum = m, final_momentum = 0.1, max_iter = 1000)$Y
colnames(t_sne) <- c("X1", "X2")
rownames(t_sne) <- GB.vMGE_14.5$cell_ID
t_sne <- data.frame(t_sne)

pdf(file = paste("./E14.5_vMGE/", Sys.Date(), "_E14.5_vMGE_tSNE_perplexity_", p, "_momentum_", 
                 m , ".pdf", sep = ""), 
    width = 12, height = 10, useDingbats = F)
ggplot(data = t_sne, aes(x = X1, y = X2)) +
  geom_point(size = 5, aes(x = X1, y = X2, color = GB.vMGE_14.5$clusters)) +
  #scale_color_brewer(palette = "Paired", name = 'Ages') +
  scale_color_discrete(name = "Cell Types") +
  #geom_text(aes(label = colnames(Y2)), hjust = 1, vjust = -1) + 
  labs(x = "tSNE 1", y = "tSNE 2") +
  ggtitle(paste("All cells, Perplexity =", p)) +
  theme_bw() +
  theme(text = element_text(size = 20),
        line = element_line(size = 1),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.key = element_blank(),
        legend.background = element_blank())
dev.off(which = dev.cur())

## Check a few marker genes

GB.vMGE_14.5$clusters <- Output_S[[1]][,9] #13 clusters
GB.vMGE_14.5$cell_type <- "T1"
GB.vMGE_14.5$cell_type[ GB.vMGE_14.5$clusters == 2 ] <- "N1"
GB.vMGE_14.5$cell_type[ GB.vMGE_14.5$clusters == 3 ] <- "T2"
GB.vMGE_14.5$cell_type[ GB.vMGE_14.5$clusters == 4 ] <- "T3"
GB.vMGE_14.5$cell_type[ GB.vMGE_14.5$clusters == 5 ] <- "N2"
GB.vMGE_14.5$cell_type[ GB.vMGE_14.5$clusters == 6 ] <- "N3"
GB.vMGE_14.5$cell_type[ GB.vMGE_14.5$clusters == 7 ] <- "P1"
GB.vMGE_14.5$cell_type[ GB.vMGE_14.5$clusters == 8 ] <- "T4"
GB.vMGE_14.5$cell_type[ GB.vMGE_14.5$clusters == 9 ] <- "T5"
GB.vMGE_14.5$cell_type[ GB.vMGE_14.5$clusters == 10 ] <- "P2"
GB.vMGE_14.5$cell_type[ GB.vMGE_14.5$clusters == 11 ] <- "P3"
GB.vMGE_14.5$cell_type[ GB.vMGE_14.5$clusters == 12 ] <- "N4"
GB.vMGE_14.5$cell_type[ GB.vMGE_14.5$clusters == 13 ] <- "P4"

GB.vMGE_14.5$cell_type <- factor(GB.vMGE_14.5$cell_type, 
                                 levels = c("P1", "P2","P3", "P4", "T1", "T2", "T3", "T4",
                                            "T5","N1", "N2", "N3", "N4"))

markers <- read.csv("20170323_Rubenstein_revised_gene_list.csv", header = F, stringsAsFactors = F)
markers <- unlist(markers[ ,1 ])
markers <- unique(sapply(markers, function(x){paste("\\b", x, "$", sep = "")}))
markers <- unlist(sapply(markers,function(x){grep(x, rownames(Y), value = T)}))

temp <- Y[ markers, vMGE_14.5]
ord <- order(GB.vMGE_14.5$cell_type, GB.vMGE_14.5$region, GB.vMGE_14.5$time_point)
ColSideColors <- c(gg_color_hue(12)[as.numeric(as.factor(GB.vMGE_14.5$cell_type))],
                 rainbow(length(unique(GB.vMGE_14.5$IFC)))[as.numeric(as.factor(GB.vMGE_14.5$IFC))])

ColSideColors <- matrix(data = ColSideColors, nrow = ncol(temp));
ColSideColors <- ColSideColors[ord, ];

colsep <- lapply(unique(GB.vMGE_14.5$cell_type[ord]), 
                 function(x){length(which(GB.vMGE_14.5$cell_type == x))});
colsep <- cumsum(unlist(colsep));  

pairs.breaks <- seq(-4, 4, length.out=101);

pdf(file = paste("./E14.5_vMGE/", Sys.Date(), "_marker_heatmap_", 
                 length(unique(GB.vMGE_14.5$cell_type)),"_vMGE_14.5_cell_type.pdf", 
                 sep = ""), 
    width = 20, height = 15, useDingbats = F);
heatmap.3(temp[, ord],
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8, 
          main="Marker gene expression",
          col = bluered(100),
          symkey = F,
          cexRow=1, cexCol = 0.6, 
          Rowv = F, 
          Colv = F,#as.dendrogram(hca), 
          ColSideColors = as.matrix(ColSideColors),
          ColSideColorsSize = 2,
          #RowSideColors = RowSideColors,
          dendrogram = "none",
          scale = "row",
          colsep = colsep,
          sepcolor = "black",
          labRow = substr(rownames(temp), start = 20, stop = 100), #c("Glast", "Tbr2"),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

write.csv(GB.vMGE_14.5, 
          file = paste(Sys.Date(), "vMGE_14.5_cluster_information.csv", sep = "_"))

## mav, pax6,  

#### Do monocle on E14.5 vMGE cells ####
expr <- Y[ , colnames(Data)]
pd <- new("AnnotatedDataFrame", data = GB.vMGE_14.5)
rownames(pd) <- GB.vMGE_14.5$cell_ID
fd <- data.frame(gene_ID = rownames(expr), 
                 gene_short_name = substr(rownames(expr), start = 20, stop = 100))
rownames(fd) <- fd$gene_ID
fd <- new("AnnotatedDataFrame", data = fd)
mm <- newCellDataSet(as(as.matrix(expr), "sparseMatrix"),
                     phenoData = pd, 
                     featureData = fd,
                     lowerDetectionLimit = 0.1, 
                     expressionFamily = tobit(Lower = 0.1))


## convert RPKM to RPC
#rpc_matrix <- relative2abs(mm)

# mm <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
#                      phenoData = pd, featureData = fd,
#                      lowerDetectionLimit = 0.5,
#                      expressionFamily = negbinomial.size())

## clean up
rm(expr, rpc_matrix, pd, fd)

## calculate size factor and dispersions
mm <- estimateSizeFactors(mm)
mm <- estimateDispersions(mm)

## filtering low-quality cells
mm <- detectGenes(mm, min_expr = 0.1)
#head(fData(mm))

expr_genes <- row.names(subset(fData(mm), num_cells_expressed >= 10))

## take a look at the distribution of mRNA totals across the cells
pData(mm)$Total_mRNAs <- colSums(exprs(mm))
mm <- mm[,pData(mm)$Total_mRNAs < 1e6]
upper_bound <- 10^(mean(log10(pData(mm)$Total_mRNAs)) + 2*sd(log10(pData(mm)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(mm)$Total_mRNAs)) - 2*sd(log10(pData(mm)$Total_mRNAs)))

pdf(file = paste("./E14.5_vMGE/", Sys.Date(), "_E14.5_vMGE_mRNA_totals.pdf", sep = ""), 
    width = 10, height = 10, useDingbats = F)
qplot(Total_mRNAs, data=pData(mm), color=time_point, geom="density") +
  geom_vline(xintercept=lower_bound) +
  geom_vline(xintercept=upper_bound)
dev.off()

mm <- mm[,pData(mm)$Total_mRNAs > lower_bound &
           pData(mm)$Total_mRNAs < upper_bound]
mm <- detectGenes(mm, min_expr = 0.1)

mm <- reduceDimension(mm, max_components = 2, num_dim = 11, reduction_method = 'tSNE',
                      verbose = T)
pData(mm)$cell_type <- GB.vMGE_14.5$cell_type[GB.vMGE_14.5$cell_ID %in% pData(mm)$cell_ID]

mm <- clusterCells(mm, num_clusters = 12)

pdf(file = paste("./E14.5_vMGE/", Sys.Date(), "_E14.5_vMGE_tSNE_cell_type.pdf", sep = ""), 
    width = 10, height = 10, useDingbats = F)
plot_cell_clusters(mm, 1, 2, color = "cell_type", cell_size = 3)
dev.off(which = dev.cur())

ordering_genes <- rownames(Data)

mm <- setOrderingFilter(mm, ordering_genes = ordering_genes)

mm <- reduceDimension(mm, max_components = 2, verbose = T)

mm <- orderCells(mm)

## plot trajectory colored by cell type
pdf(file = paste("./E14.5_vMGE/", Sys.Date(), "_E14.5_vMGE_trajectory_cell_type.pdf", sep = ""), 
    width = 11, height = 10, useDingbats = F)
plot_cell_trajectory(mm, color_by="cell_type", cell_size = 3)
dev.off()

## plot trajectory colored by pseudotime
pdf(file = paste("./E14.5_vMGE/", Sys.Date(), "_E14.5_vMGE_trajectory_pseudotime.pdf", sep = ""), 
    width = 11, height = 10, useDingbats = F)
plot_cell_trajectory(mm, color_by="Pseudotime", cell_size = 3)
dev.off()

## plot trajectory colored by individual cell type
pdf(file = paste("./E14.5_vMGE/", Sys.Date(), "_E14.5_vMGE_trajectory_cell_type_wrap.pdf", sep = ""), 
    width = 11, height = 10, useDingbats = F)
plot_cell_trajectory(mm, color_by="cell_type") + facet_wrap(~cell_type, nrow=3)
dev.off()

##
cds <- mm[markers[1:12],]
pdf(file = paste("./E14.5_vMGE/", Sys.Date(), "_E14.5_vMGE_genes_in_pseudotime.pdf", sep = ""), 
    width = 11, height = 10, useDingbats = F)
plot_genes_in_pseudotime(cds, color_by="cell_type", ncol = 3)
dev.off()

#### 7.3.2 Do spec analysis on E14.5 vMGE####
cell_type <- GB.vMGE_14.5$cell_type
Data <- scale(Y2)

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, outfile='LOG.TXT')
clusterExport(cl, c("cell_type", "Data"))

spec_scores <- parApply(cl = cl, Data, 1, function(x){
  source("~/Scripts/R/spec.R");
  gene <- as.numeric(x);
  opt_bin <- optimizebinsize(gene, header = cell_type);
  unlist(specy(gene, header = cell_type, binsize = opt_bin));
})
stopCluster(cl)

rnames <- c(unlist(lapply(levels(GB.vMGE_14.5$cell_type),function(x){paste("Spec score", x)})),
            unlist(lapply(levels(GB.vMGE_14.5$cell_type),function(x){paste("mean/bin size", x)})))
rownames(spec_scores) <- rnames
spec_scores <- data.frame(spec_scores)

write.csv(spec_scores, file = "E14.5_vMGE_spec_scores.csv")

#### 7.3.1 Do pcaReduce on E14.5 CGE cells ####
CGE_14.5 <- GB$region == "CGE" & GB$time_point == "E14.5" #& GB$year == "2016"
GB.CGE_14.5 <- GB[CGE_14.5, ]

#seed <- 3
#set.seed(seed)

#nbt = 1
#q = 15
#Output_S <- PCAreduce(t(Data[genes, ]), nbt = nbt, q = q, method = 'S')
load("2017-04-04_PCA_genes_E14.5_CGE_Output_S.Rdata")

## Check a few marker genes

GB.CGE_14.5$clusters <- Output_S[[1]][,15] #10 clusters
GB.CGE_14.5$cell_type <- "T1"
GB.CGE_14.5$cell_type[ GB.CGE_14.5$clusters == 2 ] <- "N1"
GB.CGE_14.5$cell_type[ GB.CGE_14.5$clusters == 3 ] <- "T2"
GB.CGE_14.5$cell_type[ GB.CGE_14.5$clusters == 4 ] <- "N2"
GB.CGE_14.5$cell_type[ GB.CGE_14.5$clusters == 5 ] <- "N3"
GB.CGE_14.5$cell_type[ GB.CGE_14.5$clusters == 6 ] <- "P1"
GB.CGE_14.5$cell_type[ GB.CGE_14.5$clusters == 7 ] <- "N4"


GB.CGE_14.5$cell_type <- factor(GB.CGE_14.5$cell_type, 
                                 levels = c("P1","T1", "T2", "N1", "N2", "N3", "N4"))

temp <- Y[ markers, CGE_14.5 ]
ord <- order(GB.CGE_14.5$cell_type, GB.CGE_14.5$region, GB.CGE_14.5$time_point)
n <- length(unique(GB.CGE_14.5$cell_type))
ColSideColors <- c(gg_color_hue(5)[as.numeric(as.factor(GB.CGE_14.5$IFC))],
                   rainbow(n)[as.numeric(GB.CGE_14.5$cell_type)])

ColSideColors <- matrix(data = ColSideColors, nrow = ncol(temp));
ColSideColors <- ColSideColors[ord, ];

colsep <- lapply(unique(GB.CGE_14.5$cell_type[ord]), 
                 function(x){length(which(GB.CGE_14.5$cell_type == x))});
colsep <- cumsum(unlist(colsep));  

pairs.breaks <- seq(-4, 4, length.out=101);

pdf(file = paste("./images/clusters/", Sys.Date(), "_marker_heatmap_", n,"_CGE_14.5_cell_type.pdf", sep = ""), 
    width = 20, height = 15, useDingbats = F);
heatmap.3(temp[, ord],
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8, 
          main="Marker gene expression",
          col = bluered(100),
          symkey = F,
          cexRow=1, cexCol = 0.6, 
          Rowv = F, 
          Colv = F,#as.dendrogram(hca), 
          ColSideColors = as.matrix(ColSideColors),
          ColSideColorsSize = 2,
          #RowSideColors = RowSideColors,
          dendrogram = "none",
          scale = "row",
          colsep = colsep,
          sepcolor = "black",
          labRow = substr(rownames(temp), start = 20, stop = 100), #c("Glast", "Tbr2"),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());
write.csv(GB.CGE_14.5, 
          file = paste(Sys.Date(), "CGE_14.5_cluster_information.csv", sep = "_"))

#### 7.3.2 Do spec analysis on E14.5 vMGE####
cell_type <- GB.vMGE_14.5$cell_type
Data <- scale(Y2)

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, outfile='LOG.TXT')
clusterExport(cl, c("cell_type", "Data"))

spec_scores <- parApply(cl = cl, Data, 1, function(x){
  source("~/Scripts/R/spec.R");
  gene <- as.numeric(x);
  opt_bin <- optimizebinsize(gene, header = cell_type);
  unlist(specy(gene, header = cell_type, binsize = opt_bin));
})
stopCluster(cl)

rnames <- c(unlist(lapply(levels(GB.vMGE_14.5$cell_type),function(x){paste("Spec score", x)})),
            unlist(lapply(levels(GB.vMGE_14.5$cell_type),function(x){paste("mean/bin size", x)})))
rownames(spec_scores) <- rnames
spec_scores <- data.frame(spec_scores)

write.csv(spec_scores, file = "E14.5_vMGE_spec_scores.csv")

#### Do pcaReduce on E12.5 dMGE cells ####
#### 7.3.1 Do pcaReduce on E12.5 vMGE cells ####
dMGE_12.5 <- GB$region == "dMGE" & GB$time_point == "E12.5" #& GB$year == "2016"
GB.dMGE_12.5 <- GB[dMGE_12.5, ]

#seed <- 3
#set.seed(seed)

#nbt = 1
#q = 15
#Output_S <- PCAreduce(t(Data[genes, ]), nbt = nbt, q = q, method = 'S')
load("2017-04-04_PCA_genes_E12.5_dMGE_Output_S.Rdata")

## Check a few marker genes

GB.dMGE_12.5$clusters <- Output_S[[1]][,10] #12 clusters
GB.dMGE_12.5$cell_type <- "T1"
GB.dMGE_12.5$cell_type[ GB.dMGE_12.5$clusters == 2 ] <- "N1"
GB.dMGE_12.5$cell_type[ GB.dMGE_12.5$clusters == 3 ] <- "T2"
GB.dMGE_12.5$cell_type[ GB.dMGE_12.5$clusters == 4 ] <- "P1"
GB.dMGE_12.5$cell_type[ GB.dMGE_12.5$clusters == 5 ] <- "P2"
GB.dMGE_12.5$cell_type[ GB.dMGE_12.5$clusters == 6 ] <- "T3"
GB.dMGE_12.5$cell_type[ GB.dMGE_12.5$clusters == 7 ] <- "N2"
GB.dMGE_12.5$cell_type[ GB.dMGE_12.5$clusters == 8 ] <- "N3"
GB.dMGE_12.5$cell_type[ GB.dMGE_12.5$clusters == 9 ] <- "P3"
GB.dMGE_12.5$cell_type[ GB.dMGE_12.5$clusters == 10 ] <- "N4"
GB.dMGE_12.5$cell_type[ GB.dMGE_12.5$clusters == 11 ] <- "P4"
GB.dMGE_12.5$cell_type[ GB.dMGE_12.5$clusters == 12 ] <- "T4"


GB.dMGE_12.5$cell_type <- factor(GB.dMGE_12.5$cell_type, 
                                levels = c("P1","P2", "P3", "P4","T1", "T2", "T3", "T4", 
                                           "N1", "N2", "N3", "N4"))

temp <- Y[ markers, dMGE_12.5]
ord <- order(GB.dMGE_12.5$cell_type, GB.dMGE_12.5$region, GB.dMGE_12.5$time_point)
n <- length(unique(GB.dMGE_12.5$cell_type))
ColSideColors <- c(gg_color_hue(5)[as.numeric(as.factor(GB.dMGE_12.5$IFC))],
                   rainbow(n)[as.numeric(GB.dMGE_12.5$cell_type)])

ColSideColors <- matrix(data = ColSideColors, nrow = ncol(temp));
ColSideColors <- ColSideColors[ord, ];

colsep <- lapply(unique(GB.dMGE_12.5$cell_type[ord]), 
                 function(x){length(which(GB.dMGE_12.5$cell_type == x))});
colsep <- cumsum(unlist(colsep));  

pairs.breaks <- seq(-4, 4, length.out=101);

pdf(file = paste("./images/clusters/", Sys.Date(), "_marker_heatmap_", n,"_dMGE_12.5_cell_type.pdf", sep = ""), 
    width = 20, height = 15, useDingbats = F);
heatmap.3(temp[, ord],
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8, 
          main="Marker gene expression",
          col = bluered(100),
          symkey = F,
          cexRow=1, cexCol = 0.6, 
          Rowv = F, 
          Colv = F,#as.dendrogram(hca), 
          ColSideColors = as.matrix(ColSideColors),
          ColSideColorsSize = 2,
          #RowSideColors = RowSideColors,
          dendrogram = "none",
          scale = "row",
          colsep = colsep,
          sepcolor = "black",
          labRow = substr(rownames(temp), start = 20, stop = 100), #c("Glast", "Tbr2"),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());
write.csv(GB.dMGE_12.5, 
          file = paste(Sys.Date(), "dMGE_12.5_cluster_information.csv", sep = "_"))

#### Do pcaReduce on E12.5 vMGE cells ####
#### 7.3.1 Do pcaReduce on E12.5 vMGE cells ####
vMGE_12.5 <- GB$region == "vMGE" & GB$time_point == "E12.5" #& GB$year == "2016"
GB.vMGE_12.5 <- GB[vMGE_12.5, ]

#seed <- 3
#set.seed(seed)

#nbt = 1
#q = 15
#Output_S <- PCAreduce(t(Data[genes, ]), nbt = nbt, q = q, method = 'S')
load("2017-04-04_PCA_genes_E12.5_vMGE_Output_S.Rdata")

## Check a few marker genes

GB.vMGE_12.5$clusters <- Output_S[[1]][,10] #12 clusters
GB.vMGE_12.5$cell_type <- "P3"
GB.vMGE_12.5$cell_type[ GB.vMGE_12.5$clusters == 2 ] <- "P4"
GB.vMGE_12.5$cell_type[ GB.vMGE_12.5$clusters == 3 ] <- "P5"
GB.vMGE_12.5$cell_type[ GB.vMGE_12.5$clusters == 4 ] <- "P1"
GB.vMGE_12.5$cell_type[ GB.vMGE_12.5$clusters == 5 ] <- "T1"
GB.vMGE_12.5$cell_type[ GB.vMGE_12.5$clusters == 6 ] <- "T2"
GB.vMGE_12.5$cell_type[ GB.vMGE_12.5$clusters == 7 ] <- "N1"
GB.vMGE_12.5$cell_type[ GB.vMGE_12.5$clusters == 8 ] <- "T3"
GB.vMGE_12.5$cell_type[ GB.vMGE_12.5$clusters == 9 ] <- "P2"
GB.vMGE_12.5$cell_type[ GB.vMGE_12.5$clusters == 10 ] <- "N2"
GB.vMGE_12.5$cell_type[ GB.vMGE_12.5$clusters == 11 ] <- "N3"
GB.vMGE_12.5$cell_type[ GB.vMGE_12.5$clusters == 12 ] <- "T4"


GB.vMGE_12.5$cell_type <- factor(GB.vMGE_12.5$cell_type, 
                                 levels = c("P1","P2", "P3","P4","T1", "T2", "T3", "T4", 
                                            "N1", "N2", "N3"))

temp <- Y[ markers, vMGE_12.5]
ord <- order(GB.vMGE_12.5$cell_type, GB.vMGE_12.5$region, GB.vMGE_12.5$time_point)
n <- length(unique(GB.vMGE_12.5$cell_type))
ColSideColors <- c(gg_color_hue(5)[as.numeric(as.factor(GB.vMGE_12.5$IFC))],
                   rainbow(n)[as.numeric(GB.vMGE_12.5$cell_type)])

ColSideColors <- matrix(data = ColSideColors, nrow = ncol(temp));
ColSideColors <- ColSideColors[ord, ];

colsep <- lapply(unique(GB.vMGE_12.5$cell_type[ord]), 
                 function(x){length(which(GB.vMGE_12.5$cell_type == x))});
colsep <- cumsum(unlist(colsep));  

pairs.breaks <- seq(-4, 4, length.out=101);

pdf(file = paste("./images/clusters/", Sys.Date(), "_marker_heatmap_", n,"_vMGE_12.5_cell_type.pdf", sep = ""), 
    width = 20, height = 15, useDingbats = F);
heatmap.3(temp[, ord],
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8, 
          main="Marker gene expression",
          col = bluered(100),
          symkey = F,
          cexRow=1, cexCol = 0.6, 
          Rowv = F, 
          Colv = F,#as.dendrogram(hca), 
          ColSideColors = as.matrix(ColSideColors),
          ColSideColorsSize = 2,
          #RowSideColors = RowSideColors,
          dendrogram = "none",
          scale = "row",
          colsep = colsep,
          sepcolor = "black",
          labRow = substr(rownames(temp), start = 20, stop = 100), #c("Glast", "Tbr2"),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());
write.csv(GB.vMGE_12.5, 
          file = paste(Sys.Date(), "vMGE_12.5_cluster_information.csv", sep = "_"))

#### Do pcaReduce on E12.5 CGE cells ####
CGE_12.5 <- GB$region == "CGE" & GB$time_point == "E12.5" #& GB$year == "2016"
GB.CGE_12.5 <- GB[CGE_12.5, ]

#seed <- 3
#set.seed(seed)

#nbt = 1
#q = 15
#Output_S <- PCAreduce(t(Data[genes, ]), nbt = nbt, q = q, method = 'S')
load("2017-04-04_PCA_genes_E12.5_CGE_Output_S.Rdata")

## Check a few marker genes

GB.CGE_12.5$clusters <- Output_S[[1]][,13] #9 clusters
GB.CGE_12.5$cell_type <- "T1"
GB.CGE_12.5$cell_type[ GB.CGE_12.5$clusters == 2 ] <- "N1"
GB.CGE_12.5$cell_type[ GB.CGE_12.5$clusters == 3 ] <- "N2"
GB.CGE_12.5$cell_type[ GB.CGE_12.5$clusters == 4 ] <- "T2"
GB.CGE_12.5$cell_type[ GB.CGE_12.5$clusters == 5 ] <- "P1"
GB.CGE_12.5$cell_type[ GB.CGE_12.5$clusters == 6 ] <- "P2"
GB.CGE_12.5$cell_type[ GB.CGE_12.5$clusters == 7 ] <- "T3"
GB.CGE_12.5$cell_type[ GB.CGE_12.5$clusters == 8 ] <- "T4"
GB.CGE_12.5$cell_type[ GB.CGE_12.5$clusters == 9 ] <- "N3"

GB.CGE_12.5$cell_type <- factor(GB.CGE_12.5$cell_type, 
                                 levels = c("P1","P2", "T1", "T2","T3","T4","N1","N2","N3"))

temp <- Y[ markers, CGE_12.5]
ord <- order(GB.CGE_12.5$cell_type, GB.CGE_12.5$region, GB.CGE_12.5$time_point)
n <- length(unique(GB.CGE_12.5$cell_type))
ColSideColors <- c(gg_color_hue(unique(length(GB.CGE_12.5$IFC)))[as.numeric(as.factor(GB.CGE_12.5$IFC))],
                   rainbow(n)[as.numeric(GB.CGE_12.5$cell_type)])

ColSideColors <- matrix(data = ColSideColors, nrow = ncol(temp));
ColSideColors <- ColSideColors[ord, ];

colsep <- lapply(unique(GB.CGE_12.5$cell_type[ord]), 
                 function(x){length(which(GB.CGE_12.5$cell_type == x))});
colsep <- cumsum(unlist(colsep));  

pairs.breaks <- seq(-4, 4, length.out=101);

pdf(file = paste("./images/clusters/", Sys.Date(), "_marker_heatmap_", n,"_CGE_12.5_cell_type.pdf", sep = ""), 
    width = 20, height = 15, useDingbats = F);
heatmap.3(temp[, ord],
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8, 
          main="Marker gene expression",
          col = bluered(100),
          symkey = F,
          cexRow=1, cexCol = 0.6, 
          Rowv = F, 
          Colv = F,#as.dendrogram(hca), 
          ColSideColors = as.matrix(ColSideColors),
          ColSideColorsSize = 2,
          #RowSideColors = RowSideColors,
          dendrogram = "none",
          scale = "row",
          colsep = colsep,
          sepcolor = "black",
          labRow = substr(rownames(temp), start = 20, stop = 100), #c("Glast", "Tbr2"),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());
write.csv(GB.CGE_12.5, 
          file = paste(Sys.Date(), "CGE_12.5_cluster_information.csv", sep = "_"))

## Clean up
rm(panRGC, aRGC, bRGC, IPC, EENR, pNeuron,EEN,LEN,IN,Astro,Olig,micro,vent, endo, 
   markers, x, ord, n, ColSideColors, colsep, temp, pairs.breaks)

#### 8. Do ANOVA ####
## Do ANOVA for all the cells 
data <- log2(Y + 1)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, outfile='LOG.TXT')
clusterExport(cl, c("GB"))

result <- parApply(cl, data, 1, function(x){
  df <- data.frame(x = x, time_point = GB$time_point, region = GB$region)
  fit <- aov(x ~ time_point*region, data = df)
  summary(fit)[[1]]$`Pr(>F)`[1:3]
})

stopCluster(cl)

## Do ANOVA for all the MGE cells in E12.5 period
MGE_12.5 <- GB$region != "CGE" & GB$time_point == "E12.5"
data <- log2(Y[,MGE_12.5] + 1)
GB.res <- GB[MGE_12.5,]

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, outfile='LOG.TXT')
clusterExport(cl, c("GB.res"))

result <- parApply(cl, data, 1, function(x){
  df <- data.frame(x = x, region = GB.res$region)
  fit <- aov(x ~ region, data = df)
  summary(fit)[[1]]$`Pr(>F)`[1]
})

stopCluster(cl)

## Plot the significant genes according to time and region interaction
genes <- names(sort(result)[1:100])

temp <- Y[ genes, MGE_12.5 ]

ColSideColors <- c(gg_color_hue(2)[as.numeric(as.factor(GB$time_point))],
                   rainbow(23)[as.numeric(as.factor(GB$IFC))],
                   brewer.pal(3, "Paired")[as.numeric(as.factor(GB$region))])

ColSideColors <- matrix(data = ColSideColors, nrow = ncol(Y));
ord <- order(as.factor(GB$time_point), as.factor(GB$region))
ColSideColors <- ColSideColors;

pairs.breaks <- seq(-4, 4, length.out=101);

pdf(file = paste("./images/", Sys.Date(), "_ANOVA_heatmap_MGE_E12.5.pdf", sep = ""), 
    width = 20, height = 15, useDingbats = F);
heatmap.3(temp,
          breaks = pairs.breaks,
          #symbreaks = T,
          keysize = 0.8, 
          main="Marker gene expression",
          col = bluered(100),
          symkey = F,
          cexRow=1, cexCol = 0.6, 
          Rowv = T, 
          Colv = T,#as.dendrogram(hca), 
          ColSideColors = as.matrix(ColSideColors)[ MGE_12.5, ],
          ColSideColorsSize = 2,
          #RowSideColors = RowSideColors,
          dendrogram = "none",
          scale = "row",
          #colsep = colsep,
          #sepcolor = "black",
          labRow = substr(rownames(temp), start = 20, stop = 100), #c("Glast", "Tbr2"),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

## Do ANOVA for all the MGE cells in E14.5 period
MGE_14.5 <- GB$region != "CGE" & GB$time_point == "E14.5"
data <- log2(Y[ ,MGE_14.5 ] + 1)
GB.res <- GB[ MGE_14.5, ]

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, outfile='LOG.TXT')
clusterExport(cl, c("GB.res"))

aov.result <- parApply(cl, data, 1, function(x){
  df <- data.frame(x = x, region = GB.res$region)
  fit <- aov(x ~ region, data = df)
  summary(fit)[[1]]$`Pr(>F)`[1]
})

#cor.result <- parApply(cl, data, 1, function(x){
#  cor(x, as.numeric(as.factor(GB.res$region)))
#})

stopCluster(cl)

## Plot the significant genes according to time and region interaction
genes <- names(sort(aov.result)[1:100])
#genes <- names(sort(cor.result, decreasing = T)[1:100])

temp <- log2(Y[ genes, MGE_14.5] + 1)

ColSideColors <- c(gg_color_hue(2)[as.numeric(as.factor(GB.res$time_point))],
                   rainbow(23)[as.numeric(as.factor(GB.res$IFC))],
                   brewer.pal(3, "Paired")[as.numeric(as.factor(GB.res$region))])
ColSideColors <- matrix(data = ColSideColors, nrow = ncol(temp));

ord <- order(as.factor(GB.res$region))

pairs.breaks <- seq(-4, 4, length.out=101);

pdf(file = paste("./images/", Sys.Date(), "_ANOVA_heatmap_MGE_E14.5.pdf", sep = ""), 
    width = 20, height = 15, useDingbats = F);
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
          ColSideColors = ColSideColors[ord,],
          ColSideColorsSize = 2,
          #RowSideColors = RowSideColors,
          dendrogram = "none",
          scale = "row",
          #colsep = colsep,
          #sepcolor = "black",
          labRow = substr(rownames(temp), start = 20, stop = 100), #c("Glast", "Tbr2"),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());

## Clean up
rm(MGE, cl, ord)

## Do ANOVA for all the genes in MGE ##
MGE <- GB$region != "CGE"
data <- log2(Y[,MGE] + 1)
GB.MGE <- GB[MGE,]
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, outfile='LOG.TXT')
clusterExport(cl, c("GB.MGE"))

result <- parApply(cl, data, 1, function(x){
  df <- data.frame(x = x, time_point = GB.MGE$time_point, region = GB.MGE$region)
  fit <- aov(x ~ time_point*region, data = df)
  summary(fit)[[1]]$`Pr(>F)`[1:3]
})

stopCluster(cl)

#### Estimate number of clusters with K-means ####
gap.kmeans <- data.frame(clusGap(t_sne, FUNcluster = kmeans, K.max = 30, B = 1000)$Tab)
gap_elbow <- vector()
for(i in 1:nrow(gap.kmeans)){
  gap_elbow[i] <- (gap.kmeans$gap[i] - gap.kmeans$gap[i+1] + gap.kmeans$SE.sim[i+1])
}
gap.kmeans$gap_elbow <- gap_elbow

pdf(file = paste("./images/", Sys.Date(), "_da_mi_gap_elbow.pdf", sep = ""), 
    width = 10, height = 10, useDingbats = F)
ggplot(data = gap.kmeans, aes(x = 1:30, y = gap_elbow, group = F)) + 
  geom_point() + 
  geom_line(color = "black") +
  stat_smooth() +
  ggtitle(paste('Gap')) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2),
        panel.background =  element_rect(fill = "white", colour = NA), 
        panel.border = element_rect(fill = NA, colour="grey50"), 
        panel.grid.major =  element_line(colour = "grey90", size = 0.2),
        panel.grid.minor =  element_line(colour = "grey98", size = 0.5),
        axis.text.x=element_text(angle=90, size=10),
        axis.title.x = element_blank())
dev.off(which=dev.cur())

k = 10
kmeans.clusters <- kmeans(t_sne, k, iter.max = 100000, nstart = 100)
GB$kmeans <- kmeans.clusters$cluster

#### Define cell types ####
GB$cell_type <- "T4"
GB$cell_type[GB$kmeans == 2] <- "N2"
GB$cell_type[GB$kmeans == 3] <- "P2"
GB$cell_type[GB$kmeans == "4"] <- "N3"
GB$cell_type[GB$kmeans == "5"] <- "P3"
GB$cell_type[GB$kmeans == "6"] <- "P1"
GB$cell_type[GB$kmeans == "7"] <- "T3"
GB$cell_type[GB$kmeans == "8"] <- "T2"
GB$cell_type[GB$kmeans == "9"] <- "N1"
GB$cell_type[GB$kmeans == "10"] <- "T1"
GB$cell_type <- factor(GB$cell_type, levels = c("P1", "P2", "P3", "T1", "T2", "T3", "T4","N1", "N2", "N3"))

pdf(file = paste("./images/", Sys.Date(), "_Mida_tSNE_kmeans_", p, ".pdf", sep = ""), 
    width = 10, height = 10, useDingbats = F)
ggplot(data = t_sne, aes(x = X1, y = X2)) +
  geom_point(size = 3, aes(x = X1, y = X2, color = as.character(GB$kmeans), shape = GB$region)) +
  #scale_color_brewer(palette = "Paired", name = 'Ages') +
  scale_color_discrete(name = "cell type") +
  #geom_text(aes(label = colnames(Y2)), hjust = 1, vjust = -1) + 
  labs(x = "tSNE 1", y = "tSNE 2") +
  ggtitle(paste("All cells, Perplexity =", p)) +
  theme_bw() +
  theme(text = element_text(size = 20),
        line = element_line(size = 1),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.key = element_blank(),
        legend.position = c(0.9, 0.2),
        legend.background = element_blank())
dev.off(which = dev.cur())

#### MST analysis ####
z <- data.frame(kmeans.clusters$centers)
rownames(z) <- paste0("t",1:nrow(z))
z$kmeans <- as.character(unique(GB$kmeans[order(GB$kmeans)]))
z$cell_type <- unique(GB$cell_type[order(GB$kmeans)])
m <- mst(dist(z[,1:2]))
w <- which(m!=0)
i <- as.vector(row(m))[w]
j <- as.vector(col(m))[w]

g <- data.frame(from.x = as.numeric(z[i,1]), from.y = as.numeric(z[i,2]), 
                to.x = as.numeric(z[j,1]), to.y = as.numeric(z[j,2]))

pdf(file = paste("./images/", Sys.Date(), "_Mida_tSNE_lineage_CDK6_", p, ".pdf", sep = ""), 
    width = 15, height = 15, useDingbats = F)
ggplot(data = t_sne, aes(x = X1, y = X2, colour = as.character(GB$kmeans))) +
  geom_point(alpha = 1, size = 5) +
  geom_segment(data = g, aes(x = from.x, xend = to.x, y = from.y, yend = to.y), 
               colour = "grey60", size = 2) +
  geom_point(data = data.frame(z), aes(x = X1, y = X2, colour = kmeans), 
             size = 30, alpha = 1) +
  geom_text(data = data.frame(z), aes(x = X1, y = X2, label = cell_type), color = "black") +
  
  scale_color_discrete(name = "Cell type") +
  #scale_color_continuous(high = "red", name = "TUBB3") +
  #geom_text(aes(label = colnames(Y2)), hjust = 1, vjust = -1) + 
  labs(x = "tSNE 1", y = "tSNE 2") +
  ggtitle(paste("Lineage, Perplexity =", p)) +
  theme_bw() +
  theme(text = element_text(size = 20),
        line = element_line(size = 1),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.key = element_blank(),
        legend.position = "",
        legend.background = element_blank())
dev.off(dev.cur())

#### Pseudotime reconstruction ####
z <- kmeans.clusters$centers
z <- z[order(z[,1]),]

rownames(z) <-paste0("t",rownames(z))
m <- mst(dist(z[,1:2]))

t.names <-names(which(colSums(m!=0)==1))[2] # There are two ends, then use the left most one.
t.names <- "t6"
for (i in 1:nrow(m)){
  t.names <-append(t.names,names(which(m[t.names[i],]==1))[which(names(which(m[t.names[i],]==1)) %notin% t.names)])
}

y2d <-t_sne[,1:2]
#y2d <-y2d[order(y2d[,1]),]
z2d <-z[,1:2]
z2d <-z2d[t.names,]

time_start.i <-0
updatethis.dist <-rep(Inf,nrow(y2d))
updatethis.time <-rep(0,nrow(y2d))
update.updown <-rep(0,nrow(y2d))
pseudotime.flow <-c(0)

for (i in 1:(nrow(z2d)-1)){
  
  # distance between this z2d.i and all y2d
  dot.dist.i <-apply(y2d,1,function(X){dist(rbind(X,z2d[i,]))})
  
  # distance between this z2d.i-z2d.i+1 segment and "insider" y2d
  inside_this_segment <- which(apply(y2d,1,function(X){inside_check.foo(z2d[i,],z2d[i+1,],X)}))
  seg.dist.i <-rep(Inf,nrow(y2d))
  seg.dist.i[inside_this_segment] <-apply(y2d,1,function(X){distance.foo(z2d[i,],z2d[i+1,],X)})[inside_this_segment]
  
  # intersect coordinate between this z2d.i-z2d.i+1 segment and all y2d
  intersect.i <-t(apply(y2d,1,function(X){intersect.foo(z2d[i,],z2d[i+1,],X)}))
  
  # this z2d.i-z2d.i+1 segment's unit vector
  seg_unit_vector <-unit_vector.foo(z2d[i,],z2d[i+1,])
  
  # UPDATE
  # 2. idx for the shortest distance at this round (either dot or seg)
  update.idx <-apply(cbind(dot.dist.i,seg.dist.i,updatethis.dist),1,which.min)
  # 3. update the pseudotime for y2ds with the short distance from the z2d.i
  updatethis.time[which(update.idx==1)] <-time_start.i
  # 4. update the pseudotime for y2ds with the short distance from the z2d.i-z2d.i+1 segment
  relative_cordinates <-t(apply(intersect.i[which(update.idx==2),],1,function(X){seg_unit_vector%*%(X-z2d[i,])}))
  updatethis.time[which(update.idx==2)] <-time_start.i + relative_cordinates
  # 1. update the shortest distance
  updatethis.dist <-apply(cbind(dot.dist.i,seg.dist.i,updatethis.dist),1,min)
  
  update.updown[which(update.idx==1)] <-c(apply(y2d,1,function(X){crossvec_direction(z2d[i,],z2d[i+1,],X)})*dot.dist.i)[which(update.idx==1)]
  update.updown[which(update.idx==2)] <-c(apply(y2d,1,function(X){crossvec_direction(z2d[i,],z2d[i+1,],X)})*seg.dist.i)[which(update.idx==2)]
  
  # update time for the next round
  time_start.i <-time_start.i + dist(rbind(z2d[i,],z2d[i+1,]))
  pseudotime.flow <-append(pseudotime.flow,time_start.i)
}

# For the y2ds that are closest to the starting z2d
i=1
dot.dist.i <-apply(y2d,1,function(X){dist(rbind(X,z2d[i,]))})
if (length(start.idx <-which(dot.dist.i <= updatethis.dist))>0){
  intersect.i <-t(apply(y2d,1,function(X){intersect.foo(z2d[i,],z2d[i+1,],X)}))
  seg_unit_vector <-unit_vector.foo(z2d[i,],z2d[i+1,])
  relative_cordinates <-0 + t(apply(intersect.i,1,function(X){seg_unit_vector %*% (X-z2d[i,])}))[start.idx]
  updatethis.time[start.idx] <-relative_cordinates
  seg.dist.i <-apply(y2d,1,function(X){distance.foo(z2d[i,],z2d[i+1,],X)})
  update.updown[start.idx] <-c(apply(y2d,1,function(X){crossvec_direction(z2d[i,],z2d[i+1,],X)})*seg.dist.i)[start.idx]
}
# For the y2ds that are closest to the arriving z2d
i=nrow(z2d)
dot.dist.i <-apply(y2d,1,function(X){dist(rbind(X,z2d[i,]))})
if (length(arrive.idx <-which(dot.dist.i <= updatethis.dist))>0){
  intersect.i <-t(apply(y2d,1,function(X){intersect.foo(z2d[i-1,],z2d[i,],X)}))
  seg_unit_vector <-unit_vector.foo(z2d[i-1,],z2d[i,])
  relative_cordinates <-time_start.i + as.numeric(t(apply(intersect.i,1,function(X){seg_unit_vector %*% (X-z2d[i,])})))[arrive.idx]
  updatethis.time[arrive.idx] <-relative_cordinates
  seg.dist.i <-apply(y2d,1,function(X){distance.foo(z2d[i-1,],z2d[i,],X)})
  update.updown[arrive.idx] <-c(apply(y2d,1,function(X){crossvec_direction(z2d[i-1,],z2d[i,],X)})*seg.dist.i)[arrive.idx]
}

pseudotime <-updatethis.time
pseudotime.y <-update.updown
pseudotime.flow <-pseudotime.flow

if (x.reverse==T){
  pseudotime <- -pseudotime
  pseudotime.flow <- -pseudotime.flow
}

pseudotime_range <-max(pseudotime)-min(pseudotime)

pseudotime.flow <-pseudotime.flow-min(pseudotime)
pseudotime.flow <-pseudotime.flow/pseudotime_range

pseudotime <-pseudotime-min(pseudotime)
pseudotime <-pseudotime/pseudotime_range

plot(pseudotime, pseudotime.y,cex=3,pch=20,bty="n")
points(pseudotime.flow,rep(0,length(pseudotime.flow)), col="#FF0000", pch=19,cex=1)
segments(pseudotime.flow[1], 0, pseudotime.flow[length(pseudotime.flow)], 0, col="#FF0000" ,lwd=5)
#    text(pseudotime, pseudotime.y,labels=rownames(y2d))
#box()

#### Plot genes in pseudotime ####
pseudotime.df <- data.frame(t_sne[,1:2],pseudotime)
genes <- read.csv("Pseudotemporal_ordering.csv", header = F, stringsAsFactors = F)[,1]
genes <- unique(sapply(genes, function(x){paste("\\b", x, "$", sep = "")}))
genes <-unlist(lapply(genes, function(x){grep(x,rownames(Y))}))
pseudotime.df$cluster <- GB$kmeans
pseudotime.df <- cbind(pseudotime.df, t(Y[genes,]))
pseudotime.df <- melt(pseudotime.df, id.vars = c("X1", "X2", "pseudotime", "cluster"))

pdf(file = paste("./images/", Sys.Date(), "_pseudotime_genes.pdf", sep = ""), 
    width = 36, height = 50, useDingbats = F)
ggplot(data = pseudotime.df, aes(x = pseudotime, y = value, colour = as.character(cluster), group = 1)) +
  geom_point() +
  scale_color_discrete(name = "Cell type") +
  ylab(label = "Gene expression level") +
  geom_smooth(aes(group = 1)) +
  facet_wrap(~variable, nrow = 7) +
  ggtitle("Gene expression in pseudotime") +
  theme(plot.title = element_text(size=60, face="bold", vjust=2),
        
        panel.background =  element_rect(fill = "white", colour = NA), 
        panel.border = element_rect(fill = NA, colour="grey50"), 
        panel.grid.major =  element_line(colour = "grey90", size = 0.2),
        panel.grid.minor =  element_line(colour = "grey98", size = 0.5),
        strip.text.x = element_text(size = 50),
        axis.title.x = element_blank())
dev.off()

#### Plot marker genes from Linnarsson ####
temp <- na.omit(Y[unlist(markers),])

pairs.breaks <- seq(-2, 2, length.out=101);

ColSideColors <- c(brewer.pal(3, "Paired")[as.numeric(as.factor(GB$region))],
                   gg_color_hue(length(unique(GB$IFC)))[as.numeric(as.factor(GB$IFC))])
ColSideColors <- matrix(data = ColSideColors, nrow = ncol(temp))

pdf(file = paste("./images/", Sys.Date(), "_heatmap_linnarsson_markers.pdf", sep = ""), 
    width = 20, height = 30, useDingbats = F)
heatmap.3(temp,
          breaks = pairs.breaks, 
          keysize = 0.8, 
          main="Marker gene expression",
          col = bluered(100),
          symkey = F,
          cexRow=1, cexCol = 0.6, 
          Rowv = T, 
          Colv = T,#as.dendrogram(hca), 
          ColSideColors = as.matrix(ColSideColors),
          #ColSideColorsSize = 2,
          #RowSideColors = RowSideColors,
          dendrogram = "column",
          scale = "row",
          #colsep = colsep,
          sepcolor = "black",
          labRow = substr(rownames(temp), start = 17, stop = 100), #c("Glast", "Tbr2"),
          labCol = "",
          na.rm = F);
dev.off(dev.cur())

#### Gap analysis ####
gap.kmeans <- data.frame(clusGap(t_sne, FUNcluster = kmeans, K.max = 15, B = 1000)$Tab)

gap_elbow <- vector()
for(i in 1:nrow(gap.kmeans)){
  gap_elbow[i] <- (gap.kmeans$gap[i] - gap.kmeans$gap[i+1] + gap.kmeans$SE.sim[i+1])
}
gap.kmeans$gap_elbow <- gap_elbow

pdf(file = paste("./images/", Sys.Date(), "_Mida_gap_statistics.pdf", sep = ""), width = 11, height = 10, useDingbats = F)
ggplot(data = gap.kmeans, aes(x = 1:15, y =  gap, group = F)) + 
  geom_point() + 
  geom_line(color = "black") +
  stat_smooth() +
  ggtitle(paste('Gap statistics')) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2),
        panel.background =  element_rect(fill = "white", colour = NA), 
        panel.border = element_rect(fill = NA, colour="grey50"), 
        panel.grid.major =  element_line(colour = "grey90", size = 0.2),
        panel.grid.minor =  element_line(colour = "grey98", size = 0.5),
        axis.text.x=element_text(angle=90, size=10),
        axis.title.x = element_blank())
dev.off(dev.cur())

#### WGCNA ####
datExpr <- t(Y[,keep[[k]]])

# Choose a set of soft-thresholding powers
powers = c(seq(from = 1, to=20, by=1))

# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, 
                         powerVector = powers, 
                         verbose = 0,
                         corFnc = "bicor",
                         corOptions = list(maxPOutliers =0.1),
                         networkType = "signed hybrid")

# Scale-free topology fit index as a function of the soft-thresholding power

pdf(file = paste(Sys.Date(),name, "_sft.pdf", sep = ""), width = 11, height = 10, useDingbats = F)
ggplot(data = sft$fitIndices, aes(x = Power, y = -sign(slope)*SFT.R.sq)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  geom_hline(aes(yintercept = 0.9), colour = "red")
dev.off(dev.cur())

pdf(file = paste(Sys.Date(),name, "_mean_connectivity.pdf", sep = ""), width = 11, height = 10, useDingbats = F)
ggplot(data = sft$fitIndices, aes(x = Power, y = mean.k.)) +
  geom_point(size = 3) +
  geom_line(size = 0.5) +
  ggtitle(label = "Mean connectivity")
dev.off(dev.cur())
# Mean connectivity as a function of the soft-thresholding power
softPower = 5
cor.Y <- bicor(datExpr, use = "pairwise.complete.obs", maxPOutliers = 0.1)
adj = adjacency.fromSimilarity(cor.Y, type = "signed hybrid", power = softPower)
TOM = TOMsimilarity(adj, TOMDenom = "mean", TOMType = "signed")
colnames(TOM) <- rownames(TOM) <- rownames(Y)

#load("2016-02-19_iPSC_TOM_sfp_5.Rdata")
dissTOM <- 1 - TOM

geneTree <- hclust(as.dist(dissTOM),method="average");

# Set the minimum module size
minModuleSize = 20;

# Module identification using dynamic tree cut

dynamicMods = cutreeDynamic(dendro = geneTree,
                            method="hybrid",
                            deepSplit = T,
                            #minAbsSplitHeight = 0.999,
                            minClusterSize = minModuleSize,
                            distM = dissTOM);

#Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)
#table(dynamicColors)

pdf(file = paste(Sys.Date(),name, "_WGCNA_dentro.pdf", sep = ""), width = 11, height = 10, useDingbats = F)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, 
                    hang = 0.03, 
                    addGuide = TRUE, 
                    guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")
dev.off()

restGenes= (dynamicColors != "grey")
cor.Y.rest <- bicor(datExpr[,restGenes], maxPOutliers = 0.1)
adj2 <- adjacency.fromSimilarity(cor.Y.rest, type = "signed hybrid", power = softPower)
#adj2 <- adjacency(datExpr[,restGenes], type = "signed hybrid", power = softPower)
diss1=1-TOMsimilarity(adjMat =adj2, 
                      TOMType = "signed",
                      TOMDenom = "mean")

colnames(diss1) =rownames(diss1) =rownames(Y2)[restGenes]
hier1=hclust(as.dist(diss1), method="average" )
pdf(file = paste(Sys.Date(),name, "_WGCNA_dentro2.pdf", sep = ""), width = 11, height = 10, useDingbats = F)
plotDendroAndColors(hier1, 
                    dynamicColors[restGenes], 
                    "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, 
                    guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")
dev.off()

#### Correlate genes with module eigengenes
MEs <- moduleEigengenes(datExpr[,restGenes], dynamicColors[restGenes])
MEs <- orderMEs(MEs$eigengenes)

#### Merge similar modules
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")

pdf(file = paste(Sys.Date(),name, "_WGCNA_METree.pdf", sep = ""), width = 11, height = 10, useDingbats = F)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
MEDissThres = 0.15
abline(h=MEDissThres, col = "red", lty = 2)
dev.off()

merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
MEs <- merge$newMEs
dynamicColors <- merge$colors

#### Relating modules to external factors such as clusters
moduleTraitCor <- list()
moduleTraitPvalue <- list()
header = as.character(GB.keep$predicted)
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

ColSideColors = gg_color_hue(4)

pdf(file = paste(Sys.Date(),name, "_WGCNA_module_vs_predicted.pdf", sep = ""), width = 11, height = 10, useDingbats = F)
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
          #ColSideColors=as.matrix(ColSideColors),
          #ColSideColorsSize = 1,
          main = "Correlation between module and clusters")
dev.off()

#### Find top 5 genes with highest correlation to each module ####
moduleGeneCor <- apply(Y2[restGenes,], 1, function(x){cor(MEs, x, use = "p")})
moduleGenePvalue <- apply(moduleGeneCor, 2, 
                          function(x){corPvalueStudent(x, length(which(restGenes)))})
rownames(moduleGeneCor) <- rownames(moduleGenePvalue) <- colnames(MEs)
module_colors= setdiff(unique(dynamicColors), "grey")
moduleHighCorGenes <- list()
for (i in 1:length(module_colors)){
  moduleHighCorGenes[[i]] <- moduleGeneCor[i, dynamicColors[restGenes] == module_colors[i]]
  moduleHighCorGenes[[i]] <- moduleHighCorGenes[[i]][order(moduleHighCorGenes[[i]], decreasing = T)]
}
moduleHighCorGenes <- apply(moduleGeneCor, 1, function(x){
  list(x[order(x, decreasing = T)], colnames(moduleGeneCor)[order(x, decreasing = T)])})
moduleHighCorGenes <- lapply(moduleHighCorGenes, function(x){names(x[[1]][x[[1]] > 0.6])})
moduleHighCorGenes <- unique(moduleHighCorGenes)

Colors <- substr(rownames(moduleGeneCor), 3, 100)
RowSideColors <- c()
for(i in 1:length(moduleHighCorGenes)){
  RowSideColors <- c(RowSideColors, rep(module_colors[i], times = length(moduleHighCorGenes[[i]])))
}
ColSideColors <- gg_color_hue(5)[as.numeric(as.numeric(t_sne2$IDX)+2)]
ColSideColors <- matrix(ColSideColors, ncol = 1)
colsep <- t_sne2$IDX
colsep <- lapply(unique(colsep[order(colsep)]), function(x){length(which(colsep == x))})
colsep <- cumsum(unlist(colsep))
rowsep <- cumsum(unlist(lapply(moduleHighCorGenes, length)))
rec <- Y2[unlist(moduleHighCorGenes), order(as.numeric(t_sne2$IDX))]
#rec <- Y2[names(unlist(moduleHighCorGenes)), order(t_sne$kmeans)]
pairs.breaks <- seq(-3, 3, length.out=101);
heatmap.3(rec,
          key=F, keysize = 0.8,
          col = greenred(100),
          breaks=pairs.breaks,
          density.info = "histogram", 
          symkey=T, 
          main="Marker gene expression", 
          trace="none", cexRow=1, cexCol = 0.6, 
          Rowv = F, 
          Colv = F, #as.dendrogram(hca), 
          dendrogram = "none",
          ColSideColors = as.matrix(ColSideColors[order(t_sne2$IDX),]),
          #ColSideColors = ColSideColors[order(t_sne$kmeans),],
          ColSideColorsSize = 1,
          RowSideColors = t(as.matrix(RowSideColors)),
          RowSideColorsSize = 1,
          scale = "row",
          colsep = colsep,
          rowsep = rowsep,
          labRow = substr(unlist(moduleHighCorGenes), start = 20, stop = 100),
          labCol = "",
          na.rm = F
);

#### Plot modular genes
module_colors= setdiff(unique(dynamicColors), "grey")
m <- list()
for (i in 1:length(module_colors)){
  m[[i]] <- rownames(Y2)[which(dynamicColors==module_colors[i])]
}
rec <- Y2[unlist(m), order(as.numeric(t_sne2$IDX))]
#rec <- Y2[names(unlist(moduleHighCorGenes)), order(t_sne$kmeans)]
pairs.breaks <- seq(-3, 3, length.out=101);
heatmap.3(rec,
          key=F, keysize = 0.8,
          col = greenred(100),
          breaks=pairs.breaks,
          density.info = "histogram", 
          symkey=T, 
          main="Marker gene expression", 
          trace="none", cexRow=1, cexCol = 0.6, 
          Rowv = F, 
          Colv = F, #as.dendrogram(hca), 
          dendrogram = "none",
          ColSideColors = as.matrix(ColSideColors[order(t_sne2$IDX),]),
          ColSideColorsSize = 1,
          #RowSideColors = t(as.matrix(RowSideColors)),
          #RowSideColorsSize = 1,
          scale = "row",
          colsep = colsep,
          #rowsep = rowsep,
          labRow = substr(unlist(m), start = 20, stop = 100),
          labCol = "",
          na.rm = F
);

#Extract modules
load("2016-02-19_iPSC_TOM_sfp_5.Rdata")
module_colors= setdiff(unique(dynamicColors), "grey")
for (c in module_colors){
  module <- rownames(Y2)[which(dynamicColors==c)]
  write.csv(module, file = paste0(Sys.Date(), "_module_", c, ".csv", sep=""))
}

#### Plot gene network 
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
    colnames(modTOM) <- rownames(modTOM) <- substr(rownames(modTOM), 17, 100)
    
    threshold = 0.0
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
        width = 10, height = 9, useDingbats = F)
    p <- ggplot(data = fr, aes(x = x, y = y))+
      geom_segment(data = g, aes(x = from.x, xend = to.x, y = from.y, yend = to.y), colour = "grey90") +
      geom_point(data = fr, aes(x = x, y = y, size = size), colour = "green", alpha = 0.75) +
      geom_text(data = fr, aes(x = x, y = y, label = network_name), vjust = -0.5, size = 2) +
      theme_bw() +
      theme(rect = element_blank(),
            line = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = "none") +
      ggtitle(paste("Module", module, "threshold =", threshold, sep = " ")) +
      theme(plot.title = element_text(size=20, face="bold", vjust=2))
    print(p)
    dev.off(which = dev.cur())
  }, error=function(e){cat("Module: ", module, " failed\n")})
  }
}

