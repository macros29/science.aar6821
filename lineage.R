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

df <- RunTSNE(df, dims.use = pcs[[i]], seed.use = 7)

data <- data.frame(df@dr$tsne@cell.embeddings)
# TSNEPlot(object = df, pt.size = 3)

lin <- read.csv("./raw_data/INs lineage marker 21-09-17.csv", header = F, stringsAsFactors = F)
lin[,1] <- sapply(lin[,1], function(x){paste("\\b", x, "$", sep = "")})
lin[,1] <- as.character(sapply(lin[,1],function(x){grep(x, rownames(df@scale.data), value = T)}))
lin <- lin[!(lin$V1 == "character(0)"),]

lin.gord <- read.csv("./raw_data/INs lineage marker Cord 28-09-17.csv", header = F, stringsAsFactors = F)
lin.gord[,1] <- sapply(lin.gord[,1], function(x){paste("\\b", x, "$", sep = "")})
lin.gord[,1] <- as.character(sapply(lin.gord[,1],function(x){grep(x, rownames(df@scale.data), value = T)}))
lin.gord <- lin.gord[!(lin.gord$V1 == "character(0)"),]

temp <- t(df@scale.data[unique(lin.gord$V1),])
d <- as.dist(1-cor(temp, method = "s"))
hcu <- hclust(d, method = "complete")
plot(hcu)

heatmap.3(d, Rowv = as.dendrogram(hcu),Colv = as.dendrogram(hcu))

ID2_lin <- c("\\bID2$", "\\bGPX3$", "\\bCPNE5$", "\\bCITED1$")
ID2_lin <- as.character(sapply(ID2_lin,function(x){grep(x, rownames(df@scale.data), value = T)}))

SST_lin <- c("CDK6", "ID1", "ID3", "NR2E1", "HTR3A", "SST", "RELN")
SST_lin <- sapply(SST_lin, function(x){paste("\\b", x, "$", sep = "")})
SST_lin <- as.character(sapply(SST_lin,function(x){grep(x, rownames(df@scale.data), value = T)}))

PV_lin <- c("NETO2", "TPBG", "RXRA", "MXD1", "NETO1", "CPLX1", "TRPS1", "IKZF2")
SST_lin <- sapply(SST_lin, function(x){paste("\\b", x, "$", sep = "")})
SST_lin <- as.character(sapply(SST_lin,function(x){grep(x, rownames(df@scale.data), value = T)}))

temp <- df@scale.data[ID2_lin,]
pca <- prcomp(t(temp))

ggplot(data, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(aes(color = -1*pca$x[,1], alpha = -1*pca$x[,1]), size = 3) +
  scale_color_continuous(high = "purple", low = "grey")
