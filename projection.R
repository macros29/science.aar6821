counts.allen <- read.csv("./raw_data/GSE71585_RefSeq_counts.csv", header = T, stringsAsFactors = F)
rownames(counts.allen) <- toupper(counts.allen[,1])
counts.allen <- counts.allen[,-1]

clu.allen <- read.csv("./raw_data/GSE71585_Clustering_Results.csv", header = T, stringsAsFactors = F)

keep <- clu.allen$broad_type == "GABA-ergic Neuron"

df.allen <- CreateSeuratObject(raw.data = as.vector(counts.allen[,keep]), min.cells = 5, min.genes = 1000, project = "Da_Mi")
df.allen <- NormalizeData(object = df.allen, normalization.method = "LogNormalize", scale.factor = 10000)

percent.mito <- colSums(df.allen@raw.data[mito.genes, ])/colSums(df.allen@raw.data)

# AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
keep <- GB$time_point %in% samples[[i]] & !(GB$IFC %in% c("189","180","181","182","186","190","191","192","193","194","195", "196")) 
temp <- expr[,keep]

length(which(rownames(df.allen@scale.data) %in% substr(rownames(temp), 20, 100)))

df <- CreateSeuratObject(raw.data = as.vector(expr[,keep]), min.cells = 5, min.genes = 1000, project = "Da_Mi")
df <- NormalizeData(object = df, normalization.method = "LogNormalize", scale.factor = 10000)
mito.genes <- grep(pattern = "\\bMT-", x = rownames(x = df@data), value = TRUE)
percent.mito <- colSums(df@raw.data[mito.genes, ])/colSums(df@raw.data)
df <- AddMetaData(object = df, metadata = percent.mito, col.name = "percent.mito")
df@meta.data <- cbind(df@meta.data, GB[df@cell.names,])

df <- ScaleData(df,vars.to.regress = c("nGene", "mouse"))

df.allen <- AddMetaData(object = df.allen, metadata = percent.mito, col.name = "percent.mito")
df.allen@meta.data <- cbind(df.allen@meta.data, clu.allen[keep,])
#df.allen <- ScaleData(df.allen)
df.allen <- ScaleData(df.allen)
df.allen <- FindVariableGenes(object = df.allen, do.plot = F)

hvg.E12.5 <-rownames(x = head(x = df@hvg.info, n = 2000))
hvg.E12.5 <- substr(hvg.E12.5, 20, 100)
hvg.allen <- rownames(x = head(x = df.allen@hvg.info, n = 2000))
hvg.union <- union(hvg.allen, hvg.E12.5)

df@meta.data$source <- "Da_Mi"
df.allen@meta.data$source <- "allen"

rownames(df@raw.data) <- substr(rownames(df@raw.data), 20, 100)
rownames(df@scale.data) <- substr(rownames(df@scale.data), 20, 100)
df@var.genes <- substr(df@var.genes, 20, 100)

pbmc <- RunCCA(object = df, object2 = df.allen, genes.use = hvg.union)
DimHeatmap(object = pbmc, reduction.type = "cca", dim.use = 1:9, 
           do.balanced = T)

pbmc <- CalcVarExpRatio(object = pbmc, reduction.type = "pca", grouping.var = "source",
                        dims.use = 1:5)
pbmc <- AlignSubspace(object = pbmc, reduction.type = "cca", grouping.var = "source", 
                      dims.align = 1:13)

df <- RunPCA(object = df, pc.genes=df@var.genes, pcs.compute = 30, pcs.print = 1:12,genes.print = 10)
