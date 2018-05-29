keep <- GB$region %in% samples[[i]] & GB$time_point %in% samples[[i]] & !(GB$IFC %in% c("180","181","191","192","193","194","195", "196"))
df <- CreateSeuratObject(raw.data = as.vector(expr[,keep]), min.cells = 5, min.genes = 1000, project = "Da_Mi")
df <- NormalizeData(object = df, normalization.method = "LogNormalize", scale.factor = 10000)
mito.genes <- grep(pattern = "\\bMT-", x = rownames(x = df@data), value = TRUE)
percent.mito <- colSums(df@raw.data[mito.genes, ])/colSums(df@raw.data)
df <- AddMetaData(object = df, metadata = percent.mito, col.name = "percent.mito")
df@meta.data <- cbind(df@meta.data, GB[df@cell.names,])
df <- ScaleData(df,vars.to.regress = c("nGene", "mouse"))

df <- FindVariableGenes(object = df, mean.function = ExpMean, dispersion.function = LogVMR, 
                        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = F)

df <- RunPCA(object = df, pc.genes=df@var.genes, pcs.compute = 30, pcs.print = 1:12,genes.print = 10)

pdf(file = paste( Sys.Date(),samples[[i]][2], samples[[i]][1], "PCA_norm.pdf", sep = "_"), 
    width = 11, height = 10, useDingbats = F);
PCAPlot(df, group.by = "mouse")
dev.off()

pcs <- list(c(1:3), c(1:7), c(1:4), c(1:4), c(1:4), c(1:4))
df <- RunTSNE(df, dims.use = pcs[[i]], seed.use = 7)

df <- FindClusters(object = df, 
                   reduction.type = "tsne", 
                   dims.use = 1:2, 
                   resolution = 1, # dMGE, CGE = 0.8, E14.5 vMGE = 2
                   k.scale = 25, # dMGE, CGE = 50
                   prune.SNN = 0, #dMGE, CGE = 0.5
                   plot.SNN = F, 
                   print.output = F, 
                   save.SNN = F,
                   algorithm = 2,
                   force.recalc = TRUE, 
                   random.seed = 1,
                   k.param = 30) # dMGE, CGE = 30

#df <- RunDiffusion(df, dims.use = 1:3)

pdf(file = paste(Sys.Date(), samples[[i]][2], samples[[i]][1],  
                 "tSNE.pdf", sep = "_"), width = 11, 
    height = 10, useDingbats = F);
TSNEPlot(object = df, pt.size = 3)
dev.off()

# VvD <- FindMarkers(df, ident.1 = "vMGE", ident.2 = "dMGE", only.pos = T)
df.markers <- FindAllMarkers(genes.use = df@var.genes, object = df, only.pos = TRUE, 
                             min.pct = 0.1, thresh.use = 0.25, return.thresh = 0.01)
write.csv(df.markers, file = paste(samples[[i]][2], samples[[i]][1], "_markers_clusters.csv"))

lin <- read.csv("./raw_data/INs lineage marker 21-09-17.csv", header = F, stringsAsFactors = F)
lin[,1] <- sapply(lin[,1], function(x){paste("\\b", x, "$", sep = "")})
lin[,1] <- as.character(sapply(lin[,1],function(x){grep(x, rownames(df@scale.data), value = T)}))
lin <- lin[!(lin$V1 == "character(0)"),]

