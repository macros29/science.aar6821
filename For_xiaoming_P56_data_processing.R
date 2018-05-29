RP.genes <- grep("Rp", rownames(counts.allen))

df0.allen <- CreateSeuratObject(raw.data = as.vector(counts.allen[-c(RP.genes),keep]), 
                                min.cells = 5, min.genes = 1000, project = "Da_Mi")
df0.allen@meta.data <- cbind(df0.allen@meta.data, GB.allen[keep,])
df0.allen <- NormalizeData(object = df0.allen, normalization.method = "LogNormalize", 
                           scale.factor = 1000000)

df0.allen <- FindVariableGenes(object = df0.allen, mean.function = ExpMean, 
                               dispersion.function = LogVMR, 
                               x.low.cutoff = 0.5, 
                               x.high.cutoff = Inf, 
                               y.cutoff = 0.5, 
                               y.high.cutoff = 5, 
                               do.plot = T)