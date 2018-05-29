keep <- GB$region %in% c("vMGE", "dMGE", "CGE") & GB$time_point %in% "E14.5" & !(GB$IFC %in% c("189","180","181","182","186","190","191","192","193","194","195", "196"))
df <- CreateSeuratObject(raw.data = as.vector(expr[,keep]), min.cells = 5, min.genes = 1000, project = "Da_Mi")
df <- NormalizeData(object = df, normalization.method = "LogNormalize", scale.factor = 10000)
df <- ScaleData(df)

mito.genes <- grep(pattern = "\\bMT-", x = rownames(x = df@data), value = TRUE)
percent.mito <- colSums(df@raw.data[mito.genes, ])/colSums(df@raw.data)

# AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
df <- AddMetaData(object = df, metadata = percent.mito, col.name = "percent.mito")
df@meta.data <- cbind(df@meta.data, GB[df@cell.names,])

df <- FindVariableGenes(object = df, mean.function = ExpMean, dispersion.function = LogVMR, 
                        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = F)

temp <- as.character(df@ident)
for(i in 4:6){
  temp[match(dfs[[i]]@cell.names, df@cell.names)] <- paste(names(dfs[i]), dfs[[i]]@ident, sep = "_")
}  
df@ident <- as.factor(temp)
names(df@ident) <- df@cell.names

glc <- read.csv("./raw_data/gene list for comparison.csv", header = F)
glc$V1 <- toupper(glc$V1)
glc$V1 <- sapply(glc$V1, function(x){paste("\\b", x, "$", sep = "")})
glc$V1 <- sapply(glc$V1,function(x){grep(x, rownames(df@scale.data), value = T)})
glc$V1 <- as.character(glc$V1)
glc <- glc[glc$V1 != "character(0)",]

DnvVn.glc <- FindMarkers(df,genes.use = glc$V1, ident.1 = "E14.5_dMGE_2", ident.2 = "E14.5_vMGE_1", only.pos = F)
DnvCn.glc <- FindMarkers(df,genes.use = glc$V1, ident.1 = "E14.5_dMGE_2", ident.2 = "E14.5_CGE_1", only.pos = F)
VnvCn.glc <- FindMarkers(df,genes.use = glc$V1, ident.1 = "E14.5_vMGE_1", ident.2 = "E14.5_CGE_1", only.pos = F)
DpvVp.glc <- FindMarkers(df,genes.use = glc$V1, ident.1 = "E14.5_dMGE_1", ident.2 = "E14.5_vMGE_2", only.pos = F)
DpvCp.glc <- FindMarkers(df,genes.use = glc$V1, ident.1 = "E14.5_dMGE_1", ident.2 = "E14.5_CGE_2", only.pos = F)
VpvCp.glc <- FindMarkers(df,genes.use = glc$V1, ident.1 = "E14.5_vMGE_2", ident.2 = "E14.5_CGE_2", only.pos = F)

write.csv(DnvVn.glc, file = "2017-09-20_E14.5_Dn_Vn_glc.csv")
write.csv(DnvCn.glc, file = "2017-09-20_E14.5_Dn_Cn_glc.csv")
write.csv(VnvCn.glc, file = "2017-09-20_E14.5_Vn_Cn_glc.csv")
write.csv(DpvVp.glc, file = "2017-09-20_E14.5_Dp_Vp_glc.csv")
write.csv(DpvCp.glc, file = "2017-09-20_E14.5_Dp_Cp_glc.csv")
write.csv(VpvCp.glc, file = "2017-09-20_E14.5_Vp_Cp_glc.csv")

pdf(file = paste(Sys.Date(),"Prog_E14.5_vMGE_vs_CGE_glc.pdf", sep = "_"), 
    width = 18, height = 30, useDingbats = F);
VlnPlot(object = df, features.plot = rownames(VpvCp.glc[VpvCp.glc$p_val < 0.05,]), 
        ident.include = c("E14.5_vMGE_2", "E14.5_CGE_2"),
        nCol = 3)
dev.off()

pdf(file = paste(Sys.Date(),"Neuron_E14.5_vMGE_vs_CGE_glc.pdf", sep = "_"), 
    width = 18, height = 30, useDingbats = F);
VlnPlot(object = df, features.plot = rownames(VnvCn.glc[VnvCn.glc$p_val < 0.05,]), 
        ident.include = c("E14.5_vMGE_1", "E14.5_CGE_1"),
        nCol = 3)
dev.off()

#### Lineage (Da Mi) ####
ilm <- read.csv("./raw_data/INs lineage marker 28-08-17.csv", header = F)
ilm$V1 <- toupper(ilm$V1)
ilm$V1 <- sapply(ilm$V1, function(x){paste("\\b", x, "$", sep = "")})
ilm$V1 <- sapply(ilm$V1,function(x){grep(x, rownames(df@scale.data), value = T)})
ilm$V1 <- as.character(ilm$V1)
ilm <- ilm[ilm$V1 != "character(0)",]

DnvVn.ilm <- FindMarkers(df,genes.use = ilm$V1, ident.1 = "E12.5_dMGE_1", ident.2 = "E12.5_vMGE_1", only.pos = F)
DnvCn.ilm <- FindMarkers(df,genes.use = ilm$V1, ident.1 = "E12.5_dMGE_1", ident.2 = "E12.5_CGE_2", only.pos = F)
VnvCn.ilm <- FindMarkers(df,genes.use = ilm$V1, ident.1 = "E12.5_vMGE_1", ident.2 = "E12.5_CGE_2", only.pos = F)
DpvVp.ilm <- FindMarkers(df,genes.use = ilm$V1, ident.1 = "E12.5_dMGE_2", ident.2 = "E12.5_vMGE_2", only.pos = F)
DpvCp.ilm <- FindMarkers(df,genes.use = ilm$V1, ident.1 = "E12.5_dMGE_2", ident.2 = "E12.5_CGE_1", only.pos = F)
VpvCp.ilm <- FindMarkers(df,genes.use = ilm$V1, ident.1 = "E12.5_vMGE_2", ident.2 = "E12.5_CGE_1", only.pos = F)

write.csv(DnvVn.ilm, file = "2017-09-20_E12.5_Dn_Vn_ilm.csv")
write.csv(DnvCn.ilm, file = "2017-09-20_E12.5_Dn_Cn_ilm.csv")
write.csv(VnvCn.ilm, file = "2017-09-20_E12.5_Vn_Cn_ilm.csv")
write.csv(DpvVp.ilm, file = "2017-09-20_E12.5_Dp_Vp_ilm.csv")
write.csv(DpvCp.ilm, file = "2017-09-20_E12.5_Dp_Cp_ilm.csv")
write.csv(VpvCp.ilm, file = "2017-09-20_E12.5_Vp_Cp_ilm.csv")

pdf(file = paste(Sys.Date(),"E12.5_DnvVn_ilm.pdf", sep = "_"), 
    width = 18, height = 30, useDingbats = F);
VlnPlot(object = df, features.plot = rownames(DnvVn.ilm[DnvVn.ilm$p_val < 0.05,]), 
        ident.include = c("E12.5_dMGE_1", "E12.5_vMGE_1"),
        nCol = 3)
dev.off()

pdf(file = paste(Sys.Date(),"Prog_E12.5_VpvCp_ilm.pdf", sep = "_"), 
    width = 18, height = 30, useDingbats = F);
VlnPlot(object = df, features.plot = rownames(VpvCp.ilm[VpvCp.ilm$p_val < 0.05,]), 
        ident.include = c("E12.5_vMGE_2", "E12.5_CGE_1"),
        nCol = 3)
dev.off()


DnvVn.ilm <- FindMarkers(df,genes.use = ilm$V1, ident.1 = "E14.5_dMGE_2", ident.2 = "E14.5_vMGE_1", only.pos = F)
DnvCn.ilm <- FindMarkers(df,genes.use = ilm$V1, ident.1 = "E14.5_dMGE_2", ident.2 = "E14.5_CGE_1", only.pos = F)
VnvCn.ilm <- FindMarkers(df,genes.use = ilm$V1, ident.1 = "E14.5_vMGE_1", ident.2 = "E14.5_CGE_1", only.pos = F)
DpvVp.ilm <- FindMarkers(df,genes.use = ilm$V1, ident.1 = "E14.5_dMGE_1", ident.2 = "E14.5_vMGE_2", only.pos = F)
DpvCp.ilm <- FindMarkers(df,genes.use = ilm$V1, ident.1 = "E14.5_dMGE_1", ident.2 = "E14.5_CGE_2", only.pos = F)
VpvCp.ilm <- FindMarkers(df,genes.use = ilm$V1, ident.1 = "E14.5_vMGE_2", ident.2 = "E14.5_CGE_2", only.pos = F)

write.csv(DnvVn.ilm, file = "2017-09-20_E14.5_Dn_Vn_ilm.csv")
write.csv(DnvCn.ilm, file = "2017-09-20_E14.5_Dn_Cn_ilm.csv")
write.csv(VnvCn.ilm, file = "2017-09-20_E14.5_Vn_Cn_ilm.csv")
write.csv(DpvVp.ilm, file = "2017-09-20_E14.5_Dp_Vp_ilm.csv")
write.csv(DpvCp.ilm, file = "2017-09-20_E14.5_Dp_Cp_ilm.csv")
write.csv(VpvCp.ilm, file = "2017-09-20_E14.5_Vp_Cp_ilm.csv")

pdf(file = paste(Sys.Date(),"E14.5_VnvCn_ilm.pdf", sep = "_"), 
    width = 18, height = 30, useDingbats = F);
VlnPlot(object = df, features.plot = rownames(VnvCn.ilm[VnvCn.ilm$p_val < 0.05,]), 
        ident.include = c("E14.5_vMGE_1", "E14.5_CGE_1"),
        nCol = 3)
dev.off()

pdf(file = paste(Sys.Date(),"Prog_E14.5_dMGE_vs_CGE_ilm.pdf", sep = "_"), 
    width = 18, height = 30, useDingbats = F);
VlnPlot(object = df, features.plot = rownames(DpvCp.ilm[DpvCp.ilm$p_val < 0.05,]), 
        ident.include = c("E14.5_dMGE_1", "E14.5_CGE_2"),
        nCol = 3)
dev.off()

#### Lineage (Gord) ####
ilmg <- read.csv("./raw_data/INs lineage marker Cord 18-09-17.csv", header = F)
ilmg$V1 <- toupper(ilmg$V1)
ilmg$V1 <- sapply(ilmg$V1, function(x){paste("\\b", x, "$", sep = "")})
ilmg$V1 <- sapply(ilmg$V1,function(x){grep(x, rownames(df@scale.data), value = T)})
ilmg$V1 <- as.character(ilmg$V1)
ilmg <- ilmg[ilmg$V1 != "character(0)",]

DnvVn.ilmg <- FindMarkers(df,genes.use = ilmg, ident.1 = "E12.5_dMGE_1", ident.2 = "E12.5_vMGE_1", only.pos = F)
DnvCn.ilmg <- FindMarkers(df,genes.use = ilmg, ident.1 = "E12.5_dMGE_1", ident.2 = "E12.5_CGE_2", only.pos = F)
VnvCn.ilmg <- FindMarkers(df,genes.use = ilmg, ident.1 = "E12.5_vMGE_1", ident.2 = "E12.5_CGE_2", only.pos = F)
DpvVp.ilmg <- FindMarkers(df,genes.use = ilmg, ident.1 = "E12.5_dMGE_2", ident.2 = "E12.5_vMGE_2", only.pos = F)
DpvCp.ilmg <- FindMarkers(df,genes.use = ilmg, ident.1 = "E12.5_dMGE_2", ident.2 = "E12.5_CGE_1", only.pos = F)
VpvCp.ilmg <- FindMarkers(df,genes.use = ilmg, ident.1 = "E12.5_vMGE_2", ident.2 = "E12.5_CGE_1", only.pos = F)

write.csv(DnvVn.ilmg, file = "2017-09-20_E12.5_Dn_Vn_ilmg.csv")
write.csv(DnvCn.ilmg, file = "2017-09-20_E12.5_Dn_Cn_ilmg.csv")
write.csv(VnvCn.ilmg, file = "2017-09-20_E12.5_Vn_Cn_ilmg.csv")
write.csv(DpvVp.ilmg, file = "2017-09-20_E12.5_Dp_Vp_ilmg.csv")
write.csv(DpvCp.ilmg, file = "2017-09-20_E12.5_Dp_Cp_ilmg.csv")
write.csv(VpvCp.ilmg, file = "2017-09-20_E12.5_Vp_Cp_ilmg.csv")

pdf(file = paste(Sys.Date(),"E12.5_VnvCn_ilmg.pdf", sep = "_"), 
    width = 18, height = 30, useDingbats = F);
VlnPlot(object = df, features.plot = rownames(VnvCn.ilmg[VnvCn.ilmg$p_val < 0.05,]), 
        ident.include = c("E12.5_vMGE_1", "E12.5_CGE_2"),
        nCol = 3)
dev.off()

pdf(file = paste(Sys.Date(),"Prog_E12.5_VpvCp_ilmg.pdf", sep = "_"), 
    width = 18, height = 30, useDingbats = F);
VlnPlot(object = df, features.plot = rownames(VpvCp.ilmg[VpvCp.ilmg$p_val < 0.05,]), 
        ident.include = c("E12.5_vMGE_2", "E12.5_CGE_1"),
        nCol = 3)
dev.off()

DnvVn.ilmg <- FindMarkers(df,genes.use = ilmg, ident.1 = "E14.5_dMGE_2", ident.2 = "E14.5_vMGE_1", only.pos = F)
DnvCn.ilmg <- FindMarkers(df,genes.use = ilmg, ident.1 = "E14.5_dMGE_2", ident.2 = "E14.5_CGE_1", only.pos = F)
VnvCn.ilmg <- FindMarkers(df,genes.use = ilmg, ident.1 = "E14.5_vMGE_1", ident.2 = "E14.5_CGE_1", only.pos = F)
DpvVp.ilmg <- FindMarkers(df,genes.use = ilmg, ident.1 = "E14.5_dMGE_1", ident.2 = "E14.5_vMGE_2", only.pos = F)
DpvCp.ilmg <- FindMarkers(df,genes.use = ilmg, ident.1 = "E14.5_dMGE_1", ident.2 = "E14.5_CGE_2", only.pos = F)
VpvCp.ilmg <- FindMarkers(df,genes.use = ilmg, ident.1 = "E14.5_vMGE_2", ident.2 = "E14.5_CGE_2", only.pos = F)

write.csv(DnvVn.ilmg, file = "2017-09-20_E14.5_Dn_Vn_ilmg.csv")
write.csv(DnvCn.ilmg, file = "2017-09-20_E14.5_Dn_Cn_ilmg.csv")
write.csv(VnvCn.ilmg, file = "2017-09-20_E14.5_Vn_Cn_ilmg.csv")
write.csv(DpvVp.ilmg, file = "2017-09-20_E14.5_Dp_Vp_ilmg.csv")
write.csv(DpvCp.ilmg, file = "2017-09-20_E14.5_Dp_Cp_ilmg.csv")
write.csv(VpvCp.ilmg, file = "2017-09-20_E14.5_Vp_Cp_ilmg.csv")

pdf(file = paste(Sys.Date(),"E14.5_VnvCn_ilmg.pdf", sep = "_"), 
    width = 18, height = 30, useDingbats = F);
VlnPlot(object = df, features.plot = rownames(VnvCn.ilmg[VnvCn.ilmg$p_val < 0.05,]), 
        ident.include = c("E14.5_vMGE_1", "E14.5_CGE_1"),
        nCol = 3)
dev.off()

pdf(file = paste(Sys.Date(),"Prog_E14.5_VpvCp_ilmg.pdf", sep = "_"), 
    width = 18, height = 30, useDingbats = F);
VlnPlot(object = df, features.plot = rownames(VpvCp.ilmg[VpvCp.ilmg$p_val < 0.05,]), 
        ident.include = c("E14.5_vMGE_2", "E14.5_CGE_2"),
        nCol = 3)
dev.off()





df.markers <- FindAllMarkers(genes.use = df@var.genes, object = df, only.pos = TRUE, 
                             min.pct = 0.1, thresh.use = 0.25, return.thresh = 0.01)

grep("\\bTH$", rownames(df@scale.data), value = T)
plot(df@scale.data["ENSMUSG00000000214|TH",])
top10 <- df.markers %>% group_by(cluster) %>% top_n(10, avg_diff)

pdf(file = paste(Sys.Date(),"E12.5_dMGE_vs_vMGE_top10.pdf", sep = "_"), 
    width = 10, height = 5, useDingbats = F);
DoHeatmap(object = df, 
          genes.use = top10$gene, 
          slim.col.label = TRUE,group.spacing = 0.3,
          remove.key = TRUE)
dev.off()

DoHeatmap(object = dfs[[6]], 
          genes.use = top10$gene, 
          slim.col.label = TRUE,group.spacing = 0.3,
          remove.key = TRUE)
