setwd("~/Data/Da_Mi/")
load("working.RData")
pcs <- list(c(1:3), c(1:9), c(1:7), c(1:3), c(1:5), c(1:6))
for(i in c(1:6) ){
  df <- dfs[[i]]
  
  m <- t(df@dr$pca@cell.embeddings[,pcs[[i]]]) # E12.5 dMGE = 1:3; vMGE = 1:9; CGE = 1:7; E14.5 dMGE = 1:3; vMGE = 1:5; CGE = 1:6
  hc <- hclust(as.dist(1-cor(m, method = "pearson", use = "pairwise.complete.obs")), method = "ward.D2")
  cl <- cutree(hc, k=2)
  df@ident <- as.factor(cl)
  
  df.markers <- FindAllMarkers(genes.use = df@var.genes, object = df, only.pos = TRUE, 
                              min.pct = 0.1, thresh.use = 0.25, return.thresh = 0.01)
  # 
  # write.csv(df.markers, 
  #           file = paste(paste(samples[[i]][2],"_", samples[[i]][1],"/",Sys.Date(), sep = ""),
  #                        samples[[i]][2], samples[[i]][1],"markers.csv", sep = "_"))
  # 
  # top10 <- df.markers %>% group_by(cluster) %>% top_n(10, avg_diff)
  # 
  # pdf(file = paste(paste(samples[[i]][2],"_", samples[[i]][1],"/",Sys.Date(), sep = ""), 
  #                  samples[[i]][2], samples[[i]][1],"3_clusters_heatmap_top10.pdf", sep = "_"), 
  #     width = 20, height = 20, useDingbats = F);
  # DoHeatmap(object = df, 
  #           genes.use = top10$gene, 
  #           slim.col.label = TRUE,group.spacing = 0.3,
  #           remove.key = TRUE)
  # dev.off()
  # 
  glc <- read.csv("./raw_data/gene list for comparison.csv", header = F)
  glc$V1 <- toupper(glc$V1)
  glc$V1 <- sapply(glc$V1, function(x){paste("\\b", x, "$", sep = "")})
  glc$V1 <- sapply(glc$V1,function(x){grep(x, rownames(df@scale.data), value = T)})
  glc$V1 <- as.character(glc$V1)
  # glc[glc$V1 %in% df.markers$gene,]
  # temp <- cbind(df.markers[na.omit(match(glc$V1, df.markers$gene)),], glc[sort(na.omit(match(df.markers$gene, glc$V1))),])
  # write.csv(temp, 
  #           file = paste(paste(samples[[i]][2],"_", samples[[i]][1],"/",Sys.Date(), sep = ""),
  #                        samples[[i]][2], samples[[i]][1],"Da_Mi_list.csv", sep = "_"))
  dfs[[i]] <- df
  dfs.markers[[i]] <- df.markers
  # 
  t <- df@ident

  gene_mat0 <- df@raw.data[,df@ident == "1"]
  save(gene_mat0, file = paste(samples[[i]][2],"_", samples[[i]][1],"/",
                          "Round1.C1.AllGenes_NewInput.dat", sep = ""))
  gene_mat0 <- df@raw.data[,df@ident == "2"]
  save(gene_mat0, file = paste(samples[[i]][2],"_", samples[[i]][1],"/",
                                "Round1.C2.AllGenes_NewInput.dat", sep = ""))
}