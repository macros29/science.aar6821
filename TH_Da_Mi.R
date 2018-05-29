grep("\\bTH$", df@var.genes, value = T)
hist(as.numeric(df@raw.data["ENSMUSG00000000214|TH",]), breaks = 1000, ylim = c(0,6))
hist(as.numeric(df@raw.data["ENSMUSG00000000214|TH",]))

TH.cells <- as.numeric(df@raw.data["ENSMUSG00000000214|TH",]) > 0
keep <- df@meta.data$time_point == "E14.5"

TH.cells <- df@scale.data[,TH.cells & keep]
markers <- c("TH","GAD1", "GAD2", "SST", "PVALB", "NPY", "NOS1", "CALB1", "CALB2", "VIP", "RBFOX3", "ETV1",
             "ETV5", "DBH", "DDC", "SLC18A2","SLC6A3", "SP8", "NR2F2", "NKX2.1", "SOX6", "PAX6")
markers <- paste("\\b", markers, "$", sep = "")
markers <- unlist(sapply(markers,function(x){grep(x, rownames(df@scale.data), value = T)}))

temp <- TH.cells[ markers, ]

pairs.breaks <- seq(-2, 2, length.out=101);

pdf(file = paste("./", Sys.Date(), "_marker_heatmap_TH_cells_E14.5.pdf", sep = ""), 
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
          #ColSideColors = as.matrix(ColSideColors),
          #ColSideColorsSize = 2,
          #RowSideColors = RowSideColors,
          dendrogram = "both",
          scale = "row",
          #colsep = colsep,
          sepcolor = "black",
          labRow = substr(rownames(temp), 20, 100),
          labCol = "",
          na.rm = F);
dev.off(dev.cur());