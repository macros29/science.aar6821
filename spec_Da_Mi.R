#### use Spec analysis to  ####
## Load data
args <- commandArgs(trailingOnly =T)
print(paste("Start analyze", sep = " "))

library(parallel)
library(Seurat)
library(matrixStats)
source("~/Scripts/R/spec.R")
load("./E12.5_df.Rdata")

Data <- df@scale.data[df@var.genes,]
cell_type = df@ident

n <- length(unique(df@ident))
no_cores <- detectCores() - 1
cat("\n", "Number of Cores = ", no_cores, "\n")
cl <- makeCluster(no_cores, outfile='LOG.TXT')
clusterExport(cl, c("cell_type", "Data"))


spec_scores <- parApply(cl = cl, Data, 1, function(x){
  source("~/Scripts/R/spec.R");
  gene <- as.numeric(x);
  opt_bin <- optimizebinsize(gene, header = cell_type);
  cat("=")
  return(unlist(specy(gene, header = cell_type, binsize = opt_bin)));
})
stopCluster(cl)

colnames(spec_scores) <- df@var.genes
rownames(spec_scores)[1:n] <- unique(sort(df@ident));

res <- spec_scores[1:n,]
write.csv(res, file = paste(Sys.Date(), "E12.5_spec_scores_all_clusters.csv", sep = "_"))

#save(spec_scores, file = "spec_scores.Rdata")

cat("\n", "DONE")
