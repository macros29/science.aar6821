library("reshape")
library("tsne")
source("./SNN-Cliq-Gap.R")
setwd("./")
load("./t_sne.Rdata")

path.to.Cliq.py <- "./"
B <- 1000
r <- seq(0.3, 0.8, by = 0.1)
m <- seq(0.3, 0.8, by = 0.1)
name <- "edge_temp.txt" #paste0(Sys.Date(), "_no_correction_edge_file_k_", k, ".txt", sep = "")

x <- t_sne
cat("Start SNN-cliq\n\n")
k = 4
SNN(data = x, outfile = name, k = k)

n_clust <- matrix(nrow = length(r), ncol = length(m), dimnames = list(r, m))
W.k <- matrix(nrow = length(r), ncol = length(m), dimnames = list(r, m))
E.logW <- matrix(nrow = length(r), ncol = length(m), dimnames = list(r, m))
SE.sim <- matrix(nrow = length(r), ncol = length(m), dimnames = list(r, m))

n <- nrow(x)
for(par.r in 1:length(r)){
  for(par.m in 1:length(m)){
    cat(paste0("r = ", r[par.r], ", m = ", m[par.m], "\n\n", sep = ""))
    w <- getW(dataset = x, path.to.Cliq.py = path.to.Cliq.py, r = r[par.r], m = m[par.m], name = name)
    n_clust[par.r, par.m] <- getClusterNumber(dataset = x, path.to.Cliq.py = path.to.Cliq.py, r = r[par.r], m = m[par.m], name = name)
    xs <- scale(x, center = TRUE, scale = FALSE)
    m.x <- rep(attr(xs, "scaled:center"), each = n)
    V.sx <- svd(xs)$v
    rng.x1 <- apply(xs %*% V.sx, 2, range)
    logWks <- c()
    #if (verbose) 
    #  cat("Bootstrapping, b = 1,2,..., B (= ", B, ")  [one \".\" per sample]:\n", 
    #      sep = "")
    for (b in 1:B) {
      z1 <- apply(rng.x1, 2, function(M, nn) runif(nn, min = M[1], 
                                                   max = M[2]), nn = n)
      z <- tcrossprod(z1, V.sx) + m.x
      logWks <- c(logWks, log(getW(dataset = z, path.to.Cliq.py, r[par.r], m[par.m], name)))
      cat("=")}
    #if (verbose) 
    #  cat(".", if (b%%50 == 0) 
    #    paste(b, "\n"))
    W.k[par.r, par.m] <- w
    E.logW[par.r, par.m] <- mean(logWks)
    SE.sim[par.r, par.m] <- sqrt((1 + 1/B) * var(logWks))
    cat("\n\n")
  }
  #if (verbose && (B%%50 != 0)) 
  #cat("", B, "\n")
}

logWks <- log(W.k)
gap = E.logW - logWks
gap <- melt(gap)
colnames(gap) <- c("r", "m", "gap")
gap$logWks <- as.vector(logWks)
gap$E.logW <- as.vector(E.logW)
gap$SE.sim <- as.vector(SE.sim)
gap$n_clust <- as.vector(n_clust)
gap <- gap[order(gap$gap), ]
gap_elbow <- vector()
for(i in 1:nrow(gap)){
  gap_elbow[i] <- (gap$gap[i] - gap$gap[i+1] + gap$SE.sim[i+1])
}
gap$gap_elbow <- gap_elbow

save(gap, file = "gap.Rdata")

cat("gap statistics saved\n")