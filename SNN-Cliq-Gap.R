SNN<-function(data, outfile, k, distance){
  if(missing(data)){
    stop(paste("Input data missing.",help,sep="\n"))
  }
  if(missing(outfile)){
    stop(paste("Output file name missing.",help,sep="\n"))
  }
  if(missing(k)){
    k=3
  }
  if(missing(distance)){
    distance<-"euclidean"  # other distance options refer to dist() in R
  }
  #m<-as.data.frame(data)
  numSpl<-dim(data)[1]
  m<-dist(data, distance, diag=TRUE, upper=TRUE)
  #m<-as.dist(1-cor(data, method = "p")) # calculate correlation and use as distance
  #m <- data # use distance matrix directly
  x<-as.matrix(m)
  IDX<-t(apply(x,1,order)[1:k,]) # knn list
  
  edges<-list()              # SNN graph
  for (i in 1:numSpl){
    j<-i
    while (j<numSpl){
      j<-j+1
      shared<-intersect(IDX[i,], IDX[j,])
      if(length(shared)>0){			
        s<-k-0.5*(match(shared, IDX[i,])+match(shared, IDX[j,]))
        strength<-max(s)
        if (strength>0)
          edges<-rbind(edges, c(i,j,strength))
      }				
    }
  }
  write.table(edges, outfile, quote=FALSE, sep='\t',col.names=FALSE,row.names=FALSE)
}

SNNclust <- function (path.to.cliq.py, name, k, r, m, ...) {
  system(paste0("python ", path.to.Cliq.py,"Cliq.py -i ", name, "  -o IDX_k_", k, "_temp.txt -r ",
                r, " -m ", m), ignore.stdout = T)
  clusts <- read.table(paste0("IDX_k_", k, "_temp.txt"))
  #names(clusts) <- "cluster"
  #clusts$cluster
  clusts$V1
}

getW <- function (dataset, clus) {
  n <- nrow(dataset)
  ii <- seq_len(n)
  0.5 * sum(vapply(split(ii, clus), function(I) {
    xs <- dataset[I, , drop = FALSE]
    sum(dist(xs)/nrow(xs))}, 0))
}

getClusterNumber <- function (path.to.Cliq.py, k, r, m, name) {
  clus <- SNNclust(path.to.cliq.py = path.to.Cliq.py, k = k, r = r, m = m, name = name)
  max(clus)
}