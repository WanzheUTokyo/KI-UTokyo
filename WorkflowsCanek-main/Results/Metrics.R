# Metrics

# Set up the data frames for results. 
scoresKbet <- data.frame(matrix(integer(), nrow = 1, ncol = length(datal)))
colnames(scoresKbet) <- names(datal)

scoresSilhouette <- data.frame(matrix(integer(), nrow = 1, ncol = length(datal)))
colnames(scoresSilhouette) <- names(datal)

## kBET
for(method in names(datal)){
  if(method == "Liger"){
    reduction <- "Liger"
    dims = ncol(datal$Liger[[reduction]])
  }else if(method == "Harmony"){
    reduction <- "harmony"
    dims = ncol(datal$Harmony[[reduction]])
  }else{
    reduction <- "pca"
    dims <- dimPCA
  }
  
  set.seed(seed)
  scoresKbet[method] <- RunKBET(object = datal[[method]], batch = batchKBET, reduction = reduction, dims = dims, per = per, acceptance = TRUE)
}

## Silhouette
for(method in names(datal)){
  if(method == "Liger"){
    reduction <- "Liger"
    dims = ncol(datal$Liger[[reduction]])
  }else if(method == "Harmony"){
    reduction <- "harmony"
    dims = ncol(datal$Harmony[[reduction]])
  }else{
    reduction <- "pca"
    dims <- dimPCA
  }
  
  set.seed(seed)
  scoresSilhouette[method] <- RunSilhouette(object = datal[[method]], batch = batchSilhouette, reduction = reduction, dims = dims)
}
