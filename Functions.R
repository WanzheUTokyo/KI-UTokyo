## Functions##

########
## Batch-correction
########

#' RunComBatseq
#'
#' @param object A seurat object to correct batch effects.
#' @param batch Batch labels.
#' @param runningTime Return the running time.
#' @param verbose Print verbose.
#' @param ... Arguments passed to other methods.
#'
#' @return Corrected and normalized Seurat object.
#' @export
RunComBatseq <- function(object = NULL, batch = "batch", runningTime = FALSE, verbose = FALSE, ...){
  
  features <- Seurat::VariableFeatures(object)
  if(length(features) == 0){
    warning("Variable features not defined. Running 'FindVariableFeatures' function.", call. = TRUE)
    features <- Seurat::VariableFeatures(Seurat::FindVariableFeatures(object, verbose = verbose))
  }
  
  counts <- as.matrix(Seurat::GetAssayData(object, assay = "RNA", slot = "counts")[features,])
  md <- object[[]]
  
  if(!(batch %in% colnames(md)))
    stop(paste0(batch, "not found in object's metadata. Check the batch label."))
  
  time <- system.time({
    corrCounts <- sva::ComBat_seq(counts = counts, batch = md[[batch]], full_mod = FALSE)
  })
  
  object[["integrated"]] <- Seurat::CreateAssayObject(counts = corrCounts)
  Seurat::DefaultAssay(object) <- "integrated"
  
  object <- Seurat::NormalizeData(object = object, assay = "integrated", verbose = verbose, ...)
  Seurat::VariableFeatures(object) <- features
  
  if(runningTime == FALSE)
    return(object)
  else
    return(list(object = object, time = time))
}

#' RunScMerge
#'
#' @param object A seurat object to correct batch effects.
#' @param batch Batch labels.
#' @param ks A vector indicates the kmeans's K for each batch, which length needs to be the same as the number of batches.
#' @param runningTime Return the running time.
#' @param verbose Print verbose.
#' @param ... Arguments passed to other methods.
#'
#' @return Corrected Seurat object.
#' @export
RunScMerge <- function(object = NULL, batch = "batch", ks = NULL, runningTime = FALSE, verbose = FALSE, ...){
  
  features <- Seurat::VariableFeatures(object)
  if(length(features) == 0){
    warning("Variable features not defined. Running 'FindVariableFeatures' function.", call. = TRUE)
    features <- Seurat::VariableFeatures(Seurat::FindVariableFeatures(object, verbose = verbose))
  }
  
  data <- as.matrix(Seurat::GetAssayData(object, assay = "RNA", slot = "data"))
  md <- object[[]]
  
  if(!(batch %in% colnames(md)))
    stop(paste0(batch, "not found in object's metadata. Check the batch label."))
  
  if(is.null(ks)){
    nBatches <- length(unique(md[,batch]))
    ks <- rep(5, nBatches)
  }
  
  tmp <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = data, logcounts = data), colData = md)
  seg = scMerge::scSEGIndex(exprs_mat = data)
  
  time <- system.time({
    tmp <- scMerge::scMerge(sce_combine = tmp, ctl = rownames(seg), assay_name = "scMerge",
                            kmeansK = ks, batch_name = batch, plot_igraph = FALSE, verbose = FALSE, ...)
  })
  
  # Seurat assay
  corrData <- as.matrix(SummarizedExperiment::assay(tmp, "scMerge"))
  object[["integrated"]] <- Seurat::CreateAssayObject(counts = corrData)
  Seurat::DefaultAssay(object) <- "integrated"
  Seurat::VariableFeatures(object) <- features
  
  if(runningTime == FALSE)
    return(object)
  else
    return(list(object = object, time = time))
}

#' RunMNN
#'
#' @param object A seurat object to correct batch effects.
#' @param batch Batch labels.
#' @param runningTime Return the running time.
#' @param verbose Print verbose.
#' @param ... Arguments passed to other methods.
#'
#' @return Seurat object with the corrected data in the 'integrated' assay.
#' @export
RunMNN <- function(object = NULL, batch = "batch", runningTime = FALSE, verbose = FALSE, ...){
  
  features <- Seurat::VariableFeatures(object)
  if(length(features) == 0){
    warning("Variable features not defined. Running 'FindVariableFeatures' function.", call. = TRUE)
    features <- Seurat::VariableFeatures(Seurat::FindVariableFeatures(object, verbose = verbose))
  }
  
  data <- as.matrix(Seurat::GetAssayData(object, assay = "RNA", slot = "data")[features,])
  md <- object[[]]
  
  if(!(batch %in% colnames(md)))
    stop(paste0(batch, "not found in object's metadata. Check the batch labels."))
  
  time <- system.time({
    corrData <- batchelor::mnnCorrect(data, batch = md[[batch]], ...)
  })
  
  corrData <- SummarizedExperiment::assay(corrData, "corrected")
  object[["integrated"]] <- Seurat::CreateAssayObject(counts = corrData)
  Seurat::DefaultAssay(object) <- "integrated"
  Seurat::VariableFeatures(object) <- features
  
  if(runningTime == FALSE)
    return(object)
  else
    return(list(object = object, time = time))
}

#' RunScanorama
#'
#' @param object A seurat object to correct batch effects.
#' @param batch Batch labels.
#' @param runningTime Return the running time.
#' @param verbose Print verbose.
#' @param ... Arguments passed to other methods.
#'
#' @return Seurat object with the corrected data in the 'integrated' assay.
#' @export
RunScanorama <- function(object = NULL, batch = "batch", runningTime = FALSE, verbose = FALSE, ...){
  
  Scanorama <- reticulate::import("scanorama")
  datal <- list()
  genel <- list()
  features <- Seurat::VariableFeatures(object)
  
  if(length(features) == 0){
    warning("Variable features not defined. Running 'FindVariableFeatures' function.", call. = TRUE)
    features <- Seurat::VariableFeatures(Seurat::FindVariableFeatures(object, verbose = verbose))
  }
  
  if(!(batch %in% colnames(object[[]])))
    stop(paste0(batch, "not found in object's metadata. Check the batch labels."))
  
  objectl <- Seurat::SplitObject(object, split.by = batch)
  
  for(i in seq_len(length(objectl))){
    datal[[i]] <- Seurat::GetAssayData(objectl[[i]], assay = "RNA", slot = "data")[features,] # Normalized counts
    datal[[i]] <- as.matrix(datal[[i]])
    datal[[i]] <- t(datal[[i]]) # Cell x genes
    
    genel[[i]] <- features
  }
  
  time <- system.time({
    corrDatal <- Scanorama$correct(datasets_full = datal, genes_list = genel, return_dense = TRUE)
  })
  
  corrData <- Reduce(rbind, corrDatal[[1]])
  corrData <- t(corrData)
  rownames(corrData) <- corrDatal[[2]]
  colnames(corrData) <- unlist(sapply(objectl,colnames))
  
  # Same cell names as the original object
  corrData <- corrData[,colnames(object)]
  
  ## Create Seurat assay
  object[["integrated"]] <- Seurat::CreateAssayObject(counts = corrData)
  Seurat::DefaultAssay(object) <- "integrated"
  Seurat::VariableFeatures(object) <- features
  
  if(runningTime == FALSE)
    return(object)
  else
    return(list(object = object, time = time))
}

#' RunLiger
#'
#' @param object A seurat object to correct batch effects.
#' @param batch Batch labels.
#' @param k Inner dimension of factorization (number of factors)
#' @param runningTime Return the running time.
#' @param verbose Print verbose.
#' @param ... Arguments passed to other methods.
#'
#' @return Seurat object with the corrected data in the 'Liger' reduction.
#' @export
RunLiger <- function(object = NULL, batch = "batch", k = 30, runningTime = FALSE, verbose = FALSE, ...){
  
  features <- Seurat::VariableFeatures(object)
  
  if(length(features) == 0){
    warning("Variable features not defined. Running 'FindVariableFeatures' function.", call. = TRUE)
    features <- Seurat::VariableFeatures(Seurat::FindVariableFeatures(object, verbose = verbose))
  }
  
  if(!(batch %in% colnames(object[[]])))
    stop(paste0(batch, "not found in object's metadata. Check the batch labels."))
  
  tmp <- object[features,]
  
  time <- system.time({
    tmp <- Seurat::ScaleData(tmp, split.by = "batch", do.center = FALSE, verbose = verbose, ...)
    tmp <- SeuratWrappers::RunOptimizeALS(tmp, k = k, split.by = "batch", ...)
    tmp <- SeuratWrappers::RunQuantileNorm(tmp, split.by = "batch", ...)
  })
  
  object[["Liger"]] <- tmp[["iNMF"]]
  
  if(runningTime == FALSE)
    return(object)
  else
    return(list(object = object, time = time))
}

#' RunComBat
#'
#' @param object A seurat object to correct batch effects.
#' @param batch Batch labels.
#' @param runningTime Return the running time.
#' @param verbose Print verbose.
#' @param ... Arguments passed to other methods.
#'
#' @return Seurat object with the corrected data in the 'integrated' assay.
#' @export
RunComBat <- function(object = NULL, batch = "batch", runningTime = FALSE, verbose = FALSE, ...){
  
  features <- Seurat::VariableFeatures(object)
  if(length(features) == 0){
    warning("Variable features not defined. Running 'FindVariableFeatures' function.", call. = TRUE)
    features <- Seurat::VariableFeatures(Seurat::FindVariableFeatures(object, verbose = verbose))
  }
  
  data <- as.matrix(Seurat::GetAssayData(object, assay = "RNA", slot = "data")[features,])
  md <- object[[]]
  
  if(!(batch %in% colnames(md)))
    stop(paste0(batch, "not found in object's metadata. Check the batch label."))
  
  time <- system.time({
    corrData <- sva::ComBat(dat = data, batch = md[[batch]], ...)
  })
  
  object[["integrated"]] <- Seurat::CreateAssayObject(counts = corrData)
  Seurat::DefaultAssay(object) <- "integrated"
  Seurat::VariableFeatures(object) <- features
  
  if(runningTime == FALSE)
    return(object)
  else
    return(list(object = object, time = time))
}

#' RunHarmony
#'
#' @param object A seurat object to correct batch effects.
#' @param batch Batch labels.
#' @param dims Dimensions to use in the correction.
#' @param runningTime Return the running time.
#' @param verbose Print verbose.
#' @param ... Arguments passed to other methods.
#'
#' @return Seurat object with the corrected data in the 'harmony' reduction.
#' @export
RunHarmony <- function(object = NULL, batch = "batch", dims = 10, runningTime = FALSE, verbose = FALSE, ...){
  
  if(!("pca" %in% Seurat::Reductions(object))){
    if(verbose)
      print("Running PCA.")
    time <- system.time({
      object <- GetPCA(object = object, dims = dims, verbose = verbose, ...)
      object <- harmony::RunHarmony(object = object, group.by.vars = batch, dims.use = 1:dims, verbose = verbose, ...)
    })
  }else{
    time <- system.time({
      object <- harmony::RunHarmony(object = object, group.by.vars = batch, dims.use = 1:dims, verbose = verbose, ...)
    })
  }
  
  if(runningTime == FALSE)
    return(object)
  else
    return(list(object = object, time = time))
}

#' Title
#'
#' @param object A seurat object to correct batch effects.
#' @param batch Batch labels.
#' @param reduction Reduction to use.
#' @param dims Number of dimensions to use.
#' @param per Percentages of the mean batch size.
#' @param acceptance Return the acceptance rate.
#' @param verbose Print verbose.
#'
#' @return kBET mean score.
#' @export
RunKBET <- function(object = NULL, batch = "batch", reduction = "pca", dims = 10, per = 0.1, acceptance = TRUE, verbose = FALSE){
  
  md <- object[[]]
  if(!(reduction %in% Seurat::Reductions(object)))
    stop(paste0(reduction, " not found in the object's reductions."))
  
  if(!(batch %in% colnames(md)))
    stop(paste0(batch, " not found in the object's meta data."))
  
  data <- as.data.frame(Seurat::Embeddings(object = object, reduction = reduction)[,1:dims])
  meanBatch <- mean(table(md[[batch]]))
  
  scores <- lapply(per, function(p){
    k0 = floor(p*(meanBatch))
    score <- mean(kBET::kBET(df = data, batch = md[[batch]], do.pca = FALSE,
                             heuristic = FALSE, k0 = k0,
                             plot = FALSE)$stats$kBET.observed)
    return(score)
  })
  
  scores <- unlist(scores)
  scores <- mean(scores)
  
  if(acceptance)
    scores <- 1-scores
  
  return(scores)
}

#' RunSilhouette
#'
#' @param object A seurat object to correct batch effects.
#' @param batch Batch labels.
#' @param reduction Reduction to use.
#' @param dims Number of dimensions to use.
#'
#' @return Silhouette width score.
#' @export
RunSilhouette <- function(object = NULL, batch = "celltype", reduction = "pca", dims = 10){
  
  md <- object[[]]
  
  if(!(reduction %in% Seurat::Reductions(object)))
    stop(paste0(reduction, " not found in the object's reductions."))
  
  if(!(batch %in% colnames(md)))
    stop(paste0(batch, " not found in the meta data."))
  
  batch <- factor(md[[batch]])
  
  pcaData <- as.matrix(Seurat::Embeddings(object = object, reduction = reduction)[,1:dims])
  pcaData <- list(x = pcaData)
  
  score <- kBET::batch_sil(pca.data = pcaData, batch = batch, nPCs = dims)
  
  return(score)
}

########
## Dimensionality Reduction
########


#' GetPCA
#'
#' @param object Seurat object.
#' @param dims Dimensions to obtain.
#' @param verbose Print verbose.
#' @param ... Arguments passed to other methods.
#'
#' @return PCA representation.
#' @export
GetPCA <- function(object = NULL, dims = 10, verbose = FALSE, ...){
  object <- Seurat::ScaleData(object, ...)
  object <- Seurat::RunPCA(object, npcs = dims, ...)
  
  return(object)
}

#' GetUMAP
#'
#' @param object Seurat object.
#' @param dims Dimensions to use.
#' @param reduction Reduction to use.
#' @param PCA Obtain PCA.
#' @param scale Scale the data.
#' @param verbose Print verbose.
#' @param seed Set a random seed. By default, sets the seed to 42.
#' @param ... Arguments passed to other methods.
#'
#' @return UMAP representation.
#' @export
GetUMAP <- function(object = NULL, dims = 10, reduction = "pca", PCA = TRUE, scale = TRUE, seed = 42, verbose = FALSE, ...){
  if(scale)
    object <- Seurat::ScaleData(object, verbose = verbose, ...)
  if(PCA)
    object <- Seurat::RunPCA(object, npcs = dims, verbose = verbose, ...)
  
  object <- Seurat::RunUMAP(object, reduction = reduction, dims = 1:dims, seed.use = seed, verbose = verbose, ...)
  
  return(object)
}

#' SeuratPreprocessing
#'
#' @param object Seurat object.
#' @param verbose Print verbose.
#' @param ... Arguments passed to other methods.
#'
#' @return Normalized Seurat object and its variable features.
#' @export
SeuratPreprocessing <- function(object = NULL, verbose = FALSE, ...){
  object <- Seurat::NormalizeData(object, verbose = FALSE, ...)
  object <- Seurat::FindVariableFeatures(object, verbose = FALSE, ...)
  
  return(object)
}

#' Title
#'
#' @param object Object to sample
#' @param frac Fraction of samples. Must be a value between 0 and 1.
#' @param seed Seed to use on sampling
#' @param ... Arguments passed to other methods.
#'
#' @return Sampled data
#' @export
SampleData <- function(object = NULL, frac = NULL, seed = 777, ...){
  set.seed(seed)
  nCells <- ncol(object)
  samples <- floor(frac*nCells)
  idx <- sample(x = nCells, size = samples, ... )
  
  return(object[,idx])
}


########
## Visualization
########

# Function to plot the kBET and Silhouette scores
plotMetrics <- function(scoresKbet = NULL, scoresSilhouette = NULL, methods = NULL, ...){
  df <- data.frame("Silhouette" = t(scoresSilhouette), "kBET" = t(scoresKbet), Method = methods)
  p <- ggplot(df, aes(Silhouette, kBET, color = Method, label = Method)) + 
    geom_point(size = 8, alpha = 0.5) +
    geom_text_repel(size = 9, alpha = 1.0, direction = "both", ...) + 
    ylab("kBET(acceptance rate)") + theme(legend.position = "none")
  
  return(p)
}

# Plot a correction by different groups (e.g. batch, celltype, etc.)
plotCorrection <- function(object, groups = NULL, reduction = "pca", ...){
  plotls <- list()
  for(g in groups){
    plotls[[g]] <-  DimPlot(object, reduction = reduction, group.by = g, ...) 
  }
  return(plotls)
}

# Scales the plot color gradient to a given scale (limits) and colors.
scaleGradient <- function(plotls, low = "gray", high = "purple", limits = NULL){
  if(is.null(limits)){
    stop("No defined limits. Select the limits and try again.", call. = TRUE)
  }
  
  for(i in seq_len(length(plotls))){
    plotls[[i]] <- plotls[[i]] + 
      scale_color_gradient(low = low, high = high, limits = limits)  
  }
  
  return(plotls)
}

# Plot features from a Seurat objects with a given low-high color scale
plotFeatures <- function(object, features = NULL, reduction = "pca",colorLow = "gray", colorHigh = "purple", colorLimits = c(0,1), ...){
  pFeatures <- FeaturePlot(object, features = features, reduction = reduction,
                           order = TRUE, combine = FALSE)
  pFeatures <- scaleGradient(pFeatures, low = colorLow, high = colorHigh, limits = colorLimits)
  names(pFeatures) <- features
  
  return(pFeatures)
}