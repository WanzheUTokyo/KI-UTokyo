##################################
#title: "Batch-effect corrections"
#author: "Martin Loza"
##################################

# This is the script used to correct batch-effects using different methods.
# To reproduce the Figures of each test, check the appropiate notebook file.

## Load packages
library(Canek)
#library(batchelor)
library(Seurat)
#library(sva)

## Seurat integration
set.seed(seed)
Anchors <- Seurat::FindIntegrationAnchors(object.list = xl, verbose = FALSE)
Seurat <- Seurat::IntegrateData(anchorset = Anchors, verbose = FALSE)
set.seed(seed)
Seurat <- GetUMAP(object = Seurat, dims = dimPCA, assay = "integrated")

## Set up the Uncorrected data
set.seed(seed)
features <- Seurat::VariableFeatures(Seurat, assay = "integrated")
Uncorrected <- Reduce(merge, xl)
Uncorrected <- Uncorrected[features,]
Seurat::VariableFeatures(Uncorrected) <- features
set.seed(seed)
Uncorrected <- GetUMAP(object = Uncorrected, dims = dimPCA, assay = "RNA") 
 
## Canek 
set.seed(seed)
Canek <- RunCanek(Uncorrected, batches = "batch")
set.seed(seed)
Canek <- GetUMAP(object = Canek, dims = dimPCA, assay = "Canek")

## Canek2
set.seed(seed)
Canek2 <- RunCanek(Uncorrected, batches = "batch", correctEmbeddings = TRUE)
set.seed(seed)
Canek2 <- GetUMAP(object = Canek, dims = dimPCA, assay = "Canek")

## MNN 
# 
# set.seed(seed)
# MNN <- RunMNN(object = Uncorrected, batch = "batch")
# set.seed(seed)
# MNN <- GetUMAP(object = MNN, dims = dimPCA, assay = "integrated")
# 
# ## ComBat
# set.seed(seed)
# ComBat <- RunComBat(object = Uncorrected, batch = "batch")
# set.seed(seed)
# ComBat <- GetUMAP(object = ComBat, dims = dimPCA, assay = "integrated")
# 
# ## scMerge
# set.seed(seed)
# scMerge <- RunScMerge(object = Uncorrected, batch = "batch", ks = ks)
# set.seed(seed)
# scMerge <- GetUMAP(object = scMerge, dims = dimPCA, assay = "integrated")
#  
# ## ComBat-seq
# set.seed(seed)
# ComBatseq <- RunComBatseq(object = Uncorrected, batch = "batch")
# set.seed(seed)
# ComBatseq <- GetUMAP(object = ComBatseq, dims = dimPCA, assay = "integrated")
# 
# ## Liger
# reduction <- "Liger"
# set.seed(seed)
# Liger <- RunLiger(object = Uncorrected, batch = "batch")
# dims <- ncol(Liger[[reduction]])
# set.seed(seed)
# Liger <- GetUMAP(object = Liger, dims = dims, reduction = reduction, PCA = FALSE, scale = FALSE)
# 
# ## Scanorama
# set.seed(seed)
# Scanorama <- RunScanorama(object = Uncorrected, batch = "batch")
# set.seed(seed)
# #Scanorama <- GetUMAP(object = Scanorama, dims = dimPCA, assay = "integrated")

 
## Harmony
set.seed(seed)
reduction <- "harmony"
Harmony <- RunHarmony(object = Uncorrected, batch = "batch", dims = dimPCA)
dims <- ncol(Harmony[[reduction]])
set.seed(seed)
Harmony <- GetUMAP(object = Harmony, dims = dims, reduction = reduction, PCA = FALSE, scale = FALSE)
 