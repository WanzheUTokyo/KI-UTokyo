##################################
#title: "Batch-effect corrections"
#author: "Martin Loza"
##################################

# This is the script used to correct batch-effects in the pseudo-batch test.
# To reproduce the Figures of each test, check the appropiate notebook file.

## Seurat integration
set.seed(seed) 
Anchors <- Seurat::FindIntegrationAnchors(object.list = xl, verbose = FALSE)
Seurat <- Seurat::IntegrateData(anchorset = Anchors, verbose = FALSE)
set.seed(seed)
Seurat <- RNAseqAnalysis::GetUMAP(object = Seurat, dims = dimPCA, assay = "integrated")

## Set up the Uncorrected data
set.seed(seed)
features <- Seurat::VariableFeatures(Seurat, assay = "integrated")
Uncorrected <- Reduce(merge, xl)
Uncorrected <- Uncorrected[features,]
Seurat::VariableFeatures(Uncorrected) <- features
set.seed(seed)
Uncorrected <- RNAseqAnalysis::GetUMAP(object = Uncorrected, dims = dimPCA, assay = "RNA") 
 
## Clustering the uncorrected data
set.seed(seed)
Uncorrected <- Seurat::FindNeighbors(Uncorrected, dims = 1:dimPCA, reduction = "pca", verbose = FALSE)
set.seed(seed)
Uncorrected <- Seurat::FindClusters(Uncorrected, resolution =resolution, algorithm = algorithm, verbose = FALSE)

## Add clusters to Seurat object
Seurat$seurat_clusters <- Uncorrected$seurat_clusters

## Canek 
set.seed(seed)
Canek <- RunCanek(Uncorrected, batches = "batch")
set.seed(seed)
Canek <- RNAseqAnalysis::GetUMAP(object = Canek, dims = dimPCA, assay = "Canek")
Canek$seurat_clusters <- Uncorrected$seurat_clusters

## MNN 
set.seed(seed) 
MNN <- RNAseqAnalysis::RunMNN(object = Uncorrected, batch = "batch")
set.seed(seed)
MNN <- RNAseqAnalysis::GetUMAP(object = MNN, dims = dimPCA, assay = "integrated")

## ComBat
set.seed(seed)
ComBat <- RNAseqAnalysis::RunComBat(object = Uncorrected, batch = "batch")
set.seed(seed)
ComBat <- RNAseqAnalysis::GetUMAP(object = ComBat, dims = dimPCA, assay = "integrated")

## scMerge
ks <- rep(length(unique(Uncorrected$seurat_clusters)),2) # Number of celltypes used in scMerge.
set.seed(seed)
scMerge <- RNAseqAnalysis::RunScMerge(object = Uncorrected, batch = "batch", ks = ks)
set.seed(seed)
scMerge <- RNAseqAnalysis::GetUMAP(object = scMerge, dims = dimPCA, assay = "integrated")
 
## ComBat-seq
set.seed(seed)
ComBatseq <- RNAseqAnalysis::RunComBatseq(object = Uncorrected, batch = "batch")
set.seed(seed)
ComBatseq <- RNAseqAnalysis::GetUMAP(object = ComBatseq, dims = dimPCA, assay = "integrated")

## Liger
reduction <- "Liger"
set.seed(seed)
Liger <- RNAseqAnalysis::RunLiger(object = Uncorrected, batch = "batch")
dims <- ncol(Liger[[reduction]])
set.seed(seed)
Liger <- RNAseqAnalysis::GetUMAP(object = Liger, dims = dims, reduction = reduction, PCA = FALSE, scale = FALSE)

## Scanorama
set.seed(seed)
Scanorama <- RNAseqAnalysis::RunScanorama(object = Uncorrected, batch = "batch")
set.seed(seed)
Scanorama <- RNAseqAnalysis::GetUMAP(object = Scanorama, dims = dimPCA, assay = "integrated")
 
## Harmony
set.seed(seed)
reduction <- "harmony"
Harmony <- RNAseqAnalysis::RunHarmony(object = Uncorrected, batch = "batch", dims = dimPCA)
dims <- ncol(Harmony[[reduction]])
set.seed(seed)
Harmony <- RNAseqAnalysis::GetUMAP(object = Harmony, dims = dims, reduction = reduction, PCA = FALSE, scale = FALSE)
 