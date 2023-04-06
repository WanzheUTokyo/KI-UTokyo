# ---
# title: "Figure 4. Simulations"
# author: "Martin Loza"
# ---

# This is the main workflow to reproduce the simulated data correction shown in Figure 4.

## Setup
library(here)
library(Canek)
library(Seurat)
options(future.globals.maxSize = 4e10)

dimPCA <- 10 # Number of PCA dimensions used in the analysis.
#per <- c(0.05, 0.15, 0.3) # Percentages of mean size used in kBET.
per <- c(0.05) # Percentages of mean size used in kBET for test
seed <- 666 # Luck seed.
batchKBET <- "batch" # Label used in kBET.
batchSilhouette <- "celltype" # Label used in Silhouette.

dataFile <- here("Data/Simulations") # Input data file's path.
resultsFile <- here("Data/Figure4") # Output data file's path.

##Global Functions
source(here("Functions.R"))

## Load data
xl <- readRDS(file = paste0(dataFile, "/Raw.Rds"))
GS <- readRDS(file = paste0(dataFile, "/GS_Raw.Rds"))
xl
GS

## Sampling data. For tests
set.seed(seed)
idxSample <- sample(x = ncol(GS), size = floor(0.2*ncol(GS)), replace = FALSE)

GS <- GS[,idxSample]

x <- Reduce(merge, xl)
x <- x[,idxSample]
xl <- Seurat::SplitObject(x, split.by = "batch")

GS
xl

## Data preprocessing
set.seed(seed)
xl <- lapply(xl, SeuratPreprocessing)
set.seed(seed)
GS <- SeuratPreprocessing(GS)

## Gold standard PCA and UMAP
GS <- GetUMAP(GS, dims = dimPCA)

## Corrections
ks <- c(3,2,2) # Ks for scMerge
set.seed(seed)
source(here("Results/CorrectData.R"), knitr::knit_global())

#We create a list with the correction results to ease downstream analyses.
datal <- list(GS = GS,
              Uncorrected = Uncorrected,
              Canek = Canek,
              MNN = MNN,
              Seurat = Seurat,
              scMerge = scMerge, 
              ComBat = ComBat, 
              Liger = Liger, 
              Harmony = Harmony, 
              Scanorama = Scanorama, 
              ComBatseq = ComBatseq)

## Save corrections
saveRDS(object = datal, file = paste0(resultsFile, "/datal.Rds"))

## Metrics
set.seed(seed)
source(here("Results/Metrics.R"), knitr::knit_global())

## Save scores
saveRDS(object = list(scoresKbet = scoresKbet, scoresSilhouette = scoresSilhouette), file = paste0(resultsFile, "/scores.RDS"))
        
        