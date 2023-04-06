# ---
# title: "Spleen data correction supplementary figure"
# author: "Martin Loza"
# ---

#This is the main workflow to reproduce the spleend data correction supplementary figure.

## Setup
library(here)
library(Canek)
library(Seurat)
options(future.globals.maxSize = 4e10)

dimPCA <- 10 # Number of PCA dimensions used in the analysis.
#per <- c(0.05, 0.15, 0.3) # Percentages of mean size used in kBET.
per <- c(0.05) # Percentages of mean size used in kBET for test
seed <- 777 # Luck seed.
batchKBET <- "batch" # Label used in kBET.
batchSilhouette <- "celltype" # Label used in Silhouette.

dataFile <- here("Data/Spleen_TM/Raw.Rds") # Input data file's path.
resultsFile <- here("Data/Sup_Figure_Spleen") # Output data file's path.

## Load data
xl <- readRDS(dataFile)
xl

## Sampling data. For tests
set.seed(seed)
xl <- lapply(xl, SampleData, frac = 0.1, seed = seed)
xl

# Data preprocessing
set.seed(seed)
xl <- lapply(xl, SeuratPreprocessing)

## Corrections
ks <- c(length(unique(xl[[1]]$celltype)), length(unique(xl[[2]]$celltype)))
set.seed(seed)
source(here("Results/CorrectData.R"), knitr::knit_global())

#We create a list with the correction results to ease downstream analyses.
datal <- list(Uncorrected = Uncorrected,
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
