# Load required libraries

#install.packages("Seurat")
library(Seurat)
#install.packages("devtools")
#devtools::install_github("jbergenstrahle/STUtility")
#library(STutility)
library(ggplot2)
library(patchwork)

#getwd()
#setwd("C:/Users/momou/OneDrive/Desktop/PhD_Erlangen/pipelines/Spatial")
#list.files()

main_dir <- "Visium_Output"

# List all sample directories
sample_dirs <- list.dirs(main_dir, recursive = FALSE)

# Initialize an empty list to store Seurat objects
samples_list <- list()

# Load each sample
for (sample_dir in sample_dirs) {
  sample_name <- basename(sample_dir)
  sample_seurat <- Load10X_Spatial(data.dir = sample_dir, assay = "Spatial")
  print(sample_seurat)
  sample_seurat$orig.ident <- sample_name 
  # filter zero expression spots
  sample_seurat <- sample_seurat[, unname(which(colSums(GetAssayData(sample_seurat)) != 0))]
  samples_list[[sample_name]] <- sample_seurat
}


# normalize the data
samples_list <- lapply(samples_list, function(x) {
  x <- SCTransform(x, assay = "Spatial", verbose = FALSE)
  return(x)
})


# Select integration features
features <- SelectIntegrationFeatures(object.list = samples_list, nfeatures = 3000)

# Prepare for integration
samples_list <- PrepSCTIntegration(object.list = samples_list, anchor.features = features)

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = samples_list, normalization.method = "SCT", anchor.features = features)

saveRDS(anchors, file = "./anchors.rds")
#anchors <- readRDS(anchors_file)
# Integrate data
seurat_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# Run PCA
seurat_integrated <- RunPCA(seurat_integrated, features = VariableFeatures(object = seurat_integrated))

# Find neighbors and clusters
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:30)
seurat_integrated <- FindClusters(seurat_integrated, resolution = 0.5)

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:30)

# Plot UMAP
DimPlot(seurat_integrated, reduction = "umap", label = TRUE)