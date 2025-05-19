# Script to prepare seurat objects and save the necessary metadata

library(tidyverse)
library(Seurat)
library(SeuratDisk)

###############################################################################
# SBM #########################################################################
###############################################################################

# First for SBM
load('../data/rdata_sbm/BoneMarrow_data_2025.RData') # SBM

# Load additional annotations
annots <- read.table('results/sbm_key_results.csv', sep = ';', header = T)
annots[48, 13] <- 'Correct'

# Add metadata
colnames(annots)[1] <- "seurat_clusters"

# 1. First, extract the current metadata
metadata <- Immune_norm@meta.data

# 2. Convert your annotation dataframe's cluster column to match the seurat object's type
# Check what type the seurat_clusters is in your object
print(class(metadata$seurat_clusters))
metadata$seurat_clusters <- as.integer(metadata$seurat_clusters)

# Create a named lookup vector for each annotation column
# Example for one annotation column (repeat for each column):
for(col_name in colnames(annots)[-1]) {  # Skip the first column (cluster IDs)
  # Create a named vector where names are cluster IDs and values are annotations
  annotation_vector <- setNames(
    annots[[col_name]],
    as.character(annots$seurat_clusters)
  )
  
  # Add the new column to metadata by mapping through the existing clusters
  metadata[[col_name]] <- annotation_vector[as.character(metadata$seurat_clusters)]
}

# 4. Assign the updated metadata back to the Seurat object
Immune_norm@meta.data <- metadata

# 5. Verify the new columns were added correctly
head(Immune_norm@meta.data)


# Test 
DimPlot(Immune_norm, reduction = 'umap', label = TRUE, group.by = 'SingleR_res',
        pt.size = 0.4, shuffle = T, seed = seed, label.box = T, repel = T, 
        label.size = 4, label.color = 'black', alpha = 0.75) + NoLegend() +
  ggtitle('GPTCellType corrected, GPT 4.5 Preview, detailedprompt2')


# Create table with umap coordinates
umap_reduc <- Immune_norm@reductions$umap@cell.embeddings
umap_annots <- cbind(umap_reduc, metadata)

# Change wrong col names
colnames(umap_annots)[22] <- 'GPTCelltype_new_res'

# Modify GPT to exclude numbers
GPTCelltype <- gsub("^[\\p{Z}\\s\\p{P}\\p{N}]+|[\\p{Z}\\s]+$", "", 
                    umap_annots$GPTCelltype, perl = T)
GPTCelltype <- gsub("\\x{00A0}", " ", GPTCelltype, perl = TRUE)
umap_annots$GPTCelltype <- GPTCelltype

# Save output of interest
write_csv(umap_annots, file = 'results/sbm_annotations.csv')

print(unique(GPTCelltype))
print(unique(umap_annots$Gemini))









###############################################################################
# SBM #########################################################################
###############################################################################

# First for SBM
load('../data/rdata_cp/CP_base_data.RData') # SBM

# Load additional annotations
annots <- read.table('results/cp_key_results.csv', sep = ';', header = T)
annots <- annots[,-2]

# Add metadata
colnames(annots)[1] <- "seurat_clusters"

# 1. First, extract the current metadata
metadata <- CP_Immune_norm@meta.data

# 2. Convert your annotation dataframe's cluster column to match the seurat object's type
# Check what type the seurat_clusters is in your object
print(class(metadata$seurat_clusters))
metadata$seurat_clusters <- as.integer(metadata$seurat_clusters)

# Create a named lookup vector for each annotation column
# Example for one annotation column (repeat for each column):
for(col_name in colnames(annots)[-1]) {  # Skip the first column (cluster IDs)
  # Create a named vector where names are cluster IDs and values are annotations
  annotation_vector <- setNames(
    annots[[col_name]],
    as.character(annots$seurat_clusters)
  )
  
  # Add the new column to metadata by mapping through the existing clusters
  metadata[[col_name]] <- annotation_vector[as.character(metadata$seurat_clusters)]
}

# 4. Assign the updated metadata back to the Seurat object
CP_Immune_norm@meta.data <- metadata

# 5. Verify the new columns were added correctly
head(CP_Immune_norm@meta.data)


# Test 
DimPlot(CP_Immune_norm, reduction = 'umap', label = TRUE, group.by = 'SingleR_res',
        pt.size = 0.4, shuffle = T, seed = seed, label.box = T, repel = T, 
        label.size = 4, label.color = 'black', alpha = 0.75) + NoLegend() +
  ggtitle('GPTCellType corrected, GPT 4.5 Preview, detailedprompt2')


# Create table with umap coordinates
umap_reduc <- CP_Immune_norm@reductions$umap@cell.embeddings
umap_annots <- cbind(umap_reduc, metadata)

# Save output of interest
write_csv(umap_annots, file = 'results/cp_annotations.csv')

print(unique(GPTCelltype))
print(unique(umap_annots$Gemini))






