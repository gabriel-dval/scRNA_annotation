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
# CP #########################################################################
###############################################################################

# First for CP
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

# Change all non-breakable space
umap_annots$Manual <- gsub("\\x{00A0}", " ", 
                           umap_annots$Manual, perl = TRUE)
umap_annots$SingleR <- gsub("\\x{00A0}", " ", 
                            umap_annots$SingleR, perl = TRUE)
umap_annots$scMayoMap <- gsub("\\x{00A0}", " ", 
                            umap_annots$scMayoMap, perl = TRUE)
umap_annots$GPTCelltype <- gsub("\\x{00A0}", " ", 
                            umap_annots$GPTCelltype, perl = TRUE)
umap_annots$Gemini <- gsub("\\x{00A0}", " ", 
                            umap_annots$Gemini, perl = TRUE)

# Save output of interest
write_csv(umap_annots, file = 'results/cp_annotations.csv')

print(unique(GPTCelltype))
print(unique(umap_annots$Gemini))










###############################################################################
# SBM_new #####################################################################
###############################################################################

# First for SBM
load('../data/test_datasets/SBM_new.RData') # SBM

# Load additional annotations
annots <- read.table('results/sbm_new_results.csv', sep = ';', header = T)

# Add metadata and change non-breakable space
colnames(annots)[1] <- "seurat_clusters"
Manual <- gsub("\\x{00A0}", " ", annots$Manual, perl = TRUE)
GPTCelltype <- gsub("\\x{00A0}", " ", annots$GPTCelltype, perl = TRUE)

annots$Manual <- Manual
annots$GPTCelltype <- GPTCelltype

# 1. First, extract the current metadata
metadata <- OCM_immune@meta.data

# 2. Convert your annotation dataframe's cluster column to match the seurat object's type
# Check what type the seurat_clusters is in your object
print(class(metadata$celltypes))
metadata$celltypes <- as.character(metadata$celltypes)

# Create a named lookup vector for each annotation column
# Example for one annotation column (repeat for each column):
for(col_name in colnames(annots)[-1]) {  # Skip the first column (cluster IDs)
  # Create a named vector where names are cluster IDs and values are annotations
  annotation_vector <- setNames(
    annots[[col_name]],
    as.character(annots$Manual)
  )
  
  # Add the new column to metadata by mapping through the existing clusters
  metadata[[col_name]] <- annotation_vector[as.character(metadata$celltypes)]
}

# 4. Assign the updated metadata back to the Seurat object
OCM_immune@meta.data <- metadata

# 5. Verify the new columns were added correctly
head(OCM_immune@meta.data)


# Test 
DimPlot(OCM_immune, reduction = 'umap', label = TRUE, group.by = 'GPTCelltype_res',
        pt.size = 0.4, shuffle = T, seed = seed, label.box = T, repel = T, 
        label.size = 4, label.color = 'black', alpha = 0.75) + NoLegend() +
  ggtitle('GPTCellType corrected, GPT 4.5 Preview, detailedprompt2')

# Create table with umap coordinates
umap_reduc <- OCM_immune@reductions$umap@cell.embeddings
umap_annots <- cbind(umap_reduc, metadata)

# Modify GPT to exclude numbers
Manual <- gsub("\\x{00A0}", " ", 
               umap_annots$Manual, perl = TRUE)
GPTCelltype <- gsub("\\x{00A0}", " ", 
                    umap_annots$GPTCelltype, perl = TRUE)
GPTCelltype_res <- gsub("\\x{00A0}", " ", 
                        umap_annots$GPTCelltype_res, perl = TRUE)

umap_annots$Manual <- Manual
umap_annots$GPTCelltype <- GPTCelltype
umap_annots$GPTCelltype_res <- GPTCelltype_res

# Save output of interest
write_csv(umap_annots, file = 'results/sbm_new_annotations.csv')

print(unique(GPTCelltype))











###############################################################################
# CP_new #####################################################################
###############################################################################

# First for CP
load('../data/test_datasets/CP_Amaia.RData')

DimPlot(Flex, reduction = 'umap', label = TRUE, group.by = 'CellGroups',
        pt.size = 0.4, shuffle = T, seed = seed, label.box = T, repel = T, 
        label.size = 4, label.color = 'black', alpha = 0.75) + NoLegend() +
  ggtitle('GPTCellType corrected, GPT 4.5 Preview, detailedprompt2')

DimPlot(Flex, reduction = 'umap', label = TRUE, group.by = 'seurat_clusters',
        pt.size = 0.4, shuffle = T, seed = seed, label.box = T, repel = T, 
        label.size = 4, label.color = 'black', alpha = 0.75) + NoLegend() +
  ggtitle('GPTCellType corrected, GPT 4.5 Preview, detailedprompt2')

# Load additional annotations
annots <- read.table('results/cp_new_results.csv', sep = ';', header = T)

# Add metadata and change non-breakable space
colnames(annots)[1] <- "seurat_clusters"
Manual <- gsub("\\x{00A0}", " ", annots$Manual, perl = TRUE)
GPTCelltype <- gsub("\\x{00A0}", " ", annots$GPTCelltype, perl = TRUE)
GPTCelltype_res <- gsub("\\x{00A0}", " ", annots$GPTCelltype_new, perl = TRUE)

annots$Manual <- Manual
annots$GPTCelltype <- GPTCelltype
annots$GPTCelltype_new <- GPTCelltype_res

# 1. First, extract the current metadata
metadata <- Flex@meta.data

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
Flex@meta.data <- metadata

# 5. Verify the new columns were added correctly
head(Flex@meta.data)


# Test 
DimPlot(Flex, reduction = 'umap', label = TRUE, group.by = 'GPTCelltype_new',
        pt.size = 0.4, shuffle = T, seed = seed, label.box = T, repel = T, 
        label.size = 4, label.color = 'black', alpha = 0.75) + NoLegend() +
  ggtitle('GPTCellType corrected, GPT 4.5 Preview, detailedprompt2')


# Create table with umap coordinates
umap_reduc <- Flex@reductions$umap@cell.embeddings
umap_annots <- cbind(umap_reduc, metadata)

# Change wrong col names
colnames(umap_annots)[20] <- 'GPTCelltype_res'
Manual <- gsub("\\x{00A0}", " ", 
               umap_annots$Manual, perl = TRUE)
GPTCelltype <- gsub("\\x{00A0}", " ", 
                        umap_annots$GPTCelltype, perl = TRUE)
GPTCelltype_res <- gsub("\\x{00A0}", " ", 
                        umap_annots$GPTCelltype_res, perl = TRUE)

umap_annots$Manual <- Manual
umap_annots$GPTCelltype <- GPTCelltype
umap_annots$GPTCelltype_res <- GPTCelltype_res

# Save output of interest
write_csv(umap_annots, file = 'results/cp_new_annotations.csv')

print(unique(Manual))
print(unique(umap_annots$GPTCelltype))









###############################################################################
# mATLAS Lung #################################################################
###############################################################################

# First for CP
load('../data/test_datasets/mATLAS_Lung_Droplet.RData')

# Load additional annotations
annots <- read.table('results/lung_key_results.csv', sep = ';', header = T)
annots[36, 4] <- 'Correct'

# Add metadata and change non-breakable space
Manual <- gsub("\\x{00A0}", " ", annots$Manual, perl = TRUE)
GPTCelltype <- gsub("\\x{00A0}", " ", annots$GPTCelltype, perl = TRUE)
GPTCelltype_res <- gsub("\\x{00A0}", " ", annots$GPTCelltype_res, perl = TRUE)

annots$Manual <- Manual
annots$GPTCelltype <- GPTCelltype
annots$GPTCelltype_res <- GPTCelltype_res

# 1. First, extract the current metadata
metadata <- toc@meta.data

# 2. Convert your annotation dataframe's cluster column to match the seurat object's type
# Check what type the seurat_clusters is in your object
print(class(toc$Cluster))

# Create a named lookup vector for each annotation column
# Example for one annotation column (repeat for each column):
for(col_name in colnames(annots)[-1]) {  # Skip the first column (cluster IDs)
  # Create a named vector where names are cluster IDs and values are annotations
  annotation_vector <- setNames(
    annots[[col_name]],
    as.character(annots$Cluster)
  )
  
  # Add the new column to metadata by mapping through the existing clusters
  metadata[[col_name]] <- annotation_vector[as.character(metadata$Cluster)]
}

# 4. Assign the updated metadata back to the Seurat object
toc@meta.data <- metadata

# 5. Verify the new columns were added correctly
head(toc@meta.data)


# Test 
seed <- 4
DimPlot(toc, reduction = 'umap', label = TRUE, group.by = 'GPTCelltype_res',
        pt.size = 0.4, shuffle = T, seed = seed, label.box = T, repel = T, 
        label.size = 4, label.color = 'black', alpha = 0.75) + NoLegend() +
  ggtitle('GPTCellType corrected, GPT 4.5 Preview, detailedprompt2')


# Create table with umap coordinates
umap_reduc <- toc@reductions$umap@cell.embeddings
umap_annots <- cbind(umap_reduc, metadata)

# Change wrong col names
Manual <- gsub("\\x{00A0}", " ", 
               umap_annots$Manual, perl = TRUE)
GPTCelltype <- gsub("\\x{00A0}", " ", 
                    umap_annots$GPTCelltype, perl = TRUE)
GPTCelltype_res <- gsub("\\x{00A0}", " ", 
                        umap_annots$GPTCelltype_res, perl = TRUE)

umap_annots$Manual <- Manual
umap_annots$GPTCelltype <- GPTCelltype
umap_annots$GPTCelltype_res <- GPTCelltype_res

# Save output of interest
write_csv(umap_annots, file = 'results/lung_annotations.csv')

print(unique(Manual))
print(unique(umap_annots$GPTCelltype))






