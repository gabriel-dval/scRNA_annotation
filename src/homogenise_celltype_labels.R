#' The purpose of this script is to test different methods for homogenising
#' the cell type labels of reference datasets to match the manual annotation 
#' of the Skull Bone Marrow (SBM) data. This has a number of uses, notably
#' for comparing the results of different datasets/annotation methods on the
#' SBM data.
#' 
#' NOTE : We will most likely be using the SMB data to match labels across 
#' datasets - this may cause a bias (overfitting). It is important to remain
#' cautious about this.
library(Seurat)
library(scmap)
library(SingleCellExperiment)
library(tidyverse) 
library(ggplot2)


###############################################################################
# Load Original Data ##########################################################
###############################################################################

# Load RData and corresponding clusters
load('../data/rdata_sbm/BoneMarrow_2024.05.22.RData')
norm_data <- Immune_norm@assays$SCT@scale.data
log1_data <- Immune_norm@assays$SCT@data
colnames(norm_data) <- Cells(Immune_norm)

# Load old cluster annotations
new.cluster.ids <- c("Neutrophil Mature",  #1
                     "Neutrophil Mature",  #2
                     "Mature B cells",  #3
                     "Neutrophil Immature",  #4
                     "Neutrophil Mature",  #5
                     "Pre-B cells",  #6
                     "Macrophages",  #7
                     "Neutrophil Immature",  #8
                     "Neutrophil Mature",  #9
                     "Neutrophil Mature", #10
                     "Neutrophil Mature", #11
                     "Neutrophil Pre/Pro", #12
                     "Neutrophil Immature", #13
                     "Neutrophil Mature", #14
                     "NKT", #15
                     "Monocytes", #16
                     "Immature B cells", #17
                     "HSC/CMP", #18
                     "Neutrophil Pre/Pro", #19
                     "pDCs", #20
                     "Neutrophil Mature", #21
                     "NK", #22
                     "CD4 T cells", #23
                     "Pre-B cells", #24
                     "Neutrophil Pre/Pro", #25
                     "CD8 T cells", #26
                     "Monocytes", #27
                     "Macrophages", #28
                     "GMP", #29
                     "cDCs", #30
                     "Monocytes", #31
                     "Neutrophil Immature", #32
                     "Neutrophil Mature", #33
                     "Neutrophil Immature", #34
                     "Unknown 1", #35
                     "MDP", #36
                     "Pre-B cells", #37
                     "Neutrophil Pre/Pro", #38
                     "Monocytes", #39
                     "Basophils",  #40
                     "Pro-B cells",  #41
                     "Mast cells",  #42
                     "Pre-B cells",  #43
                     "NKT",  #44
                     "Pax5 B cells",  #45
                     "Memory B cells",  #46
                     "Memory B cells",  #47
                     "Plasma cells",  #48
                     "Unknown 2",  #49
                     "Macrophages",#50
                     "Erythroblasts",#51
                     "Neutrophil Immature",#52
                     "Unknown 3"#53
)

# Add cluster names
names(new.cluster.ids) <- levels(Immune_norm)
Immune_norm[["old.ident"]] <- as.factor(Idents(Immune_norm) )
levels(Immune_norm$old.ident)
Immune_norm <- RenameIdents(Immune_norm, new.cluster.ids)
levels(Immune_norm@active.ident)

cluster_ids <- as.character(Idents(Immune_norm)) %>% as.data.frame()
colnames(cluster_ids) <- 'cluster_annot'
rownames(cluster_ids) <- colnames(Immune_norm)
Immune_norm <- AddMetaData(Immune_norm, 
                           metadata = cluster_ids)
Immune_norm[["RNA"]] <- as(object = Immune_norm[["RNA"]], Class = "Assay")
SaveH5Seurat(Immune_norm, filename = "manual_annot.h5Seurat")
Convert("manual_annot.h5Seurat", dest = "h5ad")






###############################################################################
# First test - use scMAP to lower resolution ##################################
###############################################################################

# Load general immune cell reference
help(scmap::indexCluster())








