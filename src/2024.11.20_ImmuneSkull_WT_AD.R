# Mariangeles KOVACS & Joseph JOSEPHIDES
# Institut Pasteur, Paris FR
# Pipeline created: 08 Feb 2024
# Immune cells from Skull Bone Marrow
# Conditions: Wild Type (WT) & Alzheimer's disease (AD)

# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SoupX)
library(DropletUtils)
library(scDblFinder)
library(BiocParallel)

# define seed
seed = 44
# Set seed
set.seed(seed)
# Set number of cores
cores=12

# Paths to 10X data for mac users
ADSkull_file_path='/Volumes/Biclab/Olivier/mari/single_cell_skull_dec2023/counts_231205_VH00505_176_AAF37W3M5_gex_tcr_bcr_adt/AD_Skull_immune/per_sample_outs/AD_Skull_immune/count/sample_filtered_feature_bc_matrix'
WTSkull_file_path='/Volumes/Biclab/Olivier/mari/single_cell_skull_dec2023/counts_231205_VH00505_176_AAF37W3M5_gex_tcr_bcr_adt/WT_Skull_Immune/per_sample_outs/WT_Skull_immune/count/sample_filtered_feature_bc_matrix'

# For Imad only /!\
# ADSkull_file_path='Z:/Olivier/mari/single_cell_skull_dec2023/counts_231205_VH00505_176_AAF37W3M5_gex_tcr_bcr_adt/AD_Skull_immune/per_sample_outs/AD_Skull_immune/count/sample_filtered_feature_bc_matrix'
# WTSkull_file_path='Z:/Olivier/mari/single_cell_skull_dec2023/counts_231205_VH00505_176_AAF37W3M5_gex_tcr_bcr_adt/WT_Skull_Immune/per_sample_outs/WT_Skull_immune/count/sample_filtered_feature_bc_matrix'

# Load the WTSkull dataset
WTSkull.data <- Read10X(data.dir = WTSkull_file_path)
# Load the ADSkull dataset
ADSkull.data <- Read10X(data.dir = ADSkull_file_path)

# Initialize the Seurat object with the raw (non-normalized data)
WTSkull <- CreateSeuratObject(counts = WTSkull.data, project = "WTSkull", min.cells = 15, min.features = 200)
WTSkull
ADSkull <- CreateSeuratObject(counts = ADSkull.data, project = "ADSkull", min.cells = 15, min.features = 200)
ADSkull

###############################
###### Doublet exclusion ######
###############################
# Read in data as an sce object
sce_WT = Seurat::as.SingleCellExperiment(WTSkull, assay = 'RNA')
sce_AD = Seurat::as.SingleCellExperiment(ADSkull, assay = 'RNA')
# Find Doublets
doublet.score_WT = scDblFinder(sce_WT, returnType="table", BPPARAM=MulticoreParam(cores) )
doublet.score_AD = scDblFinder(sce_AD, returnType="table", BPPARAM=MulticoreParam(cores) )
# Get results
results_WT = data.frame("Barcode" = rownames(doublet.score_WT), "scDblFinder_DropletType" = doublet.score_WT$class, "scDblFinder_Score" = doublet.score_WT$score)
results_AD = data.frame("Barcode" = rownames(doublet.score_AD), "scDblFinder_DropletType" = doublet.score_AD$class, "scDblFinder_Score" = doublet.score_AD$score)
# Add the results in the metadata
WTSkull <- AddMetaData(object = WTSkull, metadata = doublet.score_WT[colnames(WTSkull),"score"], col.name = "doublet.score")
WTSkull <- AddMetaData(object = WTSkull, metadata = doublet.score_WT[colnames(WTSkull),"class"], col.name = "doublet.class")
ADSkull <- AddMetaData(object = ADSkull, metadata = doublet.score_AD[colnames(ADSkull),"score"], col.name = "doublet.score")
ADSkull <- AddMetaData(object = ADSkull, metadata = doublet.score_AD[colnames(ADSkull),"class"], col.name = "doublet.class")
# Show number of doublets and singlets
table(WTSkull$doublet.class)
table(ADSkull$doublet.class)
# Remove cells based on doublet criteria
WTSkull = subset(WTSkull, subset = doublet.class == 'singlet')
ADSkull = subset(ADSkull, subset = doublet.class == 'singlet')

########################################################
######################## SOUP X ########################
######## Removal of cell free mRNA contamination #######
########################################################
source_names = c("WT_Skull_Immune", "AD_Skull_immune")
list_of_counts = paste0('/Volumes/Biclab/Olivier/mari/single_cell_skull_dec2023/counts_231205_VH00505_176_AAF37W3M5_gex_tcr_bcr_adt/',source_names,'/per_sample_outs/',source_names,'/count/sample_filtered_feature_bc_matrix.h5')
list_of_droplets = paste0('/Volumes/Biclab/Olivier/mari/single_cell_skull_dec2023/counts_231205_VH00505_176_AAF37W3M5_gex_tcr_bcr_adt/',source_names,'/multi/count/raw_feature_bc_matrix.h5')
list_of_Seurat_objects_for_later = c()
for (source in 1:length(source_names)){
  print (source_names[source])
  toc_file = list_of_counts[source]
  tod_file = list_of_droplets[source]
  # Load data
  toc = Read10X_h5(toc_file,use.names = T)
  tod = Read10X_h5(tod_file,use.names = T)
  # make a soup channel
  sc = SoupChannel(tod, toc)
  # Define marker genes
  toc = CreateSeuratObject(counts = toc, project = source_names[source])
  toc = SCTransform(toc, verbose = F)
  toc = RunPCA(toc, verbose = F)
  toc = RunUMAP(toc, dims = 1:40, verbose = F)
  toc = FindNeighbors(toc, dims = 1:33, verbose = F)
  toc = FindClusters(toc, verbose = F, resolution = 2)
  # Add clustering
  meta = toc@meta.data
  umap = toc@reductions$umap@cell.embeddings
  sc = setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
  sc = setDR(sc, umap)
  # head(meta)
  # Run the main SoupX function: calculate ambient RNA profile
  sc = autoEstCont(sc, doPlot = T, maxMarkers = 150)
  plotMarkerDistribution(sc)
  # plotChangeMap(sc, out, "Iglv1")
  # Save
  out = adjustCounts(sc, roundToInt = T)
  # DropletUtils:::write10xCounts(paste0("/Volumes/Biclab/Joseph/2024_skull_bone_marrow/Code shared with Mari/data/",source_names[source]), out, overwrite = T)
  # srat = CreateSeuratObject(out)
  srat = CreateSeuratObject(out)
  list_of_Seurat_objects_for_later = c(list_of_Seurat_objects_for_later, srat)
}
# The plotMarkerDistribution plot shows the distribution of log10 ratios of observed counts to 
# expected if the cell contained nothing but soup. A guess at which cells definitely express each
# gene is made and the background contamination is calculated. The red line shows the
# global estimate (i.e., assuming the same contamination fraction for all cells) of
# the contamination fraction using just that gene. This “guessing” is done using the
# quickMarkers function to find marker genes of each cluster (see “Automatic method” section).

WTSkull = list_of_Seurat_objects_for_later[[1]]
ADSkull = list_of_Seurat_objects_for_later[[2]]

WTSkull$orig.ident = 'WT Skull'
ADSkull$orig.ident = 'AD Skull'

#################################
######## Seurat Analysis ########
#################################
# Calculate mitochondrial percentage 
WTSkull$mitoPercent <- PercentageFeatureSet(WTSkull, pattern='^mt-') 
ADSkull$mitoPercent <- PercentageFeatureSet(ADSkull, pattern='^mt-')
# QC metrics
VlnPlot(WTSkull, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)
VlnPlot(ADSkull, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)
# FeatureScatter QC
plot1 <- FeatureScatter(WTSkull, feature1 = "nCount_RNA", feature2 = "mitoPercent")
plot2 <- FeatureScatter(WTSkull, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
  geom_smooth(method = 'lm')
plot3 <- FeatureScatter(ADSkull, feature1 = "nCount_RNA", feature2 = "mitoPercent")
plot4 <- FeatureScatter(ADSkull, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
  geom_smooth(method = 'lm')
plot1 + plot2 + plot3 + plot4
# Filtering based on QC
WTSkull <- subset(WTSkull, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & mitoPercent < 4 & nCount_RNA < 20000) 
ADSkull <- subset(ADSkull, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & mitoPercent < 4 & nCount_RNA < 20000) 
# QC metrics (after filtering)
VlnPlot(WTSkull, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)
VlnPlot(ADSkull, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)
# Normalise the data
WTSkull_norm <- SCTransform(WTSkull, vars.to.regress = "mitoPercent")
ADSkull_norm <- SCTransform(ADSkull, vars.to.regress = "mitoPercent")

# Identify Outliers
# WTSkull_norm <- FindVariableFeatures(WTSkull_norm, selection.method = "vst", nfeatures = 1000)
# ADSkull_norm <- FindVariableFeatures(ADSkull_norm, selection.method = "vst", nfeatures = 1000)

# Identify the 10 most highly variable genes
top10_WTSkull <- head(VariableFeatures(WTSkull_norm), 10)
top10_ADSkull <- head(VariableFeatures(ADSkull_norm), 10)
plot1 <- VariableFeaturePlot(WTSkull_norm)
plot2 <- LabelPoints(plot = plot1, points = top10_WTSkull, repel = TRUE)
plot3 <- VariableFeaturePlot(ADSkull_norm)
plot4 <- LabelPoints(plot = plot1, points = top10_ADSkull, repel = TRUE)
plot1 + plot3
plot2 + plot4

# Scale data (prepare for PCA)
all.genes <- rownames(WTSkull_norm)
WTSkull_norm <- ScaleData(WTSkull_norm, features = all.genes)
all.genes <- rownames(ADSkull_norm)
ADSkull_norm <- ScaleData(ADSkull_norm, features = all.genes)

# Run PCA
WTSkull_norm <- RunPCA(WTSkull_norm, features = VariableFeatures(object = WTSkull_norm))
ADSkull_norm <- RunPCA(ADSkull_norm, features = VariableFeatures(object = ADSkull_norm))

# Examine PCA results
print(WTSkull_norm[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(WTSkull_norm, dims = 1:2, reduction = "pca")
print(ADSkull_norm[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(ADSkull_norm, dims = 1:2, reduction = "pca")

# Heatmap to see primary sources of heterogeneity
DimHeatmap(WTSkull_norm, dims = 1:2, cells = 500, balanced = TRUE)
DimHeatmap(ADSkull_norm, dims = 1:2, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
ElbowPlot(WTSkull_norm, ndims = 50)
ElbowPlot(ADSkull_norm, ndims = 50)

# Cluster the cells at elbow point
WTSkull_norm <- FindNeighbors(WTSkull_norm, dims = 1:30)
WTSkull_norm <- FindClusters(WTSkull_norm, resolution = 1.5)
ADSkull_norm <- FindNeighbors(ADSkull_norm, dims = 1:30)
ADSkull_norm <- FindClusters(ADSkull_norm, resolution = 1.5)

# Look at cluster IDs of the first 5 cells
head(Idents(WTSkull_norm), 5)
head(Idents(ADSkull_norm), 5)

# Visualize PCA
DimPlot(WTSkull_norm, reduction = "pca")
DimPlot(ADSkull_norm, reduction = "pca")

# non-linear dimensional reduction with UMAP
WTSkull_norm <- RunUMAP(WTSkull_norm, dims = 1:40)
ADSkull_norm <- RunUMAP(ADSkull_norm, dims = 1:40)
DimPlot(WTSkull_norm, reduction = "umap")
DimPlot(ADSkull_norm, reduction = "umap")

# Optional: save 
# saveRDS(WTSkull_norm, file = "../output/WTSkull_norm.rds")
# saveRDS(ADSkull_norm, file = "../output/ADSkull_norm.rds")

# Finding differentially expressed features (cluster biomarkers)
WTSkull_norm.markers <- FindAllMarkers(WTSkull_norm, only.pos = TRUE)
WTSkull_norm.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
ADSkull_norm.markers <- FindAllMarkers(ADSkull_norm, only.pos = TRUE)
ADSkull_norm.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# find cluster 2 markers
cluster2.markers <- FindMarkers(WTSkull_norm, ident.1 = 3)
head(cluster2.markers, n = 5)
cluster2.markers <- FindMarkers(ADSkull_norm, ident.1 = 3)
head(cluster2.markers, n = 5)

# Visualizing marker expression
VlnPlot(WTSkull_norm, features = c("Chil3", "Ltf"))
VlnPlot(ADSkull_norm, features = c("Chil3", "Ltf"))

# Plot cluster 1,2,3 specific markers (2 genes for each cluster)
FeaturePlot(WTSkull_norm, features = c("Chil3", "Ltf","Ccl6","Csf3r","H2-Eb1","H2-Aa"))
FeaturePlot(ADSkull_norm, features = c("Chil3","Ltf","Mmp8","Mmp9","Ms4a6c","S100a4"))

## Neutrophils
### Pre
FeaturePlot(WTSkull_norm, features = c("Ptma", "Rpl3"), blend = TRUE)
### Pro
FeaturePlot(WTSkull_norm, features = c("Mki67", "Top2a"), blend = TRUE)
### Immature
FeaturePlot(WTSkull_norm, features = c("Ngp", "Ltf"), blend = TRUE)
### Mature
FeaturePlot(WTSkull_norm, features = c("Ccl6", "Cxcr2"), blend = TRUE)

# Generate an expression heatmap for given cells and features
top10 = WTSkull_norm.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()
DoHeatmap(ADSkull_norm, features = top10$gene) + NoLegend()
top10 = ADSkull_norm.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()
DoHeatmap(ADSkull_norm, features = top10$gene) + NoLegend()

####################################################
###### MERGE BOTH DATASETS AND REDO THE ABOVE ######
####################################################

Immune <- merge(WTSkull, y = ADSkull,
                project = 'Immune')
Immune

# Join layers (collapses the individual datasets together and recreates the original counts and data layer)
Immune <- JoinLayers(Immune)

# Normalize
Immune_norm <- SCTransform(Immune, vars.to.regress = "mitoPercent", seed.use = seed,
                           assay = 'RNA', variable.features.n = NULL,
                           variable.features.rv.th = 1.3)
# perform standard workflow steps to figure out if we see any batch effects
integrationfeatures = SelectIntegrationFeatures(c(WTSkull_norm,ADSkull_norm))
VariableFeatures(Immune_norm) = integrationfeatures

# Identify the 10 most highly variable genes
top10_immune <- head(VariableFeatures(Immune_norm), 10)
plot <- LabelPoints(plot = plot1, points = top10_immune, repel = TRUE)
plot

# Scale data (prepare for PCA)
Immune_norm <- ScaleData(Immune_norm)
# Run PCA
Immune_norm <- RunPCA(Immune_norm, features = VariableFeatures(object = Immune_norm))

# Examine PCA results
print(Immune_norm[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Immune_norm, dims = 1:2, reduction = "pca")

# Heatmap to see primary sources of heterogeneity of PC1 and 2
DimHeatmap(Immune_norm, dims = 1:2, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
ElbowPlot(Immune_norm, ndims = 50)

# Cluster the cells at elbow point
Immune_norm <- FindNeighbors(Immune_norm, dims = 1:32)
Immune_norm <- FindClusters(Immune_norm, resolution = 2.5, algorithm = 4, random.seed = seed)

# Look at cluster IDs of the first 5 cells
head(Idents(Immune_norm), 5)

# Visualize PCA
DimPlot(Immune_norm, reduction = "pca")

# non-linear dimensional reduction with UMAP
Immune_norm <- RunUMAP(Immune_norm, dims = 1:50)
DimPlot(Immune_norm, reduction = "umap", label = TRUE) #, cols = rainbow(50)

# First, prepare object to run differential expression on SCT assay with multiple models
Immune_norm <- PrepSCTFindMarkers(object = Immune_norm)

# Now, Find markers
Immune_norm.markers <- FindAllMarkers(Immune_norm, only.pos = TRUE)
markers = Immune_norm.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
head(markers)
# print top gene marker for each cluster
for (i in c(1:length(levels(Immune_norm) ) ) ){
  Cluster_markers = markers %>% filter(cluster == i)
  print(paste("Cluster ", as.character(i),":") )
  print(paste(Cluster_markers$gene[1:25]))
}
# Plot some marker genes
# Plasma
FeaturePlot(Immune_norm, features = c("Jchain", "Xbp1"), blend = TRUE)
# T cells
FeaturePlot(Immune_norm, features = c("Cd3e", "Cd3g"), blend = TRUE)
## CD8
FeaturePlot(Immune_norm, features = c("Cd8a", "Cd8b1"), blend = TRUE)
# NKT
FeaturePlot(Immune_norm, features = c("Il2rb", "Ccl5"), blend = TRUE)
# NK
FeaturePlot(Immune_norm, features = c("Gzma", "Klrb1c"), blend = TRUE)
# Macrophages + Monocytes
FeaturePlot(Immune_norm, features = c("Ly6c2", "Ms4a6c"), blend = TRUE)
# Monocytes
FeaturePlot(Immune_norm, features = c("Hopx", "Ccr2"), blend = TRUE)
FeaturePlot(Immune_norm, features = c("Bst2", "Cd14"), blend = TRUE)
# Macrophages
FeaturePlot(Immune_norm, features = c("Adgre1", "Lyz2"), blend = TRUE)
FeaturePlot(Immune_norm, features = c("Csf1r", "Fcgr1"), blend = TRUE)
# Dentritic
## Conventional 
FeaturePlot(Immune_norm, features = c("Clec9a", "Ifi205"), blend = TRUE)
## Plasmocytoid
FeaturePlot(Immune_norm, features = c("Siglech"), blend = F)
# Granulocytes
## Neutrophiles
### Pre
FeaturePlot(Immune_norm, features = c("Ptma", "Rpl3"), blend = TRUE)
### Pro
FeaturePlot(Immune_norm, features = c("Mki67", "Top2a"), blend = TRUE)
### Immature
FeaturePlot(Immune_norm, features = c("Ngp", "Ltf"), blend = TRUE)
### Mature
FeaturePlot(Immune_norm, features = c("Ccl6", "Cxcr2"), blend = TRUE)
## Eosinophils
FeaturePlot(Immune_norm, features = c("Siglecf", 'Foxp1'), blend = T)
## Basophils
FeaturePlot(Immune_norm, features = c("Fcer1a"), blend = F)
## Mast cells
FeaturePlot(Immune_norm, features = c("Kit"), blend = F)
# Erythroblasts
FeaturePlot(Immune_norm, features = c("Hba-a2", "Aqp1"), blend = TRUE)
# Myeloid Progenitor
FeaturePlot(Immune_norm, features = c("Gata2", "Itga2b"), blend = TRUE)
# G1/S progression
FeaturePlot(Immune_norm, features = c("Cdk4",'Cdk6'), blend = TRUE)
# S-phase markers
FeaturePlot(Immune_norm, features = c("Cdk2", "Rif1"), blend = TRUE)
# G2/M
FeaturePlot(Immune_norm, features = c("Cdk10"), blend = F)
# Immune marker (general)
FeaturePlot(Immune_norm, features = c("Ptprc"), blend = F)
# Fibroblasts
FeaturePlot(Immune_norm, features = c("Col1a1",'Col5a1'), blend = T)
# Endothelial cells
FeaturePlot(Immune_norm, features = c("Pecam1",'Tie1'), blend = T)
# Pericyte cells
FeaturePlot(Immune_norm, features = c("Mapk8",'Vcam1'), blend = T)
# Glutamatergic neurons
FeaturePlot(Immune_norm, features = c("Slc32a1"), blend = F)
# B cells
FeaturePlot(Immune_norm, features = c("Cd19",'Cd79a'), blend = T)
## B cells : pre-B
FeaturePlot(Immune_norm, features = c("Dnajc7"), blend = F)
## B cells : pro-B
FeaturePlot(Immune_norm, features = c("Vpreb1","Vpreb2"), blend = T)
## B cells : Immature
FeaturePlot(Immune_norm, features = c("Ms4a1"), blend = F)
## B cells : Mature
FeaturePlot(Immune_norm, features = c("Cd74","Fcer2a"), blend = T)
## Memory T
FeaturePlot(Immune_norm, features = c("Rora","Cxcr6"), blend = T)
## Memory T
FeaturePlot(Immune_norm, features = c("Rora","Cxcr6"), blend = T)
# Monocyte-derived Dendritic cells (Mo-DC)
FeaturePlot(Immune_norm, features = c("F13a1","Zeb2"), blend = T)
FeaturePlot(Immune_norm, features = c("F13a1","Zeb2"), blend = T)
FeaturePlot(Immune_norm, features = c("Slc8a1","Pid1"), blend = T)
# Apoptosis
FeaturePlot(Immune_norm, features = c("Fas","Bcl2"), blend = T)
# Apopotosis Inhibitor, promotes cell proliferation : PMC7801710 + genecard + Mari
FeaturePlot(Immune_norm, features = c("Birc5","Bax"), blend = T)
FeaturePlot(Immune_norm, features = c("Casp3","Casp8"), blend = T)
FeaturePlot(Immune_norm, features = c("Fas","Xiap"), blend = T)
FeaturePlot(Immune_norm, features = c("Tnfsf10","Cflar"), blend = T)
FeaturePlot(Immune_norm, features = c("Casp3","Bcl2"), blend = T)
# Lysis (genecard)
FeaturePlot(Immune_norm, features = c("Klrk1","Icam1"), blend = T)
FeaturePlot(Immune_norm, features = c("Ncr1","Cd55"), blend = T)
# Activated immune response
FeaturePlot(Immune_norm, features = c("Pik3cd","Stat5b"), blend = T)
FeaturePlot(Immune_norm, features = c("Rag1","Rag2"), blend = T)
FeaturePlot(Immune_norm, features = c("Mapk1","Plcg2"), blend = T)
# Mitochondrial
FeaturePlot(Immune_norm, features = c("Polg","Surf1"), blend = T)
FeaturePlot(Immune_norm, features = c("Fis1","Mfn2"), blend = T)
# T cell activation
FeaturePlot(Immune_norm, features = c("Raf1","Jun"), blend = T)
# B cell activation
FeaturePlot(Immune_norm, features = c("Map2k1","Lyn"), blend = T)
FeaturePlot(Immune_norm, features = c("Itpr1","Grb2"), blend = T)

# Asked by Mari
FeaturePlot(Immune_norm, features = c("Ptprc","Ifng","Cd8a","Itga4"), blend = F)
FeaturePlot(Immune_norm, features = c("Cd34"), blend = F)
FeaturePlot(Immune_norm, features = c("Ly6g"), blend = F)

FeaturePlot(Immune_norm, features = c("Il17ra"), blend = F)
FeaturePlot(Immune_norm, features = c("Cd19","Cd79a","Cd79b","Ighd","Ighm"), blend = F)
FeaturePlot(Immune_norm, features = c("Adgre1","Csf1r","Ly6c2"), blend = F)

FeaturePlot(Immune_norm, features = c("Adgre1","Csf1r","Ly6c2"), group_by("orig.ident"), blend = F)


# Violin Plot
VlnPlot(ADSkull_norm, features = c('Cd74','Fcer2a'))

## Define marker genes for each cell type
marker_genes = c("Ptprc",                           # General immune cell marker
                 "Jchain","Xbp1",                   # Plasma
                 "Cd3d","Cd3e","Cd3g",              # General T-Cells
                 "Cd8a","Cd8b1",                    ## T-Cells : Cd8
                 "Cd4",                             ## T-Cells : Cd4
                 "Il2rb","Ccl5",                    # NKT / NK
                 "Gzma", "Klrb1c",                  # NK
                 "Ms4a6c","Ly6c2",                  # Monocytes and Macrophages
                 "Hopx","Ccr2","Fn1","Cd14","Bst2", # Monocytes
                 "Adgre1", "Lyz2","Csf1r","Fcgr1",  # Macropgages
                 "Ifi205","Cd74",         # Dentritic Conventional
                 "Siglech",'Ccr9','Pacsin1',        # Dentritic Plasmocytoid
                 "Mmp8","S100a9","S100a8","Ly6g",   # Neutrophils
                 "Ngp",                             ## All neutrophils except Mature
                 "Ptma", "Rpl3",'Cxcr4',            ## Pre
                 "Mki67", "Top2a","Cebpe",          ## Pro
                 "Ltf",                             ## Immature
                 "Ccl6", "Cxcr2",                   ## Mature
                 "Siglecf",                         # Eosinophils
                 "Fcer1a",                          # Basophils
                 "Kit",'Itgb7','Il1rl1',            # Mast Cells
                 "Hba-a2","Hbb-bs","Aqp1",          # Erythropblasts
                 "Sox4","Gata2","Itga2b",           # Myeloid Progenitors
                 "Cd19","Cd79a","Cd79b",            # B-cells
                 "Dnajc7",                          ## Pre-B
                 "Vpreb1","Vpreb2",                 ## Pro-B
                 "Ms4a1",                           ## Immature B
                 "Fcer2a",                          ## Mature B
                 "Cdk4",'Cdk6',                     # G1/S progression
                 "Cdk2",                            # Replicating cells
                 "Pax5",
                 "Ly6a",                            # HSCs
                 "Cd34",                            # CMP
                 "Fcgr3",                           # GMP
                 "Flt3","Itgax",                    # MDP
                 "Xcr1","Itgae","Naaa","Clec9a",    # cDC1
                 "Sirpa","Cd209a"                   # cDC2
)

# DO NOT RUN:
# rm(ADSkull, WTSkull, ADSkull_norm, WTSkull_norm, ADSkull_norm.markers, WTSkull_norm.markers,
#    doublet.score_AD, doublet.score_WT, list_of_Seurat_objects_for_later, meta, out, plot,
#    plot1, plot2, plot3, plot4, results_AD, results_WT, sc, sce_AD, sce_WT, srat, toc, tod,
#    top10, umap, WTSkull.data, ADSkull.data, all.genes, i, Cluster_markers, cluster2.markers)
# save.image(file='/Volumes/Biclab/Joseph/2024_skull_bone_marrow/Code shared with Mari/BoneMarrow_2024.05.22.RData')

# Load most recent Rdata
load('../data/rdata_sbm/BoneMarrow_2024.11.15.RData')

# Immune_norm <- FindClusters(Immune_norm, resolution = 2.5, algorithm = 4, random.seed = seed)
DimPlot(Immune_norm, reduction = "umap", label = TRUE, pt.size = 0.3, cols = viridis::turbo(length(levels(Immune_norm)), begin = 0.2))
# Rename based on gene expression
new.cluster.ids <- c("Neutro. Mature",  #1
                     "Neutro. Mature",  #2
                     "Mature B cells",  #3
                     "Neutro. Immature",  #4
                     "Neutro. Mature",  #5
                     "Pre-B cells",  #6
                     "Macrophages",  #7
                     "Neutro. Immature",  #8
                     "Neutro. Mature",  #9
                     "Neutro. Mature", #10
                     "Neutro. Mature", #11
                     "Neutro. Pre/Pro", #12
                     "Neutro. Immature", #13
                     "Neutro. Mature", #14
                     "NKT", #15
                     "Monocytes", #16
                     "Immature B cells", #17
                     "HSC/CMP", #18
                     "Neutro. Pre/Pro", #19
                     "pDCs", #20
                     "Neutro. Mature", #21
                     "NK", #22
                     "CD4 T cells", #23
                     "Pre-B cells", #24
                     "Neutro. Pre/Pro", #25
                     "CD8 T cells", #26
                     "Monocytes", #27
                     "Macrophages", #28
                     "GMP", #29
                     "cDCs", #30
                     "Monocytes", #31
                     "Neutro. Immature", #32
                     "Neutro. Mature", #33
                     "Neutro. Immature", #34
                     "Unknown 1", #35
                     "MDP", #36
                     "Pre-B cells", #37
                     "Neutro. Pre/Pro", #38
                     "Monocytes", #39
                     "Basophils",  #40
                     "Pro-B cells",  #41
                     "Mast cells",  #42
                     "Pre-B cells",  #43
                     "NKT",  #44
                     "Pax5 B cells",  #45
                     "MBC",  #46
                     "MBC",  #47
                     "Plasma cells",  #48
                     "Unknown 2",  #49
                     "Macrophages",#50
                     "Erythroblasts",#51
                     "Neutro. Immature",#52
                     "Unknown 3"#53
                     )

# Add cluster names
names(new.cluster.ids) <- levels(Immune_norm)
Immune_norm[["old.ident"]] <- as.factor(Idents(Immune_norm) )
levels(Immune_norm$old.ident)
Immune_norm <- RenameIdents(Immune_norm, new.cluster.ids)
levels(Immune_norm@active.ident)

DimPlot(Immune_norm, reduction = "umap", label = TRUE, pt.size = 0.4, shuffle = T, seed = seed, label.box = T, repel = T, label.size = 4, label.color = 'black', alpha = 0.75) + NoLegend()


# For visualisation purposes, add cluster number to cell-type name
celltypes_cluster_names = paste(Immune_norm@active.ident,Immune_norm$seurat_clusters,sep = '_')
Immune_norm$CTCN = celltypes_cluster_names

# DotPlot(Immune_norm, features = marker_genes, group.by = 'CTCN') + RotatedAxis()
DimPlot(Immune_norm, reduction = "umap", label = TRUE, pt.size = 0.3, group.by = 'CTCN') + NoLegend()

DoHeatmap(subset(Immune_norm, downsample = 200), features = marker_genes, group.by = 'seurat_clusters')  + NoLegend()
# DotPlot(subset(Immune_norm, downsample = 200), features = marker_genes, split.by = 'seurat_clusters',
        # cols = sample(viridis::turbo(length(levels(Immune_norm$seurat_clusters))) )) + RotatedAxis() + NoLegend()
DimPlot(Immune_norm, reduction = "umap", label = TRUE)


# Immune_norm@active.ident = factor(Immune_norm@active.ident,levels = sort(levels(Immune_norm@active.ident)) )

cell_colors <- c(
  "Pax5 B cells" = "#22ace0",                      # Light green for B-like cells
  "Basophils" = "#a869f5",                   # Light orange for Basophils
  "CD4 T cells" = "#629c7b",                 # Light blue for CD4 T cells
  "CD8 T cells" = "#7ed9a5",                 # Light green for CD8 T cells
  "cDCs" = "#c9bc2c",                        # Light yellow for cDCs
  "Erythroblasts" = "#cb1818",               # Light pink for Erythroblasts
  "HSC/CMP" = "#abd179",                      # Light cyan for HSPCs
  "MDP" = "#82b342",  
  "GMP" = "#73bf0f",  
  "Immature B cells" = "#fb6a4a",            # Light blue for Immature B cells
  "Macrophages" = "#69afc9",                 # Light red for Macrophages
  "Mast cells" = "#1882cc",                  # Light green for Mast Cells
  "Mature B cells" = "#cb1818",
  "MBC" = "#cb3c18",                         # Dark blue for Mature B cells
  "Monocytes" = "#ced1a5",                   # Light orange for Monocytes
  "Neutro. Immature" = "#fb6a4a",            # Light red for Neutro. Immature
  "Neutro. Mature" = "#3e7e96",               # Medium red for Neutro. Mature
  "Neutro. Pre/Pro" = "#bad2db",             # Dark red for Neutro. Pre/Pro/Pro
  "NK" = "#306948",                          # Light tan for NK cells
  "NKT" = "#c5f0d7",                         # Light salmon for NKT cells
  "pDCs" = "#dceb6e",                        # Light olive for pDCs
  "Plasma cells" = "#c947a7",                      # Light olive green for Plasma cells
  "Pre-B cells" = "#ff5972",                 # Very light blue for Pre-B cells
  "Pro-B cells" = "#fcae91",
  "Unknown 1" = "#87868f",
  "Unknown 2" = "#a8a6bd",
  "Unknown 3" = "#a8a6ad"
)


# Plot cell types with new colours
DimPlot(Immune_norm, reduction = "umap", label = TRUE, pt.size = 0.4, cols = cell_colors, shuffle = T, seed = seed, label.box = T, repel = T, label.size = 4, label.color = 'black', alpha = 0.75) + NoLegend()
# Plot cluster numbers (for reference)
DimPlot(Immune_norm, reduction = "umap", label = TRUE, pt.size = 0.5, cols = viridis::turbo(length(levels(Immune_norm$old.ident)), begin = 0.1, end=0.9), shuffle = T, seed = seed, label.box = T, repel = T, label.size = 3.5, label.color = 'black', alpha = 0.75, group.by = 'old.ident') + NoLegend()
# Plot origin (Ad or WT)
DimPlot(Immune_norm, reduction = "umap", label = F, pt.size = 0.3, cols =  viridis::cividis(2, begin = 0.1, end=0.9), shuffle = T, seed = seed, label.box = F, repel = T, alpha = 0.75, group.by = 'orig.ident')

# RIDGE PLOR FOR A GENE (OR MANY)
RidgePlot(Immune_norm, features = c("Cdk4"), ncol = 1)

# Dot plot for the marker genes we used
DotPlot(Immune_norm, features = marker_genes) + RotatedAxis()
# Same but with heatmap
DoHeatmap(subset(Immune_norm, downsample = 100), features = marker_genes, size = 3) + scale_fill_gradient2( low = rev(c('#2166ac')), mid = "white", high = rev(c('#b2182b')), midpoint = 0, guide = "colourbar", aesthetics = "fill")

# Check mitochondrial percentage per cell type and cluster
VlnPlot(Immune_norm, features = "mitoPercent", split.by = "ident") + NoLegend()
VlnPlot(Immune_norm, features = "mitoPercent", group.by = "seurat_clusters")+ NoLegend()
# Check RNA counts
VlnPlot(Immune_norm, features = "nCount_SCT", split.by = "ident") + NoLegend()
VlnPlot(Immune_norm, features = "nCount_SCT", group.by = "seurat_clusters")+ NoLegend()
# Check number of genes (=features)
VlnPlot(Immune_norm, features = "nFeature_SCT", split.by = "ident") + NoLegend()
VlnPlot(Immune_norm, features = "nFeature_SCT", group.by = "seurat_clusters")+ NoLegend()

##################################
###### Barplot with counts #######
##################################
testtt = data.frame( Immune_norm$orig.ident, Immune_norm@active.ident )
colnames(testtt) = c('Condition','Cell_Type')
testtt$Cell_Type <- factor(testtt$Cell_Type, levels = sort(levels(testtt$Cell_Type)))

# Downsample
original_counts <- table(testtt$Condition)
min_count <- min(original_counts)
# Downsampling
downsampled_df <- testtt %>%
  group_by(Condition) %>%
  slice_sample(n = min_count, replace = FALSE) %>%
  ungroup()

library(ggstatsplot)

p <- ggbarstats(downsampled_df, x = Condition, y = Cell_Type,
                bf.message=F, label="both")
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=15))


p <- ggbarstats(downsampled_df, x = Condition, y = Cell_Type,
                bf.message=F, label="both")
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, size=15))  +
  scale_fill_viridis_d(option='cividis', begin=0.1, end=0.9, direction = -1)

# Subset the data for specific clusters
selected_clusters <- c(
  "Mature B cells", 
  "MBC", 
  "Plasma cells", 
  "Pre-B cells", 
  "Pro-B cells", 
  "Immature B cells", 
  "Pax5 B cells"
)
filtered_df <- testtt %>%
  filter(Cell_Type %in% selected_clusters)

# Downsampling for balanced data (optional)
original_counts <- table(filtered_df$Condition)
min_count <- min(original_counts)

downsampled_df <- filtered_df %>%
  group_by(Condition) %>%
  slice_sample(n = min_count, replace = FALSE) %>%
  ungroup()

# Calculate scaled frequencies (0-1 range)
relative_frequencies <- downsampled_df %>%
  group_by(Condition, Cell_Type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Condition) %>%
  mutate(freq_scaled = count / sum(count))  # Scale within each condition

# Ensure cell type order matches `cell_colors`
relative_frequencies$Cell_Type <- factor(
  relative_frequencies$Cell_Type,
  levels = names(cell_colors)  # Ensure factor levels match color names
)

# Plot scaled frequencies using `cell_colors`
p <- ggplot(relative_frequencies, aes(x = Condition, y = freq_scaled, fill = Cell_Type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cell_colors) +  # Use UMAP colors
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +  # Show raw scale (0-1)
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15),
    legend.title = element_blank()
  ) +
  labs(
    x = "Condition",
    y = "Scaled Frequency (0-1)",
    title = "Cell Type Frequencies Scaled to 0-1"
  )

# Display plot
print(p)

# # Create a data frame with unique combinations of Condition and Cell_Type
# unique_pairs <- testtt %>% distinct()
# # Count occurrences for each pair of Condition and Cell_Type
# counts <- testtt %>% count(Condition, Cell_Type)
# # Merge unique pairs with counts
# result <- merge(unique_pairs, counts, by = c("Condition", "Cell_Type"), all.x = TRUE)
# p <- ggbarstats(
#   data = result,
#   x = Condition,
#   y = Cell_Type,
#   counts = n,
#   bf.message=F, label="both"
# )
# p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Counts per condition
global.prop = as.data.frame( prop.table(table(Idents(Immune_norm))) )
global.prop
global.counts = as.data.frame.matrix( table(Idents(Immune_norm), Immune_norm$orig.ident) )
global.counts



# library(rstatix)
# xyz = prop_test(global.counts, detailed = TRUE)
# prop_test(t(global.counts), detailed = TRUE)
# 
# prop.by.origin = prop.table(table(Idents(Immune_norm), Immune_norm$orig.ident), margin = 2)
# prop.by.origin = as.data.frame.matrix(prop.by.origin)
# prop.by.origin
# 
# library(ggplot2)
# library(dplyr)
# 
# # Compute chi-squared test p-values for each row
# p_values <- apply(global.counts, 1, function(x) {
#   chisq.test(x)$p.value
# })
# 
# # Create a dataframe with the p-values
# p_values_df <- data.frame(row.names = rownames(global.counts), p_values = p_values)
# 
# # Combine the counts and p-values
# counts_with_p <- cbind(global.counts, p_values_df)
# library(tibble)
# library(tidyr)
# # Melt the data for ggplot2
# counts_melted <- counts_with_p %>%
#   rownames_to_column(var = "Category") %>%
#   gather(key = "Group", value = "Count", -Category, -p_values)
# 
# # Create the barplot with ggplot2
# ggplot(counts_melted, aes(x = Category, y = Count, fill = Group)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   geom_text(aes(label = cut(p_values, breaks = c(0, 0.001, 0.01, 0.05, 1),
#                             labels = c("***", "**", "*", ""))),
#             vjust = -0.5, size = 6, y=-200) +
#   labs(title = "Barplot of Counts by Category",
#        x = "Cell Type",
#        y = "Cell Count",
#        fill = "Group") +
#   theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#########################
### Cell-cycle score ####
#########################

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

library(stringr)
s.genes <- str_to_title(s.genes)
g2m.genes <- str_to_title(g2m.genes)

Immune_norm_Phase <- CellCycleScoring(Immune_norm, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

RidgePlot(Immune_norm_Phase, features = c("Mcm5", "Rad51", "Cdk1", "Top2a"), ncol = 2)

# Immune_norm_Phase <- ScaleData(Immune_norm_Phase, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Immune_norm_Phase), block.size = 2000)

DimPlot(Immune_norm_Phase, reduction = "umap", label = TRUE, pt.size = 0.5, cols = viridis::magma(3, begin = 0.15, end=0.85), shuffle = T, seed = seed, label.box = T, repel = T, label.size = 5, label.color = 'white', alpha = 0.75, group.by = 'Phase')

DimPlot(Immune_norm_Phase, reduction = "umap", label = TRUE, pt.size = 0.3, cols = viridis::turbo(length(levels(Immune_norm)), begin = 0.1, end=0.9), shuffle = T, seed = seed, label.box = F, repel = T, label.size = 4, label.color = 'black', alpha = 0.75, group.by = 'old.ident', split.by = 'Phase') + NoLegend()

DimPlot(Immune_norm, reduction = "umap", label = TRUE, pt.size = 0.4, shuffle = T, seed = 44, label.box = T, repel = T, label.size = 4, label.color = 'black', alpha = 0.75)+ NoLegend()
# Immune_norm <- RunPCA(Immune_norm, features = VariableFeatures(Immune_norm), nfeatures.print = 10)
Immune_norm$Phase = Immune_norm_Phase$Phase
rm(Immune_norm_Phase)

##########################
#### Diff. Expression ####
##########################
## Igkv, Ighv
# inverse order WT / AD ?

# Load library for the volcano plot and ktest
library(EnhancedVolcano)
# library(ktest)
# # create dedicated Python virtual environment
# reticulate::virtualenv_create("ktest")
# # activate the python environment
# reticulate::use_virtualenv(virtualenv = "ktest", required = TRUE)
# # verify python version
# reticulate::py_config()
# Function for Wilcoxon Test
DE_Wilcox <- function(object, ident1, ident2=NULL){
  markers = FindMarkers(object, ident.1 = ident1, ident.2 = ident2,
                        only.pos = F, min.pct = .025, test.use = 'wilcox', densify=T)
  return (markers)
}
# Function for MAST Test
DE_MAST <- function(object, ident1, ident2=NULL){
  markers = FindMarkers(object, ident.1 = ident1, ident.2 = ident2,
                        only.pos = F, min.pct = .025, test.use = 'MAST', densify=T)
  return (markers)
}
DE_DESeq2 <- function(object, ident1, ident2=NULL){
  markers = FindMarkers(object, ident.1 = ident1, ident.2 = ident2,
                        only.pos = F, min.pct = .025, test.use = 'DESeq2', densify=T)
  return (markers)
}
# Function for ktest
ktest <- function(data_tab, metadata_tab, comp1, comp2, markers){
  new_markers = markers[metadata_tab %in% c(comp1,comp2),]
  data_tab <- t(as.matrix(data_tab))
  metadata_tab <- as.data.frame(metadata_tab)
  # create Ktest object
  kt_1 = ktest_init(
    data = data_tab, metadata = metadata_tab, 
    sample_names = c(comp1, comp2)
  )
  # run test
  test(kt_1)
  # p values
  pvalues = get_pvalues(kt_1)
  padj = p.adjust(pvalues, method="BH")
  new_markers$p_val = as.numeric(pvalues)
  new_markers$p_val_adj = as.numeric(padj)
  return(new_markers)
  # get_statistics(kt_1, stat = 'kfda', contrib = FALSE, t_max = 50)
  # get_pvalues(kt_1, stat = 'kfda', permutation = FALSE, t_max = 50)
  # get_proj(kt_1)
  # proj <- get_proj(kt_1, contrib = FALSE, t_max = 50)
  # names(proj)
  # as_tibble(proj[[1]])
  # proj_contrib <- get_proj(kt_1, contrib = TRUE, t_max = 50)
  # names(proj_contrib)
  # as_tibble(proj_contrib[[1]])
  # as_tibble(proj_contrib[[2]])
  # 
  # plt <- reticulate::import("matplotlib.pyplot")
  # fig <- kt_1$plot_density(t = 1L)
  # fig[[1]]
  # fig[[2]]
  # plt$show()
  # fig <- kt_1$scatter_projection(t_x = 1L, t_y = 2L, proj_xy = c('kfda_contrib', 'kfda_contrib'))
  # fig[[1]]
  # fig[[2]]
  # plt$show()
}
# Function for creating and saving volcano plots
Volcano <- function(name, markers, title, comparison, strict=FALSE){
  if (strict == TRUE){logFC = 2}else{logFC = .5}
  pdf(paste0('/Volumes/Biclab/Joseph/2024_skull_bone_marrow/DE_genes_per_cell_type/Plots/', name, '_Volcano.pdf'),
      width = 12, height=9)
  print(EnhancedVolcano(markers,
                        rownames(markers),
                        x ="avg_log2FC",  arrowheads = F, max.overlaps = 30,
                        y ="p_val_adj", FCcutoff = logFC, title = title, boxedLabels = T, parseLabels = F, drawConnectors = T,
                        subtitle = comparison,
                        xlim = c(min(markers$avg_log2FC) - .5, max(markers$avg_log2FC) + .5) ))
  dev.off()
}
# Function to save the DE gene lists
Save_DE_genes <- function(dataframe, filename){
  write.table(dataframe,
              file = paste0('/Volumes/Biclab/Joseph/2024_skull_bone_marrow/DE_genes_per_cell_type/DFs/',filename,'.csv'),
              sep = ',', quote = F, row.names = T, col.names = NA, eol = '\r\n')
  
}

tmp = Immune_norm
tmp$origin_cell_type <- paste(Immune_norm$orig.ident, Immune_norm@active.ident, sep = "_")
Idents(tmp) <- "origin_cell_type"
for (i in levels(Immune_norm@active.ident)){
  print(i)
  tryCatch({
    i_2 = gsub("/", "-", i)
    # Find Markers of the cell type in question Vs all the other cells in the dataset
    d.markers = DE_Wilcox(Immune_norm, i)
    # Find markers between the same cell type of WT vs AD
    d.markers2 = DE_Wilcox(tmp, paste('AD Skull', i, sep = "_"), paste('WT Skull', i, sep = "_"))
    # Same with MAST
    mast_d.markers = DE_MAST(Immune_norm, i)
    mast_d.markers2 = DE_MAST(tmp, paste('AD Skull', i, sep = "_"), paste('WT Skull', i, sep = "_"))
    # # Same with DESeq2
    # deseq2_d.markers = DE_DESeq2(Immune_norm, i, strict=T)
    # deseq2_d.markers2 = DE_DESeq2(tmp, paste('AD Skull', i, sep = "_"), paste('WT Skull', i, sep = "_"))
    # Save the dataframes with DE info 
    Save_DE_genes(d.markers, paste0('Wilcoxon/',i_2,'_vsALL'))
    Save_DE_genes(d.markers2, paste0('Wilcoxon/',i_2,'_ADvsWT'))
    Save_DE_genes(mast_d.markers, paste0('MAST/',i_2,'_vsALL'))
    Save_DE_genes(mast_d.markers2, paste0('MAST/',i_2,'_ADvsWT'))
    # Save_DE_genes(deseq2_d.markers, paste0('DESeq2/',i_2,'_vsALL'))
    # Save_DE_genes(deseq2_d.markers2, paste0('DESeq2/',i_2,'_ADvsWT'))
    # Create and save the volcano plots
    Volcano(paste0('Wilcoxon/',i_2,'_vs_ALL'), d.markers, i_2, 'vs ALL - Wilcoxon', strict=T)
    Volcano(paste0('Wilcoxon/',i_2,'_ADvsWT'), d.markers2, i_2, 'AD vs WT - Wilcoxon')
    Volcano(paste0('MAST/',i_2,'_vs_ALL'), mast_d.markers, i_2, 'vs ALL - MAST', strict=T)
    Volcano(paste0('MAST/',i_2,'_ADvsWT'), mast_d.markers2, i_2, 'AD vs WT - MAST')
    # Volcano(paste0('DESeq2/',i_2,' vs ALL'), deseq2_d.markers, i_2, 'vs ALL - DESeq2')
    # Volcano(paste0('DESeq2/',i_2,'ADvsWT'), deseq2_d.markers2, i_2, 'AD vs WT - DESeq2')
    # AD vs WT with ktest
    # k.markers = ktest(tmp@assays$SCT$data, tmp$origin_cell_type, paste('AD Skull', i, sep = "_"), paste('WT Skull', i, sep = "_"), as.data.frame(as.numeric(d.markers2$avg_log2FC), row.names = rownames(d.markers2)))
    # Save_DE_genes(k.markers, paste0('ktest/','ADvsWT_',i_2))
    # Volcano(paste0('ktest/',i_2,' vs ALL'), k.markers, i_2, 'AD vs WT - Wilcoxon')
  })
}
## Now all cells of WT vs AD
tmp = Immune_norm
Idents(tmp) <- "orig.ident"
d.markers = DE_Wilcox(tmp, 'AD Skull', 'WT Skull')
MAST_d.markers = DE_MAST(tmp, 'AD Skull', 'WT Skull')
Save_DE_genes(d.markers, paste0('Wilcoxon/WTvsAD_ALL'))
Save_DE_genes(MAST_d.markers, paste0('MAST/WTvsAD_ALL'))
Volcano(paste0('Wilcoxon/ALL_ADvsWT'), d.markers, 'All cell types', 'AD vs WT - Wilcoxon')
Volcano(paste0('MAST/ALL_ADvsWT'), MAST_d.markers, 'All cell types', 'AD vs WT - MAST')

# Now for each cluster
tmp = Immune_norm
Idents(tmp) <- "seurat_clusters"
for (i in levels(tmp$seurat_clusters)){
  print(i)
  tryCatch({
    i_2 = gsub("/", "-", i)
    # Same with MAST
    mast_d.markers = DE_MAST(tmp, i)
    # mast_d.markers2 = DE_MAST(tmp, paste('AD Skull', i, sep = "_"), paste('WT Skull', i, sep = "_"))
    Save_DE_genes(mast_d.markers, paste0('MAST/Clusters/',i_2,'_vsALL'))
    # Save_DE_genes(mast_d.markers2, paste0('MAST/Clusters/',i_2,'_ADvsWT'))
    Volcano(paste0('MAST/Clusters/',i_2,'_vs_ALL'), mast_d.markers, i_2, 'vs ALL - MAST', strict=T)
    # Volcano(paste0('MAST/Clusters/',i_2,'_ADvsWT'), mast_d.markers2, i_2, 'AD vs WT - MAST')
  })
}


#### Here, only G1 cells #####
for (i in levels(Immune_norm@active.ident)){
  print(i)
  tryCatch({
    i_2 = gsub("/", "-", i)
    tmp = subset(Immune_norm, subset = Phase != 'S')
    # Find Markers of the cell type in question Vs all the other cells in the dataset
    # d.markers = DE_Wilcox(tmp, i)
    mast_d.markers = DE_MAST(tmp, i)
    # Find markers between the same cell type of WT vs AD
    tmp$origin_cell_type <- paste(Immune_norm$orig.ident, Immune_norm@active.ident, sep = "_")
    Idents(tmp) <- "origin_cell_type"
    # d.markers2 = DE_Wilcox(tmp, paste('AD Skull', i, sep = "_"), paste('WT Skull', i, sep = "_"))
    mast_d.markers2 = DE_MAST(tmp, paste('AD Skull', i, sep = "_"), paste('WT Skull', i, sep = "_"))
    # Save tables
    # Save_DE_genes(d.markers, paste0('Wilcoxon/G1_only/',i_2))
    # Save_DE_genes(d.markers2, paste0('Wilcoxon/G1_only/ADvsWT_',i_2))
    Save_DE_genes(mast_d.markers, paste0('MAST/G1_only/',i_2))
    Save_DE_genes(mast_d.markers2, paste0('MAST/G1_only/ADvsWT_',i_2))
    # Save Volcanos
    # Volcano(paste0('Wilcoxon/G1_only/',i_2,'_vs_ALL'), d.markers, i_2, 'vs ALL - Wilcoxon', strict=T)
    # Volcano(paste0('Wilcoxon/G1_only/',i_2,'_ADvsWT'), d.markers2, i_2, 'AD vs WT - Wilcoxon')
    Volcano(paste0('MAST/G1_only/',i_2,'_vs_ALL'), mast_d.markers, i_2, 'vs ALL - MAST', strict=T)
    Volcano(paste0('MAST/G1_only/',i_2,'_ADvsWT'), mast_d.markers2, i_2, 'AD vs WT - MAST')
  })
}
## G1 Now all cells of WT vs AD
tmp = subset(Immune_norm, subset = Phase != 'S')
Idents(tmp) <- "orig.ident"
d.markers = DE_MAST(tmp, 'AD Skull', 'WT Skull')
Save_DE_genes(d.markers, paste0('MAST/G1_only/WTvsAD_ALL'))
Volcano(paste0('MAST/G1_only/ALL_ADvsWT'), d.markers, 'No S phase cells', 'AD vs WT - MAST')
rm(tmp)

###########################
### PHATE TRAJECTORIES ####
###########################
library(phateR)
library(viridis)
library(Rmagic)

Immune_norm$cell_type = Immune_norm@active.ident

T_WT = subset(Immune_norm, subset = orig.ident == 'WT Skull' & cell_type %in% c('T-like', 'CD4 T cells', 'CD8 T cells', 'NKT'))
T_AD = subset(Immune_norm, subset = orig.ident == 'AD Skull' & cell_type %in% c('T-like', 'CD4 T cells', 'CD8 T cells', 'NKT'))
B_WT = subset(Immune_norm, subset = orig.ident == 'WT Skull' & cell_type %in% c('B-like', 'Pro-B cells', 'Pre-B cells', 'Immature B cells', 'Mature B cells', 'Pax5 B cells'))
B_AD = subset(Immune_norm, subset = orig.ident == 'AD Skull' & cell_type %in% c('B-like', 'Pro-B cells', 'Pre-B cells', 'Immature B cells', 'Mature B cells', 'Pax5 B cells'))

matrix = t(as.data.frame(GetAssayData(object = Immune_norm, assay = "RNA", layer = "counts")))

# Filtering data
## keep genes expressed in at least 10 cells
keep_cols <- colSums(matrix > 0) > 10
matrix <- matrix[,keep_cols]

## look at the distribution of library sizes
ggplot() +
  geom_histogram(aes(x=rowSums(matrix)), bins=50) +
  geom_vline(xintercept = 1000, color='red')

# keep cells with at least 1000 UMIs
keep_rows <- rowSums(matrix) > 1000
matrix <- matrix[keep_rows,]

# Normalise
matrix <- library.size.normalize(matrix)
matrix <- sqrt(matrix)
matrix = as.data.frame(matrix)

#subset
T_WT_matrix = na.omit(matrix[names(T_WT$orig.ident),])
T_AD_matrix = na.omit(matrix[names(T_AD$orig.ident),])
B_WT_matrix = na.omit(matrix[names(B_WT$orig.ident),])
B_AD_matrix = na.omit(matrix[names(B_AD$orig.ident),])

sample_numbers = sample(1:length(matrix), 1000)
# PCA
matrix_PCA <- as.data.frame(prcomp(matrix[sample_numbers,])$x)
ggplot(matrix_PCA) +
  geom_point(aes(PC1, PC2, color=matrix$Ngp[sample_numbers])) +
  labs(color="Ngp") +
  scale_color_viridis(option="B")

# Running PHATE
data_phate <- phate(matrix, n.jobs = 13, seed = seed, npca = 35, t = 100)
data_phate_TWT <- phate(T_WT_matrix, n.jobs = 13, seed = seed, npca = 35, t = 100)
data_phate_TAD <- phate(T_AD_matrix, n.jobs = 13, seed = seed, npca = 35, t = 100)
data_phate_BWT <- phate(B_WT_matrix, n.jobs = 13, seed = seed, npca = 35, t = 100)
data_phate_BAD <- phate(B_AD_matrix, n.jobs = 13, seed = seed, npca = 35, t = 100)

# use gene marker to visualise
ggplot(data_phate) +
  geom_point(aes(PHATE1, PHATE2, color=c(matrix$Ngp))) +
  labs(color="Ngp") +
  scale_color_viridis(option="B")
ggplot(data_phate) +
  geom_point(aes(PHATE1, PHATE2, color=c(matrix$Cd3g))) +
  labs(color="Cd3g") +
  scale_color_viridis(option="B")
ggplot(data_phate) +
  geom_point(aes(PHATE1, PHATE2, color=c(matrix$Cd19))) +
  labs(color="Cd19") +
  scale_color_viridis(option="B")

ggplot(data_phate) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=Immune_norm@active.ident[rownames(matrix)] ) ) + 
  scale_color_manual(values = cell_colors) + theme_classic()
ggplot(data_phate_TWT) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=T_WT@active.ident[rownames(T_WT_matrix)] ) ) + 
  scale_color_manual(values = cell_colors) + theme_classic()
ggplot(data_phate_TAD) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=T_AD@active.ident[rownames(T_AD_matrix)] ) ) + 
  scale_color_manual(values = cell_colors) + theme_classic()

# # Rerunning PHATE with new parameters
# data_phate <- phate(matrix, n.jobs = 13, knn=7, decay=100, init=data_phate)
# ggplot(data_phate) +
#   geom_point(aes(PHATE1, PHATE2, color=matrix$Ngp)) +
#   labs(color="Ngp") +
#   scale_color_viridis(option="B")

# Visualizing imputed genes on PHATE with MAGIC
MAGIC <- magic(matrix, genes='all_genes', seed = seed, n.jobs = 13, npca = 35, t = 100)

ggplot(MAGIC) +
  geom_point(aes(Cd3d, Cd3g, color=Cd3e)) +
  scale_color_viridis(option="B")
ggplot(MAGIC) +
  geom_point(aes(Ccl6, Ly6g, color=Ngp)) +
  scale_color_viridis(option="B")

# ggplot(Immune_norm@reductions$pca) +
#   geom_point(aes(x=PC_1, y=PC_2, color=MAGIC$result$Cd3e)) +
#   scale_color_viridis(option="B") +
#   labs(color="Cd3e")

# Phate + Magic
ggplot(data_phate) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=MAGIC$result$Ngp)) +
  scale_color_viridis(option="B") +
  labs(color="Ngp")+ ylim(-.02,0.03)
ggplot(data_phate) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=MAGIC$result$Cd3e)) +
  scale_color_viridis(option="B") +
  labs(color="Cd3e")+ ylim(-.02,0.03)
ggplot(data_phate) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=MAGIC$result$Cd19)) +
  scale_color_viridis(option="B") +
  labs(color="Cd19")+ ylim(-.02,0.03)

ggplot(data_phate) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=MAGIC$result$Ly6g)) +
  scale_color_viridis(option="B") +
  labs(color="Ly6g")+ ylim(-.02,0.03)


ggplot(data_phate) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=Immune_norm@active.ident[rownames(matrix)] ) ) + 
  scale_color_manual(values = cell_colors) + theme_classic() + ylim(-.02,0.03)

ggplot(data_phate, aes(x=PHATE1, y=PHATE2, color=Immune_norm$orig.ident[rownames(matrix)])) +
  geom_point() + ylim(-.02,0.03)


# For cell-types
ggplot(data_phate_TWT) +
  geom_point(aes(PHATE1, PHATE2, color=c(MAGIC$result[rownames(data_phate_TWT$embedding),]$Cd3d))) +
  labs(color="Cd3d") +
  scale_color_viridis(option="B")
ggplot(data_phate_TAD) +
  geom_point(aes(PHATE1, PHATE2, color=c(MAGIC$result[rownames(data_phate_TAD$embedding),]$Cd3d))) +
  labs(color="Cd3d") +
  scale_color_viridis(option="B")

ggplot(data_phate_BWT) +
  geom_point(aes(PHATE1, PHATE2, color=c(MAGIC$result[rownames(data_phate_BWT$embedding),]$Cd19))) +
  labs(color="Cd19") +
  scale_color_viridis(option="B")
ggplot(data_phate_BAD) +
  geom_point(aes(PHATE1, PHATE2, color=c(MAGIC$result[rownames(data_phate_BAD$embedding),]$Cd19))) +
  labs(color="Cd19") +
  scale_color_viridis(option="B")

ggplot(data_phate_TAD) +
  geom_point(aes(PHATE1, PHATE2, color=c(MAGIC$result[rownames(data_phate_TAD$embedding),]$Cd44))) +
  labs(color="Cd3d") +
  scale_color_viridis(option="B")
ggplot(data_phate_BWT) +
  geom_point(aes(PHATE1, PHATE2, color=c(MAGIC$result[rownames(data_phate_BWT$embedding),]$Bank1))) +
  labs(color="Cd19") +
  scale_color_viridis(option="B")


ggplot(data_phate_TWT) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=Immune_norm@active.ident[rownames(data_phate_TWT$embedding)] ) ) + 
  scale_color_manual(values = cell_colors) + theme_classic()
ggplot(data_phate_TAD) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=Immune_norm@active.ident[rownames(data_phate_TAD$embedding)] ) ) + 
  scale_color_manual(values = cell_colors) + theme_classic()

ggplot(data_phate_BWT) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=Immune_norm@active.ident[rownames(data_phate_BWT$embedding)] ) ) + 
  scale_color_manual(values = cell_colors) + theme_classic()
ggplot(data_phate_BAD) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=Immune_norm@active.ident[rownames(data_phate_BAD$embedding)] ) ) + 
  scale_color_manual(values = cell_colors) + theme_classic()


#######################
#### GENE ONTOLOGY ####
#######################
library(gprofiler2)

GO <- function(x){
  gostres = gost(query = rownames(x),
                 organism = "mmusculus", ordered_query = T,
                 exclude_iea = T, significant = F,
                 user_threshold = 0.05, correction_method = "g_SCS",
                 domain_scope = "annotated",
                 numeric_ns = "", sources = NULL)
  gostres$result = gostres$result[gostres$result$significant == T,]
  return(gostres)
}

for (i in levels(Immune_norm@active.ident)){
  print(i)
  i = gsub("/", "-", i)
  tryCatch({
  x = read.csv(paste0('/Volumes/Biclab/Joseph/2024_skull_bone_marrow/DE_genes_per_cell_type/DFs/WTvsAD_',i,'.csv'), sep = ';' )
  x = x[x$p_val_adj<1e-05,]
  x = x[abs(x$avg_log2FC)>.5,]
  gostres = GO(x)
  if (length(gostres$result$p_value) >1){
    head(gostres$result, 10)
    gostplot(gostres, capped = TRUE, interactive = TRUE)
    p = gostplot(gostres, capped = TRUE, interactive = F)
    pp = publish_gostplot(p, highlight_terms = unlist(head(gostres$result[order(gostres$result$p_value),], 20)$term_id), filename = NULL) +
      labs(title = i)
    ggsave(
      paste0('/Volumes/Biclab/Joseph/2024_skull_bone_marrow/DE_genes_per_cell_type/Plots/GO/WTvsAD_',i,'_GO.pdf'),
      plot = pp,
      width = 10,
      height = 4+(length(gostres$result$p_value)*.25),
      dpi = 250,
      create.dir = T
    )
  }
  })
}

## All 
x = read.csv('/Volumes/Biclab/Joseph/2024_skull_bone_marrow/DE_genes_per_cell_type/DFs/WTvsAD_ALL.csv', sep = ';' )
x = x[x$p_val_adj<1e-05,]
x = x[abs(x$avg_log2FC)>.5,]
gostres = GO(x)
head(gostres$result, 10)
gostplot(gostres, capped = TRUE, interactive = TRUE)
p = gostplot(gostres, capped = TRUE, interactive = F)
pp = publish_gostplot(p, highlight_terms = unlist((gostres$result[order(gostres$result$p_value),])$term_id), filename = NULL) +
  labs(title = 'All cell-types')
ggsave(
  paste0('/Volumes/Biclab/Joseph/2024_skull_bone_marrow/DE_genes_per_cell_type/Plots/GO/WTvsAD_ALL_GO.pdf'),
  plot = pp,
  width = 13,
  height = 5+(length(unlist((gostres$result[order(gostres$result$p_value),])$term_id))*.25),
  dpi = 250,
  create.dir = T
)
# AD up
x = read.csv('/Volumes/Biclab/Joseph/2024_skull_bone_marrow/DE_genes_per_cell_type/DFs/WTvsAD_ALL.csv', sep = ';' )
x = x[x$p_val_adj<1e-05,]
x = x[(x$avg_log2FC)>.5,]
gostres = GO(x)
head(gostres$result, 10)
gostplot(gostres, capped = TRUE, interactive = TRUE)
p = gostplot(gostres, capped = TRUE, interactive = F)
pp = publish_gostplot(p, highlight_terms = unlist((gostres$result[order(gostres$result$p_value),])$term_id), filename = NULL) +
  labs(title = 'All cell-types: AD upregulated')
ggsave(
  paste0('/Volumes/Biclab/Joseph/2024_skull_bone_marrow/DE_genes_per_cell_type/Plots/GO/AD_up_ALL_GO.pdf'),
  plot = pp,
  width = 13,
  height = 5+(length(unlist((gostres$result[order(gostres$result$p_value),])$term_id))*.25),
  dpi = 250,
  create.dir = T
)
# AD down
x = read.csv('/Volumes/Biclab/Joseph/2024_skull_bone_marrow/DE_genes_per_cell_type/DFs/WTvsAD_ALL.csv', sep = ';' )
x = x[x$p_val_adj<1e-05,]
x = x[(x$avg_log2FC)<.5,]
gostres = GO(x)
head(gostres$result, 10)
gostplot(gostres, capped = TRUE, interactive = TRUE)
p = gostplot(gostres, capped = TRUE, interactive = F)
pp = publish_gostplot(p, highlight_terms = unlist(head(gostres$result[order(gostres$result$p_value),], 25)$term_id), filename = NULL) +
  labs(title = 'All cell-types: AD downregulated')
ggsave(
  paste0('/Volumes/Biclab/Joseph/2024_skull_bone_marrow/DE_genes_per_cell_type/Plots/GO/AD_down_ALL_GO.pdf'),
  plot = pp,
  width = 13,
  height = 5+(length(unlist(head(gostres$result[order(gostres$result$p_value),], 25)$term_id))*.25),
  dpi = 250,
  create.dir = T
)
