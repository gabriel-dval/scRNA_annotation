# Samir ALI MOUSSA & Joseph JOSEPHIDES
# Institut Pasteur, Paris FR
# 07 Mar 2024

# Load packages
library(Seurat)
library(SeuratObject)
library(BiocParallel)
library(scDblFinder)
library(Matrix)
library(ggplot2)
library(tidyverse)
library(harmony)

# Set Random Seed
seed = 1111111 # this does not set the seed. It will attribute a number to a variable called seed.
set.seed(seed)

# Multiprocessing parameters
# Define number of cores
cores=6
library(future)
plan("multisession", workers = 2)
handlers(global = TRUE)

## Load Files
Ad_CP_R2 = Read10X("/Volumes/Biclab/Joseph/Development_project/CiteSeq/SoupX_matrix_Ad_CP_Immune_R2/")
Ad_CP_R2_seurat = CreateSeuratObject(counts = Ad_CP_R2, project = "Ad_R2")
Ad_Ms.R2=CreateAssay5Object(counts = Read10X('/Volumes/BICommunication/processed_data/samir/240413_VH00505_236_AAFK3HCM5__Sun_Apr_14_06h27m34_2024_CiteSeq_2024/counts_240413_VH00505_236_AAFK3HCM5/Ad3_CP_Immune/outs/per_sample_outs/Ad3_CP_Immune/count/sample_filtered_feature_bc_matrix')[["Antibody Capture"]])
Ad_CP_R2_seurat[["ADT"]] = Ad_Ms.R2

Ad_CP = Read10X("/Volumes/Biclab/Joseph/Development_project/CiteSeq/SoupX_matrix_Ad_CP_Immune/")
Ad_CP_seurat = CreateSeuratObject(counts = Ad_CP, project = "Ad")
Ad_ADT=CreateAssay5Object(counts = Read10X("/Volumes/Biclab/Samir/single\ cell/Cite-seq/Cite-seq\ P0-Ad/Ad_CP_Immune/count/sample_filtered_feature_bc_matrix")[["Antibody Capture"]])
Ad_CP_seurat[["ADT"]] = Ad_ADT

P0_CP = Read10X("/Volumes/Biclab/Joseph/Development_project/CiteSeq/SoupX_matrix_P0_CP_Immune/")
P0_CP_seurat = CreateSeuratObject(counts = P0_CP, project = "P0")
P0_ADT=CreateAssay5Object(counts = Read10X("/Volumes/Biclab/Samir/single\ cell/Cite-seq/Cite-seq\ P0-Ad/P0_CP_Immune/count/sample_filtered_feature_bc_matrix")[["Antibody Capture"]])
P0_CP_seurat[["ADT"]] = P0_ADT

P6_CP =Read10X("/Volumes/Biclab/Joseph/Development_project/CiteSeq/SoupX_matrix_P6_CP_Immune/")
P6_CP_seurat = CreateSeuratObject(counts = P6_CP, project = "P6")
P6_ADT=CreateAssay5Object(counts = Read10X("/Volumes/Biclab/Samir/single\ cell/Cite-seq/Cite-seq\ P0-Ad/P6_CP_Immune/count/sample_filtered_feature_bc_matrix")[["Antibody Capture"]])
P6_CP_seurat[["ADT"]] = P6_ADT

P10_CP=Read10X("/Volumes/Biclab/Joseph/Development_project/CiteSeq/SoupX_matrix_P10_CP_Immune/")
P10_CP_seurat = CreateSeuratObject(counts = P10_CP, project = "P10")
P10_ADT=CreateAssay5Object(counts = Read10X("/Volumes/Biclab/Samir/single\ cell/Cite-seq/Cite-seq\ P0-Ad/P10_CP_Immune/count/sample_filtered_feature_bc_matrix")[["Antibody Capture"]])
P10_CP_seurat[["ADT"]] = P10_ADT

P15_CP=Read10X("/Volumes/Biclab/Joseph/Development_project/CiteSeq/SoupX_matrix_P15_CP_Immune/")
P15_CP_seurat = CreateSeuratObject(counts = P15_CP, project = "P15")
P15_ADT=CreateAssay5Object(counts = Read10X("/Volumes/Biclab/Samir/single\ cell/Cite-seq/Cite-seq\ P0-Ad/P15_CP_Immune/count/sample_filtered_feature_bc_matrix")[["Antibody Capture"]])
P15_CP_seurat[["ADT"]] = P15_ADT

P21_CP=Read10X("/Volumes/Biclab/Joseph/Development_project/CiteSeq/SoupX_matrix_P21_CP_Immune/")
P21_CP_seurat = CreateSeuratObject(counts = P21_CP, project = "P21")
P21_ADT=CreateAssay5Object(counts = Read10X("/Volumes/Biclab/Samir/single\ cell/Cite-seq/Cite-seq\ P0-Ad/P21_CP_Immune/count/sample_filtered_feature_bc_matrix")[["Antibody Capture"]])
P21_CP_seurat[["ADT"]] = P21_ADT

###################################
### Individual analyses per age ###
###################################
List_seurat_ages = c(P0_CP_seurat, P6_CP_seurat, P10_CP_seurat, P15_CP_seurat, P21_CP_seurat, Ad_CP_seurat, Ad_CP_R2_seurat)

for (age in 1:length(List_seurat_ages)){
  print(age)
  List_seurat_ages[[age]]$percent.mt = PercentageFeatureSet(List_seurat_ages[[age]], pattern='^mt-')
}

for (age in 1:length(List_seurat_ages)){
  print(age)
  p1 = VlnPlot(List_seurat_ages[[age]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  p2 = FeatureScatter(List_seurat_ages[[age]], feature1 = "nCount_RNA", feature2 = "percent.mt") +
    theme_linedraw() 
  p3 = FeatureScatter(List_seurat_ages[[age]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_smooth(method = 'lm') + theme_linedraw() 
  print(p1+p2+p3 )
  readline(prompt="Press [enter] to continue")
}

List_seurat_ages[[1]] = subset(List_seurat_ages[[1]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5 & nCount_RNA < 20000) 
List_seurat_ages[[2]] = subset(List_seurat_ages[[2]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5 & nCount_RNA < 21000) 
List_seurat_ages[[3]] = subset(List_seurat_ages[[3]], subset = nFeature_RNA > 200 & nFeature_RNA < 4800 & percent.mt < 5 & nCount_RNA < 22000) 
List_seurat_ages[[4]] = subset(List_seurat_ages[[4]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5 & nCount_RNA < 22000) 
List_seurat_ages[[5]] = subset(List_seurat_ages[[5]], subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 5 & nCount_RNA < 18000) 
List_seurat_ages[[6]] = subset(List_seurat_ages[[6]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5 & nCount_RNA < 20000) 
List_seurat_ages[[7]] = subset(List_seurat_ages[[7]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5 & nCount_RNA < 20000) 

# QC metrics (after filtering)
for (age in 1:length(List_seurat_ages)){
  print(age)
  p1 = VlnPlot(List_seurat_ages[[age]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  p2 = FeatureScatter(List_seurat_ages[[age]], feature1 = "nCount_RNA", feature2 = "percent.mt") +
    theme_linedraw() 
  p3 = FeatureScatter(List_seurat_ages[[age]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_smooth(method = 'lm') + theme_linedraw() 
  print(p1+p2+p3 )
  readline(prompt="Press [enter] to continue")
}

# Merge the Seurat objects into one
CP_Immune = merge(List_seurat_ages[[1]], y = List_seurat_ages[2:length(List_seurat_ages)], project = 'CP_Immune')

# Join layers (collapses the individual datasets together and recreates the original counts and data layer)
CP_Immune = JoinLayers(CP_Immune)

# Add ages and batch
CP_Immune$Age = CP_Immune$orig.ident
CP_Immune$Age[CP_Immune$orig.ident =='Ad_R2'] = "Ad"
CP_Immune$Batch = "Batch1"
CP_Immune$Batch[CP_Immune$orig.ident =='Ad_R2'] =  "Batch2"

# Quality Control
## Violin plot
VlnPlot(CP_Immune, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter QC
plot1 = FeatureScatter(CP_Immune, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(CP_Immune, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
plot1 + plot2

# # Doublet exclusion
# ## Read in data as an sce object
sce = Seurat::as.SingleCellExperiment(CP_Immune, assay = 'RNA')
# ## Find Doublets
doublet.score = scDblFinder(sce, samples="orig.ident", returnType="table", BPPARAM=MulticoreParam(cores))
results = data.frame("Barcode" = rownames(doublet.score),
                     "scDblFinder_DropletType" = doublet.score$class,
                     "scDblFinder_Score" = doublet.score$score)

CP_Immune = AddMetaData(object = CP_Immune, metadata = doublet.score[colnames(CP_Immune),"score"], col.name = "doublet.score")
CP_Immune = AddMetaData(object = CP_Immune, metadata = doublet.score[colnames(CP_Immune),"class"], col.name = "doublet.class")

table(CP_Immune$doublet.class)
summary(CP_Immune$doublet.score)

# Remove cells based on QC criteria
CP_Immune = subset(CP_Immune, subset = doublet.class == 'singlet')

# Increase global max size
options(future.globals.maxSize= 5368709120) # ~1GB = 1024(Mb)*1024^2

# Normalise the data.
CP_Immune_norm = SCTransform(CP_Immune,
                             vars.to.regress = c("percent.mt"),
                             seed.use = seed, assay = 'RNA', variable.features.n = NULL,
                             variable.features.rv.th = 1.3)

# Run PCA
CP_Immune_norm = RunPCA(CP_Immune_norm, features = VariableFeatures(object = CP_Immune_norm), assay = "SCT")

# Run Harmony
CP_Immune_harmony = RunHarmony(CP_Immune_norm, assay.use="SCT", group.by.vars = c("Batch"))

# Examine PCA/harmony results
print(CP_Immune_harmony[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(CP_Immune_harmony, dims = 1:2, reduction = "pca")

print(CP_Immune_harmony[["harmony"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(CP_Immune_harmony, dims = 1:2, reduction = "harmony")

# Heatmap to see primary sources of heterogeneity
DimHeatmap(CP_Immune_harmony, dims = 1:2, cells = 500, balanced = TRUE, reduction = 'harmony')

# Determine the ‘dimensionality’ of the dataset
ElbowPlot(CP_Immune_harmony, ndims = 50, reduction = 'pca')
ElbowPlot(CP_Immune_harmony, ndims = 50, reduction = 'harmony')

# Cluster the cells at elbow point
CP_Immune_harmony = FindNeighbors(CP_Immune_harmony, dims = 1:30, reduction = 'harmony')
CP_Immune_harmony = FindClusters(CP_Immune_harmony, resolution = 2, algorithm = 4, random.seed = seed, group.singletons = F)

# Look at cluster IDs of the first 5 cells
head(Idents(CP_Immune_harmony), 5)

# Visualize PCA
DimPlot(CP_Immune_harmony, reduction = "pca", label = TRUE)
DimPlot(CP_Immune_harmony, reduction = "harmony", label = TRUE)
# DimPlot(CP_Immune_norm, reduction = "harmony", label = TRUE)

# non-linear dimensional reduction with UMAP
CP_Immune_harmony = RunUMAP(CP_Immune_harmony, dims = 1:30, seed.use = seed, reduction = 'harmony', assay = 'SCT')
# Visualise UMAP
order = c('P0','P6','P10','P15','P21','Ad','Ad_R2')
CP_Immune_harmony$orig.ident = factor(CP_Immune_harmony$orig.ident, levels = order)
order = c('P0','P6','P10','P15','P21','Ad')
CP_Immune_harmony$Age = factor(CP_Immune_harmony$Age, levels = order)
p = DimPlot(CP_Immune_harmony, reduction = "umap", label = TRUE, label.box = T, repel = T, shuffle = T) + NoLegend()
p
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/Dim_reductions/UMAP_clusters.pdf',
  plot = p,
  width = 10,
  height = 7,
  dpi = 300,
  create.dir = T
)
p = DimPlot(CP_Immune_harmony, reduction = "umap", label = FALSE, label.box = T, repel = T, shuffle = T,
            group.by  = 'Age', cols = viridis::magma(6, direction = -1, end = .9)) + 
  theme(legend.text=element_text(size=20))
p
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/Dim_reductions/UMAP_ages.pdf',
  plot = p,
  width = 10,
  height = 7,
  dpi = 300,
  create.dir = T
)
p = DimPlot(CP_Immune_harmony, reduction = "umap", label = TRUE, label.box = T, repel = T, shuffle = T,
            group.by = 'Batch', cols = viridis::cividis(2, direction = -1)) + NoLegend()
p
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/Dim_reductions/Plots/UMAP_batch.pdf',
  plot = p,
  width = 10,
  height = 7,
  dpi = 300,
  create.dir = T
)
DimPlot(CP_Immune_harmony, reduction = "umap", label = TRUE, label.box = T, repel = T, shuffle = T,
        group.by = 'Batch', split.by = 'Age') + NoLegend()
############################################

# Finding DE features (cluster biomarkers)
CP_Immune_norm.markers = FindAllMarkers(CP_Immune_harmony, only.pos = TRUE)
CP_Immune_norm.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# Print top gene marker for each cluster
for (i in c(1:length(levels(CP_Immune_harmony) ) ) ){
  Cluster_markers = CP_Immune_norm.markers %>% filter(cluster == i)  %>% filter(avg_log2FC > 1)
  print(c(paste0("Cluster",as.character(i),":"), paste0(Cluster_markers$gene[1:5]) ))
}

# reset back to immune_norm
CP_Immune_norm = CP_Immune_harmony

# Remove things we no longer need
# rm(P0_ADT, P0_CP, P0_CP_seurat,
#    P6_ADT, P6_CP, P6_CP_seurat,
#    P10_ADT, P10_CP, P10_CP_seurat,
#    P15_ADT, P15_CP, P15_CP_seurat,
#    P21_ADT, P21_CP, P21_CP_seurat,
#    Ad_ADT, Ad_CP, Ad_CP_seurat,
#    Ad_Ms.R2, Ad_CP_R2, Ad_CP_R2_seurat,
#    List_seurat_ages, CP_Immune, sce, results, plot1, plot2, doublet.score,
#    CP_Immune_norm.markers, CP_Immune_harmony)
# save.image(file='/Volumes/Biclab/Joseph/Development_project/CiteSeq/Code_shared_with_Samir/CiteSeq_2024.10.09.RData')

####################################################################
########### Samir, start from here #################################
####################################################################
load('/Volumes/Biclab/Joseph/Development_project/CiteSeq/Code_shared_with_Samir/CiteSeq_2024.10.09.RData')

# Plot some marker genes
# Mature B cells
FeaturePlot(CP_Immune_norm, features = c("Ms4a1", "Cd19"), blend = TRUE, order = T)
FeaturePlot(CP_Immune_norm, features = c("H2-Aa", "Cd19"), blend = TRUE, order = T)

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
                 "Clec9a", "Ifi205","Cd74",         # Dentritic Conventional
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
                 "Cdk2", "Rif1"                    # Replicating cells
)

DotPlot(subset(CP_Immune_norm, downsample = 300), features = marker_genes, split.by = 'seurat_clusters',
        cols = sample(viridis::turbo(length(levels(CP_Immune_norm$seurat_clusters))) )) + RotatedAxis() + NoLegend()

#Naming Cell types
Cell_type.ids = c("Macs 1",  #1
                  "Macs 1",  #2
                  "Macs 1",  #3
                  "Epiplexus",  #4
                  "Macs 1",  #5
                  "Macs 1",  #6
                  "Macs 1",  #7
                  "Epiplexus",  #8
                  "Macs 1", #9
                  "Macs 1", #10
                  "Macs 1", #11
                  "Epiplexus", #12
                  "Macs 2", #13
                  "Macs 1", #14
                  "Macs 2", #15
                  "Monocytes", #16
                  "Epithelium", #17
                  "Macs 1", #18
                  "B cells", #19
                  "Proliferating Macs", #20
                  "Epithelium", #21
                  "Neutrophils", #22
                  "Macs 1", #23
                  "Macs 1", #24
                  "Non B Lymphocytes", #25
                  "Macs 1", #26
                  "Non B Lymphocytes", #27
                  "Endothelium", #28
                  "B cells", #29
                  "cDC1", #30
                  "Macs 1", #31
                  "cDC2", #32
                  "Mesenchymal cells", #33
                  "Macs 1", #34
                  "Neutrophils", #35
                  "Epithelium", #36
                  "Macs 2", #37
                  "Proliferating Macs", #38
                  "Non B Lymphocytes", #39
                  "Epiplexus",  #40
                  "B cells",  #41
                  "Epithelium",  #42
                  "Macs 2",  #43
                  "Macs 1",  #44
                  "Epithelium",  #45
                  "Mesenchymal cells",  #46
                  "Mast cells",  #47
                  "Non B Lymphocytes", #48
                  "cDC2", #49
                  "Non B Lymphocytes", #50
                  "migDC", #51
                  "Macs 1", #52
                  "Basophils", #53
                  "Endothelium") #54

names(Cell_type.ids) = levels(CP_Immune_norm)
CP_Immune_norm = RenameIdents(CP_Immune_norm, Cell_type.ids)
CP_Immune_norm$cell_type = CP_Immune_norm@active.ident

# Save raw and normalised matrices for SCENIC
write.table(as.matrix(t(CP_Immune_norm@assays$SCT@counts)), file = '/Volumes/Biclab/Joseph/Development_project/CiteSeq/SCENIC/data/raw_matrix_2024.10.28.tsv',
            row.names = T, col.names = T, quote = F, sep = '\t')
write.table(as.matrix(t(CP_Immune_norm@assays$SCT@data)), file = '/Volumes/Biclab/Joseph/Development_project/CiteSeq/SCENIC/data/normalised_matrix_2024.10.28.tsv',
            row.names = T, col.names = T, quote = F, sep = '\t')
write.table(as.matrix(CP_Immune_norm@meta.data), file = '/Volumes/Biclab/Joseph/Development_project/CiteSeq/SCENIC/data/metadata_2024.10.28.tsv',
            row.names = T, col.names = T, quote = F, sep = '\t')

# subsetting to only have immune cells
Immune_only = subset(x = CP_Immune_norm, idents = c("Epithelium", "Endothelium", "Mesenchymal cells"), invert = TRUE)

#Define the colours
cell_types = sort(levels(CP_Immune_norm@active.ident))
cell_types

# Combine colors
color_palette = c("#7E1BC9",  # B cells - deep, distinct purple
                  "#89C2D6",  # Basophils - light, distinguishable blue
                  "#A63A3A",  # cDC1 - clear, warm red
                  "#FF6F61",  # cDC2 - bright coral red, distinct from cDC1
                  "#BBA0D7",  # Endothelium - soft lavender
                  "#2A6B44",  # Epiplexus - deep green, close to lymphoid but distinct
                  "#FFB74D",  # Epithelium - bright orange, higher contrast with Mesenchymal
                  "#66A773",  # Macs 1 - teal green, distinct yet cohesive with Macs 2
                  "#A6D89E",  # Macs 2 - light green for relationship with Macs 1
                  "#3B8DB4",  # Mast cells - deep blue, distinct from other blues
                  "#D95F0E",  # Mesenchymal cells - burnt orange, more contrast from Epithelium
                  "#8B4513",  # migDC - warm, earthy red to contrast with cDC1
                  "#C7E3CC",  # Monocytes - mint green, distinct among related immune cells
                  "#4F94B5",  # Neutrophils - mid blue, visually distinct from Mast cells
                  "#EAACC9",  # Non B Lymphocytes - muted pink for the lymphocyte group
                  "#4DBDB2"   # Proliferating Mac - clear teal, distinguishable among macrophages
)

# Create a named vector with cell types and corresponding colours
cell_colors = setNames(color_palette, cell_types)
immune_colors = cell_colors[names(cell_colors) %in% levels(Immune_only@active.ident)]

# UMAP of all cell types
p = DimPlot(CP_Immune_norm, reduction = "umap", label = TRUE, pt.size = 0.3,
            cols = cell_colors, shuffle = T, seed = seed, label.box = T, repel = T,
            label.size = 4, label.color = 'white', alpha = 0.75) + 
            theme(legend.position = "none")
p
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/Dim_reductions/UMAP_all_cells.pdf',
  plot = p,
  width = 8,
  height = 6,
  dpi = 300,
  create.dir = T
)

# UMAP of Immune cells only
p = DimPlot(Immune_only, reduction = "umap", label = TRUE, pt.size = 0.5,
            cols = immune_colors, shuffle = T, seed = seed, label.box = T, repel = T,
            label.size = 4, label.color = 'white', alpha = 0.75) + theme(legend.position = "none")
p
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/Dim_reductions/UMAP_immune_cells_subset.pdf',
  plot = p,
  width = 8,
  height = 6,
  dpi = 300,
  create.dir = T
)

# re-run UMAP on immune subset
Immune_only = RunUMAP(Immune_only, dims = 1:30, seed.use = seed, reduction = 'harmony', assay = 'SCT')
# visualuse new UMAP of Immune cells
p = DimPlot(Immune_only, reduction = "umap", label = TRUE, pt.size = 0.5,
            cols = immune_colors, shuffle = T, seed = seed, label.box = T, repel = T,
            label.size = 4, label.color = 'white', alpha = 0.75) + theme(legend.position = "none")
p
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/Dim_reductions/UMAP_immune_cells.pdf',
  plot = p,
  width = 8,
  height = 6,
  dpi = 300,
  create.dir = T
)

##################################
###### Barplot with counts #######
##################################
testtt = data.frame(Immune_only$Age, Immune_only@active.ident )
colnames(testtt) = c('Condition','Cell_Type')
testtt$Cell_Type = factor(testtt$Cell_Type, levels = sort(levels(testtt$Cell_Type)))

# Downsample
original_counts = table(testtt$Condition)
min_count = min(original_counts)
# Downsampling
downsampled_df = testtt %>%
  group_by(Condition) %>%
  slice_sample(n = min_count, replace = FALSE) %>%
  ungroup()

library(ggstatsplot)

p = ggbarstats(downsampled_df, x = Condition, y = Cell_Type,
               bf.message=F, label="both")
p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=15))
p
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/proportions/ggbarstats_downsampled_cell_type_composition_age.pdf',
  plot = p,
  width = 14,
  height = 6,
  dpi = 300,
  create.dir = T
)

#Plotting immune cell type proportions
prop = as.data.frame(Immune_only$Age)
colnames(prop) = c('Age')
prop$celltype = Immune_only@active.ident
head(prop)

prop$Age = factor(prop$Age, levels = c("P0", "P6", "P10", "P15", "P21", "Ad"))
prop$celltype = factor(prop$celltype, levels = sort(levels(Immune_only@active.ident)))

ggplot(prop, aes(x = Age, fill = celltype)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = immune_colors) +  # Assigning custom colors based on cell types
  labs(title = "Proportions of Immune Cells per Age",
       x = "Age",
       y = "Proportion",
       fill = "Cell Type") +
  theme_classic() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + theme(axis.line = element_blank(),legend.key.size = unit(1, 'cm'), legend.text=element_text(size=17))

# Subset to have an equal number of cells
cells = as.data.frame(Immune_only@meta.data)
x = min(table(cells$Age))
x
cells$barcodes = rownames(cells)
cells = cells %>% group_by(Age) %>% sample_n(size = x, replace = F)
cells = cells$barcodes

#plotting each age separately
p = DimPlot(subset(x = Immune_only, cells = cells), pt.size = 0.3, cols = immune_colors,
        shuffle = T, seed = seed, alpha = 0.75, ncol = 3, split.by = 'Age')
p
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/Dim_reductions/UMAP_immune_cells_downsampled_split_Age.pdf',
  plot = p,
  width = 12,
  height = 8,
  dpi = 300,
  create.dir = T
)

# Check mitochondrial percentage per cell type and cluster
p1 = VlnPlot(CP_Immune_norm, features = "percent.mt", split.by = "ident") + NoLegend()
p2 = VlnPlot(CP_Immune_norm, features = "percent.mt", group.by = "Age")+ NoLegend()
p1
p2
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/QC/VlnPlot_percent.mt_ident.pdf',
  plot = p1,
  width = 10,
  height = 7,
  dpi = 300,
  create.dir = T
)
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/QC/VlnPlot_percent.mt_age.pdf',
  plot = p2,
  width = 10,
  height = 7,
  dpi = 300,
  create.dir = T
)
# Check RNA counts
p1 = VlnPlot(CP_Immune_norm, features = "nCount_SCT", split.by = "ident") + NoLegend()
p2 = VlnPlot(CP_Immune_norm, features = "nCount_SCT", group.by = "Age")+ NoLegend()
p1
p2
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/QC/VlnPlot_nCount_SCT_ident.pdf',
  plot = p1,
  width = 10,
  height = 7,
  dpi = 300,
  create.dir = T
)
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/QC/VlnPlot_nCount_SCT_age.pdf',
  plot = p2,
  width = 10,
  height = 7,
  dpi = 300,
  create.dir = T
)
# Check number of genes (=features)
p1 = VlnPlot(CP_Immune_norm, features = "nFeature_SCT", split.by = "ident") + NoLegend()
p2 = VlnPlot(CP_Immune_norm, features = "nFeature_SCT", group.by = "Age")+ NoLegend()
p1
p2
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/QC/VlnPlot_nFeature_SCT_ident.pdf',
  plot = p1,
  width = 10,
  height = 7,
  dpi = 300,
  create.dir = T
)
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/QC/VlnPlot_nFeature_SCT_age.pdf',
  plot = p2,
  width = 10,
  height = 7,
  dpi = 300,
  create.dir = T
)

#########################
### Cell-cycle score ####
#########################
s.genes = cc.genes$s.genes
g2m.genes = cc.genes$g2m.genes

library(stringr)
s.genes = str_to_title(s.genes)
g2m.genes = str_to_title(g2m.genes)

Immune_norm_Phase = CellCycleScoring(CP_Immune_norm, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

p = DimPlot(Immune_norm_Phase, reduction = "umap", label = TRUE, pt.size = 0.3,
        cols = viridis::cividis(3, begin = 0.15, end=0.85), shuffle = T, seed = seed,
        label.box = T, repel = T, label.size = 5, label.color = 'white', alpha = 0.75,
        group.by = 'Phase')
p
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/Dim_reductions/UMAP_cell_cycle.pdf',
  plot = p,
  width = 12,
  height = 8,
  dpi = 300,
  create.dir = T
)

CP_Immune_norm$Phase = Immune_norm_Phase$Phase
rm(Immune_norm_Phase)

###############################################
############# SURFACE PROTEINS 💪 #############
###############################################
library(dsb)
CP_Immune_norm.prot = JoinLayers(CP_Immune_norm, assay = "ADT")
CP_Immune_norm.prot$cell_types = Idents(CP_Immune_norm.prot)
isotypes = rownames(CP_Immune_norm.prot$ADT$counts)[grep("^Iso", rownames(CP_Immune_norm.prot$ADT$counts))]
norm.adt = ModelNegativeADTnorm(cell_protein_matrix = CP_Immune_norm.prot$ADT$counts,
                                denoise.counts = TRUE,
                                use.isotype.control = TRUE,
                                isotype.control.name.vec = isotypes,
                                return.stats = F,
                                quantile.clipping = TRUE
)

# add dsb normalized matrix "cell.adt.dsb" to the "CITE" data (not counts!) slot
CP_Immune_norm.prot[["ADTnorm"]] = Seurat::CreateAssayObject(data = norm.adt)
DefaultAssay(CP_Immune_norm.prot) = "ADTnorm"

# Add RNA reductions
CP_Immune_norm.prot@reductions = CP_Immune_norm@reductions

### Run PCA, UMAP
# define proteins to use in clustering (non-isptype controls)
prots = rownames(CP_Immune_norm.prot@assays$ADTnorm@data)[!(rownames(CP_Immune_norm.prot@assays$ADTnorm@data) %in% isotypes)]
# Variable Features
VariableFeatures(CP_Immune_norm.prot) <- prots
# Run PCA
CP_Immune_norm.prot = CP_Immune_norm.prot %>%
  ScaleData() %>%
  RunPCA(assay = "ADTnorm", seed.use = seed, reduction.name = 'ADT_PCA')

# Determine the ‘dimensionality’ of the dataset
ElbowPlot(CP_Immune_norm.prot, ndims = 50, reduction = 'ADT_PCA')

# Cluster the cells at elbow point
CP_Immune_norm.prot = FindNeighbors(CP_Immune_norm.prot, dims = 1:17, reduction = 'ADT_PCA', assay = 'ADTnorm', features = prots, seed = seed)
CP_Immune_norm.prot = FindClusters(CP_Immune_norm.prot, resolution = 1.5, algorithm = 4, random.seed = seed,
                                           group.singletons = F, graph.name = 'ADTnorm_nn')
# Visualize PCA
p = DimPlot(CP_Immune_norm.prot, reduction = "ADT_PCA", label = TRUE,
            shuffle = T, seed = seed, label.box = T, repel = T,
            label.size = 4, label.color = 'white', alpha = 0.75)
p
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/Dim_reductions/PCA_ADT_clusters.pdf',
  plot = p,
  width = 12,
  height = 8,
  dpi = 300,
  create.dir = T
)

# non-linear dimensional reduction with UMAP
CP_Immune_norm.prot = RunUMAP(CP_Immune_norm.prot, dims = 1:17, seed.use = seed, reduction = 'ADT_PCA', assay = 'ADTnorm',
                              reduction.name = 'ADT_UMAP')
p = DimPlot(CP_Immune_norm.prot, reduction = "ADT_UMAP", label = TRUE,
            shuffle = T, seed = seed, label.box = T, repel = T,
            label.size = 4, label.color = 'white', alpha = 0.75)
p
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/Dim_reductions/UMAP_ADT_clusters.pdf',
  plot = p,
  width = 12,
  height = 8,
  dpi = 300,
  create.dir = T
)

p = DimPlot(CP_Immune_norm.prot, reduction = "ADT_UMAP", label = TRUE,
            shuffle = T, seed = seed, label.box = T, repel = T,
            label.size = 4, label.color = 'white', alpha = 0.75, group.by = 'cell_types') + ggtitle('Cell types on surface proteins')
p
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/Dim_reductions/UMAP_ADT_cell_types.pdf',
  plot = p,
  width = 12,
  height = 8,
  dpi = 300,
  create.dir = T
)

p = DimPlot(CP_Immune_norm.prot, reduction = "ADT_UMAP", label = TRUE,
            shuffle = T, seed = seed, repel = T, cols = viridis::magma(6, direction = -1, end = .9),
            alpha = 0.75, group.by = 'Age') + ggtitle('Ages on surface proteins')
p
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/Dim_reductions/UMAP_ADT_Ages.pdf',
  plot = p,
  width = 12,
  height = 8,
  dpi = 300,
  create.dir = T
)

# make results dataframe 
d = cbind(CP_Immune_norm.prot@meta.data, 
          as.data.frame(t(CP_Immune_norm.prot@assays$ADTnorm@data))
)

library(magrittr)
# calculate the median protein expression separately for each cluster 
adt_plot = d %>% 
  dplyr::group_by(cell_types) %>% 
  dplyr::summarize_at(.vars = prots, .funs = median) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("cell_types") 
# plot a heatmap of the average dsb normalized values for each cluster
htmap = pheatmap::pheatmap(t(adt_plot), 
                           color = viridis::viridis(25, option = "B"), 
                           fontsize_row = 8, border_color = NA, fontsize = 8, main = "Surface Proteins per Cell Type")
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/Heatmaps/ADT_Norm_Heatmap_cell_type.pdf',
  plot = htmap,
  width = 7,
  height = 13,
  dpi = 300,
  create.dir = T
)

# calculate the median protein expression separately for each AGE 
adt_plot = d %>% 
  dplyr::group_by(Age) %>% 
  dplyr::summarize_at(.vars = prots, .funs = median) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("Age") 
# plot a heatmap of the average dsb normalized values for each cluster
htmap = pheatmap::pheatmap(t(adt_plot), 
                           color = viridis::viridis(25, option = "B"), 
                           fontsize_row = 8, border_color = NA, fontsize = 8, main = "Surface Proteins per Age")
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/Heatmaps/ADT_Norm_Heatmap_age.pdf',
  plot = htmap,
  width = 7,
  height = 13,
  dpi = 300,
  create.dir = T
)

###############################################
######### Samir start from here #2 ############
###############################################
# rm(p,p1,p2,p3,adt_plot,Cluster_markers,d,downsampled_df,Immune_only,norm.adt,prop,testtt,htmap,age,i,cells)
# save.image(file='/Volumes/Biclab/Joseph/Development_project/CiteSeq/Code_shared_with_Samir/CiteSeq_2024.11.14.RData')
load('/Volumes/Biclab/Joseph/Development_project/CiteSeq/Code_shared_with_Samir/CiteSeq_2024.11.14.RData')

p1 = FeaturePlot(CP_Immune_norm.prot, "adtnorm_Ms.CD19", cols = c("lightgrey", "darkgreen"), order = T, shape.by = "Age") + ggtitle("CD19 Protein")
p2 = FeaturePlot(CP_Immune_norm.prot, "sct_Cd19", order = T, shape.by = "Age") + ggtitle("CD19 RNA")
p1 | p2
# Interactive plot
FeaturePlot(CP_Immune_norm.prot, "adtnorm_Ms.CD19", cols = c("lightgrey", "darkgreen"), order = T, shape.by = "Age", interactive = T)
VlnPlot(CP_Immune_norm.prot, "adtnorm_Ms.CD19", group.by = 'cell_types')

##########################
#### Diff. Expression ####  ### Pseudo-Bulk RNA ###
##########################
# Load library for the volcano plot and pp
library(EnhancedVolcano)
library(foreach)
library(doParallel)

# Function for Wilcoxon Test
DE_Wilcox_prot = function(object, ident1, ident2=NULL){
  markers = FindMarkers(object, ident.1 = ident1, ident.2 = ident2, max.cells.per.ident = 500, random.seed = seed,
                        only.pos = F, min.pct = .025, test.use = 'wilcox', densify=T)
  return (markers)
}
# Function for MAST Test
DE_MAST = function(object, ident1, ident2=NULL){
  markers = Seurat::FindMarkers(object, ident.1 = ident1, ident.2 = ident2, max.cells.per.ident = 500, random.seed = seed,
                                only.pos = F, min.pct = .025, test.use = 'MAST', densify=T)
  return (markers)
}
# Function for creating and saving volcano plots
Volcano = function(name, markers, title, comparison, strict=FALSE){
  if (strict == TRUE){logFC = 2}else{logFC = .5}
  pdf(paste0('/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/Volcano/', name, '_Volcano.pdf'),
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
Save_DE_genes = function(dataframe, filename){
  write.table(dataframe,
              file = paste0('/Volumes/Biclab/Joseph/Development_project/CiteSeq/DE_genes_per_cell_type/',filename,'.csv'),
              sep = ',', quote = F, row.names = T, col.names = NA, eol = '\r\n')
  
}
n.cores = 6
my.cluster = parallel::makeCluster(
  n.cores)
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)
tmp = CP_Immune_norm.prot
tmp$origin_cell_type = paste(CP_Immune_norm.prot$orig.ident, CP_Immune_norm.prot$cell_types, sep = "_")
Idents(tmp) = "origin_cell_type"
DefaultAssay(tmp) = "SCT"
DefaultAssay(CP_Immune_norm.prot) = "SCT"
Idents(CP_Immune_norm.prot) = "cell_types"
foreach(i = levels(CP_Immune_norm.prot@active.ident),
        .packages = c("Seurat","ggplot2",'EnhancedVolcano')) %dopar% {
          print(i)
          tryCatch({
            i_2 = gsub("/", "-", i)
            # Find Markers of the cell type in question Vs all the other cells in the dataset
            mast_d.markers = DE_MAST(CP_Immune_norm.prot, i)
            # Save the dataframes with DE info 
            Save_DE_genes(mast_d.markers, paste0('/RNA/',i_2,'_vsALL'))
            # Create and save the volcano plots
            Volcano(paste0('/RNA/',i_2,'_vs_ALL'), mast_d.markers, i_2, 'vs ALL - MAST', strict=T)
            # Find markers between the same cell type of Ages
            ages = levels(CP_Immune_norm.prot$Age)
            for (age1 in ages){
              for (age2 in ages){
                if (age1 != age2){
                  mast_d.markers2 = DE_MAST(tmp, paste(age1, i, sep = "_"), paste(age2, i, sep = "_"))
                  Save_DE_genes(mast_d.markers2, paste0('/RNA/',i_2,'_',age1,'_vs_',age2))
                  Volcano(paste0('/RNA/',i_2,'_',age1,'_vs_',age2), mast_d.markers2, i_2, paste0(age1,' vs ',age2,' - MAST test'))
                }
              }
              ages = ages[2:length(ages)]
            }
            
          })
        }

# differentially expressed proteins
DefaultAssay(tmp) = "ADTnorm"
DefaultAssay(CP_Immune_norm.prot) = "ADTnorm"
foreach(i = levels(CP_Immune_norm.prot@active.ident),
        .packages = c("Seurat","ggplot2",'EnhancedVolcano')) %dopar% {
          print(i)
          tryCatch({
            i_2 = gsub("/", "-", i)
            # Find Markers of the cell type in question Vs all the other cells in the dataset
            d.markers = DE_Wilcox_prot(CP_Immune_norm.prot, i)
            # Save the dataframes with DE info 
            Save_DE_genes(d.markers, paste0('/SurfProt/',i_2,'_vsALL'))
            # Create and save the volcano plots
            Volcano(paste0('/SurfProt/',i_2,'_vs_ALL'), d.markers, i_2, 'vs ALL - Wilcoxon', strict=T)
            # Find markers between the same cell type of Ages
            ages = levels(CP_Immune_norm.prot$Age)
            for (age1 in ages){
              for (age2 in ages){
                if (age1 != age2){
                  d.markers2 = DE_Wilcox_prot(tmp, paste(age1, i, sep = "_"), paste(age2, i, sep = "_"))
                  Save_DE_genes(d.markers2, paste0('/SurfProt/',i_2,'_',age1,'_vs_',age2))
                  Volcano(paste0('/SurfProt/',i_2,'_',age1,'_vs_',age2), d.markers2, i_2, paste0(age1,' vs ',age2,' - Wilcoxon'))
                }
              }
              ages = ages[2:length(ages)]
            }
            
          })
        }
parallel::stopCluster(cl = my.cluster)

rm(tmp)

DefaultAssay(CP_Immune_norm.prot) = "SCT"

###########################
### PHATE TRAJECTORIES ####
###########################
library(phateR)
library(viridis)
library(Rmagic)

# DefaultAssay(CP_Immune_norm.prot) = "RNA"
for (celltype in levels(CP_Immune_norm.prot$cell_types) ){
  print(celltype)
  subset_only = subset(x = CP_Immune_norm.prot, idents = c(celltype))
  matrix = t(as.data.frame(GetAssayData(object = subset_only, assay = "RNA", layer = "counts")))
  # Filtering data
  ## keep genes expressed in at least 10 cells
  keep_cols = colSums(matrix > 0) > 10
  matrix = matrix[,keep_cols]
  # keep cells with at least 1000 UMIs
  keep_rows = rowSums(matrix) > 1000
  matrix = matrix[keep_rows,]
  # Normalise
  matrix = library.size.normalize(matrix)
  matrix = sqrt(matrix)
  matrix = as.data.frame(matrix)
  # Running PHATE
  data_phate = phate(matrix, n.jobs = cores, seed = seed, npca = 50, t = 'auto')
  # plot cell types on PHASE
  p = ggplot(data_phate) +
    geom_point(aes(x=PHATE1, y=PHATE2) ) +
    theme_classic() #+
  # scale_color_manual(values = cell_colorsB, name="B Cells")
  p
  ggsave(
    paste0('/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/Dim_reductions/Phate_plots/',celltype,'/',celltype,'_on_Phate.pdf'),
    plot = p,
    width = 6.5,
    height = 5,
    dpi = 300,
    create.dir = T
  )
  # Plot by age
  p2 = ggplot(data_phate, aes(x=PHATE1, y=PHATE2, color=CP_Immune_norm.prot$Age[rownames(matrix)])) +
    geom_point(position = "jitter")+  
    scale_color_manual(values = viridis::magma(6, direction = -1, end = .9), name="Age") + theme_classic()
  p2
  ggsave(
    paste0('/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/Dim_reductions/Phate_plots/',celltype,'/',celltype,'Ages_on_Phate.pdf'),
    plot = p2,
    width = 6,
    height = 5,
    dpi = 300,
    create.dir = T
  )
  ggsave(
    paste0('/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/Dim_reductions/Phate_plots/',celltype,'/',celltype,'Ages_2_Phate_plots.pdf'),
    plot = (p+p2),
    width = 12,
    height = 5,
    dpi = 300,
    create.dir = T
  )
  ### ONLY B CELLS ###
  if (celltype == 'B cells'){
    B_genes = c("Cd19","Kit", "Spn", "Vpreb1","Vpreb3", "Dntt", "Lef1","Rag1", "Rag2",
                "Top2a", "H2afx","Tuba1b","Mki67", "Ccnb2", "Cd24a","Lmo4","Ighm", "Ms4a1","H2-Aa", "H2-Eb1","Cd74", "Ighd", 
                "Bank1", "Fcer2a","Bcl2","Pax5","Reln")
    # Visualizing imputed genes on PHATE with MAGIC
    MAGIC = magic(matrix, genes='all_genes', n.jobs = cores, seed = seed,
                  t = 'auto', decay = NULL)
    for (B in B_genes){
      print(B)
      p = ggplot(data_phate) +
        geom_point(aes(x=PHATE1, y=PHATE2, color=MAGIC$result[,B])) +
        scale_color_viridis(option="B") +
        labs(color=B) + theme_classic()
      # print(p)
      # readline(prompt="Press [enter] to continue")
      ggsave(
        paste0('/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/Dim_reductions/Phate_plots/B cells/Gene_plots/',B,'_B_cell_magic.pdf'),
        plot = p,
        width = 6,
        height = 5,
        dpi = 300,
        create.dir = T
      )
    }
  }
}

########################
##### Only B cells #####
########################
B_only = subset(x = CP_Immune_norm.prot, idents = c("B cells"))
B_only = JoinLayers(B_only, assay = "ADT")

# Run PCA
B_only = RunPCA(B_only, features = VariableFeatures(object = B_only), assay = "SCT", seed.use = seed, reduction.name = 'B_PCA')

# Determine the ‘dimensionality’ of the dataset
ElbowPlot(B_only, ndims = 50, reduction = 'B_PCA')

# Cluster the cells at elbow point
B_only = FindNeighbors(B_only, dims = 1:10, reduction = 'B_PCA', assay = 'SCT', graph.name = 'B_SCT_nn')
B_only = FindClusters(B_only, algorithm = 4, random.seed = seed, graph.name = 'B_SCT_nn', resolution = 1)

# non-linear dimensional reduction with UMAP
B_only = RunUMAP(B_only, dims = 1:10, seed.use = seed, reduction = 'B_PCA', assay = 'SCT')

# visualize UMAP
DimPlot(B_only, reduction = "umap", label = TRUE,
        label.box = T, repel = T,shuffle = T, label.size = 8, label.color = "#FFFFFF")

DimPlot(B_only, reduction = "umap", label = FALSE, label.box = T, repel = T, shuffle = T,
        group.by = 'Age', cols = viridis::magma(6, direction = -1, end = .9))+ 
  theme(legend.text=element_text(size=20))

############################################
# Finding DE features (cluster biomarkers)
B_only.markers = FindAllMarkers(B_only, only.pos = TRUE)

# Remove non B cells which are DCs 
# B_only = subset(x = B_only, idents = c("8", "9"), invert = TRUE)

B_only.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1, p_val_adj<0.05)


B_only.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1, p_val_adj<0.05) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
DoHeatmap(B_only, features = top5$gene, group.colors = cell_colorsB, size = 5 ) + NoLegend()

# Rename clusters
B_only = RenameIdents(B_only, '0' = 'Small Pre B')
B_only = RenameIdents(B_only, '1' = 'Immature B cells')
B_only = RenameIdents(B_only, '2' = 'Cycling Pre B cells')
B_only = RenameIdents(B_only, '3' = 'Mature B cells')
B_only = RenameIdents(B_only, '4' = 'Mature B cells')
B_only = RenameIdents(B_only, '5' = 'Pre Pro B cells')
B_only = RenameIdents(B_only, '6' = 'Memory B cells')
B_only = RenameIdents(B_only, '7' = 'Mature B cells')

DotPlot(B_only, features = B_genes,  scale.by = "size", cluster.idents = TRUE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#Define the colors
cell_typesB = c("Mature B cells", "Memory B cells", "Pre Pro B cells", "Cycling Pre B cells", 
                "Immature B cells", "Small Pre B")

# Define different colors for the remaining cell types based on their personalities
MatImmatB = "#6148ae"  
memB = "#452c9f"
PreProB =  "#d2b7eb"
CyclingB = "#b59ddb"
ImmatureB = "#7d64be"
SmallpreB = "#9980cd"


# Combine colors
color_paletteB = c(MatImmatB, memB, PreProB, CyclingB, ImmatureB, SmallpreB)

# Create a named vector with cell types and corresponding colours
cell_colorsB = setNames(color_paletteB, cell_typesB)
# cell_colorsB = setNames(viridis::plasma(4, end = 0.95), cell_typesB)


matrix = t(as.data.frame(GetAssayData(object = B_only, assay = "RNA", layer = "counts")))

# Filtering data
## keep genes expressed in at least 10 cells
keep_cols = colSums(matrix > 0) > 10
matrix = matrix[,keep_cols]

## look at the distribution of library sizes
ggplot() +
  geom_histogram(aes(x=rowSums(matrix)), bins=50) +
  geom_vline(xintercept = 1000, color='red')

# keep cells with at least 1000 UMIs
keep_rows = rowSums(matrix) > 1000
matrix = matrix[keep_rows,]

# Normalise
matrix = library.size.normalize(matrix)
matrix = sqrt(matrix)
matrix = as.data.frame(matrix)

sample_numbers = sample(1:length(matrix), 1000)
# PCA
matrix_PCA = as.data.frame (prcomp(matrix)$x)
ggplot(matrix_PCA) +
  geom_point(aes(PC1, PC2, color=matrix$Vpreb1)) +
  labs(color="Cd9") +
  scale_color_viridis(option="B")

# Running PHATE
data_phate = phate(matrix, n.jobs = 10, seed = seed, npca = 35 ,
                   t = 'auto', decay = NULL)

#data_phate = phate(matrix, n.jobs = cores, seed = seed, npca = 35, t = 60, mds.dist.method = 'cosine')


order = c("Pre Pro B cells", "Cycling Pre B cells", 
          "Small Pre B","Immature B cells","Mature B cells","Memory B cells")
B_only@active.ident = factor(B_only@active.ident, levels = order)

# plot cell types on PHASE
p = ggplot(data_phate) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=B_only@active.ident[rownames(matrix)] ) ) +
  scale_color_manual(values =cell_colorsB)+
  theme_classic() +theme(legend.text=element_text(size=20),legend.key.size = unit(1, 'cm'))

p

p1 = ggplot(data_phate) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=c(matrix[,"Cd19"]))) +
  labs(color="Cd19") +
  scale_color_viridis(option="B")+
  theme_classic() 

p1


ggsave(
  paste0('/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/Dim_reductions/Phate_plots/B_cells/B_cell_subtypes_on_Phate.pdf'),
  plot = p,
  width = 6.5,
  height = 5,
  dpi = 300,
  create.dir = T
)

# Plot by age
p2 = ggplot(data_phate, aes(x=PHATE1, y=PHATE2, color=CP_Immune_norm$Age[rownames(matrix)])) +
  geom_point(position = "jitter")+  
  scale_color_manual(values = viridis::magma(6, direction = -1, end = .9), name="Age") + theme_classic()
p2
ggsave(
  paste0('/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/Dim_reductions/Phate_plots/B_cells/B_cell_Ages_on_Phate.pdf'),
  plot = p,
  width = 6,
  height = 5,
  dpi = 300,
  create.dir = T
)
ggsave(
  paste0('/Volumes/Biclab/Joseph/Development_project/CiteSeq/plots/Dim_reductions/Phate_plots/B_cells/B_cell_Ages_Subtypes_Marriage_on_Phate.pdf'),
  plot = (p+p2),
  width = 12,
  height = 5,
  dpi = 300,
  create.dir = T
)

# use gene marker to visualise
B_genes = c("Cd19","Kit", "Spn", "Vpreb1","Vpreb3", "Dntt", "Lef1","Rag1", "Rag2",
            "Top2a", "H2afx","Tuba1b","Mki67", "Ccnb2", "Cd24a","Lmo4","Ighm", "Ms4a1","H2-Aa", "H2-Eb1","Cd74", "Ighd", 
            "Bank1", "Fcer2a","Bcl2","Pax5","Cd38",'Cd44','Tle3','Zeb2','Bcl6','Mzb1',"Cd52","Cxcr4","Cxcr5","Foxp1")
for (B in B_genes){
  print(B)
  B = as.character(B)
  p = ggplot(data_phate) +
    geom_point(aes(PHATE1, PHATE2, color=c(matrix[,B]))) +
    labs(color=B) +
    scale_color_viridis(option="B") + ggtitle(B) + theme_classic()
  # print(p)
  # readline(prompt="Press [enter] to continue")
  ggsave(
    paste0('/Users/saalimou/Desktop/BIC/single cell/Cite-seq P0-Ad/Analysis/Plots/Phate',B,'.pdf'),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300,
    create.dir = T
  )
}



# # Rerunning PHATE with new parameters
# data_phate = phate(matrix, n.jobs = 13, knn=7, decay=100, init=data_phate)
# ggplot(data_phate) +
#   geom_point(aes(PHATE1, PHATE2, color=matrix$Ngp)) +
#   labs(color="Ngp") +
#   scale_color_viridis(option="B")

# Visualizing imputed genes on PHATE with MAGIC
MAGIC = magic(matrix, n.jobs = cores, seed = seed, npca = 35,
              t = 'auto', decay = NULL)


ggplot(data_phate) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=MAGIC$result[,"Cxcr4"])) +
  scale_color_viridis(option="B") +
  labs(color="Cxcr4") + theme_classic()+ theme(legend.key.size = unit(2, 'cm'), legend.text=element_text(size=15))


# Phate + Magic
for (B in B_genes){
  print(B)
  p = ggplot(data_phate) +
    geom_point(aes(x=PHATE1, y=PHATE2, color=MAGIC$result[,B])) +
    scale_color_viridis(option="B") +
    labs(color=B) + theme_classic()
  # print(p)
  # readline(prompt="Press [enter] to continue")
  ggsave(
    paste0('/Users/saalimou/Desktop/BIC/single cell/Cite-seq P0-Ad/Analysis/Plots/Magic',B,'.pdf'),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300,
    create.dir = T
  )
}



DimPlot(CP_Immune_norm, reduction = "umap", label = F, label.box = F, repel = T, shuffle = T,
        group.by = 'orig.ident', cols =  viridis::magma(6, direction = -1, end = .9), pt.size = .5, seed = seed, alpha = 0.75)







Epithelium = subset(CP_Immune_norm, ident= "Epithelium")
cells = as.data.frame(Epithelium@meta.data)
x = min(table(cells$Age))
cells$barcodes = rownames(cells)
cells = cells %>% group_by(Age) %>% sample_n(x)
cells = cells$barcodes

Epithelium = subset(Epithelium, cells= cells)

Idents(Epithelium) = Epithelium$Age
Developing_Epithelium_Markers = FindAllMarkers(Epithelium,only.pos = TRUE, min.pct = 0.3)
Developing_Epithelium_Markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)


Developing_Epithelium_Markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 15) %>%
  ungroup() -> top15
DoHeatmap(Epithelium, features = top15$gene, size = 5 ) 

############# PHATE ON EPITHELIUM
matrix_epithelium = t(as.data.frame(GetAssayData(object = Epithelium, assay = "RNA", layer = "counts")))

# Filtering data
## keep genes expressed in at least 10 cells
keep_cols = colSums(matrix_epithelium > 0) > 10
matrix_epithelium = matrix_epithelium[,keep_cols]

## look at the distribution of library sizes
ggplot() +
  geom_histogram(aes(x=rowSums(matrix_epithelium)), bins=50) +
  geom_vline(xintercept = 1000, color='red')

# keep cells with at least 1000 UMIs
keep_rows = rowSums(matrix_epithelium) > 500
matrix_epithelium = matrix_epithelium[keep_rows,]

# Normalise
matrix_epithelium = library.size.normalize(matrix_epithelium)
matrix_epithelium = sqrt(matrix_epithelium)
matrix_epithelium = as.data.frame(matrix_epithelium)

sample_numbers = sample(1:length(matrix_epithelium), 1000)
# PCA
matrix_epithelium_PCA = as.data.frame (prcomp(matrix_epithelium)$x)
ggplot(matrix_epithelium_PCA) +
  geom_point(aes(PC1, PC2, color=matrix$Ttr)) +
  labs(color="Ttr") +
  scale_color_viridis(option="B")

# Running PHATE
data_phate = phate(matrix_epithelium, n.jobs = 10, seed = seed, npca = 100,knn = 6 ,
                   t = 'auto', decay = NULL)


ggplot(data_phate) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=Epithelium$percent.mt)) + scale_color_manual(values = viridis::magma(6, direction = -1, end = .9))+
  theme_classic() 


ggplot(data_phate) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=c(matrix_epithelium[,"Ins2"]))) +
  labs(color="Ins2") +
  scale_color_viridis(option="B")+
  theme_classic() 

ggplot(data_phate, aes(x=PHATE1, y=PHATE2, color=Epithelium$Age[rownames(matrix_epithelium)])) +
  geom_point(position = "jitter")+  
  scale_color_manual(values = viridis::magma(6, direction = -1, end = .9), name="Age") + theme_classic()


MAGIC = magic(matrix_epithelium, n.jobs = cores, seed = seed, npca = 35,
              t = 'auto', decay = NULL)

ggplot(data_phate) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=MAGIC$result[,"Fos"])) +
  scale_color_viridis(option="B") +
  labs(color="Fos") + theme_classic()





##### Split matrix per B and T cells
phate_TWT = phate(matrix[rownames(matrix) %in% names(T_AD$orig.ident),], n.jobs = 13, seed = seed, npca=30, t=50)
phate_TAD = phate(matrix[rownames(matrix) %in% names(T_AD$orig.ident),], n.jobs = 13, seed = seed, npca=30, t=50)
phate_BWT = phate(matrix[rownames(matrix) %in%  names(B_WT$orig.ident),], n.jobs = 13, seed = seed, npca=30, t=50)
phate_BAD = phate(matrix[rownames(matrix) %in%  names(B_AD$orig.ident),], n.jobs = 13, seed = seed, npca=30, t=50)

ggplot(phate_TWT) +
  geom_point(aes(PHATE1, PHATE2, color=c(matrix$Cd3e[rownames(phate_TWT$embedding)]))) +
  labs(color="Cd3e") +
  scale_color_viridis(option="B")

ggplot(phate_TAD) +
  geom_point(aes(PHATE1, PHATE2, color=c(matrix$Cd3e[rownames(phate_TAD$embedding)]))) +
  labs(color="Cd3e") +
  scale_color_viridis(option="B")

ggplot(phate_BWT) +
  geom_point(aes(PHATE1, PHATE2, color=c(matrix$Cd3e[rownames(phate_BWT$embedding)]))) +
  labs(color="Cd3e") +
  scale_color_viridis(option="B")

ggplot(phate_BAD) +
  geom_point(aes(PHATE1, PHATE2, color=c(matrix$Cd3e[rownames(phate_BAD$embedding)]))) +
  labs(color="Cd3e") +
  scale_color_viridis(option="B")


ggplot(phate_TWT) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=Immune_norm@active.ident[rownames(phate_TWT$embedding)] ) ) + 
  scale_color_manual(values = cell_colors) + theme_classic()

ggplot(phate_TAD) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=Immune_norm@active.ident[rownames(phate_TAD$embedding)] ) ) + 
  scale_color_manual(values = cell_colors) + theme_classic()

ggplot(phate_BWT) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=Immune_norm@active.ident[rownames(phate_BWT$embedding)] ) ) + 
  scale_color_manual(values = cell_colors) + theme_classic()

ggplot(phate_BAD) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=Immune_norm@active.ident[rownames(phate_BAD$embedding)] ) ) + 
  scale_color_manual(values = cell_colors) + theme_classic()



############################################################################################################


