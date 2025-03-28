# Joseph JOSEPHIDES
# Institut Pasteur, Paris FR
# 07 Mar 2024
# Joint RNA + Surface protein analysis

# Load packages
library(Seurat)
library(SeuratObject)
library(BiocParallel)
library(Matrix)
library(ggplot2)
library(tidyverse)
library(future)
library(phateR)
library(viridis)
library(Rmagic)
# library(magrittr)

# Set Random Seed
seed = 1111111 # this does not set the seed. It will attribute a number to a variable called seed.
set.seed(seed)

# Multiprocessing parameters
# Define number of cores
cores=6
plan("multisession", workers = 8)
handlers(global = TRUE)

# Load data
load('~/Documents/data/rdata_cp/CP_base_data.RData')

p1 = FeaturePlot(CP_Immune_norm.prot, "adtnorm_Ms.CD19", cols = c("lightgrey", "darkgreen"), order = T, shape.by = "Age") + ggtitle("CD19 Protein")
p2 = FeaturePlot(CP_Immune_norm.prot, "sct_Cd19", order = T, shape.by = "Age") + ggtitle("CD19 RNA")
p1 | p2
# Interactive plot
FeaturePlot(CP_Immune_norm.prot, "adtnorm_Ms.CD19", cols = c("lightgrey", "darkgreen"), order = T, shape.by = "Age", interactive = T)
VlnPlot(CP_Immune_norm.prot, "adtnorm_Ms.CD19", group.by = 'cell_type')

#################################################
###### Weighted Nearest Neighbor Analysis #######
#################################################
Idents(CP_Immune_norm.prot) = 'cell_type'
VariableFeatures(CP_Immune_norm.prot, assay = 'ADTnorm')
VariableFeatures(CP_Immune_norm.prot, assay = 'SCT')
CP_Immune_norm.prot <- FindMultiModalNeighbors(
  CP_Immune_norm.prot, reduction.list = list("harmony", "ADT_PCA"), 
  dims.list = list(1:30, 1:17), modality.weight.name = "RNA.weight",
)
CP_Immune_norm.wnn <- RunUMAP(CP_Immune_norm.prot, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
                              reduction.key = "wnnUMAP_", seed.use = seed)

CP_Immune_norm.wnn <- FindClusters(CP_Immune_norm.wnn, 
                                   method = 'igraph',
                                   graph.name = "wsnn", 
                                   algorithm = 4, 
                                   resolution = 1, 
                                   verbose = T)


Markers = FindAllMarkers(CP_Immune_norm.wnn, min.pct = 0.1, 
                         logfc.threshold = 0.5, only.pos = T, assay = "SCT")

MarkersADT = FindAllMarkers(CP_Immune_norm.wnn, min.pct = 0.1, 
                         logfc.threshold = 0.5, only.pos = T, assay = "ADTnorm")
  #Naming Cell types
Cell_type.ids = c("Macs",  #1
                  "MHC2+ Macs",  #2
                  "Choroid plexus epithelial cells",  #3EXPRESS IL18 BUT OTHER THINGS, DO GO !
                  "Zeb2hi Macs",  #4 CAUTIOOUS
                  "Iron metabolism Macs",  #5
                  "Sall1+ Epiplexus",  #6
                  "Il31ra+ Macs",  #7
                  "Macs",  #8
                  "DAM-like Epiplexus", #9
                  "Interferon-responsive Macs", #10Express CCL12 AND CCL2 AND CCL7
                  "Metabolically active Macs", #11
                  "CCR1+ Epiplexus", #12
                  "Lyve1+ Macs", #13
                  "Neutrophils", #14
                  "Proliferating Macs", #15
                  "Monocytes", #16
                  "Developing B cells", #17
                  "Mature B cells", #18
                  "Unknown", #19
                  "CD11a+ IL10+ Macs", #20
                  "Mesenchymal cells", #21
                  "cDC2", #22
                  "NK cells", #23
                  "T cells", #24
                  "Endothelium", #25
                  "cDC1", #26
                  "Other lymphocytes", #27
                  "Htr2c+ epithelial cells", #28
                  "Unknown 2", #29 pax5 and otehr weird
                  "Mast cells", #30
                  "Basophils") #31
names(Cell_type.ids) = levels(CP_Immune_norm.wnn)
CP_Immune_norm.wnn = RenameIdents(CP_Immune_norm.wnn, Cell_type.ids)




DimPlot(CP_Immune_norm.wnn, reduction = 'wnn.umap',
              label = TRUE, repel = TRUE, label.box = T, label.size = 2.5,
              label.color = 'white', shuffle = T, seed = seed) +
  ggtitle('WNN UMAP with clusters') 
p1
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/Immune_single_cell_samir/CiteSeq/plots/Dim_reductions/UMAP_WNN_clusters.pdf',
  plot = p1,
  width = 9,
  height = 7,
  dpi = 300,
  create.dir = T
)

p2 <- DimPlot(CP_Immune_norm.wnn, reduction = 'wnn.umap', group.by = 'cell_type',
              label = TRUE, repel = TRUE, label.box = T, label.size = 2.5,
              label.color = 'white', cols = cell_colors, shuffle = T, seed = seed) +
  ggtitle('WNN UMAP with cell types') + NoLegend()
p2
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/Immune_single_cell_samir/CiteSeq/plots/Dim_reductions/UMAP_WNN_celltypes.pdf',
  plot = p2,
  width = 9,
  height = 7,
  dpi = 300,
  create.dir = T
)

p3 <- DimPlot(CP_Immune_norm.prot, reduction = 'umap', group.by = 'cell_type',
              label = TRUE, repel = TRUE, label.box = T, label.size = 2.5,
              label.color = 'white', cols = cell_colors, shuffle = T, seed = seed) +
  ggtitle('RNA UMAP with cell types') + NoLegend()
p3
p4 <- DimPlot(CP_Immune_norm.prot, reduction = 'ADT_UMAP', group.by = 'cell_type',
              label = TRUE, repel = TRUE, label.box = T, label.size = 2.5,
              label.color = 'white', cols = cell_colors, shuffle = T, seed = seed) +
  ggtitle('Surface protein UMAP with cell types') + NoLegend()
p4
p3 + p4 + p2
ggsave(
  '/Volumes/Biclab/Joseph/Development_project/Immune_single_cell_samir/CiteSeq/plots/Dim_reductions/UMAP_RNA+ADT+WNN_celltypes.pdf',
  plot = p3 + p4 + p2,
  width = 17,
  height = 6,
  dpi = 300,
  create.dir = T
)

rm(p1,p2,p3,p4,p)

# Investigate Proliferating macs 2 populations
CP_Immune.wnn = CP_Immune_norm.prot
CP_Immune.wnn$cell_type = as.character(CP_Immune.wnn$cell_type)
CP_Immune.wnn$cell_type[CP_Immune.wnn$cell_type == "Proliferating Macs"] = 'Right Proliferating Macs'
CP_Immune.wnn$cell_type[CP_Immune.wnn$wsnn_res.2 == 37] = 'Left Proliferating Macs'
CP_Immune.wnn$cell_type = as.factor(CP_Immune.wnn$cell_type)
Idents(CP_Immune_norm.prot) = 'cell_type'
RNA.macs.markers.r = FindMarkers(CP_Immune.wnn, group.by = 'cell_type', assay = 'SCT', ident.1 = 'Right Proliferating Macs', ident.2 = 'Left Proliferating Macs') %>% dplyr::filter(avg_log2FC > .5)
ADT.macs.markers.r = FindMarkers(CP_Immune.wnn, group.by = 'cell_type', assay = 'ADTnorm', ident.1 = 'Right Proliferating Macs', ident.2 = 'Left Proliferating Macs') %>% dplyr::filter(avg_log2FC > .5)
RNA.macs.markers.l = FindMarkers(CP_Immune.wnn, group.by = 'cell_type', assay = 'SCT', ident.1 = 'Left Proliferating Macs', ident.2 = 'Right Proliferating Macs') %>% dplyr::filter(avg_log2FC > .5)
ADT.macs.markers.l = FindMarkers(CP_Immune.wnn, group.by = 'cell_type', assay = 'ADTnorm', ident.1 = 'Left Proliferating Macs', ident.2 = 'Right Proliferating Macs') %>% dplyr::filter(avg_log2FC > .5)

head(RNA.macs.markers.l)
head(RNA.macs.markers.r)

head(ADT.macs.markers.r)
head(ADT.macs.markers.l)

write.table(RNA.macs.markers.l, '/Volumes/Biclab/Joseph/Development_project/Immune_single_cell_samir/CiteSeq/DE_genes_per_cell_type/Left-right-macrophages/RNA.macs.markers.l.tsv', quote = F, sep = '\t', row.names = T, col.names = T)
write.table(RNA.macs.markers.r, '/Volumes/Biclab/Joseph/Development_project/Immune_single_cell_samir/CiteSeq/DE_genes_per_cell_type/Left-right-macrophages/RNA.macs.markers.r.tsv', quote = F, sep = '\t', row.names = T, col.names = T)
write.table(ADT.macs.markers.r, '/Volumes/Biclab/Joseph/Development_project/Immune_single_cell_samir/CiteSeq/DE_genes_per_cell_type/Left-right-macrophages/ADT.macs.markers.r.tsv', quote = F, sep = '\t', row.names = T, col.names = T)
write.table(ADT.macs.markers.l, '/Volumes/Biclab/Joseph/Development_project/Immune_single_cell_samir/CiteSeq/DE_genes_per_cell_type/Left-right-macrophages/ADT.macs.markers.l.tsv', quote = F, sep = '\t', row.names = T, col.names = T)

# rm(ADT.macs.markers.l, ADT.macs.markers.r, RNA.macs.markers.l, RNA.macs.markers.r, ADT.VariableFeatures, data_magic, data_phate, Macs_only, MAGIC, matrix, p1, p, p2, p3, RNA.VariableFeatures, subset_only, g2m.genes,
#    immune_colors, isotypes, min_count, order, original_counts, prots, s.genes, x, marker_genes, keep_rows, keep_cols, celltype, CP_Immune.wnn)
save.image(file='/Volumes/Biclab/Joseph/Development_project/Immune_single_cell_samir/Data_and_code/R_data/CiteSeq_WNN_2025.03.13.RData')
##########################
#### Diff. Expression ####  ### Pseudo-Bulk RNA ###
##########################
DefaultAssay(CP_Immune_norm.prot) = "SCT"
Idents(CP_Immune_norm.prot) = 'cell_type'
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
  pdf(paste0('/Volumes/Biclab/Joseph/Development_project/Immune_single_cell_samir/CiteSeq/plots/Volcano/', name, '_Volcano.pdf'),
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
              file = paste0('/Volumes/Biclab/Joseph/Development_project/Immune_single_cell_samir/CiteSeq/DE_genes_per_cell_type/',filename,'.csv'),
              sep = ',', quote = F, row.names = T, col.names = NA, eol = '\r\n')
  
}
n.cores = 6
my.cluster = parallel::makeCluster(
  n.cores)
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)
tmp = CP_Immune_norm.prot
tmp$origin_cell_type = paste(CP_Immune_norm.prot$orig.ident, CP_Immune_norm.prot$cell_type, sep = "_")
Idents(tmp) = "origin_cell_type"
DefaultAssay(tmp) = "SCT"
DefaultAssay(CP_Immune_norm.prot) = "SCT"
Idents(CP_Immune_norm.prot) = "cell_type"
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
DefaultAssay(CP_Immune_norm.prot) = "RNA"
Idents(CP_Immune_norm.prot) = 'cell_type'
for (celltype in levels(CP_Immune_norm.prot$cell_type) ){
  print(celltype)
  subset_only = subset(x = CP_Immune_norm.prot, idents = celltype)
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
    paste0('/Volumes/Biclab/Joseph/Development_project/Immune_single_cell_samir/CiteSeq/plots/Dim_reductions/Phate_plots/',celltype,'/',celltype,'_on_Phate.pdf'),
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
    paste0('/Volumes/Biclab/Joseph/Development_project/Immune_single_cell_samir/CiteSeq/plots/Dim_reductions/Phate_plots/',celltype,'/',celltype,'Ages_on_Phate.pdf'),
    plot = p2,
    width = 6,
    height = 5,
    dpi = 300,
    create.dir = T
  )
  ggsave(
    paste0('/Volumes/Biclab/Joseph/Development_project/Immune_single_cell_samir/CiteSeq/plots/Dim_reductions/Phate_plots/',celltype,'/',celltype,'Ages_2_Phate_plots.pdf'),
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
        paste0('/Volumes/Biclab/Joseph/Development_project/Immune_single_cell_samir/CiteSeq/plots/Dim_reductions/Phate_plots/B cells/Gene_plots/',B,'_B_cell_magic.pdf'),
        plot = p,
        width = 6,
        height = 5,
        dpi = 300,
        create.dir = T
      )
    }
  }
}

DefaultAssay(CP_Immune_norm.prot) = "SCT"

##### Analysis for Aleks - 13 Mar 2025
Idents(CP_Immune_norm.prot) = 'cell_type'
DefaultAssay(CP_Immune_norm.prot) = "RNA"
Macs_only = subset(x = CP_Immune_norm.prot, idents = c( "Monocytes", "Macs 1", "Macs 2", "Proliferating Macs")) #"Epiplexus",
matrix = t(as.data.frame(GetAssayData(object = Macs_only, assay = "RNA", layer = "counts")))
keep_cols = colSums(matrix > 0) > 10
matrix = matrix[,keep_cols]
keep_rows = rowSums(matrix) > 1000
matrix = matrix[keep_rows,]
matrix = library.size.normalize(matrix)
matrix = sqrt(matrix)
matrix = as.data.frame(matrix)
# Running PHATE
data_magic = magic(matrix, n.jobs = cores, seed = seed, t = 'auto')
data_phate = phate(matrix, n.jobs = cores, seed = seed, t = 'auto')
# plot cell types on PHASE
p1 = ggplot(data_phate) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=Macs_only$Age[rownames(matrix)])) +
  scale_color_manual(values = viridis::magma(6, direction = -1, end = .9), name="Age") + theme_classic()
p1

p2 = ggplot(data_phate) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=Macs_only$cell_type[rownames(matrix)])) +
  scale_color_manual(values = cell_colors, name = 'Cell type') + theme_classic() 
p2
p1+p2

############################################################################################################

write.table(capture.output(sessionInfo()),
            file = paste0('/Volumes/Biclab/Joseph/Development_project/Immune_single_cell_samir/Data_and_code/R_data/SessionInfo_WNN_',Sys.Date(),'.txt'),
            quote = F, row.names = F, col.names = F)


