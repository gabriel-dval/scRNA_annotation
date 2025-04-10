#' Even with the best cell type annotation methods, we still only get between
#' 80 and 90% accuracy - without any form of confidence score, this could be
#' a problem
#' 
#' 2 solution tested here : 
#' 
#' Establish a confidence score for the annotation
#' 
#' If dealing with large number of clusters, merge annotations of clusters
#' that closely ressemble each other


###############################################################################
# A - CONFIDENCE SCORE ########################################################
###############################################################################

library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(wordcloud)
library(grid)
library(Seurat)
library(SeuratDisk)
library(SeuratData)

# Load annotations and original data
load('../data/rdata_cp/CP_base_data.RData')  
gpt_annots <- read.table('results/cp_gptcelltype_annotations.csv', sep = ';')

# Info : 
# V1 : Manual annotation
# V2 - V6 : 40 marker genes, LogFold > 0.5
# V7 - V11 : 40 marker genes, LogFold > 0
# V12 - V16 : 30 marker genes, LogFold > 0.5
# V17 - V21 : 30 marker genes, LogFold > 0

# First step - get all words in the annotation

clean_words <- function(string_vec) {
  #' Convert a character vector into a vector of its component
  #' single words, in lower case and all punct/numbers removed.
  #' 
  #' Args
  #' ---
  #' string_vec : data.frame row
  #' 
  #' Returns
  #' --------
  #' single_words : character vector
  #' 
  
  # Convert to vector
  string_vec <- unlist(string_vec) %>%
    as.vector()
  
  # Convert to lowercase
  string_vec <- tolower(string_vec)
  
  # Replace separators (space, '/', and ',') with a single space
  string_vec <- gsub("[ /,]", " ", string_vec)
  
  # Remove any non-alphabetic characters EXCEPT spaces (preserve word boundaries)
  string_vec <- gsub("[[:digit:][:punct:]]", "", string_vec)
  
  # Trim trailing whitespace
  single_words <- trimws(string_vec, whitespace = "[\\h\\v]")
  
  # Split into individual words using whitespace
  single_words <- unlist(strsplit(single_words, "[ \u00A0]+"))
  
  # Remove plurals
  single_words <- sub("s$", "", single_words)
  
  return(single_words)
}


calculate_proportions <- function(clean_word_vector) {
  #' Function to calculate the proportion of each word appearing in
  #' a phrase and plot a corresponding pie chart
  #' 
  #' Args
  #' ----
  #' clean_word_vector : character vector
  #' 
  #' Returns
  #' -------
  #' proportions : list of two - data.frame and graph object
  
  # Remove words that are not useful (e.g : 'cells')
  clean_word_vector <- clean_word_vector[!(clean_word_vector %in% c("cell", 
                                                                    "cells"))]
  
  # Calculate proportions of the vector
  word_frequency <- table(clean_word_vector)
  word_frequency <- as.data.frame(word_frequency)
  colnames(word_frequency) <- c('Label', 'Count')
  
  # Create pie chart using ggplot2
  p <- ggplot(word_frequency, aes(x = "", y = Count, fill = Label)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +  # Convert bar chart to pie chart
    labs(title = "Word cloud", fill = 'Label') +
    theme_void() + 
    geom_text(aes(label = Count), 
              position = position_stack(vjust = 0.5), 
              color = "white")
  
  # Get heterogeneity score
  h_score <- length(word_frequency$Label)
  
  
  # Create output object
  res <- list(word_frequency, p, h_score)
  names(res) <- c('proportions', 'pie_chart', 'heterogeneity')
  
  return(res)
  
}


# Test on a single row
row_1 <- gpt_annots[37,2:21] %>%
  unlist() %>%
  as.vector()

test <- clean_words(row_1)
prop_test <- calculate_proportions(test)
prop_test


# Now test function by applying it over matrix
res <- apply(X = gpt_annots[,2:21], MARGIN = 1, FUN = clean_words)
res <- lapply(X = res, FUN = calculate_proportions)


# wordcloud function example for downstream
wordcloud(words = res[[23]]$proportions$Label, 
                freq = res[[23]]$proportions$Count, min.freq = 1, 
                scale = c(3,0.5), colors = brewer.pal(8, "Dark2"))



# Plot some examples
res[[37]]$pie_chart


# From this, we can do two things : 
# Calculate a heterogeneity score
# Set the most common word as the dominant annotation, 
# and allow for more differentials using this heterogeneity
# score ? 

select_annotation <- function(single_cluster_annot) {
  #' Function to select most appropriate annotation
  #' 
  #' Args
  #' ---
  #' list_of_cluster_annots : list of 3 elements
  #'  Output of the calculate_proportions function
  #'  
  #' Returns
  #' -------
  #' annotation : character
  #'  Annotation for the cluster
  
  # Get most common label
  props <- single_cluster_annot$proportions %>%
    arrange(desc(Count))
  most_common <- as.character(props$Label[1])
  
  if(single_cluster_annot$heterogeneity >= 8) {
    differential_1 <- props$Label[2]
    differential_2 <- props$Label[3]
      most_common <- paste(most_common, differential_1, 
                           differential_2, sep = ' + ')
  }
  
  else if(single_cluster_annot$heterogeneity >= 2) {
    differential <- props$Label[2]
    most_common <- paste(most_common, differential, sep = ' + ')
  }
  
  return(most_common)
}


# Now apply this to all clusters
res <- lapply(X = res, FUN = select_annotation)

# Get clusters
annotations <- unlist(res)







###############################################################################
# B - TEST ON CP DATA #########################################################
###############################################################################

# Load data
load('../data/rdata_cp/CP_base_data.RData')


# Make equivalency vector
cluster_annotation <- setNames(annotations, 
                               levels(Idents(CP_Immune_norm)))
ct <- cluster_annotation %>% as.data.frame()

# Alternatively
names(cluster_annotation)
names(cluster_annotation) <- levels(CP_Immune_norm)
CP_Immune_norm[["old.ident"]] <- as.factor(Idents(CP_Immune_norm) )
levels(CP_Immune_norm$old.ident)
CP_Immune_norm <- RenameIdents(CP_Immune_norm, cluster_annotation)
levels(CP_Immune_norm@active.ident)

# UMAP
DimPlot(CP_Immune_norm, reduction = 'umap', label = TRUE, 
        pt.size = 0.4, shuffle = T, seed = seed, label.box = T, repel = T, 
        label.size = 4, label.color = 'black', alpha = 0.75) + NoLegend() +
  ggtitle('GPTCellType corrected, GPT 4.5 Preview, detailedprompt2')
FeaturePlot(CP_Immune_norm, features = 'Ly6g')

ggsave('results/cp_gptcelltype_corrected_annotation.pdf')


View(as.data.frame(annotations))











###############################################################################
# C - Look at wordcloud of ref datasets #######################################
###############################################################################

# Original dataset
load('../data/rdata_sbm/BoneMarrow_data_2025.RData') # SBM
Immune_norm <- LoadH5Seurat('../data/raw_sbm/manual_annot.h5seurat')
wf_manual <- table(Immune_norm@meta.data$cluster_annot)
wf_manual <- as.data.frame(wf_manual)
colnames(wf_manual) <- c('Label', 'Count')
wf_manual$Label <- sub("s$", "", as.character(wf_manual$Label))


# Load all the different labels of reference datasets
# FACS
labels <- read_tsv(file = '../data/ref_sbm/mATLAS_Marrow_facs_10x/cluster_names.tsv',
                   col_names = F)
clean_labels_facs <- sub(".*_", "", labels$X1) %>% as.data.frame()
wf_facs <- table(clean_labels_facs)
wf_facs <- as.data.frame(wf_facs)
colnames(wf_facs) <- c('Label', 'Count')
wf_facs$Label <- sub("s$", "", as.character(wf_facs$Label))


# Droplet
labels_droplet <- read_tsv(file = paste('../data/ref_sbm/mATLAS_Marrow_droplet_10x/',
                                        'cluster_names.tsv', sep = ''),
                           col_names = F)
clean_labels_droplet <- sub(".*_", "", labels_droplet$X1) %>% as.data.frame()
wf_droplet <- table(clean_labels_droplet)
wf_droplet <- as.data.frame(wf_droplet)
colnames(wf_droplet) <- c('Label', 'Count')
wf_droplet$Label <- sub("s$", "", as.character(wf_droplet$Label))


# PangLao
ref_clusters <- read.table(paste('../data/ref_sbm/PangLao_DB_SRA653146/',
                                 'SRA653146.clusters.txt', sep = ''), sep = ' ')

ref_cluster_annotations <- c(
  '0' =	'Hematopoietic stem cells',
  '1' =	'B cells',
  '2' = 'Neutrophils',
  '3' = 'B cells',
  '4' = 'Hematopoietic stem cells',
  '5' = 'Neutrophils',
  '6' = 'B cells',
  '7' = 'Macrophages',
  '8' = 'Hematopoietic stem cells',
  '9' = 'Neutrophils',
  '10' = 	'B cells',
  '11' = 	'T memory cells',
  '12' = 	'Plasmacytoid dendritic cells',
  '13' = 	'Dendritic cells',
  '14' = 	'Unknown',
  '15' = 	'T memory cells',
  '16' = 	'Basophils',
  '17' = 	'B cells',
  '18' = 	'B cells'
)

labels_panglao <- ref_cluster_annotations[as.character(ref_clusters$V2)] %>% 
  as.data.frame()
wf_panglao <- table(labels_panglao)
wf_panglao <- as.data.frame(wf_panglao)
colnames(wf_panglao) <- c('Label', 'Count')
wf_panglao$Label <- sub("s$", "", as.character(wf_panglao$Label))


# CellXGene
labels_cellxgene <- read_tsv(file=paste('../data/ref_sbm/BoneMarrow_multi_assay_10x/',
                                          'cluster_names.tsv', sep = ''),
                             col_names = F) %>% as.data.frame()
wf_cellxgene <- table(labels_cellxgene)
wf_cellxgene <- as.data.frame(wf_cellxgene)
colnames(wf_cellxgene) <- c('Label', 'Count')
wf_cellxgene$Label <- sub("s$", "", as.character(wf_cellxgene$Label))


# Tabula Muris
annot_tm <- read_csv('../data/ref_sbm/TM_facs/annotations_FACS.csv')
annot_tm <- annot_tm[annot_tm[,22] == 'Marrow', ]
tm_labels <- annot_tm$cell_ontology_class
wf_tm <- table(tm_labels)
wf_tm <- as.data.frame(wf_tm)
colnames(wf_tm) <- c('Label', 'Count')
wf_tm$Label <- sub("s$", "", as.character(wf_tm$Label))



wordcloud(words = word_frequency$Label, freq = word_frequency$Count, min.freq = 5, 
          scale = c(3,0.5), colors = brewer.pal(8, "Dark2"))



# Find common annotations
common <- intersect(tolower(wf_manual$Label), tolower(wf_s$Label))


common <- Reduce(intersect, list(tolower(wf_manual$Label), 
                                 tolower(wf_droplet$Label),
                       tolower(wf_facs$Label), tolower(wf_panglao$Label),
                       tolower(wf_cellxgene$Label), tolower(wf_tm$Label)))





