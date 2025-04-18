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
library(SnowballC)
library(tm)
library(textstem)

# Load annotations and original data
load('../data/rdata_cp/CP_base_data.RData')  
load('../data/rdata_sbm/BoneMarrow_data_2025.RData') # SBM
gpt_annots <- read.table('output_sbm_2/gptcelltype_annots.csv', sep = ',')
gpt_annots <- gpt_annots[-1,-1]
rownames(gpt_annots) <- NULL
colnames(gpt_annots) <- NULL

# Info : 
# V1 : Manual annotation
# V2 - V6 : 30 marker genes, LogFold > 0
# V7 - V11 : 30 marker genes, LogFold > 0.5
# V12 - V16 : 40 marker genes, LogFold > 0
# V17 - V21 : 40 marker genes, LogFold > 0.5

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
  string_vec <- gsub("[[:punct:]]", "", string_vec)
  
  # This regex looks for numbers that have space or boundary before and after them
  string_vec <- gsub("(^|\\s)\\d+(\\s|$)", " ", string_vec)
  
  # Trim trailing whitespace
  single_words <- trimws(string_vec, whitespace = "[\\h\\v]")
  
  # Split into individual words using whitespace
  single_words <- unlist(strsplit(single_words, "[ \u00A0]+"))
  
  # Remove plurals
  single_words <- sub("s$", "", single_words)
  
  return(single_words)
}

# Remove useless words
remove_useless_words <- function(clean_word_vector) {
  #' Function to eliminate numbers, useless words or any other
  #' expression not useful to the cell type annotation (e.g : 
  #' words with the same root such as neuron and neuronal)
  #' 
  #' Args
  #' ----
  #' clean_word_vector : character vector
  #' 
  #' Returns 
  #' ------
  #' useful_words : character vector
  
  # Filter out words that are only numeric
  clean_word_vector <- clean_word_vector[!grepl("^\\d+$", clean_word_vector)]
  
  # Remove words that are not useful (e.g : 'cells')
  clean_word_vector <- clean_word_vector[!(clean_word_vector %in% c("cell", 
                                                                    "cells"))]
  
  # Remove common stopwords
  stopwords <- tm::stopwords("en")
  clean_word_vector <- clean_word_vector[!(clean_word_vector %in% stopwords)]
}

# Proportions
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

# Basic root
same_root <- function(word1, word2) {
  # Check if inputs are valid strings
  if (!is.character(word1) || !is.character(word2) || 
      length(word1) != 1 || length(word2) != 1) {
    stop("Both inputs must be single character strings")
  }
  
  # Load the SnowballC package
  # Install if needed with: install.packages("SnowballC")
  if (!requireNamespace("SnowballC", quietly = TRUE)) {
    stop("The SnowballC package is required. Please install it with install.packages('SnowballC')")
  }
  
  # Convert to lowercase
  word1 <- tolower(word1)
  word2 <- tolower(word2)
  
  # Use Porter stemming algorithm from SnowballC
  stem1 <- SnowballC::wordStem(word1, language = "english")
  stem2 <- SnowballC::wordStem(word2, language = "english")
  
  # Check if stems are the same
  result <- stem1 == stem2
  
  cat("Word 1:", word1, "→ Stem:", stem1, "\n")
  cat("Word 2:", word2, "→ Stem:", stem2, "\n")
  
  return(result)
}

# Lemmatization test
same_root_lemma <- function(word1, word2) {
  # Check if inputs are valid strings
  if (!is.character(word1) || !is.character(word2) || 
      length(word1) != 1 || length(word2) != 1) {
    stop("Both inputs must be single character strings")
  }
  
  # Convert to lowercase
  word1 <- tolower(word1)
  word2 <- tolower(word2)
  
  # Custom dictionary for scientific/biological terms
  bio_lemma_dict <- list(
    "epithelial" = "epithelium",
    "neuronal" = "neuron",
    "dendritic" = "dendrite",
    "axonal" = "axon",
    "glial" = "glia",
    "astrocytic" = "astrocyte",
    "myelinated" = "myelin",
    "somatic" = "soma",
    "endothelial" = "endothelium",
    "mesenchymal" = "mesenchyme",
    "cytoplasmic" = "cytoplasm",
    "nuclear" = "nucleus",
    "vascular" = "vasculature"
  )
  
  # Helper function to get lemma
  get_lemma <- function(word) {
    if (word %in% names(bio_lemma_dict)) {
      return(bio_lemma_dict[[word]])
    }
    
    # Try to lemmatize with textstem if available
    if (requireNamespace("textstem", quietly = TRUE)) {
      lemmas <- textstem::lemmatize_words(word)
      return(lemmas[1]) # Take first lemma if multiple returned
    }
    
    # Fallback to Porter stemming if textstem not available
    if (requireNamespace("SnowballC", quietly = TRUE)) {
      return(SnowballC::wordStem(word, language = "english"))
    }
    
    # If nothing else works, return the original word
    return(word)
  }
  
  # Get lemmas
  lemma1 <- get_lemma(word1)
  lemma2 <- get_lemma(word2)
  
  # Print for debugging/information
  cat("Word 1:", word1, "→ Lemma:", lemma1, "\n")
  cat("Word 2:", word2, "→ Lemma:", lemma2, "\n")
  
  # Check if lemmas are the same
  return(lemma1 == lemma2)
}


###############################################################################
# Alternatively

find_most_common_expression <- function(expressions, ignore_case = TRUE,
                                        standardize_punctuation = TRUE,
                                        standardize_plural = TRUE,
                                        standardize_conjunctions = TRUE) {
  
  # Input validation
  if(!is.vector(expressions) || !is.character(expressions)) {
    if (is.list(expressions)) {
      expressions <- as.character(expressions)
    } else {
      stop('Incorrect input format - character vector or list only')
    }
  }
  
  if(length(expressions) == 0) {
    return(NULL)
  }
  
  # Create a clean version of the expressions for normalization
  clean_expressions <- expressions
  
  # Apply case normalization if requested
  if(ignore_case) {
    clean_expressions <- tolower(clean_expressions)
  }
  
  # Standardize punctuation if requested
  if(standardize_punctuation) {
    # Replace slashes and plus signs with spaces
    clean_expressions <- gsub("/|\\+", " ", clean_expressions)
    # Remove other punctuation
    clean_expressions <- gsub("[[:punct:]]", "", clean_expressions)
  }
  
  # Replace multiple spaces with single space
  clean_expressions <- gsub("\\s+", " ", clean_expressions)
  
  # Trim leading and trailing whitespace
  clean_expressions <- trimws(clean_expressions)
  
  # Standardize plurals if requested
  if(standardize_plural) {
    # Function to selectively remove trailing 's' from words
    remove_plural <- function(str) {
      words <- unlist(strsplit(str, " "))
      words <- ifelse(nchar(words) > 2 & substr(words, nchar(words), nchar(words)) == "s", 
                      substr(words, 1, nchar(words)-1), 
                      words)
      paste(words, collapse = " ")
    }
    
    clean_expressions <- sapply(clean_expressions, remove_plural)
  }
  
  # Standardize conjunctions if requested
  if(standardize_conjunctions) {
    # Replace variations of "and" and "or" with a standard form
    clean_expressions <- gsub(" and | & ", " and ", clean_expressions)
    clean_expressions <- gsub(" or ", " or ", clean_expressions)
    
    # Optionally, completely remove conjunctions to group similar expressions
    # Uncomment the line below if you want to remove conjunctions entirely
    # clean_expressions <- gsub(" and | & | or ", " ", clean_expressions)
    
    # Replace multiple spaces again after conjunction standardization
    clean_expressions <- gsub("\\s+", " ", clean_expressions)
    clean_expressions <- trimws(clean_expressions)
  }
  
  # Count the standardized expressions
  expr_counts <- table(clean_expressions)
  
  # Find the most common expression(s)
  max_count <- max(expr_counts)
  most_common_expr <- names(expr_counts)[expr_counts == max_count]
  
  # Map back to the original expressions to maintain original formatting
  # Create a mapping from clean expressions to original expressions
  expr_mapping <- data.frame(
    original = expressions,
    clean = clean_expressions,
    stringsAsFactors = FALSE
  )
  
  # Find all original expressions that match the most common clean expression
  matching_originals <- expr_mapping$original[expr_mapping$clean == most_common_expr[1]]
  
  # Return the first one (or you could return all matches)
  return(list(
    expression = matching_originals[1],
    normalized = most_common_expr[1],
    count = max_count,
    all_matches = matching_originals
  ))
}


###############################################################################
# Test these functions

# Test alternative function
test <- find_most_common_expression(res[31, -1])



# Test on a single row
test <- clean_words(c(gpt_annots[37,2:21], 'CD8', ' and', 'with'))
print(test)
test <- remove_useless_words(test)
print(test)


prop_test <- calculate_proportions(test)
prop_test


# Now test function by applying it over matrix
res <- apply(X = gpt_annots[,2:21], MARGIN = 1, FUN = clean_words)
res <- lapply(X = res, FUN = remove_useless_words)
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
  
  if(single_cluster_annot$heterogeneity >= 2) {
    differential <- as.character(props$Label[2])
    print(most_common)
    print(differential)
    i <- 3
    while(same_root_lemma(most_common, differential) == T) {
      differential <- as.character(props$Label[i])
      i <- i + 1
    }
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



wordcloud(words = wf_manual$Label, freq = wf_manual$Count, min.freq = 5, 
          scale = c(3,0.5), colors = brewer.pal(8, "Dark2"))



# Find common annotations
common <- intersect(tolower(wf_manual$Label), tolower(wf_cellxgene$Label))


common <- Reduce(intersect, list(tolower(wf_manual$Label), 
                                 tolower(wf_droplet$Label),
                       tolower(wf_facs$Label), tolower(wf_panglao$Label),
                       tolower(wf_cellxgene$Label), tolower(wf_tm$Label)))












###############################################################################
# D - Build heatmap of results per cluster ####################################
###############################################################################
library(plotly)
library(reshape2)


# Run pipeline until calculate proportion function
post_res <- apply(X = gpt_annots, MARGIN = 1, FUN = clean_words)
post_res <- lapply(X = post_res, FUN = remove_useless_words)
post_res <- lapply(X = post_res, FUN = calculate_proportions)

# Calculate size of the matrix
# Take all clean words generated by ChatGPT by iterating through all available
# annotations, concatenate and get unique ones. Create a named vector with all
# these - this is the size of the x axis

# All labels
all_labels <- lapply(X = post_res, 
                     FUN = function(i) {as.character(i$proportions$Label)})
all_labels <- str_to_title(unique(unlist(all_labels)))
print(all_labels)


# Fill heatmap function
fill_heatmap <- function(all_labels, annotations, index) {
  
  # Initiate fill vector
  result <- rep(0, length(all_labels))
  
  # Get annots
  annots <- annotations[[index]]$proportions
  
  # Fill in vector
  for(i in 1:nrow(annots)) {
    term <- str_to_title(as.character(annots$Label[i]))
    count <- annots$Count[i]
    
    position <- which(all_labels == term)
    print(position)
    
    if(length(position) > 1) {
      stop('SOMETHING IS WRONG')
    }
    
    # Fill in position in vector
    result[position] <- count
  }
  
  return(result)
}


# Run
heatmap_filled <- lapply(
  1:nrow(heatmap_data),
  function(i) fill_heatmap(all_labels, post_res, i)
)

# Concatenate 
heatmap_data <- do.call(rbind, heatmap_filled) %>%
  as.data.frame() 
heatmap_data <- data.frame(lapply(heatmap_data, 
                                  function(x) as.numeric(x)))
colnames(heatmap_data) <- all_labels
heatmap_data$Row <- 1:nrow(heatmap_data)

# Convert matrix to long format
heatmap_long <- melt(heatmap_data, id.vars = 'Row')
names(heatmap_long) <- c("Row", "Column", "Value")

# Create with ggplot first
p <- ggplot(heatmap_long, aes(x = Column, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_viridis_c(option = 'magma') +
  scale_y_continuous(breaks = 1:nrow(heatmap_data)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Convert to plotly
ggplotly(p)

