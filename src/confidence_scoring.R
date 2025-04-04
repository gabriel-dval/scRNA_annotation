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

row_1 <- gpt_annots[1,2:21] %>%
  as.vector() %>%
  unlist()

print(row_1)


















