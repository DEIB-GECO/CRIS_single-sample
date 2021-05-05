# Description -------------------------------------------------------------

# Script for installing required libraries.

# Configuration -----------------------------------------------------------

.INSTALL <- TRUE

# Required libraries
libraries <- c("caret",              # single-label classification
               "crayon",             # coloured print 
               "glmnet",             # generalized linear model for Lasso
               "here",               # handle paths
               "fs",                 # handle paths
               "R6",                 # handle R classes
               "stringr",            # string utilities (string length, substring etc.)
               "tidyverse",          # handle data.frames (filter, group, etc.)
               "readxl",             # read excel files
               "utiml"               # apply multi-label classification
              )
  

# Installation ------------------------------------------------------------

if (.INSTALL){
  
  # BiocManager
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  # Biobase (ExpressionSet)
  if (!require("Biobase"))
    BiocManager::install("Biobase")
  
  for (lib in libraries){
    if (!require(lib))
      install.packages(lib)
  }
  
}

