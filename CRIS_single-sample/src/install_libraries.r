# Description -------------------------------------------------------------

# Script for installing required libraries.

# Configuration -----------------------------------------------------------

.INSTALL <- TRUE

# Required libraries
libraries <- c("caret",              # single-label classification
               "crayon",             # coloured print
               "e1071",              # machine learning library
               "glmnet",             # generalized linear model for Lasso
               "here",               # handle paths
               "fs",                 # handle paths
               "kernlab",            # for single-label classifiers (SVM)
               "openxlsx",           # writing multiple excel sheets
               "R6",                 # handle R classes
               "randomForest",       # execute the random forest
               "stringr",            # string utilities (string length, substring etc.)
               "tidyverse",          # handle data.frames (filter, group, etc.)
               #"RColorBrewer",       # color palette for plots
               "readxl"#,             # read excel files
               #"utiml"               # apply multi-label classification TODO: check since it is not on CRAN anymore
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

