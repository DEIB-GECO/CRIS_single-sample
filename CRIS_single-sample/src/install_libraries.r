# Description -------------------------------------------------------------

# Script for installing required libraries.

# Configuration -----------------------------------------------------------

.INSTALL <- FALSE

# Required libraries
libraries <- c("crayon",             # coloured print
               "here",               # handle paths
               "fs",                 # handle paths
               "R6",                 # handle R classes
               "stringr",            # string utilities (string length, substring etc.)
               "tidyverse",          # handle data.frames (filter, group, etc.)
               "readxl"              # read excel files
              )
  

# Installation ------------------------------------------------------------

if (.INSTALL){
  
  # BiocManager
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  # Biobase 
  if (!require("Biobase"))
    BiocManager::install("Biobase")
  
  for (lib in libraries){
    if (!require(lib))
      install.packages(lib)
  }
  
}

