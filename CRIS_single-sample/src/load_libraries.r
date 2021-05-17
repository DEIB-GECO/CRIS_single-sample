# Load the required libraries

library(Biobase)                  # ExpressionSet object
library(caret)                    # Single-label classification
library(crayon)                   # Colored prints
library(CRISclassifier)           # Library for the CRIS classifiers (NTP and TSP)
library(e1071)                    # Machine learning
library(glmnet)                   # Generalized linear model for lasso
library(here)                     # Relative paths
library(fs)                       # File paths manipulation
library(kernlab)                  # Single-label classifier (svm)
library(openxlsx)                 # writing multiple excel sheets
library(R6)                       # R classes
library(randomForest)             # Execute the random forest
library(RColorBrewer)             # Colors for legends in R plots
library(stringr)                  # String manipulation
library(tidyverse)                # Data.frames manipulation and %>% pipe
library(readxl)                   # Read excel files
library(utiml)                    # Multi-label classification