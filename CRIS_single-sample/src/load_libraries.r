# Load the required libraries

library(Biobase)                  # ExpressionSet object
library(caret)                    # Single-label classification
library(crayon)                   # Colored prints
library(glmnet)                   # Generalized linear model for lasso
library(here)                     # Relative paths
library(fs)                       # File paths manipulation
library(R6)                       # R classes
library(stringr)                  # String manipulation
library(tidyverse)                # Data.frames manipulation and %>% pipe
library(readxl)                   # Read excel files
library(utiml)                    # Multi-label classification