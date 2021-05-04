# Libraries ---------------------------------------------------------------

# Class
library(R6)

# Expression Set
library(Biobase)

# Data manipulation and %>% pipe
library(tidyverse)

# Single-label classification
library(caret)


# Class definition --------------------------------------------------------

SLData <- R6Class(
  
  'SLData',
  inherit      = Data,
  lock_objects = TRUE,
  lock_class   = TRUE,
  
  #### Public fields ####
  
  private = list(
    
    .prepare_data = function(){
      
      private$.ref <- private$.ref %>% 
        select(ALIQUOT_LABEL, CLASS_LABEL, private$.classes, CLASS_FDR_LABEL) %>%
        as.data.frame()
      private$.ref[,CLASS_LABEL] <- factor(private$.ref[,CLASS_LABEL], levels = levels(F_CRIS_CLASSES))
      rownames(private$.ref) <- private$.ref[,ALIQUOT_LABEL]
      
      
      # Get expression (samples x genes)
      expr <- as.data.frame(t(private$.data@assayData$exprs))
      expr <- expr %>% rownames_to_column(ALIQUOT_LABEL) %>% as.data.frame()
      
      # Join with labels from reference
      expr <- full_join(expr, private$.ref[,c(ALIQUOT_LABEL, CLASS_LABEL)],
                        by = ALIQUOT_LABEL)
      # Rename rows
      expr <- expr %>% column_to_rownames(ALIQUOT_LABEL) %>% as.data.frame()
      
      # Save data
      private$.data <- expr
      
    }
    
  ),
  
  
  public = list(
    
    initialize = function(cris_data, classes){
     
      if (all(class(cris_data) != 'CRISData'))
        stop('`cris_data` must be a CRISData object')
      
      super$initialize(cris_data$data,cris_data$ref,classes)
      
      private$.prepare_data()
    },
    
    stratified_split  = function(train_perc, seed){
      
      # Call superclass to check percentage and seed.
      super$stratified_split(train_perc, seed)
      
      if (train_perc == 0){
        private$.train_ <- NULL
        private$.test_  <- private$.data
        
        # Split reference
        train_samples <- c()
        test_samples  <- rownames(private$.test_)
      }else if (train_perc == 1){
        private$.train_ <- private$.data
        private$.test_  <- NULL
        
        # Split reference
        train_samples <- rownames(private$.train_)
        test_samples  <- c()
      }else{
        # Split data
        train_indices <- createDataPartition(private$.data[ ,CLASS_LABEL], 
                                             p     = train_perc,
                                             list  = FALSE,
                                             times = 1)
        
        # Set train and test data
        private$.train_ <- private$.data[train_indices, ]
        private$.test_  <- private$.data[-train_indices, ]
        
        # Split reference
        train_samples <- rownames(private$.train_)
        test_samples  <- rownames(private$.test_)
      
      }
      
      private$.split_reference(train_samples, test_samples)
    },
    
    split_by_list  = function(train_samples){
      
      if (length(train_samples) == 0){
        private$.train_ <- NULL
        private$.test_  <- private$.data
        
        train_samples <- c()
        test_samples  <- rownames(private$.data)
      }
      else{
        if (class(train_samples) != 'character')
          stop('`train_samples` must be a character vector.')
      
        # Split data
        
        train_indices <- which(rownames(private$.data) %in% train_samples)
        if (length(train_indices) > 0){
          private$.train_ <- private$.data[train_indices, ]
          private$.test_  <- private$.data[-train_indices, ]
          
          # Split reference
          train_samples <- rownames(private$.train_)
          test_samples  <- rownames(private$.test_)
        }else{
          private$.train_ <- NULL
          private$.test_  <- private$.data
          
          # Split reference
          train_samples <- c()
          test_samples  <- rownames(private$.test_)
        }
        
      }  
    
      private$.split_reference(train_samples, test_samples)
      
    }
    
  )
  
  
)






