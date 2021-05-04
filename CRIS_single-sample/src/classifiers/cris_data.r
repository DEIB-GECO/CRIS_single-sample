
# Description -------------------------------------------------------------

# Class for replication of CRIS classifiers (NTP and TSP). Contains the data, 
# the reference and the list of classes. Contains data, the reference and
# the list of classes. Moreover, it contains train and test data, together
# with reference splitted in train and test. Among the settings, the
# percentage of training data is saved, together with the seed used for
# random stratified split.
# 
# Finally, it contains the original scores and the max number of labels
# allowed (if NULL, no maximum is set).

# Libraries ---------------------------------------------------------------

# Class
library(R6)

# Expression Set
library(Biobase)

# Data manipulation and %>% pipe
library(tidyverse)


# Class definition --------------------------------------------------------
CRISData <- R6Class(
  
  'CRISData',
  lock_objects = TRUE,
  lock_class   = TRUE,
  
  #### Public fields ####
  
  private = list(
    
    .data = NULL,
    .ref  = NULL,
    .classes = NULL,
    
    .validate_data = function(data){
      
      if(!check_type(data,'ExpressionSet'))
        stop('use ExpressionSet data')
      
      if(!check_type(data@assayData$exprs,'matrix',c(1,1)) | !check_type(data@assayData$exprs,'numeric'))
        stop('Invalid expression matrix')

      
    },
    
    .validate_ref = function(ref){
      
      # Check type
      if (!check_type(ref,'data.frame'))
        stop('use data.frame as reference')
      
      # Check samples
      if (!any(all.equal(rownames(ref), ref[,ALIQUOT_LABEL]) == TRUE))
        stop(paste('reference must have non-null rownames (', ALIQUOT_LABEL, 'column)'))
      
      # Check if columns are missing
      required_columns <- c(ALIQUOT_LABEL,
                            CLASS_LABEL, 
                            CLASS_FDR_LABEL,
                            CLASS_DISTANCE_LABEL)
      
      if (!any(all.equal(sort(required_columns), sort(colnames(ref))) == TRUE))
        stop(paste(c('All and only required columns:', required_columns), collapse = '\n'))

      if (!check_type(ref[,ALIQUOT_LABEL],'character'))
        stop(paste(ALIQUOT_LABEL, 'must be character'))
      
      # Check class label column is a factor with CRIS levels
      if (!check_type(ref[,CLASS_LABEL],'factor'))
        stop(paste(CLASS_LABEL, 'must be a factor'))
      
      if (!any(all.equal(levels(ref[,CLASS_LABEL]), private$.classes)))
        stop(paste(CLASS_LABEL, 'must be a factor with', private$.classes, 'levels.'))
      
      # Check other columns are numeric
      numeric_columns <- c(CLASS_FDR_LABEL, CLASS_DISTANCE_LABEL)
      dist_fdr_modes <- lapply(numeric_columns, function(k){
        return(check_type(ref[,k], 'numeric'))
      }) %>% unlist()
      
      if (!all(dist_fdr_modes))
        stop(paste(numeric_columns, 'must be numeric'))
      
      # Check missing values
      if (any_uncomplete(ref))
        stop('NA, NULL, NaN or Inf values are not allowed')
      
    },
    
    #' Prepare reference by translating NTP class distances into correlations
    .prepare_ref = function(){
      
      # Extract class names in from distance columns
      colnames(private$.ref) <- colnames(private$.ref) %>%
        gsub(pattern = 'distCRIS.', replacement = 'CRIS-')
      
      # Correlation = 1 - distance (i.e. class score)
      private$.ref[,private$.classes] <- 1 - private$.ref[,private$.classes]
      private$.ref <- as.data.frame(private$.ref)

    }
  ),
  
  #### Public fields ####
  
  public = list(
    
    initialize = function(data, ref = NULL){
     
      # Set classes
      private$.classes <- levels(F_CRIS_CLASSES)
      
      # Check and assign data expressionset
      private$.validate_data(data)
      private$.data    <- data
      
      # Check and assign reference
      if (check_type(ref,'data.frame')) {
        
        private$.validate_ref(ref)
        
        # Check samples correspondence
        data_samples <- sort(colnames(private$.data@assayData$exprs))
        ref_samples  <- sort(rownames(ref))
        
        if (!any(all.equal(data_samples, ref_samples) == TRUE))
          stop('Not NULL reference with no samples correspondence')
        
        # If reference is valid, assign and prepare it
        private$.ref <- ref
        private$.prepare_ref()
      }
      
        
    }
    
  ),
  
  #### Active fields ####
  # Read only access to private fields
  
  active = list(
    
    # Global
    
    data    = function(v){
      if (missing(v))
        private$.data
      else
        stop('`data` is read-only.')
    },
    
    ref     = function(v){
      if (missing(v))
        private$.ref
      else
        stop('`ref` is read-only.')
    },
    
    classes = function(v){
      if (missing(v))
        private$.classes
      else
        stop('`classes` is read-only.')
    }
  )
  
  
)






