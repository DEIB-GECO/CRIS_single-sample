
# Description -------------------------------------------------------------

# Data class used in the replication of CRIS classifiers (NTP and TSP). 
# Contains the data and the reference and the list of classes. 

# Class definition --------------------------------------------------------
CRISData <- R6Class(
  
  'CRISData',
  lock_objects = TRUE,
  lock_class   = TRUE,
  
  #### Private fields ####
  
  private = list(
    
    .data = NULL,
    .ref  = NULL,
    .classes = NULL,
    
    #' Checks the type of data is an ExpressionSet with not-emtpy numeric
    #' assayData field. In case of errors, it stops the execution of the program.
    #' 
    #' @param data   The data object to check
    .validate_data = function(data){
      
      if(!check_type(data,'ExpressionSet'))
        stop('CRIData: use ExpressionSet data')
      
      if(!check_type(data@assayData$exprs,'matrix',c(1,1)) | !check_type(data@assayData$exprs,'numeric'))
        stop('CRISData: Invalid expression matrix')

      
    },
    
    #' Checks the type of reference is data.frame with aliquot id (character), predicted class
    #' (factor), distance and fdr for all CRIS classes (numeric). Rownames must
    #' be equal to aliquot IDs. In case of errors, it stops the execution of the program. 
    #' 
    #' @param ref   The reference object to check
    .validate_ref = function(ref){
      
      # Check type
      if (!check_type(ref,'data.frame'))
        stop('CRISData: use data.frame as reference')
      
      # Check samples
      if (!any(all.equal(rownames(ref), ref[,ALIQUOT_LABEL]) == TRUE))
        stop(paste('CRISData: reference must have non-null rownames (', ALIQUOT_LABEL, 'column)'))
      
      # Check if columns are missing
      required_columns <- c(ALIQUOT_LABEL,
                            CLASS_LABEL, 
                            CLASS_FDR_LABEL,
                            CLASS_DISTANCE_LABEL)
      
      if (!any(all.equal(sort(required_columns), sort(colnames(ref))) == TRUE))
        stop(paste(c('CRISData: All and only required columns:', required_columns), collapse = '\n'))

      if (!check_type(ref[,ALIQUOT_LABEL],'character'))
        stop(paste('CRISData:', ALIQUOT_LABEL, 'must be character'))
      
      # Check class label column is a factor with CRIS levels
      if (!check_type(ref[,CLASS_LABEL],'factor'))
        stop(paste('CRISData:',CLASS_LABEL, 'must be a factor'))
      
      if (!any(all.equal(levels(ref[,CLASS_LABEL]), private$.classes)))
        stop(paste('CRISData:',CLASS_LABEL, 'must be a factor with', private$.classes, 'levels.'))
      
      # Check other columns are numeric
      numeric_columns <- c(CLASS_FDR_LABEL, CLASS_DISTANCE_LABEL)
      dist_fdr_modes <- lapply(numeric_columns, function(k){
        return(check_type(ref[,k], 'numeric'))
      }) %>% unlist()
      
      if (!all(dist_fdr_modes))
        stop(paste('CRISData:',numeric_columns, 'must be numeric'))
      
      # Check missing values
      if (any_uncomplete(ref))
        stop('CRISData:NA, NULL, NaN or Inf values are not allowed')
      
    },
    
    #' Prepare reference by translating NTP class distances into correlations
    #' (correlation = 1 - distance)
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
    
    #' Create the CRISData object; it receives an ExpressionSet (data) and possibly
    #' a reference, obtained by NTP classification (ref). 
    initialize = function(data, ref = NULL){
     
      # Set classes
      private$.classes <- CRIS_CLASSES
      
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
          stop('CRISData: Not NULL reference with no samples correspondence')
        
        # If reference is valid, assign and prepare it
        private$.ref <- ref
        private$.prepare_ref()
      }
      
        
    }
    
  ),
  
  #### Active fields ####
  # Provide read-only access to private fields
  
  active = list(
    
    # Global
    
    data    = function(v){
      if (missing(v))
        private$.data
      else
        stop('CRISData: `data` is read-only.')
    },
    
    ref     = function(v){
      if (missing(v))
        private$.ref
      else
        stop('CRISData: `ref` is read-only.')
    },
    
    classes = function(v){
      if (missing(v))
        private$.classes
      else
        stop('CRISData: `classes` is read-only.')
    }
  )
  
  
)






