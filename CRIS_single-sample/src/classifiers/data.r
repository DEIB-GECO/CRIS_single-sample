
# Description -------------------------------------------------------------

# Superclass of classifier data objects (sl data and ml data). 
# Contains data, the reference and the list of classes. 
# Moreover, it contains train and test data, together
# with reference splitted in train and test. Among the settings, the
# percentage of training data is saved, together with the seed used for
# random stratified split.

# Class definition --------------------------------------------------------

Data <- R6Class(
  
  'Data',
  lock_objects = TRUE,
  lock_class   = TRUE,

  ##### Private fields #####
  
  private = list(
    
    # Global (original data and reference + class names vector)
    .data = NULL,
    .ref  = NULL,
    .classes   = NULL,
    
    # Splitted data and reference
    .train_     = NULL,
    .test_      = NULL,
    .train_ref  = NULL,
    .test_ref   = NULL,
    
    
    #' Split reference into training an testing
    #'
    #' @param train_samples Character vector with training samples (ids)
    #' @param test_samples  Character vector with testing samples (ids)
    .split_reference   = function(train_samples, test_samples){
      
      # Check input
      
      if (class(train_samples) != 'character' & length(train_samples) > 0)
        stop('`train_samples` must be a character vector.')
      
      if (class(test_samples) != 'character' & length(test_samples) > 0)
        stop('`test_samples` must be a character vector.')
      
      # Training and testing must be disjoint
      
      train_test <- intersect(train_samples, test_samples)
      if (length(train_test) > 0)
        stop('`test_samples` and `train_samples` must be disjoint')
      
      # Split extracting train and test samples
      
      f <- Filter$new()

      private$.train_ref <-
        f$filter_by_attribute(private$.ref, ALIQUOT_LABEL, train_samples) 

      private$.test_ref <-
        f$filter_by_attribute(private$.ref, ALIQUOT_LABEL, test_samples) 

    }
    
    
  ),
  
  ##### Public fields #####
  
  public = list(
    
    
    #' Constructor
    #'
    #' @param data    ExpressionSet with data
    #' @param ref     Data.frame with aliquot ids, class label, class scores
    #'                (column for each class, with class name)
    #' @param classes List of classes (character vector)
    initialize = function(data, ref, classes){
     
      # Check input 
      
      if (class(data) != 'ExpressionSet')
        stop('Data: `data` must be an ExpressionSet')
      
      if (class(ref) != 'data.frame')
        stop('Data: `ref` must be a data.frame')
      
      if (class(classes) != 'character' | length(classes) < 2)
        stop('Data: `classes` must be a character vector with at least two classes.')
      
      # Reference fields
      
      if (!all(c(ALIQUOT_LABEL, CLASS_LABEL, classes) %in% colnames(ref)))
        stop('Data: `ref` must contain aliquot_id, class label and classes scores.')
      
      # Samples correspondence
      
      data_samples <- sort(colnames(data@assayData$exprs))
      ref_samples  <- sort(ref[,ALIQUOT_LABEL])
      
      if (!all.equal(data_samples,ref_samples))
        stop('Data: data and reference samples do not correspond.')
      
      # Assign fields (all controls satisfied)
      private$.data    <- data
      private$.ref     <- ref
      private$.classes <- classes
      
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
        stop('Data: `data` is read-only.')
    },
    
    ref     = function(v){
      if (missing(v))
        private$.ref
      else
        stop('Data: `ref` is read-only.')
    },
    
    classes = function(v){
      if (missing(v))
        private$.classes
      else
        stop('Data: `classes` is read-only.')
    },
    
    # Splitted
    
    train_    = function(v){
      if (missing(v))
        private$.train_
      else
        stop('Data: `train_` is read-only.')
    },
    
    test_     = function(v){
      if (missing(v))
        private$.test_
      else
        stop('Data: `test_` is read-only.')
    },
    
    train_ref = function(v){
      if (missing(v))
        private$.train_ref
      else
        stop('Data: `train_ref` is read-only.')
    },
    
    test_ref  = function(v){
      if (missing(v))
        private$.test_ref
      else
        stop('Data: `test_ref` is read-only.')
    }
    
  )
  
)