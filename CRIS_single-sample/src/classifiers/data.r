
# Description -------------------------------------------------------------

# Superclass of classifier data objects. Contains data, the reference and
# the list of classes. Moreover, it contains train and test data, together
# with reference splitted in train and test. Among the settings, the
# percentage of training data is saved, together with the seed used for
# random stratified split.

# Libraries ---------------------------------------------------------------

# Class
library(R6)

# Expression Set
library(Biobase)

# Data manipulation and %>% pipe
library(tidyverse)



# Class definition --------------------------------------------------------

Data <- R6Class(
  
  'Data',
  lock_objects = TRUE,
  lock_class   = TRUE,

  private = list(
    
    # Global
    .data = NULL,
    .ref  = NULL,
    .classes   = NULL,
    
    # Splitted
    .train_     = NULL,
    .test_      = NULL,
    .train_ref  = NULL,
    .test_ref   = NULL,
    
    # Settings
    .perc       = NULL,
    .seed       = NULL,
    
    
    #' Split reference into training an testing
    #'
    #' @param train_samples Character vector with training samples (ids)
    #' @param test_samples  Character vector with testing samples (ids)
    .split_reference   = function(train_samples, test_samples){
      
      # Check input
      
      if (class(private$.ref) != 'data.frame')
        stop('`ref` can be splitted only if it is a data.frame')

      if (class(train_samples) != 'character' & length(train_samples) > 0)
        stop('`train_samples` must be a character vector.')
      
      if (class(test_samples) != 'character' & length(test_samples) > 0)
        stop('`test_samples` must be a character vector.')
      
      # Training and testing must be disjoint
      
      train_test <- intersect(train_samples, test_samples)
      if (length(train_test) > 0)
        stop('`test_samples` and `train_samples` must be disjoint')
      
      # Filtering
      
      f <- Filter$new()

      private$.train_ref <-
        f$filter_by_attribute(private$.ref, ALIQUOT_LABEL, train_samples) 
      # rownames(private$.train_ref) <- NULL
      # private$.train_ref <- private$.train_ref %>% column_to_rownames(ALIQUOT_LABEL)

      private$.test_ref <-
        f$filter_by_attribute(private$.ref, ALIQUOT_LABEL, test_samples) 
      # rownames(private$.test_ref) <- NULL
      # private$.test_ref <- private$.test_ref %>% column_to_rownames(ALIQUOT_LABEL)
      
    }
    
    
  ),
  
  public = list(
    
    
    #' Constructor
    #'
    #' @param data    ExpressionSet with data
    #' @param ref     Data.frame with aliquot ids, class label, class scores
    #'                (column for each class, with class name)
    #' @param classes List of classes
    #'
    #' @return        An object of class Data
    initialize = function(data, ref, classes){
     
      # Check input 
      
      if (class(data) != 'ExpressionSet')
        stop('`data` must be an ExpressionSet')
      
      if (class(ref) != 'data.frame')
        stop('`ref` must be a data.frame')
      
      if (class(classes) != 'character' | length(classes) < 2)
        stop('`classes` must be a character vector with at least two classes.')
      
      # Reference fields
      
      if (!all(c(ALIQUOT_LABEL, CLASS_LABEL, classes) %in% colnames(ref)))
        stop('`ref` must contain aliquot_id, class label and classes scores.')
      
      # Samples correspondence
      
      data_samples <- sort(colnames(data@assayData$exprs))
      ref_samples  <- sort(ref[,ALIQUOT_LABEL])
      
      if (!all.equal(data_samples,ref_samples))
        stop('data and reference samples do not correspond.')
      
      # Assign fields (all controls satified)
      private$.data    <- data
      private$.ref     <- ref
      private$.classes <- classes
      
    },
    
    
    
    #' Create split for training and testing with required percentage
    #' of training samples. The seed is used for reproducibility
    #'
    #' @param train_perc Percentage of samples to put into training (in[0,1])
    #' @param seed       Integer seed for reproducibility of sampling
    #'
    #' @return           Train and test data are saved into the object.
    stratified_split  = function(train_perc, seed){
      
      # Check input
      
      if (class(train_perc) != 'numeric' | class(seed) != 'numeric')
        stop('`train_perc` and `seed` must be numbers')
      
      if (train_perc < 0 | train_perc > 1)
        stop('`train_perc` must be between 0 and 1 (included)')
      
      if (round(seed) != seed)
        warning('Rounding the seed to integer.')
      
      # Assign and save parameters
      
      private$.perc <- train_perc
      private$.seed <- seed 
      
      # Set seed
      set.seed(private$.seed)

      # Splitting implemented in subclasses
      print('Split')
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
    },
    
    # Splitted
    
    train_    = function(v){
      if (missing(v))
        private$.train_
      else
        stop('`train_` is read-only.')
    },
    
    test_     = function(v){
      if (missing(v))
        private$.test_
      else
        stop('`test_` is read-only.')
    },
    
    train_ref = function(v){
      if (missing(v))
        private$.train_ref
      else
        stop('`train_ref` is read-only.')
    },
    
    test_ref  = function(v){
      if (missing(v))
        private$.test_ref
      else
        stop('`test_ref` is read-only.')
    },
    
    # Settings
    
    perc  = function(v){
      if (missing(v))
        private$.perc
      else
        stop('`perc` is read-only.')
    },
    
    seed  = function(v){
      if (missing(v))
        private$.seed
      else
        stop('`seed` is read-only.')
    }
    
    
  )
  
)