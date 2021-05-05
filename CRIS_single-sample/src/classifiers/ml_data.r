# Description -------------------------------------------------------------

# Class definition for multi-label data. It extends the Data class but includes
# also the thresholds for assigning multi-label reference (from NTP single-label 
# to NTP multi-label) and a copy of the original data as data.frame (used for
# creating mldr object). 

# Class definition --------------------------------------------------------

PercentileMLData <- R6Class(
  
  'PercentileMLData',
  inherit      = Data,
  lock_objects = TRUE,
  lock_class   = TRUE,

  ##### Private fields #####
  
  private = list(
    
    .class_thr = NULL,
    .data_as_df = NULL,
    
    #' Prepare the reference saving a copy of original ntp scores (attributes: ntpCRIS-A, ntpCRIS-B, ...)
    #' and renaming the scores to be binarized with class name.
    .init_ref = function(){
      
      # Select aliquot ID, primary NTP class, ntp scores (to be binarized), ntp fdr scores
      private$.ref <- private$.ref %>% 
        select(ALIQUOT_LABEL, CLASS_LABEL, private$.classes, CLASS_FDR_LABEL) %>%
        as.data.frame()
      
      private$.ref[ ,CLASS_LABEL] <- 
        factor(private$.ref[ ,CLASS_LABEL], levels = private$.classes)
      
      # Save a copy of the original ntp scores in reference
      scores           <- private$.ref[,c(ALIQUOT_LABEL,private$.classes)]
      colnames(scores) <- c(ALIQUOT_LABEL, paste('ntp', private$.classes, sep = ''))
      private$.ref <- full_join(scores, private$.ref, by = ALIQUOT_LABEL)
      
      # Rownames of reference: aliquot ids
      rownames(private$.ref) <- private$.ref[,ALIQUOT_LABEL]
    },
    
    #' Binarize the NTP class scores. For each class, in each sample, assign 1
    #' if fdr < 0.2 and class = CRIS-class OR if fdr < 0.2 and score >= class-threshold
    .prepare_ref = function(){
      
      private$.init_ref()
      
      # For each class, in each sample assign 1
      # if fdr < 0.2 and class = CRIS-class OR 
      # if fdr < 0.2 and score >= class-threshold
      for (c in seq(length(private$.classes))) {
        cl <- private$.classes[c]
        
        for (s in seq(nrow(private$.ref))) {
          cl_score <- private$.ref[s,cl]
          cl_fdr   <- private$.ref[s,CLASS_FDR_LABEL[c]]
          cl_ntp   <- private$.ref[s,CLASS_LABEL]
          
          if (cl_fdr < BH_FDR_THRESHOLD) {
            
            if (cl_ntp == cl)
               private$.ref[s,cl] <- 1
            else if (cl_ntp != cl & cl_score >= private$.class_thr[c])
               private$.ref[s,cl] <- 1
            else
               private$.ref[s,cl] <- 0
            
          }else{
            private$.ref[s,cl] <- 0
          }
            
        }
        
      }
      
      # Save as data.frame
      private$.ref <- as.data.frame(private$.ref)
      
    },
    
    
    #' Prepare data_as_df object with expression data and binarized scores, obtained
    #' from prepared reference.
    #' 
    #' @return set the data and data_as_df fields with mldr object and data.frame
    #' object representing data and binarized multi-label target, respectively.
    .prepare_data = function(){
      
       # Expression matrix (samples x genes)
       expr   <- as.data.frame(t(private$.data@assayData$exprs))

       # Join with reference to add binarized scores of each class
       expr   <- expr %>% 
         rownames_to_column(ALIQUOT_LABEL) %>% 
         as.data.frame()
       
       private$.data_as_df <- 
         full_join(expr, private$.ref[ ,c(ALIQUOT_LABEL, CRIS_CLASSES)], by = ALIQUOT_LABEL)

       # Rename rows with aliquot IDs
       private$.data_as_df <- private$.data_as_df %>% column_to_rownames(ALIQUOT_LABEL)
       
       # Get indices of columns with labels
       lab_cols <- which(colnames(private$.data_as_df) %in% private$.classes)

       # Create mldr from dataframe
       private$.data <-
         mldr_from_dataframe(private$.data_as_df, labelIndices = lab_cols, name = "ml_data")
       
       
    }
    
    
  ),
  
  public = list(
    
    #' Constructor
    #' 
    #' @param cris_data A CRISData object from which data and reference are extracted
    #' @param classes   The class names character vector
    #' @param class_thr The distances used to obtain multi-label NTP reference
    #' from NTP single-label reference
    initialize = function(cris_data, classes, class_thr){
     
      # Check input
      if (names(class_thr) != CRIS_CLASSES | mode(class_thr) != 'numeric')
        stop('PercentileMLData: `class_thr` must contain a nuemric threshold for each class')
      
      if (all(class(cris_data) != 'CRISData'))
        stop('PercentileMLData: `cris_data` must be a CRISData object')
      
      # Assign input
      private$.class_thr <- class_thr 
      super$initialize(cris_data$data, cris_data$ref, classes)
      
      # Adapt ML ref (NTP correlations >> binary class assignment)
      private$.prepare_ref()

      # Convert data into MLDR object
      private$.prepare_data()

    },
    
   
    #' Split the data and reference putting into training the provided samples.
    #' The remainder is put into testing.
    #' 
    #' @param train_samples A character vector with the training samples, provided
    #' through aliquot IDs.
    #' 
    #' @return Set the train_ and test_ fields, together with corresponding
    #' train_ref and test_ref of the current object.
    split_by_list  = function(train_samples){
      
      # If no train sample, put all data into testing
      if (length(train_samples) == 0){
        print_info('all testing')
        private$.train_ <- NULL
        private$.test_  <- private$.data
        
        train_samples <- c()
        test_samples  <- rownames(private$.test_$dataset)
      }
      else{
        if (class(train_samples) != 'character')
          stop('`train_samples` must be a character vector.')
      
        # Samples belonging to required train set
        train_indices <- which(rownames(private$.data_as_df) %in% train_samples)
        
        # At least one of required samples is present in the dataset
        if (length(train_indices) > 0){
          
          # Get indices of columns named by class labels (NTP scores)
          lab_cols <- which(colnames(private$.data_as_df) %in% private$.classes)
      
          # Create MLDR objects for train and test data
          private$.train_ <- private$.data_as_df[train_indices, ] %>% 
            mldr_from_dataframe(labelIndices = lab_cols, name = "train_")
          private$.test_  <- private$.data_as_df[-train_indices, ] %>% 
            mldr_from_dataframe(labelIndices = lab_cols, name = "test_")
          
          # Set samples for training and testing
          train_samples <- rownames(private$.train_$dataset)
          test_samples  <- rownames(private$.test_$dataset)
        
        } # Some training samples have been required, but none is present in the data
        else{
          
          private$.train_ <- NULL
          private$.test_  <- private$.data
          
          # Set samples for training and testing
          train_samples <- c()
          test_samples  <- private$.test_$dataset$aliquot_id
        }
        
      }  
      
      # Split the reference with the train and test samples set before
      private$.split_reference(train_samples, test_samples)
    
    }
    
  ),
  
  ##### Active fields #####
  # Read-only access to private fields; the missing ones are inherited from Data.
  
  active = list(
    class_thr = function(v){
      if (missing(v))
        private$.class_thr
      else
        stop('PercentileMLData: class_thr is read-only')
        
    }
  )
)
