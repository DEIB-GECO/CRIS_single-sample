# Description -------------------------------------------------------------

# Class definition for single-label data. It extends the Data class. 

# Class definition --------------------------------------------------------

SLData <- R6Class(
  
  'SLData',
  inherit      = Data,
  lock_objects = TRUE,
  lock_class   = TRUE,
  
  #### Public fields ####
  
  private = list(
    
    #' Prepare the reference saving aliquot ID, NTP class, class scores and FDRs;
    #' predicted class is transformed into a factor and aliquot ids are used as 
    #' rownames of reference.
    .prepare_ref = function(){
      
      private$.ref <- private$.ref %>% 
        select(ALIQUOT_LABEL, CLASS_LABEL, private$.classes, CLASS_FDR_LABEL) %>%
        as.data.frame()
      
      private$.ref[,CLASS_LABEL] <- 
        factor(private$.ref[,CLASS_LABEL], levels = private$.classes)
      
      rownames(private$.ref) <- private$.ref[,ALIQUOT_LABEL]
      
    },
    
    #' Get expression data and add the predicted class.
    .prepare_data = function(){
      
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
    
    #' Constructor
    #' 
    #' @param cris_data  CRISData object with data ExpressionSet and reference
    #' @param classes    The character vector with class names
    initialize = function(cris_data, classes){
     
      if (all(class(cris_data) != 'CRISData'))
        stop('`cris_data` must be a CRISData object')
      
      super$initialize(cris_data$data,cris_data$ref,classes)
      
      private$.prepare_ref()
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
      
      # If no sample is provided, put all data into testing
      if (length(train_samples) == 0){
        private$.train_ <- NULL
        private$.test_  <- private$.data
        
        # Set train and test samples for reference
        train_samples <- c()
        test_samples  <- rownames(private$.data)
      
      }# If at least one sample is provided
      else{
        
        if (class(train_samples) != 'character')
          stop('`train_samples` must be a character vector.')
      
        # Positions of training samples
        train_indices <- which(rownames(private$.data) %in% train_samples)
        
        # If at least one training sample is present
        if (length(train_indices) > 0){
          private$.train_ <- private$.data[train_indices, ]
          private$.test_  <- private$.data[-train_indices, ]
          
          # Set train and test samples for reference
          train_samples <- rownames(private$.train_)
          test_samples  <- rownames(private$.test_)
          
        }# If none of the required samples is present, put all in testing
        else{
          private$.train_ <- NULL
          private$.test_  <- private$.data
          
          # Set train and test samples for reference
          train_samples <- c()
          test_samples  <- rownames(private$.test_)
        }
        
      }  
    
      # Split the reference with training and testing samples set above
      private$.split_reference(train_samples, test_samples)
      
    }
    
  )
  
  
)






