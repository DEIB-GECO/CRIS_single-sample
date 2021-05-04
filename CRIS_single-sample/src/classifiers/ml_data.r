
PercentileMLData <- R6Class(
  
  'PercentileMLData',
  inherit      = Data,
  lock_objects = TRUE,
  lock_class   = TRUE,

  private = list(
    
    .class_thr = NULL,
    .data_as_df = NULL,
    
    .init_ref = function(){
      
      # Select aliquot, primary NTP class, ntp scores (to be binarized), ntp fdr
      private$.ref <- private$.ref %>% 
        select(ALIQUOT_LABEL, CLASS_LABEL, private$.classes, CLASS_FDR_LABEL) %>%
        as.data.frame()
      
      private$.ref[ ,CLASS_LABEL] <- factor(private$.ref[ ,CLASS_LABEL], levels = private$.classes)
      
      # Save original ntp scores in reference
      scores           <- private$.ref[,c(ALIQUOT_LABEL,private$.classes)]
      colnames(scores) <- c(ALIQUOT_LABEL, paste('ntp', private$.classes, sep = ''))
      private$.ref <- full_join(scores, private$.ref, by = ALIQUOT_LABEL)
      
      # Rownames of reference: aliquot ids
      rownames(private$.ref) <- private$.ref[,ALIQUOT_LABEL]
    },
    
    .add_ranks = function(){
      
      classes <- levels(F_CRIS_CLASSES)
      rank_labels     <- paste('rank', classes, sep = '')
      w_scores_labels <- paste('weighted', classes, sep = '')
      
      # Empty weighted scores
      w_scores <- matrix(0, nrow = nrow(private$.ref), ncol = length(classes))
      colnames(w_scores) <- w_scores_labels
      
      # Empty ranks
      ranks <- matrix(0, nrow = nrow(private$.ref), ncol = length(classes))
      colnames(ranks) <- rank_labels
      
      # Add weighted scores and ranks to result
      private$.ref     <- cbind(private$.ref, as.data.frame(w_scores))
      private$.ref     <- cbind(private$.ref, as.data.frame(ranks))
      
 
      # Fill the ranks when at least one class is assigned
      for (k in seq(nrow(private$.ref))){
        
        # Get and sort (descending) scores for current sample
        scores        <- private$.ref[k, paste('ntp',classes, sep = '')] %>% 
                         unlist() %>% 
                         as.numeric() %>%
                         signif(8)
        
        # Scores/thresholds 
        weighted_scores <- scores/private$.class_thr[classes]
        sorted_w_scores <- sort(unique(weighted_scores), decreasing = TRUE)
        
        # Get position (rank) of each score in the sorted list of unique scores
        ranking <- lapply(X   = seq(length(weighted_scores)), 
                          FUN = function(k){which(sorted_w_scores == weighted_scores[k])}) 
        
        # Assign weighted scores and ranking
        private$.ref[k,w_scores_labels] <- weighted_scores
        private$.ref[k,rank_labels]     <- ranking
        
      }
      
    },
    
    .prepare_ref = function(){
      
      private$.init_ref()
      
      # For each class, assign 1
      # if fdr < 0.2 and class = CRIS-class or 
      # if fdr < 0.2 and score > class-threshold
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
      
      # Add ranks
      private$.add_ranks()
    },
    
    #' Add one-hot encoding to data and create an mldr object
    .prepare_data = function(){
      
       # Expression matrix (samples x genes)
       expr   <- as.data.frame(t(private$.data@assayData$exprs))

       # Join with reference to get binary vector of classes
       expr   <- expr %>% rownames_to_column(ALIQUOT_LABEL) %>% as.data.frame()
       private$.data_as_df <- full_join(expr, private$.ref[ ,c(ALIQUOT_LABEL, levels(F_CRIS_CLASSES))], by = ALIQUOT_LABEL)

       # Rename rows
       private$.data_as_df <- private$.data_as_df %>% column_to_rownames(ALIQUOT_LABEL)
       
       # Get indices of columns with labels
       lab_cols <- which(colnames(private$.data_as_df) %in% private$.classes)

       # Create mldr from dataframe
       private$.data <-
         mldr_from_dataframe(private$.data_as_df, labelIndices = lab_cols, name = "ml_data")
       
       
    }
    
    
  ),
  
  public = list(
    
    # max_labels is ignored here
    initialize = function(cris_data, classes, class_thr){
     
      private$.class_thr <- class_thr 
      
      if (all(class(cris_data) != 'CRISData'))
        stop('`cris_data` must be a CRISData object')
      
      super$initialize(cris_data$data, cris_data$ref, classes)
      
      # Adapt ML ref (NTP correlations >> binary class assignment)
      private$.prepare_ref()

      # Convert data into MLDR object
      private$.prepare_data()

    },
    
    stratified_split  = function(train_perc, seed){
      
      super$stratified_split(train_perc, seed)
      if (train_perc == 0){
        private$.train_ <- NULL
        private$.test_  <- private$.data
        
        # Split reference
        train_samples <- c()
        test_samples  <- rownames(private$.test_$dataset)
      }else if (train_perc == 1){
        private$.train_ <- private$.data
        private$.test_  <- NULL
        
        # Split reference
        train_samples <- rownames(private$.train_$dataset)
        test_samples  <- c()
      }else{
        
        # Split data
        perc <- c(train = train_perc, test = 1 - train_perc)
        ds   <- create_holdout_partition(private$.data, perc, method = 'stratified')
      
        # Set train and test data
        private$.train_ <- ds$train
        private$.test_  <- ds$test
        
        # Split reference
        train_samples <- rownames(ds$train$dataset)
        test_samples  <- rownames(ds$test$dataset)
      
      }
      
      private$.split_reference(train_samples, test_samples)
      
    },
    
    split_by_list  = function(train_samples){
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
      
        # Split data
        train_indices <- which(rownames(private$.data_as_df) %in% train_samples)
        
        if (length(train_indices) > 0){
          # Get indices of columns with labels
          lab_cols <- which(colnames(private$.data_as_df) %in% private$.classes)
      
          private$.train_ <- private$.data_as_df[train_indices, ] %>% 
            mldr_from_dataframe(labelIndices = lab_cols, name = "train_")
          private$.test_  <- private$.data_as_df[-train_indices, ] %>% 
            mldr_from_dataframe(labelIndices = lab_cols, name = "test_")
          
          # Split reference
          train_samples <- rownames(private$.train_$dataset)
          test_samples  <- rownames(private$.test_$dataset)
        }else{
          private$.train_ <- NULL
          private$.test_  <- private$.data
          
          # Split reference
          train_samples <- c()
          test_samples  <- private$.test_$dataset$aliquot_id
        }
        
      }  
      
      private$.split_reference(train_samples, test_samples)
    
    }
    
  ),
  
  active = list(
    class_thr = function(v){
      if (missing(v))
        private$.class_thr
      else
        stop('class_thr is read-only')
        
    }
  )
)
