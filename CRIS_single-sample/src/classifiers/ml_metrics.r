# Description -------------------------------------------------------------

# Class for metrics of multi-label classifiers.

# Class definition --------------------------------------------------------

MLMetrics <- R6Class(
  'MLMetrics',
  inherit = Metrics,
  lock_objects = FALSE,
  lock_class = TRUE,
  
  private = list(
    
    
    #' Computes a confusion matrix with respect to a specific class (used for additional metrics,
    #' e.g. mcc score)
    #' 
    #' @param ml_conf Confusion matrix table obtained from utiml library
    #' @param class_name  name of the class with respect to which computing 
    #' true positives (tp), false positives (fp), false negatives (fn), true negatives (tn).
    #' 
    #' @return list with tp, tn, fp, fn computed with respect to the required class.
    .class_conf_mat     = function(ml_conf, class_name){
      
      if (!is.character(class_name))
        stop(paste('class name must be a character') )
      
      if (any(!class_name %in% CRIS_CLASSES))
        stop(paste('Cannot found the class',class_name,'for class confusion matrix'))
      
      conf <- list(
        tp = ml_conf$TPl[[class_name]],
        tn = ml_conf$TNl[[class_name]],
        fp = ml_conf$FPl[[class_name]],
        fn = ml_conf$FNl[[class_name]]
      )
      
      
      return(conf)
    },
    
    #' Check that the samples of reference and prediction coincide. Rownames
    #' of reference and prediction are aliquot ids.
    #' 
    #' @param ref reference target. If ref is an mldr object, pred must be a mlresult object
    #' @param pred prediction. If ref is an mldr object, pred must be a mlresult object
    #' 
    #' @return stops the execution if the correspondence is not verified.
    .check_correspondence = function(ref, pred){
      
      if (class(ref) == 'mldr' & class(pred) == 'mlresult'){

        # Check class labels correspondence
        if (!any(all.equal(rownames(ref$labels), colnames(as.matrix(pred))) == TRUE))
          stop('Label names do not coincide')
        
        
        # Check samples correspondence
        if (!any(all.equal(rownames(ref$dataset), rownames(as.matrix(pred))) == TRUE))
          stop('Samples do not coincide')
        
      }else {
        
        # Check samples correspondence
        if (!any(all.equal(rownames(ref), rownames(pred)) == TRUE))
          stop('Samples do not coincide')
        
      }
      
      
      
    },
    
    # Functions to check the input
    
    .is_binary = function(obj){
      
      classes   <- CRIS_CLASSES
      is_binary <- FALSE
      
      # Convert to df
      if (class(obj) == 'mlresult'){
        obj <- as.data.frame(as.matrix(obj))
      }
      
      if (any(class(obj) %in% c('matrix','data.frame'))){
        
        # Check labels columns exist
        if (!all(classes %in% colnames(obj)))
          stop('Cannot found columns for classes.')
      
        # Check labels columns contain only 0 and 1 values
        type_check <- lapply(classes, function(k){
          return(all(unique(obj[,k]) %in% c(0,1)))
        }) %>% unlist()

        if (all(type_check == TRUE))
          is_binary <- TRUE
        
      }
      
      return(is_binary)
    },
    
    .validate_binary = function(ref, pred){
      
      # Check reference
      if (class(ref) != 'mldr'){
        stop('reference must be a mldr object')
      }
      
      if (!private$.is_binary(ref$dataset)){
        stop('reference labels must contain only 0 and 1 values')
      }
      
      # Check prediction
      if (all(class(pred) != 'mlresult'))
        stop('prediction must be an mlresult')
      
      
      if (!private$.is_binary(attr(pred, 'classes'))){
        stop('Prediction labels must contain only 0 and 1 values')
      }
      
      # Check correspondence
      private$.check_correspondence(ref, pred)
    },
    
    .validate_df = function(ref,pred){
      
      # Check reference
      if (class(ref) != 'data.frame'){
        stop('reference must be an data.frame object')
      }
      
      if (!CLASS_LABEL %in% colnames(ref)){
        stop(paste('reference must contain ', CLASS_LABEL))
      }
      
      # Check prediction
      if (class(pred) != 'data.frame')
        stop('prediction can be only a data.frame')
      
      if (!private$.is_binary(pred))
        stop('Prediction labels must contain only 0 and 1 values')
      
      # Check missing values
      if (any_uncomplete(pred))
        stop('Prediction cannot be NA, NULL, infinite or NaN')
      
      if (any_uncomplete(ref[,CLASS_LABEL]))
        stop('primary NTP class cannot be NA, NULL, infinite or NaN')
      
      # Check correspondence of samples
      private$.check_correspondence(ref,pred)
      
    },
    
    .validate_sim = function(df){
      
      classes <- CRIS_CLASSES
      rank_cols     <- paste('rank',classes, sep = '')
      weighted_cols <- paste('weighted',classes,sep = '')
      
      # Check reference
      if (class(df) != 'data.frame'){
        stop('must be an data.frame object')
      }
      
      if (!all(rank_cols %in% colnames(df)))
        stop('must contain all class ranks')
      
      if (!all(weighted_cols %in% colnames(df)))
        stop('must contain all class weighted scores')
      
      if (!all(classes %in% colnames(df)) | !private$.is_binary(df))
        stop('must contain all class binary values')
      
      # Check missing values
      if (any_uncomplete(df))
        stop('cannot contain NA, NULL, infinite or NaN')

    },
    
    # Function for the relaxed accuracy and the rank accuracies
   
    rel_acc = function(ref,pred){
      
      # List and number of samples
      samples <- rownames(pred)
      n <- nrow(pred)
      # True positives (primary NTP class has been predicted)
      tp <- lapply(samples, function(k){
        primary_class <- ref[k,CLASS_LABEL]
        return(pred[k, primary_class] == 1)
      })
      
      # Count the true positives
      tp <- unlist(tp) %>% sum()
      
      # Compute the relaxed accuracy
      if (min(tp,n) == 0)
        return(0)
      else
        return(tp/n)
      
    },
    
    rank_accuracy = function(ref,pred){
      
      # List of classes
      classes    <- CRIS_CLASSES
      rank_cols  <- paste('rank', classes,sep = '')
      
      # Get reference ranks of samples in result
      ref_ranks <- ref[rownames(pred), rank_cols]
      
      # Initialize value of accuracy
      accuracy <- 0
      
      # No sample: return warning
      if (nrow(pred) < 1){
        warning('No data. Returning rank accuracy of -1.')
      }else{
        
        # Count accurate samples
        for (k in seq(nrow(pred))){
          
          # Find assigned cl
          assigned        <- classes[which(pred[k,classes] == 1)]
          r_assigned_cols <- paste('rank', assigned, sep = '')
          
          # If at least one class is assigned, check ranks
          if (length(assigned) >= 1) {
            
            # Get ranks for assigned classes 
            pred_ranks_k <- pred[k,r_assigned_cols]
            ref_ranks_k  <- ref_ranks[k,r_assigned_cols]
            
            # Consider accurate if ranks are exactly the same
            if (any(all.equal(pred_ranks_k, ref_ranks_k) == TRUE))
              accuracy <- accuracy + 1
            
          }
            
        }
        
        # Compute accuracy
        accuracy <- accuracy/nrow(pred)
      }
      

      return(accuracy)
      
    },
    
    rank_subset_accuracy = function(ref,pred){
      
      classes <- CRIS_CLASSES
      
      # Get samples predicted exactly as the reference
      exact_samples <- lapply(rownames(ref), function(k){
        
        assigned       <- classes[which(pred[k,classes] == 1)]
        assigned_ref   <- classes[which(ref[k,classes] == 1)]
        
        if (any(all.equal(assigned, assigned_ref) == TRUE))
          return(k)  
        else
          return(NULL)
          
      }) %>% unlist()
      
      # Return rank accuracy computed on exact samples only
      return(private$rank_accuracy(ref[exact_samples, ], pred[exact_samples, ]))
      
    },
    
    scores_correlation = function(ref, pred, method){
      
      if (!any(method %in% c('pearson','spearman')))
        stop('Type of correlation must be either pearson or spearman')
      
      # Empty vector for correlation
      correl  <- rep(NA,nrow(pred))
      samples <- rownames(pred)
      scores_col <- paste('weighted',CRIS_CLASSES, sep = '')
      names(correl) <- samples
      
      # Compute correlation of scores for each samples
      for (s in samples){
        ref_scores  <- t(ref[s, scores_col])
        pred_scores <- t(pred[s,scores_col])
        correl[s] <- cor(ref_scores, pred_scores, method = method) %>% as.numeric()
      }
      
      # Save correlations with also average value
      correlation <- list(
        values = correl,
        avg = mean(correl, na.rm = TRUE)
      )
      
      return(correlation)
    }
    
  ),
  
  ##### Public fields #####
  
  public = list(
    
    
    #' Returns the confusion matrix of multi-label classifier using the
    #' utiml library funciton (utiml::multilabel_confusion_matrix)
    #' 
    #' @param ref the reference (mldr dataset)
    #' @param pred the predicted outcome (mlresult object)
    #' 
    #' @return the confusion matrix obtained with the utiml library function.
    get_conf_mat = function(ref, pred){
      
      private$.validate_binary(ref,pred)
      
      # Compute confusion matrix
      return(multilabel_confusion_matrix(mdata = ref, mlresult = pred))
      
    },
    
    
    #' Prints the confusion matrix as a data.frame
    #' 
    #' @param ml_conf the confusion matrix to print
    #' @return the confusion matrix in a data.frame
    print_confusion_matrix = function(ml_conf){
      
      
      df_conf <- data.frame(tp = ml_conf$TPl, 
                            fp = ml_conf$FPl,
                            tn = ml_conf$TNl,
                            fn = ml_conf$FNl) %>% t()
      
      df_conf <- as.data.frame(df_conf) %>% rownames_to_column('Quantity')
    
      return(df_conf)
    },
    
    #' Returns the global metrics that can be computed
    #' @return character vector with the computed global metrics
    list_global_metrics = function(){
      measures <- c('accuracy',
                    'subset-accuracy',
                    'hamming-loss',
                    'average-precision'
                    )
      
      return(measures)
    },
    
    
    #' Returns the local metrics that can be computed
    #' @return character vector with the computed local metrics
    list_local_metrics = function(){
      
      measures <- c('precision',
                    'recall',
                    'F1',
                    'mcc')
      
      return(measures)
    },
    
    
    
    
    #' Uses the function get_conf_mat to compute a confusion matrix
    #' for multi-label classifier and returns the global metrics obtained from it
    #' 
    #' @param ref the reference with class targets. It must contain a class label column
    #' @param pred the predicted outcome. It must contain a class label column
    #' @param sel  a subset of global metrics to extract
    #' 
    #' @return a data.frame containing the global metrics and their value.
    global_metrics = function(ref,pred, sel = NULL){

      ml_conf <- self$get_conf_mat(ref,pred)
      if (is.null(sel)){
        metrics <- multilabel_evaluate(ml_conf, measures = self$list_global_metrics())
      }else if (class(sel) != 'character' | length(sel) < 1){
        stop('Selected features must be a vector of character')
      }else{
        metrics <- multilabel_evaluate(ml_conf, measures = sel)
      }
      
      metrics <- data.frame(Value = metrics) %>% rownames_to_column('Metric')
      
      return(metrics)
    },
    
    #' Computes relaxed metrics, which refer to the original ntp class (accuracy)
    #' 
    #' @param ref the reference with class targets
    #' @param pred the predicted outcome
    #' 
    #' @return a vector containing the metrics and their value.
    relaxed_metrics = function(ref,pred){
      
      private$.validate_df(ref,pred)
      
      rel_metrics = list(
        Accuracy = private$rel_acc(ref,pred),
        Kappa    = NA
      )
      
      return(unlist(rel_metrics))
      
      
    },
    
    
    #' Computes correlation metrics
    #' TODO: not working for now, to be checked
    #' 
    #' @param ref the reference
    #' @param pred the predicted outcome
    #' 
    #' @return a list containing the metrics and their value.
    similarity_metrics = function(ref, pred){
      
      # Check input
      private$.validate_sim(ref)
      private$.validate_sim(pred)
      
      # Check samples correspondence
      private$.check_correspondence(ref,pred)

      pearson  <- private$scores_correlation(ref,pred,'pearson')
      spearman <- private$scores_correlation(ref,pred,'spearman')
      
      # Compute metrics
      sim_metr <- list(
        'rank accuracy' = private$rank_accuracy(ref,pred),
        'rank subset accuracy' = private$rank_subset_accuracy(ref, pred),
        'pearson' = pearson$avg,
        'spearman' = spearman$avg
      )
      
      # Convert into data.frame
      sim_metr <- data.frame(Value = unlist(sim_metr)) %>% rownames_to_column('Metric')
      
      return(list(
        table = sim_metr,
        pearson  = pearson$values,
        spearman = spearman$values)
      )
      
    },
    
    #' Uses the function get_conf_mat to compute a confusion matrix
    #' for multi-label classifier and returns the local metrics obtained from it.
    #' Additional metrics are added (mcc)
    #' 
    #' @param ref the reference with class targets
    #' @param pred the predicted outcome
    #' 
    #' @return a data.frame containing the global metrics and their value.
    local_metrics = function(ref,pred){
      
      ml_conf <- self$get_conf_mat(ref,pred)
      local_metrics <- multilabel_evaluate(ml_conf, labels = TRUE)
      local_metrics <- t(local_metrics$labels)
      computed_metrics <- rownames(local_metrics)
      
      
      for (m in self$list_local_metrics()) {
        
        if (!m %in% computed_metrics){
          
          # Empty row for new metric
          local_metrics   <- rbind(local_metrics, NA)
          rownames(local_metrics)[nrow(local_metrics)] <- m
          
          # Fill the table with metric value for each class
          for (c in CRIS_CLASSES) {
            args  <- list(conf = private$.class_conf_mat(ml_conf, c))
            local_metrics[m,c] <- do.call(get(m, envir = private), args)
          }
        }
        
      }
      
      local_metrics <- as.data.frame(local_metrics[self$list_local_metrics(), ])
      rownames(local_metrics) <- tolower(rownames(local_metrics))
      return(local_metrics %>% rownames_to_column('Metric'))
      
    }
    
    
    
    
  )
)
