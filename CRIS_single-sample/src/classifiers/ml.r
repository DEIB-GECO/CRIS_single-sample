# Description -------------------------------------------------------------

# Class for the execution of single-label classifiers.

# Class definition --------------------------------------------------------

library(matrixStats)

MLClassifier <- R6Class(
  'MLClassifier',
  inherit  = SingleSampleClassifier,
  lock_objects = TRUE,
  lock_class   = TRUE,
  
  private = list(
    .cv_set = NULL,
    
    .exec_cv    = function(data, alg_settings, tune_vals){

      set.seed(private$.seed)

      # Partition
      kfcv <- create_kfold_partition(mdata  = data,
                                     k      = private$.cv_set$folds,
                                     method = private$.cv_set$sampling)

      # Results of each fold
      result <- lapply(seq(private$.cv_set$folds), function(k) {
         
         print_info(paste('Fold',k))
         # Divide partitions in validation (k-th fold) and training
         ds <- partition_fold(kfcv, k)

         # Merge needed arguments into a unique vector
         args <- merge_lists(seed  = private$.seed,
                             base.algorithm = alg_settings$alg,
                             alg_settings$other_params,
                             tune_vals)
         
         args[['mdata']] <- ds$train
         
         # Manager of utiml metrics
         um   <- MLMetrics$new()

         # Train on all folds (except k-th) and then predict on k-th fold
         model <- do.call(alg_settings$strat, args)
         pred  <- predict(model, ds$test)
         
         # Binarize with fixed 0.5 threshold
         binary_pred <- fixed_threshold(pred,threshold = 0.5)
         attr(binary_pred, 'classes') <- as.matrix(binary_pred)[,CRIS_CLASSES]
         # Compute metrics for 
         metr <- um$global_metrics(ds$test, binary_pred, private$.cv_set$measures)
         
         # Vectorize obtained metrics
         metr_sel        <- metr$Value %>% as.numeric()
         names(metr_sel) <- metr$Metric
         return(metr_sel)
      })
      result_cv <<- result
      # Put cv result in data.frame with parameters on rows and a column for each fold
      res_df <- as.data.frame(result)
      colnames(res_df) <- NULL
      
      # Clear the seed; TODO REFERENCE
      rm(.Random.seed, envir = globalenv())
      
      # Compute mean and sd of each metric
      # https://stackoverflow.com/questions/12861734/calculating-standard-deviation-of-each-row
    
      means  <- rowMeans(res_df)
      sds    <- rowSds(as.matrix(res_df))
      res_df <- cbind(mean = means,sd = sds, res_df)
      
      # Create a numeric vector with mean and sd of each metric
      cv_metr_vect <- numeric(0)
      for (m in rownames(res_df)){
        cv_metr_vect[paste('mean', m,sep = '_')] <- as.numeric(res_df[m,'mean'])
        cv_metr_vect[paste('sd', m,sep = '_')]   <- as.numeric(res_df[m,'sd'])
      }
      
      # Return result
      print_info('Ended CV')
      return(cv_metr_vect)
    },
   
    # Used for ROC
    .get_cl_ref = function(ref,cl){
      return(private$.get_cl_ref_ml(ref,cl))
    }
    
  ),
  
  public = list(
    
    initialize  = function(cv_set, seed, cl_thresholds = NULL){
      
      # Assign seed and class thresholds
      super$initialize(seed, cl_thresholds)
      
      # Check cv settings
      if (all(class(cv_set) != 'CVSettings'))
        stop('`cv_set` must be a CVSettings object')
      
      # Assign inputs
      private$.cv_set    <- cv_set
      private$.seed      <- seed
      
    },
    
    cv      = function(data, alg_settings){
      
     # Get parameters to tune and number of attempts
     params     <- rownames(alg_settings$tune_grid)
     n_attempts <- ncol(alg_settings$tune_grid)
     
     # Matrix for metrics
     eval_matr <- data.frame()

     print_info("Preparing cv arguments")
    
     # Tuning
     for (j in seq(n_attempts)) {

       print_info(paste("Attempt", j, "of", n_attempts, "..."))
       
       # Tune values for j-th attempt
       tune_vals <- as.list(alg_settings$tune_grid[ ,j])
       print_info('Tuning values:')
       print(tune_vals)
       
       cv_res    <- private$.exec_cv(data, alg_settings, tune_vals)
       
       # Add the CV result to matrix
       if (check_type(eval_matr,'data.frame',c(0,0),c(0,0))) {
         eval_matr <- matrix(c(alg_settings$tune_grid[,j],cv_res), nrow = 1)
         colnames(eval_matr) <- c(params, names(cv_res))
         rownames(eval_matr) <- NULL
       }else{
         eval_matr <- rbind(eval_matr, c(alg_settings$tune_grid[,j], cv_res))
       }

       print(eval_matr)
     }

     # Return matrix with metrics
     return(eval_matr)
    },
    
    train   = function(data, alg_set,...){
      
      # Check input
      if (!check_type(data,'mldr'))
        stop('`data` must be a mldr object') 
      
      if (!check_type(alg_set,'AlgSettings'))
        stop('`alg_set` must be a AlgSettings object')
      
      print_info('ML training without cv...')
      set.seed(private$.seed)
      
      # Merge parameters 
      if (length(list(...)) > 0){
        args  <- merge_lists(base.algorithm = alg_set$alg,
                             alg_set$other_params,
                             as.list(...))
      }else{
        # Merge parameters 
        args  <- merge_lists(base.algorithm = alg_set$alg,
                             alg_set$other_params)
      }
      
      print_info('Arguments used for training (+ data):')
      print(args)
      args[['mdata']] <- data
      
      # Training
      model <- do.call(what = alg_set$strat, args)
      
      return(model)
    },
    
    classify     = function(model, data, ...){
      
      return(super$classify(model, data, ...))
      
    },
    
    best_tuning  = function(tune_results, alg_settings){

      # Portion that holds the best parameters
      selection <- tune_results
      print(selection)
      params    <- rownames(alg_settings$tune_grid)
      print(params)
      selection_order <- private$.cv_set$sel_order
      
      # Counter for columns in the selection order
      c <- 1

      # Keep iterating until one row is selected or all columns have been checked
      while (nrow(selection) != 1 & c <= length(selection_order)) {

         # Search '<' to check if order is decreasing
         decreasing <- grepl(selection_order[c],
                              pattern = "<",
                              fixed   = TRUE)

         # Get parameter name by removing sign of decreasing order
         colname <- gsub(x = selection_order[c],
                             pattern     = '<',
                             replacement = '',
                             fixed       = TRUE)
         
         # Use mean of parameter to select rows
         colname <- paste('mean',colname, sep = '_')
         if (decreasing)
            val <- min(selection[,colname])
         else
            val <- max(selection[,colname])
         
         print(val)
         
         # Select rows with required value (max or min)
         selection <- subset(selection, selection[,colname] == val)
         
         # Check the next column
         c <- c + 1
      }
      
      # If columns provided are not sufficient, warn and take the first row
      selection <- as.data.frame(selection)
      if (nrow(selection) > 1){
         print_warning("Selecting first row from following remaining values:")
         print(selection)
      }
      
      # Get the chosen parameters as a list
      print_success('Selected:')
      chosen        <- selection[1,params]
      names(chosen) <- params
      chosen        <- as.list(chosen)
      
      print(chosen)

      return(chosen)

     },
    
    metrics = function(mldr_ref, df_ref, binary_pred){
      
      
      # Convert binarized data into ML result
      mm <- MLMetrics$new()
      
      ml_conf <- list(conf = mm$get_conf_mat(mldr_ref, binary_pred))
      ml_conf[['table']] <-  mm$print_confusion_matrix(ml_conf$conf)
      
      ml_metrics <- list(
        conf    = ml_conf,
        global  = mm$global_metrics(mldr_ref, binary_pred),
        local   = mm$local_metrics(mldr_ref, binary_pred),
        relaxed = mm$relaxed_metrics(df_ref, attr(binary_pred, 'classes'))#,
        #sim     = mm$similarity_metrics(df_ref, binary_pred)
      )
      
      return(ml_metrics)
    }
  ),
  
  active = list(
    cv_set = function(v){
      if (missing(v))
        private$.cv_set
      else
        stop('`cv_set` is read-only')
    }
  )
)