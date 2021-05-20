# Description -------------------------------------------------------------

# Class for the execution of single-label classifiers.

# Class definition --------------------------------------------------------

SLClassifier <- R6Class(
  'SLClassifier',
  inherit  = SingleSampleClassifier,
  lock_objects = TRUE,
  lock_class   = TRUE,

  #### Private fields ####
  
  private = list(
    .method         = NULL,  # method name for caret library (see caret available models)
    
    #' Get the data for the ROC curves, considering also the case where the dataset
    #' is a PercentileMLData object.If training flag is TRUE, get training data,
    #' otherwise get testing data.
    #' 
    #' @param dataset   the data (SLData or PercentileMLData)
    #' @param training  a boolean flag; if TRUE, get training data and ref, otherwise
    #' get testing data and ref.
    #' 
    #' @return a list with data and ref fields.
    .get_roc_data = function(dataset,training){
      
      roc_data <- super$.get_roc_data(dataset,training)  
      
      if (check_type(dataset, 'PercentileMLData'))
        roc_data$data <- roc_data$data$dataset
    
      return(roc_data)
    }
  ),
  
  #### Public fields ####
  
  public = list(
    
    #' Constructor
    #' 
    #' @param method the name of the classifier method, used by caret
    #' @param seed   an integer seed to allow replication
    #' @param cl_thresholds  class-specific numeric thresholds. Default is NULL.
    initialize = function(method, seed, cl_thresholds = NULL){
      
      super$initialize(seed, cl_thresholds)
      
      # Check method
      if(!class(method) == 'character')
        stop('`method` must be a string')
      
      if(length(method) != 1)
        stop('`method` must be a unique string')
      
      # Assign inputs
      private$.method <- method
      
 
    },
    
    #' Apply cross-validation to tune hyperparameters of the model and train again
    #' on the entire dataset with the best hyperparameter values
    #' 
    #' @param formula  The formula for the model
    #' @param data     The (training) data on which to apply the model (must be a data.frame)
    #' @param fit_control trainControl function for specifying cross-validation parameters
    #' @param tune_grid  Tune grid with the parameters to tune. 
    #' 
    #' @return the computed model
    cv       = function(formula, data, fit_control, tune_grid, ...){
      
      print_info('Training with cv')
      
      adj_data <- data
      adj_data[,CLASS_LABEL] <- make.names(adj_data[,CLASS_LABEL])
      
      set.seed(private$.seed)
      
      model <- train(formula, adj_data, 
                     method    = private$.method,
                     trControl = fit_control,
                     tuneGrid  = tune_grid)
      
      return(model) 
      
    },
    
    
    #' Train the model with the given formula on the required data. No cross-validation
    #' is applied.
    #' 
    #' @param formula  the formula of the model
    #' @param data     the training data (data.frame)
    #' 
    #' @return The trained model
    train    = function(formula, data){
      
      print_info('Training')
      
      set.seed(private$.seed)
      adj_data <- data
      adj_data[,CLASS_LABEL] <- make.names(adj_data[,CLASS_LABEL])
      
      model <- caret::train(formula, adj_data, 
                            method = private$.method,
                            trControl = trainControl(classProbs = TRUE))
      
      return(model)
    },
    
    
    #' Classify the data given the provided model
    #' 
    #' @param model the model to apply
    #' @param data  the data to apply (data.frame)
    #' 
    #' @return classification result with class scores and assigned class (the one
    #' with the highest score)
    classify = function(model, data){
      
      # Adjust class labels
      ignored_cols <- c(CRIS_CLASSES, CLASS_LABEL, '.labelcount','.SCUMBLE')
      adj_data <- data[,-which(colnames(data) %in% ignored_cols)]
      
      
      # Predict
      pred  <- predict(model, newdata = adj_data, type = 'prob') %>%
                      format(scientific = FALSE) %>%
                      as.data.frame()
      
      # Adjust columns mode
      for (c in colnames(pred))
        pred[,c] <- as.numeric(unlist(pred[,c]))
    
      
      # Samples as rownames, class labels as colnames
      rownames(pred) <- rownames(adj_data)
      colnames(pred) <- gsub(colnames(pred), 
                             pattern     = 'CRIS.', 
                             replacement = 'CRIS-', 
                             fixed = TRUE)
    
      pred <- pred[,CRIS_CLASSES]
      
      # Class with max score (factorized)
      cls      <- apply(pred[,CRIS_CLASSES], 1, which.max) %>% as.numeric()
      fact_cls <- factor(colnames(pred)[cls], CRIS_CLASSES)
      
      # Add class with max score to result
      pred <- cbind(as.data.frame(pred), fact_cls)
      colnames(pred)[ncol(pred)] <- CLASS_LABEL
      
      # Return result
      return(as.data.frame(pred))
    },
    
    
    #' Overwrite the superclass method to compute class thresholds through ROC curves.
    #' The procedure is the same described for SingleSampleClassifier class but,
    #' before computing the ROC curves, a softmax layer is removed and scores are normalized
    #' into [0-1]: each score s (for class c) is normalized with (s-min)/(max-min), where 
    #' max and min are maximum and minimum values of the scores for the required class (computed
    #' on training data only).
    #' 
    #' @param model The model of the classifier
    #' @param dataset The data (SLData/MLData) from which taking the data to be classified.
    #' Both SLData and MLData can be used since this function is used also for computing
    #' ROC curves on testing data for single-label classifier adapted to multi-label context.
    #' These curves maybe used only for verifications and are not further exploited.
    #' @param png_path The path to save the ROC curve
    #' @param training Boolean flag; if TRUE (default), compute the ROC on traning data,
    #' otherwise use testing data.
    #' @param max_cls  a vector with maximum scores for each class (computed on training), used
    #' for the normalization
    #' @param min_cls  a vector with minimum scores for each class (computed on training), used
    #' for the normalization
    #' 
    #' @return the class-specific thresholds (computed with private$.roc_thr)
    roc_thresholds  = function(model, dataset, png_path, training = TRUE, max_cls, min_cls){
      
      roc_data <- private$.get_roc_data(dataset, training)
      data <- roc_data$data
      ref  <- roc_data$ref
    
      # Predict model on data
      pred <- self$classify(model, data)
      
      # Remove the softmax layer
      for (i in 1:nrow(pred)){
        pred[i, CRIS_CLASSES] <- ohenery::inv_smax(as.matrix(pred[i,CRIS_CLASSES]))
      }
      
      
      # Move into range [0,1], removing Infinite values
      for (c in CRIS_CLASSES){
        min_not_inf <- min_cls[c] 
        max_not_inf <- max_cls[c]
        pred[pred[,c] == Inf, c]  <- max_not_inf
        pred[pred[,c] == -Inf, c] <- min_not_inf
        if (max_not_inf == min_not_inf)
          pred[,c] <- 0
        else
          pred[,c] <- as.double(pred[,c] - min_not_inf)/as.double(max_not_inf - min_not_inf)
      }
       
      # Through ROC curves (saved in PNG) get class thresholds
      png(png_path)
      thresholds  <- private$.roc_thr(pred, ref) %>% 
                     as.data.frame()
      dev.off()
      
      return(thresholds)
    },
    
    
    #' Binarize the result of single-label classifier adapted to multi-label using
    #' the superclass procedure. In addition, maintain the original predicted class,
    #' forcing to 1 the assignment even if it is not above the threshold.
    #' 
    #' @param prediction the result obtained by the classification
    #' @return the result with binary assignment and the original single-label class,
    #' together with number of assigned classes and classification flag (at least one class is assigned).
    binarize_result_as_ml = function(prediction){
      
      pred_binary <- super$binarize_result(prediction)
      for (c in CRIS_CLASSES){
        if (length(pred_binary[pred_binary[,CLASS_LABEL] == c, c]) > 0)
          pred_binary[pred_binary[,CLASS_LABEL] == c, c] <- 1
      }
      pred_binary[,'n_classes'] <- rowSums(pred_binary[,CRIS_CLASSES])
      pred_binary[,'classified'] <- TRUE
      return(pred_binary[,c(CLASS_LABEL, CRIS_CLASSES,'n_classes','classified')])
      
    },
    
    
    #' Binarize the result of single-label classifier, assigning 1 to the assigned
    #' class and 0 to the others.
    #' 
    #' @param prediction the result obtained by the classification
    #' @return the result with binary assignment and the original single-label class,
    #' together with number of assigned classes and classification flag.
    #' The number of assigned classes is always 1.
  
    binarize_result_as_sl = function(prediction){
      
      pred_binary <- super$binarize_result(prediction) %>% rownames_to_column(ALIQUOT_LABEL)
      
      # Force binarization to be 1 on best class, 0 on the others
      for (c in CRIS_CLASSES){
        if (length(pred_binary[pred_binary[,CLASS_LABEL] == c, c]) > 0)
          pred_binary[pred_binary[,CLASS_LABEL] == c, c] <- 1
        if (length(pred_binary[pred_binary[,CLASS_LABEL] != c, c]) > 0)
          pred_binary[pred_binary[,CLASS_LABEL] != c, c] <- 0
      }
      pred_binary <- pred_binary %>% column_to_rownames(ALIQUOT_LABEL)
      pred_binary[,'n_classes'] <- rowSums(pred_binary[,CRIS_CLASSES])
      pred_binary[,'classified'] <- TRUE
      return(pred_binary[,c(CLASS_LABEL, CRIS_CLASSES,'n_classes','classified')])
      
    },
    
    
    #' Compute single-label classification metrics (confusion matrix, global
    #' metrics, local metrics)
    #' 
    #' @param ref_data the reference targets
    #' @param pred_data the predicted class (CLASS_LABEL column with aliquots)
    #' @return a list with confusion matrix (conf), global metrics, namely accuracy and
    #' kappa score (global) and local metrics, namely precision, recall, f1 score and
    #' mcc for each class (local)
    metrics  = function(ref_data, pred_data){
      
      sm <- SLMetrics$new()

      sl_metrics <- list(
        conf   = sm$confusion_matrix(ref_data, pred_data),
        global = sm$global_metrics(ref_data, pred_data),
        local  = sm$local_metrics(ref_data, pred_data)
      )
      
      return(sl_metrics)
      
    }
    
  )
)


