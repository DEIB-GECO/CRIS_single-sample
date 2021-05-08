
# Description -------------------------------------------------------------

# Class for metrics of single-label classifiers (accuracy, kappa, precision,
# recall, f1 score, mcc)

# Class definition --------------------------------------------------------


SLMetrics <- R6Class(
  'SLMetrics',
  inherit = Metrics,
  lock_objects = FALSE,
  lock_class = TRUE,
  
  ##### Private fields #####
  
  private = list(
    
    #' Computes a confusion matrix with respect to a specific class (used for additional metrics,
    #' e.g. mcc score)
    #' 
    #' @param conf Confusion matrix table obtained from caret::confusionMatrix
    #' @param class_name  name of the class with respect to which computing 
    #' true positives (tp), false positives (fp), false negatives (fn), true negatives (tn).
    #' 
    #' @return list with tp, tn, fp, fn computed with respect to the required class.
    .class_conf_mat     = function(conf, class_name){
      
      if (!is.character(class_name))
        stop(paste('class name must be a character') )
      
      if (any(!class_name %in% CRIS_CLASSES))
        stop(paste('Cannot found the class',class_name,'for class confusion matrix'))
      
      conf <- conf %>% as.data.frame()
      
      tp <- conf %>% filter(Prediction == class_name & Reference == class_name)
      fp <- conf %>% filter(Prediction == class_name & Reference != class_name)
      tn <- conf %>% filter(Prediction != class_name & Reference != class_name)
      fn <- conf %>% filter(Prediction != class_name & Reference == class_name)
      
      conf <- list(
        tp = sum(tp$Freq),
        tn = sum(tn$Freq),
        fp = sum(fp$Freq),
        fn = sum(fn$Freq)
      )
      

      return(conf)
    },
    
    
    
    #' Check that the samples of reference and prediction coincide. Rownames
    #' of reference and prediction are aliquot ids.
    #' 
    #' @param ref reference target
    #' @param pred prediction
    #' 
    #' @return stops the execution if the correspondence is not verified.
    .check_correspondence = function(ref, pred){
      
      if (!any(all.equal(rownames(ref), rownames(pred)) == TRUE))
        stop('There is no sample correspondence in prediction and reference')
    
    },
    
    
    #' Check the provided object is a data.frame with class label column
    #' 
    #' @param df the object (data.frame) to check
    #' 
    #' @return stops the execution if the check is not verified.
    .check_df_format = function(df){
      
      if (!class(df) %in% c('data.frame'))
        stop('Can use only data.frame')
      
      if (!CLASS_LABEL %in% colnames(df))
        stop(paste('Must contain',CLASS_LABEL,'column'))
      
    },
    
    
    #' Returns the class label column of the provided data.frame as a factor.
    #' 
    #' @param df the object (data.frame) to use
    .get_class_factor = function(df){
      
      return(factor(df[,CLASS_LABEL], levels = CRIS_CLASSES))
    },
    
    
    #' Returns the confusion matrix of single-label classifier using the
    #' caret library funciton (caret::confusionMatrix)
    #' 
    #' @param ref the reference with class targets. It must contain a class label column
    #' @param pred the predicted outcome. It must contain a class label column
    #' 
    #' @return the confusion matrix obtained with the caret library function.
    .get_conf_mat = function(ref, pred){
      
      private$.validate(ref,pred)
      
      # Factorize classes
      ref_classes  <- private$.get_class_factor(ref)
      pred_classes <- private$.get_class_factor(pred)
      
      # Compute confusion matrix
      confMat <- caret::confusionMatrix(data = pred_classes, reference = ref_classes, mode = 'everything')
      
      return(confMat)
    },
    
    
    #' Validates the reference and the prediction objects, used to compute the metrics
    #' 
    #' @param ref the reference with class targets. It must contain a class label column
    #' @param pred the predicted outcome. It must contain a class label column
    #' 
    #' @return in case of errors, halts the execution of the program.
    .validate = function(ref,pred){
      
      # Check input
      private$.check_df_format(ref)
      private$.check_df_format(pred)
      
      # Check samples correspondence
      private$.check_correspondence(ref, pred)
      
    }
  ),
  
  ##### Public fields #####
  
  public = list(
    
    #' Returns the global metrics that can be computed
    #' @return character vector with the computed global metrics
    list_global_metrics = function(){
      return(c('Accuracy','Kappa'))
    },
    
    #' Returns the local metrics that can be computed
    #' @return character vector with the computed local metrics
    list_local_metrics = function(){
      return(c('Precision','Recall','F1','mcc'))
    },

    
    #' Uses the private function .get_conf_mat to compute a confusion matrix
    #' for single-label classifier and returns it with correct formatting.
    #' 
    #' @param ref the reference with class targets. It must contain a class label column
    #' @param pred the predicted outcome. It must contain a class label column
    #' 
    #' @return a data.frame containing the confusion matrix. Prediction is on the rows
    #' and reference is on the columns
    confusion_matrix = function(ref, pred){
      
      confMat <- private$.get_conf_mat(ref,pred)
      confMat <- confMat %>% as.matrix() %>% as.data.frame() %>% rownames_to_column('Pred x Ref')
      
      return(confMat)
      
    },
    
    
    #' Uses the private function .get_conf_mat to compute a confusion matrix
    #' for single-label classifier and returns the global metrics obtained from it
    #' 
    #' @param ref the reference with class targets. It must contain a class label column
    #' @param pred the predicted outcome. It must contain a class label column
    #' 
    #' @return a data.frame containing the global metrics and their value.
    global_metrics = function(ref,pred){
      
      confMat <- private$.get_conf_mat(ref,pred)
      
      global <- data.frame(Value = confMat$overall[c('Accuracy','Kappa')]) %>% 
        rownames_to_column('Metric')
      
      return(global)
    },
    
    
    
    #' Uses the private function .get_conf_mat to compute a confusion matrix
    #' for single-label classifier and returns the local metrics obtained from it. In
    #' addition, an MCC score for each class is added.
    #' 
    #' @param ref the reference with class targets. It must contain a class label column
    #' @param pred the predicted outcome. It must contain a class label column
    #' 
    #' @return a data.frame containing the local metrics and their value.
    local_metrics = function(ref,pred){
      
      
      # Get local measures from Confusion Matrix
      confMat <- private$.get_conf_mat(ref,pred)
      local   <- confMat$byClass %>% t() %>% as.data.frame()
      
      
      # Adjust class names
      colnames(local)  <- gsub(colnames(local), pattern = 'Class: ', replacement = '')
      computed_metrics <- rownames(local)
      
      for (m in self$list_local_metrics()) {

        # Add missing metrics (mcc score)
        if (!m %in% computed_metrics){
          
          # Empty row for new metric
          local   <- rbind(local, NA)
          
          rownames(local)[nrow(local)] <- m
          
          # Fill the table with metric value for each class
          for (c in CRIS_CLASSES) {
            
            args  <- list(conf = private$.class_conf_mat(confMat$table, c))
            local[m,c] <- do.call(get(m, envir = private), args)
          }
        }

      }
      local <- local[self$list_local_metrics(), ]
      rownames(local) <- tolower(rownames(local))
      
      return( local %>% rownames_to_column('Metric'))
    }
    
  )
)
