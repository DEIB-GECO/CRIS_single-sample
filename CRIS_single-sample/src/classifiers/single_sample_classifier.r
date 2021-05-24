
# Description -------------------------------------------------------------

# Base class for single-sample classifiers (both single-label and multi-label)


# Class definition --------------------------------------------------------


SingleSampleClassifier <- R6Class(
  
  'SingleSampleClassifier',
  lock_objects = TRUE,
  lock_class   = TRUE,
 
  ##### Private fields ###### 

  private = list(
    .seed   = NULL,
    .cl_thresholds    = NULL,
  
    #' Add ranks of classes depending on the weighted scores (class scores/class thresholds)
    #' 
    #' @param pred_binary the assigned classes with continuous score
    #' @param classes     the attributes relative to class scores
    .add_ranks  = function(pred_binary, classes){
      
      if (is.null(private$.cl_thresholds)){
        return(pred_binary)
      }
      
      rank_labels     <- paste('rank', classes, sep = '')
      w_scores_labels <- paste('weighted', classes, sep = '')
      
      # Empty weighted scores
      w_scores <- matrix(0, nrow = nrow(pred_binary), ncol = length(classes))
      colnames(w_scores) <- w_scores_labels
      
      # Empty ranks
      ranks <- matrix(0, nrow = nrow(pred_binary), ncol = length(classes))
      colnames(ranks) <- rank_labels
      
      # Add weighted scores and ranks to result
      pred_binary     <- cbind(pred_binary, as.data.frame(w_scores))
      pred_binary     <- cbind(pred_binary, as.data.frame(ranks))
      
      # Fill the ranks when at least one class is assigned
      for (k in seq(nrow(pred_binary))){
        
        # Get and sort (descending) scores for current sample
        scores        <- pred_binary[k, classes] %>% 
                         unlist() %>% 
                         as.numeric() %>%
                         signif(8)
        
        # Divide scores by thresholds
        weighted_scores <- scores/private$.cl_thresholds[classes, 'threshold']
        sorted_w_scores <- sort(unique(weighted_scores), decreasing = TRUE)
        
        # Get position (rank) of each score in the sorted list of unique scores
        ranking <- lapply(X   = seq(length(weighted_scores)), 
                          FUN = function(k){which(sorted_w_scores == weighted_scores[k])}) 
        
        # Assign weighted scores and ranking
        pred_binary[k,w_scores_labels] <- weighted_scores
        pred_binary[k,rank_labels]     <- ranking
        
      }
      
      return(pred_binary)
      
    },
    
    #' Get the data for roc curves depending on training flag. If TRUE, get training data,
    #' otherwise get testing data.
    #' 
    #' @param dataset   the data (SLData or PercentileMLData)
    #' @param training  a boolean flag; if TRUE, get training data and ref, otherwise
    #' get testing data and ref.
    #' 
    #' @return a list with data and ref fields.
    .get_roc_data = function(dataset, training){
      
      if (training){
        data <- dataset$train_
        ref  <- dataset$train_ref
      }else{
        data <- dataset$test_
        ref  <- dataset$test_ref
      }
      
      roc_data <- list(
        data = data,
        ref = ref
      )
      
      return(roc_data)
    },
    
    #' Get scores for the required class (used in ROC curve computation)
    #' 
    #' @param res  the obtained scores
    #' @param cl   the CRIS class of which extracting the scores
    #' 
    #' @return data.frame with aliquot ID and required class score, sorted by aliquot ID
    .get_cl_res = function(res,cl){
      
      cl_res <- res %>% 
          as.matrix() %>% 
          as.data.frame()
      if (!ALIQUOT_LABEL %in% colnames(cl_res))
        cl_res <- cl_res %>% rownames_to_column(ALIQUOT_LABEL)
      
      cl_res <- cl_res[,c(ALIQUOT_LABEL,cl)] %>% arrange(aliquot_id)
      
      return(cl_res)
    },
    
    #' Get binary reference assignment for the required class (used in ROC curve computation)
    #' for a single-label classifier
    #' 
    #' @param ref  the reference assignment (containing aliquot ID and predicted class)
    #' @param cl   the CRIS class of which extracting the assignment
    #' 
    #' @return data.frame with aliquot ID and required class assignment (positive column), sorted by aliquot ID.
    #' Assignment is 1 if predicted class is equal to the required class (cl), otherwise
    #' it is 0.
    .get_cl_ref_sl = function(ref,cl){
      print_debug('Binary class reference for ROC curves (single-label classifier)')
      cl_ref <- ref[ ,c(ALIQUOT_LABEL, CLASS_LABEL)] %>% mutate(positive = 0)
      cl_ref[cl_ref$predict.label2 == cl, 'positive'] <- 1
      
      return(cl_ref)
    },
    
    
    #' Get binary reference assignment for the required class (used in ROC curve computation)
    #' for a multi-label classifier. NB: assume that the original target is already binarized.
    #' 
    #' @param ref  the reference assignment (containing aliquot ID and binarized target scores)
    #' @param cl   the CRIS class of which extracting the assignment
    #' 
    #' @return data.frame with aliquot ID and required class assignment (positive column), 
    #' sorted by aliquot ID. 
    .get_cl_ref_ml = function(ref,cl){
      print_debug('Binary class reference for ROC curves (multi-label classifier)')
      cl_ref <- ref[ ,c(ALIQUOT_LABEL, cl)]
      colnames(cl_ref) <- c(ALIQUOT_LABEL, 'positive')
      return(cl_ref)
    },
    
    
    #' Get binary reference assignment for the required class (used in ROC curve computation)
    #' depending on the type of classifier. If the target scores are binary numbers, use the
    #' multi-label version (.get_cl_ref_ml), otherwise use the single-label version (.get_cl_ref_sl).
    #' 
    #' @param ref  the reference assignment (containing aliquot ID and target scores)
    #' @param cl   the CRIS class of which extracting the assignment
    #' 
    #' @return data.frame with aliquot ID and required class assignment (positive column), 
    #' sorted by aliquot ID. 
    .get_cl_ref = function(ref,cl){
      
      if (all(unique(ref[,cl]) %in% c(0,1)))
        return(private$.get_cl_ref_ml(ref,cl))
      else
        return(private$.get_cl_ref_sl(ref,cl))

    },
    
    
    #' Compute the roc curves and return the class thresholds maximizing the true positive
    #' percentage (TPP) and minimizing the false positive percentage (FPP)
    #' 
    #' @param res result of classification
    #' @param ref original target assignment
    #' 
    #' @return data.frame with, class, TPP, FPP and threshold
    .roc_thr    = function(res, ref){
      
      # Colors and classes for the plot
      colors   <- c('#FD9696','#78A4F0','#4DD049','#FDCF38','#E656C5')
      classes  <- CRIS_CLASSES
      names(colors) <- classes
      
      # Constrain plot into min and max data range (x-axis)
      par(pty = 's')
      
      # Hold roc results
      roc.info.df <- list()
      
      # Hold thresholds (have tpp = 100% and min fpp)
      thresholds <- matrix(NA, nrow = length(classes), ncol = 3)
      colnames(thresholds) <- c('tpp','fpp','threshold')
      rownames(thresholds) <- classes
      
      for (cl in classes) {
        
        # Get class scores
        cl_res <- private$.get_cl_res(res,cl)

        # Get original class (1 if belonging to cl, 0 otherwise)
        cl_ref <- private$.get_cl_ref(ref, cl) %>% 
                  filter(aliquot_id %in% cl_res$aliquot_id) %>% 
                  arrange(ALIQUOT_LABEL)

        # Debug: check order of samples
        print_debug(all.equal(cl_res$aliquot_id, cl_ref$aliquot_id))
        
        # Start a new plot or use an existing one
        if (cl == classes[1]){
          # Reset plot
          auc_col_pos <- 40
          add_to_plot <- FALSE
        }else{
          auc_col_pos <- auc_col_pos - 5
          add_to_plot <- TRUE
        }
        
        # Draw ROC curve
        roc.res <- pROC::roc(cl_ref$positive, as.numeric(cl_res[,cl]),
            # Plot customization arguments
            plot        = TRUE,
            legacy.axes = TRUE, 
            percent     = TRUE,
            xlab = '% False positives',
            ylab = '% True positives',
            col   = colors[cl],
            lwd   = 3,
            add   = add_to_plot,
            print.auc = TRUE,
            print.auc.y = auc_col_pos)
      
        # Save thresholds used in each point of the ROC curve
        roc.info.df[[cl]] <- data.frame(
          tpp = roc.res$sensitivities,
          fpp = (100 - roc.res$specificities),
          thresholds = roc.res$thresholds
        )
        
        # Find the optimal threshold for the current class
        thresholds[cl, ] <- roc.info.df[[cl]] %>%
          filter(tpp == max(tpp)) %>%
          filter(fpp == min(fpp)) %>%
          filter(thresholds == min(thresholds)) %>%
          unlist() %>%
          as.numeric()
        
      }
      
      # Add legend to the plot
      legend('bottomright', legend = names(colors), col = colors, lwd = 4)
      
      # Restore parameter for plot type (x-axis not constrained)
      par(pty = 'm')
    
      # Add column with class names
      thresholds <- cbind(class = rownames(thresholds), as.data.frame(thresholds))
      
      # Return data.frame with class, tpp, fpp and threshold
      return(thresholds)
  
    }

  ),
  
  ##### Public fields ###### 
  
  public  = list(
    
    #' Constructor
    #' 
    #' @param seed  The seed to allow replications
    #' @param cl_thresholds class thresholds (vector with a numeric threshold for each CRIS class) 
    #' (optional. Default: NULL)
    initialize = function(seed, cl_thresholds = NULL){
      
      # Check and assign seed
      if (!class(seed) == 'numeric')
        stop('`seed` must be a number')
      
      if (round(seed) != seed)
        stop('`seed` must be an integer number')
      
      if (length(seed) != 1)
        stop('`seed` must be an unique integer number')
      
      private$.seed   <- seed
      
      # Check and assign class thresholds, if provided
      if (!is.null(cl_thresholds)) {
        
        if (!class(cl_thresholds) == 'data.frame')
          stop('`cl_thresholds` must be a data.frame')
        
        if (!all.equal(rownames(cl_thresholds), CRIS_CLASSES))
          stop('`cl_thresholds` must have a row for each class')
        
        if (!all.equal(colnames(cl_thresholds), c('class','tpp','fpp','threshold')))
          stop('`cl_thresholds` must have tpp, fpp, threshold columns')
        
        if (any(is.na(cl_thresholds[,'threshold'])))
          stop('`cl_thresholds` cannot have NA thresholds')
        
        if (mode(cl_thresholds[,'threshold']) != 'numeric')
          warning('`cl_thresholds` must have a numeric threshold column')
        
        private$.cl_thresholds   <- as.data.frame(cl_thresholds)
        mode(private$.cl_thresholds$threshold) <- 'numeric'
      }
      
    
    },
    
    #' Obtain a binarized result (binary assignment of class depending on scores).
    #' class score >= class threshold means 1, 0 otherwise.
    #' 
    #' @param prediction the result of the classifier prediction
    #' @return binarized prediction with weighted scores and class ranks too.
    binarize_result = function(prediction){
      
      classes <- CRIS_CLASSES
      if (any(class(prediction) == 'mlresult')){
        pred_binary <- prediction %>% as.matrix() %>% as.data.frame()
      }else{
        pred_binary <- prediction[,classes] %>% as.matrix() %>% as.data.frame()
        pred_binary <- cbind(pred_binary, prediction[,CLASS_LABEL])
        colnames(pred_binary)[ncol(pred_binary)] <- CLASS_LABEL
      }
      
      # Add ranking basing on scores
      pred_binary <- private$.add_ranks(pred_binary, classes)
      
      # Binarize scores (class score >= class threshold means 1, 0 otherwise)
      for (cl in classes){
        thr <- private$.cl_thresholds[cl, 'threshold']

        if (length(pred_binary[pred_binary[,cl] < thr, cl]) > 0){

          pred_binary[pred_binary[,cl] < thr, cl]  <- 0
        }
        
        if (length(pred_binary[pred_binary[,cl] >= thr, cl]) > 0){

          pred_binary[pred_binary[,cl] >= thr, cl] <- 1
        }
      }
 
      # Add number of classes assigned at each sample
      n_classes   <- rowSums(pred_binary[,classes])
      pred_binary <- pred_binary %>% cbind(n_classes)
      
      # At least one class has been assigned
      classified <-  lapply(seq(nrow(pred_binary)), function(k){
        if (pred_binary[k,'n_classes'] < 1)
          return(FALSE)
        else
          return(TRUE)
      }) %>% unlist()
      
      pred_binary <- pred_binary %>% cbind(classified)
      return(pred_binary)
    },
    
    #' Public function to compute the class thresholds through roc curves
    #' 
    #' @param model The model of the classifier
    #' @param dataset The data (SLData or MLData) from which taking the data to be classified
    #' @param png_path The path to save the ROC curve
    #' @param training Boolean flag; if TRUE (default), compute the ROC on traning data,
    #' otherwise use testing data.
    #' 
    #' @return the class-specific thresholds (computed with private$.roc_thr)
    roc_thresholds  = function(model, dataset, png_path, training = TRUE){
      
      roc_data <- private$.get_roc_data(dataset, training)
      data <- roc_data$data
      ref  <- roc_data$ref
    
      # Predict model on train data
      pred <- self$classify(model, data)
       
      # Through ROC curves (saved in PNG) get class thresholds
      png(png_path)
      thresholds  <- private$.roc_thr(pred, ref) %>% 
                     as.data.frame()
      dev.off()
      
      return(thresholds)
    },
    
    classify = function(model, data, ...){
      
      return(predict(model,data))
      
    }
    
  ),
  
  ##### Active fields ###### 
  # Read-only access to seed fields, read and write access to class thresholds
  
  active  = list(
    seed   = function(v){
      if (missing(v))
        private$.seed
      else
        stop('SingleSampleClassifier: `seed` is read-only')
    },
    
    cl_thresholds = function(v){
      if (missing(v))
        private$.cl_thresholds
      else
        private$.cl_thresholds <- v
    }
  ),
)
