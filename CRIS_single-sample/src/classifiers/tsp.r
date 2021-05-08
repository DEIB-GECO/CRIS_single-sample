# Description -------------------------------------------------------------

# Execute TSP classification and compute metrics on result

# Class definition ------------------------------------------------------------

TSPClassifier <- R6Class(
  'TSPClassifier',
  lock_objects = TRUE,
  lock_class   = TRUE,
  
  #### Private fields ####
  private = list(
    
    .validate_data = function(data){
      
      if (!check_type(data,'ExpressionSet'))
        stop('`data` must be an ExpressionSet object')
        
      if (!check_type(data@assayData$exprs, 'matrix',c(1,1)) | 
          !check_type(data@assayData$exprs, 'numeric'))
        stop('data assayData Matrix must be numeric and non empty')
      
      if (check_type(colnames(data@assayData$exprs), 'null') | 
          check_type(rownames(data@assayData$exprs), 'null'))
        stop('data assayData Matrix have non null rows and non null columns')
    }
    
  ),
  
  #### Public fields ####
  
  public = list(
    

    # data is expression set
    classify = function(data){
      
      print_info('Check input')
      
      private$.validate_data(data)
      
      print('Prediction')
      tsp_pred  <- predictCRISclassKTSP(data@assayData$exprs)
      
      tsp_class <- as.data.frame(tsp_pred$tspSetClassPredsFinal)
      colnames(tsp_class) <- CLASS_LABEL
      tsp_class <- tsp_class %>% rownames_to_column(ALIQUOT_LABEL)
      
      # Factorize labels
      tsp_class[,CLASS_LABEL] <- tsp_class[,CLASS_LABEL] %>%
           gsub(pattern = 'CRIS', replacement = 'CRIS-') %>%
           factor(levels = CRIS_CLASSES)
      
      print('Scores')
  
      tsp_scores <- as.data.frame(tsp_pred$tspSetClassPercent)
      colnames(tsp_scores) <- gsub(colnames(tsp_scores),
                                   pattern = 'CRIS',
                                   replacement = 'CRIS-',
                                   fixed = TRUE)
      tsp_scores <- tsp_scores %>% rownames_to_column(ALIQUOT_LABEL)
      
      print('Results')
  
      tsp_result <- inner_join(tsp_scores, tsp_class, by = ALIQUOT_LABEL) %>%
                    column_to_rownames(ALIQUOT_LABEL)
      
      tsp_res   <- list(
        raw   = tsp_pred,
        cl    = tsp_result %>% filter(!is.na(predict.label2)) %>% as.data.frame(),
        uncl  = tsp_result %>% filter(is.na(predict.label2))  %>% as.data.frame()
      )
      
      return(tsp_res)
    },
    
    metrics = function(ref_data, pred_data){
      
      # Metrics computer
      sm  <- SLMetrics$new()
      
      # Consider only classified samples
      pred <- pred_data$cl
      ref  <- ref_data[rownames(pred), ]

      sl_metrics <- list(
        conf   = sm$confusion_matrix(ref, pred),
        global = sm$global_metrics(ref, pred),
        local  = sm$local_metrics(ref, pred)
      )
      
      return(sl_metrics)
    },
    
    prepare_for_saving = function(res, ref){
  
      # Reference class and scores
      ref <- ref[,c(ALIQUOT_LABEL, CLASS_LABEL, CRIS_CLASSES)]
      colnames(ref) <- c(ALIQUOT_LABEL, 
                         'true_class', 
                         paste('ntp', CRIS_CLASSES, sep = ''))
      
      # Classification and scores
      classif <- rbind(res$result$cl, res$result$uncl)
      result_summary <- full_join(classif, ref, by = ALIQUOT_LABEL) %>%
                        full_join(res$result$score, by = ALIQUOT_LABEL)
      
      # Summary metrics
      gl <- c(res$metrics$sl$global, res$metrics$sim$avg) %>% unlist()
      gl <- data.frame(Value = gl)
      rownames(gl) <- c(names(res$metrics$sl$global), names(res$metrics$sim$avg))
      gl <- gl %>% rownames_to_column('Metric')
      
      # Correlations
      correl  <- full_join(res$metrics$sim$pearson, 
                           res$metrics$sim$spearman, by = ALIQUOT_LABEL)
        
      # Store in list
      what.to.save <- list(
        result  = result_summary,
        confusion_matrix = res$metrics$sl$conf,
        local   = res$metrics$sl$local %>% rownames_to_column('Metric'),
        global  = gl,
        scores_correlations = correl
      )
      
      return(what.to.save)
    }
    
  )
  
)




