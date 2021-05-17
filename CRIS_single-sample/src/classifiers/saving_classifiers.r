compare_metrics <- function(res, type){
  
  # Comparison of local metrics
  loc_comp <- compare_local_metr(res) %>% 
      t() %>% as.data.frame()
  
  colnames(loc_comp) <- paste(loc_comp[1,], loc_comp[2,], sep = ' ')
  loc_comp <- loc_comp[-(1:2), ]
  
  for (c in colnames(loc_comp)[which(colnames(loc_comp)!= 'Algorithm')]){
    mode(loc_comp[,c]) <- 'numeric'
  }
  
  if (all(type == 'sl')){
    
    # Comparison of metrics (similarity, global, local)
    comparison <- list(
        gl_comp   = compare_metr(res,'sl') %>% rownames_to_column('Algorithm'),
        loc_comp  = loc_comp %>% rownames_to_column('Algorithm')
    )
    
    
  }else if (all(type == 'ml')){
    
    
    # Comparison of metrics (multilabel, similarity, global, local)
    comparison <- list(
      ml_comp   = compare_metr(res,'ml') %>% rownames_to_column('Algorithm'),
      # sim_comp  = compare_metr(res,'sim') %>% rownames_to_column('Algorithm'),
      gl_comp   = compare_metr(res,'gl') %>% rownames_to_column('Algorithm'),
      loc_comp  = loc_comp %>% rownames_to_column('Algorithm')
    )
    
    
  }else{
    stop('compare_metrics: comparison only for type `ml` and `sl`.')
  }
  
  
  return(comparison)
}


compare_metr <- function(res, type = 'gl'){       
  
  comp <- data.frame()
  for(i in seq(length(res))) {
  
    set_n     <- names(res)[i]
    
    # Global metrics (accuracy, kappa)
    if (type == 'gl'){
      metr <- res[[set_n]]$metrics$relaxed
    }
    
    else if (type == 'sl'){
      metr <- as.numeric(res[[set_n]]$metrics$global$Value)
      names(metr) <- res[[set_n]]$metrics$global$Metric
    }
    
    # Similarity metrics (pearson, spearman scores correlation)
    else if (type == 'sim') {
      metr <-  as.numeric(res[[set_n]]$metrics$sim$table$Value)
      names(metr) <- res[[set_n]]$metrics$sim$table$Metric
    }
    
    # Multilabel metrics
    else if (type == 'ml'){
      metr <-  as.numeric(res[[set_n]]$metrics$global$Value)
      names(metr) <- res[[set_n]]$metrics$global$Metric
    }else{
      stop('Supported types: gl, sim, ml')
    }
    
    if (any(dim(comp) == 0)){
      comp <- rownames(metr)
      comp <- as.data.frame(t(metr))
    }else{
      comp <- rbind(comp, t(metr))
    }
    
    rownames(comp)[i] <- set_n
  
  }
  
  return(comp)
}


compare_local_metr <- function(res){
  
  loc_metr <- SLMetrics$new()$list_local_metrics() %>% tolower()

  # Result matrix
  overall_comp <- data.frame()
  classes  <- levels(F_CRIS_CLASSES)
  
  for (m in loc_metr){
    
    # Matrix for current metric (classes x methods)
    loc_comp <- matrix(0, nrow = length(classes), ncol = length(res))
    rownames(loc_comp) <- classes
    colnames(loc_comp) <- names(res)
    
    # For each method, get local metrics
    
    for (n in names(res)){
      metr_vals <- res[[n]]$metrics$local %>% column_to_rownames('Metric')
      loc_comp[classes,n] <- metr_vals[m, ] %>% t()
    }
      
    # Add column for class and metric
    loc_comp <- as.data.frame(loc_comp) %>% 
                rownames_to_column('Class')
    loc_comp <- cbind(Metric = m, loc_comp)
    
    # Add metric results to comparison matrix
    overall_comp <- rbind(overall_comp, loc_comp)
  }
  
  return(overall_comp)
}


prepare_excel_res_sl <- function(res,ref){
  
  pred <- res$pred %>% rownames_to_column(ALIQUOT_LABEL)
  
  ref <- ref 
  colnames(ref)[which(colnames(ref) == CLASS_LABEL)] <- 'true_class'
  for (c in CRIS_CLASSES){
    colnames(ref)[which(colnames(ref) == c)] <- paste('ntp',c, sep = '')
  }
  
  # Prediction
  result_summary <- full_join(pred, ref, by = ALIQUOT_LABEL)
  
  # Binarization
  result_binary  <- res$binary_res %>% rownames_to_column(ALIQUOT_LABEL)
  
  # Confusion matrix
  conf_mat <- res$metrics$conf
  
  # Global metrics
  global <- res$metrics$global
  
  # Local metrics
  local  <- res$metrics$local
  
  excel_res <- list(
    result = result_summary,
    binarization = result_binary,
    confusion_matrix = conf_mat,
    global = global,
    local  = local
  )
  
  return(excel_res)
}


prepare_excel_res_ml <- function(res, ref){

  # Prediction
  classes <- levels(F_CRIS_CLASSES)
  scores <- as.matrix(res$pred) %>%
            as.data.frame() %>%
            rownames_to_column(ALIQUOT_LABEL)
  
  cls_indices <- which(colnames(scores) %in% CRIS_CLASSES)
  colnames(scores)[cls_indices] <- paste('score', colnames(scores)[cls_indices], sep = '')
  
  for (c in paste('score', classes, sep = '')){
    mode(scores[,c]) <- 'numeric'
  }
  
  # Assigned labels
  result_binary <- res$binary_res %>% as.data.frame() %>% rownames_to_column(ALIQUOT_LABEL)
  
  # NTP REFERENCE
  ntp_scores <- ref 
  colnames(ntp_scores)[which(colnames(ntp_scores) == CLASS_LABEL)] <- 'true_class'
  
  # Confusion matrix
  conf_mat_ml    <- as.matrix(res$metrics$conf$table) %>% as.data.frame()
  
  for (c in CRIS_CLASSES){
     mode(conf_mat_ml[,c]) <- 'numeric'
  }
  
  # Global metrics (ML)
  global_ml      <- as.data.frame(res$metrics$global)
  mode(global_ml[ ,'Value']) <- 'numeric'
  
  # Relaxed metrics
  global_relaxed <- data.frame(Value = res$metrics$relaxed) %>% rownames_to_column('Metric')
  mode(global_relaxed[ ,'Value']) <- 'numeric'
  
  # Local metrics
  local_ml  <- res$metrics$local %>% as.data.frame()
  for (c in CRIS_CLASSES){
     mode(local_ml[,c]) <- 'numeric'
  }
  
  excel_res <- list(
    result_scores = scores,
    result_binary = result_binary,
    ntp_scores = ntp_scores,
    conf_mat = conf_mat_ml,
    global_relaxed = global_relaxed,
    global_ml = global_ml,
    local  = local_ml
  )
  
  return(excel_res)
}

