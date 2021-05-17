
# NTP replication ---------------------------------------------------------

ntp_pipeline <- function(crisdata, features, repl_id, ntp_ref = NULL, n_resempl = NULL){
  
  # Check n resamplings input
  if (any_uncomplete(n_resempl) & any(!is.null(n_resempl)))
    stop('NTP pipeline: cannot accept infinite, NaN, NA as number of resamplings.')
  else if (check_type(n_resempl, 'null', length_min = 1, length_max = 1))
    ntp <- NTPClassifier$new(features)
  else if (all(check_type(n_resempl, 'numeric', length_min= 1, length_max = 1) &
           n_resempl >= 1))
    ntp <- NTPClassifier$new(features, n_resempl = n_resempl)
  else
    stop('NTP pipeline: invalid number of resamplings.')

  # Replication
  ntp_pred   <- ntp$classify(crisdata$data, id = repl_id)
  conf_count <- SummaryMaker$new()$ntp_summary_by_class(ntp_pred)
  ntp_summ   <- SummaryMaker$new()$count_by_attribute(ntp_pred, CLASS_LABEL)

  # Comparison with reference, if any
  if (!is.null(ntp_ref)){

    # Compare results sample by sample
    comparison <- NTPComparator$new()$compare_results(ntp_ref, ntp_pred)

    # Metrics
    ref_id_col <- grep(colnames(ntp_ref), pattern = 'id', fixed = TRUE)
    ref_data   <- subset(ntp_ref, ntp_ref[,ref_id_col] %in% comparison$id_1)
    rownames(ref_data) <- ref_data[,ref_id_col]

    pred_data  <- subset(ntp_pred, ntp_pred$aliquot_id %in% comparison$id_2)
    rownames(pred_data) <- pred_data$aliquot_id

    slm <- SLMetrics$new()
    metrics <- list(
      conf = slm$confusion_matrix(ref_data, pred_data),
      global = slm$global_metrics(ref_data, pred_data),
      local = slm$local_metrics(ref_data, pred_data)
    )

  }else{
    comparison  <- NULL
    metrics     <- NULL
  }

  ntp_result <- list(
    data = crisdata$data,
    result  = ntp_pred,
    summary = ntp_summ,
    confident  = conf_count,
    comparison = comparison,
    metrics    = metrics#,
    # cl_distr   = plot_distr(ntp_summ, percentage = TRUE, repl_id)
  )

  result_file <- paste('NTP_',repl_id, '.rds', sep = '')
  saveRDS(ntp_result,paste(ntp$output_folder, result_file, sep = '/'))
  return(ntp_result)
}



# ML classifier -----------------------------------------------------------

ml_pipeline_train <- function(mldata, seed, cv_set, alg_set, tune = FALSE){

  # Declare classifier
  if (is.null(cv_set))
    mlc <- MLClassifier$new(CVSettings$new(), seed)
  else
    mlc <- MLClassifier$new(cv_set, seed)

  # Training and (if required) tuning
  if (tune) {
    print_info('Tuning...')
    cv_res  <- mlc$cv(mldata$train_, alg_set)
    print_info('Best hyperparameters...')
    best_hp <- mlc$best_tuning(cv_res, alg_set)
    print_info('Training with obtained HP on whole training set...')
    model   <- mlc$train(mldata$train_, alg_set, best_hp)
  }else{
    cv_res  <- NULL
    best_hp <- NULL
    print_info('Training with default parameters...')
    model   <- mlc$train(mldata$train_, alg_set)
  }

  ml_result <- list(
    cv_set  = cv_set,
    alg_set = alg_set,
    cv_res  = cv_res,
    best_hp = best_hp,
    model = model
  )
  return(ml_result)
}


ml_pipeline_test <- function(mldata, seed, cv_set, model, cl_thresholds){

  # Declare classifier
  mlc <- MLClassifier$new(cv_set, seed, cl_thresholds)

  print_info('Prediction on testing data...')
  # Predict and compute all metrics
  pred           <- mlc$classify(model,mldata$test_) %>% as.mlresult()

  print_info('Binarization...')
  ml_pred_binary <- mlc$binarize_result(pred)
  attr(pred, 'classes') <- ml_pred_binary[,CRIS_CLASSES]
  print_info('Metrics...')
  metrics        <- mlc$metrics(mldr_ref = mldata$test_,
                                df_ref   = mldata$test_ref,
                                binary_pred = pred)

  ml_result <- list(
    pred = pred,
    metrics = metrics,
    thresholds = mlc$cl_thresholds,
    binary_res = ml_pred_binary
  )

  return(ml_result)
}


ml_pipeline_thresholds <- function(mldata, seed, cv_set, alg_set, model, png_path){

  # Declare classifier
  mlc <- MLClassifier$new(cv_set, seed)

  print_info('Best threshold with ROC curves...')
  thresholds <- mlc$roc_thresholds(model, mldata, png_path)

  print_info('Setting classifier thresholds...')
  mlc$cl_thresholds <- thresholds

  print_info('Prediction on training data...')
  pred <- mlc$classify(model, mldata$train_)

  print_info('Binarization...')
  ml_pred_binary <- mlc$binarize_result(pred)

  print_info('Metrics...')
  attr(pred, 'classes') <- ml_pred_binary[,CRIS_CLASSES]
  
  metrics <- mlc$metrics(mldr_ref = mldata$train_,
                         df_ref   = mldata$train_ref,
                         binary_pred = pred)

  ml_result <- list(
    pred = pred,
    metrics = metrics,
    thresholds = mlc$cl_thresholds,
    binary_res = ml_pred_binary
  )

  return(ml_result)

}

# SL classifier -----------------------------------------------------------


sl_pipeline_train <- function(sldata, method, sl_formula, seed,
                        tune = FALSE, fit_control = NULL, tune_grid = NULL){

  # Declare classifier
  slc <- SLClassifier$new(method,seed)

  # Training and (if required) tuning
  if (tune) {
    print_info('Tuning hyperparameters...')
    model  <- slc$cv(sl_formula, sldata$train_, fit_control, tune_grid)
  }else{
    print_info('Training with default hyperparameters...')
    model  <- slc$train(sl_formula, sldata$train_)
  }

  # Return result
  sl_result <- list(fit_control = fit_control,
                    tune_grid   = tune_grid,
                    cv_res      = model$results,
                    best_hp     = model$bestTune,
                    model       = model)

  return(sl_result)
}


sl_pipeline_test <- function(sldata, method, seed, model, cl_thresholds = NULL, mldata = NULL, max_cls = NULL, min_cls = NULL){

  # Declare classifier
  slc <- SLClassifier$new(method,seed, cl_thresholds)

  # Predict and binarize result
  print_info('Predicting on testing data...')
  sl_pred        <- slc$classify(model, sldata$test_)

  if (!check_type(mldata, 'PercentileMLData')){

    print_info('Binarizing results (for SL interpretation)...')

    sl_pred_binary <- slc$binarize_result_as_sl(sl_pred)

    print_info('Metrics...')
    metr <- slc$metrics(sldata$test_ref,sl_pred)

  }else{

    if (is.null(max_cls) | is.null(min_cls))
      stop('Missing max and min values for normalization into 0-1')

    # Multi label algorithm adaptation requires
    pred <- ohenery::inv_smax(as.matrix(sl_pred[,CRIS_CLASSES]))

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
    pred <- pred %>% as.data.frame() %>% rownames_to_column(ALIQUOT_LABEL) %>%
      inner_join(sl_pred %>% rownames_to_column(ALIQUOT_LABEL) %>% select(ALIQUOT_LABEL, CLASS_LABEL), by = ALIQUOT_LABEL)

    print_info('Binarizing results (for ML interpretation)...')
    sl_pred_binary <- slc$binarize_result_as_ml(pred %>% column_to_rownames(ALIQUOT_LABEL))
    
    print_info('Metrics...')
    mlm <- MLMetrics$new()
    dfres <- sl_pred_binary[,CRIS_CLASSES]
    mlres <- as.mlresult(sl_pred[,CRIS_CLASSES])
    attr(mlres, 'classes') <- dfres

    ml_conf <- list(conf = mlm$get_conf_mat(mldata$test_, mlres))
    ml_conf[['table']] <-  mlm$print_confusion_matrix(ml_conf$conf)

    metr <- list(
      conf    = ml_conf,
      global  = mlm$global_metrics(ref = mldata$test_,pred = mlres),
      local   = mlm$local_metrics(ref = mldata$test_,pred = mlres),
      relaxed = mlm$relaxed_metrics(ref = mldata$test_ref,pred = attr(mlres, 'classes'))#,
      #sim     = mlm$similarity_metrics(ref = mldata$test_ref,pred = sl_pred_binary)
    )

  }

  sl_result <- list(model   = model,
                    pred    = sl_pred,
                    metrics = metr,
                    thresholds  = slc$cl_thresholds,
                    binary_res = sl_pred_binary
                    )

  return(sl_result)

}


sl_pipeline_thresholds <- function(sldata, method, seed, model, png_path, mldata){#, roc_mldata = FALSE){

  # Declare classifier
  slc <- SLClassifier$new(method,seed)

  min_cls <- c()
  max_cls <- c()

  print_info('ROC curves and best thresholds...')
  # if (check_type(roc_mldata, 'logical',1,1) & roc_mldata & check_type(mldata, 'PercentileMLData')){
  # 
  #   pred <- slc$classify(model, mldata$train_$dataset)
  #   pred <- ohenery::inv_smax(as.matrix(pred[,CRIS_CLASSES]))
  # 
  #   for (c in CRIS_CLASSES){
  #     min_cls[c] <- min(pred[!is.infinite(pred[,c]), c], na.rm = TRUE)
  #     max_cls[c] <- max(pred[!is.infinite(pred[,c]), c], na.rm = TRUE)
  #   }
  # 
  #   thresholds <- slc$roc_thresholds(model, mldata, png_path, max_cls = max_cls, min_cls = min_cls)
  # }else{

  pred <- slc$classify(model, sldata$train_)
  pred <- ohenery::inv_smax(as.matrix(pred[,CRIS_CLASSES]))
  for (c in CRIS_CLASSES){
    min_cls[c] <- min(pred[!is.infinite(pred[,c]), c], na.rm = TRUE)
    max_cls[c] <- max(pred[!is.infinite(pred[,c]), c], na.rm = TRUE)
  }

  thresholds <- slc$roc_thresholds(model, sldata, png_path, max_cls = max_cls, min_cls = min_cls)
  # }

  print_info('Setting classifier threhsolds...')
  slc$cl_thresholds <- thresholds

  print_info('Predicting on training data...')
  sl_pred <- slc$classify(model, sldata$train_) %>% as.data.frame()

  print_info('Binarizing result (for ML interpretaion)...')
  sl_pred_binary <- slc$binarize_result_as_ml(sl_pred)

  print_info('ML Metrics...')
  mlm <- MLMetrics$new()
  dfres <- sl_pred_binary[,CRIS_CLASSES]
  mlres <- as.mlresult(sl_pred[,CRIS_CLASSES])
  attr(mlres, 'classes') <- dfres

  ml_conf <- list(conf = mlm$get_conf_mat(mldata$train_, mlres))
  ml_conf[['table']] <-  mlm$print_confusion_matrix(ml_conf$conf)

  metr <- list(
    conf    = ml_conf,
    global  = mlm$global_metrics(ref = mldata$train_,pred = mlres),
    local   = mlm$local_metrics(ref = mldata$train_,pred = mlres),
    relaxed = mlm$relaxed_metrics(ref = mldata$train_ref,pred = attr(mlres, 'classes'))#,
      #sim     = mlm$similarity_metrics(ref = mldata$test_ref,pred = sl_pred_binary)
  )

  sl_result <- list(pred    = sl_pred,
                    max_cls = max_cls,
                    min_cls = min_cls,
                    metrics = metr,
                    thresholds  = slc$cl_thresholds,
                    binary_res  = sl_pred_binary)

  return(sl_result)

}

# TSP classifier ----------------------------------------------------------

#' Apply the TSP classifier on using the given CRISData object
#' 
#' @param crisdata CRISData object with data and reference
#' @return a list with TSP prediction and computed metrics
tsp_pipeline <- function(crisdata){

  tsp <- TSPClassifier$new()
  tsp_pred <- tsp$classify(crisdata$data)
  tsp_metr <- tsp$metrics(crisdata$ref, tsp_pred)

  # Binarization
  binary_res <- tsp_pred$cl 
  binary_res[,CLASS_LABEL] <- as.character(binary_res[,CLASS_LABEL])
  for (c in levels(F_CRIS_CLASSES)){
    binary_res[binary_res[,CLASS_LABEL] == c,c] <- 1
    binary_res[binary_res[,CLASS_LABEL] != c,c] <- 0
  }
  binary_res <- rbind(binary_res, tsp_pred$uncl)
  binary_res[rownames(tsp_pred$uncl), CRIS_CLASSES] <- 0
  binary_res[rownames(tsp_pred$uncl), CLASS_LABEL] <- ''
  
  tsp_result <- list(
    pred = rbind(tsp_pred$cl, tsp_pred$uncl) %>% rownames_to_column(ALIQUOT_LABEL),
    ref  = crisdata$ref,
    metrics = tsp_metr,
    binary_res = binary_res %>% as.data.frame()
  )
  
  return(tsp_result)
}



# # Comparing pipelines results ---------------------------------------------
# 
# compare_metr <- function(res, type = 'gl'){
# 
#   comp <- data.frame()
#   for(i in seq(length(res))) {
# 
#     set_n     <- names(res)[i]
# 
#     # Global metrics (accuracy, kappa)
#     if (type == 'gl'){
#       metr <- res[[set_n]]$metrics$relaxed
#     }
# 
#     else if (type == 'sl'){
#       metr <- res[[set_n]]$metrics$global$Value
#       names(metr) <- res[[set_n]]$metrics$global$Metric
#     }
# 
#     # Similarity metrics (pearson, spearman scores correlation)
#     else if (type == 'sim') {
#       metr <- res[[set_n]]$metrics$sim$table$Value
#       names(metr) <- res[[set_n]]$metrics$sim$table$Metric
#     }
# 
#     # Multilabel metrics
#     else if (type == 'ml'){
#       metr <- res[[set_n]]$metrics$global$Value
#       names(metr) <- res[[set_n]]$metrics$global$Metric
#     }else{
#       stop('Supported types: gl, sim, ml')
#     }
# 
#     if (any(dim(comp) == 0)){
#       comp <- rownames(metr)
#       comp <- as.data.frame(t(metr))
#     }else{
#       comp <- rbind(comp, t(metr))
#     }
# 
#     rownames(comp)[i] <- set_n
# 
#   }
# 
#   return(comp)
# }
# 
# compare_local_metr <- function(res){
# 
#   loc_metr <- SLMetrics$new()$list_local_metrics() %>% tolower()
# 
#   # Result matrix
#   overall_comp <- data.frame()
#   classes  <- CRIS_CLASSES
# 
#   for (m in loc_metr){
# 
#     # Matrix for current metric (classes x methods)
#     loc_comp <- matrix(0, nrow = length(classes), ncol = length(res))
#     rownames(loc_comp) <- classes
#     colnames(loc_comp) <- names(res)
#     # For each method, get local metrics
# 
#     for (n in names(res)){
#       metr_vals <- res[[n]]$metrics$local %>% column_to_rownames('Metric')
#       loc_comp[classes,n] <- metr_vals[m, ] %>% t()
#     }
# 
#     # Add column for class and metric
#     loc_comp <- as.data.frame(loc_comp) %>%
#                 rownames_to_column('Class')
#     loc_comp <- cbind(Metric = m, loc_comp)
# 
#     # Add metric results to comparison matrix
#     overall_comp <- rbind(overall_comp, loc_comp)
#   }
# 
# 
#   return(overall_comp)
# }
# 
# plot_local_metr_comp <- function(comp, type){
# 
#   metrics <- unique(comp$metric)
#   plots   <- list()
#   for (m in metrics){
#     data <- comp %>% filter(metric == m) %>% as.data.frame()
#     plot <- ggplot(data, aes(fill=class, y=score, x=alg)) +
#             geom_bar(position="dodge", stat="identity") +
#             labs(title = paste(m, ' (',type,')', sep = '')) +
#             xlab('algorithm') +
#             ylab('value')
#     plots[[m]] <- ggplotly(plot)
#     print(plots[[m]])
#   }
# 
#   return(plots)
# 
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
