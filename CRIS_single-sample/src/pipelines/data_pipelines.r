
# CRIS data ---------------------------------------------------------------

cris_no_ref_data_pipeline <- function(expset, log2_tr = FALSE){

  # Apply log2 transformation if required
  if (log2_tr){
    print_info('Log2 (1 + expr) applied')
    expset <- ExpressionSet(assayData = log2(1 + expset@assayData$exprs))
  }

  # Create crisdata
  crisdata <- CRISData$new(expset, NULL)

  return(crisdata)

}


cris_data_pipeline <- function(expset, ref, conf = 'all', log2_tr = FALSE){

  # Extract (non) confident samples only, if required
  if (any(conf == c('conf','non_conf'))) {

    sf  <- SamplesFilter$new()

    if (conf == 'conf'){
      print_info('Confident samples only')
      ref      <- subset(ref, ref$BH.FDR < BH_FDR_THRESHOLD)
    }else{
      print_info('Non confident samples only')
      ref      <- subset(ref, ref$BH.FDR >= BH_FDR_THRESHOLD)
    }

    filt_exp <- sf$filter_by_list(expset@assayData$exprs, ref[,ALIQUOT_LABEL])
    expset   <- ExpressionSet(assayData = filt_exp)

  }else if (conf == 'all'){
    print_info('All samples')
  }else{
    stop('Available values for confidence: conf, non_conf, all')
  }

  # Apply log2 transformation if required
  if (log2_tr){
    print_info('Log2 (1 + expr) applied')
    expset <- ExpressionSet(assayData = log2(1 + expset@assayData$exprs))
  }

  # Aliquot ids put as rownames of reference
  rownames(ref) <- ref[ ,ALIQUOT_LABEL]
  ref <- ref[colnames(expset@assayData$exprs), ] %>%
    select(ALIQUOT_LABEL, CLASS_LABEL, CLASS_DISTANCE_LABEL, CLASS_FDR_LABEL)

  ref[,CLASS_LABEL] <- factor(ref[,CLASS_LABEL], levels = CRIS_CLASSES)

  # Create crisdata
  crisdata <- CRISData$new(expset, ref)

  return(crisdata)
}


# ML data -----------------------------------------------------------------


ml_data_pipeline <- function(expset, ref, dist_thr, conf = 'all', train_samples = c(), log2_tr = TRUE){

  crisdata <- cris_data_pipeline(expset, ref, conf, log2_tr)
  # Create mldata
  mldata   <- PercentileMLData$new(crisdata, CRIS_CLASSES, dist_thr)
  mldata$split_by_list(train_samples)

  return(mldata)
}



# SL data -----------------------------------------------------------------

sl_data_pipeline <- function(expset, ref, conf = 'all', train_samples = c(), log2_tr = TRUE){

  crisdata <- cris_data_pipeline(expset, ref, conf, log2_tr)
  # Create sldata
  sldata   <- SLData$new(crisdata, CRIS_CLASSES)

  sldata$split_by_list(train_samples)

  return(sldata)
}

