# DescriptionS ------------------------------------------------------------

# Exposes methods to load specific datasets

# Dependencies ------------------------------------------------------------

# Base functions to load data
source(here('src', 'loader_writer', 'load_functions.r'))

# Utilities (paths, print, constants etc.)
source(here('src', 'utils', 'source_utils.r'))


# Datasets ----------------------------------------------------------------


# TODO: reference to Isella et. al.
#' Load the original TCGA dataset aligned to HG19 in the global environment
load_candiolo_hg19 <- function(){
  
  # Load data in global environment
  assign(x     = "CANDIOLO_HG19",
         value = load_file(path = path_loader$get_path("CANDIOLO_HG19")), 
         envir = .GlobalEnv)
  
  # Load adapted format in global environment
  # assign(x     = "adapted_CANDIOLO_HG19",
  #        value = adapt_candiolo_tcga(CANDIOLO_HG19@assayData[["exprs"]]), 
  #        envir = .GlobalEnv)
}





#' Load the TCGA dataset aligned to GRCh38 in the global environment. Logical flags help
#' to control which version must be loaded. Flags are all put to FALSE by default.
#' 
#' @param original Logical flag; load the original version of the dataset.
#' @param filtered Logical flag; load the pre-processed version of the dataset.
#' @param metadata Logical flag; load metadata and clinical annotations of the dataset
#' @param uniformed Logical flag; if both uniformed and filtered flags are true,
#' load the preprocessed version with genes uniformed to PDX_GRCh38 reference. Only genes
#' in common with PDX_GRCh38 are considered.
load_gmql_grch38 <- function(original = FALSE, filtered = FALSE, metadata = FALSE, uniformed = FALSE){
  
  # Load filtered data
  if (any(filtered == TRUE) & any(uniformed == FALSE))
    assign(x     = "GMQL_GRCH38_FILTERED",
           value = load_file(path = path_loader$get_path("GMQL_GRCH38_FILT")),
           envir = .GlobalEnv)
  else if (any(filtered == TRUE) & any(uniformed == TRUE))
    assign(x     = "GMQL_GRCH38_FILTERED",
           value = load_file(path = path_loader$get_path("GMQL_GRCH38_FILT_UNIF")),
           envir = .GlobalEnv)
  
  # Load original data
  if (any(original == TRUE))
    assign(x     = "GMQL_GRCH38_LIST",
           value = load_file(path = path_loader$get_path("GMQL_GRCH38")),
           envir = .GlobalEnv)
  
  # Load adapted format
  if (any(metadata == TRUE)){
    
    assign(x     = "GMQL_GRCH38_metadata",
           value = load_file(path = path_loader$get_path("GMQL_GRCH38_META")) %>% as.data.frame(),
           envir = .GlobalEnv)
    
    assign(x     = "GMQL_GRCH38_annot",
           value = load_file(path = path_loader$get_path("GMQL_GRCH38_ANNOT")) %>% as.data.frame(),
           envir = .GlobalEnv)
    
    # assign(x     = "adapted_GMQL_GRCH38",
    #        value = adapt_gmql_grch38(GMQL_GRCH38_metadata),
    #        envir = .GlobalEnv)
  }
  
}





#' Load the PDX batches aligned to GRCh38, possibly merged, in the global environment. Flags help
#' to control which version must be loaded. Logical flags are all put to FALSE by default.
#' 
#' @param original Logical flag; load the original version of the batches enumerated in batches parameter.
#' @param filtered Logical flag; load the pre-processed version the batches enumerated in batches parameter or
#' of the merged batches, if merged is set to TRUE.
#' @param merged Logical flag; load the pre-processed merged version the batches batches.
#' @param batches Integer vector; number of batches to load. Any combination of values from 1 to 6 is accepted.
#' By default, it is empty.
#' @param metadata Logical flag; load metadata and clinical annotations of the merged batches.
#' @param uniformed Logical flag; if both uniformed and merged flags are true,
#' load the preprocessed version with genes uniformed to PDX_GRCh38 reference. Only genes
#' in common with TCGA_GRCh38 are considered.
load_pdx <- function(original = FALSE, filtered = FALSE, 
                     batches = integer(0), merged = FALSE,
                     uniformed = FALSE, metadata = FALSE){
  
  # Check batch numbers
  if (check_type(batches, "numeric",1)){
    batches <- as.integer(batches)
  }
  
  # Load required batches, possibly filtered
  for (i in batches){
    
    if (any(original == TRUE)){
      dataset_name  <- paste("Hbiod",i,"_LMX", sep = "")
      dataset_label <- paste("PDX",i, sep = "_")
      assign(x     = dataset_name,
             value = load_file(path = path_loader$get_path(dataset_label)),
             envir = .GlobalEnv)
    }
    if (any(filtered == TRUE)){
      dataset_name  <- paste("Hbiod",i,"_LMX_FILTERED", sep = "")
      dataset_label <- paste("PDX",i, "FILT",sep = "_")
      assign(x     = dataset_name,
             value = load_file(path = path_loader$get_path(dataset_label)),
             envir = .GlobalEnv)
    }
  }
  
  # Load merged PDX, possibly filtered and uniformed with TCGA_GRCh38 genes
  if (any(merged == TRUE) & any(uniformed == FALSE)){
    dataset_name <- 'Hbiod_MERGED_FILTERED'
    assign(x     = dataset_name,
           value = load_file(path = path_loader$get_path("PDX_MERGED_FILT")),
           envir = .GlobalEnv)
    
  }else if (any(merged == TRUE) & any(uniformed == TRUE)){
    dataset_name <- 'Hbiod_MERGED_FILTERED'
    assign(x     = dataset_name,
           value = load_file(path = path_loader$get_path("PDX_MERGED_FILT_UNIF")),
           envir = .GlobalEnv)
  }
  
  # Load annotations
  if (any(metadata == TRUE)){
    dataset_name <- 'pdx_annotations'
    assign(x     = dataset_name,
           value = load_file(path = path_loader$get_path("PDX_MERGED_ANNOT")),
           envir = .GlobalEnv)
  }
  
}



#' Load the PDX batches aligned to GRCh38, merged, in the global environment, 
#' excluding the samples used in bio-driven feature selection
#' @param uniformed Logical flag; if both uniformed and merged flags are true,
load_unused_pdx <- function(uniformed = TRUE){
  
  load_pdx(merged = TRUE, uniformed = uniformed)

  # Remove samples used in bio-driven feature selection
  used_pdx_samples <-
      load_file(path_loader$get_path("HD_PDX")) %>%
      unlist() %>% 
      as.character()
    
  sf <- SamplesFilter$new()
  expset_list <- list()
  for (n in names(Hbiod_MERGED_FILTERED)){
    print_info(n)
    expr <- Hbiod_MERGED_FILTERED[[n]]@assayData$exprs
    expr <- sf$remove_by_list(expr, used_pdx_samples)
    expset_list[[n]] <- ExpressionSet(assayData = expr)
  }
  
  assign('Hbiod_MERGED_FILTERED', expset_list, .GlobalEnv)
  
  print_success("Unused PDX loaded in the environment.")
}
# References --------------------------------------------------------------

#' Load the samples that have been assigned to training
#' @return a character vector with the training samples
load_training_samples <- function(){
  
  path    <- path_loader$get_path("TCGA_SPLITTING")
  samples <- load_file(path, sheet = 'training')
  return(samples[,ALIQUOT_LABEL])
  
}

#' Load the TCGA samples that have been assigned to testing
#' @return a character vector with the testing samples
load_tcga_testing_samples <- function(){
  
  path    <- path_loader$get_path("TCGA_SPLITTING")
  samples <- load_file(path, sheet = 'testing')
  return(samples[,ALIQUOT_LABEL])
  
}


#' Load the published NTP classification results into a data.frame
load_published_ntp <- function(){
  
  # Read the published results
  data    <- load_file(path = path_loader$get_path("PUB_NTP"))
  
  if (is.null(data) | any(dim(data) == 0))
    stop("load_published_ntp: invalid read data.")
  
  # Take only the TCGA results
  adapted <- subset(data, data[,1] == "TCGA")
  
  # Remove the first column (dataset) and set the colnames
  adapted <- data.frame(adapted[,2:ncol(adapted)])
  colnames(adapted) <- colnames <- c(PATIENT_LABEL,
                                     CLASS_LABEL, 
                                     BEST_DISTANCE_LABEL,
                                     BEST_FDR_LABEL)
  
  # Adjust the class labels
  adapted[,CLASS_LABEL] <- gsub(x = adapted[,CLASS_LABEL], pattern = "CRIS", replacement = "CRIS-")

  assign(x     = "published_ntp",
         value = adapted,
         envir = .GlobalEnv)
}



#' Load the published TSP classification results
load_published_tsp <- function(){
  
  # Read the data
  data <- load_file(path = path_loader$get_path("PUB_TSP"))
  
  if (is.null(data) | any(dim(data) == 0))
    stop("load_published_tsp: invalid read data.")
  
  # Take only the TCGA results
  adapted <- subset(data, data[,2] == "TCGA")
  
  # Extract the results of TSP
  adapted  <-  data.frame(adapted[,1],  
                          adapted[,5])
  
  # Set the colnames
  colnames(adapted) <- colnames <- c(PATIENT_LABEL, 
                                     CLASS_LABEL)
  
  # Adjust the class labels
  adapted[,CLASS_LABEL] <- gsub(x = adapted[,CLASS_LABEL], pattern = "CRIS", replacement = "CRIS-")
  adapted[which(adapted[,CLASS_LABEL] == "NA"),CLASS_LABEL] <- NA
  
  # Save result into global environment
  assign(x     = "published_tsp",
         value = adapted %>% arrange(predict.label2),
         envir = .GlobalEnv)
}




# Features  ---------------------------------------------------------------

load_features_hg19 <- function(original = FALSE, filtered = FALSE){
  
  # Load original features
  if (any(original == TRUE))
    assign(x     = "features_hg19_original",
           value = load_file(path_loader$get_path("FEATURES_ORIGINAL")),
           envir = .GlobalEnv)
  
  # Load filtered features
  if (any(filtered == TRUE))
    assign(x     = "features_hg19_updated",
           value = load_file(path_loader$get_path("FEATURES_HG19")),
           envir = .GlobalEnv)
  
}


load_features_grch38 <- function(dataset = 'tcga'){
  
  if (any(dataset == 'tcga'))
    assign(x     = "features_grch38",
           value = load_file(path_loader$get_path("FEATURES_GRCH38")),
           envir = .GlobalEnv)
  
  if (any(dataset == 'pdx')) 
    assign(x     = "features_pdx",
           value = load_file(path_loader$get_path("FEATURES_PDX")),
           envir = .GlobalEnv)
  
}


# Prepared datasets -------------------------------------------------------

# OK
.get_feature_manager <- function(uniformed, dataset, fs_type, lasso_type){
  
  # Feature selection
  if (uniformed | dataset == 'pdx'){
    print_info('Feature Selection (PDX aliases)...')
    fs <- FeatureManagerPDX$new(fs_type,lasso_type)
  }else if (!uniformed & dataset == 'tcga'){
    print_info('Feature Selection (TCGA aliases)...')
    fs <- FeatureManagerTCGA$new(fs_type,lasso_type)
  }
  
  return(fs)
}

# OK
.get_training_samples <- function(load_training){
  
  if (any(load_training == TRUE)){
    print_info('Loading training samples...')
    train_samples <- load_training_samples()
  }else {
    print_info('No training sample..')
    train_samples <- c()
  }
  
}

# OK
.get_ref <- function(data_type, dataset){
  
  if (data_type == 'cris_no_ref' | !dataset %in% c('tcga','pdx')){
    print_info('No reference...')
    ref_data <- NULL
  }else if (data_type != 'cris_no_ref' & dataset %in% c('tcga','pdx')){
    
    print_info('Loading reference...')
    label <- toupper(paste('ntp','ref',dataset, sep = '_'))
    ref   <- load_file(path_loader$get_path(label))
    ref_data <- ref$result
  }

  return(ref_data)
}

# OK
.filt_samples_and_ref <- function(samples_filter, data, ref){
  
  # If filter is null, return original data
  if (is.null(samples_filter))
    return(list(data = data, ref = ref))
  
  # Extract samples from data
   if (!is.null(data)){
    filt_data <- SamplesFilter$new()$filter_by_list(data, samples_filter)
  }else{
    filt_data <- data
  }

  # Extract samples from ref
  if (!is.null(ref)){
    filt_ref <- Filter$new()$filter_by_attribute(ref, ALIQUOT_LABEL, samples_filter)
  }else{
    filt_ref <- ref
  }
  
  return(list(data = filt_data, ref = filt_ref))
}


# OK
.apply_data_pipeline <- function(data_type, filt_data, filt_ref, confident, train_samples, bin_thr){
  
  print_info(paste('Selecting pipeline:', data_type))
  
  data <- switch(data_type,
                 'ml'   = ml_data_pipeline(filt_data, filt_ref, bin_thr, confident, train_samples),
                 'sl'   = sl_data_pipeline(filt_data, filt_ref, confident, train_samples),
                 'cris' = cris_data_pipeline(filt_data, filt_ref, confident, log2_tr = FALSE),
                 'cris_no_ref' = cris_no_ref_data_pipeline(filt_data, log2_tr = FALSE)
                 )
  
  if (is.null(data))
    stop('data pipeline returned NULL: check type is either `ml` or `sl` or `cris` or `cris_no_ref`')
  
  return(data)
}

# TODO: CHECK
load_testing_tcga_tsp <- function(confident,uniformed,fs_type, lasso_type = '', type = 'train'){
  
  # Data loading
  load_gmql_grch38(original = FALSE, filtered = TRUE, uniformed)
   
  filt_data     <- GMQL_GRCH38_FILTERED$cpm@assayData$exprs
  train_samples <- load_training_samples()
  
  if (type == 'train')
    samp_filt <- train_samples
  else if (type == 'test')
    samp_filt <- dplyr::setdiff(colnames(filt_data), train_samples)
  else
    samp_filt <- NULL
  
  # Apply samples filter
  filt_res  <- .filt_samples_and_ref(samp_filt, filt_data, .get_ref('cris', 'tcga'))
  filt_data <- ExpressionSet(assayData = filt_res$data)
  filt_ref  <- filt_res$ref

  # Feature selection
  fs        <- .get_feature_manager(uniformed, 'tcga', fs_type, lasso_type)
  filt_data <- fs$feature_selection(filt_data, fs_type)
   
  return(.apply_data_pipeline('cris', filt_data, filt_ref, confident, train_samples, bin_thr))

}

# TODO: CHECK
load_prepared_tcga_data <- function(confident, uniformed, fs_type, type = 'sl', lasso_type = '', samples_filter = NULL, load_training = TRUE){
  
  # Data loading
  load_gmql_grch38(original = FALSE, filtered = TRUE, uniformed = uniformed)
  
  filt_data     <- GMQL_GRCH38_FILTERED$cpm@assayData$exprs
  filt_ref      <- .get_ref(type, 'tcga')
  train_samples <- .get_training_samples(load_training)
  bin_thr       <- load_file(path_loader$get_path('NTP_THR'))$tcga
  

  # Apply samples filter
  filt_res  <- .filt_samples_and_ref(samples_filter, filt_data, filt_ref)
  filt_data <- ExpressionSet(assayData = filt_res$data)
  filt_ref  <- filt_res$ref

  # Feature selection
  fs        <- .get_feature_manager(uniformed, 'tcga', fs_type, lasso_type)
  filt_data <- fs$feature_selection(filt_data, fs_type)
  
  
  return(.apply_data_pipeline(type, filt_data, filt_ref, confident, train_samples, bin_thr))

}

# TODO: CHECK
load_prepared_pdx_data <- function(confident, uniformed, fs_type, type = 'sl', lasso_type = '', samples_filter = NULL){
  
  # Loading data
  load_unused_pdx(uniformed)
  filt_data     <-  Hbiod_MERGED_FILTERED$cpm@assayData$exprs
  ref_data      <- .get_ref(type, 'pdx')
  train_samples <- .get_training_samples(FALSE)

  bin_thr       <- load_file(path_loader$get_path('NTP_THR'))$pdx
  
  # Apply samples filter
  filt_res  <- .filt_samples_and_ref(samples_filter, filt_data, ref_data)
  
  filt_data <- ExpressionSet(assayData = filt_res$data)
  filt_ref  <- filt_res$ref
  
  # Feature selection
  fs   <- .get_feature_manager(uniformed, 'pdx', fs_type, lasso_type)
  filt_data <- fs$feature_selection(filt_data, fs_type)
  
  return(.apply_data_pipeline(type, filt_data, filt_ref, confident, train_samples, bin_thr))
  
}
