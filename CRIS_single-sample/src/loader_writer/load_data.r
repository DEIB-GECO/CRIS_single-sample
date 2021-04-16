# DescriptionS ------------------------------------------------------------

# Exposes methods to load specific datasets

# Dependencies ------------------------------------------------------------

# Base functions to load data
source(here('src', 'loader_writer', 'load_functions.r'))

# Utilities (paths, print, constants etc.)
source(here('src', 'utils', 'source_utils.r'))


# DATASETS ----------------------------------------------------------------





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




