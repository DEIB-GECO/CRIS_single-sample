# Description --------------------------------------------------------------

# Exposes methods to adapt read data and transform them into a common format:
# 1. aliquot_id: TCGA-XX-XXXX-XXX-XXX-XXXX-XX (28 characters)
# 2. sample_id:  TCGA-XX-XXXX-XXX (16 characters)
# 3. patient_id: TCGA-XX-XXXX (12 characters)
# 4. tumoral_flag: TRUE if the sample is tumoral, FALSE otherwise

# The main functions used to adapt the data are
#   adapt_gmql_grch38(dataset)
#   adapt_candiolo_tcga(dataset)

# Exposes methods to filter the data before the classification

# Dependencies -------------------------------------------------------------

# Strings and regular expression
library(stringr)

# Handling of relative paths
library(here)

# Utilities (paths, print, constants etc.)
source(here('src', 'utils', 'source_utils.r'))

# ADAPT DATA --------------------------------------------------------------


#' Create an adapted dataframe from the metadata of GRCH38 (GQML)
#' 
#' @param dataset_metadata  The GMQL_GRC38 dataset (metadata list)
#' @return         Dataset in adapted form 
adapt_gmql_grch38 <- function(dataset_metadata){
  
  if (is.null(dataset_metadata) | any(dim(dataset_metadata) == 0)) {
    stop("Invalid input for adaptation of GMQL")
    return(data.frame())
  }
  
  aliquot_id   <- str_to_upper(dataset_metadata[,ALIQUOT_META]);
  patient_id   <- str_to_upper(dataset_metadata[,PATIENT_META]);
  sample_id    <- str_to_upper(dataset_metadata[,SAMPLE_META]);
  tumoral_flag <- .build_tumoral_flag(sample_id)
  
  #platform <- dataset[,"gdc__input_files__platform"];
  #normalization <- dataset[,"gcm_curated__pipeline"];
  
  adapted <- data.frame(aliquot_id, patient_id, sample_id, tumoral_flag)
  # adapted <- .addCodeMeaning(adapted, load_file(get_source_path("CODE_REF"), 
  #                                               sheet = SHEETS_CODE_REF))
  return(adapted)
  
}




#' Create an adapted dataframe from the expr matrix of TCGA data provided by Candiolo
#' 
#' @param expr_matrix  The expr matrix to adapt
#' @return             Dataset in adapted form 
adapt_candiolo_tcga <- function(expr_matrix){
  
  if (is.null(expr_matrix) | any(dim(expr_matrix) == 0)) {
    stop("Invalid input for adaptation of TCGA")
    return(data.frame())
  }
  
  names         <- colnames(expr_matrix)
  aliquot_id    <- extract_tcga_id(names, ALIQUOT_LABEL_LENGTH)
  patient_id    <- extract_tcga_id(names, PATIENT_LABEL_LENGTH)
  sample_id     <- extract_tcga_id(names, SAMPLE_LABEL_LENGTH)
  tumoral_flag  <- .build_tumoral_flag(sample_id)
  
  adapted <- data.frame(aliquot_id, patient_id, sample_id, tumoral_flag)
  # adapted <- .addCodeMeaning(adapted, load_file(get_source_path("CODE_REF"),
  #                                               sheet = SHEETS_CODE_REF))
  return(adapted)
}


#' Create an adapted dataframe from the expr matrix of PDX data provided by Candiolo
#' 
#' @param expr_matrix  The expr matrix to adapt
#' @return             Dataset in adapted form 
adapt_pdx <- function(expr_matrix){
  
  if (is.null(expr_matrix) | any(dim(expr_matrix) == 0)) {
    stop("Invalid input for adaptation pf PDX")
    return(data.frame())
  }
  # complete id model: CRC xxxx yyy tt h
  aliquot_id    <- colnames(expr_matrix)
  
  # extract the xxxx part
  patient_id    <- str_sub(aliquot_id, 4, 7)
  
  # do not distinguish between aliquot and sample
  sample_id     <- aliquot_id
  
  # samples are always tumoral
  tumoral_flag  <- TRUE
  
  # source site
  source_site <- 'Istituto di Candiolo - IRCCS'
  
  # tumor type
  tumor_type <- str_sub(aliquot_id, 8, 10)
  
  # study name
  study_name <- str_sub(aliquot_id, 1, 3)
  
  # analyte
  analyte <- 'RNA'
  
  # center
  center <- NA
  adapted <- data.frame(aliquot_id, 
                        patient_id, 
                        sample_id, 
                        tumoral_flag, 
                        source_site,
                        study_name,
                        tumor_type, 
                        analyte, 
                        center)
  
  colnames(adapted) <-  c(ALIQUOT_LABEL,
                          PATIENT_LABEL,
                          SAMPLE_LABEL,
                          TUMORAL_FLAG,
                          SOURCE_SITE_LABEL,
                          STUDY_NAME_LABEL,
                          TUMOR_TYPE_LABEL,
                          ANALYTE_LABEL,
                          CENTER_LABEL)
  return(adapted)
}


# SUPPORT FUNCTIONS (PRIVATE) ---------------------------------------------




#' Add meaning of portions of TCGA ids, given the adapted dataset and the
#' reference tables
#' 
#' @param adapted_dataset  adapted dataset with aliquot, sample, patient id 
#'                         and tumoral flag
#' @param code_ref         reference of id portion meanings 
#' @return                 Dataset in complete adapted form 
.addCodeMeaning <- function(adapted_dataset, code_ref){

  # Add tissue source site and study name (colon/rectum adenocarcinoma)
  adapted_dataset <- .addCodeColumn(adapted_dataset, 
                                   new_col_names = c(SOURCE_SITE_LABEL, STUDY_NAME_LABEL),
                                   id_portion_start = 6,
                                   id_portion_end = 7,
                                   code_ref_df = code_ref[[SHEETS_CODE_REF[1]]],
                                   ref_cols = c(2,3))
  
  # Add tumor type acronym
  adapted_dataset <- .addCodeColumn(adapted_dataset,
                                   new_col_names = TUMOR_TYPE_LABEL,
                                   id_portion_start = 14,
                                   id_portion_end = 15,
                                   code_ref_df = code_ref[[SHEETS_CODE_REF[2]]],
                                   ref_cols = 3,
                                   is_numeric_code = TRUE)

  # Check the length of the ID, before adding the last two columns
  aliquot_length <- str_length(adapted_dataset[,ALIQUOT_LABEL])
  if (length(which(aliquot_length == ALIQUOT_LABEL_LENGTH)) > 0){
    
    # Add analyte meaning
    adapted_dataset <- .addCodeColumn(adapted_dataset,
                                     new_col_names = ANALYTE_LABEL,
                                     id_portion_start = 20,
                                     id_portion_end = 20,
                                     code_ref_df = code_ref[[SHEETS_CODE_REF[3]]],
                                     ref_cols = 2)
  
    # Add center meaning
    adapted_dataset <- .addCodeColumn(adapted_dataset,
                                     new_col_names = CENTER_LABEL,
                                     id_portion_start = 27,
                                     id_portion_end = 28,
                                     code_ref_df = code_ref[[SHEETS_CODE_REF[4]]],
                                     ref_cols = 5,
                                     is_numeric_code = TRUE)
  }else{
  # Add NA columns if the length of ids is not sufficient  
    adapted_dataset <- cbind(adapted_dataset, NA, NA)
    n <- ncol(adapted_dataset)
    colnames(adapted_dataset)[(n - 1):n] <- c(ANALYTE_LABEL, CENTER_LABEL)
  }
  
  return(adapted_dataset)
}





#' Add meaning of a single portions of TCGA ids
#' 
#' @param adapted_dataset  adapted dataset with aliquot, sample, patient id 
#'                         and tumoral flag
#' @param new_col_names    names of the columns that are added
#' @param id_portion_start the position at which the portion to decode starts
#' @param id_portion_end   the position at which the portion to decode end
#' @param code_ref_df      reference of id portion meanings
#' @param ref_cols         columns of the reference to consider
#' @param is_numeric_code  flag that specifies if the portion of id to translate
#'                         is numeric or not
#' @return                 Dataset in complete adapted form 
.addCodeColumn <- function(adapted_dataset, new_col_names, 
                          id_portion_start, id_portion_end,
                          code_ref_df, ref_cols, is_numeric_code = FALSE){
  
  # Add empty columns with given colnames to the dataset
  for (i in 1:length(new_col_names)) {
    adapted_dataset <- cbind(adapted_dataset, "")
    colnames(adapted_dataset)[ncol(adapted_dataset)] <- new_col_names[i]
  }
  
  # For each aliquot 
  for (i in 1:nrow(adapted_dataset)) {
    # Extract the requested id portion and cast it
    code <- str_sub(adapted_dataset[i,ALIQUOT_LABEL], 
                    start = id_portion_start, end = id_portion_end)
    if (is_numeric_code && !is.na(as.numeric(code)))
      code <- as.numeric(code)
    else
      code <- as.character(code)
    
    # Search the row with the given code and bind the requested columns
    adapted_dataset[i, new_col_names] <- as.character(code_ref_df[code_ref_df[,1] == code,ref_cols])
  }
  
  return(adapted_dataset)
}





#' For each sample in the given samples vector, create a flag set to TRUE 
#' if the sample is tumoral, set to FALSE otherwise.
#' 
#' @param samples_list  List of samples on which building the flag
#' @return A logical vector of same length of samples_list, or with "-"
#'         if the samples list is null.
.build_tumoral_flag <- function(samples_list){
  
  if (is.null(samples_list)){
    stop("Provide a non-null samples vector")
    return(replicate(length(samples_list), "-"))
  }
  # Regular expression for tumoral sample is TCGA-XX-XXXX-0AB
  # X = alphanumeric character, A =  number, B = uppercase letter).
  # $ is the end of the string
  regexp_tumoral <- paste("TCGA[-]",               #project
                          "[A-Z 0-9]{2}[-]",       #TSS
                          "[A-Z 0-9]{4}[-]",       #participant
                          "0{1}[0-9]{1}[A-Z]{1}$", #tumoral sample (start is 0)
                          sep = "")
  
  # Grepl returns a logical value for each sample in the list if it
  # matches the given regexp
  tumoral <- grepl(regexp_tumoral, str_to_upper(samples_list), fixed = FALSE)
  
  return(tumoral)
  
}





#' Adjust the labels of the classes so that they coincide with the 
#' defined factor
#' 
#' @param prediction_vector   The vector with labels to be adjusted
#' @return the prediction vector transformed into a factor with correct
#'         class labeling
.adjust_class_label <- function(prediction_vector){
  
  if (!is.character(prediction_vector)) {
    stop("Please put a vector of character label as input.")
    return(predictionVector);
  }
  
  if (is.null(prediction_vector) | length(prediction_vector) == 0) {
    warning("Null or empty vector.")
    return(c())
  }
  
  for (i in 1:length(prediction_vector)) {
    prediction_vector[i] <- switch( as.character(prediction_vector[i]),
                                   "CRISA" =  as.character(F_CRIS_CLASSES[1]),
                                   "CRISB" =  as.character(F_CRIS_CLASSES[2]),
                                   "CRISC" =  as.character(F_CRIS_CLASSES[3]),
                                   "CRISD" =  as.character(F_CRIS_CLASSES[4]),
                                   "CRISE" =  as.character(F_CRIS_CLASSES[5]),
                                   "NA" = NA
    )
  }
  print_success("Adaptation of class labels: done.")
  return(factor(prediction_vector, levels = levels(F_CRIS_CLASSES)))
}

















