# Description -------------------------------------------------------------

# Set of miscellaneous functions that perform small heterogeneous tasks and that 
# can be used by other scripts.

# VALIDITY ---------------------------------------------------------------

# Validity functions checked
check_length <- function(object, length_min, length_max){
  
  if(!any(class(object) %in% c('logical','integer','double','numeric','character','vector', 'list')))
    return(FALSE)
  if (is.null(object) | is.array(object))
    return(FALSE)
  
  # Min length check
  if (any(is.null(length_min)) | any(is.na(length_min))){
    min_length_check <- TRUE
  }else if (any(is.nan(length_min))){
    min_length_check <- FALSE
  }else if (any(!is.null(length_min)) & any(!is.na(length_min)))
    min_length_check <- length(object) >= length_min
  else
    min_length_check <- FALSE

  # Max length check
  if (any(is.null(length_max)) | any(is.na(length_max)))
    max_length_check <- TRUE
  else if (any(is.nan(length_max)))
    max_length_check <- FALSE
  else if (any(!is.null(length_max)) & any(!is.na(length_max)))
    max_length_check <- length(object) <= length_max
  else
    max_length_check <- FALSE

  # Both conditions must be satisfied 
  outcome <- min_length_check & max_length_check

  if (any(is.na(outcome) | is.null(outcome) | is.nan(outcome) | is.infinite(outcome)))
    return(FALSE)
  else
    return(outcome)

}

check_dim <- function(object, dim_min, dim_max){
  
  if(is.null(object))
    return(FALSE)
  
  if(!is.data.frame(object) & !is.array(object))
    return(FALSE)
  
  # Type check skipped if all dims are null
  if (is.null(dim_min) & is.null(dim_max))
    return(TRUE)
  
  if (!is.null(dim_min) & (!is.numeric(dim_min) | length(dim_min) < 1))
    stop('min dim, if specified, must be a vector of 1 or more numbers')
  
  if (!is.null(dim_max) & (!is.numeric(dim_max) | length(dim_max) < 1))
    stop('max dim, if specified, must be a vector of 1 or more numbers')

  # Coherence of dimensions
  if (!is.null(dim_min) & length(dim(object)) != length(dim_min))
    stop('min dim and object dim must be vectors of the same length')
  
  if (!is.null(dim_max) & length(dim(object)) != length(dim_max))
    stop('max dim and object dim must be vectors of the same length')
  
  # Vector with outcome of check for each dimension
  n_dim        <- length(dim(object))
  
  length_check <- lapply(seq(n_dim), function(i){
    
    # Min check
    if (any(is.na(dim_min[i]) | is.null(dim_min[i])))
      min_check <- TRUE
    else if (any(is.nan(dim_min[i])))
      min_check <- FALSE
    else
      min_check  <- dim(object)[i] >= dim_min[i]
      
    # Max check
    if (any(is.na(dim_max[i]) | is.null(dim_max[i])))
      max_check <- TRUE
    else if (any(is.nan(dim_max[i])))
      max_check <- FALSE
    else
      max_check  <- dim(object)[i] <= dim_max[i]
      
    # Both checks must be verified
    return(min_check & max_check)
  }) %>% unlist() %>% as.logical()
  
  # Global check
  outcome <- all(length_check)
  
  if (is.na(outcome) | is.null(outcome) | is.nan(outcome) | is.infinite(outcome))
    return(FALSE)
  else
    return(outcome)
  
}

#' Check that a given object is not NULL nor NA
#' 
#' @param object An object to check
#' @param type An optional string if type checking is required
#' @return TRUE id the object is not NULL and not NA (possibly type is respected), 
#'         FALSE otherwise
check_type <- function(object, type, length_min = NULL, length_max = NULL, dim_min = NULL, dim_max = NULL){
  
  # Allowed types
  allowed_vectors    <- c('logical','integer','double','numeric','character','vector')
  allowed_list       <- c('list','data.frame')
  allowed_array      <- c('array','matrix')
  allowed_particular <- c('null','na','nan','inf','factor','ordered','environment')
  allowed_classes    <- c('ExpressionSet','Filter','GenesFilter','SamplesFilter',
                          'MLMetrics','SLMetrics','Metrics',
                          'Data','CRISData','PercentileMLData','SLData',
                          'AlgSettings','CVSettings',
                          'FeatureManager','FeatureManagerPDX','FeatureManagerTCGA',
                          'SingleSampleClassifier','mldr',
                          'MLClassifier','SLClassifier','NTPClassifier','TSPClassifier'
                          )
  
  allowed_types <- c(allowed_vectors,
                     allowed_list,
                     allowed_array,
                     allowed_particular,
                     allowed_classes)
  
  # Check type of object
  if (!is.character(type) | length(type) != 1)
    stop('type must a non null empty string')
  
  if (!type %in% allowed_types)
    stop(paste(c('Allowed types:', allowed_types), collapse = '\n'))
  
  # Check for atomic and standard types
  if (!type %in% allowed_classes){
    type_check <- switch(type,
                         "null"          = is.null(object),
                         "na"            = is.na(object),
                         "nan"           = is.nan(object),
                         "inf"           = is.infinite(object),
                         "environment"   = mode(object) == 'environment',
                         "factor"        = is.factor(object) | is.ordered(object),
                         "ordered"       = is.ordered(object),
                         "logical"       = is.logical(object),
                         "integer"       = is.integer(object),
                         "double"        = is.double(object),
                         "numeric"       = is.numeric(object),
                         "character"     = is.character(object),
                         "vector"        = is.vector(object),
                         "list"          = is.list(object),
                         "data.frame"    = is.data.frame(object),
                         "array"         = is.array(object),
                         "matrix"        = is.matrix(object),
                         
                         )
    type_check <- any(type_check == TRUE)
  # Check for class objects
  }else {
    type_check <- any(class(object) == type)
  }
  
  # Dimension check for atomic vectors and lists (excluding arrays and data.frames)
  if (all(type_check &  type %in% c(allowed_vectors, 'list') & !is.array(object))){
    type_check <- type_check & any(check_length(object, length_min, length_max) == TRUE)
  
    # Dimension check for arrays (matrices) and data.frames
  }else if (all(type_check & (is.data.frame(object) | is.array(object)))){
    type_check <- type_check & any(check_dim(object, dim_min, dim_max) == TRUE)
  }
    
  return(type_check)
  
}


#' Check if any undesired value (NA, NaN, NULL, Inf) is present in a dataframe.
#' @param df    Dataframe to check
#' 
#' @return      TRUE if any undersired value is present, FALSE otherwise.
#'  
any_uncomplete <- function(df){
  
  if (check_type(df,'null'))
    return(TRUE)
  
  data <- unlist(df, recursive = TRUE)
  
  check <- lapply(data, function(k){
    missing   <- any(is.na(k) | k == 'NA')  | any(is.null(k) | k == 'NULL')
    undefined <- any(is.nan(k) | k == 'NaN') | any(is.infinite(k) | k == 'Inf' | k == '-Inf')
    return(missing | undefined)
  }) %>% unlist()
 
  return(any(check))
}



merge_lists <- function(...){
  
  # Put arguments into list
  args <- list(...)
  merged <- list()

  for (n in seq(length(args))){
    if (check_type(args[[n]], 'list', 0)){
      for (m in names(args[[n]]))
        merged[[m]] <- args[[n]][[m]]
    }else if (!check_type(args[[n]], 'list') & !check_type(names(args)[n], 'null')){
      merged[[names(args)[n]]] <- args[[n]]
    }else{
      i <- length(merged)
      merged[[i + 1]] <- args[[n]]
    }
  }
  
  if (length(merged) > 0)
    return(simplify_list(merged))
  else
    return(merged)
}


simplify_list <- function(l){
  if (class(l) == 'list' & length(l) == 1 & class(l[[1]]) == 'list'){
    return(l[[1]])
  }
  
  if (class(l) == 'list' & length(l) >= 1){
    for (i in seq(length(l))){
      while (class(l[[i]]) == 'list' & length(l[[i]]) == 1){
          n <- ''
          if (!is.null(names(l[i])) & any(names(l[i]) != ''))
            n <- names(l[i])
          if (!is.null(names(l[[i]])) & any(names(l[[i]]) != ''))
            n <- names(l[[i]])
    
          l[i] <- l[[i]]
          names(l)[i] <- n
        }
    }
    return(l)
  }
  
  
}
