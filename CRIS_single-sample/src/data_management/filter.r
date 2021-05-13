# Description -------------------------------------------------------------

# Superclass of all filters

# Dependencies ------------------------------------------------------------

library(R6)

# Class definition --------------------------------------------------------

Filter  <- R6Class('Filter',              
  lock_objects = FALSE,
  lock_class   = TRUE,
  
  public = list(
    
    #' Extract the samples on the given dataset that satisfy the given condition
    #' 
    #' @param adapted_dataset  A dataset in the common format on which the samples 
    #'                         are counted
    #' @param feature          name of the feature to search for
    #' @param feature_value    value of the feature to search for
    filter_by_attribute = function(data, feature_name, feature_value){
      
      # Check data
      if (is.null(data) | is.null(colnames(data))) {
        stop("Provide non-null and non-NA data to be filtered by feature.")
        return(NULL)
      }

      feature_name <- as.character(feature_name)

      # Check feature name
      stopifnot(check_type(feature_name, 'character',1,1))

      if (length(feature_name) != 1 | !feature_name %in% colnames(data)) {
        stop("Provide a single valid feature name.")
        return(NULL)
      }
      
      # Convert to character and compare
      data[,feature_name] <- as.character(data[,feature_name])
      feature_value <- as.character(feature_value)
      return(subset(data, data[,feature_name] %in% feature_value))
    },
    
    
    
    
    
    #' Select only the samples that appear in the given list.
    #' 
    #' @param expr_matrix   expression matrix to be filtered. Samples
    #'                      are column names
    #' @param list          list of elements that can be selected (character vector)
    #' @return expression matrix with reduced number of columns/rows, corresponding
    #'         to the selected elements
    filter_by_list = function(expr_matrix, list){
      # Define filter in class extensions
      
      # Check parameters
      stopifnot(check_type(expr_matrix,'matrix',c(1,1)), check_type(list,'character',0))

      if (any(check_type(rownames(as.matrix(expr_matrix)),'null') | check_type(colnames(as.matrix(expr_matrix)),'null'))) {
        warning("The matrix must have rownames and colnames.")
      }
      
      # Define filter in class extensions
      print_info("Filter by list")
      
    }
    
    
    
  )

)
