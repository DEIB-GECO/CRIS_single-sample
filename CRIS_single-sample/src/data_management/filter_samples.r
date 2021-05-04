# Description -------------------------------------------------------------

# Superclass of all filters

# Dependencies ------------------------------------------------------------

library(R6)

# Class definition --------------------------------------------------------

SamplesFilter  <- R6Class('SamplesFilter',
     inherit      = Filter,
     lock_objects = FALSE,
     lock_class   = TRUE,
     
     public = list(
       
       #' Select only the samples that appear in the given samples list.
       #' 
       #' @param expr_matrix   expression matrix to be filtered. Samples
       #'                      are column names
       #' @param list          list of samples that can be selected (character vector)
       #' @return expression matrix with reduced number of columns, corresponding
       #'         to the selected samples.
       filter_by_list = function(expr_matrix, list){
         
         super$filter_by_list(expr_matrix, list)
        
         samples_list <- list %>% unlist() %>% as.character()
        
         # Find samples of matrix that appear in the list
         dataset_samples  <- data.frame(sample = colnames(expr_matrix))
         selected_samples <- dataset_samples %>% filter(sample %in% samples_list) %>%
           unlist() %>% as.character()
        
         # Apply filter on the matrix
         expr_matrix <- expr_matrix[ ,sort(selected_samples)]
        
         # Print result info
         print_info(paste("Selected samples:", ncol(expr_matrix)))
        
         # Return expression of requested genes
         return(expr_matrix)
         
       },
       
       
       remove_by_list =  function(expr_matrix, list){
         # Check parameters
         super$filter_by_list(expr_matrix, list)
         
         # Extract matrix samples
         dataset_samples  <- colnames(expr_matrix)
         selected_samples <- vector("logical", length(dataset_samples))
         
         # Find samples of matrix that appear in the list
         for (i in seq(length(dataset_samples))) 
           selected_samples[i] <- ! dataset_samples[i] %in% list
         
         # Print result info
         print_info(paste("Selected samples:", ncol(expr_matrix[,selected_samples])))
         
         # Return expression of requested samples
         return(expr_matrix[,selected_samples])
       }
     )
                   
)
