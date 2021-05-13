# Description -------------------------------------------------------------

# Class to compare datasets/results of classification

# Dependencies ------------------------------------------------------------

library(R6)

# Class definition --------------------------------------------------------

Comparator  <- R6Class('Comparator',              
      lock_objects = FALSE,
      lock_class   = TRUE,
      
      private = list(
        
        
        #' Checks if an id is present in the dataset
        #' @param id                The id to search
        #' @param adapted_dataset   A dataset in the common format
        #' @param id_type           Type of id (sample_id, aliquot_id, patient_id), factor value
        .is_id_in_dataset = function(id, adapted_dataset, id_type){
          
          if (!check_type(id,'character',1,1) | check_type(adapted_dataset,'null') | !check_type(id_type,'factor')) {
            stop("Provide valid input paramters. Returning FALSE.")
            return(FALSE)
          }
          
          if (!id_type %in% F_ID_TYPES) {
            stop("Id type must be chosen inside the factor F_ID_TYPES")
            return(FALSE)
          }
          
          return(id %in% adapted_dataset[,id_type])
          
        },
        
        
        
        
        
        #' Performs the intersection on two adapted datasets (PRIVATE)
        #' 
        #' @param adapted_dataset_1   A dataset in the common format
        #' @param adapted_dataset_2   A dataset in the common format
        #' @param id_type             Type of id (sample_id, aliquot_id, patient_id), factor value
        .find_intersection = function(adapted_dataset_1, adapted_dataset_2, id_type) {
          
          difference <- data.frame()
          samples <- unique(adapted_dataset_1[, id_type])
          
          for (i in seq(length(samples))) {
            if (private$.is_id_in_dataset(samples[i], adapted_dataset_2, id_type)) 
              difference <- rbind(difference, adapted_dataset_1[i,])
          }
          
          return(difference)
        },
        
        
        
        
        
        #' Performs the difference of two adapted datasets (PRIVATE)
        #' 
        #' @param adapted_dataset_1   A dataset in the common format
        #' @param adapted_dataset_2   A dataset in the common format
        #' @param id_type             Type of id (sample_id, aliquot_id, patient_id), factor value
        .find_difference = function(adapted_dataset_1, adapted_dataset_2, id_type) {
          
          difference <- data.frame()
          samples <- unique(adapted_dataset_1[, id_type])
          
          for (i in seq(length(samples))) {
            if (!private$.is_id_in_dataset(samples[i], adapted_dataset_2, id_type)) 
              difference <- rbind(difference, adapted_dataset_1[i,])
          }
          
          return(difference)
        }
        
      ),
      
      public = list(
        
        #' Performs a given operation (intersection/difference) on two adapted datasets
        #' 
        #' @param adapted_dataset_1   A dataset in the common format
        #' @param adapted_dataset_2   A dataset in the common format
        #' @param id_type             Type of id (sample_id, aliquot_id, patient_id) 
        #'                            on which datasets are compared,factor value
        #' @param operation           Type of operation (intersection/difference) to perform,
        #'                            factor value
        #' @return A data.frame with comparison of datasets
        compare_datasets = function(adapted_dataset_1, adapted_dataset_2, id_type, operation){
          
          stopifnot(check_type(adapted_dataset_1,'data.frame',c(1,1),c(Inf,Inf)))
          stopifnot(check_type(adapted_dataset_2,'data.frame',c(1,1),c(Inf,Inf)))
          stopifnot(id_type %in% F_ID_TYPES)
          stopifnot(operation %in% F_OPERATION)
          
          # Intersection
          if (operation == F_OPERATION[1])   
            return(private$.find_intersection(adapted_dataset_1, adapted_dataset_2, id_type))
          else                            
          # Operation  
            return(private$.find_difference(adapted_dataset_1, adapted_dataset_2, id_type)) 
          
        },
        
        
        
        #' Performs on two classification results
        #' 
        #' @param result_1   A classification result
        #' @param result_2   A classification result
        #' @return Comparison of two results. Implementation done on subclasses
        compare_results = function(result_1, result_2){
          
          stopifnot(!check_type(result_1,'null'), !check_type(result_2,'null'))
          print_info("Compare results")
        },
        
        
        
        #' Extracts the samples that differ in two classifications because of the given condition
        #' 
        #' @param compared_ntp_result   The comparison of the result of two ntp classifications
        #' @param condition             The condition to be considered 
        #' @param sensibility           Minimum sensibility (if required)
        check_comparison = function(compared_ntp_result, condition, sensibility = NULL){}
        
        
      )
)
