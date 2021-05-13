library(R6)


SummaryMakerTCGA  <- R6Class('SummaryMakerTCGA',
                            
  inherit = SummaryMaker,
  lock_objects = FALSE,
  lock_class   = TRUE,
  
  ####### Public fields #######
  
  public = list(
    
    
    #' Computes the number of normal, tumoral and overall samples for each patient,
    #' other than the overall number of normal and tumoral samples in the entire
    #' dataset.
    #' 
    #' @param  adapted_dataset   one of the datasets in the common format
    #' @return Dataframe that specifies, for each patient, the number of tumoral 
    #'         samples, of normal samples and of overall samples.
    samples_summary = function(adapted_dataset){
      
      if (is.null(adapted_dataset)) {
        stop("Provide valid adapted dataset.")
        return(data.frame())
      }
      
      normal  <- self$count_by_attribute_value(adapted_dataset, TUMORAL_FLAG, 'FALSE')
      tumoral <- self$count_by_attribute_value(adapted_dataset, TUMORAL_FLAG, 'TRUE')
      
      patient_counts           <- self$count_patient_samples(adapted_dataset, TUMORAL_FLAG, c('TRUE','FALSE'))
      colnames(patient_counts) <- c(PATIENT_LABEL, 
                                    COUNT_TUMORAL_LABEL, 
                                    COUNT_NORMAL_LABEL, 
                                    COUNT_SAMPLES_LABEL)
      
      samples_distribution     <- self$summarize_samples_distribution(patient_counts)
      tumor_type_distribution  <- self$summarize_tumor_types(adapted_dataset)
      
      analysis_result <- list()
      analysis_result[[COUNT_NORMAL_LABEL]]      <- normal
      analysis_result[[COUNT_TUMORAL_LABEL]]     <- tumoral
      analysis_result[[COUNT_SAMPLES_LABEL]]     <- normal + tumoral
      analysis_result[[COUNT_PATIENTS_DF]]       <- patient_counts
      analysis_result[[SAMPLES_DISTRIB_DF]]      <- samples_distribution
      analysis_result[[TUMOR_TYPES_DISTRIB_DF]]  <- tumor_type_distribution
      
      
      return(analysis_result)
    
    },
    
    
    
    
    
    #' Create a dataframe that shows the number of samples for each type of tumor, 
    #' distinguished by rectum and colon cancer
    #' 
    #' @param adapted_dataset  The dataset (in common format) on which the count is 
    #'                         performed
    #' @return A dataframe with the count of each type of tumor, distinguished
    #' by colon and rectum adenocarcinoma
    summarize_tumor_types = function(adapted_dataset){
      
      sf <- SamplesFilter$new()
      # Samples that belong to colon or rectum carcinoma
      colon_samples  <- sf$filter_by_attribute(adapted_dataset, 
                                          STUDY_NAME_LABEL, 
                                          "Colon adenocarcinoma")
      rectum_samples <- sf$filter_by_attribute(adapted_dataset, 
                                          STUDY_NAME_LABEL, 
                                          "Rectum adenocarcinoma")
      
      tumor_types   <- levels(factor(adapted_dataset[,TUMOR_TYPE_LABEL]))
      n_tumor_types <- length(tumor_types)
      
      # Count the tumor types for both colon and rectum
      colon_count  <- vector("integer", n_tumor_types)
      rectum_count <- vector("integer", n_tumor_types)
      
      for (i in seq(n_tumor_types)) {
        colon_count[i]  <- self$count_by_attribute_value(colon_samples, TUMOR_TYPE_LABEL, tumor_types[i])
        rectum_count[i] <- self$count_by_attribute_value(rectum_samples,TUMOR_TYPE_LABEL, tumor_types[i])
      }
      
      # Create and return a result dataframe
      tumor_type_counts <- data.frame(colon  = colon_count,
                                      rectum = rectum_count,
                                      total  = colon_count + rectum_count)
      rownames(tumor_type_counts) <- tumor_types
      
      return(tumor_type_counts)
    },
    
    
    
    
    
    
    #' Create a dataframe describing, for each combination of normal and tumoral samples,
    #' the number of patients that belong to it.
    #' 
    #' https://stackoverflow.com/questions/33710240/how-to-attach-a-title-to-a-data-frame-in-r
    #' for how to add title to dataframe with comments.
    #' 
    #' @param  counts_dataset  A dataframe containing, for each sample, the count
    #'                         of normal, total and tumoral samples
    #' @return A dataframe with the counts of all possible combinations of normal
    #'         and tumoral samples.
    summarize_samples_distribution = function(counts_dataset){
      
      # The number of tumoral/normal samples (per patient) ranges from 
      # 0 to the max total number of samples
      max_samples_per_patient <- max(counts_dataset[,COUNT_SAMPLES_LABEL])
      n_cases <- max_samples_per_patient + 1
      
      # The results contains also the total, thus add an additional column
      summarized_results <- matrix(0, n_cases + 1, n_cases + 1)
      
      sf <- SamplesFilter$new()
      
      for (i in seq(n_cases + 1)) {				
        for (j in seq(n_cases)) {
          # Compute the counts of the current combination
          if (i <= n_cases) {
            # The count starts from 0, not from 1
            n_normal  <- i - 1
            n_tumoral <- j - 1
            normal_samples  <- sf$filter_by_attribute(counts_dataset, 
                                                 COUNT_NORMAL_LABEL, as.integer(n_normal))
            tumoral_samples <- sf$filter_by_attribute(counts_dataset, 
                                                 COUNT_TUMORAL_LABEL, as.integer(n_tumoral))
            
            summarized_results[i,j] <- length(intersect(normal_samples [,PATIENT_LABEL],
                                                        tumoral_samples[,PATIENT_LABEL]))
          }else{
            # Sum the values of the row (number of patients with i-1 normal samples)
            summarized_results[i,j] <- sum(summarized_results[1:n_cases,j]) 
          }
        }
        # Sum the values of the column (number of patients with j-1 tumoral samples)
        summarized_results[i,n_cases + 1] <- sum(summarized_results[i,1:n_cases]) 
      }
      
      summarized_results <- data.frame(summarized_results)
      for (j in 1:n_cases) {
        colnames(summarized_results)[j] <- paste(j - 1, "tumoral samples", sep = " ")
        rownames(summarized_results)[j] <- paste(j - 1, "normal samples",  sep = " ")
      }
      
      
      colnames(summarized_results)[n_cases + 1] <- "Total (row)"
      rownames(summarized_results)[n_cases + 1] <- "Total (column)"
      
      return(summarized_results)
      
    },
    
    
    
    
    
    
    #' Creates a dataframe that summarizes the cardinality of intersection/diffrence
    #' between all possible combinations inside a list of datasets
    #' 
    #' @param adapted_dataset_list   A list of datasets in the common format
    #' @param id_type_list           A list of type of ids (one for each dataset), list of factor values
    #' @param adapted_dataset_names  A list of dataset names (one for each dataset)
    #' @param operation              Type of operation (intersection/difference), factor value
    create_summary_for_operation = function(adapted_dataset_list,
                                            id_type_list,
                                            adapted_dataset_names,
                                            operation) {
        
      # Check of input dataset list and names
      stopifnot(check_type(adapted_dataset_list,'list',1))
      stopifnot(check_type(adapted_dataset_names,'character',length(adapted_dataset_list),length(adapted_dataset_list)))

      # Check id types
      stopifnot(check_type(id_type_list,'factor',length(adapted_dataset_list),length(adapted_dataset_list)), all(id_type_list %in% F_ID_TYPES))
    
      
      # Check operation
      stopifnot(check_type(operation,'factor',1,1), operation %in% F_OPERATION)

      # Check that each dataset is associated to an id_type and a dataset name.
      n_datasets <- length(adapted_dataset_list)
      

      # Comparison within all datasets of the list
      comp <- Comparator$new()
      summary    <- matrix(0, n_datasets, n_datasets)
      for (i in 1:n_datasets) {
        for (j in 1:n_datasets) {
          # Take the id with highest order (shorter, less specific)
          id_type    <- max(id_type_list[[i]], id_type_list[[j]])
          comparison <- comp$compare_datasets(adapted_dataset_list[[i]], 
                                              adapted_dataset_list[[j]], 
                                              id_type, 
                                              operation)
          # If the comparison is successful, a non-empty comparison is obtained  
          if (all(dim(comparison)) > 0) {
            summary[i,j] <- nrow(comparison)
          }
        }
      }
      
      
      colnames(summary) <- adapted_dataset_names
      rownames(summary) <- adapted_dataset_names
      
      return(data.frame(summary))
      
    }
    
  )
)
