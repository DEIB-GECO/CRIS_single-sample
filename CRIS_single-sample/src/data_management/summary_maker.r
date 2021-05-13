# Description -------------------------------------------------------------

# Class that contains functions to create tabular summaries.


# Class definition --------------------------------------------------------



SummaryMaker  <- R6Class('SummaryMaker',
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
    samples_summary = function(adapted_dataset){},
    
    
    
    #' Create a dataframe that shows the number of samples for each type of tumor, 
    #' distinguished by rectum and colon cancer
    #' 
    #' @param adapted_dataset  The dataset (in common format) on which the count is 
    #'                         performed
    #' @return A dataframe with the count of each type of tumor, distinguished
    #' by colon and rectum adenocarcinoma
    summarize_tumor_types = function(adapted_dataset){},
    
    
    
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
    summarize_samples_distribution = function(counts_dataset){},
    
    
    
    
    
    #' Creates a dataframe that summarizes the cardinality of intersection/difference
    #' between all possible combinations inside a list of datasets
    #' 
    #' @param adapted_dataset_list   A list of datasets in the common format
    #' @param id_type_list           A list of type of ids (one for each dataset), list of factor values
    #' @param adapted_dataset_names  A list of dataset names (one for each dataset)
    #' @param operation              Type of operation (intersection/difference), factor value
    create_summary_for_operation =  function(adapted_dataset_list, 
                                             id_type_list, 
                                             adapted_dataset_names, 
                                             operation) {},
    
    
    
    #' Creates a dataframe with count of samples with the given attribute, grouped by
    #' the attribute values. If not specified, the values are taken from the data.
    #' 
    #' @param data       Data to count
    #' @param attribute  Column on which group data and then count
    #' @param values     The values to be used (if any)
    #' @return A dataframe with 1 column and a row for the count of each attribute value,
    #'         plus a row that sums all the others.
    count_by_attribute = function(data, attribute, values = NULL){
      
      # Check input
      stopifnot(check_type(data,'data.frame',c(1,1)))
      stopifnot(check_type(attribute, 'character',1,Inf))
      stopifnot(attribute %in% colnames(data))
      
      # Identify type of samples in the dataset
      if (is.null(values))
        values <- sort(unique(data[ ,attribute]))
      
      # Prepare an empty result
      count <- matrix(0, nrow = length(values) + 1, ncol = 1)
      rownames(count) <- c(as.character(values), 'total')
      colnames(count) <- 'count'
      
      # Count the data by type
      s_filt <- SamplesFilter$new()
      for(t in values) {
        selection <- s_filt$filter_by_attribute(data, attribute, t)
        count[as.character(t),'count'] <- nrow(selection)
      }
      
      count['total','count'] <- sum(count[as.character(values),'count'])
      
      return(count)
      
    },
    
    
    
    
    
    #' Returns the total number of samples on the given dataset with the 
    #' specified feature value. It is actually an element of count_by_attribute 
    #' 
    #' @param adapted_dataset  A dataset in the common format on which the samples 
    #'                         are counted
    #' @param attribute        name of the attribute to search for
    #' @param value            value of the feature to search for
    count_by_attribute_value = function(data, attribute, value){
      
      return(self$count_by_attribute(data, attribute, value)[value, ])
    
    },
    
    
    
    
    #' Creates a dataframe that specifies, for each patient, the number of samples
    #' for each value of the given attribute. If only some values are needed,
    #' they can be specified as character vector.
    #' 
    #' @param adapted_dataset  the dataset where the samples are counted
    #' @param attribute        the attribute for which values must be counted
    #' @param values           the values searched, otherwise NULL
    #' 
    #' @return A dataframe with the count of samples for each attribute value
    count_patient_samples = function(adapted_dataset, attribute, values = NULL){
      
      # Check input
      stopifnot(check_type(adapted_dataset,'data.frame',c(1,1)))
      stopifnot(PATIENT_LABEL %in% colnames(adapted_dataset))
      stopifnot(check_type(attribute, 'character',1,1), attribute %in% colnames(adapted_dataset))
      
      # Get samples by type for each patient
      sm         <- SamplesManager$new()
      samples    <- sm$get_patient_samples_by_type(adapted_dataset, attribute, values)
      
      # Get distinct attribute values
      if (is.null(values))
        types <- sort(unique(adapted_dataset[,attribute]))
      else
        types <- values
      
      # Overall count result
      counts_tot <- data.frame(names(samples))
      n_patients <- length(names(samples))
      # Keep the count for each attribute value
      
      for (t in types){
        
        count  <- vector(mode = "integer", length = n_patients)
        # Get the count for each patient
        for (i in seq(n_patients)){ 
          p <- names(samples)[i]
          count[i]  <- length(samples[[p]][[t]])
        }
        # Add count to dataframe
        counts_tot <- cbind(counts_tot, count)
      }    
      
      # Set colnames
      colnames(counts_tot) <- c(PATIENT_LABEL, as.character(types))
      if (length(types) > 1)
        counts_tot  <- cbind(counts_tot, total = rowSums(counts_tot[,-1]))
      else
        counts_tot  <- cbind(counts_tot, total = counts_tot[,types])
      return(counts_tot)
    },
    
    
    
    
    
    
    #' Returns a matrix with statistical measures for each gene of the given
    #' expression matrix. The computed measures (min, max, median, mean, 1st
    #' quartile, 3rd quartile, standard deviation) aggregate the measures on
    #' all samples, for each gene.
    #' 
    #' @param expression_matrix The matrix on which to compute the statistics
    #' @return A matrix with genes on the rows and statistical measures
    #' on the columns. A summary of this matrix is printed too.
    #' 
    #' References
    #' https://stackoverflow.com/questions/10945703/calculate-row-means-on-subset-of-columns
    #' http://www.r-tutor.com/elementary-statistics/numerical-measures/quartile
    get_genes_statistics = function(expression_matrix){
      
      # Check input
      if (!check_type(expression_matrix, "matrix",c(1,1)) | !check_type(expression_matrix,'numeric')){
        stop("Provide a matrix as input for gene statistics.")
      }
      
      statistics <- data.frame(GeneID = rownames(expression_matrix)) %>% 
        cbind(Min       = apply(X = data.frame(expression_matrix), 1, min    )) %>% 
        cbind('1st Qu.' = apply(X = data.frame(expression_matrix), 1, quantile, probs = 0.25)) %>%
        cbind(Mean      = apply(X = data.frame(expression_matrix), 1, mean   )) %>% 
        cbind(Median    = apply(X = data.frame(expression_matrix), 1, median )) %>%
        cbind('3rd Qu.' = apply(X = data.frame(expression_matrix), 1, quantile, probs = 0.75)) %>%
        cbind('80% Qu.' = apply(X = data.frame(expression_matrix), 1, quantile, probs = 0.80)) %>%
        cbind(Max       = apply(X = data.frame(expression_matrix), 1, max    )) %>%
        cbind(Sd        = apply(X = data.frame(expression_matrix), 1, sd     )) %>%
        format(scientific = FALSE) 
      statistics[,-1] <- as.numeric(as.matrix(statistics[,-1]))
      
      return(statistics)
      
    },
    
    
    #' Summarize the result of the NTP classification by counting, for each class:
    #' 1. the samples confidentially assigned to it (BH.FDR low) and 
    #' 2. the samples not confidentially assigned to it (BH.FDR high) and
    #' 3. the overall samples assigned to it (any BH.FDR)
    #' 
    #' @param classification    the result of an NTP classification
    #' @result A dataframe with classes on the rows and above categories on columns
    ntp_summary_by_class = function(ntp_pred){
      
      if (is.null(ntp_pred) | any(dim(ntp_pred) == 0)) {
        stop("Provide a non-null and non-empty NTP classification result.")
        return(data.frame())
      }
      
      # Separate the samples confidently and non-confidently classified
      best_fdr <- as.numeric(ntp_pred[,BEST_FDR_LABEL])
      conf     <- ntp_pred[best_fdr <  BH_FDR_THRESHOLD, ]
      non_conf <- ntp_pred[best_fdr >= BH_FDR_THRESHOLD, ]                           
    
      # Count the samples for each class, considering confidently classified,
      # non confidently classified and overall classified samples
      summary <- data.frame(
        self$count_by_attribute(conf,     CLASS_LABEL, levels(F_CRIS_CLASSES)),
        self$count_by_attribute(non_conf, CLASS_LABEL, levels(F_CRIS_CLASSES)),
        self$count_by_attribute(ntp_pred, CLASS_LABEL, levels(F_CRIS_CLASSES))
      )
      
      # Percentage of samples accurately classified
      summary <- cbind(summary, summary[,1]/summary[3])
      
      # Set rows and column names
      colnames(summary) <- c(paste("BH.DFR <" , BH_FDR_THRESHOLD, sep = " "),
                             paste("BH.DFR >=", BH_FDR_THRESHOLD, sep = " "),
                             "Overall classified",
                             "Percentage confidently classified")
    
      rownames(summary) <- c(levels(F_CRIS_CLASSES),"Classified")
      
      # Round numbers at 4 significant digits
      summary <- round(summary,4)
      
      return(summary)
      
    },
    
    
    #' Summarize the result of the TSP classification by saying:
    #' 1. the samples confidentially assigned to each class (not NA)
    #' 2. the samples that remained unclassified
    #' 3. the overall samples analyzed (classified + unclassified)
    #' 
    #' @param classification  The result of a TSP classification
    #' @return A dataframe with the above counting for each class
    tsp_summary_by_class = function(tsp_pred){
      
      # Check the input is non-null and non-empty
      if (is.null(tsp_pred) | any(dim(tsp_pred) == 0)) {
        stop("Provide a non-null and non-empty TSP prediction.")
        return(data.frame())
      }
      
      # Extract the classified samples
      classified   <- 
        Filter$new()$filter_by_attribute(tsp_pred, CLASS_LABEL, F_CRIS_CLASSES)
      
      # Count the classified samples for each class
      classified_count   <- self$count_by_attribute(classified, CLASS_LABEL)
     
      # Total of classified, unclassified samples
      tot_unclass  <- self$count_by_attribute_value(tsp_pred, CLASS_LABEL, NA)
      tot_classif  <- classified_count['total', 'count']
      
      # Overall number of samples
      total   <- sum(tot_classif + tot_unclass)
      
      # Put counts in a single dataframe
      summary <- data.frame(c(classified_count[,'count'], tot_unclass, total))
      
      rownames(summary) <-
        c(levels(F_CRIS_CLASSES), "Classified", "Unclassified", "Total")
      colnames(summary) <- c("Count of samples")
      
      return(summary)
      
    }
    
    
    
  )
  
)
