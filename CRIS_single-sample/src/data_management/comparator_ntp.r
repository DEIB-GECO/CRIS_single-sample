# Description -------------------------------------------------------------

# Class to compare datasets/results of classification

# Dependencies ------------------------------------------------------------

library(R6)

# Class definition --------------------------------------------------------

NTPComparator  <- R6Class('NTPComparator',              
   
   inherit      = Comparator,
   lock_objects = FALSE,
   lock_class   = TRUE,
   
   public = list(
     
     
     
     #' Create a dataframe comparing corresponding samples of two 
     #' NTP classification results.The first dataframe to be given must
     #' have shorter (or comparable) ids because the correspondence if found 
     #' checking if the ids inside the second result dataframe contain the 
     #' ids inside the first result dataframe. Simple equality of ids cannot be
     #' done as it is not ensured that the ids have the same length.
     #'   
     #' 
     #' @param result_1   A result of an NTP classification (the one with shortest ids)
     #' @param result_2   A result of an NTP classification (the one with longest or 
     #'                   comparable ids)
     #' @return Comparison of two results. A dataframe with ids, labels, distances and 
     #'         FDR of the assigned class for both classifications
     compare_results = function(result_1, result_2){
       
       # Check input
       super$compare_results(result_1, result_2)
       
       # Check ids length
       id_1 <- as.character(result_1[1,1])
       id_2 <- as.character(result_2[1,1])
       
       if (compare_id_length(id_1,id_2) == 1) {
         warning("Ids' lengths are not compatible. Try to swap the datasets.")
         return(NULL)
       }
       
       # Build result
       result_comparison <- data.frame()
       n_samples_df1 <- nrow(result_1)
       
       # Search correspondences of ids of the first datasets in the second
       # and create a row in the result
       for (i in seq(n_samples_df1)) {
         
         # Extract ID, class, FDR and distance for sample of first dataframe
         id_1    <- as.character(result_1[i,1])
         class_1 <- as.character(result_1[i,CLASS_LABEL])
         dist_1  <- as.numeric(result_1[i,BEST_DISTANCE_LABEL])
         FDR_1   <- as.numeric(result_1[i,BEST_FDR_LABEL])
         
         # IDs in the second dataset that correspond to id_1
         corr_samples <- result_2[startsWith(result_2[,1], prefix = id_1), ]
         n_corr_samples <- nrow(corr_samples)
         
         if (!is.null(corr_samples) & n_corr_samples > 0)
           # For each corresponding sample, extract ID, class, FDR and distance
           for (k in seq(n_corr_samples)) {
             id_2    <- as.character(corr_samples[k,1])
             class_2 <- as.character(corr_samples[k,CLASS_LABEL])
             dist_2  <- as.numeric(corr_samples[k,BEST_DISTANCE_LABEL])
             FDR_2   <- as.numeric(corr_samples[k,BEST_FDR_LABEL])
             
             # Add a row with the extracted fields to the result
             comparison_row <- c(id_1,  id_2,  class_1, class_2, 
                                 FDR_1, FDR_2, dist_1,  dist_2 )
             result_comparison <- rbind(result_comparison, comparison_row)
           }
         
       }
       
       # If something went wrong, the result is empty
       if (ncol(result_comparison) == 0) 
         warning("Empty comparison dataset.")
       else{
         
         colnames(result_comparison) <- COMPARISON_NTP_LABELS
         
         # Factorize class labels
         result_comparison[, COMPARISON_NTP_LABELS[3]] <-
           factor(result_comparison[, COMPARISON_NTP_LABELS[3]],
                  levels = sort(unique(result_comparison[, COMPARISON_NTP_LABELS[3]])))
         
         result_comparison[, COMPARISON_NTP_LABELS[4]] <-
           factor(result_comparison[, COMPARISON_NTP_LABELS[4]],
                  levels = sort(unique(result_comparison[, COMPARISON_NTP_LABELS[4]])))
       }
       return(result_comparison);
     },
     
     
     
     
     
     #' Extracts the samples that differ in two classifications because of the given condition
     #' 
     #' @param compared_ntp_result   The comparison of the result of two ntp classifications
     #' @param condition             The condition to be considered (CLASS_LABEL, 
     #'                              BEST_DISTANCE_LABEL, BEST_FDR_LABEL, DIFF_CONFIDENCE, 
     #'                              UNCLASSIFIED)
     #' @param sensibility           Minimum difference of distance/fdr (optional). Used
     #'                              only with conditions BEST_DISTANCE_LABEL, BEST_FDR_LABEL
     check_comparison = function(compared_ntp_result, condition, sensibility = 0.0001){
       
       if (is.null(compared_ntp_result)) {
         stop("Given comparing dataframe cannot be NULL.")
         return(NULL)
       }
       
       if (!condition %in% c(CLASS_LABEL,BEST_DISTANCE_LABEL, BEST_FDR_LABEL, 
                             DIFF_CONFIDENCE, UNCLASSIFIED)) {
         stop(paste("The attribute on which to check differences must be either:", 
                           CLASS_LABEL,
                           BEST_DISTANCE_LABEL,
                           BEST_FDR_LABEL,
                           DIFF_CONFIDENCE,
                           UNCLASSIFIED, sep = "\n"))
         return(NULL)
       }
       # Get and cast the required fields
       cls1  <- as.character(compared_ntp_result[,COMPARISON_NTP_LABELS[3]])   
       cls2  <- as.character(compared_ntp_result[,COMPARISON_NTP_LABELS[4]])
       dist1 <- as.numeric(compared_ntp_result[,COMPARISON_NTP_LABELS[7]])    
       dist2 <- as.numeric(compared_ntp_result[,COMPARISON_NTP_LABELS[8]])
       FDR1  <- as.numeric(compared_ntp_result[,COMPARISON_NTP_LABELS[5]])    
       FDR2  <- as.numeric(compared_ntp_result[,COMPARISON_NTP_LABELS[6]]) 
       
       # Select the filter to check according to the given attribute
       # Avoid comparison with "==" that may lead to numerical errors
       
       # Misclassification
       if (condition == CLASS_LABEL)                
         filter <- cls1 != cls2
       
       # distance diff >= sensibility
       else if (condition == BEST_DISTANCE_LABEL)   
         filter <- !(abs(dist1 - dist2) < sensibility)  
       
       # fdr diff >= sensibility 
       else if (condition == BEST_FDR_LABEL)       
         filter <- !(abs(FDR1 - FDR2) < sensibility)     
       
       # not confidently classified in one case 
       else if (condition == DIFF_CONFIDENCE)      
         
         filter <- ((!FDR1  < BH_FDR_THRESHOLD) &  FDR2 <  BH_FDR_THRESHOLD  ) |
         (FDR1   < BH_FDR_THRESHOLD  &  !(FDR2 < BH_FDR_THRESHOLD) )
       
       # unclassified in one case  
       else if (condition == UNCLASSIFIED)         
         filter <- is.na(cls1) || is.na(cls2)
       
       
       return(compared_ntp_result[filter, ])
     }
     
   )

)
                         
