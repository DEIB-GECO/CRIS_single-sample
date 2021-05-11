library(R6)
library(tidyverse)

SummaryMakerPDX  <- R6Class('SummaryMakerPDX',
                             
     inherit = SummaryMakerTCGA,
     lock_objects = FALSE,
     lock_class   = TRUE,
     
     ####### Public fields #######
     public = list(
       
             
             
       #' Create a dataframe that shows the number of samples for each type of tumor
       #' for CRC samples (colorectal cancer)
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
                                                  "CRC")
         
         tumor_types   <- levels(factor(adapted_dataset[,TUMOR_TYPE_LABEL]))
         n_tumor_types <- length(tumor_types)
         
         # Count the tumor types for both colon and rectum
         colon_count  <- vector("integer", n_tumor_types)
         
         for (i in seq(n_tumor_types)) {
           colon_count[i]  <- self$count_by_attribute_value(colon_samples, 
                                                              TUMOR_TYPE_LABEL, 
                                                              tumor_types[i])
         }
         
         # Create and return a result dataframe
         tumor_type_counts <- data.frame(CRC  = colon_count)
         rownames(tumor_type_counts) <- tumor_types
         
         return(tumor_type_counts)
       },
       
       
       batch_distribution = function(pdx_list){
         
         # Create dataframe with samples of each batch
         distr <- data.frame()
         
         for (b in seq(length(pdx_list))) {
           
           # Get expression matrix
           expr <- pdx_list[[b]]@assayData$exprs
           
           for (id in colnames(expr)) {
             
             # Total reads of sample
             reads <- sum(expr[, id])
             
             # Insert in (not) empty dataframe
             if (any(dim(distr)) == 0) {
               distr <- data.frame(id, b, reads)
               colnames(distr) <- c(ALIQUOT_LABEL,
                                    'batch',
                                    'reads')
             }else {
               distr <- rbind(distr, c(id, b, reads))
             }
           }
           
         }
         
         # Extract patient ids
         patient_ids <- extract_pdx_id(distr[,ALIQUOT_LABEL], PATIENT_LABEL)
         distr       <- cbind(patient_ids, distr)
         colnames(distr)[1] <- PATIENT_LABEL
         
         # Sort rows
         distr <- distr %>% arrange(patient_id, aliquot_id, batch, reads)
         
         # Return the distribution
         return(distr) 
       },
       
       
       
       technical_replicas = function(batch_distr){
         
         n_batches <- as.numeric(max(batch_distr[,'batch']))
         
         tech_replicas <- as.data.frame(batch_distr) %>% 
          group_by(aliquot_id) %>%
          count(name = 'technical') %>%
          as.data.frame()
                                      
         tech_plot <- 
           ggplot(tech_replicas, aes(x = tech_replicas[,ALIQUOT_LABEL],
                                     y = tech_replicas[,'technical'])) + 
           geom_col() + 
           labs(title = 'Technical replicas for each sample') + 
           xlab('aliquot id') + 
           ylab('count')
         
         counts_l <- list()
         for (b in seq(n_batches)){
           b_tech      <- batch_distr %>% group_by(aliquot_id) %>%
                          mutate(technical = n()) %>%
                          filter(batch == b) %>%
                          as.data.frame()
           
           counts_l[[b]] <- self$count_by_attribute(b_tech, 'technical')
         }
         
         
         return(list(
           data   = tech_replicas,
           plot   = tech_plot,
           global_count = self$count_by_attribute(tech_replicas, 'technical'),
           batch_count  = counts_l
         ))
            
       },
       
       
       
       select_technical_replicas = function(batch_distr, by = 'batch'){
         
         if(!any(by == c('batch', 'reads'))) {
           stop('Can filter only by batch or by max number of reads.')
         }
         
         # Apply filter on everly aliquot id
         filt_distr <- as.data.frame(batch_distr) %>% group_by(aliquot_id)
         
         # If necessary, select the replicas with highest number of reads
         if (by == 'reads')
           filt_distr <- filt_distr %>% filter(reads == max(reads))
        
         # Select, in remaining replicas, the one in highest batch
         filt_distr <- filt_distr %>% filter(batch == max(batch)) %>% 
            select(patient_id, aliquot_id, batch, reads) %>%
            as.data.frame()
         
         return(filt_distr)
       
       },
       
       
       
       biological_replicas = function(batch_distr) {
         
         n_batches <- as.integer(max(batch_distr$batch))
         
         # Distribution of biological replicas by batch
         batch_bio_distr  <- batch_distr %>% 
                             group_by(patient_id, batch) %>%
                             count(name = 'biological')

         bio_list   <- list()
         bio_plot   <- list()
         count_list <- list()
         for (b in seq(n_batches)) {
           
           # Select data of current batch
           data <- batch_bio_distr %>% filter(batch == b) %>%
                   as.data.frame()

           # Plot count of replicas for current batch
           bio_plot[[b]] <-
             ggplot(data, aes(x = patient_id, y = biological)) +
             geom_col() +
             labs(title = paste('Biological replicas by patient, batch', b)) +
             xlab('patient id') +
             ylab('count')
          
          bio_list[[b]]   <- data
          count_list[[b]] <- self$count_by_attribute(data, 'biological')
           
         }

         # Global distribution of biological replicas
         glob_bio_distr <- batch_distr %>%
                           group_by(patient_id) %>%
                           count(name = 'biological') %>%
                           as.data.frame()

         glob_bio_count <- self$count_by_attribute(glob_bio_distr, 'biological')
         
         glob_bio_plot <-
           ggplot(glob_bio_distr, aes(x = glob_bio_distr[,PATIENT_LABEL],
                                      y = glob_bio_distr[,'biological'])) +
           geom_col() +
           labs(title = 'Biological replicas by patient, globally') +
           xlab('patient id') +
           ylab('count')

         # Return global/by-batch results
         return(list(
           byBatch = list(
             data  = bio_list,
             counts = count_list,
             plot  = bio_plot
           ),
           global  = list(
             data  = glob_bio_distr,
             counts = glob_bio_count,
             plot  = glob_bio_plot
           ))
         )
         
       }
      
     )
)
