# Description -------------------------------------------------------------

# Manager class for samples

# Dependencies ------------------------------------------------------------

library(R6)

# Class definition --------------------------------------------------------

SamplesManager  <- R6Class('SamplesManager',            
   lock_objects = FALSE,
   lock_class   = TRUE,
   
   public = list(
     
     #' Function to read the region data of a single sample, whose file path is given.
     #' Extension check is preformed too, together with column names for the resulting
     #' dataframe
     #' 
     #' @param path            The path of sample file
     #' @param file_extension  The extension of the file
     #' @param region_colnames The names of the region data fields
     get_single_sample = function(path, file_extension, region_colnames){
       
       if (!is_path_valid(path, 'INPUT', extension = file_extension)){
         stop("Path not valid")
         return(data.frame())
       }
       
       sample <- read.table(path, sep = "\t", header = FALSE)
       colnames(sample) <- region_colnames
       return(sample)
     },
     
     
     
     
     
     #' Extract the type of samples for each patient
     #' 
     #' @param adapted_dataset  the dataset where the samples are retrieved
     #' 
     #' @return List with, for each patient, a list of aliquots, subdivided
     #' by attribute type
     get_patient_samples_by_type = function(adapted_dataset, attribute, values = NULL){
       
       # Check input
       stopifnot(check_type(adapted_dataset,'data.frame',c(1,1)))
       stopifnot(check_type(attribute, 'character',1,1))
       stopifnot(attribute %in% colnames(adapted_dataset))
       
       # Identify type of samples in the dataset
       if (is.null(values))
         types <- sort(unique(adapted_dataset[ ,attribute]))
       else
         types <- values
       
       # Select the samples relative to the given patient
       patients <- unique(adapted_dataset[,PATIENT_LABEL])
       
       # Distinguish the samples by type
       result <- list()
       for (p in patients) {
         # Find patient samples
         samples  <- adapted_dataset[adapted_dataset[,PATIENT_LABEL] == p, ]
         
         # Divide aliquots by type
         patient_data <- list()
         for (t in types){
           if (t == TRUE)
            patient_data[['TRUE']] <- samples[samples[,attribute] == t, ALIQUOT_LABEL]
           else if (t == FALSE)
            patient_data[['FALSE']] <- samples[samples[,attribute] == t, ALIQUOT_LABEL]
           else
            patient_data[[t]] <- samples[samples[,attribute] == t, ALIQUOT_LABEL]
         }
         
         # Add aliquots to patient list
         result[[p]] <- patient_data
       }
       
       return(result)
     },
     
     
     
     
     
     
     #' Extract all the samples for each patient
     #' 
     #' @param adapted_dataset  the dataset where the samples are retrieved
     #' 
     #' @return List with, for each patient, a list of all its samples
     get_all_patient_samples = function(adapted_dataset){
       
       stopifnot(check_type(adapted_dataset,'data.frame',c(1,1)))
       
       # Select the samples relative to the given patient
       patients <- unique(adapted_dataset[,PATIENT_LABEL])
       
       # Distinguish the samples by type
       result <- list()
       for (p in patients) {
         # Find patient samples
         samples     <- adapted_dataset[adapted_dataset[,PATIENT_LABEL] == p, ]
         result[[p]] <- samples[, ALIQUOT_LABEL]
       }
       
       return(result)
     }
     
     
     
   )
)
