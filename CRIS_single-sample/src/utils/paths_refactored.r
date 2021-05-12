# Description -------------------------------------------------------------

# Script that handles the paths

# Libraries ---------------------------------------------------------------

source(here('src','utils','misc.r'))

# Class definition --------------------------------------------------------

PathLoader <- R6Class("PathLoader",
    lock_objects = TRUE,
    lock_class = TRUE,
    
    private = list(
      .path_db = NULL,
      .allowed_extensions = c(".gct",".csv", ".rda", ".xlsx", ".rds")
    ),
    
    public = list(
      
      #' Get the path of the required source file given its label.
      #' 
      #' @param  file_label    Label of the file (string)
      #' @return Returns the path corresponding of the requested file, if it exists,
      #' otherwise returns NULL.
      get_path = function(file_label){
        
        if (!check_type(file_label, 'character',1,1))
          return(NULL)
        
        selectedPath <- fs::path_rel(private$.path_db[file_label, 'path'])
	
      	if (is.null(selectedPath) | is.na(selectedPath)) {
      		stop(paste("PathLoader: No matching file label has been found for", file_label))
      	  return(NULL)
      	}
      
      	return(here(selectedPath))
        
      },
      
      #' Initialize the PathLoader reading the file with the list of paths
      #' 
      #' @param  files_db_path    Path of the file containing the db of all loadable sources.
      #' @return Creates the instance of the PathLoader
      initialize = function(files_db_path){
        
        # Check the file paths dataset path
        if (!check_type(files_db_path, 'character',1,1))
          stop("PathLoader: path of the files paths database is not valid.")
        
        if (!fs::file_exists(files_db_path))
          stop("PathLoader: files paths database not found")
        
        if (!grepl(x = files_db_path, pattern = '.xlsx',fixed = TRUE))
          stop("PathLoader: files paths database must be an xlsx file.")
        
        # Read the file
        private$.path_db <- openxlsx::read.xlsx(files_db_path)
        
        if (!is.data.frame(private$.path_db))
          stop("PathLoader: could not initialize the files paths database.")
        
        # Check the file structure
        
        if (ncol(private$.path_db) != 2)
          stop("PathLoader: file paths database must contain 2 character columns (file label, file path).")
        
        if (nrow(private$.path_db) < 1)
          warning("PathLoader: file paths database is empty.")
        
        if (mode(private$.path_db[,1]) != 'character' |mode(private$.path_db[,2]) != 'character')
          stop("PathLoader: file label and file path column must be characters")
        
        # Set colnames of the file paths database
        colnames(private$.path_db) <- c('label','path')
        
        # Remove rows with NA path or NA label and index with label
        private$.path_db <- private$.path_db %>% 
          dplyr::filter(!is.na(label) & !is.na(path)) %>%
          as.data.frame() %>% 
          column_to_rownames('label')
        
      },
      
      #' Checks that the source path has been built correctly
      #' 
      #' @param path path to check
      #' @param extension the type of extensions (.csv, .gct, .rda, .rds, .xlsx)
      #' @return TRUE or FALSE
      has_extension = function(path, extension = ""){
        
        # Check input parameter validity
        if (!check_type(path,'character',1,1) | !check_type(extension,'character',1,1))
          stop("PathLoader: invalid input.")
      
        # Check the extension of the file
        extension_check <- FALSE
        if (any(tolower(extension) == private$.allowed_extensions))
          extension_check <- grepl(extension, path, fixed = TRUE)
        else
          warning(paste("PathLoader: extension", extension, "not found in the path."))
          
        return(extension_check)
      },
      
      get_classifier_file_path = function(type, fs_type, tuned, path_type, testing_folder = NULL){
        
        # Check input parameter validity
        if (!check_type(type,'character',1,1) | !check_type(fs_type,'character',1) |
            !check_type(tuned,'logical',1,1)  | !check_type(path_type,'character',1,1))
          stop("PathLoader: invalid input format.")
      
        if (any_uncomplete(c(type, fs_type, tuned, path_type)))
          stop("PathLoader: some input is NULL, NA, NaN or Infinite")
        
        # Check range of values
        
        if (!path_type %in% c('models','train_settings', 'thresholds', 'roc', 'testing'))
          stop('PathLoader: path_type can be only models or train_settings or thresholds or roc or testing.')
          
        if (any(!fs_type %in% c('bio_driven','ntp_only','lasso','none')))
          stop('PathLoader: type can include only bio_driven,ntp_only,lasso,none')
        
        if (any(!tuned %in% c(TRUE, FALSE)))
          stop('PathLoader: tuned can be only true or false')
        
        # Set the subfolder (to be placed in feature set folder)
        if (path_type == 'testing'){
          if (!check_type(testing_folder,'character',1,1) | any_uncomplete(testing_folder))
            stop("PathLoader: testing folder must be specified")
          
          if (any(!testing_folder %in% c('tcga','pdx','comparison')))
            stop('PathLoader: testing_folder can be only tcga or pdx or comparison.')
          
          subfolder <- paste(path_type, testing_folder, sep = '/')
        }else{
          subfolder <- path_type
        }
        
        # Create sub-directory for feature selection and ROC curves and set extension
        base_folder        <- self$get_path('OUT_FOLDER_MODELS')
        feature_set_folder <- paste(sort(fs_type), collapse = '_')
        folder <- dir_create(paste(base_folder, feature_set_folder, subfolder, sep = '/'))
        
        # Build the file name
        if (tuned)
          tuned_label <- 'custom_tuned'
        else
          tuned_label <- 'default'
        
        
        if (path_type == 'roc'){
          extension <- '.png'
        }else{
          extension <- '.rds'
        }
      
        filename <- paste(type, path_type, tuned_label, sep = '_') %>%
                    paste(extension, sep = '')
        
        # return the complete path
        return(paste(folder, filename, sep = '/'))

      }
      
    ),
    
    active = list(
      path_db = function(v){
        if (missing(v))
          return(private$.path_db)
        else
          stop('PathLoader: path_db is read-only.')
      },
      
      allowed_extensions = function(v){
        if (missing(v))
          return(private$.allowed_extensions)
        else
          stop('PathLoader: allowed_extensions is read-only.')
      }
    )
    

)
