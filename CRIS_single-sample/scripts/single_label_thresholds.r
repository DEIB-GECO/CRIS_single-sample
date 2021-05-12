# Description -------------------------------------------------------------

# Training on TCGA for SL classification

# Dependencies ------------------------------------------------------------

source(here('src','load_libraries.r'))
source(here('src','utils','source_utils.r'))
source(here('src','loader_writer','load_data.r'))
source(here('src','data_management','source_data_management.r'))
source(here('src','classifiers','source_classifiers.r'))
source(here('src','pipelines','source_pipelines.r'))

# Configuration constants -------------------------------------------------

# Path for saving a file with all single-label classifier models
.model_file       <- path_loader$get_classifier_file_path('sl', .FS_TYPE, .TUNE, path_type = 'models')

# Path for saving a file with settings used in the training of the models
.thresholds_file  <- path_loader$get_classifier_file_path('sl', .FS_TYPE, .TUNE, path_type = 'thresholds')

# Flag to decide if saving the results on file system or not
.SAVE <- TRUE

# Classifier settings -----------------------------------------------------

# Lists for saving models and train settings used for the computation
models <- load_file(path = .model_file)
thresholds  <- list()

# Load data with settings specified in the config.r file
sldata <- load_prepared_tcga_data(confident = .CONFIDENT_ONLY, 
                                  uniformed = .UNIFORMED, 
                                  fs_type   = .FS_TYPE, 
                                  type      = 'sl',
                                  load_training = TRUE)

mldata <- load_prepared_tcga_data(confident = .CONFIDENT_ONLY, 
                                  uniformed = .UNIFORMED, 
                                  fs_type   = .FS_TYPE, 
                                  type      = 'ml',
                                  load_training = TRUE)

# Hold the result of the training pipeline
threshold_res <- list()

# Compute the thresholds for the methods appearing in the models file
for (m in intersect(methods, names(models))){
    
  print_debug(m)
  
  # File for saving the ROC curves plot
  method <- paste('sl',m, 'train', sep = '_')
  roc_file <- path_loader$get_classifier_file_path(method, .FS_TYPE, .TUNE, path_type = 'roc')
  
  # Compute the thresholds
  threshold_res[[m]] <- sl_pipeline_thresholds(
    sldata = sldata,
    method = m,
    seed     = .SEED,
    model    = models[[m]],
    png_path = roc_file,
    mldata   = mldata)
  
  # Save a time flag to specify last time thresholds have been updated
  threshold_res[['last_update']] <- Sys.time()
  
  # If requested, save the models and the settings
  if (.SAVE){
    saveRDS(threshold_res, .thresholds_file)
  }
}