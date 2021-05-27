# Description -------------------------------------------------------------

# Thresholds for assignment of binary target of multi-label problem transformation
# classifiers

# Dependencies ------------------------------------------------------------
library(here)
source(here('src','load_libraries.r'))
source(here('src','utils','source_utils.r'))
source(here('src','loader_writer','load_data.r'))
source(here('src','data_management','source_data_management.r'))
source(here('src','classifiers','source_classifiers.r'))
source(here('src','pipelines','source_pipelines.r'))

# Configuration constants -------------------------------------------------

# Path for saving a file with all multi-label classifier models
.model_file <- path_loader$get_classifier_file_path('ml', .FS_TYPE, .TUNE, path_type = 'models')

# Path for saving a file with settings used in the training of the models
.thresholds_file  <- path_loader$get_classifier_file_path('ml', .FS_TYPE, .TUNE, path_type = 'thresholds')


# Classifier settings -----------------------------------------------------

# Lists for saving models and train settings used for the computation
models <- load_file(path = .model_file)
thresholds  <- list()

# Load data with settings specified in the config.r file
mldata <- load_prepared_tcga_data(confident = .CONFIDENT_ONLY, 
                                  uniformed = .UNIFORMED, 
                                  fs_type   = .FS_TYPE, 
                                  type      = 'ml',
                                  load_training = TRUE)
# Hold the result of the training pipeline
threshold_res <- list()

# Compute the thresholds for the classifiers appearing in the models file
for (m in intersect(names(settings), names(models))){
    
  print_debug(m)
  
  # File for saving the ROC curves plot
  method <- paste('ml',m, 'train', sep = '_')
  roc_file <- path_loader$get_classifier_file_path(method, .FS_TYPE, .TUNE, path_type = 'roc')
  
  # Compute the thresholds
  res <- ml_pipeline_thresholds(
    mldata = mldata,
    seed = .SEED,
    cv_set = CVSettings$new(),
    alg_set = settings[[m]],
    model = models[[m]],
    png_path = roc_file
    )
  
  threshold_res[[m]] <- res$thresholds
  
  # Save a time flag to specify last time thresholds have been updated
  threshold_res[['last_update']] <- Sys.time()
  threshold_res[['models_last_update']] <- models[['last_update']]
  
  
  # If requested, save the models and the settings
  if (.SAVE){
    saveRDS(threshold_res, .thresholds_file)
  }
}