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

# Flag to decide if saving the results on file system or not
.SAVE <- TRUE

# Classifier settings -----------------------------------------------------

# Lists for saving models that have been trained
models <- load_file(.model_file)

# Load data with settings specified in the config.r file
sldata <- load_prepared_tcga_data(confident = .CONFIDENT_ONLY, 
                                  uniformed = .UNIFORMED, 
                                  fs_type   = .FS_TYPE, 
                                  type      = 'sl',
                                  load_training = TRUE)

# Hold the result of the testing_pipeline
testing_res <- list()

for (m in intersect(methods, names(models))){
    
  print_debug(m)
  method <- paste(.DATA, 'sl', sep = '_')
  testing_file <- path_loader$get_classifier_file_path(method, .FS_TYPE, .TUNE, path_type = 'testing', testing_folder = .DATA)
  testing_res[[m]] <- sl_pipeline_test(
      sldata = sldata,
      method = m,
      seed = .SEED,
      model = models[[m]]
    )
  
  # Save the time flag that identifies the model version used
  testing_res[['models_last_update']] <-  models[['last_update']]
  
  # Save the time of last testing
  testing_res[['last_update']] <-  Sys.time()
  
  # If requested, save the models and the settings
  if (.SAVE){
    saveRDS(object = testing_res, testing_file)
  }
}