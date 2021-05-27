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
.settings_file    <- path_loader$get_classifier_file_path('sl', .FS_TYPE, .TUNE, path_type = 'train_settings')


# Classifier settings -----------------------------------------------------

# Lists for saving models and train settings used for the computation
models <- list()
train_settings  <- list()

# Load data with settings specified in the config.r file
sldata <- load_prepared_tcga_data(confident = .CONFIDENT_ONLY, 
                                  uniformed = .UNIFORMED, 
                                  fs_type   = .FS_TYPE, 
                                  type      = 'sl',
                                  load_training = TRUE)

# Hold the result of the training pipeline
train_res <- list()

for (m in methods){
    
  print_debug(m)
  
  train_res[[m]] <- sl_pipeline_train(
      sldata  = sldata,
      method  = m,
      sl_formula = .SL_FORMULA,
      seed    = .SEED,
      tune    = .TUNE,
      fit_control = fit_control,
      tune_grid = tune_grid[[m]]
    )
  
  # Save the results
  models[[m]] <- train_res[[m]]$model
  
  
  # Save the settings of the training
  train_settings[[m]] <- list(
    fit_control = train_res[[m]]$fit_control,
    tune_grid   = train_res[[m]]$tune_grid,
    cv_res      = train_res[[m]]$cv_res,
    best_hp     = train_res[[m]]$best_hp
  )
  
  # Save a time flag to specify last time both models and settings have been updated
  models[['last_update']] <- Sys.time()
  train_settings[['last_update']] <- models[['last_update']]
  
  # If requested, save the models and the settings
  if (.SAVE){
    saveRDS(train_settings, .settings_file)
    saveRDS(models, .model_file)
  }
}