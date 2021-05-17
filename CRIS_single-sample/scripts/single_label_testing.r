# Description -------------------------------------------------------------

# Training on TCGA for SL classification

# Dependencies ------------------------------------------------------------
library(here)
source(here('src','load_libraries.r'))
source(here('src','utils','source_utils.r'))
source(here('src','loader_writer','load_data.r'))
source(here('src','data_management','source_data_management.r'))
source(here('src','classifiers','source_classifiers.r'))
source(here('src','pipelines','source_pipelines.r'))

# Configuration constants -------------------------------------------------

if(any(.PUBLISHED_MODELS)){
  .model_file  <- path_loader$get_path('NTP_ONLY_SL_MODELS')
}else{
  .model_file <- path_loader$get_classifier_file_path('sl', .FS_TYPE, .TUNE, path_type = 'models')
}

method <- paste(.DATA, 'sl', sep = '_')
testing_file <- path_loader$get_classifier_file_path(method, .FS_TYPE, .TUNE, path_type = 'testing', testing_folder = .DATA)

# Classifier settings -----------------------------------------------------

# Lists for saving models that have been trained
models <- load_file(.model_file)

# Load data with settings specified in the config.r file

if (.DATA == 'tcga'){
  sldata <- load_prepared_tcga_data(confident = .CONFIDENT_ONLY, 
                                    uniformed = .UNIFORMED, 
                                    fs_type   = .FS_TYPE, 
                                    type      = 'sl',
                                    load_training = TRUE)
}else if (.DATA == 'pdx'){
  
  sldata <- load_prepared_pdx_data(confident = .CONFIDENT_ONLY, 
                                    uniformed = .UNIFORMED, 
                                    fs_type   = .FS_TYPE, 
                                    type      = 'sl')

}else {
  stop('.DATA flag must be either tcga or pdx')
}

# Hold the result of the testing_pipeline
testing_res <- list(
  results = list()
)

for (m in intersect(methods, names(models))){
    
  print_debug(m)
  
  testing_res$results[[m]] <- sl_pipeline_test(
      sldata = sldata,
      method = m,
      seed = .SEED,
      model = models[[m]]
    )
  
  # Save the time flag that identifies the model version used
  testing_res[['models_last_update']] <-  models[['last_update']]
  
  # Save the time of last testing
  testing_res[['last_update']] <-  Sys.time()
  
  # If requested, save the result on rds and excel
  if (.SAVE){
    saveRDS(object = testing_res, testing_file)
    
    excel_path <- path_loader$get_classifier_file_path(
                    type = paste(m, 'sl', sep = '_'),
                    fs_type = .FS_TYPE,
                    tuned = .TUNE,
                    path_type = 'testing',
                    testing_folder = .DATA,
                    extension = '.xlsx'
                  )
    res_to_save <- prepare_excel_res_sl(testing_res$results[[m]], sldata$test_ref)
    save_data_list(res_to_save, excel_path, names(res_to_save))
  }
}


# Comparison of metrics ---------------------------------------------------

comparison_path <- path_loader$get_classifier_file_path(
                      type = paste(.DATA, .CONFIDENT_ONLY, 'sl_comparison', sep = '_'),
                      fs_type = .FS_TYPE,
                      tuned = .TUNE,
                      path_type = 'testing',
                      testing_folder = 'comparison'
                    )

comparison <- compare_metrics(res = testing_res$results, type = 'sl')

if (.SAVE){
  
  # Save RDS file with comparison and info on published models (if used or not)
  saveRDS(object = list(comparison = comparison, 
                        published_models = .PUBLISHED_MODELS), 
          file = comparison_path)
  
  # Save excel version
  save_data_list(data_list = comparison, 
                 path_xlsx = gsub(x = comparison_path,pattern = '.rds', replacement = '.xlsx', fixed = TRUE),
                 sheet_names = names(comparison))
}  




