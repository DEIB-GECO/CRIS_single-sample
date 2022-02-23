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

# Path file with all single-label classifier models
if (any(.PUBLISHED_MODELS)) {
  .model_file <- path_loader$get_path('NTP_ONLY_SL_MODELS')
  .thresholds_file  <- path_loader$get_path('NTP_ONLY_ML_AA_THR')
}else{
  .model_file      <- path_loader$get_classifier_file_path('sl', .FS_TYPE, .TUNE, path_type = 'models')
  .thresholds_file <- path_loader$get_classifier_file_path('sl', .FS_TYPE, .TUNE, path_type = 'thresholds')
}

method <- paste(.DATA, 'ml_alg_adapted', sep = '_')
testing_file <- path_loader$get_classifier_file_path(method,
                                                     .FS_TYPE,
                                                     .TUNE,
                                                     path_type = 'testing',
                                                     testing_folder = .DATA)

# Classifier settings -----------------------------------------------------

# Lists for saving models that have been trained
models <- load_file(.model_file)
thresholds <- load_file(.thresholds_file)

# Load data with settings specified in the config.r file


if (.DATA == 'tcga'){
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
}else if (.DATA == 'pdx'){
  
  sldata <- load_prepared_pdx_data(confident = .CONFIDENT_ONLY, 
                                    uniformed = .UNIFORMED, 
                                    fs_type   = .FS_TYPE, 
                                    type      = 'sl')
  
  mldata <- load_prepared_pdx_data(confident = .CONFIDENT_ONLY, 
                                    uniformed = .UNIFORMED, 
                                    fs_type   = .FS_TYPE, 
                                    type      = 'ml')
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
      model = models[[m]],
      cl_thresholds = thresholds[[m]]$thresholds,
      mldata = mldata,
      max_cls = thresholds[[m]]$max_cls,
      min_cls = thresholds[[m]]$min_cls
      
    )
  
  # Save the time flag that identifies the model version used
  testing_res[['models_last_update']] <-  models[['last_update']]
  
  # Save the time of last testing
  testing_res[['last_update']] <-  Sys.time()
  
  # Save the time of last testing
  testing_res[['published_models']] <- .PUBLISHED_MODELS
  
  # If requested, save the results 
  if (.SAVE){
    
    # Save all the results in rds file
    saveRDS(object = testing_res, testing_file)
    
    # Save model-specific results in excel file
    res_excel <- prepare_excel_res_ml(testing_res$results[[m]], mldata$test_ref)
    res_excel_path <- path_loader$get_classifier_file_path(paste(m, 'ml_alg_adapted', sep = '_'),
                                                           .FS_TYPE,
                                                           .TUNE,
                                                           path_type = 'testing',
                                                           testing_folder = .DATA,
                                                           extension = '.xlsx')
    save_data_list(data_list = res_excel, 
                   path_xlsx = res_excel_path,
                   sheet_names = names(res_excel))
  }
}



# Comparison of metrics ---------------------------------------------------

# RDS path
comparison_path <- path_loader$get_classifier_file_path(
                      type = paste(.DATA, .CONFIDENT_ONLY, 'ml_alg_adapted_comparison', sep = '_'),
                      fs_type = .FS_TYPE,
                      tuned = .TUNE,
                      path_type = 'testing',
                      testing_folder = 'comparison'
                    )

# Compare the results
comparison <- compare_metrics(res = testing_res$results, type = 'ml')

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

  
  

