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

# Path files
if (any(.PUBLISHED_MODELS)) {
  .model_file <- path_loader$get_path('NTP_ONLY_ML_MODELS')
  .thresholds_file  <- path_loader$get_path('NTP_ONLY_ML_PT_THR')
}else{
  .model_file       <- path_loader$get_classifier_file_path('ml', .FS_TYPE, .TUNE, path_type = 'models')
  .thresholds_file  <- path_loader$get_classifier_file_path('ml', .FS_TYPE, .TUNE, path_type = 'thresholds')
}

method <- paste(.DATA, 'ml_problem_transformation', sep = '_')
testing_file <- path_loader$get_classifier_file_path(method, .FS_TYPE, .TUNE, path_type = 'testing', testing_folder = .DATA)
  
# Classifier settings -----------------------------------------------------

# Lists for saving models that have been trained
models <- load_file(.model_file)
thresholds <- load_file(.thresholds_file)

# Load data with settings specified in the config.r file
if (.DATA == 'tcga'){
  mldata <- load_prepared_tcga_data(confident = .CONFIDENT_ONLY, 
                                    uniformed = .UNIFORMED, 
                                    fs_type   = .FS_TYPE, 
                                    type      = 'ml',
                                    load_training = TRUE)
}else if (.DATA == 'pdx'){
  mldata <- load_prepared_pdx_data(confident = .CONFIDENT_ONLY, 
                                    uniformed = .UNIFORMED, 
                                    fs_type   = .FS_TYPE, 
                                    type      = 'ml')
}else {
  stop('.DATA flag must be either tcga or pdx')
}

# Hold the result of the testing_pipeline
testing_res <- list(
  results = list()
)

for (m in intersect(names(settings), names(models))){
    
  print_debug(m)
  testing_res$results[[m]] <- ml_pipeline_test(
      mldata = mldata,
      seed = .SEED,
      cv_set = CVSettings$new(),
      model = models[[m]],
      cl_thresholds = thresholds[[m]]
    )
  
  # Save the time flag that identifies the model version used
  testing_res[['models_last_update']] <-  models[['last_update']]
  
  # Save the time of last testing
  testing_res[['last_update']] <-  Sys.time()
  
  # If requested, save the models and the settings
  if (.SAVE){ 
    
    # Save all the results in rds file
    saveRDS(object = testing_res, testing_file)
    
    # Save model-specific results in excel file
    res_excel <- prepare_excel_res_ml(testing_res$results[[m]], mldata$test_ref)
    res_excel_path <- path_loader$get_classifier_file_path(m,
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
                      type = paste(.DATA, .CONFIDENT_ONLY, 'ml_problem_transf_comparison', sep = '_'),
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