# Description -------------------------------------------------------------

# Check classifiers are actually single-sample

# Dependencies ------------------------------------------------------------
library(here)
source(here('src','load_libraries.r'))
source(here('src','utils','source_utils.r'))
source(here('src','loader_writer','load_data.r'))
source(here('src','data_management','source_data_management.r'))
source(here('src','classifiers','source_classifiers.r'))
source(here('src','pipelines','source_pipelines.r'))


# Data Loading -----------------------------------------------------------------

# Models and thresholds
sl_models <- load_file(path_loader$get_path('NTP_ONLY_SL_MODELS'))
sl_thr    <- load_file(path_loader$get_path('NTP_ONLY_ML_AA_THR'))

ml_models <- load_file(path_loader$get_path('NTP_ONLY_ML_MODELS'))
ml_thr    <- load_file(path_loader$get_path('NTP_ONLY_ML_PT_THR'))

# Complete testing for problem transformation
method       <- 'tcga_ml_problem_transformation'
testing_file <- path_loader$get_classifier_file_path(method, .FS_TYPE, .TUNE, path_type = 'testing', testing_folder = 'tcga')
ml_pt_testing_ref <- load_file(testing_file)

# Complete testing for algorithm_adaptation
method       <- 'tcga_ml_alg_adapted'
testing_file <- path_loader$get_classifier_file_path(method, .FS_TYPE, .TUNE, path_type = 'testing', testing_folder = 'tcga')
ml_aa_testing_ref <- load_file(testing_file)

# Complete testing for single-label
method       <- 'tcga_sl'
testing_file <- path_loader$get_classifier_file_path(method, .FS_TYPE, .TUNE, path_type = 'testing', testing_folder = 'tcga')
sl_testing_ref <- load_file(testing_file)


testing_samples <- load_tcga_testing_samples()
s_filt <- list()

for (i in 1:10){
  s_filt[[i]] <- testing_samples[sample.int(n = length(testing_samples),size = 3)]
}


# Function for checking
is_single_sample <- function(res, reference_res, model_name, label, score_tolerance = 1e-7){
  
  samples <- rownames(res$binary_res)
  
  check_scores <- all.equal.numeric(res$pred[samples, CRIS_CLASSES], 
                                    reference_res$results[[model_name]]$pred[samples, CRIS_CLASSES], 
                                    tolerance = score_tolerance)
  
  check_binary <- all.equal(res$binary_res[samples, ],
                            reference_res$results[[model_name]]$binary_res[samples, ])
  
  if (any(check_binary == TRUE) & any(check_scores == TRUE)){
    print_success(label)
    return(TRUE)
  }else{
    warning(paste(label, 'results are not the same'))
    return(FALSE)
  }
  
}
# Testing -----------------------------------------------------------------

i <- 1
sl_m <- methods[[1]]
ml_m <- 'ecc_lsvm'


ss_test <- list()

for (i in 1:length(s_filt)){
  
  sldata <- load_prepared_tcga_data(confident = .CONFIDENT_ONLY, 
                                    uniformed = .UNIFORMED, 
                                    fs_type   = .FS_TYPE, 
                                    type      = 'sl',
                                    samples_filter = s_filt[[1]])

  mldata <- load_prepared_tcga_data(confident = .CONFIDENT_ONLY, 
                                   uniformed = .UNIFORMED, 
                                   fs_type   = .FS_TYPE, 
                                   type      = 'ml',
                                   samples_filter = s_filt[[1]],
                                   load_training = FALSE)


  sl_testing <- sl_pipeline_test(
                sldata = sldata,
                method = methods[sl_m],
                seed = .SEED,
                model = sl_models[[sl_m]]
              )


  ml_aa_testing <- sl_pipeline_test(
                  sldata = sldata,
                  method = sl_m,
                  seed = .SEED,
                  model = sl_models[[sl_m]],
                  cl_thresholds = sl_thr[[sl_m]]$thresholds,
                  mldata = mldata,
                  max_cls = sl_thr[[sl_m]]$max_cls,
                  min_cls = sl_thr[[sl_m]]$min_cls
                )
  
  ml_pt_testing <- ml_pipeline_test(
        mldata = mldata,
        seed = .SEED,
        cv_set = CVSettings$new(),
        model = ml_models[[ml_m]],
        cl_thresholds = ml_thr[[ml_m]]
      )

  
  
  ss_test[[i]] <- list(
    sl = is_single_sample(sl_testing, sl_testing_ref, sl_m, 'lsvm_sl'),
    ml_aa = is_single_sample(ml_aa_testing, ml_aa_testing_ref, sl_m, 'lsvm_ml_aa'),
    ml_pt = is_single_sample(ml_pt_testing, ml_pt_testing_ref, ml_m, 'ecc_lsvm')
  )
  
}

