library(here)
source(here('src','load_libraries.r'))
source(here('src','utils','source_utils.r'))
source(here('src','loader_writer','load_data.r'))
source(here('src','data_management','source_data_management.r'))
source(here('src','classifiers','source_classifiers.r'))
source(here('src','pipelines','source_pipelines.r'))
source(here('src','biological_validation','source_biological_validation.r'))


# Data reading ------------------------------------------------------------

# Annotation and time window
annot_tcga <- load_file(path_loader$get_path('GMQL_GRCH38_ANNOT'))
annot_pdx  <- load_file(path_loader$get_path('PDX_MERGED_ANNOT'))

# Objects that allows to draw kaplan meier
biopl_tcga <- BiologicalPlots$new(annotations = annot_tcga, xlim = .XLIM)
biopl_pdx <- BiologicalPlots$new(annotations = annot_pdx, xlim = .XLIM)

ann_tcga <- AnnotationTCGA$new()
ann_pdx  <- AnnotationPDX$new()

# Data reading ----------------------------------------------------------------
sldata <- list(
  tcga = load_prepared_tcga_data(confident = 'conf',
                                   uniformed = TRUE,
                                   fs_type = 'ntp_only',
                                   type = 'sl'),
  pdx = load_prepared_pdx_data(confident = 'conf',
                                   uniformed = TRUE,
                                   fs_type = 'ntp_only',
                                   type = 'sl')
)

mldata <-  list(
  tcga = load_prepared_tcga_data(confident = 'conf',
                                   uniformed = TRUE,
                                   fs_type = 'ntp_only',
                                   type = 'ml'),
  pdx = load_prepared_pdx_data(confident = 'conf',
                                   uniformed = TRUE,
                                   fs_type = 'ntp_only',
                                   type = 'ml')
)

# Results reading ------------------------------------------------------------

# Single-label
sl_testing_tcga <- load_file(
  path = path_loader$get_classifier_file_path(
            type    = 'tcga_sl', 
            fs_type = .FS_TYPE, 
            tuned =  .TUNE, 
            path_type = 'testing', 
            testing_folder = 'tcga')
  )

sl_testing_pdx <- load_file(
  path = path_loader$get_classifier_file_path(
            type    = 'pdx_sl', 
            fs_type = .FS_TYPE, 
            tuned =  .TUNE, 
            path_type = 'testing', 
            testing_folder = 'pdx')
  )

# Multi-label algorithm adaptation
ml_adapted_testing_tcga <- load_file(
  path = path_loader$get_classifier_file_path(
            type    = 'tcga_ml_alg_adapted', 
            fs_type = .FS_TYPE, 
            tuned =  .TUNE, 
            path_type = 'testing', 
            testing_folder = 'tcga')
  )

ml_adapted_testing_pdx <- load_file(
  path = path_loader$get_classifier_file_path(
            type    = 'pdx_ml_alg_adapted',
            fs_type = .FS_TYPE, 
            tuned =  .TUNE, 
            path_type = 'testing', 
            testing_folder = 'pdx')
  )

# Statistical tests for NTP -----------------------------------------------

tests_ntp_sl_tcga <- ann_tcga$perform_stat_tests(ann_tcga$get_count_annot_sl(sldata$tcga$test_ref, annot_tcga, 'ntp_sl'))
tests_ntp_sl_pdx  <- ann_pdx$perform_stat_tests(ann_pdx$get_count_annot_sl(sldata$pdx$test_ref, annot_pdx, 'ntp_sl'))

tests_ntp_ml_tcga <- ann_tcga$perform_stat_tests(ann_tcga$get_count_annot_ml(res = NULL, annot_tcga, mldata$tcga$ref, type = 'ntp_ml'))
tests_ntp_ml_pdx  <- ann_pdx$perform_stat_tests(ann_pdx$get_count_annot_ml(res = NULL, annot_pdx, mldata$pdx$ref, type = 'ntp_ml'))


ntp_sl_forest_plots <- biopl_tcga$compute_all_forest_plots(.FOREST_ATTRIBUTES$tcga, 
                                                           .FOREST_LABEL_ATTRIBUTES$tcga, 
                                                           .FOREST_CLASSES, 
                                                           tests_ntp_sl_tcga)

ntp_sl_forest_plots <- ntp_sl_forest_plots %>% 
  rbind(biopl_pdx$compute_all_forest_plots(.FOREST_ATTRIBUTES$pdx, 
                                          .FOREST_LABEL_ATTRIBUTES$pdx, 
                                          .FOREST_CLASSES, 
                                          tests_ntp_sl_pdx,
                                          max(ntp_sl_forest_plots$Y) + 1))

ntp_ml_forest_plots <- biopl_tcga$compute_all_forest_plots(.FOREST_ATTRIBUTES$tcga, 
                                                           .FOREST_LABEL_ATTRIBUTES$tcga, 
                                                           .FOREST_CLASSES, 
                                                           tests_ntp_ml_tcga)

ntp_ml_forest_plots <- ntp_ml_forest_plots %>% 
  rbind(biopl_pdx$compute_all_forest_plots(.FOREST_ATTRIBUTES$pdx, 
                                          .FOREST_LABEL_ATTRIBUTES$pdx, 
                                          .FOREST_CLASSES, 
                                          tests_ntp_ml_pdx,
                                          max(ntp_ml_forest_plots$Y) + 1))

# Statistical tests for single-label ------------------------------------------

tests_sl_tcga <- ann_tcga$perform_stat_tests(ann_tcga$get_count_annot_sl(sl_testing_tcga$results$svmLinear2, annot_tcga, 'sl'))
tests_sl_pdx  <- ann_pdx$perform_stat_tests(ann_pdx$get_count_annot_sl(sl_testing_pdx$results$svmLinear2, annot_pdx, 'sl'))


sl_forest_plots <- biopl_tcga$compute_all_forest_plots(.FOREST_ATTRIBUTES$tcga, 
                                                           .FOREST_LABEL_ATTRIBUTES$tcga, 
                                                           .FOREST_CLASSES, 
                                                           tests_sl_tcga)

sl_forest_plots <- sl_forest_plots %>% 
  rbind(biopl_pdx$compute_all_forest_plots(.FOREST_ATTRIBUTES$pdx, 
                                          .FOREST_LABEL_ATTRIBUTES$pdx, 
                                          .FOREST_CLASSES, 
                                          tests_sl_pdx,
                                          max(sl_forest_plots$Y) + 1))

# Statistical tests for multi-label ------------------------------------------

tests_ml_tcga <- ann_tcga$perform_stat_tests(ann_tcga$get_count_annot_ml(ml_adapted_testing_tcga$results$svmLinear2, annot_tcga, mldata$tcga$ref, type = 'aa_ml'))
tests_ml_pdx  <- ann_pdx$perform_stat_tests(ann_pdx$get_count_annot_ml(ml_adapted_testing_pdx$results$svmLinear2, annot_pdx, mldata$pdx$ref, type = 'aa_ml'))



ml_forest_plots <- biopl_tcga$compute_all_forest_plots(.FOREST_ATTRIBUTES$tcga, 
                                                           .FOREST_LABEL_ATTRIBUTES$tcga, 
                                                           .FOREST_CLASSES, 
                                                           tests_ml_tcga)

ml_forest_plots <- ml_forest_plots %>% 
  rbind(biopl_pdx $compute_all_forest_plots(.FOREST_ATTRIBUTES$pdx, 
                                          .FOREST_LABEL_ATTRIBUTES$pdx, 
                                          .FOREST_CLASSES, 
                                          tests_ml_pdx,
                                          max(ml_forest_plots$Y) + 1))

# openxlsx::write.xlsx(sl_forest_plots, here('sl_forest_plots.xlsx'))
# openxlsx::write.xlsx(ntp_sl_forest_plots, here('ntp_sl_forest_plots.xlsx'))
# openxlsx::write.xlsx(ntp_ml_forest_plots[which(ntp_ml_forest_plots$Y %% 2 != 0), ], here('ntp_ml_primary_forest_plots.xlsx'))
# openxlsx::write.xlsx(ntp_ml_forest_plots[which(ntp_ml_forest_plots$Y %% 2 == 0), ], here('ntp_ml_all_forest_plots.xlsx'))
# openxlsx::write.xlsx(ml_forest_plots[which(ml_forest_plots$Y %% 2 == 0), ], here('ml_aa_all_forest_plots.xlsx'))
# openxlsx::write.xlsx(ml_forest_plots[which(ml_forest_plots$Y %% 2 != 0), ], here('ml_aa_primary_forest_plots.xlsx'))
# 





