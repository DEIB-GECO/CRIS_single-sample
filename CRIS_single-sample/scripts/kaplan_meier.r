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

# Objects that allows to draw kaplan meier
biopl_tcga <- BiologicalPlots$new(annotations = annot_tcga, xlim = .XLIM)

.DATA <- 'tcga'

# Data reading ----------------------------------------------------------------
sldata <- load_prepared_tcga_data(confident = 'conf',
                                   uniformed = TRUE,
                                   fs_type = 'ntp_only',
                                   type = 'sl')

mldata <- load_prepared_tcga_data(confident = 'conf',
                                 uniformed = TRUE,
                                 fs_type = 'ntp_only',
                                 type = 'ml')

# Results reading ------------------------------------------------------------

sl_testing <- load_file(
  path = path_loader$get_classifier_file_path(
            type    = paste(.DATA, 'sl', sep = '_'), 
            fs_type = .FS_TYPE, 
            tuned =  .TUNE, 
            path_type = 'testing', 
            testing_folder = .DATA)
  )

ml_adapted_testing <- load_file(
  path = path_loader$get_classifier_file_path(
            type    = paste(.DATA, 'ml_alg_adapted', sep = '_'), 
            fs_type = .FS_TYPE, 
            tuned =  .TUNE, 
            path_type = 'testing', 
            testing_folder = .DATA)
  )

ml_pt_testing <- load_file(
  path = path_loader$get_classifier_file_path(
            type    = paste(.DATA, 'ml_problem_transformation', sep = '_'), 
            fs_type = .FS_TYPE, 
            tuned =  .TUNE, 
            path_type = 'testing', 
            testing_folder = .DATA)
  )

# References for kaplan meier -----------------------------------------------

km_ntp_sl_ref <- biopl_tcga$ref_kaplan_meier(sldata$test_ref, type = 'ntp_sl', cl = 'CRIS-B')
km_ntp_ml_ref <- biopl_tcga$ref_kaplan_meier(mldata$test_ref, type = 'ntp_ml', cl = 'CRIS-B')


# Computation of the Kaplan Meier ------------------------------------------

# NTP
km_ntp_sl <- biopl_tcga$compute_kaplan_meier(km_ntp_sl_ref, group = 'ref',cl = 'CRIS-B')
km_ntp_ml <- biopl_tcga$compute_kaplan_meier(km_ntp_ml_ref, group = 'ref',cl = 'CRIS-B')
km_ntp_ml_primary <- biopl_tcga$compute_kaplan_meier(km_ntp_ml_ref, group = 'is_primary',cl = 'CRIS-B')

km_ntp_ml_secondary <- biopl_tcga$compute_kaplan_meier_additional(km_ntp_ml_ref, group = 'is_secondary',cl = 'CRIS-B')


# Single label
km_sl <- list()
for (m in names(sl_testing$results)){
  
  .sl_res <- sl_testing$results[[m]]$binary_res %>% rownames_to_column(ALIQUOT_LABEL)
  .km_ref <- biopl_tcga$ref_kaplan_meier(.sl_res, type = 'sl', cl = 'CRIS-B')
  km_sl[[m]] <- biopl_tcga$compute_kaplan_meier(.km_ref, group = 'ref', cl = 'CRIS-B')
}


# Multi-label algorithm adaptation
km_ml_ref <- list()
km_ml_primary <- list()

km_ml_secondary <- list()


for (m in names(ml_adapted_testing$results)){
  
  .ml_res <- ml_adapted_testing$results[[m]]$binary_res %>% rownames_to_column(ALIQUOT_LABEL)
  .km_ref <- biopl_tcga$ref_kaplan_meier(.ml_res, type = 'aa_ml', cl = 'CRIS-B')
  km_ml_ref[[m]] <- biopl_tcga$compute_kaplan_meier(.km_ref, group = 'ref', cl = 'CRIS-B')
  km_ml_primary[[m]] <- biopl_tcga$compute_kaplan_meier(.km_ref, group = 'is_primary', cl = 'CRIS-B')
  
  km_ml_secondary[[m]] <- biopl_tcga$compute_kaplan_meier_additional(.km_ref, group = 'is_secondary', cl = 'CRIS-B')

  
}

# Multi-label problem transformation
for (m in names(ml_pt_testing$results)){
  
  .ml_res <- ml_pt_testing$results[[m]]$binary_res %>% rownames_to_column(ALIQUOT_LABEL)
  .km_ref <- biopl_tcga$ref_kaplan_meier(.ml_res, type = 'pt_ml', cl = 'CRIS-B')
  km_ml_ref[[m]] <- biopl_tcga$compute_kaplan_meier(.km_ref, group = 'ref', cl = 'CRIS-B')
  km_ml_primary[[m]] <- biopl_tcga$compute_kaplan_meier(.km_ref, group = 'is_primary', cl = 'CRIS-B')
}


# Draw plots ---------------------------------------------------------------

km_fold <- paste(path_loader$get_path('OUT_FOLDER_MODELS'), .FS_TYPE, 'kaplan_meier', sep = '/')
dir_create(km_fold)

fname <- paste('ADDITIONAL_kaplan_meier', .XLIM, '.jpeg', sep = '')
jpeg(paste(km_fold, fname, sep = '/'), width = 800, height = 800,quality = 100)

# 6 graphs on a 3x2 grid
par(mfrow=c(4,2))

# # Margins of the plots
par(mar=c(4.1, 4.1, 4.1, 2.1))

# First row
biopl_tcga$show_km_plot(km_ntp_sl, alg_name = 'a) NTP single-label')
biopl_tcga$show_km_plot(km_sl$svmLinear2, alg_name = 'b) LSVM single-label')

# Second row
biopl_tcga$show_km_plot(km_ntp_ml_primary, alg_name = 'c) NTP multi-label - primary')
biopl_tcga$show_km_plot(km_ml_primary$svmLinear2, alg_name = 'd) LSVM multi-label - primary')

biopl_tcga$show_km_plot(km_ntp_ml_secondary, alg_name = 'e) NTP multi-label - secondary')
biopl_tcga$show_km_plot(km_ml_secondary$svmLinear2, alg_name = 'f) LSVM multi-label - secondary')


# Third row
biopl_tcga$show_km_plot(km_ntp_ml, alg_name = 'g) NTP multi-label - all')
biopl_tcga$show_km_plot(km_ml_ref$svmLinear2, alg_name = 'h) LSVM multi-label - all')


#biopl_tcga$show_km_plot(km_ml_ref$ecc_lsvm, alg_name = 'g) ecc_lsvm multi-label - all')
#biopl_tcga$show_km_plot(km_ml_primary$ecc_lsvm, alg_name = 'h) ecc_lsvm multi-label - primary')

dev.off()

fname <- paste('LSVM single-label', .XLIM, '.jpeg', sep = '')
jpeg(paste(km_fold, fname, sep = '/'), width = 800, height = 800,quality = 100)
par(1,1)
biopl_tcga$show_km_plot_and_table(km_sl$svmLinear2, alg_name = 'a) LSVM single-label')
dev.off()

fname <- paste('LSVM multi-label primary', .XLIM, '.jpeg', sep = '')
jpeg(paste(km_fold, fname, sep = '/'), width = 800, height = 800,quality = 100)
par(1,1)
biopl_tcga$show_km_plot_and_table(km_ml_primary$svmLinear2, alg_name = 'c) LSVM multi-label - primary')
dev.off()

fname <- paste('LSVM multi-label secondary', .XLIM, '.jpeg', sep = '')
jpeg(paste(km_fold, fname, sep = '/'), width = 800, height = 800,quality = 100)
par(1,1)
biopl_tcga$show_km_plot_and_table(km_ml_secondary$svmLinear2, alg_name = 'd) LSVM multi-label - secondary')
dev.off()

fname <- paste('LSVM multi-label all', .XLIM, '.jpeg', sep = '')
jpeg(paste(km_fold, fname, sep = '/'), width = 800, height = 800,quality = 100)
par(1,1)
biopl_tcga$show_km_plot_and_table(km_ml_ref$svmLinear2, alg_name = 'b) LSVM multi-label - all')
dev.off()
