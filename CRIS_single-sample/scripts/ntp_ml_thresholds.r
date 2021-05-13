# Description -------------------------------------------------------------

# Computation of thresholds for assignemnt of multi-label NTP

# Dependencies ------------------------------------------------------------

library(here)
source(here('src','load_libraries.r'))
source(here('src','utils','source_utils.r'))
source(here('src','loader_writer','load_data.r'))
source(here('src','data_management','source_data_management.r'))
source(here('src','classifiers','source_classifiers.r'))
source(here('src','pipelines','source_pipelines.r'))


# Configuration -----------------------------------------------------------

# The threshold is the fifth percentile of the distribution of NTP scores 
# (NTP score = 1 - NTP distance)
.EXCLUDED_PERC <- 0.05


compute_ntp_corr_thr <- function(ref_data,excluded_perc){
  
  dist_thr <- numeric(length(CRIS_CLASSES))

  for (c in seq(length(CRIS_CLASSES))){
    
    # Current class
    cl <- levels(F_CRIS_CLASSES)[c]
    
    # Get confident samples of current class
    cl_conf_samp <- ref_data %>% filter(predict.label2 == cl) %>% 
                    filter(BH.FDR < BH_FDR_THRESHOLD)
    
    dist_thr[c] <- quantile(1 - cl_conf_samp[,CLASS_DISTANCE_LABEL[c]], excluded_perc)
    
  }
  
  names(dist_thr) <- levels(F_CRIS_CLASSES)
  
  return(dist_thr)
}

# Get data
tcga_ref <- readRDS(path_loader$get_path('NTP_REF_TCGA'))
pdx_ref  <- readRDS(path_loader$get_path('NTP_REF_PDX'))

# Threshold computation
thresholds <- list(
  tcga = compute_ntp_corr_thr(tcga_ref$result, .EXCLUDED_PERC),
  pdx  = compute_ntp_corr_thr(pdx_ref$result, .EXCLUDED_PERC)
)


