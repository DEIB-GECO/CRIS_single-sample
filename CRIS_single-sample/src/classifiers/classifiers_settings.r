# Description -------------------------------------------------------------

# Script with constants and settings for the proposed models

# General settings --------------------------------------------------------

# Replication seed
.SEED       <- 598

# Supported models --------------------------------------------------------

# Multi label strategies and base algorithms
STRATEGIES <- c('br','cc','ecc')
ALGORITHMS <- c('lsvm','rsvm','psvm','rf')

# Single Label models
methods    <- c('svmLinear2', 'svmRadial', 'svmPoly','rf')

# Single Label settings ----------------------------------------------------

# Formula for the generic classifier
.SL_FORMULA <- as.formula(paste(CLASS_LABEL,'~ .', sep = ''))

# Cross-Validation settings for caret
fit_control  <- trainControl(
  method = "cv", 
  number = 10,
  classProbs = TRUE
)

# Tuning values for each single-label classifier
tune_grid <- list(
  
  svmLinear2 = data.frame(cost   = c(1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3)),
  svmRadial  = data.frame(C      = c(1e-2,1e-1,1e0,1e1,1e2),
                          sigma  = c(1e-2,1e-1,1e0,1e1,1e2)),
  # https://r-stat-sc-donnees.github.io/code_html/SVM.nb.html
  # https://datascience.stackexchange.com/questions/1074/polynomial-kernel-parameters-in-svms
  svmPoly    = data.frame(C      = c(1e-2,1e0,1e2),
                          degree = c(1,2,3),
                          scale  = 1),
  logreg     = data.frame(),
  rf         = data.frame(mtry  = c(10,20,50,100,200))
  
)

# Multi-label settings ---------------------------------------------------------

# Base algorithms parameters ---------------------------------------------

# tuning: hyperparameters subjected to tuning
# fixed:  fixed hyperparameters
.p_lsvm <- list(
  tuning = list(cost = c(1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3)),
  fixed  = list(
    kernel = 'linear'
  )
)

.p_rsvm <- list(
  tuning = list(cost  = c(1e-2,1e-1,1e0,1e1,1e2),
                gamma = c(1e-2,1e-1,1e0,1e1,1e2)),
  fixed = list(
    kernel = 'radial'
  )
)


.p_psvm <- list(
  tuning = list(cost = c(1e-2,1e-1,1e0,1e1,1e2),
                gamma = c(1e-2,1e-1,1e0,1e1,1e2),
                degree = c(2,3)),
  fixed = list(
    kernel = 'polynomial'
  )
)


.p_rf <- list(
  tuning = list(ntree = c(250, 500, 800)),
  fixed  = list()
)


# Strategy parameters -----------------------------------------------------

.p_cc <- list(
  tuning = list(),
  fixed  = list(chain = levels(F_CRIS_CLASSES)[c(3,1,5,2,4)])
)

.p_ecc <- list(
  tuning = list(m = c(60,120)),
  fixed = list(attr.space = 1)
)

.p_br <- list(
  tuning = list(),
  fixed = list()
)

# Settings definition -----------------------------------------------------

settings <- list()

for (s in STRATEGIES){
  
  str_p_var <- paste('.p',s, sep = '_')
  for (a in ALGORITHMS){
    
    # Name of setting and of base algorithm parameters
    name  <- paste(s,a, sep = '_') %>% tolower()
    
    alg_p_var <- paste('.p',a, sep = '_')
    # psvm, lsvm, rsvm are associated to SVM base algorithm
    if (grepl(a, pattern = 'svm', fixed = TRUE))
      alg <- 'SVM'
    else
      alg <- toupper(a)
    
    # Get parameters of base algorithm, if any
    if (exists(alg_p_var)){
      tuning_p <- get(alg_p_var)$tuning
      fixed_p  <- get(alg_p_var)$fixed
    }else{
      tuning_p <- list()
      fixed_p  <- list()
    }
    
    # Add parameters of strategy, if any
    if (exists(str_p_var)){
      tuning_p <- merge_lists(tuning_p, get(str_p_var)$tuning)
      fixed_p  <- merge_lists(fixed_p, get(str_p_var)$fixed)
    }
    
    if (!name %in% c('ecc_rsvm', 'ecc_psvm', 'ecc_rf')){
      # Create and add the setting to the list
      settings[[name]] <- AlgSettings$new(
        strat = s,
        alg = alg,
        params = tuning_p,
        fixed_p
      )
    }
  }
}


