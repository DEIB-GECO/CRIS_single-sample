
# Description -------------------------------------------------------------

# Classes for slgorithm settings (problem transformation) and Cross Validation Settings.

# Problem Transformation Settings  ----------------------------------------

AlgSettings <- R6Class(
  'AlgSettings',
  lock_objects = FALSE,
  lock_class   = TRUE,
  
  ##### Private fields #####
  
  private = list(
    .strat        = NULL,                 # chosen strategy
    .alg          = NULL,                 # chosen base algorithm
    .params       = NULL,                 # names of parameters to tune
    .tune_mode    = NULL,                 # tuning mode
    .other_params = NULL,                 # additional parameters (not tuned)
    .tune_grid    = NULL,                 # tune grid with combinations to attempt during tuning
    .supp_str     = c('br','cc','ecc'),   # supported strategies
    .supp_alg     = c('SVM','RF','KNN'),  # supported base algorithms
    
    #' Grid search algorithm applied to define the set of combination parameters
    #' to be attempted during the tuning (i.e. all possible combination of parameter
    #' values).
    #' 
    #' @param params List of parameters; each element of the list has a name (
    #' the parameter name) and a vector (the set of values to attempt).
    #' @return data.frame with grid of parameters (parmaeter names on the rows,
    #' one combination for each column)
    .grid_search = function(params){
      
      grid <- data.frame()
      
      for (i in seq(length(params))) {
        
        p    <- sort(params[[i]])
        grid <- full_join(grid, data.frame(p), by = character())
        colnames(grid) <- names(params)[1:i]
      }
      
      return(t(grid))
      
    }
    
  ),
  
  public = list(
    
    #' Constructor
    #' 
    #' @param strat Character string with a single strategy
    #' @param alg   Character string with a single base algorithm
    #' @param params List with parameter values to be attempted during the tuning
    initialize = function(strat, alg, params, ...) {
      
      # Check strategy
      if (!any(strat == private$.supp_str))
        stop('AlgSettings: `strat` must be a supported strategy')
      
      if (length(strat) != 1)
        stop('AlgSettings: `strat` must be have length 1')
      
      # Check algorithm
      if (!any(alg == private$.supp_alg))
        stop('AlgSettings: `alg` must be a character')
      
      if (length(alg) != 1)
        stop('AlgSettings: `alg` must be have length 1')
      
      
      # Check parameters
      if (all(class(params) != 'list'))
        stop('AlgSettings: `params` must be a list.')
      
      private$.strat        <- strat
      private$.alg          <- alg
      private$.params       <- names(params)

      # Build the tune_grid
      if (length(params) > 0)
        private$.tune_grid <- private$.grid_search(params)

      # Optional parameters of the algorithm, not subjected to tuning
      if (length(list(...)) > 0)
        private$.other_params <- list(...) %>% simplify_list()
      else
        private$.other_params <- list()
    }
  
  ),
  
  ##### Active fields ##### 
  # Read-only access to private fields
  
  active = list(
    alg          = function(v){
      if (missing(v))
        private$.alg
      else
        stop('AlgSettings: `alg` is read-only')
    },
    
    params       = function(v){
      if (missing(v))
        private$.params
      else
        stop('AlgSettings:`params` is read-only')
    },
    
    tune_mode    = function(v){
      if (missing(v))
        private$.tune_mode
      else
        stop('AlgSettings:`tune_mode` is read-only')
    },
    
    other_params = function(v){
      if (missing(v))
        private$.other_params
      else
        stop('AlgSettings:`other_params` is read-only')
    },
    
    tune_grid    = function(v){
      if (missing(v))
        private$.tune_grid
      else
        stop('AlgSettings:`tune_grid` is read-only')
    },
    
    supp_str    = function(v){
      if (missing(v))
        private$.supp_str
      else
        stop('AlgSettings:`supp_str` is read-only')
    },
    
    supp_alg    = function(v){
      if (missing(v))
        private$.supp_alg
      else
        stop('AlgSettings:`supp_alg` is read-only')
    },
    
    strat    = function(v){
      if (missing(v))
        private$.strat
      else
        stop('AlgSettings:`strat` is read-only')
    }
    
  )
)

# CV settings -------------------------------------------------------------

# Class that stores the Cross Validation settings (number of folds, type of
# sampling, metrics to be used to choose best tuning value)

CVSettings <- R6Class(

  'CVSettings',
  lock_objects = FALSE,
  lock_class   = TRUE,
  
  ##### Private fields #####
  
  private = list(
    .folds     = NULL,   # Number of folds for CV
    .sampling  = NULL,   # Type of sampling (stratified is better)
    .measures  = NULL,   # metrics to be used for best parameter selections
    .sel_order = NULL    # Order of selection of metrics; a '<' before the name means 
                         # the metric is good when its value is low
  ),
  
  ##### Public fields #####

  public = list(
    
    #' Constructor (TODO: CHECK)
    #' 
    #' @param folds The number of folds (default: 10)
    #' @param sampling The sampling type (default: stratified) for the fold definition
    #' @param measures The metrics to be used during tuning
    #' @param sel_order  Order of selection of metrics; a '<' before the name means 
    #' the metric is good when its value is low
    #' 
    #' @example  CVSettings$new(folds = 10, sampling = 'stratified', measures = 'hamming-loss',
    #' sel_order = '<hamming-loss')
    initialize = function(folds = 10,
                          sampling = 'stratified',
                          measures = c('hamming-loss', 'average-precision'),
                          sel_order = c('<hamming-loss', 'average-precision')) {
      
      # Check fold number
      if (class(folds) != 'numeric' | any_uncomplete(folds))
        stop('CVSettings: `folds` must be a valid number')
      
      if (length(folds) != 1 | folds < 0)
        stop('CVSettings: `folds` vector must be a vector with length 1')
      
      # Check sampling type
      if (class(sampling) != 'character' | any_uncomplete(folds))
        stop('CVSettings: `sampling` must be a valid string')
      
      if (length(sampling) != 1)
        stop('CVSettings: `sampling` must be a vector with length 1')
      
      # Check measures
      if (class(measures) != 'character' | any_uncomplete(folds))
        stop('CVSettings: `measures` must contain valid string')
      
      if (length(measures) < 1)
        stop('CVSettings: `measures` must contain at least one metric')
      
      # Check selection order
      if (class(sel_order) != 'character' | any_uncomplete(folds))
        stop('CVSettings: `sel_order` must contain valid string')
      
      if (length(sel_order) < 1)
        stop('CVSettings: `sel_order` must contain at least one metric')
      
      # Check measures and selection order coherence
      .tmp_order <- sort(gsub(x = sel_order, pattern = '<', replacement = ''))
      if(!all.equal(.tmp_order, sort(measures)))
         stop('CVSettings: `sel_order` and `measures` must contain the same metrics')
      
      # Assign the input
      private$.folds     <- as.integer(folds)
      private$.sampling  <- sampling
      private$.sel_order <- sel_order
      private$.measures  <- measures
      
    }
  ),
  
  ##### Active fields #####
  # Read-only access for private fields
  
  active = list(
    folds        = function(v){
      if (missing(v))
        private$.folds
      else
        stop('CVSettings: `folds` is read-only')
    },
    
    sampling     = function(v){
      if (missing(v))
        private$.sampling
      else
        stop('CVSettings: `sampling` is read-only')
    },
    
    sel_order    = function(v){
      if (missing(v))
        private$.sel_order
      else
        stop('CVSettings: `sel_order` is read-only')
    },
    measures    = function(v){
      if (missing(v))
        private$.measures
      else
        stop('CVSettings: `measures` is read-only')
    }
  )
  
)




