# Description -------------------------------------------------------------

# Class for the execution of NTP classifier. 

# Class definition ------------------------------------------------------------

.DEF_N_RESEMPL <- 100000
.DEF_OUTPUT    <- here('output/Replication')

NTPClassifier <- R6Class(
  'NTPClassifier',
  lock_objects = TRUE,
  lock_class   = TRUE,

  #### Private fields ####
  # Classifier parameters
  private = list(
    
    .output_folder = NULL,
    .features      = NULL,
    .n_resempl     = NULL,
    
    # TODO: check format of path is valid os path
    .validate_path = function(path){
      
      # Output folder path
      
      if (missing(path))
        stop('`path` cannot be empty')
      
      if (!any(class(path) == 'character'))
        stop('`path` must be a character vector')
      
      if (length(path) != 1)
        stop('`path` must be a character vector of length 1')
      
      
    },
    
    .validate_id   = function(id){
      
      # Id of the replication
      
      if (missing(id))
        stop('`id` cannot be empty')
      
      if (!any(class(id) %in% c('character','numeric')))
        stop('`id` must be a numeric or character vector')
      
      if (length(id) != 1)
        stop('`id` must be a numeric or character vector of length 1')
      
    },
    
    .validate_data = function(data){
      
      if (class(data) != 'ExpressionSet')
        stop('`data` must be an ExpressionSet object')
    },
    
    # Assume assayData$exprs has a matrix with samples on columns and
    # genes on the rows (with correct names)
    .exprset_to_gct = function(expr_set){
  
      # Get expression matrix
      data    <- expr_set@assayData$exprs %>% as.data.frame()
      genes   <- rownames(data)
      samples <- colnames(data)
      
      # Format for GCT (genes in first 2 columns)
      gct_data           <- cbind(genes, genes, data)
      colnames(gct_data) <- c('NAMES','NAMES', samples) 
      
      # Heading
      version <- c('#1.0')
      gct_dim <- paste(dim(gct_data), collapse = '\t')
      heading <- c(version, gct_dim)
      
      # Colnames
      gct_cols <- paste(colnames(gct_data), collapse = '\t')
      
      return(list(
        heading = heading,
        cols    = gct_cols,
        data    = gct_data
      ))
      
    },
    
    .save_gct_data = function(gct_data, gct_path){
      
      # Write gct heading (version + dim)
      writeLines(gct_data$heading, gct_path)
      
      # Append colnames of gct data
      write(x      = gct_data$cols,
            file   = gct_path,
            append = TRUE)
      
      # Append data
      write.table(
        gct_data$data,
        gct_path,
        append = TRUE,
        col.names = FALSE,
        row.names = FALSE,
        quote     = FALSE,
        sep       = '\t'
      )
       
    },
    
    .ntp_classification = function(in_path, out_path){
      
      # Paths
      private$.validate_path(in_path)
      private$.validate_path(out_path)

      # Classify
      NTPez(input.exp.filename = in_path, 
            output.name        = out_path, 
            input.features     = self$features, 
            nresmpl            = self$n_resempl)
        
      # Result
      res_path <- paste(out_path, "_prediction_result.xls", sep = "")
      ntp_res  <- read.table(res_path, header = TRUE)
      
      if(!is.null(ntp_res))
        colnames(ntp_res)[2] <- CLASS_LABEL
      else
        warning('Reading of NTP result returned NULL')
      
      return(ntp_res)

    }
    
  ),
  
  #### Public fields ####
  
  public = list(
    
    initialize = function(features,
                          n_resempl = .DEF_N_RESEMPL,
                          output_folder = .DEF_OUTPUT){
      
      # Check input
      
      if (all(class(features) != 'data.frame'))
        stop('`features` must be a data.frame')
      if (ncol(features) != 3)
        stop('`features` must have 3 columns (gene, class, class order).')
      
      if (all(class(n_resempl) != 'numeric'))
        stop('`n_resempl` must be an integer number')
      if (n_resempl < 1 | round(n_resempl) != n_resempl)
        stop('`n_resempl` must be a positive integer')
      

      private$.validate_path(output_folder)
      dir_create(output_folder)
      
      # Assign values
      private$.output_folder <- output_folder
      private$.n_resempl     <- n_resempl
      private$.features      <- features
      
      
    },
    
    classify = function(data, id){
      
      print_info('Check input')
      
      # Check input
      private$.validate_id(id)
      private$.validate_data(data)
      
      print_info('Result folder')
      
      # Create folder for NTP result
      out_f_name <- paste('NTP', id, sep = '_')
      out_path   <- paste(self$output_folder, '/', out_f_name, '/', sep = '')
      dir_create(path_abs(out_path))
      
      print_info('Data saved into GCT')
      
      # Input adaptation
      gct_data <- private$.exprset_to_gct(data)
      in_path  <- paste(out_path,'data.gct', sep = '')
      private$.save_gct_data(gct_data, in_path)
      
      print_info('NTP classification')
      
      # Classification
      pred <- private$.ntp_classification(in_path, out_path) %>% as.data.frame()
      colnames(pred)[colnames(pred) == 'sample.names'] <- ALIQUOT_LABEL
      
      return(as.data.frame(pred))
    },
    
    metrics = function(ref_data, pred_data){
      
      # Adapt prediction
      labels    <- c(ALIQUOT_LABEL, CLASS_LABEL, CLASS_DISTANCE_LABEL)
      pred_data <- pred_data %>% select(labels) %>% as.data.frame()
      
      colnames(pred_data) <- colnames(pred_data) %>%
                             gsub(pattern = 'distCRIS.', replacement = 'CRIS-')
      
      # Correlations
      pred_data[,CRIS_CLASSES] <- 1 - pred_data[,CRIS_CLASSES]
      
      # Metrics
      sl  <- SLMetrics$new()$metrics(ref_data, pred_data)
      sim <- SimMetrics$new()$metrics(ref_data, pred_data)
      
      metrics <- list(
        sl = sl,
        sim = sim,
        avg = list(
          pearson  = mean(sim$pearson, na.rm = TRUE),
          spearman = mean(sim$spearman, na.rm = TRUE)
        )
      )
      
      return(metrics)
      
    }
    
  ),
  
  #### Active fields ####
  # Read only access to private fields
  
  active = list(
    
    output_folder = function(v){
      if (missing(v))
        private$.output_folder
      else
        stop('`output_folder` is read-only.')
    },
    
    features  = function(v){
      if (missing(v))
        private$.features
      else
        stop('`features` is read-only.')
    },
    
    n_resempl = function(v){
      if (missing(v))
        private$.n_resempl
      else
        stop('`n_resempl` is read-only.')
    }
  )
  
  
)




