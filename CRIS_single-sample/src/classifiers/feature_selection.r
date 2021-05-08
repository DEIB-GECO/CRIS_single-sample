
# Description -------------------------------------------------------------

# Class with set of functions for feature selection. All functions receive
# data as expression set and return data as expression set (filtered).


# Feature Selection --------------------------------------------------------

FeatureManager <- R6Class(
  
  'FeatureManager',
  lock_objects = TRUE,
  lock_class   = TRUE,

  ##### Private fields #####
  
  private = list(
    .bio_features = NULL,
    .ntp_features = NULL,
    .lasso_features = NULL,
    
    .fs_type = NULL,
    
    .fs_by_list = function(data, fs_list){
      
      gf   <- GenesFilter$new()
      filt_expr <- data@assayData$exprs
      filt_expr <- gf$filter_by_list(filt_expr, fs_list)
      
      return(ExpressionSet(assayData = filt_expr))
      
    },
    
    .prepare_lasso_data = function(data, ref){
      
      # Features data (samples on rows, ordered by sample)
      fs_data <- data@assayData$exprs %>% t() %>% as.data.frame()
      fs_data <- fs_data[order(rownames(fs_data)), ] %>% as.matrix()
      
      # Target label (ordered by sample)
      target  <- ref %>%  arrange(aliquot_id) %>%
                 select(CLASS_LABEL) %>%
                 unlist() %>%
                 factor(CRIS_CLASSES)
                
      # Debug
      # if (!all.equal(rownames(target), rownames(fs_data)))
      #   stop('Samples of data and ref do not match.') 
      
      res <- list(
        fs_data = fs_data,
        target  = target
      )
      
      return(res)
      
    },
    
    .f_imp_lasso = function(fs_data, target){
      
      # Feature importance depending on lambda
      fs_imp <- glmnet(x = fs_data, 
                       y = target, 
                       family  = "multinomial", 
                       alpha   = 1, 
                       standardize = TRUE)
      
      fs_imp_plots <- plot(fs_imp, xvar = 'lambda')
      
      res <- list(
        fs_imp = fs_imp,
        plots  = fs_imp_plots
      )
    },
    
    .cv_lambda_lasso = function(fs_data, target){
      
      # Cross validation to select optimal lambda value
      cv_lambda <- cv.glmnet(x = fs_data, 
                            y = target, 
                            family  = "multinomial", 
                            alpha   = 1, 
                            standardize = TRUE,
                            type.measure = 'class')
      
      # Plot of misclassification error depending on lambda used
      cv_lambda_plot <- plot(cv_lambda)
      
      # Optimal lambda value
      lambda_opt <- cv_lambda$lambda.min
      
      res <- list(
        cv_lambda  = cv_lambda,
        plot       = cv_lambda_plot,
        lambda_opt = lambda_opt
      )
    },
    
    .extract_lasso_features = function(cv_lambda, lambda_opt){
      
      features  <- data.frame()
      for (c in F_CRIS_CLASSES) {
        
        # Weights of features
        weights <- as.matrix(coef(cv_lambda, lambda_opt)[[as.character(c)]])
        
        # Select features with non-zero weight, ignoring the intercept
        interc_index  <- which(rownames(weights) == '(Intercept)')
        cl_features   <- subset(weights, weights[-interc_index,1] > 0) %>% 
                         rownames() %>% 
                         data.frame()
        
        # Add class label and number
        cl_features   <- cl_features %>% 
                         mutate(class_label = as.character(c)) %>%
                         data.frame()
        
        # Set colnames
        colnames(cl_features) <- c('Gene ID', CLASS_LABEL)
        
        # Add to feature list result
        features <- rbind(features, cl_features) 
        
      }
      
      features[,CLASS_LABEL] <- factor(features[,CLASS_LABEL], CRIS_CLASSES)
      features <- features %>% mutate(Class = as.integer(features[,CLASS_LABEL]))
      return(features)
      
    }
    
  ),
  
  ##### Public fields #####
  
  public = list(
    
    initialize = function(ntp_features = NULL, bio_features = NULL, lasso_features = NULL){
      
      private$.bio_features   <- bio_features
      private$.ntp_features   <- ntp_features
      private$.lasso_features <- lasso_features

    },
    
    get_features = function(fs_type){
      
      sel_genes <- switch(paste(fs_type, collapse = '_'),
                   'bio_driven_lasso' = 
                     intersect(fs$bio_features$human, fs$lasso_features$`Gene ID`),
                   'bio_driven' = 
                     fs$bio_features,
                   'bio_driven_ntp_only' = 
                     intersect(fs$bio_features$human, fs$ntp_features$`Gene ID`),
                   'lasso' = fs$lasso_features,
                   'lasso_ntp_only' = 
                     intersect(fs$lasso_features$`Gene ID`, fs$ntp_features$`Gene ID`),
                   'bio_driven_lasso_ntp_only' = 
                     intersect(fs$bio_features$human, fs$lasso_features$`Gene ID`) %>% 
                     intersect(fs$ntp_features$`Gene ID`)
                   )

      return(sel_genes)
    },
    
    feature_selection = function(data, fs_type, ref = NULL){
      
      private$.fs_type <- fs_type
      
      # Check type
      if (!any(fs_type %in% c('ntp_only', 'bio_driven', 'none','lasso')))
        stop('Supported types: ntp_only, bio_driven, none, lasso')
      
      
      # Get expression matrix and filter
      filt_expr <- data

      # Filtering (towards higher selectivity)
      if (any(fs_type == 'none'))
        return(filt_expr)
      
      if (any(fs_type == 'lasso') & check_type(private$.lasso_features,'data.frame',c(1,1))){
        print_debug('lasso')
        filt_expr <- private$.fs_by_list(filt_expr, unlist(private$.lasso_features$`Gene ID`))

      }else if (any(fs_type == 'bio_driven') & check_type(private$.bio_features,'data.frame',c(1,1))){
        print_debug('bio_driven')
        filt_expr <- private$.fs_by_list(filt_expr, unlist(private$.bio_features$human))
      }
      
      if (any(fs_type == 'ntp_only')  & check_type(private$.ntp_features,'data.frame',c(1,1))){
        print_debug('ntp_only')
        filt_expr <- private$.fs_by_list(filt_expr, unlist(private$.ntp_features$`Gene ID`))
      }
      
      return(filt_expr)
    },
    
    apply_lasso = function(data, ref){
    
      print_info('Prepare data')
      data   <- private$.prepare_lasso_data(data,ref)
      
      print_info('Variable importance')
      fs_imp <- private$.f_imp_lasso(data$fs_data, data$target)
      
      print_info('Select lambda with CV')
      cv_lambda  <- private$.cv_lambda_lasso(data$fs_data, data$target)
      lambda_opt <- cv_lambda$lambda_opt
      
      glm_mod <- glmnet(x = data$fs_data, 
                        y = data$target, 
                        family  = "multinomial",
                        lambda  = lambda_opt,
                        alpha   = 1, 
                        standardize = TRUE)
      
      # Store remaining features by class
      features <- private$.extract_lasso_features(glm_mod, lambda_opt)
      
      # Return results
      lasso_res <- list(
        
        fs_imp         = fs_imp,
        glm_mod        = glm_mod,
        cv_lambda      = cv_lambda,
        features       = features
      )
      
      return(lasso_res)
      
    },
    
    integrate_ref = function(original_features, entrez_corr, genecol = 'Gene ID'){
  
      updated_features <- original_features
      
      for (i in seq(nrow(updated_features))) {
      
        gene <- updated_features[i, genecol] %>% as.character()
        if ( gene %in% entrez_corr[,'ORIGINAL_SYMBOL']) {
          corresp <- entrez_corr[entrez_corr[,'ORIGINAL_SYMBOL'] == gene, 'FOUND_SYMBOL']
          updated_features[i, genecol] <- corresp
        }
      
      }
       
      return(as.data.frame(updated_features)[,genecol])
      
    },
    
    restore_alias = function(features, entrez_corr, genecol = 'Gene ID'){
      
      updated_features <- features
      
      for (i in seq(nrow(updated_features))) {
      
        gene <- updated_features[i,genecol] %>% as.character()
        if ( gene %in% entrez_corr[,'FOUND_SYMBOL']) {
          corresp <- entrez_corr[entrez_corr[,'FOUND_SYMBOL'] == gene, 'ORIGINAL_SYMBOL']
          updated_features[i,genecol] <- corresp
        }
      
      }
       
      return(as.data.frame(updated_features)[,genecol])
      
    }
  ),
  
  ##### Active fields #####
  
  active = list(
    
    ntp_features = function(v){
      if (missing(v))
        private$.ntp_features
      else
        stop('ntp_features is read-only')
    },
    
    bio_features = function(v){
      if (missing(v))
        private$.bio_features
      else
        stop('bio_features is read-only')
    },
    
    lasso_features = function(v){
      if (missing(v))
        private$.lasso_features
      else
        stop('lasso_features is read-only')
    },
    
    fs_type = function(v){
      if (missing(v))
        private$.fs_type
      else
        stop('fs_type is read_only')
    }
    
  ),
  
)


# FS for TCGA -------------------------------------------------------------

FeatureManagerTCGA <- R6Class(
  
  'FeatureManagerTCGA',
  inherit = FeatureManager,
  lock_objects = TRUE,
  lock_class   = TRUE,
  
  private = list(
    .load_lasso = function(features_type = NULL){
      
      lasso <- NULL
      
      if (all(c('bio_driven','lasso') %in% features_type)){
        print_info('bio lasso')
        lasso  <- load_file(path_loader$get_path("BIO_LASSO_TCGA"))
      }
      
      return(lasso)
    }
  ),
  
  public = list(
    initialize = function(lasso_features = NULL){
      
      # Features remaining after biological feature selection
      bio <- load_file(path_loader$get_path("BIO_DRIVEN_TCGA"))
      # NTP features
      load_features_grch38('tcga')
      ntp <- features_grch38
      
      # Init with superclass
      super$initialize(ntp, bio, private$.load_lasso(lasso_features))
      
    },
    
    tcga_alias_table = function(tcga_cpm){
     
        entrez_corr_pdx  <- load_file(path_loader$get_path("ENTREZ_CORR_PDX"))
        entrez_corr_tcga <- load_file(path_loader$get_path("ENTREZ_CORR_TCGA"))
        
        # Get TCGA grch38 alias
        alias_table <- data.frame(tcga = rownames(tcga_cpm))
        rownames(alias_table) <- alias_table[,'tcga']
        
        # Add TCGA hg19 alias
        alias_table <- alias_table %>%
          cbind(self$restore_alias(alias_table, entrez_corr_tcga, 'tcga')) %>%
          as.data.frame()
        colnames(alias_table) <- c('tcga','original')

        # Add PDX alias
        alias_table <- alias_table %>%
          cbind(self$integrate_ref(alias_table, entrez_corr_pdx, 'original')) %>%
          as.data.frame()
        colnames(alias_table) <- c('tcga','original','pdx')
        
        # Sort rows
        alias_table <- alias_table[rownames(alias_table), ]
        return(alias_table)
    }
    
  )
)


# FS for pdx --------------------------------------------------------------

FeatureManagerPDX <- R6Class(
  
  'FeatureManagerPDX',
  inherit = FeatureManager,
  lock_objects = TRUE,
  lock_class   = TRUE,
   private = list(
    .load_lasso = function(features_type = NULL){
      
      lasso <- NULL
      
      if (all(c('bio_driven','lasso') %in% features_type)){
        print_info('bio lasso')
        lasso  <- load_file(path_loader$get_path("BIO_LASSO_PDX"))
      }
      
      return(lasso)
    }
  ),
  
  public = list(
    initialize = function(lasso_features = NULL){
      
      # Features remaining after biological feature selection
      bio <- load_file(path_loader$get_path("BIO_DRIVEN_PDX"))
      # NTP features
      load_features_grch38('pdx')
      ntp <- features_pdx
      
      # Init with superclass
      super$initialize(ntp, bio, private$.load_lasso(lasso_features))
      
    },
    
    tcga_alias_table = function(pdx_cpm){
        
        # Gene aliases correspondences
        entrez_corr_pdx  <- load_file(path_loader$get_path("ENTREZ_CORR_PDX"))
        entrez_corr_tcga <- load_file(path_loader$get_path("ENTREZ_CORR_TCGA"))
        
        # Get TCGA grch38 alias
        alias_table <- data.frame(pdx = rownames(pdx_cpm))
        rownames(alias_table) <- alias_table[,'pdx']
        
        # Add TCGA hg19 alias
        alias_table <- alias_table %>%
          cbind(self$restore_alias(alias_table, entrez_corr_pdx, 'pdx')) %>%
          as.data.frame()
        colnames(alias_table) <- c('pdx','original')

        # Add PDX alias
        alias_table <- alias_table %>%
          cbind(self$integrate_ref(alias_table, entrez_corr_tcga, 'original')) %>%
          as.data.frame()
        colnames(alias_table) <- c('pdx','original','tcga')
        
        # Sort rows
        alias_table <- alias_table[rownames(alias_table), ]
        return(alias_table)
    }
      
  )
)
