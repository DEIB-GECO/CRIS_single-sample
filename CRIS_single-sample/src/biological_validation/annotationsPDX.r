
# Description -------------------------------------------------------------

# Class for biological validation of PDX data

# Class definition --------------------------------------------------------

AnnotationPDX <- R6Class(
  
  'AnnotationPDX',
  lock_objects = TRUE,
  lock_class   = TRUE,
  inherit      = Annotation,
  
  ##### Public fields #####
  
  public = list(
    
    count_annotations = function(test_annot){
  
      # Class cardinality
      n_samples  <- test_annot %>% group_by(CRIS_class) %>% dplyr::count()
      colnames(n_samples) <- c('Class', 'n_samples')
      
      # Responsiveness cardinality
      n_sensitive <- test_annot %>% group_by(CRIS_class, responsiveness)  %>% dplyr::count() %>%
        filter(responsiveness == 'SENS') %>% ungroup() %>% select(CRIS_class, 'n')
      colnames(n_sensitive) <- c('Class', 'n_sensitive')
      
      # # Stable samples cardinality
      n_stable <- test_annot %>% group_by(CRIS_class, responsiveness)  %>% dplyr::count() %>%
        filter(responsiveness == 'SD') %>% ungroup() %>% select(CRIS_class, 'n')
      colnames(n_stable) <- c('Class', 'n_stable')
      
      n_not_resistant <- full_join(n_sensitive, n_stable, by = 'Class')
     
      for (c in colnames(n_not_resistant)[which(colnames(n_not_resistant) != 'Class')]){
        n_not_resistant[is.na(n_not_resistant[,c]),c] <- 0
      }
      n_not_resistant[,'n_sensitive'] <- n_not_resistant[,'n_sensitive'] + n_not_resistant[,'n_stable']
      n_not_resistant <- n_not_resistant[,c('Class','n_sensitive')]
      colnames(n_not_resistant)  <- c('Class', 'n_not_resistant')
      
      # Responsiveness cardinality
      n_resistant <- test_annot %>% group_by(CRIS_class, responsiveness)  %>% dplyr::count() %>%
        filter(responsiveness == 'RES') %>% ungroup() %>% select(CRIS_class, 'n')
      colnames(n_resistant) <- c('Class', 'n_resistant')
      
      # Build comparison
      
      annot_count <- n_samples %>% 
        full_join(n_sensitive, by = 'Class') %>%
        full_join(n_stable, by = 'Class') %>%
        full_join(n_resistant, by = 'Class') %>%
        full_join(n_not_resistant, by = 'Class') %>%
        as.data.frame()
      
      
      
      # NAs are turned into 0  
      for (c in colnames(annot_count)[which(colnames(annot_count)!='Class')]){
        annot_count[is.na(annot_count[,c]),c] <- 0
      }
      
      return(annot_count)
    },

    count_not_annotated = function(test_annot){
      
      # Class cardinality
      n_samples  <- test_annot %>% group_by(CRIS_class) %>% dplyr::count()
      colnames(n_samples) <- c('Class', 'n_samples')
      
      # Responsiveness cardinality
      n_sensitive <- test_annot %>% group_by(CRIS_class, responsiveness)  %>% dplyr::count() %>%
        filter(responsiveness == 'ND') %>% ungroup() %>% select(CRIS_class, 'n')
      colnames(n_sensitive) <- c('Class', 'n_sensitive')
      
      n_stable    <- n_sensitive
      colnames(n_stable) <- c('Class', 'n_stable')
      
      n_resistant <- n_sensitive
      colnames(n_resistant) <- c('Class', 'n_resistant')
      
      n_not_resistant <- n_sensitive
      colnames(n_not_resistant)  <- c('Class', 'n_not_resistant')
      
      annot_count <- n_samples %>% 
        full_join(n_sensitive, by = 'Class') %>%
        full_join(n_stable, by = 'Class') %>%
        full_join(n_resistant, by = 'Class') %>%
        full_join(n_not_resistant, by = 'Class') %>%
        as.data.frame()
      
      # NAs are turned into 0  
      for (c in colnames(annot_count)[which(colnames(annot_count)!='Class')]){
        annot_count[is.na(annot_count[,c]),c] <- 0
      }
      
      return(as.data.frame(annot_count))
      
    },

    perform_stat_tests = function(count_annot){
      
      tests <- list()

      for (c in names(count_annot)){
        print_info(paste(c, '-------------------------------------------------------'))
        tests[[c]] <- list()
        tests[[c]][['n_sensitive']]   <-  list()
        for (cl in count_annot[[c]]$Class){
          
          if (!grepl(cl, pattern = '^not CRIS-[A|B|C|D|E]$', fixed = FALSE)){
              tests[[c]][['n_sensitive']][[cl]]   <-  self$validate_bio(count_annot[[c]], cl, 'n_sensitive')
              tests[[c]][['n_not_resistant']][[cl]]   <-  self$validate_bio(count_annot[[c]], cl, 'n_not_resistant')
          }
          
        }
        
      }
      
      return(tests) 
    }

  )

)
