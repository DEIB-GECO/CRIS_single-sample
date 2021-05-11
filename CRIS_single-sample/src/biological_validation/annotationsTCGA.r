
# Description -------------------------------------------------------------

# Class for biological validation of TCGA data

# Class definition --------------------------------------------------------

AnnotationTCGA <- R6Class(
  
  'AnnotationTCGA',
  lock_objects = TRUE,
  lock_class   = TRUE,
  inherit      = Annotation,
  
  ##### Public fields #####
  
  public = list(
    
    count_annotations = function(test_annot){
  
      # Class cardinality
      n_samples  <- test_annot %>% group_by(CRIS_class) %>% dplyr::count()
      colnames(n_samples) <- c('Class', 'n_samples')
      
      # Mucinous cardinality
      n_mucinous <- test_annot %>% group_by(CRIS_class, is_mucinous)  %>% dplyr::count() %>%
        filter(is_mucinous == 'YES') %>% ungroup() %>% select(CRIS_class, 'n')
      colnames(n_mucinous) <- c('Class', 'n_mucinous')
      
      # Micro Satellite Instability cardinality
      n_msi_h    <- test_annot %>% group_by(CRIS_class, msi_status)  %>% dplyr::count() %>%
        filter(msi_status == 'MSI-H') %>% ungroup() %>% select(CRIS_class, 'n')
      colnames(n_msi_h) <- c('Class', 'n_msi_h')
      
      # Micro satellite stability cardinality
      n_msi_low  <- test_annot %>% group_by(CRIS_class, msi_status)  %>% dplyr::count() %>%
        filter(msi_status == 'MSI-L') %>% ungroup() %>% select(CRIS_class, 'n')
      colnames(n_msi_low) <- c('Class', 'n_msi_l')
      
      n_mss      <- test_annot %>% group_by(CRIS_class, msi_status)  %>% dplyr::count() %>%
        filter(msi_status == 'MSS') %>% ungroup() %>% select(CRIS_class, 'n')
      colnames(n_mss) <- c('Class', 'n_mss')
      
      # KRAS mutated cardinality
      n_kras_mut <- test_annot %>% group_by(CRIS_class, KRAS_mutated)  %>% dplyr::count() %>%
        filter(KRAS_mutated == 'YES') %>% ungroup() %>% select(CRIS_class, 'n')
      colnames(n_kras_mut) <- c('Class', 'kras_mutated')
      
      # BRAF mutated cardinality
      n_braf_mut <- test_annot %>% group_by(CRIS_class, BRAF_mutated)  %>% dplyr::count() %>%
        filter(BRAF_mutated == 'YES') %>% ungroup() %>% select(CRIS_class, 'n')
      colnames(n_braf_mut) <- c('Class', 'braf_mutated')
      
      # NRAS mutated cardinality
      n_nras_mut <- test_annot %>% group_by(CRIS_class, NRAS_mutated)  %>% dplyr::count() %>%
        filter(NRAS_mutated == 'YES') %>% ungroup() %>% select(CRIS_class, 'n')
      colnames(n_nras_mut) <- c('Class', 'nras_mutated')
      
      # NRAS mutated cardinality
      n_tp53_mut <- test_annot %>% group_by(CRIS_class, TP53_mutated)  %>% dplyr::count() %>%
        filter(TP53_mutated == 'YES') %>% ungroup() %>% select(CRIS_class, 'n')
      colnames(n_tp53_mut) <- c('Class', 'tp53_mutated')
      
      # Samples with recurrence before 24 months
      n_bad_prog_36 <- test_annot %>% 
        group_by(CRIS_class)  %>% 
        filter(disease_free_months_consensus < 36 & disease_free_status_consensus_number == 1) %>% 
        dplyr::count() %>% ungroup() %>% 
        select(CRIS_class, 'n')
      
      colnames(n_bad_prog_36) <- c('Class', 'recurred_before_36')
      
      # Build comparison
      
      annot_count <- n_samples %>% 
        full_join(n_msi_low, by = 'Class') %>%
        full_join(n_mss, by = 'Class') %>%
        full_join(n_msi_h, by = 'Class') %>%
        full_join(n_mucinous, by = 'Class') %>%
        full_join(n_kras_mut, by = 'Class') %>%
        full_join(n_braf_mut, by = 'Class') %>%
        full_join(n_nras_mut, by = 'Class') %>%
        full_join(n_tp53_mut, by = 'Class') %>%
        full_join(n_bad_prog_36, by = 'Class') %>%
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
      
      # Mucinous cardinality
      n_mucinous <- test_annot %>% group_by(CRIS_class, is_mucinous)  %>% dplyr::count() %>%
        filter(is.na(is_mucinous)) %>% ungroup() %>% select(CRIS_class, 'n')
      colnames(n_mucinous) <- c('Class', 'n_mucinous')
      
      # Micro satellite stability cardinality
      n_msi_low  <- test_annot %>% group_by(CRIS_class, msi_status)  %>% dplyr::count() %>%
        filter(is.na(msi_status)) %>% ungroup() %>% select(CRIS_class, 'n')
      colnames(n_msi_low) <- c('Class', 'n_msi_l')
      
      n_msi_h    <- n_msi_low
      colnames(n_msi_h) <- c('Class', 'n_msi_h')
      
      n_mss      <- n_msi_low
      colnames(n_mss) <- c('Class', 'n_mss')
      
      # KRAS mutated cardinality
      n_kras_mut <- test_annot %>% group_by(CRIS_class, KRAS_mutated)  %>% dplyr::count() %>%
        filter(is.na(KRAS_mutated)) %>% ungroup() %>% select(CRIS_class, 'n')
      colnames(n_kras_mut) <- c('Class', 'kras_mutated')
      
      # BRAF mutated cardinality
      n_braf_mut <- test_annot %>% group_by(CRIS_class, BRAF_mutated)  %>% dplyr::count() %>%
        filter(is.na(BRAF_mutated)) %>% ungroup() %>% select(CRIS_class, 'n')
      colnames(n_braf_mut) <- c('Class', 'braf_mutated')
      
      # NRAS mutated cardinality
      n_nras_mut <- test_annot %>% group_by(CRIS_class, NRAS_mutated)  %>% dplyr::count() %>%
        filter(is.na(NRAS_mutated)) %>% ungroup() %>% select(CRIS_class, 'n')
      colnames(n_nras_mut) <- c('Class', 'nras_mutated')
      
      # NRAS mutated cardinality
      n_tp53_mut <- test_annot %>% group_by(CRIS_class, TP53_mutated)  %>% dplyr::count() %>%
        filter(is.na(TP53_mutated)) %>% ungroup() %>% select(CRIS_class, 'n')
      colnames(n_tp53_mut) <- c('Class', 'tp53_mutated')
      
      # Samples with recurrence before 24 months
       n_bad_prog_36 <- test_annot %>% 
        group_by(CRIS_class)  %>% 
        filter(is.na(disease_free_months_consensus) | is.na(disease_free_status_consensus_number)) %>% 
        dplyr::count() %>% ungroup() %>% 
        select(CRIS_class, 'n')
      
      colnames(n_bad_prog_36) <- c('Class', 'recurred_before_36')
      
      annot_count <- n_samples %>% 
        full_join(n_msi_low, by = 'Class') %>%
        full_join(n_mss, by = 'Class') %>%
        full_join(n_msi_h, by = 'Class') %>%
        full_join(n_mucinous, by = 'Class') %>%
        full_join(n_kras_mut, by = 'Class') %>%
        full_join(n_braf_mut, by = 'Class') %>%
        full_join(n_nras_mut, by = 'Class') %>%
        full_join(n_tp53_mut, by = 'Class') %>%
        full_join(n_bad_prog_36, by = 'Class') %>%
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
        tests[[c]][['n_mucinous']]   <-  list()
        tests[[c]][['kras_mutated']] <-  list()
        tests[[c]][['braf_mutated']] <-  list()
        tests[[c]][['nras_mutated']] <-  list()
        tests[[c]][['tp53_mutated']] <-  list()
        tests[[c]][['n_msi_h']]      <-  list()
        tests[[c]][['recurred_before_36']]      <-  list()
        for (cl in count_annot[[c]]$Class){
          
          if (!grepl(cl, pattern = '^not CRIS-[A|B|C|D|E]$', fixed = FALSE)){
              tests[[c]][['n_mucinous']][[cl]]   <-  self$validate_bio(count_annot[[c]], cl, 'n_mucinous')
              
              tests[[c]][['kras_mutated']][[cl]] <-  self$validate_bio(count_annot[[c]], cl, 'kras_mutated')
              tests[[c]][['braf_mutated']][[cl]] <-  self$validate_bio(count_annot[[c]], cl, 'braf_mutated')
              tests[[c]][['nras_mutated']][[cl]] <-  self$validate_bio(count_annot[[c]], cl, 'nras_mutated')
              tests[[c]][['tp53_mutated']][[cl]] <-  self$validate_bio(count_annot[[c]], cl, 'tp53_mutated')
              tests[[c]][['n_msi_h']][[cl]]      <-  self$validate_bio(count_annot[[c]], cl, 'n_msi_h')
              tests[[c]][['recurred_before_36']][[cl]]  <-  self$validate_bio(count_annot[[c]], cl, 'recurred_before_36')
          }
          
        }
        
      }
      
      return(tests) 
    }

  )

)





