
# Description -------------------------------------------------------------

# Generic class to handle the biological validation of classifiers


# Constants definitions ----------------------------------------------------

.GROUPS <- c('CRIS_class','ref','is_primary','is_secondary','is_other')
      
# Class definition --------------------------------------------------------
Annotation <- R6Class(
  
  'Annotation',
  lock_objects = TRUE,
  lock_class   = TRUE,
  
  private = list(
    
    .is_invalid_number = function(number) {

      not_numeric <- !check_type(number, 'numeric', 1, 1)
      is_missing  <- is.null(number) | is.na(number)
      is_invalid  <- is.nan(number)  | is.infinite(number)
      
      return(any(is_missing | is_invalid | not_numeric))
      
    }
    
  ),
  
  public = list(
 

# Measures ----------------------------------------------------------------


    effect_size = function(contingency, class, test_type){
      
      contingency <- as.matrix(contingency)
      effect_size <- 0
      
      # Numerators (samples annotated for required test in class or in whole dataset)
      class_succ  <- contingency[class, test_type]
      total_succ  <- sum(contingency[ , test_type])
      
      # Denominators (samples in class or in whole dataset)
      class_samps <- sum(contingency[class, ])
      total_samps <- sum(contingency)
      
      # Class frequency 
      if (!private$.is_invalid_number(class_succ) & !private$.is_invalid_number(class_samps)) {
        class_ratio <- class_succ / class_samps
      }else {
        class_ratio <- 0
      }
      # Expected frequency
      if (!private$.is_invalid_number(total_succ) & !private$.is_invalid_number(total_samps)) {
        expected_ratio <- total_succ / total_samps
      }else {
        expected_ratio <- 0
      }
      
      # Compute effect size
      if (!private$.is_invalid_number(expected_ratio)){
        effect_size <- class_ratio/expected_ratio
      }else{
        # Uncomputable value
        effect_size <- -1
      }
      
      return(effect_size)
    },


# Common ------------------------------------------------------------------

    count_annotations = function(test_annot){
  
      stop('To be implemented')
    },
    
    count_not_annotated = function(test_annot){
  
      stop('To be implemented')
    },
    
    merge_annot_counts = function(annot_count, not_annotated){

      colnames(not_annotated)[-(1:2)] <- paste('na', colnames(not_annotated)[-(1:2)], sep = '_')
      
      annot_count <- full_join(annot_count, not_annotated, by = c('Class','n_samples'))
      
      return(annot_count)
    },
    
    add_annot_proportions = function(annot_count, not_annotated){
      
      cls <- annot_count$Class
      not_annotated <- not_annotated %>% column_to_rownames('Class')
      annot_count <- annot_count %>% column_to_rownames('Class')
      
      # Count proportions of above quantities on each class
      for (c in colnames(annot_count)[-1]){
        prop <- annot_count[cls,c]/(not_annotated[cls,'n_samples'] - not_annotated[cls, c]) %>% round(4)
        prop[is.nan(prop)] <- 0
        annot_count <- cbind(annot_count, prop)
        colnames(annot_count)[ncol(annot_count)] <- paste('prop',c, sep = '_')
      }
      
      annot_count <- annot_count %>% rownames_to_column('Class')
      return(annot_count)
      
    },
    
    get_test_value = function(test, type, test_type, class){
      test <- test[[test_type]][[class]]$test_out
      
      value <- switch(type,
                      'p_value'     = test$p.value,
                      'odds_ratio'  = test$odds.ratio,
                      'effect_size' = test$eff.size)
      
      if (private$.is_invalid_number(value)){
        return(NA)
      }else{
        return(value %>% format(scientific = TRUE))
      }
    },

# SL annotations ----------------------------------------------------------

    prepare_annotations_sl = function(res, all_annot, type){
      
      if ('CRIS_class' %in% colnames(all_annot)){
        all_annot <- all_annot[,-(which(colnames(all_annot) == 'CRIS_class'))]
      }
      
      test_annot <- list()
      biopl <- BiologicalPlots$new(annotations = all_annot)
      
      for (c in levels(F_CRIS_CLASSES)){
        if (type == 'ntp_sl'){
          test_annot[[c]] <- biopl$ref_forest_plots(res$result, type = 'ntp_sl', cl = c)
        }else if (type == 'sl'){
          test_annot[[c]] <- biopl$ref_forest_plots(
            res$binary_res %>% rownames_to_column(ALIQUOT_LABEL),
            type = 'sl',
            cl = c
          )
        }
      }
      return(test_annot)
    },


    get_count_annot_sl = function(res, annot, type){
      
      test_annot  <- self$prepare_annotations_sl(res, annot, type)
      
      count_annot <- list()
      
      for (c in levels(F_CRIS_CLASSES)){
        if (is.null(test_annot)){
          count_annot[[c]] <- NULL
        }else{
          
          not_annotated <- self$count_not_annotated(test_annot[[c]])
          
          count_annot[[c]] <- self$count_annotations(test_annot[[c]]) %>%
            self$add_annot_proportions(not_annotated) %>%
            self$merge_annot_counts(not_annotated)
          
          
          count_annot[[c]][,'Class'] <- gsub(x = count_annot[[c]][,'Class'], pattern = ' primary', replacement = '')  
        }
      }
      
      
      return(count_annot)
    },


# ML annotations ----------------------------------------------------------


    prepare_annotations_ml = function(ref_data, res, all_annot, type){
        

        if ('CRIS_class' %in% colnames(all_annot)){
          all_annot <- all_annot[,-(which(colnames(all_annot) == 'CRIS_class'))]
        }
        # Read ml result and select assigned class
     
        res <- res$binary_res %>% rownames_to_column(ALIQUOT_LABEL)
        
        test_annot <- list()
        biopl <- BiologicalPlots$new(annotations = all_annot)
        if (nrow(res) == 0)
          return(NULL)
        
        for (c in levels(F_CRIS_CLASSES)){
          test_annot[[c]] <- biopl$ref_forest_plots(res, type, cl = c)
          
        }
        
        return(test_annot)
    },
    
    prepare_annotations_ntp_ml = function(ref_data, all_annot){
      
      if ('CRIS_class' %in% colnames(all_annot)){
          all_annot <- all_annot[,-(which(colnames(all_annot) == 'CRIS_class'))]
      } 
      if (nrow(ref_data) == 0)
          return(NULL)
      
      test_annot <- list()
      biopl <- BiologicalPlots$new(annotations = all_annot)
      for (c in levels(F_CRIS_CLASSES)){
          
        ref <-  biopl$ref_forest_plots(ref_data, type = 'ntp_ml', cl = c)
        
        test_annot[[c]] <- inner_join(ref[,c(ALIQUOT_LABEL, .GROUPS)], all_annot, by = ALIQUOT_LABEL) %>% 
          arrange(CRIS_class, aliquot_id)
      }
      
      return(test_annot)
    },
    
    get_count_annot_ml = function(res, annot, ref_data, type){
  
      if (check_type(res, 'null')){
        test_annot <- self$prepare_annotations_ntp_ml(ref_data, annot)
      }else{
        test_annot <- self$prepare_annotations_ml(ref_data, res, annot, type)
      }
      
      count_annot <- NULL
      
      if (!is.null(test_annot)) {
        
        count_annot <- list()
        for (c in names(test_annot)){
          
          if (is.null(test_annot[[c]])){
            count_annot[[c]] <- NULL
          }else{
            not_annotated <- self$count_not_annotated(test_annot[[c]])
            # Summ all classes counts
            cum_classes   <- not_annotated %>% filter(Class != paste('not', c)) %>% column_to_rownames('Class')
            
            # Add cumulative count
            not_annotated <- rbind(not_annotated, as.numeric(c(NA, colSums(cum_classes))))
            not_annotated[nrow(not_annotated), 'Class'] <- c
            
            count_annot[[c]] <- self$count_annotations(test_annot[[c]])
            cum_classes   <- count_annot[[c]] %>% filter(Class != paste('not', c)) %>% column_to_rownames('Class')
            count_annot[[c]] <- rbind(count_annot[[c]], as.numeric(c(NA, colSums(cum_classes))))
            count_annot[[c]][nrow(count_annot[[c]]), 'Class'] <- c
            
            count_annot[[c]] <- count_annot[[c]] %>%
              self$add_annot_proportions(not_annotated) %>%
              self$merge_annot_counts(not_annotated)
          }
        }
        
      }
      
      return(count_annot)
    },
    

# Tests -------------------------------------------------------------------

    #' @references # https://towardsdatascience.com/fishers-exact-test-in-r-independence-test-for-a-small-sample-56965db48e87
    validate_bio = function(count_annot, class, attr){
  
      na_attr <- paste('na', attr, sep = '_')
      
      # Build contingency table with required attribute and class
      contingency <- count_annot[,c('Class','n_samples',na_attr, attr)] %>%
        column_to_rownames('Class')
      
      # Find not annotated samples
      n_def <- contingency$n_samples - contingency[,na_attr] - contingency[,attr]
      contingency <- cbind(contingency, n_def) %>%  select(attr,n_def)
      
      # Get the data for the row 'not CRIS-X'
      control_group <- count_annot[grepl(count_annot[,'Class',], 
                                     pattern = '^not CRIS-[A|B|C|D|E]$', 
                                     fixed = FALSE), 
                                  'Class']
      
      # Contingency Table
      contingency <- contingency[c(class,control_group), ]
      test <- NULL
      type <- NULL
      test_out <- list()
      
      # Perform Fisher Test
      if (TRUE){
        test <- fisher.test(contingency)
        type <- 'Fisher'
        test_out <- list(
          odds.ratio    = test$estimate %>% as.numeric(),
          conf.int.low  = test$conf.int[1],
          conf.int.high = test$conf.int[2],
          p.value       = test$p.value,
          eff.size      = self$effect_size(contingency, class, attr)
        )
      }# Never executed. Change the If condition to execute this too.
      else{ 
        test <- chisq.test(contingency)
        type <- 'Chisq'
        test_out <- list(
          odds.ratio    = 1,
          conf.int.low  = 1,
          conf.int.high = 1,
          p.value       = test$p.value,
          eff.size      = self$effect_size(contingency, class, attr)
        )
      }
      
      
      if (test$p.value <= 0.05){
        print_success(paste(test$method, ' (', attr, ', ', class, '). p-value:',test$p.value , sep = ''))
      }else{
        print_info(paste(test$method, ' (', attr, ', ', class, '). p-value:',test$p.value , sep = ''))
      }
      
      # Return result
      return(list(
        table = contingency,
        test  = test,
        test_out = test_out,
        type  = type
      ))
    
    },
    
    
    stat_tests_class_comparison = function(class, test_ntp_sl, test_ntp_ml, test_sl, test_ml, test_type, value = 'p_value'){
      
      ntp_sl  <- c(algorithm = 'ntp_sl', type = 'ref')
      ntp_ml  <- c(algorithm = 'ntp_ml', type = 'ref')
      sl   <- list()
      ml   <- list()
      cl_types   <- c(paste(class, 'primary'),
                      paste(class, 'not primary'))
      
      sl_methods <- names(test_sl)
      ml_methods <- names(test_ml)
      
      for (c in cl_types){
        
        ntp_sl[c] <- self$get_test_value(test_ntp_sl[[class]], value, test_type, c)
        ntp_ml[c] <- self$get_test_value(test_ntp_ml[[class]], value, test_type, c)
        
        sl[[c]] <- list()
        for (m in sl_methods){
            sl[[c]][[m]] <- self$get_test_value(test_sl[[m]][[class]], value, test_type, c)
        }
        
        sl[[c]] <- unlist(sl[[c]])
        
        ml[[c]] <- list()
        for (m in ml_methods){
          ml[[c]][[m]] <- self$get_test_value(test_ml[[m]][[class]], value, test_type, c)
        }
        
        ml[[c]] <- unlist(ml[[c]])
        
      }  
      
      sl <- cbind(type = 'sl', as.data.frame(sl, check.names = FALSE)) %>% rownames_to_column('algorithm')
      ml <- cbind(type = 'ml', as.data.frame(ml, check.names = FALSE)) %>% rownames_to_column('algorithm')
      
      comparison <- rbind(ntp_sl, ntp_ml, sl, ml)  
    
      return(comparison)
    }

  )
)

