library(survival)
library(R6)

BiologicalPlots <- R6Class(
  'BiologicalPlots',
  lock_objects = TRUE,
  lock_class = TRUE,
  
  ##### Private fields #####
  
  private = list(
    .xlim = 36,
    .time = NULL,
    .status = NULL,
    .annotations = NULL,
    .groups      = c('is_primary','is_secondary','is_other'),
    .prepare_km_annotations = function(){
      
      # Check the time and status attribute appear in the annotations
      if (!all(c(private$.time, private$.status) %in% colnames(private$.annotations)))
          stop('BiologicalPlots: cannot find time and status attribute in annotaions for Kaplan Meier data')
      
      # Extract aliquot ID, time and status attribute from annotations and rename them
      km_annot <- private$.annotations[!is.na(private$.annotations[,private$.time]) &
                                       !is.na(private$.annotations[,private$.status]), ]
      km_annot <- km_annot[,c(ALIQUOT_LABEL, private$.time, private$.status)]
      colnames(km_annot) <- c(ALIQUOT_LABEL, 'df_time','df_status')
      

      # Adjust the data: samples with recurrence (status = 1) after xlim months must
      # be ignored and thus set as disease-free (status = 0)
      df_samples <- km_annot %>% 
        filter(df_time >= private$.xlim) %>% 
        select(aliquot_id) %>%
        unlist()
      
      km_annot[km_annot$aliquot_id %in% df_samples, 'df_status'] <- 0
      km_annot[km_annot$aliquot_id %in% df_samples, 'df_time']   <- private$.xlim

      # Return the prepared annotations for KM plot
      return(km_annot)
      
    },
    
    .km_data_sl = function(cl_ref, cl){
      
      # Set as primary/ref if assigned class is equal to the provided one (cl)
      cl_ref[cl_ref[,CLASS_LABEL] == cl, 'is_primary'] <- TRUE
      cl_ref[cl_ref[,CLASS_LABEL] == cl, 'ref'] <- TRUE
      
      # Set as primary if assigned class is different from the provided one (cl)
      cl_ref[cl_ref[,CLASS_LABEL] != cl, 'is_other'] <- TRUE
      
      # The is_secondary group is never set as the predicted class is always 
      # primary in single-label classifiers
      cl_ref[, 'is_secondary'] <- FALSE
      
      # Assign class label basing to primary
      for (k in seq(nrow(cl_ref))){
        
        if (cl_ref[k,'is_primary'] == TRUE){
          cl_ref[k,'CRIS_class'] <- paste(cl, 'primary')
        }else if (cl_ref[k,'is_secondary'] == TRUE){
          cl_ref[k,'CRIS_class'] <- paste(cl, 'not primary')
        }else if (cl_ref[k,'is_other'] == TRUE){
          cl_ref[k,'CRIS_class'] <- paste('not', cl)
        }else{
          stop('Missing group')
        }
          
      }
      
      return(cl_ref)
    },
    
    .km_data_ml_alg_adapt = function(cl_ref, cl, required_cols){
      
      # Set ref on samples assigned to provided class (cl)
      cl_ref[cl_ref[,cl] == 1, 'ref'] <- TRUE

      # Set as primary if the class with highest score (class_label) is equal to the provided one (cl)
      cl_ref[cl_ref[,CLASS_LABEL] == cl, 'is_primary'] <- TRUE
      
      # Set as secondary on samples assigned to cl but the class with highest score is different
      cl_ref[cl_ref[,CLASS_LABEL] != cl & cl_ref[,cl] == 1, 'is_secondary'] <- TRUE

      # Set as other if the provided class is never assigned
      cl_ref[cl_ref[,cl] == 0, 'is_other'] <- TRUE
      
      # Assign class label basing on groups
      for (k in seq(nrow(cl_ref))){
        
        if (cl_ref[k,'is_primary'] == TRUE){
          cl_ref[k,'CRIS_class'] <- paste(cl, 'primary')
        }else if (cl_ref[k,'is_secondary'] == TRUE){
          cl_ref[k,'CRIS_class'] <- paste(cl, 'not primary')
        }else if (cl_ref[k,'is_other'] == TRUE){
          cl_ref[k,'CRIS_class'] <- paste('not', cl)
        }else{
          stop('Missing group')
        }
          
      }
      
      return(cl_ref)
    },
    
    .km_data_ml = function(cl_ref, cl){
      
      # Set ref on samples assigned to provided class (cl)
      cl_ref[cl_ref[,cl] == 1, 'ref'] <- TRUE
  
      # Fill groups
      for (k in seq(nrow(cl_ref))) {
  
        # Get assigned classes
        assigned_classes <- cl_ref[k,c(ALIQUOT_LABEL, CRIS_CLASSES)]
        assigned_classes <- colnames(assigned_classes)[which(assigned_classes[1,] == 1)]
  
        # No class has been assigned >> no primary class
        if (length(assigned_classes) == 0){
          primary_class <- ''
        }
        # One class only has been assigned >> it is also primary
        else if (length(assigned_classes) == 1){
          primary_class <- assigned_classes[1]
        }
        # At least two classes assigned >> primary class is the one with minimum rank
        else{
          ranks <- cl_ref[k, paste('rank',assigned_classes, sep = '')]
          primary_class <- which.min(ranks) %>%                      
                           names() %>%                               
                           gsub(pattern = 'rank', replacement = '')
        }
        
        # Assign groups and class label (in words)
        
        # Required class cl has not been assigned
        if (cl_ref[k,cl] == 0){
          cl_ref[k,'is_other'] <- TRUE
          cl_ref[k,'CRIS_class'] <- paste('not', cl)
        }
        # Required class cl has been assigned as primary
        else if (primary_class == cl & cl_ref[k,cl] == 1){
          cl_ref[k,'is_primary'] <- TRUE
          cl_ref[k,'CRIS_class'] <- paste(cl, 'primary')
        }
        # Required class cl has been assigned but is not primary
        else if (primary_class != cl & cl_ref[k,cl] == 1){
          cl_ref[k,'is_secondary'] <- TRUE
          cl_ref[k,'CRIS_class'] <- paste(cl, 'not primary')
        }
        # Error
        else {
          stop('Error in group assignment of ML')
        }
  
      }
  
      return(cl_ref)
      
    },
      
    .bio_plot_ref = function(ref, type,cl, annot){
      
      # Check the type parameter in input
      if (!check_type(type, 'character', length_min = 1, length_max = 1) |
          !type %in% c('ntp_sl','sl', 'ntp_ml','aa_ml','pt_ml'))
        stop('BiologicalPlots: type parameter for kaplan meier ref is invalid.')
        
      # Check the class
      if (!check_type(cl, 'character', length_min = 1, length_max = 1) |
          !cl %in% CRIS_CLASSES)
        stop('BiologicalPlots: cl parameter for kaplan meier ref is invalid.')
      
      # Set required columns in ref depending on the type
      if (any(type %in% c('ntp_sl','sl'))){
        required_cols <- c(ALIQUOT_LABEL, CLASS_LABEL)
      }else if (any(type %in% c('ntp_ml','aa_ml'))){
        required_cols <- c(ALIQUOT_LABEL, CLASS_LABEL, cl)
      }else{
        ranks <- paste('rank', CRIS_CLASSES, sep = '')
        required_cols <- c(ALIQUOT_LABEL, CRIS_CLASSES, ranks) 
      }
     
      # Check the reference
      if (!check_type(ref, 'data.frame', dim_min = c(1,1))|
          !all(required_cols %in% colnames(ref)))
        stop('BiologicalPlots: ref parameter for kaplan meier ref is invalid (missing columns)')
      
      # Prepare the reference specific for the required class
      cl_ref <- ref %>% 
        mutate(ref = FALSE) %>%           # Belongs to CRIS-X (primary/not primary)
        mutate(is_primary = FALSE) %>%    # Assigned to CRIS-X as primary 
        mutate(is_secondary = FALSE) %>%  # Assigned to CRIS-X as not primary 
        mutate(is_other = FALSE) %>%      # Not assigned to CRIS-X 
        mutate(CRIS_class = 'none') %>%   # Final class label (in words)
        as.data.frame()
      
      # Call the right function depending on the type
      if (any(type %in% c('ntp_sl','sl'))){
        cl_ref <- private$.km_data_sl(cl_ref, cl)
      }else if (any(type %in% c('ntp_ml','aa_ml'))){
        cl_ref <- private$.km_data_ml_alg_adapt(cl_ref, cl)
      }else {
        cl_ref <- private$.km_data_ml(cl_ref, cl)
      } 
      
      # Exit if two or more groups are assigned contemporary on some samples
      # (each sample can be either primary, secondary or not assigned)
      if (any(rowSums(cl_ref[,private$.groups]) > 1))
        stop(paste('BiologicalPlots (KM_SL): Inconsistent groups:', private$.groups))
      
      # Select only columns that are strictly needed
      cl_ref <- cl_ref[ ,c(ALIQUOT_LABEL, 'CRIS_class', 'ref', private$.groups)]
     
      # Add annotations and return
      cl_ref <- inner_join(cl_ref, annot, by = ALIQUOT_LABEL)
      
      return(cl_ref)
      
    }
    
  ),
  
  ##### Public fields #####
  
  public = list(
    
    initialize = function(annotations,
                          xlim = 60,
                          time   = 'disease_free_months_consensus',
                          status = 'disease_free_status_consensus_number') {
      
      # Set annotations
      if (check_type(annotations, type = 'data.frame', dim_min = c(1,1)) &
          ALIQUOT_LABEL %in% colnames(annotations))
        private$.annotations <- annotations
      else
        warning(paste('BiologicalPlots: invalid annotations. Setting default value: ', 
                      private$.annotations))
      
      # Set xlim for plot (time window in months)
      if (check_type(xlim, type = 'numeric', length_min = 1, length_max = 1) &
          xlim > 0)
        private$.xlim <- xlim
      else
        warning(paste('BiologicalPlots: invalid xlim. Setting default value: ', 
                      private$.xlim))
      
      # Set time attribute 
      if (check_type(time, type = 'character', length_min = 1, length_max = 1))
        private$.time <- time
      else
        warning(paste('BiologicalPlots: invalid time. Setting default value: ', 
                      private$.time))
      
      # Set status attribute 
      if (check_type(status, type = 'character', length_min = 1, length_max = 1))
        private$.status <- status
      else
        warning(paste('BiologicalPlots: invalid status. Setting default value: ', 
                      private$.status))
      
    },
    
    ref_forest_plots = function(ref, type,cl){
      
      # Create the reference for kaplan meier using the entire annotation data
      return(private$.bio_plot_ref(ref, type, cl, private$.annotations))
      
    },
    
    ref_kaplan_meier = function(ref, type, cl){
      
      # Prepare the annotation for Kaplan Meier
      km_annot <- private$.prepare_km_annotations()

      # Create the reference for kaplan meier using the above annotations
      return(private$.bio_plot_ref(ref, type, cl, km_annot))
      
    },
    
    #' @references issue: https://github.com/kassambara/survminer/issues/443
    compute_kaplan_meier = function(cl_ref, group, cl){
      
      # Check group input
      if (!check_type(group, 'character',1,1) |
          !group %in% c('ref','is_primary'))
        stop('use a single group among ref and is_primary')
      
      # Formula for the Kaplan-meier
      f_req_fit <- as.formula(paste('Surv(df_time, df_status) ~', group))
  
      # Get samples that belong to the current group and that do not belong to the class only.
      cl_ref <- cl_ref[cl_ref[,group] | cl_ref[,'is_other'], ]  
      
      # In case of empty data, exit
      if (nrow(cl_ref) < 1){
        warning(paste('BiologicalPlots: empty data for kaplan meier -', group, ', ', cl))
        return(list(
          data = cl_ref,
          group = group,
          fit  = NA,
          cox  = NA,
          plot_data = NA
        ))
      }
      
      
      # If everything is ok, use the aliquot ID to set the rownames
      rownames(cl_ref) <- NULL  
      cl_ref <- cl_ref %>% 
          column_to_rownames(ALIQUOT_LABEL) 
      
      # Fit the survival curve
      req_fit <- survfit(f_req_fit, data = cl_ref)
      
      # Compute COX model and extract p-value of log rank test
      req_cox <- tryCatch(coxph(f_req_fit, data = cl_ref), error = function(e){
        stop(paste('BiologicalPlots: error in cox model for kaplan meier -', group, ', ', cl))
      })
      
      pvalue  <- summary(req_cox)$sctest['pvalue']
  
      # Define the plot data 
  
      # Classes for the plot legend
      if (group == 'ref')
        classes <- c(paste('not', cl), cl)
      else
        classes <- c(paste('not', cl), paste(cl, 'primary'))
      
      # Title of the plot
      plot_title_pvalue <- paste('(p-value = ', signif(pvalue,3), ')',sep = '')
      
      # Put together all plot data
      plot_data  <- list(
        xlim = private$.xlim,
        ylim = c(0,1),
        xlab = 'disease free months',
        ylab = 'disease free probability',
        classes = classes,
        title = plot_title_pvalue
      )
      
      # Return result (data, survfit, cox model, plot data)
      return(list(
        data = cl_ref,
        group = group,
        fit  = req_fit,
        cox  = req_cox,
        plot_data = plot_data)
      )
      
    },
    
    show_km_plot = function(km_res,alg_name){
  
      # Select colour for CRIS-X as the fourth tone in OrRd palette (orange tones)
      cris_col <- brewer.pal(n=4,name="OrRd")[4]
      
      # Select colour for not CRIS-X as the fourth tone in GnBu palette (green tones)
      not_cris_col <- brewer.pal(n=4,name="GnBu")[4]
    
      # Create palette with the two colours
      cols <- c(not_cris_col,cris_col)
      
      # Create the plot
      plot(km_res$fit,                                # survfit object
           col = cols,                                # curve colours
           mark.time = TRUE,                          # show ticks for censored data
           xlab = km_res$plot_data$xlab,              # label x axis
           ylab = km_res$plot_data$ylab,              # label y axis
           cex.lab  = 1.4,                            # expansion of axes labels text
           cex.axis = 1.4                             # expansion of tick labels text
           )
      
      # Add the legend
      legend("bottomleft",                             # legend position
             legend = km_res$plot_data$classes,        # legend names
             col = cols,                               # legend colours
             pch = 8,                                  # size of coloured symbol in the legend
             cex = 1.4)                                # text size

      # Add the title
      title(main = paste(alg_name, km_res$plot_data$title), # title text
            cex.main = 1.4)                                 # text size
     
    },
    
    
    # Prepare data for the forest plot
    single_forest_plot = function(stats, test_type, test_attribute, test_number, classes = CRIS_CLASSES){
      test_out <- data.frame()
      
      for (class in classes){
        cl_listed   <- names(stats[[class]][[test_type]])
        test_out_cl <- data.frame()
        for (c in cl_listed) {
          if (!grepl(x = c, pattern = 'not primary',fixed = TRUE)){
            
            # Build vector with test name
            test_name <- c(paste(c,test_attribute), NA, NA, NA, NA, NA)
            
            # Get required data
            test_data <- stats[[class]][[test_type]][[c]]$test_out
            Description <- c('Samples assigned to class','P-value','Effect size','Odds ratio','CI low','CI high')
            samples     <- sum(stats[[class]][[test_type]][[c]]$table[c, ])
            p.value     <- test_data$p.value       %>% signif(digits = 3)
            eff.size    <- test_data$eff.size      %>% round(digits = 3)
            odds.ratio  <- test_data$odds.ratio  %>% round(digits = 3)
            conflow     <- test_data$conf.int.low  %>% round(digits = 3)
            confhigh    <- test_data$conf.int.high %>% round(digits = 3)
            test_values <- c(samples, p.value, eff.size, odds.ratio, conflow, confhigh)
            
            # Create data.frame with result
            res <- data.frame(Test = test_name, Description, X = test_values, Y = test_number)
            
            # Add the result to global data.frame of the required class
            if (check_type(test_out_cl, 'data.frame', dim_min = c(0,0), dim_max = c(0,0)))
              test_out_cl <- res
            else
              test_out_cl <- rbind(test_out_cl, res)
            
            # Increment test number
            test_number <- test_number + 1
          }
          
        }
        # Add class results to the method result
        if (check_type(test_out, 'data.frame', dim_min = c(0,0), dim_max = c(0,0)))
          test_out <- test_out_cl
        else
          test_out <- rbind(test_out, test_out_cl)
        
        
      }
      
      return(test_out)
      
    },

    compute_all_forest_plots = function(forest_attributes, label_attributes, classes, tests, test_number_init = NULL){
      
      all_forest_plots <- data.frame()
      if (is.null(test_number_init))
        test_number      <- 1
      else
        test_number <- test_number_init
      for (i in seq(length(forest_attributes))){
        test_attr  <- forest_attributes[[i]]
        label_attr <- label_attributes[[i]]
        
        # Forest plot for SL 
        fplot_sl  <- self$single_forest_plot(stats = tests, 
                                             test_type      = test_attr, 
                                             test_attribute = label_attr,
                                             test_number    = test_number,
                                             classes = classes[[test_attr]])
      
        test_number <- max(fplot_sl$Y) + 1
        
        if (check_type(all_forest_plots, 'data.frame', dim_min = c(0,0), dim_max = c(0,0)))
          all_forest_plots <- fplot_sl
        else
          all_forest_plots <- rbind(all_forest_plots, fplot_sl)
        
      }
      
      return(all_forest_plots)
    }
  ),
  
  ##### Active fields #####
  
  active = list(
    
    xlim = function(v){
      if (missing(v))
        return(private$.xlim)
      else
        stop('KaplanMeier: xlim is read-only')
    },
    
    time = function(v){
      if (missing(v))
        return(private$.time)
      else
        stop('KaplanMeier: time is read-only')
    },
    
    status = function(v){
      if (missing(v))
        return(private$.status)
      else
        stop('KaplanMeier: status is read-only')
    },
    
    annotations = function(v){
      if (missing(v))
        return(private$.annotations)
      else
        stop('KaplanMeier: annotations is read-only')
    }
    
    
  )
   
  ##### End of class ####
)

