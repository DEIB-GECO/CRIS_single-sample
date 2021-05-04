# Description -------------------------------------------------------------

# Superclass of all filters

# Dependencies ------------------------------------------------------------

library(R6)

# Class definition --------------------------------------------------------

GenesFilter  <- R6Class('GenesFilter',
    inherit      = Filter,
    lock_objects = FALSE,
    lock_class   = TRUE,
    
    public = list(
      
      
      #' Select only the samples that appear in the given genes list.
      #' 
      #' @param expr_matrix   expression matrix to be filtered. Genes
      #'                      are row names
      #' @param list          list of genes that can be selected (character vector)
      #' @return expression matrix with reduced number of rows, corresponding
      #'         to the selected genes
      filter_by_list = function(expr_matrix, list){
        
        # Check parameters
        super$filter_by_list(expr_matrix, list)
        
        genes_list <- list %>% unlist() %>% as.character()
        
        # Find genes of matrix that appear in the list
        dataset_genes  <- data.frame(gene = rownames(expr_matrix))
        selected_genes <- dataset_genes %>% filter(gene %in% genes_list) %>%
          unlist() %>% as.character()
        
        # Apply filter on the matrix
        expr_matrix <- expr_matrix[sort(selected_genes), ]
        
        # Print result info
        print_info(paste("Selected genes:", nrow(expr_matrix)))
        
        # Return expression of requested genes
        return(expr_matrix)
        
      },
      
      
      
      
      
      #' Keep only the genes for which there is a consistent portion of samples (above
      #' a minimum percentage) that has an expression above the threshold.
      #' 
      #' @param expr_matrix  the expression matrix to be filtered
      #' @param expr_thr     the threshold to be used for expression (numeric)
      #' @param samples_thr  the threshold to be used for samples percentage (numeric in [0;1])
      #' @param features     an optional list of genes that must not be removed
      #' 
      #' @return A subset of the initial expression matrix
      filter_by_threshold = function(expr_matrix, expr_thr, samples_thr, features = NULL){
        
        print_info("Filter by expression threshold:")
        # Check input
        stopifnot(check_type(expr_matrix,'matrix',c(1,1)), check_type(expr_matrix,'numeric'))
        stopifnot(check_type(expr_thr, "numeric",1,1), check_type(samples_thr, "numeric",1,1))
        
        # Check percentage is correct
        stopifnot(samples_thr >= 0, samples_thr <= 1)
        
        # Build the filter
        
        gene_selected <- vector("logical", nrow(expr_matrix))
        n_samples     <- ncol(expr_matrix)
        
        for (i in seq(nrow(expr_matrix))) {
          
          gene <- rownames(expr_matrix)[i]
          
          # Percentage of samples with current gene's expression > expr_thr
          high_expr_samples    <- expr_matrix[gene, expr_matrix[gene,] > expr_thr]
          percentage_high_expr <- length(high_expr_samples)/n_samples
          
          # Keep current gene only if percentage of samples >= samples_thr
          gene_selected[i] <- percentage_high_expr >= samples_thr
          
          
          # If requested, keep the genes that belong to the template (even if
          # they should have been discarded by the above filter). 
          if (!is.null(features) & any(features[,1] == gene))
            gene_selected[i] <- TRUE
          
        }
        
        # Apply filter
        filtered_expr_matrix <- expr_matrix[gene_selected,]
        
        # Print status information
        n_selected_genes <- nrow(filtered_expr_matrix)
        ntot_genes       <- nrow(expr_matrix)
        
        print_info("Selected genes:")
        print_info(paste(n_selected_genes, 
                         "over", 
                         ntot_genes, 
                         "genes", sep = " ")
        ) 
        return(filtered_expr_matrix)
        
      }
    )
                          
)
