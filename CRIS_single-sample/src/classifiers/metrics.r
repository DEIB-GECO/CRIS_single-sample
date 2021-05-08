
# Description -------------------------------------------------------------

# Abstract superclass for single-label and multi-label metrics classes.
# It defines the mcc function for computing the mcc score given the 
# confusion matrix with respect to a specific class. Moreover, it defines
# the public methods that should be implemented in its subclasses.


# Class definition --------------------------------------------------------


Metrics <- R6Class(
  'Metrics',
  lock_objects = FALSE,
  lock_class = TRUE,
  
  ##### Private fields #####
  private = list(
    
    
    #' Binary MCC (computed for each label)
    #' 
    #' @param conf  The confusion matrix, defined as a list with tp, tn, fp, fn.
    #' (true positives, true negatives, false positives, false negatives).
    #' 
    #' @return the MCC score value
    #' @references https://arxiv.org/pdf/2008.05756.pdf
    mcc           = function(conf){
      
      num   <- conf$tp*conf$tn - conf$fp*conf$fn
      denom <- as.double((conf$tp + conf$fn))*
               as.double((conf$tp + conf$fp))*
               as.double((conf$tn + conf$fn))*
               as.double((conf$tn + conf$fp))
      
      if(is.na(denom) | is.na(num)){
        warning('overflow')
        return(NA)
      }
      
      if (min(num,denom) == 0)
        mcc <- 0
      else
        mcc <- num/sqrt(denom)
      

      return(mcc)
      
      
    }
  ),
  
  ##### Public fields #####
  # Here there are abstract methods to be implemented, whose documentation
  # is described in subclasses.
  
  public = list(
    
    confusion_matrix = function(ref,pred){},
    
    list_local_metrics = function(){},
    
    list_global_metrics = function(){},
    
    global_metrics = function(ref,pred){},
    
    #relaxed_metrics = function(ref,pred){},
    
    local_metrics = function(ref,pred){}
    
    #similarity_metrics = function(ref,pred){}
    
  )

)