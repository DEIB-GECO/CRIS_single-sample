
# Description -------------------------------------------------------------


# Exposes four methods to print errors, warnings, success and info messages; 
# The format_message method adds color to the text and the label the message. 
# At the end of the message, an end-of-line "\n" character is added.


# DEPENDENCIES ------------------------------------------------------------


# CONFIGURATION FLAGS -----------------------------------------------------


COLOUR_ACTIVE  <- TRUE                  

PRINT_SUCCESS  <- TRUE
PRINT_DEBUG    <- TRUE
PRINT_INFO     <- TRUE


# FUNCTIONS ---------------------------------------------------------------


#' Print "SUCCESS: ...." message
#' 
#' @param message      the message to print (character)
print_success <- function(message){
  
  if (PRINT_SUCCESS) {
    .format_message(paste("SUCCESS:", message, sep = " "),"SUCCESS")
  }
  
}




#' Print "DEBUG: ...." message
#' 
#' @param message      the message to print (character)
print_debug <- function(message){
  
  if (PRINT_DEBUG) {
    .format_message(paste("DEBUG:", message, sep = " "),"DEBUG")
  }
  
}





#' Print "INFO: ...." message
#' 
#' @param message      the message to print (character)
print_info <- function(message){
  
  if (PRINT_INFO) {
    .format_message(paste("INFO:", message, sep = " "),"INFO")
  }
  
}




#' Print time info message with the time difference required
#' 
#' @param message      the message to print (character)
print_time <- function(operation, start_time, end_time){
  
  print_info(paste("Time required for:", as.character(operation),":", end_time - start_time, "s"))
  
}


#' Print info with outcome of test
#' 
#' @param test_outcome   the test result (true = passed, false = not passed)
#' @param test_name      the name of the test, specifying the objective to be verified
print_test <- function(test_outcome, test_name){
  if (!exists('.N_TESTS'))
    .N_TESTS <<- 0
  
  if (!exists('.N_TESTS_FAILED'))
    .N_TESTS_FAILED <<- 0
  
  if (!exists('.N_TESTS_SUCCEDED'))
    .N_TESTS_SUCCEDED <<- 0
  
  .N_TESTS <<- .N_TESTS + 1
  if (!all(test_outcome == TRUE)){
    print_error(paste('Test failed. Test objective: ', test_name))
    .N_TESTS_FAILED   <<- .N_TESTS_FAILED + 1
  }else{
    print_success(paste('Test succeded. Test objective:', test_name))
    .N_TESTS_SUCCEDED <<- .N_TESTS_SUCCEDED + 1  
  }
}

# SUPPORT FUNCTIONS (PRIVATE) ---------------------------------------------




#' Add color and \n to a message, according to the type
#' 
#' @param message      the message to format (character)
#' @param message_type the type of message (character)
#' @return formatted message (character)
.format_message <- function(message, message_type){
  
  formatted_message <- paste(message,"\n", sep = "")
  
  if (COLOUR_ACTIVE) {
    return(switch(message_type,
         "SUCCESS" = cat(green(formatted_message)),
         "DEBUG" =   cat(silver(formatted_message)),
         "INFO" =    cat(white(formatted_message))
         )
    )
  }else{
    return(cat(formatted_message))
  }
  
}


