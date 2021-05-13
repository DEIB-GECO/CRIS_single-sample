
# DescriptionS ------------------------------------------------------------

# Exposes methods to load data from GCT, RDA, XLSX, CSV files.
# The load_file method finds the extension of the file and calls the appropriate
# method.

# Dependencies ------------------------------------------------------------

# Utilities
source(here('src', 'utils', 'source_utils.r'))


# FUNCTIONS ---------------------------------------------------------------


# Environment assignment
# https://stuff.mit.edu/afs/athena/software/r/current/lib/R/library/base/html/assignOps.html

#' Read a generic data file, calling the right function according to the type
#' 
#' @param path       The path of the file to read (single string)
#' @param separator  The separator for CSV files (single string)
#' @param colnames   Colnames parameter for reading from xlsx files (logical)
#' @param sheet      Names of sheets to read for xlsx files (character vector)
#' 
#' @return           The loaded data or NULL if something goes wrong
load_file <- function(path, separator = ",", colnames = TRUE, sheet = NULL) {
  
  # Check the input
  if (!all(check_type(path, 'character',1,1),fs::file_exists(path)))
    stop("load file: invalid or not existing path.")

  if (!check_type(colnames, "logical",1,1)) {
    warning("load file: provide valid colnames paramter. Defaulting to TRUE.")
    colnames <- TRUE
  }
  

  # Read according to extension of the file
  if (path_loader$has_extension(path, extension = ".gct")) {
    return(.load_from_gct(path))
    
  }else if (path_loader$has_extension(path, extension = ".csv")) {
    return(.load_from_csv(path, separator))
    
  }else if (path_loader$has_extension(path, extension = ".rda")) {
    return(.load_from_rda(path))
    
  }else if (path_loader$has_extension(path, extension = ".xlsx")) {
    if (is.null(sheet) | any(length(sheet) == 1))
      return(.load_from_excel(path, colnames, sheet))
    else
      return(.load_excel_sheets(path, sheet))
    
  }else if (path_loader$has_extension(path, extension = ".rds")) {
    return(.load_from_rds(path))
  }
  
  else{
    msg <- paste(
      c(
        "load file: The file format must be one of",
        path_loader$allowed_extensions
      ),
      collapse = " "
    )
    stop(msg)
    return(NULL)
  }
  
  
}




# SUPPORT FUNCTIONS (PRIVATE) ---------------------------------------------




#' Read a csv file
#' 
#' @param path    The path of the file to read
#' @return        The loaded data 
.load_from_csv <- function(path, separator) {
  

	data <- read.csv(file = path, sep = separator, check.names = FALSE);
	
	if (is.null(data)) {
	  stop(paste("The csv file", path, "could not be read.", sep = " "));
	}else {
		print_success(paste("The csv file", path, "has been read.", sep = " "));
	}

	return(data);
	
}




#' Read a gct file
#' 
#' @param path    The path of the file to read
#' @return        The loaded data 
.load_from_gct <- function(path){
	
	data <- phantasus::read.gct(path)
	
	if (is.null(data)) {
	  stop(paste("The gct file", path, "could not be read.", sep = " "));
	}else{
		print_success(paste("The gct file", path, "has been read.", sep = " "));
	}
	
	return(data);
	
}




#' Read an rda file
#' 
#' @param path    The path of the file to read
#' @return        The loaded data 
.load_from_rda <- function(path) {
	
	data <- load(path, envir = globalenv())
	
	if (is.null(data)) {
	  stop(paste("The rda file", path, "could not be read.", sep = " "));
	}else{
		print_success(paste("The rda file", path, "has been read.", sep = " "));
	}
	
	return(data);
}




#' Read an rds file
#' 
#' @param path    The path of the file to read
#' @return        The loaded data 
.load_from_rds <- function(path) {
  
  data <- readRDS(path)
  
  if (is.null(data)) {
    stop(paste("The rds file", path, "could not be read.", sep = " "));
  }else{
    print_success(paste("The rds file", path, "has been read.", sep = " "));
  }
  
  return(data);
}





#' Read an excel file
#' 
#' @param path      The path of the file to read
#' @param colnames  flag that specifies if colnames are included (TRUE/FALSE)
#' @return          A dataframe with loaded data 
.load_from_excel <- function(path, colnames, sheet_name = NULL){
  
  if (check_type(sheet_name, "character",1,1)){
    print_info(paste("Sheet name:", sheet_name))
    data <- read_excel(path, col_names = colnames, sheet = sheet_name)
  }else{
    data <- read_excel(path, col_names = colnames)
  }

  if (is.null(data)) {
    stop(paste("The excel file", path, "could not be read.", sep = " "));
  }else{
    print_success(paste("The excel file", path, "has been read.", sep = " "));
  }
  
  return(as.data.frame(data));
  
}





#' Read an excel file with a several sheets in it
#' 
#' @param path    The path of the file to read
#' @param sheets_list  The list of sheets to read
#' @return A list with loaded sheets (dataframes)
.load_excel_sheets <- function(path, sheets_list, colnames = TRUE){
  
  # Check colnames parameter validity
  if (length(colnames) == 1) {
    use_colnames <- replicate(length(sheets_list), colnames)
  }
  else if (length(colnames) != length(sheets_list)) {
    warning("Colnames length does not match sheets length. Using TRUE.")
    use_colnames <- vector("logical", length(sheets_list))
  }
  
  # Reading the excel sheets  
  data <- list()
  
  for (i in seq(sheets_list)) {
    sheet <- sheets_list[i]
    data[[sheet]] <- read_excel(path, sheet = sheet, col_names = use_colnames[i])
    if (is.null(data[[sheet]])){
      stop(paste("The sheet", sheet, "could not be read.", sep = " "));
    }
  }
  return(data)
}




