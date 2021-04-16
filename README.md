# CRIS_single-sample

The work presented in this repository allows to apply single-sample classifiers for CRIS on colorectal cancer data. The work has been organized in the CRIS_single-sample Rproject folder. In order to use the functionalities provided by the source code, some preliminary steps must be executed. Furthermore, we recall that, since we have based the research inside an R project, any considered path is relative to the working directory of the project itself. 

The installation required is summarized in the next steps:

1. Download and install R 4.0
2. Download and install RStudio
3. Download and install Rtools.
4. Run the command
	        
  >  <code>    writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron") </code>

5. Restart R and check that path of command is returned with
  
  >  <code>    Sys.which("make") </code>
  
6. In RStudio, open the project with File > Open project and open the CRIS_single-sample project folder (it must contain the .Rproj file).

7. Install the required libraries by running the script in the src folder of the project:

  >  <code>    source("src/install_libraries.r") </code>
  
8. Load the required libraries by running the script in the src folder of the project:

  >  <code>    source("src/load_libraries.r") </code>
  
*NB*: The step 7 must be repeated whenever a new R session starts, otherwise errors of missing libraries are reported. In order to understand how to apply the implemented functionalities, we suggest to look at the R notebooks (.Rmd files) organized within the *vignette* folder.


