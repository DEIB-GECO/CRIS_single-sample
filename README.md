# CRIS_single-sample

The work presented in this repository allows to apply single-sample gene-expression-based classifiers for CRIS stratification of colorectal cancer patients. The work has been organized in the [**CRIS_single-sample**](https://github.com/DEIB-GECO/CRIS_single-sample/tree/main/CRIS_single-sample) Rproject folder, which includes the source code as well as several vignettes showing how to use such code and different scripts useful to replicate the training, testing and biological validations of all the models under evaluation. Further details about the content of the [**CRIS_single-sample**](https://github.com/DEIB-GECO/CRIS_single-sample/tree/main/CRIS_single-sample) Rproject folder are reported on the [**Project structure.md**](https://github.com/DEIB-GECO/CRIS_single-sample/tree/main/CRIS_single-sample/Project_structure.md) file provided within the same folder, as to be easily accessible also on R environment.

In order to use the functionalities provided by the source code, some preliminary steps must be executed. Furthermore, we recall that, since we have based the research inside an R project, any considered path is relative to the working directory of the project itself. 

The installation required is summarized in the next steps:

1. Download and install R 4.0
2. Download and install RStudio
3. Download and install Rtools.
4. Run the command
	        
  > <code>    writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron") </code>

5. Restart R and check that path of command is returned with
  
  > <code>    Sys.which("make") </code>

6. Download the zip file of the repository; in RStudio, open the project with File > Open project and open the **CRIS_single-sample.Rproj** file located in the folder [*CRIS_single-sample*](https://github.com/DEIB-GECO/CRIS_single-sample/tree/main/CRIS_single-sample).

7. To install the original CRIS classifier, download the tar.gz file of the library (here provided in [CRISclassifier](https://github.com/DEIB-GECO/CRIS_single-sample/blob/main/CRISclassifier_1.0.0.tar.gz) and then run the command

  > <code>	install.packages("path_to_CRIS-classifier.tar.gz") </code>
  
8. To install the *utiml* package, download the tar.gz file of the library (from [utiml](https://cran.r-project.org/src/contrib/Archive/utiml/utiml_0.1.6.tar.gz) and then run the command

  > <code>	install.packages("path_to_utiml.tar.gz") </code>

9. Install all the other required libraries by running the script in the src folder of the project:

  > <code>	source("src/install_libraries.r") </code>
  
10. Load the required libraries by running the script in the src folder of the project:

  > <code>    source("src/load_libraries.r") </code>

11. Install any further dependecies if required 

12. Download the **data.zip** file from **XX**, containing all the needed data and extract it as it within the project folder [*CRIS_single-sample*](https://github.com/DEIB-GECO/CRIS_single-sample/tree/main/CRIS_single-sample). Notice that such data are required to run all the available scripts and vignettes. 

*NB*: The step 10 must be repeated whenever a new R session starts, otherwise errors of missing libraries are reported. In order to better understand how to apply the implemented functionalities, we suggest to look at the R notebooks (.Rmd files) organized within the [*vignette*](https://github.com/DEIB-GECO/CRIS_single-sample/tree/main/CRIS_single-sample/vignette) folder.

The reference work (NTP classifier, TSP classifier, published results of NTP and TSP, CRIS signature) is *Isella C, Brundu F, Bellomo SE, Galimi F, Zanella E, Porporato R, Petti C, Fiori A, Orzan F, Senetta R, Boccaccio C, Ficarra E, Marchionni L, Trusolino L, Medico E, Bertotti A. Selective analysis of cancer-cell intrinsic transcriptional traits defines novel clinically relevant subtypes of colorectal cancer. Nat Commun. 2017 May 31;8:15107. doi: 10.1038/ncomms15107. PMID: 28561063; PMCID: PMC5499209.*
