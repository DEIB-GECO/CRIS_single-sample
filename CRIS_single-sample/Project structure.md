# Structure of CRIS_single-sample project

The present file is used to describe the structure of the repository.
In particular:

1. *additional files* folder: contains additional files indirectly related to the project. In particular, there is an excel forest plot template that can be used to draw forest plots starting from the results obtained by the source code.
It is sufficient to copy the obtained dataframe in the template file (following its structure) to see how the forest plots update automatically. **NB**: the labels present in the forest plots (representing p-values) must be updated manually.
2. *scripts* folder: contains scripts that can be run to obtain training/testing/biological validation of all models.
In particular, the logical order in which they should be executed is:
	1. *ntp_ml_thresholds.r* : generation of the thresholds for the assignment of multiple classes in NTP.
	2. *single_label_training.r* and *ml_problem_transformation_training.r*: training of single-label and problem-transformation multi-label classifiers.
	3. *single_label_thresholds.r*  and *ml_problem_transformation_thresholds.r*: computation of thresholds for assignment of multiple classes starting from scores of single-label and problem-transformation multi-label classifiers.
	4. *single_label_testing.r*, *ml_algorithm_adaptation_testing.r*  and *ml_problem_transformation_testing.r*: testing of single-label and multi-label classifiers.
	5. *kaplan_meier.r* and *forest_plots.r*: computation of Kaplan-Meier curves and forest plots for all the main methods.
3. *src* folder: source code (divided in several subfolders. The *utils* subfolder contains the *paths_db.xlsx* excel file with the paths of all accessible data.
4. *vignette* folder: contains R notebooks that explain how to use the code. In particular, the logical order in which they should be read is:
	1. *data_management/load_datasets.rmd* and *data_management/load_references.rmd*: vignettes that explain how to read the datasets and the references (e.g. NTP reference classifcation).
	2. *data_management/filter_expression.rmd*: vignette explaining how to filter the datasets by samples/genes.
	3. *data_management/data_pipelines.rmd*: vignette explaining the classes used to manage the data and the functions to load the data ready for training/testing of all classifiers (both for PDX and TCGA)
	4. *classification/ntp_replication.rmd*: vignette explanining how to apply NTP on the datasets (replication and biological validation).
	5. *classification/tsp_replication.rmd*: vignette explanining how to apply TSP on the datasets (replication and biological validation).
	6. *classification/sl_replication.rmd*: vignette explanining how to apply single-label classifiers on the datasets (training, testing and biological validation).
	7. *classification/sl_as_ml_replication.rmd*: vignette explanining how to apply single-label classifiers adapted to multi-label on the datasets (training, threshold computation for multi-label assignment, testing and biological validation).
	8. *classification/ml_pt_replication.rmd*: vignette explanining how to apply multi-label problem transformation classifiers on the datasets (training, threshold computation for multi-label assignment, testing and biological validation).

**NB**: the execution of the R notebooks requires the *data* folder, which must be required separately and extracted at the same level of *src*, *vignette*, *scripts* and *additional files* folders.