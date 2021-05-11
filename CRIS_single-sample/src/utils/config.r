
# Description ------------------------------------------------------

# The file contains constants useful for the configuration of the project

# Constants --------------------------------------------------------

# File with database of paths
.FILE_PATHS_DB_PATH <- here('src','utils','paths_db.xlsx')

# Configuration for NTP replication  -------------------------------

# Number of resamplings for NTP classification
#.DEF_N_RESEMPL <- 100000
.DEF_N_RESEMPL <- 10

# Output folder (where NTP results are stored)
.DEF_OUTPUT    <- here('output/replication')
dir_create(.DEF_OUTPUT)


# Configuration for the forest plots -------------------------------

# Attributes to be tested
.FOREST_ATTRIBUTES <- list(
  tcga = c('n_mucinous','kras_mutated','n_msi_h', 'recurred_before_36'),
  pdx  = c('n_sensible')
)

# Name of each test
.FOREST_LABEL_ATTRIBUTES <- list(
  tcga = c('Mucinous','KRAS mutated', 'MSI-High', 'Recurred within 36 months'),
  pdx  = c('Cetuximab sensitivity')
)

# For each attribute, classes for which the forest plot is required
.FOREST_CLASSES <- list(
  n_mucinous = c('CRIS-A'),
  kras_mutated = c('CRIS-A','CRIS-C'),
  n_msi_h = c('CRIS-A'),
  recurred_before_36 = c('CRIS-B'),
  n_sensible = c('CRIS-C')
)