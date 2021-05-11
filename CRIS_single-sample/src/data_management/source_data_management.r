# SOURCE SAMPLES HANDLING FOLDER
library(here)

# Comparators and adapters
source(here('src', 'data_management', 'comparator.r'))
source(here('src', 'data_management', 'comparator_ntp.r'))
source(here('src', 'data_management', 'adapter.r'))

# Filters
source(here('src', 'data_management', 'filter.r'))
source(here('src', 'data_management', 'filter_genes.r'))
source(here('src', 'data_management', 'filter_samples.r'))

# Filters

source(here('src', 'data_management', 'samples_manager.r'))
source(here('src', 'data_management', 'summary_maker.r'))
source(here('src', 'data_management', 'summary_maker_tcga.r'))
source(here('src', 'data_management', 'summary_maker_pdx.r'))