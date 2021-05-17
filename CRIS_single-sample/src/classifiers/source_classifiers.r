# SOURCE CLASSIFIERS HANDLING FOLDER
library(here)

# Data
source(here('src', 'classifiers', 'data.r'))
source(here('src', 'classifiers', 'cris_data.r'))
source(here('src', 'classifiers', 'ml_data.r'))
source(here('src', 'classifiers', 'sl_data.r'))

# Feature selection
source(here('src', 'classifiers', 'feature_selection.r'))

# Algorithms settings
source(here('src', 'classifiers', 'settings.r'))

# Classifiers
source(here('src', 'classifiers', 'ntp_classifier_code.r'))
source(here('src', 'classifiers', 'ntp.r'))
source(here('src', 'classifiers', 'tsp.r'))
source(here('src', 'classifiers', 'single_sample_classifier.r'))  # Base class for sl and ml
source(here('src', 'classifiers', 'sl.r'))
source(here('src', 'classifiers', 'ml.r'))

# Metrics
source(here('src', 'classifiers', 'metrics.r'))
source(here('src', 'classifiers', 'sl_metrics.r'))
source(here('src', 'classifiers', 'ml_metrics.r'))

# Classifiers settings
source(here('src', 'classifiers', 'classifiers_settings.r'))
source(here('src', 'classifiers', 'saving_classifiers.r'))
