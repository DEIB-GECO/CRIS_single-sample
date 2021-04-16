# SOURCE OF UTILS SCRIPTS
source(here('src', 'utils','prints.r'))
source(here('src', 'utils','config.r'))
source(here('src', 'utils', 'paths_refactored.r'))

# Instantiate the path loader
path_loader <- PathLoader$new(.FILE_PATHS_DB_PATH)
