# Description -----------------------------------------------------------------

# Useful constants used throughout the code

# Constants ------------------------------------------------------------------

# Number of pdx bathces to consider
N_PDX_BATCHES        <- as.integer(6)

# IDs attributes
ALIQUOT_LABEL        <- "aliquot_id"
PATIENT_LABEL        <- "patient_id"
SAMPLE_LABEL         <- "sample_id"

# Metadata from which the ids are taken in GMQL_GRCH38
ALIQUOT_META         <- "gdc__aliquots__submitter_id"
PATIENT_META         <- "biospecimen__shared__bcr_patient_barcode"
SAMPLE_META          <- "biospecimen__bio__bcr_sample_barcode"

# Lengths of the different types of ids in TCGA
ALIQUOT_LABEL_LENGTH <- as.integer(28)
SAMPLE_LABEL_LENGTH  <- as.integer(16)
PATIENT_LABEL_LENGTH <- as.integer(12)


# CRIS Class attributes
CLASS_LABEL          <- "predict.label2"
CRIS_CLASSES         <- c("CRIS-A", "CRIS-B", "CRIS-C", "CRIS-D", "CRIS-E")
F_CRIS_CLASSES       <- ordered(CRIS_CLASSES)
N_CLASSES            <- nlevels(F_CRIS_CLASSES)

# Distance attributes for NTP
BEST_DISTANCE_LABEL  <- "dist.to.template"
CLASS_DISTANCE_LABEL <- paste("dist", gsub("-",".",CRIS_CLASSES), sep = "")

# BH.FDR attributes for NTP
BEST_FDR_LABEL       <- "BH.FDR"
CLASS_FDR_LABEL      <- paste("BH.FDR", gsub("-",".",CRIS_CLASSES), sep = ".")
BH_FDR_THRESHOLD     <- as.numeric(0.2)
