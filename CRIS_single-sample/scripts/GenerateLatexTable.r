library(here)
library(tidyverse)

# Read the comparison file
path <- here("additional files/Compare_ntp_only.xlsx")
ref  <- readxl::read_excel(path, sheet = 2)

# Save and then remove row with data labels (TCGA/PDX)
data_labels <- ref[1,]
ref <- ref[-1, ]

# Round the numbers to 3 digits
for (c in colnames(ref)[-1]){
  ref[,c] <- round(as.numeric(unlist(ref[,c])),3)
}
ref <- as.matrix(ref)

# Put in bold numbers above 0.8
for (i in seq(1,nrow(ref))){
  # m <- max(ref[i,-1], na.rm = TRUE) %>% as.numeric()
  # print(m)
  for (c in colnames(ref)[-c(1,6)])
    if (!is.na(ref[i,c]) & any(as.numeric(ref[i,c]) >= 0.80))
      ref[i,c] <- paste("\\textbf{", ref[i,c], "}", sep = '')

}

# Re-add the data labels and add the algorithm names
ref <- rbind(data_labels, ref)
ref <- rbind(colnames(ref),ref)

# Create the tables
row_str <- ''
for (n in seq(nrow(ref))){
  row_str <- row_str %>% paste(ref[n, ], collapse = ' & ') %>% cat('\\\\\n',"\\" ,'hline',' \n', sep = '')
  
}



