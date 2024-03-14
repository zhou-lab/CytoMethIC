
library(tidyverse)

metadata_path <- system.file("extdata", "metadata.csv", package = "CytoMethIC")

metadata <- read.csv(metadata_path)

print_model_summary <- function(model_data) {
  cat("Title:", model_data[["Title"]], "\n")
  cat("Description:", model_data[["Description"]], "\n")
  cat("BiocViews:", model_data[["BiocViews"]], "\n")
  cat("Source URL:", model_data[["SourceUrl"]], "\n")
  cat("R Data Path:", model_data[["RDataPath"]], "\n")
  cat("Maintainer:", model_data[["Maintainer"]], "\n")
  cat("---------\n")
}

# Print summary for each model
apply(as.data.frame(metadata), 1, print_model_summary)


print(paste("All models for CNS tumor classification were trained on publically available samples" ,
" from Capper, D. et al. DNA methylation-based classification of central nervous system tumours. These models use
features from the probabilistic feature importance report of an initial RandomForest model.",
"All other models were trained using publically available samples from The Cancer Genome Atlas,"
 ,"and use the most variable CpG probes in the dataset as features."))







