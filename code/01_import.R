install.packages("BiocManager")


BiocManager::install(
  c("GEOquery", "Biobase", "limma"),
  ask = FALSE,
  update = FALSE
)

library(GEOquery)

matrix_file <- "data/GSE56659_series_matrix.txt.gz"

if (!file.exists(matrix_file)) {
  stop("Place the downloaded .txt.gz file in the data folder under ", getwd())
}

dataset <- GEOquery::getGEO(filename = matrix_file, getGPL = FALSE)

# 02_explore.R expects a list of datasetas.
datasets <- list(GSE56659 = dataset)
saveRDS(datasets, "data/GSE56659.rds")

# Print the accession and title for each sample.
for (dataset in datasets) {
  metadata <- Biobase::pData(dataset)
  print(metadata[, c("geo_accession", "title")], row.names = FALSE)
}