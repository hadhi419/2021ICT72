install.packages("BiocManager")


BiocManager::install(
  c("GEOquery", "Biobase", "limma"),
  ask = FALSE,
  update = FALSE
)

library(GEOquery)

dir.create("data", showWarnings = FALSE)

# Download the processed expression matrices from NCBI GEO.
datasets <- getGEO("GSE56659", GSEMatrix = TRUE, destdir = "data")

# Show each downloaded matrix and how many samples it contains.
print(data.frame(
  matrix = names(datasets),
  samples = vapply(datasets, ncol, integer(1))
))

# Print the accession and title for each sample.
for (dataset in datasets) {
  metadata <- Biobase::pData(dataset)
  print(metadata[, c("geo_accession", "title")], row.names = FALSE)
}