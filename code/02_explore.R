# Read the data saved by 01_import.R.
datasets <- readRDS("data/GSE56659.rds")

# Match the sample IDs to their four groups.
samples <- data.frame(
  id = paste0("GSM", c(1366461:1366470, 1366481:1366490)),
  tissue = rep(c("Shoots", "Roots"), each = 10),
  condition = rep(rep(c("Flight", "Ground"), each = 5), times = 2)
)

# Find the expression matrix containing all 20 samples.
matches <- vapply(datasets, function(dataset) {
  all(samples$id %in% colnames(Biobase::exprs(dataset)))
}, logical(1))
stopifnot(sum(matches) == 1)

expression <- Biobase::exprs(datasets[[which(matches)]])[, samples$id]

# Keep probes with a finite measurement in every selected sample.
complete_probes <- rowSums(is.finite(expression)) == ncol(expression)
expression <- expression[complete_probes, , drop = FALSE]
stopifnot(nrow(expression) > 1)

# Transpose the matrix so each row represents one sample.
# Log2-transform the expression values before comparing samples.
stopifnot(all(expression > 0))
log_expression <- log2(expression)
pca <- prcomp(t(log_expression))
variance_pct <- 100 * pca$sdev^2 / sum(pca$sdev^2)

colours <- c(Roots = "#0072B2", Shoots = "#D55E00")
symbols <- c(Flight = 19, Ground = 17)

dir.create("figures", showWarnings = FALSE)
png("figures/pca_log2.png", width = 1000, height = 700, res = 125)
plot(pca$x[, 1], pca$x[, 2],
     col = colours[samples$tissue],
     pch = symbols[samples$condition],
     cex = 1.5,
     xlab = sprintf("PC1 (%.1f%%)", variance_pct[1]),
     ylab = sprintf("PC2 (%.1f%%)", variance_pct[2]),
     main = "Gene expression in roots and shoots")
legend("topright",
       legend = c("Roots / Flight", "Roots / Ground",
                  "Shoots / Flight", "Shoots / Ground"),
       col = rep(colours, each = 2),
       pch = rep(symbols, times = 2), bty = "n")
invisible(dev.off())
print(with(samples, table(tissue, condition)))