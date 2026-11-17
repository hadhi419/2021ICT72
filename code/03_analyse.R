# Load the selected samples and log2 expression from the exploration script.
source("code/02_explore.R")
stopifnot(identical(colnames(log_expression), samples$id))

# The GEO sample IDs were selected in plate order within each group.
samples$plate <- rep(c(2, 4, 4, 6, 6), times = 4)

# Each row describes a plate in one condition.
plates <- expand.grid(
  plate = c(2, 4, 6),
  condition = c("Ground", "Flight"),
  stringsAsFactors = FALSE
)

# Calculate root minus shoot log2 expression within each plate.
plate_differences <- vapply(seq_len(nrow(plates)), function(i) {
  on_plate <- samples$plate == plates$plate[i] &
    samples$condition == plates$condition[i]
  
  roots <- on_plate & samples$tissue == "Roots"
  shoots <- on_plate & samples$tissue == "Shoots"
  stopifnot(sum(roots) > 0, sum(roots) == sum(shoots))
  
  rowMeans(log_expression[, roots, drop = FALSE]) -
    rowMeans(log_expression[, shoots, drop = FALSE])
}, numeric(nrow(log_expression)))

rownames(plate_differences) <- rownames(log_expression)
plates$condition <- factor(plates$condition, levels = c("Ground", "Flight"))

# Test how much the root–shoot difference changes from ground to flight.
design <- model.matrix(~ factor(plate) + condition, data = plates)
fit <- limma::lmFit(plate_differences, design)
fit <- limma::eBayes(fit)

results <- limma::topTable(
  fit, coef = "conditionFlight", number = Inf, adjust.method = "BH"
)
results <- data.frame(probe_id = rownames(results), results, row.names = NULL)

dir.create("results", showWarnings = FALSE)
write.csv(results, "results/tissue_by_flight.csv", row.names = FALSE)

cat("Probes tested:", nrow(results), "\n")
cat("Probes with FDR < 0.05:", sum(results$adj.P.Val < 0.05, na.rm = TRUE), "\n")
print(head(results[, c("probe_id", "logFC", "P.Value", "adj.P.Val")], 10))