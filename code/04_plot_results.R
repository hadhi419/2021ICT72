source("code/03_analyse.R")

top_probes <- head(results$probe_id, 3)
stopifnot(length(top_probes) == 3L)

# Average repeated measurements within each tissue, condition, and plate.
plate_groups <- with(samples, interaction(tissue, condition, plate, drop = TRUE))
plate_data <- do.call(rbind, lapply(split(seq_len(nrow(samples)), plate_groups),
                                    function(indices) {
                                      sample_info <- samples[indices[1], ]
                                      data.frame(
                                        probe_id = top_probes,
                                        tissue = sample_info$tissue,
                                        condition = sample_info$condition,
                                        plate = sample_info$plate,
                                        expression = rowMeans(log_expression[top_probes, indices, drop = FALSE])
                                      )
                                    }))

conditions <- c("Ground", "Flight")
colours <- c(Roots = "#0072B2", Shoots = "#D55E00")
offsets <- c(Roots = -0.05, Shoots = 0.05)

dir.create("figures", showWarnings = FALSE)
png("figures/top_probe_responses.png", width = 1500, height = 650, res = 150)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), oma = c(2, 0, 0, 0))

for (probe in top_probes) {
  probe_data <- subset(plate_data, probe_id == probe)
  fdr <- results$adj.P.Val[match(probe, results$probe_id)]
  
  plot(
    NA, xlim = c(0.8, 2.2), ylim = range(probe_data$expression),
    xaxt = "n", xlab = "", ylab = "Log2 expression",
    main = paste0(probe, "\nFDR = ", format.pval(fdr, digits = 2))
  )
  axis(1, at = 1:2, labels = conditions)
  
  for (tissue_name in names(colours)) {
    tissue_data <- subset(probe_data, tissue == tissue_name)
    x <- match(tissue_data$condition, conditions) + offsets[tissue_name]
    
    # Open circles show individual plates.
    points(x, tissue_data$expression, col = colours[tissue_name], pch = 1)
    
    # Filled circles and lines show the tissue mean in each condition.
    means <- tapply(tissue_data$expression, tissue_data$condition, mean)
    lines(1:2 + offsets[tissue_name], means[conditions],
          col = colours[tissue_name], type = "b", pch = 19, lwd = 2)
  }
}

mtext("Blue = Roots   |   Orange = Shoots   |   Open circles = plates",
      side = 1, outer = TRUE, line = 0.5, cex = 0.9)
invisible(dev.off())