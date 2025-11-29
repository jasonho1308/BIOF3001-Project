# Plot distribution of Horvath age predictions

library(ggplot2)
library(dplyr)

# Load combined predictions
predictions <- readRDS("results/clock/gdc_pan/gdc_pancan_predictions.rds")

# Filter to samples with valid Horvath residuals
horvath_data <- predictions %>%
  filter(!is.na(Horvath_residuals))

message(paste("Plotting distribution for", nrow(horvath_data), "samples"))

# Calculate quartiles

q1 <- quantile(horvath_data$Horvath_residuals, 0.25, na.rm = TRUE)
q3 <- quantile(horvath_data$Horvath_residuals, 0.75, na.rm = TRUE)
horvath_mean <- mean(horvath_data$Horvath_residuals, na.rm = TRUE)
horvath_median <- median(horvath_data$Horvath_residuals, na.rm = TRUE)

# Create distribution plot
dist_plot <- ggplot(horvath_data, aes(x = Horvath_residuals)) +
  geom_histogram(bins = 50, fill = "#1875B9", color = "white", alpha = 0.7) +
  geom_vline(aes(xintercept = horvath_mean),
    color = "#E75C21", linetype = "solid", linewidth = 2
  ) +
  geom_vline(aes(xintercept = horvath_median),
    color = "#F49E2D", linetype = "dashed", linewidth = 2
  ) +
  geom_vline(aes(xintercept = q1),
    color = "#D079B1", linetype = "dashed", linewidth = 2
  ) +
  geom_vline(aes(xintercept = q3),
    color = "#D079B1", linetype = "dashed", linewidth = 2
  ) +
  labs(
    title = "Distribution of Horvath Age Acceleration (Residuals)",
    subtitle = sprintf(
      "Mean = %.1f (red), Median = %.1f (orange), Q1 = %.1f, Q3 = %.1f (blue)",
      horvath_mean, horvath_median, q1, q3
    ),
    x = "Horvath Age Acceleration (years)",
    y = "Count"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 24, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18)
  )

# Save plot
ggsave(
  filename = "results/clock/gdc_pan/horvath_residuals_distribution.png",
  plot = dist_plot,
  width = 16,
  height = 9,
  dpi = 300
)


# Print summary statistics
cat("\n=== Horvath Age Acceleration (Residuals) Summary ===\n")
cat(sprintf("Sample size: %d\n", nrow(horvath_data)))
cat(sprintf("Mean: %.2f years\n", mean(horvath_data$Horvath_residuals, na.rm = TRUE)))
cat(sprintf("Median: %.2f years\n", median(horvath_data$Horvath_residuals, na.rm = TRUE)))
cat(sprintf("SD: %.2f years\n", sd(horvath_data$Horvath_residuals, na.rm = TRUE)))
cat(sprintf("Q1: %.2f years\n", q1))
cat(sprintf("Q3: %.2f years\n", q3))
cat(sprintf("IQR: %.2f years\n", q3 - q1))
cat(sprintf(
  "Range: %.2f - %.2f years\n",
  min(horvath_data$Horvath_residuals, na.rm = TRUE),
  max(horvath_data$Horvath_residuals, na.rm = TRUE)
))

message("Distribution plots saved!")
