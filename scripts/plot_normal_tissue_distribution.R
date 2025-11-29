#!/usr/bin/env Rscript
# Plot distribution of normal tissue samples across TCGA project types

library(ggplot2)
library(dplyr)

# Load pre-computed project counts
project_counts <- read.csv("data/processed/gdc_pancan/methyl/normal_counts_by_project.csv",
  row.names = 1
)

# Convert to data frame with proper column names
project_counts <- data.frame(
  project = rownames(project_counts),
  n_samples = project_counts$sample_counts
) %>%
  arrange(desc(n_samples))

message(sprintf("Loaded counts for %d projects", nrow(project_counts)))

message("\n=== Normal Tissue Sample Counts by Project ===")
print(project_counts)

# Create output directory
if (!dir.exists("results/tissue_analysis")) {
  dir.create("results/tissue_analysis", recursive = TRUE)
}

# Bar plot of sample counts per project
bar_plot <- ggplot(project_counts, aes(x = reorder(project, n_samples), y = n_samples)) +
  geom_bar(stat = "identity", fill = "#2c7bb6", alpha = 0.8) +
  geom_text(aes(label = n_samples), hjust = -0.2, size = 5) +
  coord_flip() +
  labs(
    title = "Distribution of Normal Tissue Samples Across TCGA Projects",
    x = "TCGA Project",
    y = "Number of Normal Samples"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    panel.grid.major.y = element_blank()
  )

ggsave(
  filename = "results/tissue_analysis/normal_tissue_by_project.png",
  plot = bar_plot,
  width = 16,
  height = 9,
  dpi = 300
)

message("\nBar plot saved to: results/tissue_analysis/normal_tissue_by_project.png")

# 3. Summary statistics
cat("\n=== Summary Statistics ===\n")
cat(sprintf("Total normal samples: %d\n", sum(project_counts$n_samples)))
cat(sprintf("Number of projects: %d\n", nrow(project_counts)))
cat(sprintf("Mean samples per project: %.1f\n", mean(project_counts$n_samples)))
cat(sprintf("Median samples per project: %.1f\n", median(project_counts$n_samples)))
cat(sprintf("Range: %d - %d\n", min(project_counts$n_samples), max(project_counts$n_samples)))
cat(sprintf("\nTop 5 projects by sample count:\n"))
print(head(project_counts %>% select(project, n_samples), 5))

message("\n=== Analysis complete! ===")
