library(tidyverse)
library(RColorBrewer)

# 1. Read and process data
neo <- read_tsv("Neojunction/neo_junctions.annotated.bed", col_types = cols())

# 2. Create a better data structure
junction_data <- neo %>%
  # Create junction labels with chromosome and position info
  mutate(
    junction_label = paste0("chr", chr, ":", start, "-", end),
    pos = (start + end) / 2
  ) %>%
  # Expand samples
  separate_rows(samples, sep = ",") %>%
  rename(sample = samples) %>%
  select(junction_label, chr, pos, sample, total_reads, samples_expressed, frequency)

# 3. Create a presence/absence matrix but keep read counts
plot_data <- junction_data %>%
  # Order chromosomes properly
  mutate(chr = factor(chr, levels = c(1:22, "X", "Y"))) %>%
  arrange(chr, pos)

# 4. Create an informative heatmap-style plot
p1 <- ggplot(plot_data, aes(x = reorder(junction_label, pos), y = sample)) +
  geom_tile(aes(fill = log10(total_reads + 1)), color = "white", size = 0.1) +
  scale_fill_gradient(low = "white", high = "darkred", 
                      name = "log10(reads+1)") +
  facet_wrap(~chr, scales = "free_x", nrow = 2) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.1, "lines")
  ) +
  labs(
    title = "Neo-junctions: Read Support Across Samples by Chromosome",
    subtitle = paste0("32 novel junctions across ", length(unique(plot_data$sample)), " samples"),
    x = "Junction position (ordered by chromosome)",
    y = "Sample ID"
  )

# 5. Create a summary bar plot
summary_data <- neo %>%
  mutate(chr = factor(chr, levels = c(1:22, "X", "Y"))) %>%
  arrange(chr, start)

p2 <- ggplot(summary_data, aes(x = reorder(paste0("chr", chr, ":", start, "-", end), start))) +
  geom_col(aes(y = samples_expressed), fill = "steelblue", alpha = 0.7) +
  geom_text(aes(y = samples_expressed + 2, label = samples_expressed), 
            size = 3, angle = 90) +
  facet_wrap(~chr, scales = "free_x", nrow = 1) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.1, "lines")
  ) +
  labs(
    title = "Number of Samples per Neo-junction",
    x = "Junction position (ordered by chromosome)",
    y = "Number of samples"
  )

# 6. Create a frequency distribution plot
p3 <- ggplot(neo, aes(x = samples_expressed)) +
  geom_histogram(bins = 15, fill = "darkgreen", alpha = 0.7, color = "white") +
  geom_vline(xintercept = median(neo$samples_expressed), 
             color = "red", linetype = "dashed", size = 1) +
  theme_classic() +
  labs(
    title = "Distribution of Neo-junction Recurrence",
    subtitle = paste0("Median: ", median(neo$samples_expressed), " samples per junction"),
    x = "Number of samples with junction",
    y = "Number of junctions"
  )


# 7. Expanded Neo_juncton per sample count to one row.
junction_samples <- neo %>%
  select(key, samples) %>%
  separate_rows(samples, sep = ",") %>%
  rename(sample = samples)

# Count neo-junctions per sample
sample_counts <- junction_samples %>%
  count(sample, name = "neo_junction_count") %>%
  arrange(desc(neo_junction_count))

# Make bar plot with sample names on x-axis
p <- ggplot(sample_counts, aes(x = factor(sample, levels = sample), y = neo_junction_count)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = neo_junction_count), vjust = -0.5, size = 3) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Number of Neo-Junctions per Sample",
    y = "Neo-Junction Count"
  )

# 8. Fixed Circos Plot with Chr link sampple counts  

library(circlize)
library(readr)
library(dplyr)
library(RColorBrewer)

# Process data: sample-chromosome connections with counts
sample_chr_data <- neo %>%
  select(chr, samples, total_reads) %>%
  separate_rows(samples, sep = ",") %>%
  rename(sample = samples) %>%
  group_by(sample, chr) %>%
  summarise(
    junction_count = n(),
    total_reads = sum(total_reads, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_reads))

# Create hierarchical sample groups
sample_summary <- sample_chr_data %>%
  group_by(sample) %>%
  summarise(
    total_reads_sample = sum(total_reads),
    .groups = "drop"
  ) %>%
  mutate(
    activity_group = case_when(
      total_reads_sample >= quantile(total_reads_sample, 0.8) ~ "High",
      total_reads_sample >= quantile(total_reads_sample, 0.5) ~ "Medium", 
      total_reads_sample >= quantile(total_reads_sample, 0.2) ~ "Low",
      TRUE ~ "Very Low"
    )
  )

# Create color schemes
chr_colors <- setNames(
  rainbow(length(unique(sample_chr_data$chr))), 
  sort(unique(sample_chr_data$chr))
)

activity_colors <- c("High" = "#D32F2F", "Medium" = "#FF9800", 
                     "Low" = "#4CAF50", "Very Low" = "#2196F3")

sample_colors <- sample_summary %>%
  mutate(base_color = activity_colors[activity_group]) %>%
  select(sample, base_color) %>%
  deframe()

# Create adjacency matrix
adj_matrix <- sample_chr_data %>%
  select(sample, chr, total_reads) %>%
  pivot_wider(names_from = chr, values_from = total_reads, values_fill = 0) %>%
  column_to_rownames("sample") %>%
  as.matrix()

# Sort by total activity
adj_matrix <- adj_matrix[order(rowSums(adj_matrix), decreasing = TRUE), ]

cat("Matrix dimensions:", nrow(adj_matrix), "x", ncol(adj_matrix), "\n")

# Calculate appropriate gaps - CRITICAL FIX
n_chroms <- ncol(adj_matrix)
n_samples <- nrow(adj_matrix) 
total_sectors <- n_chroms + n_samples

# Make sure total gaps < 360 degrees
small_gap <- 0.5
big_gap <- 3
total_gap <- (total_sectors - 2) * small_gap + 2 * big_gap

cat("Total sectors:", total_sectors, "\n")
cat("Total gap degrees:", total_gap, "\n")

# If gaps too large, reduce them
if(total_gap > 300) {
  small_gap <- 0.2
  big_gap <- 2
  cat("Reduced gaps to fit\n")
}

# Create the chord diagram
png("fixed_sample_chr_circos.png", width = 3000, height = 3000, res = 300)

circos.clear()

# FIXED: Calculate exact gaps
gap_vector <- c(
  rep(small_gap, n_chroms-1), big_gap,      # chromosomes + separator
  rep(small_gap, n_samples-1), big_gap      # samples + final separator  
)

circos.par(
  gap.after = gap_vector,
  start.degree = 90
)

# Create colors for all sectors
grid_colors <- c(
  chr_colors,  
  sample_colors[rownames(adj_matrix)]
)

# Draw simplified chord diagram (remove problematic parameters)
chordDiagram(
  adj_matrix,
  grid.col = grid_colors,
  transparency = 0.5,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.1)
)

# Add labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.name = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  
  if(sector.name %in% names(chr_colors)) {
    # Chromosome labels - bold and larger
    circos.text(mean(xlim), ylim[1] + 0.2, paste0("Chr", sector.name), 
                cex = 0.8, facing = "clockwise", niceFacing = TRUE, 
                col = "black", font = 2)
  } else {
    # Sample labels - smaller
    circos.text(mean(xlim), ylim[1] + 0.05, sector.name, 
                cex = 0.3, facing = "clockwise", niceFacing = TRUE, 
                col = "navy")
  }
}, bg.border = NA)

# Add legends
legend("topright", 
       legend = names(activity_colors), 
       fill = activity_colors,
       title = "Sample frequency", 
       cex = 0.8, bg = "white")

# Add title
text(x = 0, y = 1.0, 
     labels = "Sample-Chromosome Neo-Junction Network", 
     cex = 1.4, font = 2, col = "black")
text(x = 0, y = 1.15, 
     labels = "Link Width = Read Counts", 
     cex = 1.0, col = "darkblue")

dev.off()
circos.clear()

cat("\nSuccessfully created: fixed_sample_chr_circos.png\n")


# 9. Save all plots
ggsave("neo_junctions_heatmap.png", plot = p1, width = 16, height = 29, dpi = 300)
ggsave("neo_junctions_summary.png", plot = p2, width = 16, height = 6, dpi = 300)
ggsave("neo_junctions_distribution.png", plot = p3, width = 8, height = 6, dpi = 300)
ggsave("neo_junctions_per_sample.png", plot = p, width = 20, height = 6, dpi = 300)

# 10. Print summary statistics
cat("Summary of Neo-junction Results:\n")
cat("Total junctions found:", nrow(neo), "\n")
cat("Total samples analyzed:", length(unique(junction_data$sample)), "\n")
cat("Chromosomes involved:", paste(sort(unique(neo$chr)), collapse = ", "), "\n")
cat("Median samples per junction:", median(neo$samples_expressed), "\n")
cat("Range of read support:", min(neo$total_reads), "-", max(neo$total_reads), "\n")


