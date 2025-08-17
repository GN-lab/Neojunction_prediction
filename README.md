
####################################################################
# To generate Annotation files
####################################################################

#!/usr/bin/env bash

# For canonical_junctions.txt
zcat ../../../Ref/genome.gtf.gz | awk -F'\t' '
$3 == "exon" && $9 ~ "gene_biotype \"protein_coding\"" {
    # Extract gene_name
    gene_name = "";
    split($9, attrs, ";");
    for (i in attrs) {
        if (attrs[i] ~ /gene_name/) {
            split(attrs[i], gn, "\"");
            gene_name = gn[2];
            break;
        }
    }
    # Print in BED format
    if (gene_name != "") {
        print $1 "\t" $4-1 "\t" $5 "\t" gene_name "\t0\t" $7;
    }
}' | sort -k1,1 -k2,2n | uniq > canonical_junctions.txt

# For protein_coding.bed
zcat ../../../Ref/genome.gtf.gz | awk -F'\t' '
$3 == "gene" && $9 ~ "gene_biotype \"protein_coding\"" {
    # Extract gene_name
    gene_name = "";
    split($9, attrs, ";");
    for (i in attrs) {
        if (attrs[i] ~ /gene_name/) {
            split(attrs[i], gn, "\"");
            gene_name = gn[2];
            break;
        }
    }
    # Print in BED format
    if (gene_name != "") {
        print $1 "\t" $4-1 "\t" $5 "\t" gene_name "\t0\t" $7;
    }
}' | sort -k1,1 -k2,2n | uniq > protein_coding.bed

# For Gtex junction file
zcat GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz \
  | tail -n +4 \
  | cut -f1 \
  | awk -F'_' '{
      chrom=$1; sub(/^chr/, "", chrom);
      strand = ($4 == "1") ? "+" : "-";
      print chrom ":" $2 ":" $3 ":" strand;
    }' \
  > gtex_junctions_fixed.txt


####################################################################
# Neo-Junction Pipeline - Optimized for STAR 2.7.10a
# Parameters: --alignSJoverhangMin 8 --alignIntronMin 20 --outSJfilterOverhangMin 15 20 20 20
####################################################################


#!/usr/bin/env bash
set -euo pipefail

# ------------------ Define and activate venv ------------------

WORKDIR=$(pwd)
source "${WORKDIR}/venv/bin/activate"
PYTHON="${WORKDIR}/venv/bin/python3"

module load bedtools/2.29.2

# ------------------------- Configuration --------------------------

export INPUTDIR="/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/nextflow_rnaseq/output_hartwig/star_salmon/log/"
export OUTDIR="${WORKDIR}/results"
export ANNOT_DIR="${WORKDIR}/annotations"
export SAMPLES_TXT="${WORKDIR}/samples.txt"
export CANON_SET="${ANNOT_DIR}/canonical_junctions.txt"
export GTEX_SET="${ANNOT_DIR}/gtex_junctions.txt"
export PROTEIN_BED="${ANNOT_DIR}/protein_coding.bed"  # Must be defined!

# STAR-aligned filters
export MIN_JUNC_OVERHANG=8        # --alignSJoverhangMin
export MIN_INTRON_SIZE=20         # --alignIntronMin
export MAX_INTRON_SIZE=1000000    # --alignIntronMax

# Neo-junction filters
export MIN_READS_PER_JUNC=10      # Per-sample read support
export MIN_TOTAL_READS=20         # Cohort-wide depth
export MIN_FREQUENCY=0.01         # 1% spliced frequency
export MIN_SAMPLES_PCT=0.10       # 10% samples for "public" junctions
export MAX_GTEX_FREQUENCY=0.01    # 1% GTEx occurrence

# ------------------------- Sample List Generation --------------------------
if [[ ! -s "${SAMPLES_TXT}" ]]; then
  echo "[INFO] Generating samples.txt from *.SJ.out.tab files..."
  find "${INPUTDIR}" -type f -name "*.SJ.out.tab" \
    | awk -F'/' '{fname=$NF; sub(/\.SJ\.out\.tab$/, "", fname); print fname}' \
    | sort -u > "${SAMPLES_TXT}"
  echo "[INFO] Created ${SAMPLES_TXT} with $(wc -l < "${SAMPLES_TXT}") samples."
fi

# ------------------------- Sample Processing ----------------------

process_sample() {
  local sample=$1
  local sj_file="${INPUTDIR}/${sample}.SJ.out.tab"
  local norm_bed="${OUTDIR}/${sample}.juncs.bed"

  echo "[PROC] Processing ${sample}"
  awk -v MIN_OVERHANG=$MIN_JUNC_OVERHANG \
      -v MIN_INTRON=$MIN_INTRON_SIZE \
      -v MAX_INTRON=$MAX_INTRON_SIZE \
      'BEGIN {OFS="\t"} 
       $6 == 0 && $9 >= MIN_OVERHANG && ($3-$2) >= MIN_INTRON && ($3-$2) <= MAX_INTRON && ($7+$8) >= ENVIRON["MIN_READS_PER_JUNC"] {
         strand = ($4 == 1) ? "+" : ($4 == 2) ? "-" : ".";
         print $1, $2-1, $3, "junction_"NR, $7+$8, strand
       }' "${sj_file}" > "${norm_bed}"
}

while read -r sample; do
  process_sample "$sample"
done < "$SAMPLES_TXT"

# ------------------------- Aggregate & Filter ---------------------
echo "[AGG] Aggregating junctions..."

"${PYTHON}" << 'PYCODE'
import pandas as pd
import os

# Load config
outdir = os.environ['OUTDIR']
samples = [s.strip() for s in open(os.environ['SAMPLES_TXT'])]
total_samples = len(samples)
min_samples = int(max(2, float(os.environ['MIN_SAMPLES_PCT']) * total_samples))

# Load all junctions
dfs = []
for sample in samples:
    bed_file = f"{outdir}/{sample}.juncs.bed"
    if os.path.exists(bed_file):
        df = pd.read_csv(bed_file, sep='\t', 
                         names=['chr','start','end','name','reads','strand'])
        df['sample'] = sample
        dfs.append(df)

# Aggregate (key uses 1-based start)
agg = (pd.concat(dfs)
       .groupby(['chr','start','end','strand'])
       .agg(total_reads=('reads','sum'),
            samples_expressed=('sample','nunique'))
       .reset_index())
agg['key'] = agg.apply(lambda x: f"{x['chr']}:{x['start']+1}:{x['end']}:{x['strand']}", axis=1)

# Load canonical junctions and mark status
with open(os.environ['CANON_SET']) as f:
    canon_juncs = set(line.strip() for line in f)
agg['is_canonical'] = agg['key'].isin(canon_juncs)

# Calculate canonical reads for frequency calculation
canon_df = agg[agg['is_canonical']].set_index('key')
def get_canonical_reads(row):
    return canon_df['total_reads'][row['key']] if row['key'] in canon_df.index else 0
agg['canonical_reads'] = agg.apply(get_canonical_reads, axis=1)
agg['frequency'] = agg['total_reads'] / (agg['total_reads'] + agg['canonical_reads'].replace(0, 1))

# Save temp file for bedtools
agg[['chr','start','end','key','total_reads','strand']].to_csv(
    f"{outdir}/temp_junctions.bed", sep='\t', index=False, header=False)

# Protein-coding filter (requires 0-based BED)
os.system(f"bedtools intersect -a {outdir}/temp_junctions.bed -b {os.environ['PROTEIN_BED']} -wa -u > {outdir}/protein_filtered.bed")
protein_keys = set()
protein_filtered_path = f"{outdir}/protein_filtered.bed"
if os.path.exists(protein_filtered_path) and os.path.getsize(protein_filtered_path) > 0:
    protein_keys = set(pd.read_csv(protein_filtered_path, sep='\t', header=None)[3])

# GTEx filter
gtex_juncs = set()
if os.path.exists(os.environ['GTEX_SET']) and os.path.getsize(os.environ['GTEX_SET']) > 0:
    with open(os.environ['GTEX_SET']) as f:
        gtex_juncs = set(line.strip() for line in f if line.strip())

# Apply final filters
neo_juncs = agg[
    (agg['total_reads'] >= int(os.environ['MIN_TOTAL_READS'])) &
    (~agg['is_canonical']) &  # Exclude canonical
    (agg['key'].isin(protein_keys)) &  # Protein-coding only
    (~agg['key'].isin(gtex_juncs)) &  # Not in GTEx normals
    (agg['frequency'] > float(os.environ['MIN_FREQUENCY'])) &  # Frequency >1%
    (agg['samples_expressed'] >= min_samples)  # Public junctions: ≥10% samples
]

# Save results (1-based coordinates) with frequency
neo_juncs[['chr','start','end','strand','total_reads','samples_expressed','frequency']]\
    .to_csv(f"{outdir}/neo_junctions.bed", sep='\t', index=False, header=True)

print(f"[RESULT] Found {len(neo_juncs)} neo-junctions passing filters")
PYCODE
echo "[DONE] Results saved to ${OUTDIR}/neo_junctions.bed"


####################################################################
# To generate Plots using R stuido 
####################################################################


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


<img width="4800" height="8700" alt="neo_junctions_heatmap" src="https://github.com/user-attachments/assets/6dc8adbb-1b58-4688-b1b3-cb1e86c997ec" />



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


<img width="4800" height="1800" alt="neo_junctions_summary" src="https://github.com/user-attachments/assets/f5394c8c-f46c-449c-b9bf-090689ac6239" />




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

<img width="2400" height="1800" alt="neo_junctions_distribution" src="https://github.com/user-attachments/assets/6effd4da-dcad-40cc-8550-bc0fc4770c5f" />




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



<img width="6000" height="1800" alt="neo_junctions_per_sample" src="https://github.com/user-attachments/assets/a7de233a-9675-4e7d-ab00-a97fd6846537" />



# 8. Fixed Circos Plot Chr with link counts  
library(circlize)
library(readr)
library(dplyr)
library(RColorBrewer)

# Read the data
neo <- read_tsv("Neojunction/neo_junctions.annotated.bed", col_types = cols())

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
       title = "Sample Activity", 
       cex = 0.8, bg = "white")

# Add title
title("Sample-Chromosome Neo-Junction Network\nLink Width = Read Counts", 
      cex.main = 1.2, font.main = 2)

dev.off()
circos.clear()

cat("\nSuccessfully created: fixed_sample_chr_circos.png\n")


<img width="3000" height="3000" alt="fixed_sample_chr_circos" src="https://github.com/user-attachments/assets/aebff2ce-4701-4ffb-ab66-9a1fa8fa035b" />




# 9. Save all plots
ggsave("neo_junctions_heatmap.png", plot = p1, width = 16, height = 29, dpi = 300)
ggsave("neo_junctions_summary.png", plot = p2, width = 16, height = 6, dpi = 300)
ggsave("neo_junctions_distribution.png", plot = p3, width = 8, height = 6, dpi = 300)
ggsave("neo_junctions_per_sample.png", plot = p, width = 20, height = 6, dpi = 300)

# 10. Print summary statistics
> cat("Summary of Neo-junction Results:\n")
Summary of Neo-junction Results:
> cat("Total junctions found:", nrow(neo), "\n")
Total junctions found: 32 
> cat("Total samples analyzed:", length(unique(junction_data$sample)), "\n")
Total samples analyzed: 282 
> cat("Chromosomes involved:", paste(sort(unique(neo$chr)), collapse = ", "), "\n")
Chromosomes involved: 1, 10, 11, 14, 16, 17, 19, 4, 6, 7, 8, X 
> cat("Median samples per junction:", median(neo$samples_expressed), "\n")
Median samples per junction: 41.5 
> cat("Range of read support:", min(neo$total_reads), "-", max(neo$total_reads), "\n")
Range of read support: 517 - 9233 

