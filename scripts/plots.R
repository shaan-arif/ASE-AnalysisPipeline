library(data.table)
library(ggplot2)

combined_file <- "output/ase/MBASED_results/gene_level_ASE_MBASED_all_samples.tsv"
dt <- fread(combined_file)
ggplot(dt, aes(x = majorAlleleProportion, fill = sample)) +
  geom_density(alpha = 0.3) +
  theme_bw() +
  labs(title = "Distribution of Major Allele Proportion per Sample",
       x = "Major Allele Proportion (MAF)",
       y = "Density")
ggplot(dt, aes(x = majorAlleleProportion, y = -log10(pValue))) +
  geom_point(alpha = 0.4) +
  theme_bw() +
  labs(title = "ASE Volcano Plot (All samples combined)",
       x = "Major Allele Proportion",
       y = "-log10(p-value)")
ggplot(dt, aes(x = majorAlleleProportion, y = -log10(pValue))) +
  geom_point(alpha = 0.4) +
  theme_bw() +
  facet_wrap(~ sample) +
  labs(title = "ASE Volcano Plots per Sample",
       x = "Major Allele Proportion",
       y = "-log10(p-value)")
library(reshape2)

wide <- dcast(dt, gene ~ sample, value.var = "majorAlleleProportion", fill = NA)
library(pheatmap)
library(data.table)
library(pheatmap)
library(dplyr)
library(tidyr)

# Load combined MBASED results (gene-level)
combined_file <- "output/ase/MBASED_results/gene_level_ASE_MBASED_all_samples.tsv"
combined <- fread(combined_file)

# Optional: filter significant ASE genes (p < 0.05 & majorAlleleProportion away from 0.5)
combined_sig <- combined %>%
  filter(pValue < 0.05 & abs(majorAlleleProportion - 0.5) > 0.1)

# Pivot to wide format: genes x samples
wide <- combined_sig %>%
  select(sample, gene, majorAlleleProportion) %>%
  pivot_wider(names_from = sample, values_from = majorAlleleProportion)

# Convert to matrix for pheatmap
mat <- as.matrix(wide[,-1])      # remove 'gene' column
rownames(mat) <- wide$gene

# Fill NA values with row mean
mat_filled <- t(apply(mat, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}))

# Select top 50 genes by median majorAlleleProportion across samples
top_genes <- names(sort(apply(mat_filled, 1, median, na.rm = TRUE), decreasing = TRUE))[1:50]
mat_top <- mat_filled[top_genes, ]

# Plot heatmap
pheatmap(mat_top,
         scale = "none",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "Top 50 ASE Genes (Major Allele Proportion)",
         fontsize_row = 6,
         fontsize_col = 8,
         color = colorRampPalette(c("skyblue", "white", "pink"))(100))



########################################################33


sig <- dt[majorAlleleProportion > 0.7 & pValue < 0.05]

ggplot(sig, aes(x = sample, fill = sample)) +
  geom_bar() +
  theme_bw() +
  labs(title = "Number of Significant ASE Genes per Sample",
       x = "Sample", y = "Count")
top <- dt[order(pValue)][1:30]   # top 30 genes

ggplot(top, aes(x = reorder(gene, -log10(pValue)),
                y = -log10(pValue),
                fill = sample)) +
  geom_bar(stat="identity") +
  coord_flip() +
  theme_bw() +
  labs(title = "Top ASE Genes (across samples)",
       x = "Gene", y = "-log10(pValue)")
####################################################################

library(data.table)
library(ggplot2)

# Load combined MBASED results
combined_file <- "output/ase/MBASED_results/gene_level_ASE_MBASED_all_samples.tsv"
combined <- fread(combined_file)

# Optional: filter significant ASE genes (p < 0.05 & majorAlleleProportion away from 0.5)
combined_sig <- combined[pValue < 0.05 & abs(majorAlleleProportion - 0.5) > 0.1]

# Violin + boxplot per sample
ggplot(combined_sig, aes(x = sample, y = majorAlleleProportion)) +
  geom_violin(fill = "maroon", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Distribution of Major Allele Proportion per Sample",
       x = "Sample",
       y = "Major Allele Proportion")


#####################################################
library(pheatmap)
library(reshape2)

# Prepare matrix: rows = genes, columns = samples
wide <- dcast(combined_sig, gene ~ sample, value.var = "majorAlleleProportion", fill = NA)
mat <- as.matrix(wide[, -1, with = FALSE])
rownames(mat) <- wide$gene

# Replace NAs with row mean
mat <- t(apply(mat, 1, function(x) { x[is.na(x)] <- mean(x, na.rm=TRUE); x }))

pheatmap(mat, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = FALSE,
         color = colorRampPalette(c("blue","white","red"))(100),
         main = "Heatmap of Significant ASE Genes")
#########################################
#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(pheatmap)
library(reshape2)

# ================================
# Paths
# ================================
ase_dir <- "output/ase"
mba_dir <- "output/ase/MBASED_results"
combined_file <- file.path(mba_dir, "gene_level_ASE_MBASED_all_samples.tsv")
map_csv <- "config/patient_sample_map.csv"
gene_map_dir <- "data/annotation/gene_mapping"
out_dir <- file.path(mba_dir, "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ================================
# Load data
# ================================
combined <- fread(combined_file)  # gene-level MBASED
ase_all <- fread(file.path(ase_dir, "combined/ase_all_samples.tsv"))  # variant-level ASE
map_samples <- fread(map_csv)  # Patient → WES_ID

# Merge ASE table with WES_ID
ase_all <- merge(
  ase_all,
  map_samples[, .(Patient, WES_ID)],
  by.x = "sample",
  by.y = "Patient",
  all.x = TRUE
)

# Load SNP → gene maps
map_files <- list.files(gene_map_dir, pattern = "*.snp_gene_map.tsv$", full.names = TRUE)
snp_maps <- rbindlist(lapply(map_files, function(f) {
  dt <- fread(f, header = FALSE, col.names = c("chr","pos","gene"))
  dt[, WES_ID := gsub(".snp_gene_map.tsv","",basename(f))]
  return(dt)
}))

# Merge ASE variants with gene mapping
ase_all[, pos := as.integer(position)]
ase_gene <- merge(
  ase_all,
  snp_maps,
  by.x = c("WES_ID","contig","pos"),
  by.y = c("WES_ID","chr","pos"),
  all.x = TRUE
)

# ================================
# Calculate allele ratios
# ================================
ase_gene[, majorAlleleProportion := refCount / (refCount + altCount)]
ase_gene[, abs_diff := abs(majorAlleleProportion - 0.5)]

# Filter significant ASE
ase_sig <- ase_gene[totalCount >= 10 & pmin(refCount, altCount) >= 3 & abs_diff > 0.1]

# ================================
# Manhattan-style plot
# ================================
pdf(file.path(out_dir,"ASE_manhattan_plot.pdf"), width = 12, height = 6)
ggplot(ase_sig, aes(x = pos, y = abs_diff, color = contig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  facet_wrap(~contig, scales = "free_x", nrow = 2) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(title = "Genome-wide ASE imbalance", x = "Genomic Position", y = "|Major Allele Freq - 0.5|")
dev.off()

# ================================
# Violin plot per gene (top genes with ASE)
# ================================
# Select top 20 genes with highest mean abs_diff
top_genes <- ase_sig[, .(mean_abs = mean(abs_diff, na.rm=TRUE)), by = gene][order(-mean_abs)][1:20, gene]
ase_top <- ase_sig[gene %in% top_genes]

pdf(file.path(out_dir,"ASE_violin_top_genes.pdf"), width = 10, height = 6)
ggplot(ase_top, aes(x = gene, y = abs_diff)) +
  geom_violin(fill = "skyblue", alpha = 0.7) +
  geom_jitter(height = 0, width = 0.2, size = 1, alpha = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Top Genes ASE Distribution", y = "|Major Allele Freq - 0.5|", x = "Gene")
dev.off()

# ================================
# Heatmap of ASE per gene x sample
# ================================
# Compute mean ASE per gene per sample
heat_data <- ase_sig[, .(mean_abs = mean(abs_diff, na.rm=TRUE)), by = .(gene, sample)]

# Cast to wide format
wide_heat <- dcast(heat_data, gene ~ sample, value.var = "mean_abs", fill = 0)

# Convert to matrix for pheatmap
mat <- as.matrix(wide_heat[, -1, drop = FALSE])
rownames(mat) <- wide_heat$gene

# Plot heatmap
pdf(file.path(out_dir,"ASE_heatmap.pdf"), width = 10, height = 8)
pheatmap(mat,
         scale = "none",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "ASE Heatmap (|Major Allele Freq - 0.5|)")
dev.off()


