# ================================================================
# MINUTE_3 — Heatmap of significant peaks + significant-region BED export
# Runs standalone after MINUTE_1, or in sequence via run_MINUTE.R
# ================================================================
source("config.R")

# DESeq2 + annotation results produced by MINUTE_1
annotated_results <- readRDS(annotated_rds)


###### HEATMAPS

###### ONE LARGE HEATMAP OF ALL SIG CHANGES

# Combine all annotated results
combined_sig <- do.call(rbind, lapply(names(annotated_results), function(mark) {
  df <- annotated_results[[mark]]
  df$ChIP <- mark
  df$peak_id <- paste0(df$chr, ":", df$start, "-", df$end)
  df
}))


# Keep only significant peaks (thresholds defined in config.R)
combined_sig <- combined_sig[is_significant(combined_sig), ]

combined_sig <- combined_sig[unique(combined_sig$peak_id),]

# Build GRanges
sig_gr <- GRanges(
  seqnames = combined_sig$chr,
  ranges = IRanges(combined_sig$start, combined_sig$end),
  strand = "*",
  peak_id = combined_sig$peak_id,
  ChIP = combined_sig$ChIP
)

# Standardize seqnames
seqlevelsStyle(sig_gr) <- "NCBI"
sig_gr <- keepStandardChromosomes(sig_gr, pruning.mode = "coarse")


library(GenomicRanges)
library(parallel)

# Create master bigwig list across all ChIPs
all_bigwigs <- unlist(lapply(chips, function(prefix) {
  get_bigwig_files(prefix, bigwigDir)
}))

# Metadata for columns (samples)
col_annot <- data.frame(
  Sample = names(all_bigwigs),
  File = unname(all_bigwigs),
  stringsAsFactors = FALSE
)
col_annot$ChIP <- rep(names(chips), each = 11)  # 11 samples per ChIP
col_annot$Genotype <- ifelse(grepl("HP1gKO", col_annot$Sample), "HP1gKO", "WT")
col_annot$Replicate <- sub("^[^.]+\\.", "", col_annot$Sample)



# Extract signal
library(parallel)

signal_matrix <- mclapply(seq_len(nrow(col_annot)), function(i) {
  bw <- col_annot$File[i]
  if (!file.exists(bw)) return(rep(NA_real_, length(sig_gr)))
  
  vals <- tryCatch(
    bw_loci(bw, sig_gr),
    error = function(e) {
      message(sprintf("❌ Failed on %s: %s", bw, e$message))
      return(rep(NA_real_, length(sig_gr)))
    }
  )
  
  mcols(vals)[[1]]
}, mc.cores = detectCores())


# Format into matrix
signal_mat <- do.call(cbind, signal_matrix)
colnames(signal_mat) <- col_annot$Sample
rownames(signal_mat) <- sig_gr$peak_id
signal_mat <- as.matrix(signal_mat)

sig_df <- as.data.frame(mcols(sig_gr))
sig_df <- sig_df[match(rownames(signal_mat), sig_df$peak_id), ]

# Join with combined_sig to get annotations
sig_df <- left_join(sig_df, combined_sig, by = c("peak_id","ChIP"))
rownames(sig_df) <- sig_df$peak_id

#sig_df$chr <- factor(sig_df$chr, levels = c(paste(seq(1,19, by=1),sep = ","), "X","Y"))
# calculate domain size
sig_df$peak_size_kb <- (sig_df$end - sig_df$start) / 1000  # Convert to kb

# COL ANNOTATIONS
# Fix replicate names to standard form
col_annot <- droplevels(col_annot)

# Ensure all columns are properly factored
col_annot$ChIP <- factor(col_annot$ChIP, levels = c("H3K4me3","H3K36me3","H3K9me2","H3K9me3","H4K20me3"))
col_annot$Genotype <- factor(col_annot$Genotype, levels = c("WT", "HP1gKO"))
col_annot$Replicate <- factor(col_annot$Replicate, levels = unique(col_annot$Replicate))

# Recalculate color palettes
chip_colors <- structure(scales::hue_pal()(length(levels(col_annot$ChIP))), names = levels(col_annot$ChIP))
geno_colors <- c(WT = "#1f78b4", HP1gKO = "#e31a1c")

# Column annotation
col_ha <- HeatmapAnnotation(
  ChIP = col_annot$ChIP,
  Genotype = col_annot$Genotype,
  Replicate = col_annot$Replicate,
  col = list(
    ChIP = chip_colors,
    Genotype = geno_colors
  ),
  show_annotation_name = TRUE,
  annotation_name_side = "right"
)


# ROW ANNOTATIONS

# Defensive handling in case of missing or empty annotation
unique_regions <- unique(na.omit(sig_df$genomic_region))
region_colors <- if (length(unique_regions) > 0) {
  structure(scales::hue_pal()(length(unique_regions)), names = unique_regions)
} else {
  NULL
}

# Get unique chromosomes (safely handle NAs)
chrom_levels <- unique(na.omit(sig_df$chr))

# Assign colors using hue palette
chrom_colors <- structure(
  scales::hue_pal()(length(chrom_levels)),
  names = chrom_levels
)

# domain size
# Clip extreme peak sizes for color scaling (but keep full values in annotation)
sig_df$peak_size_kb_capped <- pmin(sig_df$peak_size_kb, 2)

# Define color gradient for peak size (adjust color/scale if needed)
size_col_fun <- circlize::colorRamp2(
  c(0.3, 1, 6),   # Adjust based on your distribution
  c("lightyellow", "orange", "red")
)


# Rebuild row annotation
row_ha <- rowAnnotation(
  Chromosome = sig_df$chr,
  Region = sig_df$genomic_region,
  Repeat = sig_df$repeat_class,
  TAD_Boundary = sig_df$overlaps_with_Tad_boundary,
  PeakSize_kb = sig_df$peak_size_kb_capped,
  col = list(
    Chromosome = chrom_colors,
    Region = region_colors,
    Repeat = c(LINE = "#7570b3", SINE = "#d95f02", LTR = "#1b9e77", complex = "#e7298a", none = "grey80"),
    TAD_Boundary = c(`TRUE` = "firebrick", `FALSE` = "white"),
    PeakSize_kb = size_col_fun
  ),
  show_annotation_name = TRUE,
  annotation_name_side = "top"
)

signal_scaled <- t(scale(t(signal_mat)))  # row Z-score

# ---------------------------
# 4. Plot
# ---------------------------

ht <- Heatmap(
  signal_scaled,
  name = "Z-score",
  top_annotation = col_ha,
  left_annotation = row_ha,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  row_km = 5,
  #row_split = sig_df$chr, 
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_split = col_annot$ChIP,
  column_title_gp = gpar(fontsize = 10),
  row_title = "Significant Peaks",
  heatmap_legend_param = list(title = "Signal (Z-score)")
)

#draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")

pdf("2000maxgap_indsignificance_with_TAD.pdf", width = 10, height = 16)  # adjust size as needed
ht_drawn <- draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()


# ---------------------------
# 5. Extract k-means clusters and export significant-region BED files
#     (rehomed from Minute_count_and_annotate_from_scaled_bw_v2.R)
# ---------------------------

# Read the row k-means assignment back off the heatmap we just drew
main_ht        <- ht_drawn@ht_list[[1]]
km_assignments <- ComplexHeatmap::row_order(main_ht)   # list of row indices, one element per cluster

cluster_vector <- rep(NA_integer_, nrow(signal_mat))
names(cluster_vector) <- rownames(signal_mat)
for (k in seq_along(km_assignments)) {
  cluster_vector[km_assignments[[k]]] <- k
}
sig_df$Cluster <- cluster_vector[rownames(sig_df)]

# --- Combined BED of all marks (score = log2FoldChange) ---
sig_gr_bed <- GRanges(
  seqnames = sig_df$chr,
  ranges   = IRanges(start = sig_df$start, end = sig_df$end),
  strand   = "*"
)
mcols(sig_gr_bed)$name    <- sig_df$peak_id
mcols(sig_gr_bed)$score   <- round(sig_df$log2FoldChange, 3)
mcols(sig_gr_bed)$ChIP    <- sig_df$ChIP
mcols(sig_gr_bed)$Cluster <- sig_df$Cluster
export(sig_gr_bed, "significant_peaks_with_metadata.bed", format = "BED")

# --- One BED per mark (NCBI seqnames; score = -log10(pvalue)) ---
bed_style <- "NCBI"   # switch to "UCSC" for chr-prefixed names
for (chip in unique(sig_df$ChIP)) {
  df_chip <- subset(sig_df, ChIP == chip)
  gr <- GRanges(
    seqnames = df_chip$chr,
    ranges   = IRanges(start = df_chip$start, end = df_chip$end),
    strand   = "*"
  )
  seqlevelsStyle(gr) <- bed_style
  mcols(gr)$name    <- df_chip$peak_id
  mcols(gr)$score   <- -log10(df_chip$pvalue)
  mcols(gr)$cluster <- df_chip$Cluster
  mcols(gr)$chip    <- df_chip$ChIP
  out_file <- paste0("significant_peaks_", chip, ".bed")
  export(gr, con = out_file, format = "BED")
  message("Exported: ", out_file)
}
