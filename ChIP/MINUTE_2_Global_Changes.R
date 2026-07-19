# ================================================================
# MINUTE_2 - Global changes: the genome-wide picture of the HP1gKO response
# ----------------------------------------------------------------
# Everything about the OVERALL change, in one place:
#   * main clustered heatmap of significant peaks (k-means, seeded)
#   * per-chromosome change plots + per-mark log2FC distributions (global-shift
#     diagnostic)
#   * relationships between marks (H3K9me3 vs H4K20me3) and size x signal x log2FC
#   * k-means cluster CHARACTERISATION (signal/loss profiles, ChromHMM coverage,
#     repeat/region composition, per-cluster ChromHMM enrichment)
#   * TAD-boundary + ChromHMM-state enrichment of the significant set
# Repeat enrichment lives in MINUTE_4_Repeats; the differential H4K20me3-vs-
# H3K9me3 split lives in MINUTE_3_Differential_loss.
#
# Persists results/rds/cluster_analysis_inputs.rds (used within this stage for the
# cluster characterisation) and significant_peaks_clusters.tsv.
# geno_colors / theme_m / mark_levels + hyper_row() come from config.R.
# Runs after MINUTE_1 (needs the annotated .rds).
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

# Column metadata (samples) — straight from the sample sheet, all marks
col_annot <- data.frame(
  Sample    = samples$sample_id,
  File      = samples$bigwig_path,
  ChIP      = as.character(samples$mark),
  Genotype  = as.character(samples$genotype),
  Replicate = samples$replicate,
  stringsAsFactors = FALSE
)



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

set.seed(42)   # reproducible k-means cluster assignment/numbering across runs
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

# Draw once to capture ht_drawn (used for cluster extraction below), then re-draw
# the SAME laid-out object for the other format so png/pdf share one clustering.
hm_dir <- file.path(fig_dir, "heatmap")
if (!dir.exists(hm_dir)) dir.create(hm_dir, recursive = TRUE, showWarnings = FALSE)
hm_base <- "2000maxgap_indsignificance_with_TAD"
pdf(file.path(hm_dir, paste0(hm_base, ".pdf")), width = 10, height = 16)
ht_drawn <- draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
png(file.path(hm_dir, paste0(hm_base, ".png")), width = 10, height = 16, units = "in", res = 200)
draw(ht_drawn); dev.off()
message("Saved: heatmap/", hm_base, ".{png,pdf}")


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

# --- Persist cluster table + signal matrix (self-consumed below for the cluster --
# rtracklayer BED export keeps only standard columns (it DROPS Cluster/ChIP), so
# the metadata BED below is browser-only. Save the real cluster table + the
# per-sample signal matrix here so the cluster characterisation runs without re-reading.
cluster_cols <- intersect(
  c("peak_id", "chr", "start", "end", "ChIP", "Cluster", "log2FoldChange",
    "pvalue", "baseMean", "peak_size_kb", "genomic_region", "repeat_class",
    "repeat_family", "SYMBOL", "label", "overlaps_with_Tad_boundary",
    "chromHMM_state", "chromHMM_purity"),
  names(sig_df))
saveRDS(list(sig_df = sig_df, signal_mat = signal_mat, col_annot = col_annot),
        file.path(rds_dir, "cluster_analysis_inputs.rds"))
fwrite(sig_df[, cluster_cols],
       file.path(tables_dir, "significant_peaks_clusters.tsv"), sep = "\t")
message("Saved: cluster_analysis_inputs.rds + significant_peaks_clusters.tsv")

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
export(sig_gr_bed, file.path(bed_dir, "significant_peaks_with_metadata.bed"), format = "BED")

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
  out_file <- file.path(bed_dir, paste0("significant_peaks_", chip, ".bed"))
  export(gr, con = out_file, format = "BED")
  message("Exported: ", out_file)
}


# ================================================================
# 6. Per-chromosome "changes" plots
#    Domain size (kb, log10) vs log2FoldChange, faceted by chromosome,
#    point size = ChIP signal (baseMean), colour = genomic region or repeat.
#    Highlighted peaks are coloured; the rest are grey. Ported from
#    Minute_count_and_annotate_from_scaled_bw_v2.R — operates purely on
#    annotated_results, so no bigWigs are re-read here.
#    Deps (ggplot2, ggrepel, dplyr) come from config.R.
#
#    HIGHLIGHT CRITERION IS PER-MARK (plot_criterion below):
#      "pvalue" = is_significant() (|log2FC|>0.5 AND p<thr). Correct for marks
#                 whose log2FC distribution is CENTRED at 0 (per-peak testing
#                 answers the right question): H3K4me3, H3K36me3, H3K9me2.
#      "effect" = |log2FC| > effect_lfc only. For marks with a GLOBAL directional
#                 shift (median log2FC far below 0) where per-peak p is
#                 underpowered — the change is genome-wide, not a testable
#                 minority, and counts are tiny (baseMean ~1-3): H3K9me3,
#                 H4K20me3. See the log2FC-distribution diagnostic (section 7).
#    NOTE: this criterion governs ONLY these plots. The heatmap, BED exports and
#    repeat hypergeometric test (MINUTE_4) still use is_significant() (p-based) throughout.
# ================================================================

# --- Appearance / export knobs -----------------------------------
# Two independent levers:
#  * TEXT/POINT size is set by the canvas size (ggplot sizes are absolute mm/pt,
#    so a BIGGER canvas makes them look SMALLER). If text/points are too big,
#    raise plot_w/plot_h; too small, lower them. Fine-tune with label_size /
#    point_range / base_size.
#  * PANEL SHAPE is set by panel_aspect (height/width per chromosome facet),
#    independent of canvas size: 1 = square, <1 = wider, >1 = taller.
plot_w       <- 14           # export width  (inches)
plot_h       <- 7            # export height (inches)
plot_dpi     <- 300
panel_aspect <- 1            # per-facet aspect ratio (1 = square panels)
base_size    <- 11           # theme base font size
label_size   <- 2.5          # gene-label text size
point_range  <- c(0, 6)      # min/max point size — as in the original script
size_var     <- "peak_size_kb"  # variable mapped to point SIZE:
                             #  "peak_size_kb" = domain size — reproduces the
                             #    original plots (large domains render large);
                             #    note this duplicates the x-axis.
                             #  "baseMean" = ChIP signal (independent info, but
                             #    H4K20me3 signal is uniformly low ~1-4, so points
                             #    look uniformly small — cor(domain,baseMean)~0).

# --- Highlight criterion per mark (see section header) -----------
# "pvalue" (centred marks) or "effect" (globally-shifted marks). Flip a mark
# here to change how its change-plot is coloured; nothing else is affected.
plot_criterion <- list(
  H3K4me3  = "pvalue",
  H3K36me3 = "pvalue",
  H3K9me2  = "pvalue",   # borderline (median log2FC -0.10); flip to "effect" if desired
  H3K9me3  = "effect",
  H4K20me3 = "effect"
)
effect_lfc <- 0.5            # |log2FC| cutoff used in "effect" mode

# Chromosome facet order (autosomes 1-19 then X, Y; drops scaffolds)
chr_levels <- c(as.character(1:19), "X", "Y")

# Genes to always label when highlighted (curated in the original script)
keep_labels <- c("Eif1ad2", "Scgb2b13-ps", "Scgb1b29", "Scgb1b7", "Zfp980",
                 "Zfp781b", "Gm6756", "Muc17", "Slc16a14", "Hnf4g", "Cyp3a11",
                 "Gm31138", "Vmn1r238", "Eif1ad4", "Tlr4", "Vmn2r50", "Vmn2r44",
                 "Vmn2r40", "Pira1", "Mrgprb8")

# Reference genes: ALWAYS labelled (open ring + italic navy label) even when NOT
# highlighted, so they're shown for context without implying significance.
# Cbx3 encodes HP1gamma — the knockout target — so it belongs on every plot.
ref_labels <- c("Cbx3")

# Repeat-class palette (matches the heatmap row annotation above)
repeat_palette <- c(LINE = "#7570b3", SINE = "#d95f02", LTR = "#1b9e77",
                    complex = "#e7298a", none = "grey80")

# Title: "<mark> Changes in HP1gamma^fl/fl Emx1^Cre Cortex"
changes_title <- function(mark) {
  bquote(.(mark) * " Changes in " *
           italic("HP1" * gamma^"fl/fl" * " Emx1"^"Cre") * " Cortex")
}

# Which peaks are highlighted (coloured), given the per-mark criterion.
highlight_peaks <- function(df, mark, criterion) {
  if (criterion == "effect") {
    !is.na(df$log2FoldChange) & abs(df$log2FoldChange) > effect_lfc
  } else {
    df$ChIP <- mark
    is_significant(df)          # |log2FC|>0.5 AND p<thr
  }
}

# Honest caption describing what "coloured" means for this mark.
highlight_caption <- function(mark, criterion) {
  if (criterion == "effect") {
    bquote("coloured: |log2FC| >" ~ .(effect_lfc) ~
             "(effect size — global shift; per-peak p underpowered)")
  } else {
    p <- if (!is.null(sig_thresholds[[mark]])) sig_thresholds[[mark]]$p else NA
    bquote("coloured: significant (|log2FC| > 0.5 and p <" ~ .(p) * ")")
  }
}

# Build one faceted plot for a mark, coloured by `color_by`.
changes_by_chr_plot <- function(df, mark, color_by = c("genomic_region", "repeat_class"),
                                criterion = c("pvalue", "effect")) {
  color_by  <- match.arg(color_by)
  criterion <- match.arg(criterion)

  df$ChIP         <- mark
  df$highlight    <- highlight_peaks(df, mark, criterion)
  df$peak_size_kb <- (df$end - df$start) / 1000
  df$size_val     <- df[[size_var]]                  # chosen size variable
  df <- df[df$chr %in% chr_levels, ]
  df$chr <- factor(df$chr, levels = chr_levels)

  # Reference genes (e.g. Cbx3): always shown, regardless of highlight status
  ref_df <- df %>%
    filter(label %in% ref_labels) %>%
    arrange(desc(baseMean)) %>%
    group_by(label) %>% slice_head(n = 1) %>% ungroup()

  # Label set: top 25 highlighted by domain size, plus the curated keep_labels.
  # Reference genes are excluded here — they get their own distinct styling.
  base_set <- df %>%
    filter(highlight, baseMean > 2, peak_size_kb > 3) %>%
    arrange(desc(peak_size_kb)) %>%
    slice_head(n = 25)
  forced_set <- df %>%
    filter(highlight, label %in% keep_labels) %>%
    arrange(desc(peak_size_kb)) %>%
    group_by(label) %>% slice_head(n = 1) %>% ungroup()
  lab_df <- bind_rows(base_set, forced_set) %>%
    distinct(peak_id, .keep_all = TRUE) %>%
    filter(!(label %in% ref_labels))

  # Size legend: honest label for the chosen variable
  size_breaks <- c(1, 10, 100, 1000)
  size_name   <- if (size_var == "peak_size_kb") "Domain Size (kb)" else "ChIP Signal baseMean"

  p <- ggplot(df, aes(x = peak_size_kb, y = log2FoldChange)) +
    # not highlighted: faint grey
    geom_point(data = subset(df, !highlight),
               aes(size = size_val), colour = "grey70", alpha = 0.15) +
    # highlighted: coloured
    geom_point(data = subset(df, highlight),
               aes(size = size_val, colour = .data[[color_by]]), alpha = 0.6) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_text_repel(data = lab_df, aes(label = label),
                    size = label_size, max.overlaps = Inf, box.padding = 0.15,
                    point.padding = 0.1, min.segment.length = 0,
                    segment.size = 0.2, segment.alpha = 0.4, colour = "black") +
    # reference genes: open ring + italic navy label (not implying significance)
    geom_point(data = ref_df, aes(size = size_val),
               shape = 21, fill = NA, colour = "navy", stroke = 0.6) +
    geom_text_repel(data = ref_df, aes(label = label),
                    size = label_size, fontface = "italic", colour = "navy",
                    box.padding = 0.4, min.segment.length = 0,
                    segment.colour = "navy", segment.size = 0.3) +
    scale_x_log10() +
    # Linear size scale, range/breaks as in the original script.
    scale_size_continuous(range = point_range, breaks = size_breaks,
                          name = size_name) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3)),
           size   = guide_legend(override.aes = list(colour = "grey30", alpha = 1))) +
    labs(title = changes_title(mark), subtitle = highlight_caption(mark, criterion),
         x = "Domain Size (kb)", y = "log2FoldChange", colour = color_by,
         caption = if (nrow(ref_df)) "navy ring = reference gene (shown regardless of significance)" else NULL) +
    facet_wrap(~chr, ncol = 7) +
    theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          aspect.ratio = panel_aspect)

  if (color_by == "repeat_class") {
    p <- p + scale_colour_manual(values = repeat_palette)
  }
  p
}

# Generate both colour variants for every mark, using its per-mark criterion
for (mark in names(annotated_results)) {
  crit <- plot_criterion[[mark]]
  if (is.null(crit)) crit <- "pvalue"          # default for any unlisted mark
  cat(sprintf("\n=== Change plots for: %s (criterion: %s) ===\n", mark, crit))
  df <- annotated_results[[mark]]
  for (cb in c("genomic_region", "repeat_class")) {
    suffix <- if (cb == "genomic_region") "coloured_by_genomic_region" else "coloured_by_repeat"
    save_fig(changes_by_chr_plot(df, mark, color_by = cb, criterion = crit),
             paste0(mark, "_changes_by_chr_", suffix), "change_plots",
             width = plot_w, height = plot_h)
  }
}


# ================================================================
# 7. Global-shift diagnostic
#    Per-mark density of per-peak log2FC (HP1gKO vs WT). This is the evidence
#    behind the per-mark criterion choice above: a distribution CENTRED on 0
#    means per-peak testing is valid ("pvalue"); a distribution shifted well
#    below 0 is a genome-wide change that per-peak p can't capture ("effect").
#    Because DESeq2 runs with sizeFactors = 1 (bigWigs pre-scaled), a real
#    global shift is preserved here rather than normalised away.
# ================================================================

diag_df <- do.call(rbind, lapply(names(annotated_results), function(m) {
  data.frame(ChIP = m, log2FoldChange = annotated_results[[m]]$log2FoldChange)
}))
diag_df <- diag_df[is.finite(diag_df$log2FoldChange), ]
diag_df$ChIP <- factor(diag_df$ChIP, levels = names(annotated_results))

# Per-mark median + % of domains reduced, annotated on each facet
diag_summary <- do.call(rbind, lapply(levels(diag_df$ChIP), function(m) {
  l <- diag_df$log2FoldChange[diag_df$ChIP == m]
  data.frame(ChIP = m, med = median(l), pct_down = 100 * mean(l < 0))
}))
diag_summary$ChIP <- factor(diag_summary$ChIP, levels = names(annotated_results))
diag_summary$lab  <- sprintf("median %.2f | %.0f%% down", diag_summary$med, diag_summary$pct_down)

g_diag <- ggplot(diag_df, aes(log2FoldChange)) +
  geom_density(aes(fill = ChIP), colour = NA, alpha = 0.6) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_vline(data = diag_summary, aes(xintercept = med),
             linetype = "dashed", colour = "red", linewidth = 0.4) +
  geom_text(data = diag_summary, aes(x = -Inf, y = Inf, label = lab),
            hjust = -0.05, vjust = 1.4, size = 3, colour = "red") +
  facet_wrap(~ChIP, ncol = 1, scales = "free_y") +
  coord_cartesian(xlim = c(-2, 2)) +
  labs(title = "Per-peak log2FC distribution by mark (HP1gKO vs WT)",
       subtitle = "solid = 0, dashed red = median. Median far below 0 = global loss (per-peak p underpowered → effect-size framing)",
       x = "log2FoldChange", y = "density") +
  theme_minimal(base_size = base_size) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

save_fig(g_diag, "log2FC_distribution_by_mark", "change_plots", width = 8, height = 11)


# ================================================================
# 8. Relationship: domain size x ChIP signal (baseMean) x log2FC  [TOTAL]
#    Genome-wide ("total") across all peaks; the per-cluster counterparts
#    counterparts. The change-plots use x=size and size=size (redundant). These
#    three views put domain size, baseMean and log2FC on separate channels to
#    reveal how the change depends jointly on domain length and signal intensity.
#    (baseMean is length-normalised mean coverage, so size and baseMean are
#    ~independent axes: cor ~ 0.)
#    All peaks are used (no highlight filter) — binning/averaging also tames
#    the per-peak noise from the tiny counts.
# ================================================================

rel_df <- do.call(rbind, lapply(names(annotated_results), function(m) {
  d <- annotated_results[[m]]
  data.frame(ChIP = m, kb = (d$end - d$start) / 1000,
             baseMean = d$baseMean, log2FoldChange = d$log2FoldChange)
}))
rel_df <- rel_df[is.finite(rel_df$log2FoldChange) & is.finite(rel_df$baseMean) &
                 rel_df$kb > 0 & rel_df$baseMean > 0, ]
rel_df$ChIP <- factor(rel_df$ChIP, levels = names(annotated_results))

# --- 8a. Binned heatmap: mean log2FC over the size x baseMean plane ----------
# Rectangular bins (stat_summary_2d, core ggplot2 — no hexbin dependency).
g_rel_heat <- ggplot(rel_df, aes(kb, baseMean, z = log2FoldChange)) +
  stat_summary_2d(fun = mean, bins = 30) +
  scale_x_log10() + scale_y_log10() +
  scale_fill_gradient2(low = "#2166ac", mid = "grey95", high = "#b2182b",
                       midpoint = 0, limits = c(-1, 1),
                       oob = scales::squish, name = "mean\nlog2FC") +
  facet_wrap(~ChIP, ncol = 3) +
  labs(title = "log2FC across domain size x ChIP signal (mean per bin)",
       subtitle = "blue = loss in HP1gKO; fill clipped to [-1, 1]",
       x = "Domain Size (kb)", y = "baseMean (ChIP signal)") +
  theme_minimal(base_size = base_size) +
  theme(panel.grid.minor = element_blank())
save_fig(g_rel_heat, "relationship_total_size_baseMean_log2FC_binned", "relationships",
         width = 13, height = 8)

# --- 8b. Scatter: log2FC vs domain size, coloured by baseMean ----------------
g_rel_scatter <- ggplot(rel_df, aes(kb, log2FoldChange, colour = baseMean)) +
  geom_point(size = 0.35, alpha = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_x_log10() +
  scale_colour_viridis_c(trans = "log10", name = "baseMean") +
  facet_wrap(~ChIP, ncol = 3) +
  labs(title = "log2FC vs domain size, coloured by ChIP signal (baseMean)",
       x = "Domain Size (kb)", y = "log2FoldChange") +
  theme_minimal(base_size = base_size) +
  theme(panel.grid.minor = element_blank())
save_fig(g_rel_scatter, "relationship_total_log2FC_vs_size_by_baseMean", "relationships", width = 13, height = 8)

# --- 8c. Interaction: mean log2FC vs size bin, stratified by baseMean tertile -
# Direct visual of the size x baseMean interaction: a line per baseMean tertile.
# Lines that diverge/cross = the effect of size on log2FC depends on signal.
size_breaks_kb <- c(0, 1, 2, 5, 10, 20, 50, 100, 500, Inf)
inter <- rel_df %>%
  group_by(ChIP) %>%
  mutate(base_tertile = factor(ntile(baseMean, 3), labels = c("low", "mid", "high")),
         size_bin     = cut(kb, size_breaks_kb, dig.lab = 4)) %>%
  ungroup() %>%
  filter(!is.na(size_bin)) %>%
  group_by(ChIP, base_tertile, size_bin) %>%
  summarise(mean_l2fc = mean(log2FoldChange), n = n(), .groups = "drop") %>%
  filter(n >= 20)                                   # drop sparse bins

g_inter <- ggplot(inter, aes(size_bin, mean_l2fc, colour = base_tertile, group = base_tertile)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_line() + geom_point(size = 1) +
  scale_colour_viridis_d(option = "C", end = 0.85, name = "baseMean\ntertile") +
  facet_wrap(~ChIP, ncol = 3) +
  labs(title = "Interaction: effect of domain size on log2FC, by ChIP-signal tertile",
       subtitle = "diverging/crossing lines = size effect depends on signal strength (bins with n < 20 dropped)",
       x = "Domain size bin (kb)", y = "mean log2FoldChange") +
  theme_minimal(base_size = base_size) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = base_size - 3))
save_fig(g_inter, "relationship_total_interaction_size_x_baseMean", "relationships", width = 13, height = 8)


# ================================================================
# 9. Relationship BY K-MEANS CLUSTER (size x baseMean x log2FC)
#    The same three views as Section 8, but faceted by cluster (using sig_df,
#    the clustered significant peaks) instead of genome-wide. Complements the
#    _total versions; shows whether the size/signal/change structure differs
#    between the chromatin-signature clusters.
# ================================================================

clus_df <- sig_df[is.finite(sig_df$log2FoldChange) & is.finite(sig_df$baseMean) &
                  sig_df$peak_size_kb > 0 & sig_df$baseMean > 0 & !is.na(sig_df$Cluster), ]
clus_df$Cluster <- factor(paste("Cluster", clus_df$Cluster),
                          levels = paste("Cluster", sort(unique(sig_df$Cluster))))

# --- 9a. Binned heatmap per cluster ---
g_c_heat <- ggplot(clus_df, aes(peak_size_kb, baseMean, z = log2FoldChange)) +
  stat_summary_2d(fun = mean, bins = 25) +
  scale_x_log10() + scale_y_log10() +
  scale_fill_gradient2(low = "#2166ac", mid = "grey95", high = "#b2182b",
                       midpoint = 0, limits = c(-1, 1), oob = scales::squish,
                       name = "mean\nlog2FC") +
  facet_wrap(~Cluster, ncol = 3) +
  labs(title = "log2FC across domain size x ChIP signal, by cluster",
       x = "Domain Size (kb)", y = "baseMean (ChIP signal)") +
  theme_minimal(base_size = base_size) +
  theme(panel.grid.minor = element_blank())
save_fig(g_c_heat, "relationship_bycluster_size_baseMean_log2FC_binned", "relationships", width = 12, height = 8)

# --- 9b. Scatter per cluster ---
g_c_scatter <- ggplot(clus_df, aes(peak_size_kb, log2FoldChange, colour = baseMean)) +
  geom_point(size = 0.5, alpha = 0.4) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_x_log10() +
  scale_colour_viridis_c(trans = "log10", name = "baseMean") +
  facet_wrap(~Cluster, ncol = 3) +
  labs(title = "log2FC vs domain size by cluster, coloured by baseMean",
       x = "Domain Size (kb)", y = "log2FoldChange") +
  theme_minimal(base_size = base_size) +
  theme(panel.grid.minor = element_blank())
save_fig(g_c_scatter, "relationship_bycluster_log2FC_vs_size_by_baseMean", "relationships", width = 12, height = 8)

# --- 9c. Interaction per cluster (baseMean tertiles within each cluster) ---
inter_c <- clus_df %>%
  group_by(Cluster) %>%
  mutate(base_tertile = factor(ntile(baseMean, 3), labels = c("low", "mid", "high")),
         size_bin     = cut(peak_size_kb, size_breaks_kb, dig.lab = 4)) %>%
  ungroup() %>%
  filter(!is.na(size_bin)) %>%
  group_by(Cluster, base_tertile, size_bin) %>%
  summarise(mean_l2fc = mean(log2FoldChange), n = n(), .groups = "drop") %>%
  filter(n >= 10)

g_c_inter <- ggplot(inter_c, aes(size_bin, mean_l2fc, colour = base_tertile, group = base_tertile)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_line() + geom_point(size = 1) +
  scale_colour_viridis_d(option = "C", end = 0.85, name = "baseMean\ntertile") +
  facet_wrap(~Cluster, ncol = 3) +
  labs(title = "Interaction: domain size x signal on log2FC, by cluster",
       subtitle = "bins with n < 10 dropped",
       x = "Domain size bin (kb)", y = "mean log2FoldChange") +
  theme_minimal(base_size = base_size) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = base_size - 3))
save_fig(g_c_inter, "relationship_bycluster_interaction_size_x_baseMean", "relationships", width = 12, height = 8)

# --- 9d/9e. H3K9me3 vs H4K20me3 CHANGE (DESeq2 log2FoldChange) per cluster ---
# H3K9me3 and H4K20me3 share the 2 kb-merged peak set, so their DESeq2
# log2FoldChange align by peak_id. Join both to the clustered peaks. 9d: are the
# two marks' modelled changes coordinated per peak? 9e: how is each mark's
# change distributed across clusters (uniform vs cluster-specific)?
l2fc_by_id <- function(mk) {
  d <- annotated_results[[mk]]
  setNames(d$log2FoldChange, paste0(d$chr, ":", d$start, "-", d$end))
}
k9  <- l2fc_by_id("H3K9me3")
k20 <- l2fc_by_id("H4K20me3")
clus_lv <- paste("Cluster", sort(unique(sig_df$Cluster)))
chg_df <- data.frame(
  Cluster      = factor(paste("Cluster", sig_df$Cluster), levels = clus_lv),
  H3K9me3      = unname(k9[sig_df$peak_id]),
  H4K20me3     = unname(k20[sig_df$peak_id]),
  peak_size_kb = sig_df$peak_size_kb)
d9 <- chg_df[is.finite(chg_df$H3K9me3) & is.finite(chg_df$H4K20me3) &
             is.finite(chg_df$peak_size_kb) & chg_df$peak_size_kb > 0 &
             !is.na(chg_df$Cluster), ]

# 9d: BINPLOT of the change plane, coloured by domain size.
# x = H3K9me3 log2FC, y = H4K20me3 log2FC, fill = median domain size per bin.
# Fixes the scatter overplotting and shows where the large domains sit in the
# change plane (e.g. big domains losing H4K20me3 but retaining H3K9me3?).
rlab <- d9 %>% group_by(Cluster) %>%
  summarise(r = cor(H3K9me3, H4K20me3, method = "spearman"), n = n(), .groups = "drop")
g_bin <- ggplot(d9, aes(H3K9me3, H4K20me3)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55") +
  stat_summary_2d(aes(z = peak_size_kb), fun = median, bins = 22) +
  scale_fill_viridis_c(trans = "log10", name = "median\ndomain\nsize (kb)") +
  geom_text(data = rlab, aes(x = Inf, y = Inf, label = sprintf("rho=%.2f (n=%d)", r, n)),
            inherit.aes = FALSE, hjust = 1.05, vjust = 1.5, size = 3, colour = "black") +
  facet_wrap(~Cluster, nrow = 1) +
  labs(title = "H3K9me3 vs H4K20me3 change (DESeq2 log2FC), binned by domain size",
       subtitle = "fill = median domain size of peaks in each bin; where do the large domains sit in the change plane?",
       x = "H3K9me3 log2FoldChange", y = "H4K20me3 log2FoldChange") +
  theme_minimal(base_size = base_size) + theme(panel.grid.minor = element_blank())
save_fig(g_bin, "relationship_bycluster_change_H3K9me3_vs_H4K20me3_binned", "relationships", width = 14, height = 3.6)

# 9e: DESeq2 log2FC distribution per cluster, one panel per mark (shared y — both
# are log2FoldChange, so magnitudes are directly comparable)
e9 <- melt(as.data.table(chg_df), id.vars = "Cluster",
           variable.name = "mark", value.name = "l2fc")
e9 <- e9[is.finite(l2fc) & !is.na(Cluster)]
g_chg <- ggplot(e9, aes(Cluster, l2fc, fill = Cluster)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3) +      # hide outliers; zoom below
  facet_wrap(~mark, nrow = 1) +
  scale_fill_viridis_d(option = "D", end = 0.9, guide = "none") +
  # zoom to the box/whisker range so the medians are readable (outliers span +-4)
  coord_cartesian(ylim = quantile(e9$l2fc, c(0.01, 0.99), na.rm = TRUE)) +
  labs(title = "H3K9me3 / H4K20me3 change (DESeq2 log2FoldChange) across clusters",
       subtitle = "uniform across clusters = constant effect; spread = cluster-specific change (y zoomed to 1-99%)",
       x = "Cluster", y = "log2FoldChange (HP1gKO vs WT)") +
  theme_minimal(base_size = base_size) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = base_size - 3))
save_fig(g_chg, "relationship_bycluster_mark_change_distribution", "relationships", width = 9, height = 4)

# 9f: ALL merged-set regions (both marks measured), significant coloured by
# cluster. Intersect the two marks' peak_ids for every region with both
# log2FCs; join sig_df$Cluster (NA where not significant). Shows where each
# cluster sits in the change plane relative to the whole population.
common  <- intersect(names(k9), names(k20))
all_pts <- data.frame(
  H3K9me3  = unname(k9[common]),
  H4K20me3 = unname(k20[common]),
  Cluster  = sig_df$Cluster[match(common, sig_df$peak_id)])
all_pts <- all_pts[is.finite(all_pts$H3K9me3) & is.finite(all_pts$H4K20me3), ]
all_pts$clustered <- !is.na(all_pts$Cluster)

g_all <- ggplot(all_pts, aes(H3K9me3, H4K20me3)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_point(data = subset(all_pts, !clustered), colour = "grey82",
             size = 0.25, alpha = 0.3) +
  geom_point(data = subset(all_pts, clustered), aes(colour = factor(Cluster)),
             size = 0.5, alpha = 0.6) +
  scale_colour_viridis_d(option = "D", end = 0.9, name = "Cluster") +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(title = "H3K9me3 vs H4K20me3 change: all regions, significant coloured by cluster",
       subtitle = sprintf("grey = all %d merged-set regions with both marks; coloured = %d significant/clustered",
                           nrow(all_pts), sum(all_pts$clustered)),
       x = "H3K9me3 log2FoldChange", y = "H4K20me3 log2FoldChange") +
  theme_minimal(base_size = base_size) + theme(panel.grid.minor = element_blank())
save_fig(g_all, "relationship_all_regions_H3K9me3_vs_H4K20me3_by_cluster", "relationships", width = 8, height = 7)

# ================================================================
# k-means CLUSTER CHARACTERISATION (reads the cluster inputs just written above)
# ================================================================
inp        <- readRDS(file.path(rds_dir, "cluster_analysis_inputs.rds"))
sig_df     <- inp$sig_df
signal_mat <- inp$signal_mat
col_annot  <- inp$col_annot

sig_df <- sig_df[!is.na(sig_df$Cluster), ]
sig_df$Cluster <- factor(sig_df$Cluster, levels = sort(unique(sig_df$Cluster)))
clusters <- levels(sig_df$Cluster)
cat("Clusters:", paste(clusters, collapse = ", "),
    "| sizes:", paste(as.integer(table(sig_df$Cluster)), collapse = ", "), "\n")

geno_colors <- c(WT = "#1f78b4", HP1gKO = "#e31a1c")
theme_m     <- theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank())

# ================================================================
# 1. Signal profile - mean raw signal per cluster x mark x genotype.
#    Defines each cluster: which marks are present, and how WT->KO changes them.
# ================================================================
sm <- as.data.table(signal_mat, keep.rownames = "peak_id")
sm <- melt(sm, id.vars = "peak_id", variable.name = "Sample", value.name = "signal")
sm[, Sample   := as.character(Sample)]
sm[, Cluster  := sig_df$Cluster[match(peak_id, sig_df$peak_id)]]
sm[, ChIP     := as.character(col_annot$ChIP[match(Sample, col_annot$Sample)])]
sm[, Genotype := as.character(col_annot$Genotype[match(Sample, col_annot$Sample)])]

prof <- sm[!is.na(Cluster),
           .(mean_signal = mean(signal, na.rm = TRUE)),
           by = .(Cluster, ChIP, Genotype)]
prof[, Genotype := factor(Genotype, levels = c("WT", "HP1gKO"))]
prof[, ChIP := factor(ChIP, levels = c("H3K4me3", "H3K36me3", "H3K9me2", "H3K9me3", "H4K20me3"))]

# per-peak signal distribution (mean across each genotype's replicates), so the
# profile shows spread across peaks rather than a single mean bar
peak_prof <- sm[!is.na(Cluster),
                .(signal = mean(signal, na.rm = TRUE)),
                by = .(peak_id, Cluster, ChIP, Genotype)]
peak_prof[, Genotype := factor(Genotype, levels = c("WT", "HP1gKO"))]
peak_prof[, ChIP := factor(ChIP, levels = c("H3K4me3", "H3K36me3", "H3K9me2", "H3K9me3", "H4K20me3"))]
peak_prof[, Cluster := factor(Cluster)]
# winsorize per mark at 99.5% for display so extreme peaks don't dwarf the boxes
# (per-mark cap keeps free_y's per-mark scaling; only the top whisker is clipped)
peak_prof[, signal_disp := pmin(signal, quantile(signal, 0.995, na.rm = TRUE)), by = ChIP]

g_prof <- ggplot(peak_prof, aes(Cluster, signal_disp, fill = Genotype)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3, position = position_dodge(0.8)) +
  facet_wrap(~ChIP, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = geno_colors) +
  labs(title = "ChIP signal per cluster, per peak (WT vs HP1gKO)",
       subtitle = "boxes = across-peak distribution; which marks define each cluster and which are lost",
       x = "Cluster", y = "signal") +
  theme_m
save_fig(g_prof, "cluster_signal_profile", "clusters", width = 14, height = 4)

# ================================================================
# 2. Loss heatmap - log2(HP1gKO / WT) mean signal per cluster x mark.
# ================================================================
w <- dcast(prof, Cluster + ChIP ~ Genotype, value.var = "mean_signal")
eps <- 1e-3
w[, log2_KO_WT := log2((HP1gKO + eps) / (WT + eps))]

g_loss <- ggplot(w, aes(ChIP, Cluster, fill = log2_KO_WT)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", log2_KO_WT)), size = 3) +
  scale_fill_gradient2(low = "#2166ac", mid = "grey95", high = "#b2182b",
                       midpoint = 0, name = "log2 KO/WT") +
  labs(title = "Signal change per cluster (log2 HP1gKO / WT)",
       x = "Mark", y = "Cluster") +
  theme_m
save_fig(g_loss, "cluster_loss_heatmap", "clusters", width = 7, height = 4.5)

# ================================================================
# 3. ChromHMM profile - mean per-state coverage per cluster.
#    Shown two ways: absolute mean coverage, and z-scored per state (relative
#    enrichment across clusters, so genome-abundant Quies doesn't dominate).
# ================================================================
hmm_cols <- grep("^hmm_", names(sig_df), value = TRUE)
hmm_prof <- sig_df %>% group_by(Cluster) %>%
  summarise(across(all_of(hmm_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
hmm_mat <- as.matrix(hmm_prof[, hmm_cols]); rownames(hmm_mat) <- hmm_prof$Cluster
colnames(hmm_mat) <- sub("^hmm_", "", colnames(hmm_mat))

mk_hmm_long <- function(mat, value_name) {
  d <- as.data.table(mat, keep.rownames = "Cluster")
  d <- melt(d, id.vars = "Cluster", variable.name = "state", value.name = value_name)
  d$Cluster <- factor(d$Cluster, levels = clusters)
  d
}
hmm_abs <- mk_hmm_long(hmm_mat, "cov")
hmm_z   <- mk_hmm_long(scale(hmm_mat), "z")   # per-state z across clusters

g_hmm_abs <- ggplot(hmm_abs, aes(state, Cluster, fill = cov)) +
  geom_tile(colour = "white") +
  scale_fill_viridis_c(name = "mean\ncoverage") +
  labs(title = "Mean ChromHMM state coverage per cluster (absolute)",
       x = "ChromHMM state", y = "Cluster") +
  theme_m + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g_hmm_z <- ggplot(hmm_z, aes(state, Cluster, fill = z)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(low = "#2166ac", mid = "grey95", high = "#b2182b",
                       midpoint = 0, name = "z\n(per state)") +
  labs(title = "ChromHMM state enrichment per cluster (z-scored per state)",
       subtitle = "which cluster is relatively high/low in each state",
       x = "ChromHMM state", y = "Cluster") +
  theme_m + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(g_hmm_abs, "cluster_chromHMM_coverage_abs", "clusters", width = 9, height = 4)
save_fig(g_hmm_z, "cluster_chromHMM_coverage_zscore", "clusters", width = 9, height = 4)

# ================================================================
# 4. Genomic composition per cluster: (a) repeat_class %  (b) genomic_region %
#    (per-cluster gene-family fractions removed - superseded by the dedicated
#     family exon analysis now in MINUTE_5_Clustered_Gene_Families)
# ================================================================
comp_bar <- function(df, var, title) {
  d <- df %>% mutate(v = ifelse(is.na(.data[[var]]), "NA", as.character(.data[[var]]))) %>%
    count(Cluster, v) %>% group_by(Cluster) %>% mutate(frac = n / sum(n)) %>% ungroup()
  ggplot(d, aes(Cluster, frac, fill = v)) +
    geom_col() +
    labs(title = title, x = "Cluster", y = "fraction of peaks", fill = var) +
    theme_m
}
save_fig(comp_bar(sig_df, "repeat_class", "Repeat-class composition per cluster"),
         "cluster_repeat_composition", "clusters", width = 7, height = 4.5)
save_fig(comp_bar(sig_df, "genomic_region", "Genomic-region composition per cluster"),
         "cluster_region_composition", "clusters", width = 8, height = 4.5)

# ================================================================
# 5. Cluster summary table
# ================================================================
summ <- sig_df %>% group_by(Cluster) %>%
  summarise(
    n_peaks       = n(),
    median_kb     = round(median(peak_size_kb, na.rm = TRUE), 1),
    median_l2fc   = round(median(log2FoldChange, na.rm = TRUE), 2),
    top_state     = names(sort(table(chromHMM_state), decreasing = TRUE))[1],
    top_repeat    = names(sort(table(repeat_class),  decreasing = TRUE))[1],
    pct_pericentro = round(100 * mean(start < 3e6), 1),
    pct_Pcdh      = round(100 * mean(grepl("^Pcdh", ifelse(is.na(SYMBOL), "", SYMBOL))), 1),
    pct_Zfp       = round(100 * mean(grepl("^Zfp",  ifelse(is.na(SYMBOL), "", SYMBOL))), 1),
    pct_Vmn       = round(100 * mean(grepl("^Vmn",  ifelse(is.na(SYMBOL), "", SYMBOL))), 1),
    .groups = "drop")
fwrite(summ, file.path(tables_dir, "cluster_summary.tsv"), sep = "\t")
cat("\n=== Cluster summary ===\n"); print(as.data.frame(summ))
message("Saved: cluster_summary.tsv")

# ================================================================
# 6. Per-cluster ChromHMM enrichment TEST
#    Formalises the part-3 coverage heatmaps: for each cluster x state, compare
#    the state's coverage in that cluster vs all OTHER clustered peaks (Wilcoxon
#    + BH FDR, log2 mean-coverage ratio). Tells which states each cluster is
#    distinctively enriched/depleted for, relative to the other clusters.
# ================================================================
clus_enrich <- do.call(rbind, lapply(clusters, function(cl) {
  inc <- sig_df$Cluster == cl
  # per-region (unweighted) Wilcoxon: this cluster vs all other clustered peaks
  pr <- do.call(rbind, lapply(hmm_cols, function(col) {
    a <- sig_df[[col]][inc]; b <- sig_df[[col]][!inc]
    data.frame(Cluster = cl, state = sub("^hmm_", "", col), n = sum(inc),
               mean_cov_in = mean(a, na.rm = TRUE), mean_cov_out = mean(b, na.rm = TRUE),
               log2_ratio = log2((mean(a, na.rm = TRUE) + 1e-4) / (mean(b, na.rm = TRUE) + 1e-4)),
               p_value = tryCatch(suppressWarnings(wilcox.test(a, b)$p.value),
                                  error = function(e) NA_real_),
               stringsAsFactors = FALSE)
  }))
  # size-weighted (by domain length) with permutation test
  sw <- chromHMM_size_weighted(sig_df[, hmm_cols], sig_df$peak_size_kb, inc)
  merge(pr, sw, by = "state")
}))
clus_enrich <- clus_enrich %>% group_by(Cluster) %>%
  mutate(p_adj_BH   = p.adjust(p_value, method = "BH"),
         perm_p_adj = p.adjust(perm_p,  method = "BH")) %>% ungroup()
fwrite(clus_enrich, file.path(tables_dir, "cluster_chromHMM_enrichment.tsv"), sep = "\t")

# Two dot plots: per-region (Wilcoxon) and size-weighted (permutation)
clus_dotplot <- function(lr, fdr, ttl, sub) {
  clus_enrich %>%
    filter(is.finite(.data[[lr]]), !is.na(.data[[fdr]])) %>%
    mutate(Cluster = factor(paste("Cluster", Cluster), levels = paste("Cluster", clusters)),
           neglog10FDR = -log10(pmax(.data[[fdr]], .Machine$double.xmin)),
           state = factor(state, levels = sort(unique(state)))) %>%
    ggplot(aes(Cluster, state, size = neglog10FDR, colour = .data[[lr]])) +
    geom_point() +
    scale_colour_gradient2(low = "steelblue3", mid = "grey90", high = "firebrick2",
                           midpoint = 0, name = expression(log[2]~"ratio")) +
    scale_size_continuous(name = expression(-log[10]("FDR"))) +
    labs(title = ttl, subtitle = sub, x = NULL, y = "ChromHMM state") +
    theme_m + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}
save_fig(clus_dotplot("log2_ratio", "p_adj_BH",
                    "Per-cluster ChromHMM enrichment (cluster vs other clusters) - PER-REGION",
                    "colour = log2 coverage ratio; each region equal; size = FDR"),
         "cluster_chromHMM_enrichment", "clusters", width = 8, height = 6)
save_fig(clus_dotplot("w_log2_ratio", "perm_p_adj",
                    "Per-cluster ChromHMM enrichment (cluster vs other clusters) - SIZE-WEIGHTED",
                    "colour = log2 territory ratio (weighted by domain length); size = perm FDR"),
         "cluster_chromHMM_enrichment_sizeweighted", "clusters", width = 8, height = 6)


# ================================================================
# Structural + chromatin-state enrichment of the significant set
#   TAD-boundary (hypergeometric) + ChromHMM-state (coverage) enrichment.
#   Repeat-class enrichment lives in MINUTE_4_Repeats. hyper_row() from config.R.
# ================================================================
sig_all <- do.call(rbind, lapply(names(annotated_results), function(mark) {
  df <- annotated_results[[mark]]; df$ChIP <- mark
  df$peak_id <- paste0(df$chr, ":", df$start, "-", df$end); df
}))
sig_all <- sig_all[is_significant(sig_all), ]

tad_enrich <- do.call(rbind, lapply(names(annotated_results), function(mark) {
  all_df <- annotated_results[[mark]]; all_df$ChIP <- mark
  all_df$peak_id <- paste0(all_df$chr, ":", all_df$start, "-", all_df$end)
  flag_TAD <- all_df$overlaps_with_Tad_boundary %in% TRUE
  hyper_row(mark, "TAD_boundary", all_df, sig_all[sig_all$ChIP == mark, ], flag_TAD)
}))
tad_enrich$p_adj_BH <- p.adjust(tad_enrich$p_value, "BH")
fwrite(tad_enrich, file.path(tables_dir, "enrichment_TAD.tsv"), sep = "\t")

g_tad <- tad_enrich %>%
  filter(is.finite(odds_ratio), !is.na(p_adj_BH)) %>%
  mutate(neglog10FDR = -log10(pmax(p_adj_BH, .Machine$double.xmin)),
         ChIP = factor(ChIP, levels = mark_levels)) %>%
  ggplot(aes(ChIP, odds_ratio, size = neglog10FDR, colour = odds_ratio > 1)) +
  geom_hline(yintercept = 1, linetype = "dashed") + geom_point() +
  scale_colour_manual(values = c(`TRUE` = "firebrick2", `FALSE` = "steelblue3"), guide = "none") +
  scale_y_log10() +
  labs(title = "TAD-boundary enrichment of significant peaks", x = NULL,
       y = "odds ratio (log10)", size = expression(-log[10]("FDR"))) + theme_m
save_fig(g_tad, "enrichment_TAD_dotplot", "enrichment", width = 6, height = 4)
message("Saved: enrichment_TAD.tsv + enrichment/enrichment_TAD_dotplot.{png,pdf}")

# ================================================================
# ChromHMM chromatin-state enrichment (coverage-based)
# ----------------------------------------------------------------
# ChromHMM is annotated per peak in MINUTE_1 as per-state COVERAGE FRACTIONS
# (hmm_<state> columns) rather than a single size-biased label, because the
# H4K20me3/H3K9me3 domains span many states. For each mark x state we compare
# the mean coverage fraction in the significant set vs the background
# (non-significant peaks): log2 ratio of means = effect size, Wilcoxon = test.
# ================================================================

hmm_cols   <- grep("^hmm_", names(annotated_results[[1]]), value = TRUE)
hmm_states <- sub("^hmm_", "", hmm_cols)

chromhmm_enrich <- do.call(rbind, lapply(names(annotated_results), function(mark) {
  d   <- annotated_results[[mark]]
  d$ChIP <- mark
  sig <- is_significant(d)
  # per-region (unweighted) Wilcoxon
  pr <- do.call(rbind, lapply(hmm_states, function(st) {
    col <- paste0("hmm_", st)
    a <- d[[col]][sig]; b <- d[[col]][!sig]           # background = non-significant
    data.frame(ChIP = mark, state = st, n_sig = sum(sig),
               mean_cov_sig = mean(a, na.rm = TRUE), mean_cov_bg = mean(b, na.rm = TRUE),
               log2_ratio = log2((mean(a, na.rm = TRUE) + 1e-4) / (mean(b, na.rm = TRUE) + 1e-4)),
               p_value = tryCatch(suppressWarnings(wilcox.test(a, b)$p.value),
                                  error = function(e) NA_real_),
               stringsAsFactors = FALSE)
  }))
  # size-weighted (by domain length) with permutation test
  sw <- chromHMM_size_weighted(d[, hmm_cols], d$end - d$start, sig)
  merge(pr, sw, by = "state")
}))

chromhmm_enrich <- chromhmm_enrich %>%
  dplyr::group_by(ChIP) %>%
  dplyr::mutate(p_adj_BH   = p.adjust(p_value, method = "BH"),
                perm_p_adj = p.adjust(perm_p,  method = "BH")) %>%
  dplyr::ungroup()

fwrite(chromhmm_enrich, file.path(tables_dir, "enrichment_chromHMM.tsv"), sep = "\t")

# Two dot plots: PER-REGION (Wilcoxon) and SIZE-WEIGHTED (permutation) - report both
hmm_dotplot <- function(df, lr, fdr, ttl, sub) {
  df %>% dplyr::filter(is.finite(.data[[lr]]), !is.na(.data[[fdr]])) %>%
    dplyr::mutate(neglog10FDR = -log10(pmax(.data[[fdr]], .Machine$double.xmin)),
                  state = factor(state, levels = sort(unique(state)))) %>%
    ggplot(aes(x = ChIP, y = state, size = neglog10FDR, colour = .data[[lr]])) +
    geom_point() +
    scale_colour_gradient2(low = "steelblue3", mid = "grey90", high = "firebrick2",
                           midpoint = 0, name = expression(log[2]~"ratio")) +
    scale_size_continuous(name = expression(-log[10]("FDR"))) +
    labs(title = ttl, subtitle = sub, x = NULL, y = "ChromHMM state") +
    theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank())
}
save_fig(hmm_dotplot(chromhmm_enrich, "log2_ratio", "p_adj_BH",
                     "ChromHMM enrichment of significant peaks - PER-REGION",
                     "each region equal; colour = log2 mean-coverage ratio (sig vs bg); size = FDR"),
         "enrichment_chromHMM_dotplot", "enrichment", width = 8, height = 7)
save_fig(hmm_dotplot(chromhmm_enrich, "w_log2_ratio", "perm_p_adj",
                     "ChromHMM enrichment of significant peaks - SIZE-WEIGHTED",
                     "weighted by domain length; colour = log2 territory ratio (sig vs bg); size = perm FDR"),
         "enrichment_chromHMM_sizeweighted_dotplot", "enrichment", width = 8, height = 7)

# NOTE: the H3K9me3-loss-vs-H3K9me3-retained comparison within H4K20me3-lost
# regions lives in MINUTE_3_Differential_loss (co_loss vs H4K20me3_only), which
# does it for ChromHMM + repeats + regions + genes, per-region and size-weighted.
