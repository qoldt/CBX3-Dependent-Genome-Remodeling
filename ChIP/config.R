# ================================================================
# MINUTE pipeline - shared configuration & parameters
# ----------------------------------------------------------------
# Sourced by MINUTE_1, MINUTE_2, MINUTE_3 and by run_MINUTE.R.
# SINGLE SOURCE OF TRUTH for parameters - edit values here only;
# the stage scripts must not redefine them.
#
# Run the pipeline FROM the ChIP/ directory (the folder holding
# this file); all repo-relative paths below assume that:
#     cd ChIP && Rscript run_MINUTE.R
#
# I/O layout
#   Inputs (read-only):
#     samples.tsv                sample sheet - SINGLE SOURCE OF TRUTH for samples
#     data/peaks/                consensus / master peak BEDs      (committed, small)
#     $MINUTE_DATA/bigwig/       scaled mm39 bigWigs               (external, large)
#     $MINUTE_DATA/annotation/   LINE / SINE / LTR / TAD BEDs      (external)
#   Outputs (generated, git-ignored except results/rds/):
#     results/counts/  results/rds/  results/tables/  results/figures/  results/bed/
# ================================================================

# --- Roots -------------------------------------------------------
# Large raw inputs (bigWigs + repeat/TAD annotation) are NOT committed. Download
# them from the Google Drive links in the README into the repo-relative folders
# ChIP/data/bigwig and ChIP/data/annotation (the defaults below; scripts run from
# ChIP/). To use a copy staged elsewhere (e.g. a Synology/external mount), set the
# MINUTE_DATA env var:  MINUTE_DATA=/path/to/store Rscript run_MINUTE.R
data_ext   <- path.expand(Sys.getenv("MINUTE_DATA", unset = "data"))
bigwig_dir <- file.path(data_ext, "bigwig")
annot_dir  <- file.path(data_ext, "annotation")

# Repo-relative input / output dirs (scripts are run from ChIP/)
# results_dir is overridable so a sensitivity analysis can be written to a
# SEPARATE tree without disturbing the committed run:
#   MINUTE_RESULTS=results_norep2 MINUTE_EXCLUDE=... Rscript run_MINUTE.R
peaks_dir   <- "data/peaks"
results_dir <- path.expand(Sys.getenv("MINUTE_RESULTS", unset = "results"))
counts_dir  <- file.path(results_dir, "counts")
rds_dir     <- file.path(results_dir, "rds")
tables_dir  <- file.path(results_dir, "tables")
fig_dir     <- file.path(results_dir, "figures")
bed_dir     <- file.path(results_dir, "bed")
for (d in c(counts_dir, rds_dir, tables_dir, fig_dir, bed_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# --- Figure savers: write BOTH png (quick view) and pdf (vector for papers), ---
# --- organised into per-analysis subfolders under results/figures/<subdir>.  ---
# save_fig(): for ggplot objects.  save_base_fig(): for base-graphics / draw()
# renderers (e.g. ComplexHeatmap) - pass a zero-arg function that draws.
save_fig <- function(plot, name, subdir = "", width = 8, height = 6, dpi = 300) {
  d <- if (nzchar(subdir)) file.path(fig_dir, subdir) else fig_dir
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(d, paste0(name, ".png")), plot, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(file.path(d, paste0(name, ".pdf")), plot, width = width, height = height)
  message("Saved: ", file.path(subdir, paste0(name, ".{png,pdf}")))
  invisible(file.path(d, name))
}
save_base_fig <- function(draw_fn, name, subdir = "", width = 10, height = 12, dpi = 300) {
  d <- if (nzchar(subdir)) file.path(fig_dir, subdir) else fig_dir
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  png(file.path(d, paste0(name, ".png")), width = width, height = height, units = "in", res = dpi)
  draw_fn(); dev.off()
  pdf(file.path(d, paste0(name, ".pdf")), width = width, height = height)
  draw_fn(); dev.off()
  message("Saved: ", file.path(subdir, paste0(name, ".{png,pdf}")))
  invisible(file.path(d, name))
}

# File written by MINUTE_1, read by MINUTE_2 and MINUTE_3
annotated_rds <- file.path(rds_dir, "annotated_results_H3K9me3_H4K20me3_2000bp_merged.rds")

# --- Libraries (union used across the pipeline) ---
# wigglescout is sensitive to its future/furrr stack. This is the ONLY combination
# confirmed working here; restart R after installing so newer loaded versions unload:
#   remotes::install_version("furrr",   version = "0.2.3")
#   remotes::install_version("future",  version = "1.23.0")
#   remotes::install_version("globals", version = "0.14.0")
#   remotes::install_github("cnluzon/wigglescout")
# If bw_loci() fails with "values must be length 1, but FUN(X[[i]]) result is
# length N", the pins have drifted (or a stray future::plan() is set) - reinstall
# the versions above / run future::plan("sequential"); do NOT edit the R scripts.
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(AnnotationDbi)
  library(TxDb.Mmusculus.UCSC.mm39.knownGene)
  library(org.Mm.eg.db)
  library(GenomicRanges)
  library(rtracklayer)
  library(wigglescout)
  library(parallel)
  library(data.table)
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(ltc)
})

# --- Colour palettes: SINGLE SOURCE OF TRUTH for figure colours ---------------
# Everything comes from the `ltc` package (install.packages("ltc")):
#   discrete / categorical -> "gaby"      (5 colours)
#   continuous / heatmaps  -> "heatmap0"  (9 colours, dark blue -> cream -> dark red)
# Stage scripts must use the helpers below, not raw hex codes. Greys are kept for
# "neutral / background / not-a-category" levels on purpose - they are meant to
# recede, so they are deliberately NOT drawn from the palette.

# gaby, in palette order: pale yellow, salmon, pale blue, blue, plum
gaby_cols <- unname(ltc("gaby"))

# gaby has only 5 colours. INTERPOLATING past that is not viable: gaby_pal(7)
# put 5' UTR (#F4BE99) next to Distal Intergenic (#D8B1A3) on the genomic_region
# change plots, which are indistinguishable in print - and Distal Intergenic is
# the bulk of the points. So beyond 5 levels we switch palette entirely rather
# than blend, using ltc's "hat" (10 colours), reordered so the first n are as
# far apart in hue as possible: blue, red, green, purple, yellow, teal,
# crimson, orange, green2, black.
hat_cols <- unname(ltc("hat"))[c(6, 3, 8, 5, 1, 7, 4, 2, 9, 10)]

# The project's categorical palette: gaby <=5 levels, hat 6-10, interpolated
# hat beyond that (>10 categories is a labelling problem, not a colour one).
disc_pal <- function(n) {
  if (n <= length(gaby_cols))     gaby_cols[seq_len(n)]
  else if (n <= length(hat_cols)) hat_cols[seq_len(n)]
  else                            unname(ltc("hat", n, "continuous"))
}
scale_fill_disc   <- function(...) discrete_scale(aesthetics = "fill",   palette = disc_pal, ...)
scale_colour_disc <- function(...) discrete_scale(aesthetics = "colour", palette = disc_pal, ...)

# heatmap0, in palette order (index 5 = "#E9D8A6" is the pale midpoint)
heat_cols <- unname(ltc("heatmap0"))
heat_low  <- heat_cols[2]   # teal-blue
heat_mid  <- heat_cols[5]   # cream
heat_high <- heat_cols[8]   # dark red
# Sequential (0 -> max): the full ramp.
scale_fill_heat0   <- function(...) scale_fill_gradientn(colours = heat_cols, ...)
scale_colour_heat0 <- function(...) scale_colour_gradientn(colours = heat_cols, ...)
# Diverging (log2FC, z-score): the ramp's own ends around its cream midpoint.
scale_fill_heat0_div   <- function(midpoint = 0, ...) scale_fill_gradient2(
  low = heat_low, mid = heat_mid, high = heat_high, midpoint = midpoint, ...)
scale_colour_heat0_div <- function(midpoint = 0, ...) scale_colour_gradient2(
  low = heat_low, mid = heat_mid, high = heat_high, midpoint = midpoint, ...)
# ComplexHeatmap / circlize: pass the break points, get a matching colorRamp2.
#   symmetric z-score : heat_col_fun(seq(-2, 2, length.out = length(heat_cols)))
#   sequential        : heat_col_fun(c(0.3, 1, 6))
heat_col_fun <- function(breaks) circlize::colorRamp2(
  breaks, grDevices::colorRampPalette(heat_cols)(length(breaks)))

# --- Shared plotting constants (used across stages) ---
geno_colors <- c(WT = gaby_cols[4], HP1gKO = gaby_cols[5])
theme_m     <- theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank())
mark_levels <- c("H3K4me3", "H3K36me3", "H3K9me2", "H3K9me3", "H4K20me3")

# Repeat classes (shared by MINUTE_2's heatmap annotation + scatter plots and
# MINUTE_3's composition bars).
# USES "hat", NOT gaby, deliberately. All four repeat classes are abundant on
# the change plots and need strong mutual contrast, and gaby cannot supply four
# such colours: two of its five entries are pale (#fceaab, #a8c4cc) and vanish
# against white, and its two blues are adjacent in hue. Earlier attempts here
# put LINE/LTR on the two adjacent blues (the two classes carrying the result -
# young L1 and IAP - became the hardest pair to distinguish), then put the
# abundant "complex" class on pale yellow, where it was invisible.
# hat's first four are blue / red / green / purple - all strong, all distinct.
# LINE and LTR take blue and red, the most separated pair.
repeat_palette <- c(LINE = hat_cols[1],    # blue
                    LTR  = hat_cols[2],    # red
                    SINE = hat_cols[3],    # green
                    complex = hat_cols[4], # purple
                    none = "grey80")

# --- Sample sheet: SINGLE SOURCE OF TRUTH for samples ------------
# Columns: sample_id, mark, genotype (WT/HP1gKO), replicate, bigwig (filename).
# Genotype and replicate are read from here - NOT hard-coded in the stages.
samples <- as.data.frame(data.table::fread("samples.tsv"))

# --- Sample exclusion --------------------------------------------------------
# DEFAULT: HP1gKO rep2 is EXCLUDED from every analysis.
#
# WHY. Replicates are NESTED, not independent: WT rep1-3 and rep4-6 are technical
# replicates of biological replicates 1 and 2; likewise HP1gKO rep2-3 (rep1's
# library failed) and rep4-6. HP1gKO rep2 is the highest-signal library in ALL
# FIVE marks, and it is badly discordant with rep3 - its own technical twin
# (H4K20me3 median coverage 1.92 vs 1.14, while the bio-2 trio sits at
# 0.99-1.32). Technical replicates should agree; this one does not, so it is
# treated as a technical failure rather than biological variation.
#
# EFFECT OF EXCLUDING IT: the global loss deepens (H4K20me3 median log2FC
# -0.474 -> -0.678) while the control marks stay near zero, and the apparent
# "gain" set roughly halves. Note ~0.04-0.07 of the deepening is a general
# offset from dropping the highest KO library, so it is not all correction.
#
# Override with MINUTE_EXCLUDE: a comma-separated list of sample_id values or of
# "<genotype>:<replicate>" pairs that expand across every mark. Use
# MINUTE_EXCLUDE=none to keep ALL samples (reproduces the pre-exclusion run):
#   MINUTE_EXCLUDE=none MINUTE_RESULTS=results_allreps Rscript run_MINUTE.R
#
# CAVEAT that exclusion does NOT fix: the DESeq2 design ~condition still treats
# technical replicates as independent, so the effective n is 2 biological
# replicates vs 2, not 6 vs 4. Per-peak p-values are anticonservative on that
# count - though the p-value distributions (MINUTE_2 section 7) show they are
# also conservatively miscalibrated by the global shift, which dominates.
exclude_samples <- trimws(strsplit(Sys.getenv("MINUTE_EXCLUDE", unset = "HP1gKO:2"), ",")[[1]])
if (identical(tolower(exclude_samples), "none")) exclude_samples <- character(0)
exclude_samples <- exclude_samples[nzchar(exclude_samples)]
if (length(exclude_samples)) {
  pairs <- grepl(":", exclude_samples, fixed = TRUE)
  drop  <- samples$sample_id %in% exclude_samples[!pairs]
  for (p in exclude_samples[pairs]) {
    gr <- strsplit(p, ":", fixed = TRUE)[[1]]
    drop <- drop | (as.character(samples$genotype) == gr[1] &
                    as.character(samples$replicate) == gr[2])
  }
  message("Excluding ", sum(drop), " sample(s): ",
          paste(samples$sample_id[drop], collapse = ", "))
  samples <- samples[!drop, , drop = FALSE]
}

samples$genotype    <- factor(samples$genotype, levels = c("WT", "HP1gKO"))
samples$bigwig_path <- file.path(bigwig_dir, samples$bigwig)
marks <- unique(as.character(samples$mark))   # preserves sample-sheet order

# --- Associated peak files per mark (under data/peaks) ---
regions <- list(
  "H3K4me3"  = file.path(peaks_dir, "master_peaks_consensus/H3K4me3_consensus_masterPeak.bed"),
  "H3K9me2"  = file.path(peaks_dir, "merged_H3K9me2_masterPeak.bed"),
  "H3K9me3"  = file.path(peaks_dir, "2000bp_merge_H3K9me3_H4K20me3.bed"),
  "H4K20me3" = file.path(peaks_dir, "2000bp_merge_H3K9me3_H4K20me3.bed"),
  "H3K36me3" = file.path(peaks_dir, "merged_H3K36me3_masterPeak.bed")
)

# --- Repeat / TAD annotation files (external, under $MINUTE_DATA/annotation) ---
repeat_bed_files <- list(
  LINE = file.path(annot_dir, "LINE.mm39.bed"),
  SINE = file.path(annot_dir, "SINE.mm39.bed"),
  LTR  = file.path(annot_dir, "LTR.mm39.bed")
)
tad_bed_file <- file.path(annot_dir, "TAD_boundaries_mm39.bed")

# ChromHMM 18-state cortex (P0) segmentation, mm39. Annotated per-peak in
# MINUTE_1 (column chromHMM_state), enriched in MINUTE_2, used in MINUTE_4.
chromhmm_bed <- file.path(annot_dir, "ChromHMM_18state_CortexP0_mm39.bed")

# KAP1/TRIM28 ChIP-seq peaks (mm39, Ensembl seqnames), prepared by liftOver from
# ChIP-Atlas mm10 (see README methods). Used by MINUTE_6; skipped if absent.
kap1_neural_bed  <- file.path(annot_dir, "KAP1_Neural_mm39.bed")   # Neuro-2a (GSE110032)
kap1_allcell_bed <- file.path(annot_dir, "KAP1_allCell_mm39.bed")  # ChIP-Atlas all-cell aggregate

# --- Helper: load the ChromHMM segmentation as a GRanges (state in $state) ---
# The file has a few overlapping same-state segments, which would double-count
# coverage; merge overlapping segments within each state so per-state coverage
# is well defined (fraction of a peak in any one state cannot exceed 1).
load_chromHMM <- function() {
  gr <- rtracklayer::import(chromhmm_bed, format = "bed")   # 4th BED col -> $name
  mcols(gr)$state <- as.character(mcols(gr)$name)
  red <- unlist(reduce(split(gr, mcols(gr)$state)))
  mcols(red)$state <- names(red)
  names(red) <- NULL
  red
}

# --- Helper: per-peak coverage FRACTION of each ChromHMM state ---
# Domains span many states, so we keep the full breakdown rather than a single
# size-biased label. gr: peaks GRanges (any seqnames style). Returns a matrix
# [peaks x states]; each value = fraction of that peak's width covered by the
# state (rows sum to <= 1). Downstream: dominant label + coverage-based enrichment.
chromHMM_coverage <- function(gr, chromhmm_gr) {
  suppressWarnings(seqlevelsStyle(chromhmm_gr) <- seqlevelsStyle(gr)[1])
  states <- sort(unique(as.character(mcols(chromhmm_gr)$state)))
  m <- matrix(0, nrow = length(gr), ncol = length(states),
              dimnames = list(NULL, states))
  hits <- findOverlaps(gr, chromhmm_gr, ignore.strand = TRUE)
  if (length(hits) > 0) {
    ov  <- pintersect(gr[queryHits(hits)], chromhmm_gr[subjectHits(hits)],
                      ignore.strand = TRUE)
    agg <- data.table(q = queryHits(hits),
                      s = as.character(mcols(chromhmm_gr)$state[subjectHits(hits)]),
                      w = width(ov))[, .(w = sum(w)), by = .(q, s)]
    m[cbind(agg$q, match(agg$s, states))] <- agg$w
    m <- pmin(m / width(gr), 1)   # fraction; cap at 1 (rare cross-state overlaps)
  }
  m
}

# --- Helper: SIZE-WEIGHTED per-state coverage enrichment + permutation test ---
# Complements the per-region (unweighted) Wilcoxon tests by weighting each region
# by its length, so the result reflects the changed genomic TERRITORY rather than
# treating a 2 Mb domain and a 2 kb peak equally.
#   cov_mat : regions x states coverage fractions (the hmm_<state> columns)
#   size    : region length (bp or kb - only relative weights matter)
#   group   : logical, TRUE = foreground (e.g. significant / cluster / H3K9me3-loss)
# Returns per state: weighted mean coverage in fg vs bg (= fraction of that
# group's territory in the state), log2 ratio, and a permutation p-value that
# shuffles region LABELS (domain-level; a bp-level test would pseudo-replicate).
chromHMM_size_weighted <- function(cov_mat, size, group, nperm = 1000, seed = 42) {
  cov_mat <- as.matrix(cov_mat)
  keep    <- is.finite(size) & !is.na(group)
  cov_mat <- cov_mat[keep, , drop = FALSE]; size <- size[keep]; group <- as.logical(group)[keep]
  covW <- cov_mat * size                       # region x state, bp-in-state
  totW <- sum(size)
  colW <- colSums(covW, na.rm = TRUE)          # total state-bp over all regions
  gv   <- as.numeric(group)
  fgbp <- as.numeric(crossprod(gv, covW))      # state-bp in foreground
  sf   <- sum(gv * size); sb <- totW - sf
  fg   <- fgbp / sf; bg <- (colW - fgbp) / sb
  obsd <- fg - bg
  set.seed(seed)
  n <- length(group); k <- sum(group); cnt <- numeric(ncol(cov_mat))
  for (i in seq_len(nperm)) {
    g <- numeric(n); g[sample.int(n, k)] <- 1
    fb <- as.numeric(crossprod(g, covW)); sff <- sum(g * size)
    cnt <- cnt + (abs(fb / sff - (colW - fb) / (totW - sff)) >= abs(obsd))
  }
  data.frame(state = sub("^hmm_", "", colnames(cov_mat)),
             w_cov_fg = fg, w_cov_bg = bg,
             w_log2_ratio = log2((fg + 1e-4) / (bg + 1e-4)),
             perm_p = (cnt + 1) / (nperm + 1),
             stringsAsFactors = FALSE)
}

# --- Helper: dominant state + its purity from a coverage matrix ---
# Convenience single label = the state covering the largest fraction of the peak,
# plus that fraction (purity). NA where a peak overlaps no ChromHMM segment.
chromHMM_dominant <- function(cov) {
  rs <- rowSums(cov)
  mc <- max.col(cov, ties.method = "first")
  data.frame(
    chromHMM_state  = ifelse(rs > 0, colnames(cov)[mc], NA_character_),
    chromHMM_purity = ifelse(rs > 0, cov[cbind(seq_len(nrow(cov)), mc)], NA_real_),
    stringsAsFactors = FALSE
  )
}

# --- Helper: UN-ROUNDED mean coverage per peak (a PLOTTING covariate) --------
# MINUTE_1 rounds the per-peak mean coverage to integers to satisfy DESeq2, which
# collapses `baseMean` onto a coarse k/n_samples lattice - for H4K20me3 that is
# only ~214 distinct values across 116k peaks, with one value shared by >9,000
# peaks. That is far too quantised to encode as a point size.
#
# The counts TSV is written BEFORE that rounding, so it keeps full precision.
# Because sizeFactors = 1, the row mean of that file is EXACTLY what `baseMean`
# would have been had we never rounded - the same quantity, un-degraded.
#
# Use this for AESTHETICS ONLY (point size, colour). All statistics must keep
# using DESeq2's baseMean/log2FoldChange/pvalue, which are computed from the
# rounded integer counts.
# Returns a named numeric vector (peak_id -> mean scaled coverage), or NULL if
# the counts file is missing.
mean_coverage_by_peak <- function(mark) {
  f <- file.path(counts_dir, paste0(mark, "_bigwig_counts.tsv"))
  if (!file.exists(f)) {
    warning("mean_coverage_by_peak(): no counts file for ", mark, " - ", f)
    return(NULL)
  }
  dt  <- data.table::fread(f)
  ids <- as.character(dt[[1]])                       # peak_id column
  m   <- as.matrix(dt[, -1, with = FALSE])
  stats::setNames(rowMeans(m, na.rm = TRUE), ids)
}

# --- Helper: bigWig paths for one mark, named by sample_id ---
# Order follows the sample sheet (WT rep1-6, then HP1gKO rep2-6).
get_bigwig_files <- function(mark) {
  s <- samples[as.character(samples$mark) == mark, , drop = FALSE]
  setNames(s$bigwig_path, s$sample_id)
}

# --- Helper: load a UCSC-style repeat BED into a GRanges ---
loadRepeatBED <- function(filepath) {
  expected_colnames <- c("bin", "swScore", "milliDiv", "milliDel", "milliIns",
                         "genoName", "genoStart", "genoEnd", "genoLeft", "strand",
                         "repName", "repClass", "repFamily", "repStart", "repEnd",
                         "repLeft", "id")
  df <- read.table(filepath, sep = "\t", header = FALSE, comment.char = "#", quote = "",
                   fill = TRUE, stringsAsFactors = FALSE,
                   col.names = expected_colnames, colClasses = rep("character", 17))
  df$genoStart <- suppressWarnings(as.numeric(df$genoStart))
  df$genoEnd   <- suppressWarnings(as.numeric(df$genoEnd))
  df <- df[!is.na(df$genoStart) & !is.na(df$genoEnd), ]
  GRanges(
    seqnames = df$genoName,
    ranges   = IRanges(start = df$genoStart + 1, end = df$genoEnd),
    strand   = df$strand,
    repName  = df$repName,
    repClass = df$repClass,
    repFamily = df$repFamily
  )
}

# --- Helper: load repeat + TAD annotation used by MINUTE_1 and MINUTE_2 ---
# Returns a list; the stage scripts assign to line_gr/sine_gr/ltr_gr/tad_gr.
load_annotation <- function() {
  ann <- list(
    line = loadRepeatBED(repeat_bed_files$LINE),
    sine = loadRepeatBED(repeat_bed_files$SINE),
    ltr  = loadRepeatBED(repeat_bed_files$LTR),
    tad  = import(tad_bed_file, format = "bed")
  )
  seqlevelsStyle(ann$tad) <- "UCSC"
  ann
}

# --- Helper: one hypergeometric enrichment row (OR + direction + p) for a
#     boolean annotation flag over ALL measured peaks vs the significant subset.
#     Shared by MINUTE_2 (TAD/ChromHMM) and MINUTE_4 (repeats). ---
hyper_row <- function(chip, annotation, all_df, sig_df, flag_all) {
  N <- nrow(all_df); K <- sum(flag_all, na.rm = TRUE); n <- nrow(sig_df)
  flag_sig <- flag_all[match(sig_df$peak_id, all_df$peak_id)]
  k <- sum(flag_sig, na.rm = TRUE)
  if (is.na(N) || N == 0 || is.na(n) || is.na(K) || K < 0) {
    pval <- expected <- enrich <- NA_real_
  } else {
    expected <- if (N > 0) n * (K / N) else NA_real_
    pval <- if (n > 0 && K > 0 && N > K) phyper(k - 1, K, N - K, n, lower.tail = FALSE) else NA_real_
    enrich <- if (!is.na(expected) && expected > 0) k / expected else NA_real_
  }
  a <- k; b <- n - k; c <- K - k; d <- N - K - b
  if (any(c(a, b, c, d) <= 0)) { a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d <- d + 0.5 }
  or_val <- (a * d) / (b * c); if (!is.finite(or_val)) or_val <- NA_real_
  if (!is.null(sig_df$log2FoldChange) && k > 0) {
    l2fc_med <- stats::median(sig_df$log2FoldChange[which(flag_sig)], na.rm = TRUE)
    direction <- if (is.finite(l2fc_med) && l2fc_med > 0) "enriched" else "depleted"
  } else { l2fc_med <- NA_real_; direction <- NA_character_ }
  data.frame(ChIP = chip, annotation = annotation, N_total = N, K_in_annotation = K,
             n_signif = n, k_signif_in_annotation = k, expected_overlap = expected,
             enrichment = enrich, odds_ratio = or_val,
             l2fc_median_sig_in_annotation = l2fc_med, direction = direction,
             p_value = pval, stringsAsFactors = FALSE)
}

# --- Gene families silenced by H3K9me3 / HUSH (CBX3-dependent) ---
# Quantified directly over their exons in MINUTE_4 (bw_loci), rather than relying
# on individual genes passing per-peak significance. Symbol-prefix definitions:
gene_families <- list(
  Protocadherin = "^Pcdh(a|b|g)",     # clustered protocadherins (chr18 a/b/g clusters)
  `KRAB-ZFP`    = "^Zfp",             # zinc-finger proxy for KRAB-ZFPs (large tandem clusters)
  Vomeronasal   = "^Vmn",            # Vmn1r / Vmn2r vomeronasal receptors
  Olfactory     = "^(Olfr|Or[0-9])"  # olfactory receptors: current (Or*) + legacy (Olfr*) names
)

# Which exon carries the H3K9me3/HUSH silencing, PER FAMILY — so we quantify the
# biologically relevant exon per gene rather than diluting it across all exons:
#   KRAB-ZFP      -> "last"    (3'-most exon: the ZNF array)
#   Protocadherin -> "first"   (5'-most exon: the large variable exon)
#   Vomeronasal / Olfactory -> "largest" (the single coding exon)
family_exon_rule <- list(
  Protocadherin = "first",
  `KRAB-ZFP`    = "last",
  Vomeronasal   = "largest",
  Olfactory     = "largest"
)

# --- Helper: one silencing-relevant exon PER GENE for each family (mm39 TxDb) ---
# Returns a GRanges with one interval per family member gene (the exon selected by
# family_exon_rule, strand-aware), with `family` and `gene` metadata. Seqnames are
# UCSC; set the style at query time.
family_exons <- function() {
  txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
  exg  <- exonsBy(txdb, by = "gene")                       # by Entrez ID
  sym  <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, keys = names(exg),
            column = "SYMBOL", keytype = "ENTREZID"))
  sym  <- ifelse(is.na(sym), "", sym)
  pick <- function(ex, rule) {                             # ex = one gene's exons
    if (rule == "largest") return(ex[which.max(width(ex))])
    ex <- ex[if (as.character(strand(ex))[1] == "-") order(-start(ex)) else order(start(ex))]
    if (rule == "first") ex[1] else ex[length(ex)]         # "first" = 5', else "last" = 3'
  }
  parts <- lapply(names(gene_families), function(fn) {
    idx  <- which(grepl(gene_families[[fn]], sym))
    if (!length(idx)) return(NULL)
    rule <- family_exon_rule[[fn]]
    sel  <- lapply(idx, function(i) {
      g <- pick(exg[[i]], rule); mcols(g) <- NULL
      mcols(g)$family <- fn; mcols(g)$gene <- sym[i]; g
    })
    do.call(c, sel)
  })
  do.call(c, Filter(Negate(is.null), parts))
}

# --- Significance thresholds: SINGLE SOURCE OF TRUTH ---
# A peak is significant if |log2FoldChange| > lfc AND pvalue < p, per mark.
sig_thresholds <- list(
  "H3K4me3"  = list(lfc = 0.5, p = 0.05),
  "H3K9me2"  = list(lfc = 0.5, p = 0.05),
  "H3K9me3"  = list(lfc = 0.5, p = 0.10),
  "H4K20me3" = list(lfc = 0.5, p = 0.20),
  "H3K36me3" = list(lfc = 0.5, p = 0.20)
)

# Vectorised significance flag for a data.frame that has columns
# ChIP, log2FoldChange, pvalue. Marks not listed above are never significant.
# NA log2FC/pvalue -> FALSE (matches the previous dplyr::filter behaviour).
is_significant <- function(df) {
  lfc_cut <- vapply(as.character(df$ChIP),
                    function(m) if (!is.null(sig_thresholds[[m]])) sig_thresholds[[m]]$lfc else Inf,
                    numeric(1))
  p_cut   <- vapply(as.character(df$ChIP),
                    function(m) if (!is.null(sig_thresholds[[m]])) sig_thresholds[[m]]$p else -Inf,
                    numeric(1))
  ok <- !is.na(df$log2FoldChange) & !is.na(df$pvalue) &
        abs(df$log2FoldChange) > lfc_cut & df$pvalue < p_cut
  unname(ok)
}
