# ================================================================
# SHARED - cross-assay helpers (ChIP, GAM, ...)
# ----------------------------------------------------------------
# Sourced by every assay's config.R, right after shared/setup.R has ensured the
# dependencies. Holds only what is genuinely assay-INDEPENDENT:
#
#   - the core library attachments every assay needs
#   - figure colour palettes + scales  (SINGLE SOURCE OF TRUTH for colours)
#   - save_fig() / save_base_fig()
#   - ht_raster_device
#   - hyper_row()   - one hypergeometric enrichment row
#   - size_perm()   - size-matched permutation overlap test
#   - seqname style flips, mm39 binning, GAM bin-id parsing
#
# Assay-specific parameters (paths, thresholds, sample sheets, peak files) do
# NOT belong here - they live in that assay's config.R, which remains the single
# source of truth for its own analysis.
#
# NOTE ON SOURCING: this file is source()d from a config.R which is itself
# source()d at top level, so everything below is defined in the GLOBAL
# environment. That matters for save_fig()/save_base_fig() - see the comment
# there. Do not source this file inside a function.
# ================================================================

# --- Core libraries: the set every assay uses --------------------
# Assay-specific extras (DESeq2, ChIPseeker, ComplexHeatmap, TxDb, ...) are
# attached by that assay's config.R, not here.
suppressPackageStartupMessages({
  library(GenomicRanges)
  # GenomeInfoDb: seqlevelsStyle()<- / seqlevels(). GenomicRanges re-exported
  # these up to Bioc 3.22 but no longer does (3.23+), so load it explicitly -
  # the NCBI<->UCSC seqname flips throughout depend on it.
  library(GenomeInfoDb)
  library(rtracklayer)
  library(data.table)
  library(ggplot2)
  library(circlize)
  library(ltc)
})

# --- Figure savers: write BOTH png (quick view) and pdf (vector for papers), ---
# --- organised into per-analysis subfolders under <fig_dir>/<subdir>.        ---
# save_fig(): for ggplot objects.  save_base_fig(): for base-graphics / draw()
# renderers (e.g. ComplexHeatmap) - pass a zero-arg function that draws.
#
# WHY `fig_dir` IS NOT AN ARGUMENT: both functions look `fig_dir` up by name at
# CALL time. Because this file and each config.R are source()d into the global
# environment, these closures resolve `fig_dir` from there - so each assay gets
# its OWN fig_dir (ChIP/results/figures, GAM/results/figures) with no change at
# any call site, and MINUTE_RESULTS / GAM_RESULTS overrides keep working.
# Do NOT "fix" this by capturing fig_dir in a closure at definition time: this
# file is sourced BEFORE the config defines fig_dir, so that would capture
# nothing and every figure would land in the wrong place (or error).
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

# --- Heatmap rasterisation backend -------------------------------------------
# ComplexHeatmap rasterises any heatmap body over 2000 rows: it renders the body
# to a temporary image, reads it back, and embeds it as a single raster in the
# figure. That is a FIGURE-SIZE optimisation, not data I/O - the alternative is
# ~half a million vector rectangles in the PDF - but it does put an image device
# on the pipeline's critical path.
#
# Its default device is grDevices::png(), and it forces type = "cairo" whenever
# capabilities("cairo") is TRUE. On macOS without XQuartz that check LIES: R is
# built with cairo support, so the flag is TRUE, but cairo.so cannot load (it
# links to /opt/X11/lib). The device then fails to write the temp image with only
# a warning, and the read-back dies with:
#   Error in ... (temp_image) : unable to open ...heatmap_body_....png
# ragg's agg_png is self-contained (no X11/cairo), and ComplexHeatmap does not
# force a type on it. Pass this to every Heatmap() drawn.
ht_raster_device <- if (requireNamespace("ragg", quietly = TRUE)) "agg_png" else "png"

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

# --- Shared plotting constants (used across assays) ---
# Genotype labels are CANONICAL here: WT / HP1gKO. The GAM files on disk say
# "cKO" for HP1gKO; that translation happens in GAM/config.R and nowhere else,
# so no analysis code carries two names for one genotype.
geno_colors <- c(WT = gaby_cols[4], HP1gKO = gaby_cols[5])
theme_m     <- theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank())
mark_levels <- c("H3K4me3", "H3K36me3", "H3K9me2", "H3K9me3", "H4K20me3")

# Repeat classes (ChIP: MINUTE_2 heatmap annotation + scatters, MINUTE_3
# composition bars; GAM: per-bin repeat content).
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

# --- Helper: one hypergeometric enrichment row (OR + direction + p) for a
#     boolean annotation flag over ALL measured regions vs the significant
#     subset. Shared by ChIP MINUTE_2 (TAD/ChromHMM) and MINUTE_4 (repeats). ---
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

# --- Helper: SIZE-MATCHED permutation overlap test ---------------------------
# Is `query`'s overlap rate higher than a random set of regions drawn from
# `universe` with the SAME width distribution? Needed because our domains span
# kb-Mb while the features tested against them (KAP1 peaks, TSSs) are narrow, so
# a naive overlap test is confounded by region length: a big domain overlaps
# anything. Draws a null matched on width DECILES (nbins quantiles) rather than
# unmatched, which is what removes the length confound.
#
#   obs     : observed overlap fraction in the query set
#   query_w : query region widths
#   uni_w   : universe region widths
#   uni_ov  : logical/0-1 overlap indicator per universe region
#
# Originally defined inside ChIP/MINUTE_6_Intersection.R; hoisted here because
# the GAM compartment and HUSH-target stages need exactly the same null.
size_perm <- function(obs, query_w, uni_w, uni_ov, nbins = 20, nperm = 2000, seed = 42) {
  brks <- unique(quantile(uni_w, seq(0, 1, length.out = nbins + 1)))
  ub   <- cut(uni_w,   brks, include.lowest = TRUE)
  qb   <- cut(query_w, brks, include.lowest = TRUE)
  need <- table(qb); bins <- split(seq_along(uni_w), ub)
  set.seed(seed); perm <- numeric(nperm)
  for (i in seq_len(nperm)) {
    smp <- unlist(lapply(names(need), function(b)
      if (!is.na(b) && need[[b]] > 0 && length(bins[[b]]) > 0)
        sample(bins[[b]], need[[b]], replace = TRUE)))
    perm[i] <- mean(uni_ov[smp])
  }
  data.frame(obs_frac = obs, exp_frac = mean(perm), fold = obs / mean(perm),
             z = (obs - mean(perm)) / sd(perm),
             perm_p = (sum(perm >= obs) + 1) / (nperm + 1))
}

# --- Seqname style: the most common silent bug in this repo ------------------
# The ChIP counts/DESeq2 path is NCBI ("1", "2", ...); the annotation BEDs and
# every GAM file are UCSC ("chr1"). A mismatched flip yields ZERO overlaps,
# which reads downstream as a clean null result rather than an error - so any
# cross-assay join must assert its overlap rate rather than trust it.
# These wrappers also drop scaffolds, so the style flip and the standard-
# chromosome restriction always happen together (forgetting the latter is how
# unplaced contigs leak into bin grids).
to_ucsc <- function(gr) {
  suppressWarnings(GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC")
  GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = "coarse")
}
to_ncbi <- function(gr) {
  suppressWarnings(GenomeInfoDb::seqlevelsStyle(gr) <- "NCBI")
  GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = "coarse")
}

# --- Helper: fixed-width bin grid over mm39 ----------------------------------
# The shared coordinate system for cross-assay work: GAM is binned at 20 kb, so
# putting ChIP signal on the SAME grid is what makes the two joinable.
# Seqlengths come from the mm39 TxDb (requireNamespace, not library(), so this
# file does not force a TxDb attach on assays that never bin).
# Returns a UCSC-style GRanges; `chroms` defaults to the 19 autosomes, because
# the GAM differential calls are autosome-only.
mm39_bins <- function(binsize = 20000L, chroms = NULL) {
  if (!requireNamespace("TxDb.Mmusculus.UCSC.mm39.knownGene", quietly = TRUE))
    stop("mm39_bins() needs TxDb.Mmusculus.UCSC.mm39.knownGene")
  si <- GenomeInfoDb::seqinfo(TxDb.Mmusculus.UCSC.mm39.knownGene::TxDb.Mmusculus.UCSC.mm39.knownGene)
  if (is.null(chroms)) chroms <- paste0("chr", 1:19)
  si  <- si[intersect(chroms, GenomeInfoDb::seqlevels(si))]
  gr  <- GenomicRanges::tileGenome(si, tilewidth = binsize, cut.last.tile.in.chrom = TRUE)
  GenomeInfoDb::seqlevels(gr) <- GenomeInfoDb::seqlevelsInUse(gr)
  gr
}

# --- Helper: GAM bin ids ("chr19:0-20000") -----------------------------------
# The NPMI matrices label their rows/columns this way, so these convert between
# that label and a GRanges. NOTE the coordinate convention: the ids are
# 0-based half-open (BED-like), GRanges is 1-based inclusive, so parse adds 1 to
# the start and format subtracts it. Getting this wrong shifts every bin by one
# base - invisible in a plot, wrong in an overlap.
gam_bin_id <- function(gr) {
  paste0(as.character(GenomicRanges::seqnames(gr)), ":",
         GenomicRanges::start(gr) - 1L, "-", GenomicRanges::end(gr))
}
parse_gam_bin <- function(ids) {
  m <- regmatches(ids, regexec("^(.+):([0-9]+)-([0-9]+)$", ids))
  ok <- lengths(m) == 4L
  if (!all(ok)) stop("unparseable GAM bin id(s): ",
                     paste(utils::head(ids[!ok], 3), collapse = ", "))
  GenomicRanges::GRanges(
    vapply(m, `[`, "", 2L),
    IRanges::IRanges(start = as.numeric(vapply(m, `[`, "", 3L)) + 1,
                     end   = as.numeric(vapply(m, `[`, "", 4L))))
}
