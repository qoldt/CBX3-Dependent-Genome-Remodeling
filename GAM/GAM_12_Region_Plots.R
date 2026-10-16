# ================================================================
# GAM_12 - Region plots: contact map + hotspots + histone signal
# ----------------------------------------------------------------
#     cd GAM && Rscript GAM_12_Region_Plots.R          (after GAM_1)
#
# The composite figure this project actually needs: a differential contact map
# over a region, with the differential-contact hotspots and the ChIP log2FC
# tracks aligned underneath, in the style of the collaborators' own slides
# (CA1 HP1g WT/cKO, 12.10.2026 - chr3 30-85 Mb, chr11 100-125 Mb, chr16 full).
#
# The contact panel is drawn as a ROTATED TRIANGLE rather than a square matrix:
# that is what lets the x axis be genomic coordinate, so the map, the hotspots
# and the marks share one axis and can be read off each other. Height is contact
# separation, capped at gam12_maxdist - 5 Mb by default, matching the window the
# delivered hotspot track integrates over, so the two panels describe the same
# contacts.
#
# PARAMETERS FROM THE COLLABORATORS (slides, 12.10.2026):
#   20 kb resolution; full-chromosome distance for the maps; differential
#   contacts |delta-z| > 2 with max NPMI > 0.1 (justified by their saturation
#   analysis - NPMI distributions differ in the left tail but are concordant
#   above 0.1); hotspots over 0-5 Mb; 50-75 differential contacts per window as
#   the hotspot cutoff. delta-z = z(WT) - z(cKO), so positive = stronger in WT =
#   LOST in the KO.
#
# TISSUE: the ChIP is Emx1-Cre cortex AND hippocampus, so predominantly
# pyramidal neurons; the GAM is CA1 pyramidal cells specifically. Not an exact
# match, but the same broad cell type - the best available pairing.
# ================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(patchwork)
  library(TxDb.Mmusculus.UCSC.mm39.knownGene); library(org.Mm.eg.db) })

# Full contact distance by default, as the collaborators do ("Full chrom
# distance for contact maps"). Capping at 5 Mb on a 90 Mb chromosome forces a
# 40:1 sliver that no panel height can draw honestly - that was the distortion.
# The 5 Mb hotspot ceiling is drawn as a reference line instead.
gam12_maxdist   <- NULL
gam12_hotline   <- 5e6     # where the hotspot integration window ends

# RAW NPMI panels: blue - white - yellow - orange - red. NPMI runs [-1, 1] with
# 0 meaning independence, so white belongs at zero rather than at the bottom of
# a sequential ramp.
npmi_cols <- c("#2166ac", "#67a9cf", "#d1e5f0", "#ffffff",
               "#fff0b3", "#fdae61", "#e8602c", "#a50026")
# DIFFERENCE panel: a plain two-sided blue/red. A difference has no meaningful
# intermediate colour the way NPMI's own scale does, and the yellow/orange band
# reads as a third category rather than as "more".
diff_cols <- c("#2166ac", "#6ba3cd", "#ffffff", "#d6604d", "#b2182b")
gam12_hot_cut   <- 50      # slides: 50-75 differential contacts per window

# Gene models for the annotation track, loaded once.
GENES <- local({
  g <- genes(TxDb.Mmusculus.UCSC.mm39.knownGene, single.strand.genes.only = TRUE)
  sym <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, keys = names(g),
           column = "SYMBOL", keytype = "ENTREZID"))
  mcols(g)$sym <- ifelse(is.na(sym), names(g), sym)
  g
})

g1 <- file.path(rds_dir, "gam_bins_annotated.rds")
if (!file.exists(g1)) stop("GAM_1 output not found:\n  ", g1,
  "\nRun it first:  cd GAM && Rscript GAM_1_Bins_and_Signal.R", call. = FALSE)
X <- readRDS(g1)
A <- as.data.table(X$annot); setorder(A, chr, start)
for (mk in colnames(X$l2fc)) A[[mk]] <- X$l2fc[match(A$bin_id, rownames(X$l2fc)), mk]

# z-score within each genomic distance (each diagonal). This is the
# collaborators' normalisation, and it must be computed on the WHOLE chromosome
# matrix, not on a displayed sub-region: a 2 Mb window leaves ~90 values per
# diagonal, so a region-local z-score is a different, much noisier quantity than
# the chromosome-wide one their delta-z column carries.
# Walking diagonals is O(n^2); masking the matrix once per distance is O(n^3)
# and does not finish on a full chromosome.
zbydist <- function(M) {
  n <- nrow(M)
  # One grouped pass over every cell, keyed by diagonal offset. A loop over
  # diagonals re-allocates an index matrix thousands of times and does not
  # finish on a full chromosome; this is a single data.table grouping.
  dt <- data.table(v = as.vector(M),
                   d = as.vector(abs(.row(dim(M)) - .col(dim(M)))))
  dt[, z := NA_real_]
  dt[is.finite(v), z := {
    sdv <- stats::sd(v)
    if (.N < 10L || !is.finite(sdv) || sdv < 1e-9) NA_real_ else (v - mean(v)) / sdv
  }, by = d]
  matrix(dt$z, n, n, dimnames = dimnames(M))
}

# ---------------------------------------------------------------
# gam_region_plot(chrom, from, to, marks)
#   chrom  "chr11"; from/to in bp (NULL = whole chromosome)
#   marks  which ChIP tracks to stack, in order
# Returns a patchwork object; save with save_fig().
# ---------------------------------------------------------------
gam_region_plot <- function(chrom, from = NULL, to = NULL,
                            marks = mark_levels,     # all five, in the project order
                            maxdist = gam12_maxdist, hot_cut = gam12_hot_cut,
                            smooth_map = TRUE, show_raw = TRUE,
                            diff_mode = c("deltaz", "raw")) {
  diff_mode <- match.arg(diff_mode)

  reg <- A[chr == chrom]
  if (is.null(from)) from <- min(reg$start)
  if (is.null(to))   to   <- max(reg$end)
  reg <- reg[start >= from & end <= to]
  if (nrow(reg) < 10) stop("region too small: ", chrom, ":", from, "-", to, call. = FALSE)
  if (is.null(maxdist)) maxdist <- to - from      # full distance: a true triangle

  # --- contact difference over the region ------------------------------------
  raw_tri <- function(M, pos, ttl) {
    # Raw NPMI on a sequential scale, 0 to the 99th percentile of THIS panel -
    # the convention in Winick-Ng Fig. 1c. Shown per genotype so the structure
    # is visible before any differencing.
    mm <- nrow(M)
    q <- CJ(i = seq_len(mm), j = seq_len(mm))[j > i]
    q[, d := pos[j] - pos[i]]
    q <- q[d <= maxdist]
    q[, `:=`(x = (pos[i] + pos[j]) / 2e6, y = d / 2e6, v = M[cbind(i, j)])]
    q <- q[is.finite(v)]
    cap <- as.numeric(quantile(abs(q$v), .99, na.rm = TRUE))
    ggplot(q, aes(x, y, fill = pmax(-cap, pmin(v, cap)))) +
      geom_tile(width = tw, height = tw) +
      scale_fill_gradientn(colours = npmi_cols, limits = c(-cap, cap),
                           name = "NPMI", guide = guide_colourbar(barheight = 3)) +
      scale_y_continuous(expand = c(0, 0), name = "sep (Mb)") +
      scale_x_continuous(expand = c(0, 0)) +
      coord_cartesian(xlim = c(from, to) / 1e6, expand = FALSE) +
      labs(subtitle = ttl) +
      theme_m + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                      plot.subtitle = element_text(size = 9),
                      plot.margin = margin(2, 4, 1, 4))
  }

  mats <- lapply(c(WT = "WT", KO = "HP1gKO"), function(g) {
    m <- load_npmi(g, chrom)
    zm <- if (identical(diff_mode, "deltaz")) zbydist(m) else NULL
    ids <- intersect(reg$bin_id, rownames(m))
    list(raw = m[ids, ids, drop = FALSE],
         z   = if (is.null(zm)) NULL else zm[ids, ids, drop = FALSE])
  })
  ids <- intersect(rownames(mats$WT$raw), rownames(mats$KO$raw))
  W0 <- mats$WT$raw[ids, ids]; K0 <- mats$KO$raw[ids, ids]
  # DIFFERENCE MODE.
  #   "deltaz" reproduces the collaborators' pipeline - each matrix z-scored
  #     within each genomic distance (chromosome-wide, see zbydist above), then
  #     subtracted. Raw NPMI falls ~6x from 20 kb to 20 Mb, so a raw difference
  #     is dominated by the near-diagonal and long-range structure is invisible;
  #     delta-z puts every separation on the same footing.
  #   "raw" keeps KO - WT in NPMI units, the right scale when the question is
  #     absolute contact frequency.
  # SIGN: they compute z(WT) - z(cKO), so their positive means LOST in the KO.
  # This repo uses positive = MORE in the KO throughout, so our map is their map
  # inverted. Stated on the panel so it is not read backwards.
  D <- if (identical(diff_mode, "deltaz"))
         mats$KO$z[ids, ids] - mats$WT$z[ids, ids]
       else K0 - W0
  pos <- as.numeric(sub("-.*", "", sub(".*:", "", ids)))
  rm(mats); gc(FALSE)

  # AGGREGATE before drawing. Subsampling every k-th bin aliases: the rotated
  # lattice beats against the raster grid and the map reads as diagonal streaks
  # rather than domains. Block-averaging to a display bin keeps every value and
  # produces a smooth map. Target ~350 display columns across the region.
  n <- length(ids)
  fac <- max(1L, round(n / 350))
  dbin <- fac * gam_binsize
  grp <- ((seq_len(n) - 1L) %/% fac) + 1L
  if (fac > 1L) {
    Z <- D; nn <- is.finite(Z); Z[!nn] <- 0
    S1 <- rowsum(Z, grp); N1 <- rowsum(nn * 1, grp)
    S2 <- t(rowsum(t(S1), grp)); N2 <- t(rowsum(t(N1), grp))
    D <- S2 / pmax(N2, 1); D[N2 == 0] <- NA
    agg1 <- function(M) { Zz <- M; nz <- is.finite(Zz); Zz[!nz] <- 0
      a1 <- rowsum(Zz, grp); b1 <- rowsum(nz * 1, grp)
      a2 <- t(rowsum(t(a1), grp)); b2 <- t(rowsum(t(b1), grp))
      R <- a2 / pmax(b2, 1); R[b2 == 0] <- NA; R }
    W0 <- agg1(W0); K0 <- agg1(K0)
    pos <- as.numeric(tapply(pos, grp, mean))
    rm(Z, nn, S1, N1, S2, N2); gc(FALSE)
  }
  # Light 2D smooth of the DIFFERENCE map. Differencing two noisy matrices
  # leaves a diagonal streak texture that reads as structure but is not; a 3x3
  # box filter suppresses it while preserving domain-scale features. Display
  # only - every statistic elsewhere uses the unsmoothed values. (The Pombo
  # 2026 preprint applies a Gaussian filter for the same reason.)
  if (smooth_map && nrow(D) > 9) {
    Z <- D; nn <- is.finite(Z); Z[!nn] <- 0
    k <- 1L; mm <- nrow(Z)
    cs <- function(M) { M <- apply(M, 2, cumsum); rbind(0, M) }
    box <- function(M) {
      P <- cs(M); idx <- seq_len(mm)
      lo <- pmax(1L, idx - k); hi <- pmin(mm, idx + k)
      R <- (P[hi + 1L, , drop = FALSE] - P[lo, , drop = FALSE])
      t(apply(R, 1, function(r) { Q <- c(0, cumsum(r))
        (Q[hi + 1L] - Q[lo]) }))
    }
    num <- box(Z); den <- box(nn * 1)
    D <- num / pmax(den, 1); D[den == 0] <- NA
    rm(Z, nn, num, den); gc(FALSE)
  }
  m <- nrow(D)
  ij <- CJ(i = seq_len(m), j = seq_len(m))[j > i]
  ij[, d := pos[j] - pos[i]]
  ij <- ij[d <= maxdist]
  ij[, `:=`(x = (pos[i] + pos[j]) / 2e6, y = d / 2e6, v = D[cbind(i, j)])]
  ij <- ij[is.finite(v)]
  lim <- as.numeric(quantile(abs(ij$v), 0.99, na.rm = TRUE))
  tw <- dbin / 1e6   # tiles overlap 2x on the rotated lattice, which fills it

  p_map <- ggplot(ij, aes(x, y, fill = pmax(-lim, pmin(lim, v)))) +
    geom_tile(width = tw, height = tw) +
    scale_fill_gradientn(colours = diff_cols, limits = c(-lim, lim), name = if (identical(diff_mode, "deltaz")) "KO - WT\ndelta-z" else "KO - WT\nNPMI",
                         guide = guide_colourbar(barheight = 4)) +
    { if (gam12_hotline / 2 < maxdist / 2)
        geom_hline(yintercept = gam12_hotline / 2e6, linetype = 2, linewidth = .3,
                   colour = "grey30") } +
    { if (gam12_hotline / 2 < maxdist / 2)
        annotate("text", x = from / 1e6, y = gam12_hotline / 2e6, hjust = -0.04, vjust = -0.5,
                 size = 2.3, colour = "grey30", label = "max hotspot distance") } +
    scale_y_continuous(expand = c(0, 0), name = "separation (Mb)") +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(from, to) / 1e6, expand = FALSE) +
    labs(subtitle = if (identical(diff_mode, "deltaz"))
           "delta-z, cKO - WT  (red = more contact in cKO; note their slides use the opposite sign)"
         else "cKO - WT NPMI  (red = more contact in cKO)") +
    theme_m + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                    plot.subtitle = element_text(size = 9),
                    plot.margin = margin(2, 4, 1, 4))

  # --- hotspots as a heatmap: one row per direction ---------------------------
  # Two rows read better than back-to-back bars at this width: the eye picks out
  # contiguous patches of colour, which is what a hotspot IS, rather than trying
  # to integrate a spiky bar profile.
  hm <- melt(reg[, .(mid = (start + end) / 2e6,
                     `lost in cKO` = n_contacts_lost,
                     `gained in cKO` = n_contacts_gained)],
             id.vars = "mid", variable.name = "dirn", value.name = "n")
  hm[, dirn := factor(dirn, levels = c("gained in cKO", "lost in cKO"))]
  hcap <- as.numeric(quantile(c(reg$n_contacts_lost, reg$n_contacts_gained), 0.99))
  p_hot <- ggplot(hm, aes(mid, dirn, fill = pmin(n, hcap))) +
    geom_tile(width = gam_binsize / 1e6, height = 0.92) +
    scale_fill_gradientn(colours = c("white", "#f7c6b0", "#e06b4a", "#9e2a1f"),
                         name = "diff.\ncontacts",
                         guide = guide_colourbar(barheight = 3.2)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_cartesian(xlim = c(from, to) / 1e6) +
    labs(y = NULL) +
    theme_m + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                    panel.grid = element_blank(), plot.margin = margin(1, 4, 1, 4))

  # --- ChIP tracks ------------------------------------------------------------
  # Mark tracks: a smoothed profile, not one bar per 20 kb bin. Fill is by the
  # SIGN of log2FC - red below zero is signal lost in the KO, blue above zero is
  # gained. (The previous version coloured by "above or below the genome median",
  # which put two different greys/reds on values that were both negative and was
  # unreadable.) The dashed line is that mark's genome-wide median, so you can
  # still see where a region is worse than typical.
  # Lighter smoothing than before: enough to read as a profile, not so much
  # that real features are flattened. ~1/600 of the region, floor of 3 bins.
  smw <- max(3L, round(nrow(reg) / 600))
  mk_panel <- function(mk, last = FALSE, show_key = FALSE) {
    # Keep every bin, including those with no value: dropping them lets
    # geom_area draw a straight wedge across the gap, which looks like real
    # signal. Smoothing ignores the gaps; plotting leaves them as breaks.
    d <- reg[, .(mid = (start + end) / 2e6, v = get(mk))]
    ok <- is.finite(d$v)
    sm <- rep(NA_real_, nrow(d))
    if (sum(ok) > smw) {
      f <- as.numeric(stats::filter(ifelse(ok, d$v, 0), rep(1, smw), sides = 2))
      k <- as.numeric(stats::filter(ok * 1, rep(1, smw), sides = 2))
      sm <- ifelse(is.finite(k) & k > 0, f / k, NA_real_)
    }
    d[, v := ifelse(ok, ifelse(is.na(sm), v, sm), NA_real_)]
    gm <- median(A[[mk]], na.rm = TRUE)
    rng <- range(c(d$v, 0, gm), na.rm = TRUE)
    ggplot(d, aes(mid, v)) +
      # Two ribbons clipped at zero, NOT geom_area(fill = v < 0): mapping fill to
      # the sign splits the data into two groups, and each group's polygon then
      # joins its own points straight across the other's territory - which drew
      # long straight wedges wherever positive values were sparse.
      # The fill diverges from THAT MARK'S GENOME MEDIAN, not from zero. With a
      # near-genome-wide loss (H4K20me3 median -0.497, ~100% of bins down) a
      # zero baseline paints the whole track solid red and hides everything; the
      # informative quantity is deviation from the typical loss. The y axis
      # still shows absolute log2FC, and the solid line marks the median, so
      # nothing is hidden - only the fill reference changes.
      geom_ribbon(aes(ymin = gm, ymax = pmin(v, gm), fill = "below"), na.rm = TRUE) +
      geom_ribbon(aes(ymin = gm, ymax = pmax(v, gm), fill = "above"), na.rm = TRUE) +
      geom_hline(yintercept = gm, linewidth = .35, colour = "grey30") +
      { if (rng[1] < 0 && rng[2] > 0)
          geom_hline(yintercept = 0, linewidth = .3, colour = "grey65", linetype = 3) } +
      scale_fill_manual(values = c(below = "#c0392b", above = "#2166ac"),
                        breaks = c("below", "above"),
                        labels = c("below mark's genome median", "above it"),
                        name = NULL) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(n.breaks = 3, expand = expansion(mult = 0.04)) +
      coord_cartesian(xlim = c(from, to) / 1e6, ylim = rng, expand = FALSE) +
      labs(y = mk, x = if (last) sprintf("%s (Mb)", chrom) else NULL) +
      theme_m + theme(legend.position = if (show_key) "right" else "none",
                      legend.key.size = unit(9, "pt"), legend.text = element_text(size = 7),
                      plot.margin = margin(1, 4, if (last) 2 else 1, 4),
                      axis.title.y = element_text(size = 8),
                      axis.text.y = element_text(size = 6.5),
                      axis.title.x = if (last) element_text() else element_blank(),
                      axis.text.x  = if (last) element_text() else element_blank())
  }
  panels <- lapply(seq_along(marks), function(i)
                   mk_panel(marks[i], FALSE, show_key = (i == 1L)))

  # --- gene track -------------------------------------------------------------
  gg <- GENES[as.character(seqnames(GENES)) == chrom &
              end(GENES) > from & start(GENES) < to]
  gd <- data.table(x1 = pmax(start(gg), from) / 1e6, x2 = pmin(end(gg), to) / 1e6,
                   sym = mcols(gg)$sym, strand = as.character(strand(gg)))
  if (nrow(gd)) {
    setorder(gd, x1)
    gd[, row := (seq_len(.N) - 1L) %% 4L]        # stack to avoid overlap
    # label only the longer genes, or a cluster gets unreadable
    gd[, w := x2 - x1]
    span <- (to - from) / 1e6
    # label the widest few per stacking row, and only genes wide enough to carry
    # a label - the Pcdh cluster otherwise prints 40 overlapping names
    gd[, lab := NA_character_]
    gd[w >= span * 0.015, lab := sym]
    gd[, rk := frank(-w, ties.method = "first")]
    gd[rk > 10, lab := NA_character_]
  } else gd <- data.table(x1 = numeric(0), x2 = numeric(0), row = integer(0), lab = character(0))
  p_gene <- ggplot(gd) +
    geom_segment(aes(x = x1, xend = x2, y = row, yend = row), linewidth = 2.2,
                 colour = "#3b5566", lineend = "butt") +
    { if (any(!is.na(gd$lab)))
        geom_text(data = gd[!is.na(lab)], aes((x1 + x2) / 2, row, label = lab),
                  size = 2.1, vjust = -0.85, colour = "grey20") } +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_reverse(expand = expansion(add = 0.8)) +
    coord_cartesian(xlim = c(from, to) / 1e6) +
    labs(y = "genes", x = sprintf("%s (Mb)", chrom)) +
    theme_m + theme(axis.text.y = element_blank(), panel.grid = element_blank(),
                    axis.title.y = element_text(size = 8),
                    plot.margin = margin(1, 4, 4, 4))

  # The map's panel height is set from the region's own aspect ratio, so the
  # triangle is drawn to scale. coord_fixed() would also do that, but patchwork
  # still reserves the full row height for the panel and leaves a band of
  # whitespace under a wide region. A floor keeps very wide views legible.
  map_h <- max(0.85, min(3.2, (maxdist / 2) / (to - from) * 11))
  tri <- if (show_raw) list(raw_tri(W0, pos, "WT"), raw_tri(K0, pos, "HP1gKO"), p_map)
         else list(p_map)
  # The raw panels are reference, the difference is the point, and the mark
  # panels need real height or their axis titles collapse into each other.
  tri_h <- if (show_raw) c(map_h * 0.72, map_h * 0.72, map_h) else map_h
  wrap_plots(c(tri, list(p_hot), panels, list(p_gene)), ncol = 1,
             heights = c(tri_h, 0.42, rep(0.78, length(panels)), 0.85)) +
    plot_annotation(
      title = sprintf("%s:%.2f-%.2f Mb", chrom, from / 1e6, to / 1e6),
      theme = theme(plot.caption = element_text(hjust = 0, size = 7),
                    plot.title = element_text(face = "bold", size = 13)),
      caption = paste0(
      "GAM: CA1 pyramidal cells, HP1g WT/cKO, 20 kb, |delta-z|>2, max NPMI>0.1; hotspots over 0-5 Mb. ",
      "ChIP: Emx1-Cre cortex+hippocampus (pyramidal neurons), per-bin log2FC (KO/WT).\n",
      "Contact map block-averaged and 3x3 smoothed for display; mark tracks running-mean smoothed. ",
      "All statistics elsewhere use unsmoothed values. Mark fill diverges from each mark's genome median (solid line); ",
      "dotted line is zero where in range."))
}

# --- the collaborators' own example regions, plus our blocks -----------------
regions <- list(
  list("chr3",  30e6,  85e6, "chr3_30-85Mb"),
  list("chr11", 100e6, 125e6, "chr11_100-125Mb"),
  list("chr11", 1e6,   25e6,  "chr11_1-25Mb_lostblock"),
  list("chr16", NULL,  NULL,  "chr16_full"),
  list("chr18", 36.6e6, 38.4e6, "chr18_Pcdh"),
  list("chr18", 35.5e6, 39.5e6, "chr18_Pcdh_context")
)
for (r in regions) {
  nm <- r[[4]]
  message("Plotting ", nm)
  p <- tryCatch(gam_region_plot(r[[1]], r[[2]], r[[3]]),
                error = function(e) { message("  skipped: ", conditionMessage(e)); NULL })
  if (!is.null(p)) save_fig(p, nm, subdir = "regions", width = 11, height = 14)
}
cat("\nGAM_12 complete. gam_region_plot() is reusable for any region.\n")
