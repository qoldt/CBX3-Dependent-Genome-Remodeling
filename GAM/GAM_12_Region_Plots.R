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

# RAW NPMI panels: blue (low) - white - yellow - orange - red (high), spanning
# 0 to roughly the 99th percentile. Matching the collaborators' own figures:
# their scale bar runs 0.0 to ~0.35 with blue at the BOTTOM, i.e. it is a
# sequential ramp over the positive range, NOT a diverging scale centred on
# zero. Centring at zero (an earlier attempt here) pushes nearly every value
# into the warm half and flattens the structure, because almost all NPMI in a
# cis matrix is positive.
npmi_cols <- c("#3a6daa", "#73a8cd", "#c8e0ee", "#ffffff",
               "#ffe9a8", "#fdc266", "#f08a3c", "#cf4b28")
# DIFFERENCE panel: plain two-sided blue/red on a FIXED symmetric limit.
diff_cols <- c("#2166ac", "#6ba3cd", "#ffffff", "#d6604d", "#b2182b")
# Fixed delta-z limit, again matching their figures (scale bar -4 to +4). A
# data-driven limit (we used the 99th percentile of |delta-z|) rescales every
# panel to its own extremes, saturates the whole map and makes regions
# incomparable; a fixed limit leaves the bulk pale so the strong contacts show.
gam12_dz_lim    <- 4
gam12_npmi_qcap <- 0.99
gam12_hot_cut   <- 50      # slides: 50-75 differential contacts per window
gam12_count_mode <- "net"  # diff-contact COUNT panels: "net" / "split" / "both"

# Row heights in inches. Absolute, so the layout is the same whatever the region.
gam12_h_hot    <- 0.75
gam12_h_hotmap <- 0.55
gam12_h_hotcall<- 0.50
gam12_h_lad    <- 0.22
gam12_h_mark   <- 0.70
gam12_h_gene   <- 0.80
gam12_h_chrome <- 1.70   # title, caption, x axis, inter-panel margins

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
                            fig_width = 11, side_allow = 2.6,
                            diff_mode = c("deltaz", "raw")) {
  diff_mode <- match.arg(diff_mode)

  reg <- A[chr == chrom]
  if (is.null(from)) from <- min(reg$start)
  if (is.null(to))   to   <- max(reg$end)
  reg <- reg[start >= from & end <= to]
  if (nrow(reg) < 10) stop("region too small: ", chrom, ":", from, "-", to, call. = FALSE)
  if (is.null(maxdist)) maxdist <- to - from      # full distance: a true triangle

  # --- contact difference over the region ------------------------------------
  raw_tri <- function(M, pos, ttl, cap = NULL) {
    # Raw NPMI on a sequential scale, 0 to the 99th percentile - the convention
    # in Winick-Ng Fig. 1c. Shown per genotype so the structure is visible
    # before any differencing. `cap` is passed in and SHARED between the two
    # genotypes: a per-panel cap put WT and HP1gKO on different scales, which
    # defeats the whole point of showing them one above the other.
    mm <- nrow(M)
    q <- CJ(i = seq_len(mm), j = seq_len(mm))[j > i]
    q[, d := pos[j] - pos[i]]
    q <- q[d <= maxdist]
    q[, `:=`(x = (pos[i] + pos[j]) / 2e6, y = d / 2e6, v = M[cbind(i, j)])]
    q <- q[is.finite(v)]
    if (is.null(cap)) cap <- as.numeric(quantile(q$v, gam12_npmi_qcap, na.rm = TRUE))
    ggplot(q, aes(x, y, fill = pmin(pmax(v, 0), cap))) +
      geom_tile(width = tw, height = tw) +
      scale_fill_gradientn(colours = npmi_cols, limits = c(0, cap),
                           name = "NPMI", guide = guide_colourbar(barheight = grid::unit(2.4, "cm"))) +
      scale_y_continuous(expand = c(0, 0), name = "sep (Mb)") +
      scale_x_continuous(expand = c(0, 0)) +
      coord_fixed(ratio = 1, xlim = c(from, to) / 1e6, expand = FALSE) +
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
  lim <- if (identical(diff_mode, "deltaz")) gam12_dz_lim
         else as.numeric(quantile(abs(ij$v), 0.99, na.rm = TRUE))
  tw <- dbin / 1e6   # tiles overlap 2x on the rotated lattice, which fills it

  p_map <- ggplot(ij, aes(x, y, fill = pmax(-lim, pmin(lim, v)))) +
    geom_tile(width = tw, height = tw) +
    scale_fill_gradientn(colours = diff_cols, limits = c(-lim, lim), name = if (identical(diff_mode, "deltaz")) "KO - WT\ndelta-z" else "KO - WT\nNPMI",
                         guide = guide_colourbar(barheight = grid::unit(2.8, "cm"))) +
    { if (gam12_hotline / 2 < maxdist / 2)
        geom_hline(yintercept = gam12_hotline / 2e6, linetype = 2, linewidth = .3,
                   colour = "grey30") } +
    { if (gam12_hotline / 2 < maxdist / 2)
        annotate("text", x = from / 1e6, y = gam12_hotline / 2e6, hjust = -0.04, vjust = -0.5,
                 size = 2.3, colour = "grey30", label = "diff. contact counting window") } +
    scale_y_continuous(expand = c(0, 0), name = "sep (Mb)") +
    scale_x_continuous(expand = c(0, 0)) +
    coord_fixed(ratio = 1, xlim = c(from, to) / 1e6, expand = FALSE) +
    labs(subtitle = if (identical(diff_mode, "deltaz"))
           "delta-z, cKO - WT  (red = more contact in cKO; note their slides use the opposite sign)"
         else "cKO - WT NPMI  (red = more contact in cKO)") +
    theme_m + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                    plot.subtitle = element_text(size = 9),
                    plot.margin = margin(2, 4, 1, 4))

  # --- hotspots: ONE net track, signed like the difference map above it -------
  # These counts ARE contacts: per 20 kb bin, how many significant differential
  # contacts (|delta-z| > 2, max NPMI > 0.1) that bin takes part in within 5 Mb,
  # counted separately for the two directions. Plotting the two directions as
  # separate rows asked the reader to subtract them by eye, so the default here
  # is the net:
  #
  #     net = gained in cKO - lost in cKO
  #
  # which carries exactly the sign convention of the delta-z triangle (KO - WT),
  # and takes its red and blue from the same palette. A bin that is red in the
  # map is red in this track. Magnitude is the bar height, so a bin with 60 in
  # each direction reads as "busy but balanced" (bar near zero), which a pair of
  # heatmap rows showed as two hot patches and a single heatmap row could not
  # show at all.
  #
  # gam12_count_mode = "split" restores the two-row version: the net hides
  # reciprocal bins, and that is occasionally the thing being asked about.
  # The per-direction heatmap, kept as its own track: the net above hides a bin
  # that is busy in BOTH directions (60 lost and 60 gained nets to zero), and
  # those reciprocal bins are exactly where contacts are being redistributed
  # rather than simply removed. Net answers "which way"; the heatmap answers
  # "how much traffic, each way".
  p_count_hm <- local({
    hm <- melt(reg[, .(mid = (start + end) / 2e6,
                       `lost in cKO` = n_contacts_lost,
                       `gained in cKO` = n_contacts_gained)],
               id.vars = "mid", variable.name = "dirn", value.name = "n")
    hm[, dirn := factor(dirn, levels = c("gained in cKO", "lost in cKO"))]
    hcap <- as.numeric(quantile(c(reg$n_contacts_lost, reg$n_contacts_gained), 0.99))
    if (!is.finite(hcap) || hcap <= 0) hcap <- 1
    ggplot(hm, aes(mid, dirn, fill = pmin(n, hcap))) +
      geom_tile(width = gam_binsize / 1e6, height = 0.92) +
      scale_fill_gradientn(colours = c("white", "#f7c6b0", "#e06b4a", "#9e2a1f"),
                           name = "diff. contacts\nper bin",
                           guide = guide_colourbar(barheight = grid::unit(2.4, "cm"))) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      coord_cartesian(xlim = c(from, to) / 1e6) +
      labs(y = NULL) +
      theme_m + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                      axis.text.y = element_text(size = 6.5),
                      panel.grid = element_blank(), plot.margin = margin(1, 4, 1, 4))
  })

  p_count_net <- local({
    hd <- reg[, .(mid = (start + end) / 2e6,
                  net = n_contacts_gained - n_contacts_lost)]
    # Symmetric limit so the two directions are drawn on one scale; without it a
    # region that mostly loses contacts makes its few gains look enormous.
    ncap <- as.numeric(quantile(abs(hd$net), 0.995, na.rm = TRUE))
    if (!is.finite(ncap) || ncap <= 0) ncap <- 1
    hd[, nc := pmax(-ncap, pmin(net, ncap))]
    ggplot(hd, aes(mid, nc, fill = nc > 0)) +
      geom_col(width = gam_binsize / 1e6) +
      geom_hline(yintercept = 0, linewidth = .25, colour = "grey40") +
      scale_fill_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "#2166ac"),
                        guide = "none") +
      scale_x_continuous(expand = c(0, 0)) +
      # Breaks at +/-ncap sit exactly on the panel edges (expand = 0) and collide
      # with the rows above and below, so the scale goes in the axis title and
      # only the zero line is labelled.
      scale_y_continuous(limits = c(-ncap, ncap), expand = c(0, 0), breaks = 0) +
      coord_cartesian(xlim = c(from, to) / 1e6) +
      labs(y = sprintf("net diff.\ncontacts (+/-%d)", round(ncap))) +
      theme_m + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                      axis.title.y = element_text(size = 7),
                      panel.grid.major.x = element_blank(),
                      plot.margin = margin(1, 4, 1, 4))
  })

  # --- hotspots: the CALL, not the count --------------------------------------
  # These are two different things and the earlier version of this figure
  # conflated them. The delivered Hotspots/ files are a dense per-bin COUNT -
  # every mappable bin has one, median 14 - so they are a differential-contact
  # density track, not a hit list. A hotspot is what you get after thresholding
  # that count: the collaborators' slides give 50-75 differential contacts per
  # window as the cutoff. gam12_hot_cut carries that value and until now was
  # accepted by this function and never used.
  #
  # Drawn as discrete blocks so it reads as a call. The repo also has a second,
  # smoothed definition (gam_hotspot_domains() in config.R: q90 + gap merging,
  # used by GAM_1 to emit domains); that one is deliberately NOT what is drawn
  # here, because these panels are showing the collaborators' own track under
  # the collaborators' own cutoff. Whichever is meant has to be said out loud.
  p_hotspot <- local({
    # ONE row, direction carried by colour. Only 26 bins genome-wide (0.1% of
    # hotspot bins) are called in both directions, so a single row loses almost
    # nothing - but those 26 get their own colour rather than letting whichever
    # is drawn second silently paint over the other.
    hs <- reg[, .(mid = (start + end) / 2e6,
                  l = n_contacts_lost   >= hot_cut,
                  g = n_contacts_gained >= hot_cut)]
    hs <- hs[l | g]
    hs[, call := fifelse(l & g, "both", fifelse(g, "gained in cKO", "lost in cKO"))]
    hs[, call := factor(call, levels = c("gained in cKO", "lost in cKO", "both"))]
    ggplot(hs, aes(mid, 1, fill = call)) +
      geom_tile(width = gam_binsize / 1e6, height = 1) +
      scale_fill_manual(values = c(`gained in cKO` = "#b2182b",
                                   `lost in cKO`   = "#2166ac",
                                   both            = "#5d3a7a"),
                        # drop unused levels: "both" is 0.1% of hotspot bins
                        # genome-wide, so on most regions keeping it would put an
                        # empty key in the legend of every figure.
                        name = sprintf("hotspot (>=%d)", hot_cut), drop = TRUE) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_cartesian(xlim = c(from, to) / 1e6) +
      labs(y = "hotspot") +
      theme_m + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                      axis.title.y = element_text(size = 7),
                      panel.grid = element_blank(), plot.margin = margin(1, 4, 1, 4))
  })

  # "both" (the default) stacks net count over per-direction count; either way
  # the hotspot CALL is a separate track below them.
  hot_panels <- switch(gam12_count_mode,
                       net   = list(p_count_net),
                       split = list(p_count_hm),
                       list(p_count_net, p_count_hm))
  hot_h <- switch(gam12_count_mode,
                  net   = gam12_h_hot,
                  split = gam12_h_hotmap,
                  c(gam12_h_hot, gam12_h_hotmap))
  hot_panels <- c(hot_panels, list(p_hotspot))
  hot_h      <- c(hot_h, gam12_h_hotcall)

  # --- LAD ribbon -------------------------------------------------------------
  # Constitutive LADs (Peric-Hupkes 2010, lifted to mm39). External to both our
  # genotypes and all five marks, so it is the one axis in this figure that
  # cannot be circular. Drawn as a thin ribbon rather than a track because the
  # only question asked of it is "is this stretch peripheral or not".
  lads <- tryCatch(load_lads(), error = function(e) NULL)
  p_lad <- if (is.null(lads)) NULL else {
    lg <- as.data.table(as.data.frame(
      IRanges::subsetByOverlaps(lads, GRanges(chrom, IRanges(from, to)))))
    lg <- lg[, .(s = pmax(start, from) / 1e6, e = pmin(end, to) / 1e6, class)]
    ggplot(lg) +
      geom_rect(aes(xmin = s, xmax = e, ymin = 0, ymax = 1, fill = class)) +
      scale_fill_manual(values = c(cLAD = "#53426b", ciLAD = "#cbbddc"),
                        name = NULL, drop = FALSE) +
      scale_x_continuous(expand = c(0, 0), limits = c(from, to) / 1e6) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(y = "LAD") +
      theme_m + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                      axis.title.y = element_text(size = 7),
                      panel.grid = element_blank(), plot.margin = margin(1, 4, 1, 4))
  }

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
                   mk_panel(marks[i], FALSE, show_key = (i == length(marks))))

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

  # SQUARE DATA POINTS. In a rotated triangle every matrix cell is a 45 degree
  # diamond, and it is square only when one Mb on x covers the same physical
  # distance as one Mb on y. coord_fixed(ratio = 1) is the only thing that
  # guarantees that; any layout that stretches the panel to fill an allocated
  # height shears the cells.
  #
  # coord_fixed alone is not enough inside patchwork: patchwork still reserves
  # the row height it was given, and the fixed panel then floats in it with a
  # band of whitespace. So the map rows are given an ABSOLUTE height in inches,
  # computed from the aspect the fixed coordinate system will adopt:
  #
  #     panel height = panel width * (y data range) / (x data range)
  #
  # where the y range is maxdist/2 (the triangle apex) and the x range is the
  # region width. With full contact distance (maxdist = region width) that ratio
  # is always 0.5, so a full-distance triangle is a true isoceles right triangle
  # at any scale - which is why the collaborators' full-chromosome maps are
  # undistorted however wide the region is.
  panel_w  <- max(2, fig_width - side_allow)      # inches of plotting area
  map_h_in <- panel_w * ((maxdist / 2) / (to - from))
  # One cap from the pooled values of both genotypes, so the two raw panels are
  # directly comparable and patchwork can collect their now-identical scales
  # into a single legend instead of drawing two.
  npmi_cap <- if (!show_raw) NULL else local({
    kk <- CJ(i = seq_len(nrow(W0)), j = seq_len(nrow(W0)))[j > i]
    kk <- kk[pos[j] - pos[i] <= maxdist]
    v  <- c(W0[cbind(kk$i, kk$j)], K0[cbind(kk$i, kk$j)])
    as.numeric(quantile(v[is.finite(v)], gam12_npmi_qcap, na.rm = TRUE))
  })
  tri <- if (show_raw) list(raw_tri(W0, pos, "WT", npmi_cap),
                            raw_tri(K0, pos, "HP1gKO", npmi_cap), p_map)
         else list(p_map)
  # The raw panels are reference, the difference is the point, and the mark
  # panels need real height or their axis titles collapse into each other.
  # Absolute inches for every map row (so coord_fixed is honoured exactly),
  # relative "null" units for the tracks (so they share what is left).
  n_tri   <- length(tri)
  # EVERY row gets an absolute height. Mixing "in" map rows with "null" track
  # rows made the layout non-deterministic: patchwork hands the null rows
  # whatever is left of the figure height, so shrinking the maps (a wide region,
  # or a band) inflated the tracks to fill the page and the maps came out as
  # slivers. With all rows in inches the figure height is just their sum, and a
  # map row is the same physical size whatever else is in the figure.
  lad_h   <- if (is.null(p_lad)) numeric(0) else gam12_h_lad
  track_h <- c(hot_h, lad_h, rep(gam12_h_mark, length(panels)), gam12_h_gene)
  heights <- grid::unit(c(rep(map_h_in, n_tri), track_h), "in")
  # guides = "collect" is essential here, not cosmetic. Without it every panel
  # keeps its own legend, and patchwork centres each one on ITS OWN ROW - so the
  # LAD legend (two keys) was being centred on a 0.22 in row and spilled over
  # the legends above and below it. Collected and top-justified, they stack.
  wrap_plots(c(tri, hot_panels, if (is.null(p_lad)) NULL else list(p_lad),
               panels, list(p_gene)), ncol = 1, heights = heights) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = sprintf("%s:%.2f-%.2f Mb", chrom, from / 1e6, to / 1e6),
      theme = theme(plot.caption = element_text(hjust = 0, size = 7),
                    plot.title = element_text(face = "bold", size = 13),
                    legend.spacing.y = grid::unit(12, "pt"),
                    legend.box.spacing = grid::unit(10, "pt")),
      caption = paste0(
      "GAM: CA1 pyramidal cells, HP1g WT/cKO, 20 kb, |delta-z|>2, max NPMI>0.1. ",
      "ChIP: Emx1-Cre cortex+hippocampus (pyramidal neurons), per-bin log2FC (KO/WT).\n",
      "Contact map block-averaged and 3x3 smoothed for display; mark tracks running-mean smoothed. ",
      "All statistics elsewhere use unsmoothed values. Mark fill diverges from each mark's genome median (solid line); ",
      "dotted line is zero where in range.\n",
      "Diff. contacts = count of significant differential contacts per 20 kb bin within 5 Mb (the dashed line on the map); ",
      "net = gained - lost, same sign convention as the difference map. Those are COUNTS, every bin has one. ",
      "Hotspot is the separate CALL on that count: bins at or above the collaborators' 50-75 per-window cutoff. ",
      "LADs: constitutive cLAD/ciLAD, Peric-Hupkes 2010 (GSE17051), lifted mm9->mm39.")) &
    theme(legend.justification = "top") -> out
  # Exact figure height: every row plus the title/caption/axis furniture.
  attr(out, "gam_map_in") <- list(n = n_tri, h = map_h_in,
                                  total = n_tri * map_h_in + sum(track_h) + gam12_h_chrome)
  out
}

# --- regions -----------------------------------------------------------------
# Each is rendered twice:
#   <name>        full contact distance - a true isoceles triangle, every scale
#                 comparable, but tall (three stacked maps at aspect 0.5).
#   <name>_band   distance capped at gam12_band, which keeps the cells square
#                 (the aspect just gets flatter) and fits on a screen.
#
# WHY 5 Mb FOR THE BAND. It is the window the delivered hotspot track integrates
# over, so the map and the hotspot rows below it then describe exactly the same
# contacts - anything visible in the band is in the track, and vice versa. It
# also keeps every band plot on one distance scale, so regions are comparable.
# The caveat is that several of our effects run past it: the chr16 block loses
# most at 5-14 Mb and chr11 PROX at 5-20 Mb. For those, a wider band is the
# right call; gam_region_plot(..., maxdist = 20e6) does it without touching
# anything else, and the cells stay square because coord_fixed does not care
# how flat the band is.
# BAND SIZING. A fixed band distance cannot work across region sizes, because
# what the eye reads is the panel's ASPECT, (maxdist / 2) / region width. At a
# fixed 5 Mb that aspect is 0.045 over a 55 Mb region (an unreadable sliver) and
# 1.39 over the 1.8 Mb Pcdh window - taller than the full triangle, and asking
# for 5 Mb separations that a 1.8 Mb window does not contain.
#
# So the band is specified the way it is actually judged: by how tall the map
# panel should be on the page. The distance follows from the geometry,
#
#     maxdist = 2 * (target height / panel width) * region width
#
# which keeps every band panel the same physical size. On the 24 Mb chr11 region
# this lands at ~5 Mb, i.e. it reproduces the hand-picked value that looked
# right, and it scales that same appearance to chr3 (~12 Mb) and chr16 (~21 Mb).
# Those larger bands are also the honest ones for those regions: the chr16 block
# operates at 5-14 Mb separations, which a 5 Mb band would have cut off.
gam12_band_h   <- 0.95   # inches: target height of a band map panel
gam12_band_min <- 2e6    # below this a band shows too few diagonals to read

# The same NULL-means-whole-chromosome rule gam_region_plot() uses, hoisted so
# the driver can size the band before building the plot.
region_span <- function(chrom, from = NULL, to = NULL) {
  reg <- A[chr == chrom]
  c(if (is.null(from)) min(reg$start) else from,
    if (is.null(to))   max(reg$end)   else to)
}

band_maxdist <- function(from, to, panel_w) {
  md <- 2 * (gam12_band_h / panel_w) * (to - from)
  max(md, gam12_band_min)
}

regions <- list(
  list("chr3",  30e6,  85e6, "chr3_30-85Mb"),
  list("chr11", 100e6, 125e6, "chr11_100-125Mb"),
  list("chr11", 1e6,   25e6,  "chr11_1-25Mb_lostblock"),
  list("chr16", NULL,  NULL,  "chr16_full"),
  list("chr16", 70e6,  90e6,  "chr16_73-87Mb_lostblock"),
  list("chr18", 36.6e6, 38.4e6, "chr18_Pcdh"),
  list("chr18", 35.5e6, 39.5e6, "chr18_Pcdh_context")
)
FIGW <- 11; SIDE <- 2.6
for (r in regions) {
  for (mode in c("full", "band")) {
    nm <- if (mode == "band") paste0(r[[4]], "_band") else r[[4]]
    md <- NULL
    if (mode == "band") {
      lim <- region_span(r[[1]], r[[2]], r[[3]])
      md  <- band_maxdist(lim[1], lim[2], max(2, FIGW - SIDE))
      # A "band" as tall as the full triangle is the full triangle. Emitting it
      # under a _band name would just be a duplicate file with a misleading name.
      if (md >= 0.9 * (lim[2] - lim[1])) {
        message("Skipping ", nm, ": region is ", round((lim[2] - lim[1]) / 1e6, 1),
                " Mb, the full triangle is already compact")
        next
      }
    }
    message("Plotting ", nm,
            if (is.null(md)) "" else sprintf(" (band %.1f Mb)", md / 1e6))
    p <- tryCatch(gam_region_plot(r[[1]], r[[2]], r[[3]], maxdist = md,
                                  fig_width = FIGW, side_allow = SIDE),
                  error = function(e) { message("  skipped: ", conditionMessage(e)); NULL })
    if (!is.null(p)) {
      a <- attr(p, "gam_map_in")
      h <- if (is.null(a)) 14 else a$total
      save_fig(p, nm, subdir = "regions", width = FIGW, height = h)
    }
  }
}
cat("\nGAM_12 complete. gam_region_plot() is reusable for any region;",
    "pass maxdist to set the band height.\n")
