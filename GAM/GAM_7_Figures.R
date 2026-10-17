# ================================================================
# GAM_7 - Figures: region plots and per-chromosome overviews
# ----------------------------------------------------------------
#     cd GAM && Rscript GAM_7_Figures.R              (after GAM_1, GAM_2, GAM_4)
#
# Merges the two stages that only RENDER - they compute no new statistics that
# anything else depends on, and both read the same two .rds handoffs.
# Previously GAM_12 and GAM_15.
#
#   A (was GAM_12)  gam_region_plot(): contact map, differential contacts,
#                   hotspot calls, LADs, the five marks and genes on one axis.
#                   Reusable for any region; the `regions` list at the end of the
#                   section is what gets rendered.
#   B (was GAM_15)  Per-chromosome characterisation of the five chromosomes with
#                   large effects, plus the all-19 ranking any single-chromosome
#                   claim has to be read against.
#
# export_igv.R stays a standalone utility - it is documented as being run on its
# own to refresh the IGV session, and folding it in would break that.
#
# Output file names and figure subdirectories are unchanged.
# ================================================================
source("config.R")
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(GenomicRanges); library(IRanges)
  library(patchwork); library(TxDb.Mmusculus.UCSC.mm39.knownGene); library(org.Mm.eg.db)
})





# ==============================================================================
# SECTION A - region plots                                 (was GAM_12)
# ==============================================================================

local({
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

  # HORIZONTAL LEFT-HAND LABELS. Rotated 90-degree axis titles force the reader to
  # tilt their head once per track, and with ten stacked tracks that is the main
  # friction in reading this figure. angle = 0 with hjust = 1 right-aligns them
  # into a clean column; labels are pre-wrapped so the column stays narrow.
  # patchwork aligns panels on the widest label, so this is set once and applies
  # to every row.
  ylab_h <- theme(axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5,
                                              size = 7, lineheight = 0.95))

  # Row heights in inches. Absolute, so the layout is the same whatever the region.
  gam12_h_hot    <- 0.75
  gam12_h_hotmap <- 0.55
  # The hotspot CALL is a categorical ribbon exactly like the LAD ribbon below
  # it - a presence/absence strip, not a magnitude - so it gets the same height.
  # Only the differential-contact COUNT track above needs room, because there the
  # bar height carries the number.
  # 0.26, not 0.22: with horizontal axis titles each ribbon row has to be tall
  # enough to seat its own label, and at 0.22 in the hotspot / CTCF / LAD labels
  # collided with each other.
  gam12_h_lad    <- 0.26
  gam12_h_hotcall<- 0.26
  gam12_h_mark   <- 0.70
  gam12_h_gene   <- 0.80
  gam12_h_ctcf   <- 0.26   # CTCF strip, same height as the LAD ribbon
  # Caption is ~10 wrapped lines now that the two mark-track references are stated
  # in full, so this has to cover more than it used to.
  gam12_h_chrome <- 2.30

  # MINIMUM Y SPAN of a mark panel, in log2FC. Without a floor the axis fits
  # whatever the region happens to contain, so a mark with almost no dynamic range
  # there has its noise stretched to the full panel height and reads as dramatic
  # structure. Marks with real signal already exceed this and are unaffected.
  gam12_mark_minspan <- 0.5

  # DISPLAY BIN FOR THE MARK TRACKS. The analysis grid is 20 kb, but a per-bin
  # log2FC at 20 kb is noise-dominated: after the x10 count scaling DESeq2 sees
  # integers of roughly 1-11 per bin, and the resulting unshrunk ratio has a
  # per-bin sd of 0.14-0.16 for the broad marks and up to 0.57 for H3K4me3 -
  # several times the local effects people try to read off these tracks.
  # Aggregating the DISPLAY to 100 kb gives ~5x the counts and roughly halves that
  # noise. Every statistic elsewhere still uses the 20 kb values.
  #
  # Shrinkage is deliberately NOT used for this: ashr collapses H3K36me3 and
  # H3K9me2 to exactly zero and pulls H4K20me3 from -0.497 to -0.212, because a
  # shrinkage prior assumes most effects are null - the same assumption as
  # median-of-ratios normalisation, which sizeFactors = 1 exists here to avoid.
  gam12_mark_binsize <- 100000
  gam12_mark_npoints <- 1200   # columns the mark track is resampled onto
  gam12_cov_warn     <- 0.75   # below this fraction of finite NPMI cells, warn on the figure

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
  # baseMean too: the mark tracks need it both to weight the aggregation and to
  # fade bins where there is no signal to take a ratio of.
  local({
    S <- as.data.table(as.data.frame(X$stats))
    for (mk in colnames(X$l2fc)) {
      d <- S[mark == mk]
      A[[paste0("bm_", mk)]] <<- d$baseMean[match(A$bin_id, d$bin_id)]
    }
  })

  # PER-GENOTYPE COVERAGE, which is what the mark tracks now draw. $signal is the
  # un-rounded per-bin, per-library coverage, so a genotype mean is just the mean
  # over that genotype's libraries - no weighting, because coverage IS the
  # quantity rather than a ratio built from one.
  # min/max across libraries give the replicate-range band: with 6 WT and 4 cKO
  # libraries (cKO rep2 excluded at QC) and no within-genotype test available,
  # the spread between replicates is the only honest indication of whether a gap
  # between the two means means anything.
  local({
    i <- match(A$bin_id, rownames(X$signal))
    for (mk in colnames(X$l2fc)) {
      for (gt in c("WT", "HP1gKO")) {
        cols <- grep(paste0("^", mk, "_", gt, "_rep"), colnames(X$signal))
        if (!length(cols)) next
        M <- X$signal[i, cols, drop = FALSE]
        A[[paste0("sig_", mk, "_", gt)]]     <<- rowMeans(M, na.rm = TRUE)
        A[[paste0("sig_", mk, "_", gt, "_lo")]] <<- matrixStats::rowMins(M, na.rm = TRUE)
        A[[paste0("sig_", mk, "_", gt, "_hi")]] <<- matrixStats::rowMaxs(M, na.rm = TRUE)
      }
    }
  })

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
        scale_y_continuous(expand = c(0, 0), name = "Genomic\ndistance\n(Mb)") +
        scale_x_continuous(expand = c(0, 0)) +
        coord_fixed(ratio = 1, xlim = c(from, to) / 1e6, expand = FALSE) +
        labs(subtitle = ttl) +
        theme_m + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                        axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5,
                                                    size = 7, lineheight = 0.95),
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

    # DATA COVERAGE, stated on the figure. A region whose bins are detected in few
    # nuclear profiles has mostly NA in the NPMI matrix, and an NA map renders as
    # white - indistinguishable from "no contacts here". That is this repo's
    # documented worst failure mode: an absence of measurement that reads as a
    # clean biological result. The vomeronasal cluster on chr7 is the live
    # example: 11% finite cells against 97% in its own flanks, so its apparently
    # empty contact map is missing data, not absent contact.
    npmi_cov <- mean(is.finite(W0[upper.tri(W0)]))
    if (npmi_cov < gam12_cov_warn)
      message(sprintf("  WARNING: only %.0f%% of NPMI cells are finite here - the map is mostly missing data",
                      100 * npmi_cov))

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
      scale_y_continuous(expand = c(0, 0), name = "Genomic\ndistance\n(Mb)") +
      scale_x_continuous(expand = c(0, 0)) +
      coord_fixed(ratio = 1, xlim = c(from, to) / 1e6, expand = FALSE) +
      labs(subtitle = if (identical(diff_mode, "deltaz"))
             "delta-z, cKO - WT  (red = more contact in cKO; note their slides use the opposite sign)"
           else "cKO - WT NPMI  (red = more contact in cKO)") +
      theme_m + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                      axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5,
                                                  size = 7, lineheight = 0.95),
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
        labs(y = sprintf("net diff.\ncontacts\n(+/-%d)", round(ncap))) +
        theme_m + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                        axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5,
                                                    size = 7, lineheight = 0.95),
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
                        axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 7),
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
                        axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 7),
                        panel.grid = element_blank(), plot.margin = margin(1, 4, 1, 4))
    }

    # --- ChIP tracks ------------------------------------------------------------
    # Mark tracks: a smoothed profile, not one bar per 20 kb bin.
    #
    # READ THE FILL CAREFULLY. It diverges from that mark's GENOME MEDIAN, not
    # from zero, so blue means "above the median" and NOT "increased in the KO".
    # For a mark whose median is negative, blue is still a loss - just a smaller
    # one than typical. This has already been misread once: H3K36me3 looks blue
    # over the Pcdh cluster, but its value there is -0.049 against a genome median
    # of -0.052, i.e. still down, and +0.004 off the genome is p = 0.93.
    #
    # The fill also carries no magnitude - +0.004 and +0.4 are the same blue - and
    # the y axis auto-scales per panel, so a mark with a narrow dynamic range has
    # its noise stretched to fill the panel. Hence the median value printed in the
    # axis label below: the reference is never left implicit.
    # Lighter smoothing than before: enough to read as a profile, not so much
    # that real features are flattened. ~1/600 of the region, floor of 3 bins.
    smw <- max(3L, round(nrow(reg) / 600))
    mk_panel <- function(mk, last = FALSE, show_key = FALSE) {
      # WT AND cKO COVERAGE OVERLAID, not log2 fold change.
      #
      # Why this replaced the log2FC panels. Those drew a ratio on a per-mark
      # symmetric axis, so every mark filled its panel to the same extent
      # whatever the size of its effect: H3K9me2 (genome median -0.18) looked
      # exactly as dramatic as H4K20me3 (-0.50). The encoding erased the one
      # comparison the figure most needs to support. A ratio also hides depth -
      # a 2x change on nothing is still nothing - and these are read as browser
      # tracks, where the expectation is coverage.
      #
      # What is drawn, per mark:
      #   pale band per genotype   min-max across that genotype's libraries
      #   solid line per genotype  mean across its libraries
      #   shaded gap between means red where cKO > WT, blue where cKO < WT,
      #                            the same convention as the contact maps
      # So height is signal level, the gap is the effect, its colour is the
      # direction, and the bands say whether the gap is larger than the spread
      # between replicates. A mark with a weak effect now looks weak: its two
      # lines sit on top of each other inside overlapping bands.
      #
      # Y SCALE IS PER MARK, not shared. H3K4me3 is punctate and high, H3K9me2
      # broad and low; one scale across marks would flatten the punctate ones to
      # a flat line. Comparison across marks is of the GAP relative to the bands,
      # which is scale-free, not of absolute height.
      cols <- paste0("sig_", mk, "_", rep(c("WT", "HP1gKO"), each = 3),
                     c("", "_lo", "_hi"))
      if (!all(cols %in% names(reg))) return(NULL)
      # c("start","end", cols) with with = FALSE selects COLUMNS by name. mget()
      # inside a data.table [ resolves names as variables in the calling frame,
      # not as columns, which is why it failed here.
      d <- reg[, c("start", "end", cols), with = FALSE]
      setnames(d, cols, c("wt", "wt_lo", "wt_hi", "ko", "ko_lo", "ko_hi"))

      # Aggregate to the display bin. Plain means: coverage is additive, and
      # unlike the log2FC version there is no near-empty-denominator problem to
      # weight away.
      fac <- max(1L, round(gam12_mark_binsize / gam_binsize))
      d[, grp := ((seq_len(.N) - 1L) %/% fac) + 1L]
      d <- d[, .(x0 = min(start) / 1e6, x1 = max(end) / 1e6,
                 wt = mean(wt, na.rm = TRUE), ko = mean(ko, na.rm = TRUE),
                 wt_lo = mean(wt_lo, na.rm = TRUE), wt_hi = mean(wt_hi, na.rm = TRUE),
                 ko_lo = mean(ko_lo, na.rm = TRUE), ko_hi = mean(ko_hi, na.rm = TRUE)),
             by = grp]

      fine <- local({
        mid <- (d$x0 + d$x1) / 2
        fx  <- seq(from / 1e6, to / 1e6, length.out = gam12_mark_npoints)
        itp <- function(y) { ok <- is.finite(y)
          if (sum(ok) < 2) return(rep(NA_real_, length(fx)))
          stats::approx(mid[ok], y[ok], xout = fx, rule = 1)$y }
        ppb  <- length(fx) / max(1, nrow(d))
        smw2 <- max(3L, as.integer(round(ppb / 2)))
        roll <- function(y) { g <- is.finite(y)
          num <- as.numeric(stats::filter(ifelse(g, y, 0), rep(1, smw2), sides = 2))
          den <- as.numeric(stats::filter(g * 1, rep(1, smw2), sides = 2))
          out <- ifelse(is.finite(den) & den > 0, num / den, NA_real_)
          ifelse(g & !is.finite(out), y, out) }
        data.table(x = fx, wt = roll(itp(d$wt)), ko = roll(itp(d$ko)),
                   wt_lo = roll(itp(d$wt_lo)), wt_hi = roll(itp(d$wt_hi)),
                   ko_lo = roll(itp(d$ko_lo)), ko_hi = roll(itp(d$ko_hi)))
      })
      ymax <- max(c(fine$wt_hi, fine$ko_hi, fine$wt, fine$ko), na.rm = TRUE)
      if (!is.finite(ymax) || ymax <= 0) ymax <- 1
      rng <- c(0, ymax * 1.04)

      ggplot(fine) +
        # replicate-range bands first, so the means read on top of them
        geom_ribbon(aes(x = x, ymin = wt_lo, ymax = wt_hi), fill = geno_colors[["WT"]],
                    alpha = .20, na.rm = TRUE) +
        geom_ribbon(aes(x = x, ymin = ko_lo, ymax = ko_hi), fill = geno_colors[["HP1gKO"]],
                    alpha = .20, na.rm = TRUE) +
        # the gap between the means, coloured by direction. Two ribbons each
        # clipped at the other mean, for the same reason the old panels used two:
        # one ribbon with a mapped fill joins its own points straight across
        # wherever the sign flips.
        geom_ribbon(aes(x = x, ymin = wt, ymax = pmax(ko, wt), fill = "up"), na.rm = TRUE) +
        geom_ribbon(aes(x = x, ymin = pmin(ko, wt), ymax = wt, fill = "down"), na.rm = TRUE) +
        geom_line(aes(x = x, y = wt, colour = "WT"), linewidth = .45, na.rm = TRUE) +
        geom_line(aes(x = x, y = ko, colour = "HP1gKO"), linewidth = .45, na.rm = TRUE) +
        scale_fill_manual(values = c(up = "#b2182b", down = "#2166ac"),
                          breaks = c("up", "down"),
                          labels = c("more in cKO", "less in cKO"), name = NULL) +
        scale_colour_manual(values = geno_colors, breaks = c("WT", "HP1gKO"),
                            labels = c("WT mean", "cKO mean"), name = NULL) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(n.breaks = 3, expand = c(0, 0)) +
        coord_cartesian(xlim = c(from, to) / 1e6, ylim = rng, expand = FALSE) +
        labs(y = sprintf("%s\ncoverage", mk),
             x = if (last) sprintf("%s (Mb)", chrom) else NULL) +
        theme_m + theme(legend.position = if (show_key) "right" else "none",
                        legend.key.size = unit(9, "pt"), legend.text = element_text(size = 7),
                        plot.margin = margin(1, 4, if (last) 2 else 1, 4),
                        axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5,
                                                    size = 7.5, lineheight = 0.95),
                        axis.text.y = element_text(size = 6.5),
                        axis.title.x = if (last) element_text() else element_blank(),
                        axis.text.x  = if (last) element_text() else element_blank())
    }

    # One panel per mark; the legend is collected onto the last one only, so the
    # key appears once for the whole stack rather than beside every mark.
    panels <- lapply(seq_along(marks), function(i)
                     mk_panel(marks[i], FALSE, show_key = (i == length(marks))))
    panels <- Filter(Negate(is.null), panels)

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
                      axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 8),
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
    # The coverage key strip was removed: it took a row on every figure to label
    # a fade that the caption already describes, and it read as another data
    # track rather than a legend. The fade itself is unchanged.

    # --- CTCF strip -------------------------------------------------------------
    # External to both assays and to all five marks, so unlike insulation it
    # cannot move with the KO. BULK hippocampus, NOT sorted CA1 - see
    # prepare_CTCF_annotation.R. Drawn as segments rather than filled rects
    # because at a whole-chromosome width ~90k peaks would otherwise merge into a
    # solid bar; one tick per peak keeps the density readable.
    ctcf <- tryCatch(load_ctcf(), error = function(e) NULL)
    p_ctcf <- if (is.null(ctcf)) NULL else {
      cg <- as.data.table(as.data.frame(IRanges::subsetByOverlaps(
              ctcf[ctcf$tissue == "hippocampus"], GRanges(chrom, IRanges(from, to)))))
      if (!nrow(cg)) NULL else {
        cg <- cg[, .(x = pmax(pmin((start + end) / 2, to), from) / 1e6)]
        ggplot(cg) +
          geom_segment(aes(x = x, xend = x, y = 0, yend = 1),
                       linewidth = .25, colour = "#3f6b52") +
          scale_x_continuous(expand = c(0, 0), limits = c(from, to) / 1e6) +
          scale_y_continuous(expand = c(0, 0)) +
          labs(y = "CTCF") +
          theme_m + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                          axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 7),
                          panel.grid = element_blank(), plot.margin = margin(1, 4, 1, 4))
      }
    }
    ctcf_h  <- if (is.null(p_ctcf)) numeric(0) else gam12_h_ctcf
    lad_h   <- if (is.null(p_lad)) numeric(0) else gam12_h_lad
    track_h <- c(hot_h, ctcf_h, lad_h, rep(gam12_h_mark, length(panels)), gam12_h_gene)
    heights <- grid::unit(c(rep(map_h_in, n_tri), track_h), "in")
    # guides = "collect" is essential here, not cosmetic. Without it every panel
    # keeps its own legend, and patchwork centres each one on ITS OWN ROW - so the
    # LAD legend (two keys) was being centred on a 0.22 in row and spilled over
    # the legends above and below it. Collected and top-justified, they stack.
    wrap_plots(c(tri, hot_panels,
                 if (is.null(p_ctcf)) NULL else list(p_ctcf),
                 if (is.null(p_lad))  NULL else list(p_lad),
                 panels, list(p_gene)), ncol = 1, heights = heights) +
      plot_layout(guides = "collect") +
      plot_annotation(
        title = sprintf("%s:%.2f-%.2f Mb", chrom, from / 1e6, to / 1e6),
        theme = theme(plot.caption = element_text(hjust = 0, size = 7),
                      plot.title = element_text(face = "bold", size = 13),
                      legend.spacing.y = grid::unit(12, "pt"),
                      legend.box.spacing = grid::unit(10, "pt")),
        caption = paste0(
        if (npmi_cov < gam12_cov_warn)
          sprintf("DATA COVERAGE WARNING: only %.0f%% of NPMI cells in this region are measured. White areas of the contact maps are MISSING DATA, not absent contact.\n", 100 * npmi_cov)
        else sprintf("NPMI cells measured in this region: %.0f%%. ", 100 * npmi_cov),
        "GAM: CA1 pyramidal cells, HP1g WT/cKO, 20 kb, |delta-z|>2, max NPMI>0.1. ",
        "ChIP: Emx1-Cre cortex+hippocampus (pyramidal neurons), per-bin log2FC (KO/WT).\n",
        "Contact map block-averaged and 3x3 smoothed for display. ",
        "MARK TRACKS, how they are built: per-bin log2FC(cKO/WT) at 20 kb is aggregated to 100 kb display bins, each 20 kb bin ",
        "weighted by its baseMean, then resampled to a smooth outline. Aggregation is necessary because per-bin log2FC at 20 kb is ",
        "noise-dominated (sd 0.14-0.16 for the broad marks, up to 0.57 for H3K4me3). All statistics elsewhere use unsmoothed values.\n",
        "MARK TRACKS ARE COVERAGE, NOT A RATIO. Each panel overlays the two genotypes: the solid lines are the mean per-bin coverage ",
        "across that genotype's libraries (6 WT, 4 cKO - cKO rep2 excluded at QC), and the pale band behind each is the min-max range ",
        "across those libraries. The shaded gap between the two means is the effect, coloured red where the cKO is higher and blue where ",
        "it is lower - the same convention as the contact maps above. Read it as: height is how much signal is there, the gap is the size ",
        "of the change, and the gap matters only if it is wider than the replicate bands.\n",
        "This replaced a log2 fold-change panel, which drew a ratio on a per-mark symmetric axis so every mark filled its panel equally ",
        "whatever the size of its effect - H3K9me2 (genome-wide median -0.18) read as dramatically as H4K20me3 (-0.50). The genome-wide ",
        "medians are in the ChIP section of the report; this figure is for seeing where signal actually is.\n",
        "THE Y SCALE IS PER MARK and starts at zero. H3K4me3 is punctate and high, H3K9me2 broad and low, so one shared scale would ",
        "flatten the punctate marks to a line. Compare marks by the gap relative to their own bands, which is scale-free, not by height.\n",
        "Diff. contacts = count of significant differential contacts per 20 kb bin within 5 Mb (the dashed line on the map); ",
        "net = gained - lost, same sign convention as the difference map. Those are COUNTS, every bin has one. ",
        "Hotspot is the separate CALL on that count: bins at or above the collaborators' 50-75 per-window cutoff. ",
        "LADs: constitutive cLAD/ciLAD, Peric-Hupkes 2010 (GSE17051), lifted mm9->mm39. ",
        "CTCF: ENCODE ENCSR877MSN conservative IDR peaks, BULK hippocampus (not sorted CA1), lifted mm10->mm39; one tick per peak.")) &
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
    list("chr18", 35.5e6, 39.5e6, "chr18_Pcdh_context"),
    # The other derepressed clustered families, each the largest cluster of its
    # family by gene count, padded ~1 Mb so the flanks are visible. Located from
    # the ChIP-side family regexes (^Zfp, ^Vmn, ^(Olfr|Or[0-9])) rather than
    # hand-picked coordinates, so they move if the annotation does.
    list("chr13", 65.0e6, 69.2e6,  "chr13_KRABZFP"),    # 28 Zfp genes, 66.01-68.23
    list("chr17", 20.0e6, 25.1e6,  "chr17_KRABZFP"),    # 23 Zfp genes, 20.98-24.07
    list("chr7",  18.6e6, 25.0e6,  "chr7_Vmn"),         # 87 Vmn genes, 19.58-24.01
    list("chr2",  83.9e6, 91.5e6,  "chr2_Olfr")         # 224 Olfr genes, 84.91-90.48
  )
  FIGW <- 11; SIDE <- 2.6
  # Spot-render one region without re-rendering all of them:
  #   GAM7_REGIONS=chr18_Pcdh Rscript GAM_7_Figures.R
  # Each region loads NPMI matrices, so rendering the whole list to check one
  # layout change is an hour of I/O. Comma-separated, matched against the region
  # name; unset means all of them.
  only <- Sys.getenv("GAM7_REGIONS", "")
  if (nzchar(only)) {
    want <- trimws(strsplit(only, ",")[[1]])
    regions <- Filter(function(r) r[[4]] %in% want, regions)
    message("GAM7_REGIONS set: rendering ", length(regions), " of the region list")
  }
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
})


# ==============================================================================
# SECTION B0 - all 19 autosomes, by chromatin class and separation
# ==============================================================================
# Reads gam10_perchrom_sums.tsv (sums and counts per chromosome x genotype x
# class x band, written by GAM_5) rather than the NPMI matrices, so this covers
# every autosome at no I/O cost. The contact-budget figure in section B below
# loads matrices and is therefore limited to the five FOCUS chromosomes; this one
# is the all-19 context any single-chromosome claim has to be read against.
#
# TWO PANELS, and the second is the one to read. NPMI is a relative measure, so
# anything that shifts a chromosome's detection frequency between genotypes moves
# every class on that chromosome in the same direction - a whole-chromosome
# offset. chr7 (up in all classes at all bands) and chr11 (down in all classes
# past 1 Mb) are the clearest cases. Such an offset CANCELS in the difference
# between two classes on the same chromosome, so the het-het minus eu-eu contrast
# in panel B is offset-invariant and the per-class panel A is not.
local({
  fig_sub <- "chromosomes"
  msg <- function(...) cat(sprintf(...), "\n", sep = "")
  f <- file.path(tables_dir, "gam10_perchrom_sums.tsv")
  if (!file.exists(f)) { msg("gam10_perchrom_sums.tsv missing - section B0 skipped"); return(invisible()) }
  P <- fread(f)
  BANDS_O <- c("<0.2Mb", "0.2-1Mb", "1-5Mb", "5-10Mb", ">10Mb")
  CLS_LAB <- c(intra_dom = "het-het, same domain", inter_het = "het-het, different domains",
               het_eu = "het-euchromatin", eu_eu = "euchromatin-euchromatin")
  P[, mean_npmi := s / n]
  W <- dcast(P, chr + cls + band ~ genotype, value.var = c("mean_npmi", "n"))
  W[, pct := 100 * (mean_npmi_HP1gKO / mean_npmi_WT - 1)]
  W[, band := factor(band, BANDS_O)]
  W[, chr  := factor(chr, paste0("chr", 1:19))]
  W[, cls_lab := factor(CLS_LAB[cls], CLS_LAB)]
  fwrite(W[order(chr, cls, band)], file.path(tables_dir, "gam15_perchrom_class_band.tsv"), sep = "\t")

  # --- panel A: every class on every autosome --------------------------------
  pA <- ggplot(W[is.finite(pct)], aes(band, pct, colour = cls_lab, group = cls_lab)) +
    geom_hline(yintercept = 0, linewidth = .3, colour = "grey50") +
    geom_line(linewidth = .5) + geom_point(size = 1.1) +
    facet_wrap(~ chr, ncol = 5) +
    scale_colour_manual(values = setNames(disc_pal(4), CLS_LAB), name = NULL) +
    labs(title = "Change in contact frequency by chromatin class and separation, all 19 autosomes",
         subtitle = paste0("Percent change in mean NPMI, HP1gKO vs WT. Read this panel for the ",
                           "SHAPE across separations within a chromosome;\nits overall level is ",
                           "confounded by a per-chromosome offset (chr7 up throughout, chr11 down ",
                           "throughout) - see panel B."),
         x = "genomic separation", y = "% change in mean NPMI") +
    theme_m + theme(legend.position = "top",
                    axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  save_fig(pA, "perchrom_class_band_all19", subdir = fig_sub, width = 12, height = 10)

  # --- panel B: the offset-invariant contrast --------------------------------
  CC <- dcast(W[cls %in% c("inter_het", "eu_eu")], chr + band ~ cls, value.var = "pct")
  CC[, contrast := inter_het - eu_eu]
  CC <- CC[is.finite(contrast)]
  lim <- max(abs(CC$contrast), na.rm = TRUE)
  npos <- CC[, .(n_pos = sum(contrast > 0), n = .N), by = band][order(band)]
  pB <- ggplot(CC, aes(band, chr, fill = contrast)) +
    geom_tile(colour = "white", linewidth = .4) +
    geom_text(aes(label = sprintf("%+.0f", contrast)), size = 2.7,
              colour = ifelse(abs(CC$contrast) > 0.6 * lim, "white", "grey15")) +
    scale_fill_heat0_div(limits = c(-lim, lim), name = "het-het minus\neu-eu (% pts)") +
    scale_y_discrete(limits = rev(paste0("chr", 1:19))) +
    labs(title = "Heterochromatin's gain relative to euchromatin, per autosome",
         subtitle = paste0("het-het (different domains) minus euchromatin-euchromatin, in percentage ",
                           "points. A whole-chromosome\noffset cancels in this difference, so it is the ",
                           "interpretable quantity. Positive on ",
                           paste(sprintf("%d/%d", npos$n_pos, npos$n), collapse = ", "),
                           " autosomes\nacross the five bands in order."),
         x = "genomic separation", y = NULL) +
    theme_m + theme(panel.grid = element_blank())
  save_fig(pB, "perchrom_class_contrast_all19", subdir = fig_sub, width = 8, height = 7)

  # --- panel C: ONE figure carrying all four variables ------------------------
  # chromosome (y), separation (x), chromatin class (facet) and percent change
  # (fill). 19 x 5 x 4 = 380 cells, which is the whole per-chromosome result in a
  # single view - the previous contact-budget figure showed 5 chromosomes as line
  # panels and could not be scanned down a chromosome or across a class.
  # Grey cells are combinations with no bin pairs: intra-domain pairs only exist
  # on chromosomes that have a large heterochromatin domain, and only out to the
  # separations that domain spans.
  lim2 <- max(abs(W$pct), na.rm = TRUE)
  pC <- ggplot(W[is.finite(pct)], aes(band, chr, fill = pct)) +
    geom_tile(colour = "white", linewidth = .3) +
    geom_text(aes(label = sprintf("%+.0f", pct)), size = 2.1,
              colour = ifelse(abs(W[is.finite(pct)]$pct) > 0.6 * lim2, "white", "grey15")) +
    facet_grid(~ cls_lab) +
    scale_fill_heat0_div(limits = c(-lim2, lim2), name = "% change in\nmean NPMI") +
    scale_y_discrete(limits = rev(paste0("chr", 1:19))) +
    labs(title = "Contact budget: every autosome, chromatin class and separation in one view",
         subtitle = paste0("Percent change in mean NPMI, HP1gKO vs WT. Grey = no bin pairs in that ",
                           "combination (intra-domain pairs exist only where a\nlarge heterochromatin ",
                           "domain does). Read ACROSS a row for one chromosome's distance profile, DOWN ",
                           "a column for how\nconsistent a class and separation are across the genome. ",
                           "A whole row shifting one way is a per-chromosome offset, not biology."),
         x = "genomic separation", y = NULL) +
    theme_m + theme(panel.grid = element_blank(), strip.text = element_text(size = 8, face = "bold"),
                    axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  save_fig(pC, "contact_budget_all19_unified", subdir = fig_sub, width = 13.5, height = 7.5)

  msg("\n=== het-het minus eu-eu, per autosome (percentage points) ===")
  print(dcast(CC, chr ~ band, value.var = "contrast"))
  msg("positive on: %s of %d autosomes, by band (%s)",
      paste(npos$n_pos, collapse = "/"), npos$n[1], paste(npos$band, collapse = ", "))
})

# ==============================================================================
# SECTION B - per-chromosome characterisation              (was GAM_15)
# ==============================================================================

local({
  set.seed(1)
  FOCUS   <- c("chr3", "chr7", "chr11", "chr16", "chr18")
  BANDS   <- list(c(2e4, 1e6), c(1e6, 5e6), c(5e6, 1e7), c(1e7, Inf))
  BANDNM  <- c("0.02-1Mb", "1-5Mb", "5-10Mb", ">10Mb")
  WIN     <- 2e6                      # smoothing window for the along-chromosome view
  fig_sub <- "chromosomes"
  msg <- function(...) cat(sprintf(...), "\n", sep = "")

  X <- readRDS(file.path(rds_dir, "gam_bins_annotated.rds"))
  A <- as.data.table(as.data.frame(X$annot)); setorder(A, chr, start)
  for (mk in colnames(X$l2fc)) A[[mk]] <- X$l2fc[match(A$bin_id, rownames(X$l2fc)), mk]
  A <- A[chr %in% gam_autosomes]
  A[, net := as.numeric(n_contacts_gained) - as.numeric(n_contacts_lost)]

  FE <- tryCatch(as.data.table(readRDS(file.path(rds_dir, "npmi_bin_features.rds"))$features),
                 error = function(e) NULL)
  if (!is.null(FE)) A <- merge(A, FE[, .(bin_id, insulation_WT, insulation_HP1gKO,
                                         insulation_delta, npmi_short_WT, npmi_short_HP1gKO)],
                               by = "bin_id", all.x = TRUE)
  setorder(A, chr, start)

  blocks <- tryCatch(fread(file.path(tables_dir, "gam6_large_blocks.tsv")),
                     error = function(e) NULL)
  lads <- tryCatch(load_lads(), error = function(e) NULL)

  # --- heterochromatin membership (same definition as GAM_13 bin-level) ---------
  A[, het_frac := {
    dm <- fread(het_domain_bed, header = FALSE, select = 1:3,
                col.names = c("chr", "start", "end"))
    dm[, chr := ifelse(grepl("^chr", chr), chr, paste0("chr", chr))]
    dm <- dm[chr %in% gam_autosomes]
    red <- reduce(GRanges(dm$chr, IRanges(dm$start + 1, dm$end)), ignore.strand = TRUE)
    b <- GRanges(chr, IRanges(start + 1, end))
    h <- findOverlaps(b, red, ignore.strand = TRUE)
    w <- width(pintersect(b[queryHits(h)], red[subjectHits(h)], ignore.strand = TRUE))
    v <- numeric(length(b)); agg <- tapply(w, queryHits(h), sum)
    v[as.integer(names(agg))] <- as.numeric(agg)
    pmin(v / width(b), 1)
  }]
  A[, is_het := het_frac >= 0.5]
  stopifnot(sum(A$is_het) > 1000)

  if (!is.null(lads)) {
    for (cl in c("cLAD", "ciLAD")) {
      L <- reduce(lads[lads$class == cl], ignore.strand = TRUE)
      b <- GRanges(A$chr, IRanges(A$start + 1, A$end))
      A[[paste0("frac_", cl)]] <- {
        h <- findOverlaps(b, L, ignore.strand = TRUE)
        w <- width(pintersect(b[queryHits(h)], L[subjectHits(h)], ignore.strand = TRUE))
        v <- numeric(length(b)); agg <- tapply(w, queryHits(h), sum)
        v[as.integer(names(agg))] <- as.numeric(agg); pmin(v / width(b), 1)
      }
    }
  }

  # ==============================================================================
  # 1. The genome-wide ranking these five must be read against
  # ==============================================================================
  rk <- A[, .(rho_H4K20me3 = as.numeric(suppressWarnings(
                cor(H4K20me3, net, method = "spearman", use = "complete.obs"))),
              rho_H3K9me3 = as.numeric(suppressWarnings(
                cor(H3K9me3, net, method = "spearman", use = "complete.obs"))),
              med_H4K20me3 = as.numeric(median(H4K20me3, na.rm = TRUE)),
              frac_het = mean(is_het), n_bins = .N), by = chr]
  setorder(rk, -rho_H4K20me3); rk[, rank := seq_len(.N)]
  fwrite(rk, file.path(tables_dir, "gam15_chromosome_ranking.tsv"), sep = "\t")
  msg("=== per-chromosome rho(H4K20me3 log2FC, net contact change), all 19 ===")
  print(rk[, .(chr, rho_H4K20me3 = round(rho_H4K20me3, 3),
               rho_H3K9me3 = round(rho_H3K9me3, 3), rank, focus = chr %in% FOCUS)])
  msg("genome-wide rho = %.4f; spread across chromosomes %.3f to %.3f - small and near-symmetric,",
      as.numeric(cor(A$H4K20me3, A$net, method = "spearman", use = "complete.obs")),
      min(rk$rho_H4K20me3), max(rk$rho_H4K20me3))
  msg("so a chromosome's RANK is weak evidence on its own.")

  # A tempting over-read, tested and rejected. chr7 is by far the most
  # heterochromatic chromosome (20.2% of bins vs 7.2% genome-wide) AND gains
  # contacts at every separation, while chr11 is the least (3.1%) and loses most -
  # which invites "heterochromatin-rich chromosomes gain". Across all 19 autosomes
  # that relationship is absent: rho = -0.014, p = 0.96. The two are simply the
  # extremes of both rankings at once. Recorded so it is not re-derived.
  chr_net <- A[, .(mean_net = mean(net, na.rm = TRUE)), by = chr]
  rk <- merge(rk, chr_net, by = "chr")
  ct <- suppressWarnings(cor.test(rk$frac_het, rk$mean_net, method = "spearman"))
  msg("\nheterochromatin content vs net contact change, across the 19 autosomes: rho = %.3f, p = %.3f",
      ct$estimate, ct$p.value)
  msg("  (chr7 is 20.2%% het and gains, chr11 is 3.1%% and loses - but across all 19 there is no relationship)")
  fwrite(rk, file.path(tables_dir, "gam15_chromosome_ranking.tsv"), sep = "\t")

  p <- ggplot(rk, aes(frac_het, mean_net, label = chr)) +
    geom_hline(yintercept = 0, linewidth = .3, colour = "grey50") +
    geom_point(aes(colour = chr %in% FOCUS), size = 2) +
    ggrepel::geom_text_repel(size = 2.8, colour = "grey30", max.overlaps = 20) +
    scale_colour_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "grey65"), guide = "none") +
    labs(title = "Heterochromatin content does not predict a chromosome's contact change",
         subtitle = sprintf(paste0("Spearman rho = %.3f, p = %.2f over the 19 autosomes. chr7 (most heterochromatic, gains) ",
                                   "and chr11 (least, loses most)\nare the extremes of both rankings at once, which is not ",
                                   "the same as a relationship."), ct$estimate, ct$p.value),
         x = "fraction of bins heterochromatic", y = "mean net contact change per bin") +
    theme_m
  save_fig(p, "het_content_vs_net", subdir = fig_sub, width = 7.5, height = 5.5)

  p <- ggplot(rk, aes(reorder(chr, rho_H4K20me3), rho_H4K20me3, fill = chr %in% FOCUS)) +
    geom_hline(yintercept = 0, linewidth = .3, colour = "grey50") +
    geom_col(width = .7) + coord_flip() +
    scale_fill_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "grey75"), guide = "none") +
    labs(title = "Mark-contact coupling by chromosome",
         subtitle = paste0("Spearman rho between per-bin H4K20me3 log2FC and net contact change. ",
                           "Red = the five examined here.\nThe whole spread is -0.15 to +0.13 around a ",
                           "genome-wide 0.012, so the ordering is weak evidence."),
         x = NULL, y = "Spearman rho") +
    theme_m
  save_fig(p, "chromosome_ranking", subdir = fig_sub, width = 7, height = 6)

  # ==============================================================================
  # 2. Per-chromosome: contact budget by distance and chromatin class
  # ==============================================================================
  band_of <- function(d) {
    out <- rep(NA_integer_, length(d))
    for (i in seq_along(BANDS))
      out[is.na(out) & d >= BANDS[[i]][1] & d < BANDS[[i]][2]] <- i
    out
  }

  budget <- list()
  for (ch in FOCUS) {
    sub <- A[chr == ch]
    M <- tryCatch(list(WT = load_npmi("WT", ch), KO = load_npmi("HP1gKO", ch)),
                  error = function(e) NULL)
    if (is.null(M)) { msg("%s: no NPMI, skipped", ch); next }
    ids <- intersect(sub$bin_id, rownames(M$WT)); ids <- intersect(ids, rownames(M$KO))
    sub <- sub[match(ids, bin_id)]
    W <- M$WT[ids, ids]; K <- M$KO[ids, ids]
    iu <- which(upper.tri(W), arr.ind = TRUE)
    d  <- abs(sub$start[iu[, 2]] - sub$start[iu[, 1]])
    bd <- band_of(d)
    hi <- sub$is_het[iu[, 1]]; hj <- sub$is_het[iu[, 2]]
    cls <- fifelse(hi & hj, "het-het", fifelse(!hi & !hj, "eu-eu", "het-eu"))
    wv <- W[iu]; kv <- K[iu]
    ok <- is.finite(wv) & is.finite(kv) & !is.na(bd)
    dt <- data.table(band = BANDNM[bd[ok]], cls = cls[ok], wv = wv[ok], kv = kv[ok])
    budget[[ch]] <- dt[, .(WT = mean(wv), KO = mean(kv), n = .N,
                           pct = (mean(kv) - mean(wv)) / mean(wv) * 100),
                       by = .(band, cls)][, chr := ch][]
    rm(M, W, K, iu, wv, kv, dt); gc(FALSE)
    msg("%s: contact budget done (%d bins)", ch, length(ids))
  }
  BUD <- rbindlist(budget)
  fwrite(BUD, file.path(tables_dir, "gam15_contact_budget.tsv"), sep = "\t")
  msg("\n=== contact budget by chromosome, class and distance ===")
  print(dcast(BUD, chr + cls ~ band, value.var = "pct"))

  BUD[, band := factor(band, levels = BANDNM)]
  p <- ggplot(BUD, aes(band, pct, colour = cls, group = cls)) +
    geom_hline(yintercept = 0, linewidth = .3, colour = "grey50") +
    geom_line() + geom_point(size = 1.6) +
    facet_wrap(~ chr, nrow = 1) +
    scale_colour_manual(values = c(`het-het` = "#b2182b", `eu-eu` = "#2166ac",
                                   `het-eu` = "grey60"), name = NULL) +
    labs(title = "Contact budget by chromatin class and separation",
         subtitle = paste0("Percent change in mean NPMI, HP1gKO vs WT. Genome-wide the pattern is het-het up ",
                           "and eu-eu down;\na chromosome that inverts that is the interesting case."),
         x = "genomic separation", y = "% change in mean NPMI") +
    theme_m
  save_fig(p, "contact_budget_by_chromosome", subdir = fig_sub, width = 13, height = 4.4)

  # ==============================================================================
  # 3. Along-chromosome overview + a findings file per chromosome
  # ==============================================================================
  roll <- function(y, k) {
    g <- is.finite(y)
    num <- as.numeric(stats::filter(ifelse(g, y, 0), rep(1, k), sides = 2))
    den <- as.numeric(stats::filter(g * 1, rep(1, k), sides = 2))
    ifelse(is.finite(den) & den > 0, num / den, NA_real_)
  }
  kw <- max(3L, as.integer(WIN / gam_binsize))

  dir.create(file.path(results_dir, "chromosomes"), showWarnings = FALSE, recursive = TRUE)

  for (ch in FOCUS) {
    sub <- copy(A[chr == ch]); setorder(sub, start)
    sub[, `:=`(net_s = roll(net, kw), h4_s = roll(H4K20me3, kw),
               ins_s = if ("insulation_delta" %in% names(sub)) roll(insulation_delta, kw) else NA_real_)]
    sub[, mid := (start + end) / 2e6]

    tracks <- list(
      ggplot(sub, aes(mid, net_s)) +
        geom_hline(yintercept = 0, linewidth = .3, colour = "grey50") +
        geom_area(aes(fill = net_s > 0), na.rm = TRUE) +
        scale_fill_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "#2166ac"), guide = "none") +
        labs(y = "net diff.\ncontacts", x = NULL) + theme_m +
        theme(axis.text.x = element_blank(), axis.title.y = element_text(size = 7)),
      ggplot(sub, aes(mid, h4_s)) +
        geom_hline(yintercept = median(A$H4K20me3, na.rm = TRUE), linewidth = .3, colour = "grey40") +
        geom_line(colour = "#53426b", na.rm = TRUE) +
        labs(y = "H4K20me3\nlog2FC", x = NULL) + theme_m +
        theme(axis.text.x = element_blank(), axis.title.y = element_text(size = 7)),
      ggplot(sub, aes(mid, het_frac)) +
        geom_area(fill = "grey55", na.rm = TRUE) +
        labs(y = "het\nfraction", x = sprintf("%s (Mb)", ch)) + theme_m +
        theme(axis.title.y = element_text(size = 7)))
    if ("frac_cLAD" %in% names(sub)) {
      tracks <- append(tracks, list(
        ggplot(sub, aes(mid, frac_cLAD)) + geom_area(fill = "#53426b", na.rm = TRUE) +
          labs(y = "cLAD\nfraction", x = NULL) + theme_m +
          theme(axis.text.x = element_blank(), axis.title.y = element_text(size = 7))), after = 2)
    }
    pp <- wrap_plots(tracks, ncol = 1, heights = rep(1, length(tracks))) +
      plot_annotation(title = sprintf("%s overview", ch),
        subtitle = sprintf("All tracks smoothed over %.0f Mb. Net diff. contacts = gained - lost per 20 kb bin; red = more contact in the KO.", WIN / 1e6),
        theme = theme(plot.title = element_text(face = "bold", size = 13)))
    save_fig(pp, sprintf("%s_overview", ch), subdir = fig_sub, width = 11, height = 2 + 1.5 * length(tracks))

    # --- the findings file ---
    bl <- if (!is.null(blocks)) blocks[chr == ch] else NULL
    bu <- BUD[chr == ch]
    f <- file.path(results_dir, "chromosomes", sprintf("%s.md", ch))
    con <- file(f, "w")
    writeLines(c(
      sprintf("# %s", ch), "",
      sprintf("Generated %s by `GAM_7_Figures.R` (section B). Numbers here are regenerable; do not hand-edit.", Sys.Date()), "",
      "## Mark-contact coupling, in context", "",
      sprintf("- rho(H4K20me3 log2FC, net contact change) = **%.3f**, rank %d of %d.",
              rk[chr == ch, rho_H4K20me3], rk[chr == ch, rank], nrow(rk)),
      sprintf("- The spread across all 19 autosomes is %.3f to %.3f around a genome-wide %.4f, so this rank is weak evidence on its own.",
              min(rk$rho_H4K20me3), max(rk$rho_H4K20me3),
              as.numeric(cor(A$H4K20me3, A$net, method = "spearman", use = "complete.obs"))),
      sprintf("- Median H4K20me3 log2FC here %.3f, against %.3f genome-wide.",
              rk[chr == ch, med_H4K20me3], median(A$H4K20me3, na.rm = TRUE)),
      sprintf("- %.1f%% of bins are heterochromatic (>=50%% coverage by the merged H3K9me3/H4K20me3 set), against %.1f%% genome-wide.",
              100 * rk[chr == ch, frac_het], 100 * mean(A$is_het)), "",
      "## Contact budget", "",
      "Percent change in mean NPMI, HP1gKO vs WT, by chromatin class and separation.", "",
      "| class | 0.02-1Mb | 1-5Mb | 5-10Mb | >10Mb |", "|---|---|---|---|---|"),
      con)
    for (cl in c("het-het", "het-eu", "eu-eu")) {
      r <- sapply(BANDNM, function(b) { v <- bu[cls == cl & band == b, pct]
                                        if (length(v)) sprintf("%+.1f%%", v) else "-" })
      writeLines(sprintf("| %s | %s |", cl, paste(r, collapse = " | ")), con)
    }
    # A single averaged verdict is misleading here: chr11 averages to "follows"
    # while losing 33.8% of het-het contacts at 5-10 Mb, which is the opposite of
    # the genome-wide gain. State the contrast per band and let it speak.
    dif <- sapply(BANDNM, function(b) {
      h <- bu[cls == "het-het" & band == b, pct]; e <- bu[cls == "eu-eu" & band == b, pct]
      if (length(h) && length(e)) h - e else NA_real_
    })
    writeLines(c("",
      "het-het minus eu-eu, per band (genome-wide this is positive: heterochromatin gains while euchromatin loses):", "",
      paste0("| ", paste(BANDNM, collapse = " | "), " |"),
      paste0("|", paste(rep("---", length(BANDNM)), collapse = "|"), "|"),
      paste0("| ", paste(ifelse(is.finite(dif), sprintf("%+.1f", dif), "-"), collapse = " | "), " |"),
      "",
      sprintf("Bands where this chromosome runs against the genome-wide direction: %s.",
              { neg <- BANDNM[is.finite(dif) & dif < 0]
                if (length(neg)) paste(neg, collapse = ", ") else "none" }), ""), con)
    if (!is.null(bl) && nrow(bl)) {
      writeLines(c("## Large one-directional blocks (GAM_6)", "",
        sprintf("%d block(s) on this chromosome, %.1f Mb in total.", nrow(bl),
                sum(bl$end - bl$start) / 1e6), ""), con)
      writeLines("| start (Mb) | end (Mb) | width (Mb) | direction | mean net |", con)
      writeLines("|---|---|---|---|---|", con)
      setorder(bl, start)
      for (i in seq_len(nrow(bl)))
        writeLines(sprintf("| %.2f | %.2f | %.2f | %s | %+.2f |", bl$start[i] / 1e6, bl$end[i] / 1e6,
                           (bl$end[i] - bl$start[i]) / 1e6, as.character(bl$dirn[i]),
                           bl$mean_net[i]), con)
      writeLines("", con)
    }
    writeLines(c("## Caveats carried from the pipeline", "",
      "- GAM is 2 biological replicates per genotype, already merged, so there is no within-genotype variance; every claim here is an effect size against a matched or bootstrap null, never a per-contact p-value.",
      "- Per-bin ChIP log2FC at 20 kb is noise-dominated (per-bin sd 0.14-0.16 for the broad marks). Chromosome-level medians are sound; do not read local structure off a single bin.",
      "- The A/B compartment axis is weakly determined on most chromosomes and is not used here.",
      ""), con)
    close(con)
    msg("wrote %s", f)
  }

  msg("\nGAM_15 complete.")
})
