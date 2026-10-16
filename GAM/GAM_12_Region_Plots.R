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
# NOTE ON TISSUE: the GAM is CA1; the ChIP is cortex (Emx1-Cre). Both are
# Emx1-lineage but they are different tissues, and this figure overlays them.
# ================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(patchwork) })

gam12_maxdist   <- 5e6     # matches the hotspot integration window
gam12_hot_cut   <- 50      # slides: 50-75 differential contacts per window

g1 <- file.path(rds_dir, "gam_bins_annotated.rds")
if (!file.exists(g1)) stop("GAM_1 output not found:\n  ", g1,
  "\nRun it first:  cd GAM && Rscript GAM_1_Bins_and_Signal.R", call. = FALSE)
X <- readRDS(g1)
A <- as.data.table(X$annot); setorder(A, chr, start)
for (mk in colnames(X$l2fc)) A[[mk]] <- X$l2fc[match(A$bin_id, rownames(X$l2fc)), mk]

# ---------------------------------------------------------------
# gam_region_plot(chrom, from, to, marks)
#   chrom  "chr11"; from/to in bp (NULL = whole chromosome)
#   marks  which ChIP tracks to stack, in order
# Returns a patchwork object; save with save_fig().
# ---------------------------------------------------------------
gam_region_plot <- function(chrom, from = NULL, to = NULL,
                            marks = c("H4K20me3", "H3K9me3", "H3K4me3"),
                            maxdist = gam12_maxdist, hot_cut = gam12_hot_cut,
                            smooth_map = TRUE) {

  reg <- A[chr == chrom]
  if (is.null(from)) from <- min(reg$start)
  if (is.null(to))   to   <- max(reg$end)
  reg <- reg[start >= from & end <= to]
  if (nrow(reg) < 10) stop("region too small: ", chrom, ":", from, "-", to, call. = FALSE)

  # --- contact difference over the region ------------------------------------
  mats <- lapply(c(WT = "WT", KO = "HP1gKO"), function(g) {
    m <- load_npmi(g, chrom)
    ids <- intersect(reg$bin_id, rownames(m))
    m[ids, ids, drop = FALSE]
  })
  ids <- intersect(rownames(mats$WT), rownames(mats$KO))
  D <- mats$KO[ids, ids] - mats$WT[ids, ids]
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
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "KO - WT\nNPMI",
                         guide = guide_colourbar(barheight = 4)) +
    scale_y_continuous(expand = c(0, 0), name = "separation (Mb)") +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(from, to) / 1e6) +
    labs(title = sprintf("%s:%.1f-%.1f Mb", chrom, from / 1e6, to / 1e6),
         subtitle = sprintf("Contact difference to %.0f Mb separation. Red = more contact in cKO, blue = more in WT.",
                            maxdist / 1e6)) +
    theme_m + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                    plot.margin = margin(4, 4, 1, 4))

  # --- hotspots: the two directions, back to back ----------------------------
  h <- reg[, .(mid = (start + end) / 2e6, lost = n_contacts_lost, gained = -n_contacts_gained)]
  hm <- melt(h, id.vars = "mid", variable.name = "dirn", value.name = "n")
  p_hot <- ggplot(hm, aes(mid, n, fill = dirn)) +
    geom_col(width = gam_binsize / 1e6) +
    geom_hline(yintercept = c(-hot_cut, hot_cut), linetype = 2, linewidth = .3, colour = "grey40") +
    geom_hline(yintercept = 0, linewidth = .3) +
    scale_fill_manual(values = c(lost = "#B2182B", gained = "#2166AC"),
                      labels = c("lost in cKO", "gained in cKO"), name = NULL) +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(from, to) / 1e6) +
    labs(y = "diff. contacts\n(0-5 Mb)") +
    theme_m + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                    legend.position = "right", plot.margin = margin(1, 4, 1, 4))

  # --- ChIP tracks ------------------------------------------------------------
  # Display smoothing for the mark tracks. At 20 kb a 25 Mb region is 1,250 bars
  # and reads as noise; a running mean over the same window used for the contact
  # display bin keeps the profile honest and legible. The underlying values are
  # untouched - this is a rendering choice, stated in the caption.
  smw <- max(1L, round(nrow(reg) / 350))
  mk_panel <- function(mk, last = FALSE) {
    d <- reg[, .(mid = (start + end) / 2e6, v = get(mk))][is.finite(v)]
    if (smw > 1L) {
      sm <- as.numeric(stats::filter(d$v, rep(1 / smw, smw), sides = 2))
      d[, v := ifelse(is.na(sm), v, sm)]
    }
    gm <- median(A[[mk]], na.rm = TRUE)
    ggplot(d, aes(mid, v)) +
      geom_col(aes(fill = v > gm), width = gam_binsize / 1e6 * max(1, smw), show.legend = FALSE) +
      geom_hline(yintercept = gm, linewidth = .3, colour = "grey30", linetype = 2) +
      scale_fill_manual(values = c(`TRUE` = "grey55", `FALSE` = "#B2182B")) +
      scale_x_continuous(expand = c(0, 0)) +
      coord_cartesian(xlim = c(from, to) / 1e6) +
      labs(y = mk, x = if (last) sprintf("%s (Mb)", chrom) else NULL) +
      theme_m + theme(plot.margin = margin(1, 4, if (last) 4 else 1, 4),
                      axis.title.x = if (last) element_text() else element_blank(),
                      axis.text.x  = if (last) element_text() else element_blank())
  }
  panels <- lapply(seq_along(marks), function(i) mk_panel(marks[i], i == length(marks)))

  wrap_plots(c(list(p_map, p_hot), panels), ncol = 1,
             heights = c(3.1, 1.25, rep(0.85, length(panels)))) +
    plot_annotation(theme = theme(plot.caption = element_text(hjust = 0, size = 7)), caption = paste0(
      "GAM: CA1 HP1g WT/cKO, 20 kb, |delta-z|>2, max NPMI>0.1. Dashed line on hotspots = ",
      hot_cut, " contacts/window. ChIP: Emx1-Cre cortex, per-bin log2FC (KO/WT); ",
      "dashed = genome median, red = below it.\n",
      "Contact map block-averaged and 3x3 smoothed for display; mark tracks running-mean smoothed. Statistics use unsmoothed values."))
}

# --- the collaborators' own example regions, plus our blocks -----------------
regions <- list(
  list("chr3",  30e6,  85e6, "chr3_30-85Mb"),
  list("chr11", 100e6, 125e6, "chr11_100-125Mb"),
  list("chr11", 1e6,   25e6,  "chr11_1-25Mb_lostblock"),
  list("chr16", NULL,  NULL,  "chr16_full"),
  list("chr18", 34e6,  40e6,  "chr18_Pcdh")
)
for (r in regions) {
  nm <- r[[4]]
  message("Plotting ", nm)
  p <- tryCatch(gam_region_plot(r[[1]], r[[2]], r[[3]]),
                error = function(e) { message("  skipped: ", conditionMessage(e)); NULL })
  if (!is.null(p)) save_fig(p, nm, subdir = "regions", width = 11, height = 9)
}
cat("\nGAM_12 complete. gam_region_plot() is reusable for any region.\n")
