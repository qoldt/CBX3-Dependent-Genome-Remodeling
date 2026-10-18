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

  # PALETTES: THE COLLABORATORS' OWN, RECOVERED FROM THEIR FIGURE.
  #
  # These were hand-matched by eye until 2026-10-18, when their chr16 full-length
  # slide was sampled directly. The screenshot is a macOS capture and therefore
  # carries a Display P3 profile, so its raw RGB has to be converted to sRGB
  # before anything is compared - read uncorrected, every colour looks washed out
  # and the match is meaningless. After that conversion both of their colourbars
  # fit a 9-class ColorBrewer ramp almost exactly (RMS over the whole bar, in
  # [0,1] RGB):
  #
  #   their delta-z bar   rev(RdBu, 9)     0.051     the old hand-made ramp: 0.138
  #   their NPMI bar      rev(RdYlBu, 9)   0.044     the old hand-made ramp: 0.056
  #
  # Both being 9-class is the giveaway that one routine drew both.
  #
  # WHY THE DIFFERENCE PANEL LOOKED "STRONGER" THAN THEIRS. The old 5-stop ramp
  # put #d6604d at 75% of the way out, where RdBu-9 is still at the pale salmon
  # #F4A582 and does not reach #d6604d until 87.5%. On the +-4 limit both use,
  # that is the difference between delta-z = +2 reading as a firm red and reading
  # as a faint blush - and most of a difference map lives between -2 and +2. Same
  # numbers, same limits, a map that looks about twice as changed.
  #
  # These are literal ColorBrewer hex rather than an `ltc` palette, which is the
  # repo rule elsewhere. The rule is deliberately set aside for these two: the
  # point of these panels is to be comparable cell-for-cell with the figures the
  # collaborators circulate, so the palette is a fact about their figure, not a
  # design choice of ours.
  #
  # NOT COPIED: the sign. Their delta-z is set1 - set2 = WT - cKO, so their red
  # means lost in the cKO. This repo is positive = MORE in the cKO throughout and
  # stays that way, so our maps are colour-inverted relative to theirs by design.
  # The panel subtitle says "cKO - WT" on every figure for exactly this reason.
  #
  # RAW NPMI is a sequential ramp over the positive range, NOT a diverging scale
  # centred on zero. Centring at zero (an earlier attempt here) pushes nearly
  # every value into the warm half and flattens the structure, because almost all
  # NPMI in a cis matrix is positive.
  npmi_cols <- c("#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF",
                 "#FEE090", "#FDAE61", "#F46D43", "#D73027")   # rev(RdYlBu, 9)
  diff_cols <- c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7",
                 "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")   # rev(RdBu, 9)
  # Fixed delta-z limit, again matching their figures (scale bar -4 to +4). A
  # data-driven limit (we used the 99th percentile of |delta-z|) rescales every
  # panel to its own extremes, saturates the whole map and makes regions
  # incomparable; a fixed limit leaves the bulk pale so the strong contacts show.
  gam12_dz_lim    <- 4
  # q97, not q99. The top of this distribution is the near-diagonal wedge, which
  # is a property of contact decay rather than of anything a reader is looking
  # for, and spending the last 2% of the ramp on it leaves the rest of the map
  # blue - on chr16 q99 is 0.413 against a median of 0.098. q97 is 0.299 there,
  # close to the 0.25 the collaborators use, and keeps the domain structure
  # legible without clipping 6% of cells flat the way a hard 0.25 does.
  # Rendered side by side at 0.25 / 0.30 / 0.42 before choosing.
  gam12_npmi_qcap <- 0.97
  gam12_dpi          <- 300    # must match save_fig()'s dpi in shared/common.R
  gam12_map_cols_max <- 1600L  # tile budget: ~1.3 M rects in the full triangle
  # RAW NPMI CAP: DATA-DRIVEN. Forcing theirs was tried and reverted.
  #
  # Their slides use a fixed 0 to 0.25 bar, and a fixed limit is normally the
  # better choice - it is why gam12_dz_lim is fixed. Copying it made the panels
  # read as solid warm with a thick saturated near-diagonal wedge, because 5.9%
  # of chr16 cells are above 0.25 and all of them clip to the same red.
  #
  # Two things were tried against it and both are recorded because neither
  # worked, and the second disproves the first:
  #   - "our display bin makes the wedge too thick". Measured the other way:
  #     going from a 280 kb display bin to 60 kb took the saturated-pixel share
  #     from 8.6% UP to 11.0%, because coarse blocks average the near-diagonal
  #     peak DOWN. Aggregation was hiding the clipping, not causing it.
  #   - "match their figure pixel for pixel". Their panel is a screenshot of a
  #     slide at ~42 kb per pixel, so its wedge is blended away by downscaling
  #     before we ever measure it. That number cannot be a target, and tuning a
  #     cap against a resampled JPEG is how the first wrong answer was reached.
  #
  # So the cap is the q99 of what each figure actually draws: 0.421 on chr16
  # full length, 0.694 inside the 2.4 Mb Pcdh window. A short region draws only
  # short separations, where NPMI is genuinely high, and no single number is
  # right for both. npmi_cap stays an argument for a caller that wants one
  # shared limit across two figures it is comparing; NULL means "q99 of this
  # figure". The PALETTE is still theirs - that part was a real match.
  gam12_npmi_cap <- NULL
  # THE RIBBON CUTOFF. The collaborators' slides propose 50-75 differential
  # contacts per window and we deliberately keep the BOTTOM of that range, which
  # is a choice and not a default: >=50 calls 17,527 bins, >=60 calls 10,374 and
  # >=75 calls 4,757, so moving to the top of the range would drop 73% of the
  # called bins. As a percentage that is 14.9/8.8/4.1% of the 117,273 autosomal
  # bins the delivered hotspot files cover, or 14.2/8.4/3.9% of the whole 123,133
  # -bin grid - say which denominator is meant, because both get quoted. The ribbon is read as "where is
  # this busy", so the permissive end is the honest one for a density display.
  # Anything that needs a conservative CALL should use the domains instead.
  gam12_hot_cut   <- 50
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
  # 25% shorter than before. Five coverage tracks stacked under three contact
  # maps made the figure tall enough that the maps shrank to fit the page, and
  # the coverage shapes are legible at this height.
  gam12_h_mark   <- 0.525
  gam12_h_gene   <- 0.80
  gam12_h_ctcf   <- 0.26   # CTCF strip, same height as the LAD ribbon
  # REPEAT PANEL. Taller than a ribbon because it carries three lines and a
  # numbered y axis, shorter than a mark panel because it needs only two breaks.
  gam12_h_repeat <- 0.42

  # INSULATION-SCORE MAP STRIPS - the Pombo lab's own display of melting.
  # One strip per genotype plus their difference, drawn directly under the
  # difference contact map: x is position, y is the insulation SQUARE SIZE, fill
  # is the insulation score (GAM_2c). It is the only view that shows at WHICH
  # SPATIAL SCALE a region changes - a single collapsed melting score cannot
  # distinguish a change confined to 240 kb squares from one present at every
  # size, and those are different events.
  #
  # The strips and the map above them are the same quantity at two different
  # readings: the map is one cell pair at a time, each strip row is the mean over
  # a whole square crossing that position. So they SHOULD agree, and a strip that
  # contradicts the map above it is a bug, not a finding.
  # PANEL height only. Each plot's plot.margin occupies its OWN gtable rows, so
  # the 5 pt the strips need above and below their square-size tick labels is
  # counted in the chrome below, not here.
  #
  # 0.45 in: 11 rows at ~13 px each at 300 dpi. The strips and the contact map
  # share one fixed height on a slide, so every 0.1 in of strip is 0.1 in off the
  # map - which is why this was 0.33 while three strips were drawn and can be
  # taller now that only one is.
  gam12_h_ismap  <- 0.45
  gam12_ismap_cols <- gam12_map_cols_max # display columns, matched to the contact map
  # FIXED symmetric limits, not per-region quantiles, for the same reason
  # gam12_dz_lim is fixed: a data-driven limit rescales every panel to its own
  # extremes and makes regions incomparable.
  #
  # The two limits are deliberately NOT the same percentile. A per-genotype
  # strip carries a large REGIONAL offset - the score is centred on the
  # chromosome mean, so a region sitting well below its chromosome runs hard
  # negative at every square size. At the genome-wide 99th percentile of |IS|
  # (0.145), 19.6% of the chr18 Pcdh region saturates and both strips go solid
  # blue over their flanks, which is exactly the information the display exists
  # to show. 0.20 is just above the genome-wide 99.9th (0.180) and saturates
  # 0.6% there, 0.0% on the other three slide regions.
  # THE DIFFERENCE STRIP IS ON THE SAME SCALE, AND THIS IS NOT NEGOTIABLE.
  # It was drawn at 0.10 first, on the reasoning that the difference has no
  # regional offset and so deserves its own tighter range. That made the panel
  # lie. Three strips stacked under one another are read BY COMPARISON, and at
  # half the limit a change of +0.034 - a sixth of the per-genotype range - came
  # out as a 34% red band under two panels that were both solid blue. Measured
  # at chr18:38.10 Mb: IS_WT -0.161, IS_cKO -0.127, difference +0.034. Nothing
  # was wrong with the number; the reader was simply given two rulers and no
  # warning. Genome-wide, 4.8% of the cells where BOTH genotypes are below -0.05
  # were drawn more than two-thirds saturated in the difference; on one shared
  # limit that is 0.0%.
  #
  # The difference is not washed out by sharing the scale: |IS(cKO) - IS(WT)|
  # reaches 52% of the bar at its genome-wide 99th percentile, and 29-51% within
  # each of the four slide regions. It is drawn smaller because it IS smaller.
  gam12_is_lim   <- 0.20
  gam12_is_dlim  <- gam12_is_lim

  # HOW MANY STRIPS. "diff" draws only cKO - WT, which is the default and what
  # the slide figures use; "all" draws WT, cKO and the difference (the Pombo
  # SCHEMATIC); "none" omits them.
  #
  # The default is one strip because three were measured and did not earn their
  # rows. Within the four regions these figures cover, the WT and cKO strips
  # correlate at rho 0.75-0.93 (chr18 0.925), and the SD of their difference is
  # only 0.48-0.70 of the SD of WT itself - so the two per-genotype strips are
  # close to one picture drawn twice, and what separates them is exactly what
  # the third strip already shows. Dropping them returns 0.66 in to the contact
  # map, which on a 16:9 slide is a band of 0.63 Mb instead of 0.40 at Pcdh and
  # 8.4 Mb instead of 5.5 on chr11.
  #
  # "all" is still the right choice when the question is specifically about one
  # genotype's own insulation landscape rather than the change in it.
  gam12_ismap_mode  <- "diff"
  gam12_ismap_modes <- c("none", "diff", "all")
  ismap_n <- function(mode) switch(match.arg(mode, gam12_ismap_modes),
                                   none = 0L, diff = 1L, all = 3L)

  # REPLICATE BANDS ON THE MARK TRACKS. Off by default: with five marks stacked,
  # two bands plus two mean lines per panel is twenty elements and the figure
  # reads as a thicket. The bands answer a real question - is the gap between the
  # genotypes wider than the spread between libraries - but that question is
  # answered once, genome-wide, in the report's ChIP section (gap/spread is 0.84
  # for H4K20me3, 0.54 H3K9me3, 0.30 H3K9me2, 0.12 H3K36me3, 0.04 H3K4me3), and
  # it does not need re-answering in every region. Set TRUE to bring them back.
  gam12_mark_bands <- FALSE

  # Marks that must be drawn on ONE shared y-axis. Each element is a group.
  gam12_mark_scale_groups <- list(c("H3K9me2", "H3K9me3"))

  # ROBUST CEILING FOR THE MARK PANELS. A single bin can be two orders of
  # magnitude above everything else in a region and then sets an axis nothing
  # else reaches. chr11:3.04-3.16 Mb is the case that forced this: it is the
  # first assembled sequence on chr11 (mm39 pads chr11 with 3 Mb of N), it sits
  # against the pericentromeric major satellite, and multi-mapping heterochromatin
  # reads pile up there - H3K9me3 reaches 67.8 against a region 99th percentile
  # of 3.7 and a genome-wide 99.9th of 11.1. The same pileup sits at exactly
  # 3.0-3.2 Mb on chr9, chr12 and chr14. Drawn to its own maximum it flattened
  # the whole 24 Mb region to the axis floor.
  # The ceiling is k x the 99th percentile of the values actually drawn, and is
  # only applied when it buys at least `clip_min` of height - so on 9 of the 11
  # regions it is inert (ceiling == maximum, nothing clipped) and it only bites
  # where there is a genuine off-scale outlier. Clipped spans are marked with a
  # dark bar along the top of the panel and the true maximum is printed in the
  # axis title, so nothing is hidden. Set gam12_mark_clip_k <- Inf to disable.
  gam12_mark_clip_q   <- 0.99
  gam12_mark_clip_k   <- 2
  gam12_mark_clip_min <- 1.25

  # MASK THE DELTA-Z MAP TO SIGNIFICANT CONTACTS? YES, by default, and this is
  # now the primary difference map rather than a variant.
  #
  # The unmasked map spends its colour range on cells near zero and the
  # structure washes out. That is exactly the contrast in the collaborators'
  # own figure - their pale "Delta z-score" panel against their "Significant"
  # panel, where the compartment-scale checkerboard is immediately legible -
  # and it reproduces here.
  #
  # The old reasoning for keeping it off was that at ~350 display columns a
  # block is ~14x14 native cells, so almost every block contains a significant
  # cell and the mask fills in like the unmasked map. That was true, and it was
  # an argument about the display bin, not about the mask: the bin is now driven
  # by the output raster (see the aggregation block below), so blocks are small
  # enough for the mask to stay sparse and the alpha channel to mean something.
  #
  # GAM7_DZ_SIG=0 renders the unmasked map instead, into *_alldz files, so both
  # can sit side by side. Note the suffix moved: the masked map is now the one
  # with the plain name.
  gam12_dz_significant_only <- !identical(Sys.getenv("GAM7_DZ_SIG", "1"), "0")
  # Caption is ~10 wrapped lines now that the two mark-track references are stated
  # in full, so this has to cover more than it used to.
  gam12_h_chrome <- 2.30
  # SLIDE CHROME, MEASURED NOT GUESSED. Every panel row in this figure is an
  # absolute height, so a figure saved shorter than its content does not shrink
  # to fit - it overflows the canvas and the outer rows are silently cut off.
  # At the old 0.55 the title and the caption were being lost on every slide
  # figure and nothing reported it.
  #
  # The real number comes from laying the figure out and adding up the rows that
  # are NOT panels:
  #   gt <- patchwork:::plot_table(p, "auto")
  #   sum(grid::convertHeight(gt$heights, "in", valueOnly = TRUE)) - sum(panels)
  # which for the slide panel set (one map, the count panel, three ribbons, two
  # marks, genes) is 1.487 in, PLUS 0.1389 for each insulation strip drawn.
  #
  # The per-strip term is not fitted: each strip asks for 5 pt of top and bottom
  # margin so its square-size tick labels are not drawn into the row above, and
  # a plot's margins occupy their own gtable rows, so 10 pt = 0.1389 in of chrome
  # per strip. Measured at 1.902 with three strips and 1.625 with one, a
  # difference of 0.277 for two - 0.1385 each.
  #
  # MEASURE IT ON patchworkGrob, NOT plot_table. plot_table returns the INNER
  # table - the panels and their own per-panel titles, axes and margins - and
  # leaves out the patchwork annotation that wraps it: the figure title (0.24),
  # the figure caption (0.17) and the outer margins (0.15). Measuring the inner
  # table gave 1.35, the figure was saved 0.55 in too short, and the title and
  # caption were pushed off the canvas with nothing reporting it. That is the
  # bug the slide figures shipped with before this. Every row here is absolute,
  # so a figure saved shorter than its content does not shrink - it overflows.
  gam12_h_chrome_slide <- function(n_ismap = 1L)
    1.487 + 0.1389 * as.numeric(n_ismap)

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
  # ...all of which applied to the log2FC tracks. THE MARK TRACKS NOW DRAW
  # COVERAGE, which is not a ratio of small counts and is not noise-dominated, so
  # the 100 kb aggregation that fixed the ratio is pure over-smoothing here: on a
  # 1.8 Mb region it left 18 points to draw a curve through, and the result read
  # as a signature rather than a track.
  #
  # The display bin is therefore chosen from the region width to hit roughly
  # gam12_mark_target_bins points, floored at the native 20 kb - fine regions get
  # native resolution, whole chromosomes still get aggregated enough to draw.
  gam12_mark_target_bins <- 400L
  gam12_mark_npoints <- 1200   # columns the mark track is resampled onto
  gam12_cov_warn     <- 0.75   # below this fraction of finite NPMI cells, warn on the figure

  # Insulation-score maps from GAM_2c, loaded ONCE and cached per chromosome.
  # The file is 19 MB and covers all 19 autosomes at 20 kb viewpoints x 11 square
  # sizes; reading it per region would dominate the render. Returns NULL when
  # GAM_2c has not been run, so the region figures still render without it.
  gam12_is_env <- new.env(parent = emptyenv())
  load_meltron_is <- function(chrom) {
    if (!is.null(gam12_is_env$miss)) return(NULL)
    if (is.null(gam12_is_env$all)) {
      f <- file.path(rds_dir, "meltron_is.rds")
      if (!file.exists(f)) {
        message("  insulation maps: ", f, " not found - run GAM_2c; strips omitted")
        gam12_is_env$miss <- TRUE
        return(NULL)
      }
      o <- readRDS(f)
      gam12_is_env$all   <- o$is
      gam12_is_env$sizes <- sort(unique(o$is$is_size))
    }
    gam12_is_env$all[chr == chrom]
  }

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
                              diff_mode = c("deltaz", "raw"),
                              npmi_cap = gam12_npmi_cap,
                              # Mask the difference map to the significant
                              # contacts. An ARGUMENT, not a global read from the
                              # environment: the driver renders both variants in
                              # one pass, and a global meant each needed its own
                              # run - which is how the masked figures silently
                              # kept an old palette for half a day.
                              dz_sig = gam12_dz_significant_only,
                              # Repeat composition, on every variant. It costs a
                              # 0.42 in row; on a slide that comes out of the
                              # contact band, on the full figure it just makes the
                              # figure taller, since those size by content.
                              repeats = TRUE,
                              # SLIDE MODE. The full figure is ~1.09:1 - essentially
                              # square - because three map rows, five mark tracks and
                              # a ten-line caption stack up to ~10.1 in at 11 in wide.
                              # A 16:9 slide at that width is 6.19 in, so this is a
                              # re-layout, not a resize: one map row instead of three,
                              # two marks instead of five, and the caption reduced to a
                              # single line (on a slide it belongs in the speaker notes).
                              slide = FALSE,
                              # INSULATION-SCORE MAP STRIPS under the difference
                              # map (WT, cKO, cKO - WT). On by default in slide
                              # mode, where the point of the figure is the change
                              # itself; off in the report figures, which already
                              # carry three contact maps and ten tracks.
                              # INSULATION-SCORE MAP STRIPS: "none", "diff"
                              # (the default on slides) or "all". See
                              # gam12_ismap_mode in the config block for why one
                              # strip rather than three.
                              ismap = if (slide) gam12_ismap_mode else "none") {
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
    # SIGNIFICANT CONTACTS ONLY, from the delivered Diff_contacts files rather
    # than from a local recomputation.
    #
    # Why. A full delta-z map spends its colour range on cells that are near zero,
    # so the structure washes out - that is the top panel of the collaborators'
    # own figure, against the bottom panel which masks to significant cells and in
    # which the compartment-scale checkerboard is immediately legible. The files
    # already carry their `deltaz` for exactly the pairs that passed their
    # criterion (|delta-z| > 2 AND max NPMI > 0.1), so this uses their number for
    # their pairs instead of reproducing both.
    #
    # SIGN: their `deltaz` is set1 - set2 = WT - cKO, so positive means LOST in
    # the cKO. This repo is positive = MORE in the cKO throughout, so it is
    # negated here and nowhere else.
    #
    # THREE STATES, not two. NA in this map would otherwise conflate "measured and
    # not significant" with "never measured", and an absence of measurement that
    # reads as a result is this repo's documented worst failure mode. So:
    #   coloured  measured and significant
    #   grey      measured, not significant
    #   white     not measured in one or both genotypes
    D <- if (identical(diff_mode, "deltaz") && dz_sig) {
      dz <- matrix(NA_real_, length(ids), length(ids), dimnames = list(ids, ids))
      for (sel in c("selWT", "selcKO")) {
        dc <- tryCatch(load_diff_contacts(chrom, sel,
                 select = c("start_x", "end_x", "start_y", "end_y", "deltaz")),
               error = function(e) NULL)
        if (is.null(dc) || !nrow(dc)) next
        dc <- dc[start_x >= from & end_y <= to]
        if (!nrow(dc)) next
        i <- match(sprintf("%s:%d-%d", chrom, dc$start_x, dc$end_x), ids)
        j <- match(sprintf("%s:%d-%d", chrom, dc$start_y, dc$end_y), ids)
        ok <- is.finite(i) & is.finite(j)
        if (!any(ok)) next
        dz[cbind(i[ok], j[ok])] <- -dc$deltaz[ok]   # to cKO - WT
        dz[cbind(j[ok], i[ok])] <- -dc$deltaz[ok]
      }
      dz
    } else if (identical(diff_mode, "deltaz")) {
      mats$KO$z[ids, ids] - mats$WT$z[ids, ids]
    } else K0 - W0
    # measured = finite raw NPMI in both genotypes; used to separate grey from white
    MEAS <- if (identical(diff_mode, "deltaz") && dz_sig)
              is.finite(W0) & is.finite(K0) else NULL
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

    # AGGREGATE before drawing - but to the OUTPUT RASTER, not to a fixed 350.
    #
    # Subsampling every k-th bin aliases: the rotated lattice beats against the
    # raster grid and the map reads as diagonal streaks rather than domains.
    # Block-averaging keeps every value and produces a smooth map. That part was
    # always right. The target was not.
    #
    # 350 columns is far coarser than our own PNG. These figures save at 300 dpi
    # over an 8.4 in panel, so the panel is ~2520 px wide; 350 columns made every
    # display block SEVEN PIXELS across and threw away a factor of 7 in
    # resolution that the file was already paying for. On chr16 that is a 280 kb
    # display bin against a 20 kb native grid.
    #
    # It is not only blockiness. The collaborators draw chr16 at roughly native
    # resolution, so their saturated near-diagonal wedge is one native bin thick;
    # ours was 14 bins thick, and at their 0 to 0.25 colour cap the whole wedge
    # clips flat. Measured on pixels inside the triangle, against their panel:
    #
    #                                   warm     cool     deep saturated red
    #     their Set1                    41.0%    26.5%        2.4%
    #     ours, 350 cols, cap 0.25      45.4%    41.5%        8.6%
    #     ours, 350 cols, cap q99       14.3%    76.1%        2.7%
    #
    # Widening the cap to q99 fixes the saturation and breaks the balance - it
    # pushes the ramp midpoint to 0.21 against a chr16 median NPMI of 0.098, so
    # three quarters of the map goes blue. Their 0.25 is the RIGHT cap for this
    # data; the display bin was the defect, and the cap was being used to paper
    # over it. Fix the resolution and the cap can stay theirs.
    #
    # gam12_map_cols_max is a tile budget, not a resolution choice: geom_tile
    # draws one rect per cell in the PDF as well as the PNG, so an uncapped
    # native triangle on a large chromosome produces a vector file nothing will
    # open. Raise it if the PDFs stay manageable.
    n <- length(ids)
    FRAC <- NULL
    target <- min(max(2, fig_width - side_allow) * gam12_dpi, gam12_map_cols_max)
    fac <- max(1L, round(n / target))
    dbin <- fac * gam_binsize
    grp <- ((seq_len(n) - 1L) %/% fac) + 1L
    if (fac > 1L) {
      # COLOUR FROM THE SIGNIFICANT CELLS, TRANSPARENCY FROM HOW MANY THERE WERE.
      # Two wrong ways to aggregate a sparse significant map were tried first.
      # Mean over significant cells alone fills every block at whole-chromosome
      # scale, because a 27x27 block almost always contains one. Setting
      # non-significant cells to 0 and averaging over all measured cells divides
      # a handful of delta-z by ~729 and washes the map to near-white.
      # So: the block keeps the mean delta-z of its significant cells, which is
      # the direction of change and should read at full strength, and the
      # FRACTION significant is carried as alpha. Sparse blocks fade; blocks with
      # a consistent, dense change are solid.
      if (dz_sig && !is.null(MEAS)) {
        sg <- rowsum(is.finite(D) * 1, grp); sg <- t(rowsum(t(sg), grp))
        ms <- rowsum(MEAS * 1, grp);         ms <- t(rowsum(t(ms), grp))
        FRAC <- sg / pmax(ms, 1); FRAC[ms == 0] <- NA
      } else FRAC <- NULL
      Z <- D; nn <- is.finite(Z); Z[!nn] <- 0
      S1 <- rowsum(Z, grp); N1 <- rowsum(nn * 1, grp)
      S2 <- t(rowsum(t(S1), grp)); N2 <- t(rowsum(t(N1), grp))
      D <- S2 / pmax(N2, 1); D[N2 == 0] <- NA
      agg1 <- function(M) { Zz <- M; nz <- is.finite(Zz); Zz[!nz] <- 0
        a1 <- rowsum(Zz, grp); b1 <- rowsum(nz * 1, grp)
        a2 <- t(rowsum(t(a1), grp)); b2 <- t(rowsum(t(b1), grp))
        R <- a2 / pmax(b2, 1); R[b2 == 0] <- NA; R }
      W0 <- agg1(W0); K0 <- agg1(K0)
      if (!is.null(MEAS)) {
        mz <- rowsum(MEAS * 1, grp); mz <- t(rowsum(t(mz), grp))
        MEAS <- mz > 0            # a display block counts as measured if any cell was
      }
      pos <- as.numeric(tapply(pos, grp, mean))
      rm(Z, nn, S1, N1, S2, N2); gc(FALSE)
    }
    # Light 2D smooth of the DIFFERENCE map. Differencing two noisy matrices
    # leaves a diagonal streak texture that reads as structure but is not; a 3x3
    # box filter suppresses it while preserving domain-scale features. Display
    # only - every statistic elsewhere uses the unsmoothed values. (The Pombo
    # 2026 preprint applies a Gaussian filter for the same reason.)
    if (smooth_map && nrow(D) > 9 && !dz_sig) {
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
    if (dz_sig && !is.null(MEAS)) {
      ij[, meas := MEAS[cbind(i, j)]]
      ij <- ij[meas %in% TRUE]            # drop unmeasured: those stay white
      # alpha only where blocks were aggregated; at native resolution every
      # drawn cell is one significant contact and there is no fraction to show.
      ij[, al := if (is.null(FRAC)) 1 else FRAC[cbind(i, j)]]
      if (!is.null(FRAC)) {
        q <- stats::quantile(ij$al[ij$al > 0], 0.90, na.rm = TRUE)
        ij[, al := pmax(0.12, pmin(1, al / max(q, 1e-9)))]
      }
      ij[!is.finite(al), al := 0.12]
    } else { ij <- ij[is.finite(v)]; ij[, al := 1] }
    lim <- if (identical(diff_mode, "deltaz")) gam12_dz_lim
           else as.numeric(quantile(abs(ij$v[is.finite(ij$v)]), 0.99, na.rm = TRUE))
    tw <- dbin / 1e6   # tiles overlap 2x on the rotated lattice, which fills it

    p_map <- ggplot(ij, aes(x, y, fill = pmax(-lim, pmin(lim, v)), alpha = al)) +
      geom_tile(width = tw, height = tw) +
      scale_alpha_identity() +
      scale_fill_gradientn(colours = diff_cols, limits = c(-lim, lim),
                           na.value = if (dz_sig) "grey88" else "transparent",
                           name = if (identical(diff_mode, "deltaz")) "KO - WT\ndelta-z" else "KO - WT\nNPMI",
                           guide = guide_colourbar(barheight = grid::unit(2.8, "cm"))) +
      { if (gam12_hotline / 2 < maxdist / 2)
          geom_hline(yintercept = gam12_hotline / 2e6, linetype = 2, linewidth = .3,
                     colour = "grey30") } +
      { if (gam12_hotline / 2 < maxdist / 2)
          annotate("text", x = from / 1e6, y = gam12_hotline / 2e6, hjust = -0.04, vjust = -0.5,
                   size = 2.3, colour = "grey30", label = "max hotspot distance") } +
      scale_y_continuous(expand = c(0, 0), name = "Genomic\ndistance\n(Mb)") +
      scale_x_continuous(expand = c(0, 0)) +
      coord_fixed(ratio = 1, xlim = c(from, to) / 1e6, expand = FALSE) +
      # THE FILTER IS NAMED ON THE PANEL, not only in the caption. A masked map is
      # a map of a THRESHOLD as much as of the data, and a reader who cannot see
      # which threshold cannot tell whether a white cell is "not different",
      # "not different enough" or "too little contact to call". The numbers come
      # from gam_z / gam_npmi in config.R, which are also what build the
      # Diff_contacts filenames, so the panel cannot claim a cutoff other than
      # the one that produced the file it just read. The criterion is the
      # collaborators', not ours - stated in the caption, which has room for the
      # provenance this line does not.
      labs(subtitle = if (identical(diff_mode, "deltaz"))
             if (!dz_sig) "delta-z, cKO - WT - every measured cell drawn"
             else sprintf("delta-z, cKO - WT, significant contacts only: |delta-z| > %s and max NPMI > %s%s",
                     format(gam_z), format(gam_npmi),
                     if (fac > 1L) "  (colour = mean delta-z of the significant cells; faded = few of them)"
                     else "  (grey = measured, not significant)")
           else "cKO - WT NPMI") +
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
                          # NAMED "bins", NOT "hotspot", because this repo has a
                          # SECOND thing called a hotspot - the q90 smoothed,
                          # gap-merged DOMAINS that carry every hotspot result in
                          # the results tables. The two agree at a Jaccard of only
                          # 0.39 (lost) / 0.40 (gained): the q90 cut lands at a raw
                          # count of ~44/42, close to 50, so the divergence is the
                          # smoothing and the 5-bin size floor, not the height of
                          # the cut. A reader seeing "hotspot" on the figure will
                          # assume these are the domains from the tables, and
                          # about 60% of them are not.
                          name = sprintf("hotspot bins\n(per-bin count >= %d)", hot_cut),
                          drop = TRUE) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        coord_cartesian(xlim = c(from, to) / 1e6) +
        # One line: a 0.26 in ribbon row can only seat one line of size-7 text
        # before the label collides with the ribbons above and below it, which is
        # what 0.22 in used to do. The "not domains" distinction lives in the
        # legend title and the caption, where there is room to say it properly.
        labs(y = "hotspot bins") +
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

    # --- repeat composition -----------------------------------------------------
    # LINE, SINE and LTR as the fraction of each bin they cover, from
    # gam_bins_annotated.rds (RepeatMasker, already on the 20 kb grid) - so this
    # is a display of an annotation that was always there, not a new measurement.
    #
    # ONE PANEL ON ONE SHARED Y AXIS, not three stacked tracks. The three are
    # fractions of the same bin, and the comparison a reader wants is which class
    # dominates where - the hotspot result is that domains losing contacts are
    # LINE-rich (1.25x) and SINE-DEPLETED (0.80x), which is a statement about
    # LINE against SINE and is invisible if each gets its own auto-scaled axis.
    # This repo has already made that mistake once and recorded it: the mark
    # tracks used to draw a ratio on a per-mark symmetric axis and every mark
    # filled its panel equally whatever the size of its effect. Same trap.
    #
    # The axis ceiling is the largest drawn value across all three classes, so
    # the panel is scaled by what is actually in this region. SINE is genome-wide
    # the smallest (median 0.055 against 0.120 LINE, 0.091 LTR), so a fixed 0-1
    # axis would flatten it to the baseline everywhere.
    p_rep <- if (!isTRUE(repeats)) NULL else local({
      cols <- c("frac_LINE", "frac_SINE", "frac_LTR")
      if (!all(cols %in% names(reg))) {
        message("  repeat tracks: frac_LINE/SINE/LTR not in gam_bins_annotated.rds - omitted")
        return(NULL)
      }
      d <- reg[, c("start", "end", cols), with = FALSE]
      # Same display bin as the mark tracks, so the three panels line up column
      # for column and a feature at one x is at the same x in all of them.
      fac <- max(1L, round(((to - from) / gam12_mark_target_bins) / gam_binsize))
      d[, grp := ((seq_len(.N) - 1L) %/% fac) + 1L]
      d <- d[, .(x = (min(start) + max(end)) / 2e6,
                 LINE = mean(frac_LINE, na.rm = TRUE),
                 SINE = mean(frac_SINE, na.rm = TRUE),
                 LTR  = mean(frac_LTR,  na.rm = TRUE)), by = grp]
      m <- melt(d, id.vars = "x", measure.vars = c("LINE", "SINE", "LTR"),
                variable.name = "class", value.name = "frac")
      m <- m[is.finite(frac)]
      if (!nrow(m)) return(NULL)
      # ROBUST CEILING, the same rule as the mark panels and for the same reason:
      # a single repeat-pileup bin must not flatten the whole region. Twice the
      # 99th percentile of the values actually drawn, never above the maximum
      # (and a fraction cannot exceed 1). The first version of this panel used
      # the plain maximum and Pcdh came out scaled to 0.94 by one LINE spike,
      # with every other feature squashed into the bottom third - which is the
      # failure the mark-panel ceiling exists to prevent.
      top <- min(max(m$frac, na.rm = TRUE),
                 2 * stats::quantile(m$frac, 0.99, na.rm = TRUE), 1)
      if (!is.finite(top) || top <= 0) return(NULL)
      m[, frac := pmin(frac, top)]
      ggplot(m, aes(x, frac, colour = class)) +
        geom_line(linewidth = 0.38) +
        # High-chroma and far apart in hue (~40 deg, 165 deg, 310 deg), because
        # three thin lines in a 0.42 in panel need separation that survives being
        # projected. The ColorBrewer dark endpoints used first (#8c510a, #01665e,
        # #762a83) are built for FILLED areas and read as three similar browns at
        # this linewidth. Deliberately no red and no blue: those two mean
        # direction of change everywhere else in the figure.
        scale_colour_manual(values = c(LINE = "#E69F00", SINE = "#009E73",
                                       LTR = "#AA4499"), name = NULL) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, top * 1.04),
                           breaks = c(0, signif(top, 2))) +
        coord_cartesian(xlim = c(from, to) / 1e6) +
        labs(y = "repeats\n(frac. of bin)") +
        theme_m + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                        axis.text.y = element_text(size = 5),
                        axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5,
                                                    size = 7, lineheight = 0.95),
                        panel.grid = element_blank(), plot.margin = margin(1, 4, 1, 4))
    })

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
      fac <- max(1L, round(((to - from) / gam12_mark_target_bins) / gam_binsize))
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
        # A light touch only. The old window was half a display bin, which at 18
        # display bins per region meant averaging over ~33 of 1200 columns; with
        # the display bin now set from the region width there are enough points
        # that interpolation leaves no corners to round, and anything more than a
        # 3-column mean removes real structure rather than drawing artefacts.
        smw2 <- 3L
        roll <- function(y) { g <- is.finite(y)
          num <- as.numeric(stats::filter(ifelse(g, y, 0), rep(1, smw2), sides = 2))
          den <- as.numeric(stats::filter(g * 1, rep(1, smw2), sides = 2))
          out <- ifelse(is.finite(den) & den > 0, num / den, NA_real_)
          ifelse(g & !is.finite(out), y, out) }
        data.table(x = fx, wt = roll(itp(d$wt)), ko = roll(itp(d$ko)),
                   wt_lo = roll(itp(d$wt_lo)), wt_hi = roll(itp(d$wt_hi)),
                   ko_lo = roll(itp(d$ko_lo)), ko_hi = roll(itp(d$ko_hi)))
      })
      # SHARED Y-SCALE WITHIN A GROUP. H3K9me2 and H3K9me3 are the same kind of
      # mark at different depths, and the whole point of looking at them together
      # is which one moves more - which a per-mark axis destroys by zooming each
      # to fill its own panel. They now share one scale, so the comparison
      # between them is a comparison of the drawing. Marks outside a group keep
      # their own axis, because H3K4me3 is punctate and ten times higher than
      # H3K9me2 and one scale across all five would flatten the broad marks.
      grp_marks <- if (mk %in% gam12_mark_scale_groups[[1]]) gam12_mark_scale_groups[[1]] else mk
      # The max has to come from the AGGREGATED values, not the raw 20 kb bins:
      # what is drawn is a 100 kb mean, and a raw per-bin maximum is several
      # times higher, which would set an axis the curve never reaches and squash
      # it into the floor (H3K4me3 went to 40 against a drawn maximum of ~8).
      grp_vals <- unlist(lapply(grp_marks, function(m) {
        cc <- paste0("sig_", m, "_", c("WT", "HP1gKO"))
        if (!all(cc %in% names(reg))) return(NA_real_)
        g <- ((seq_len(nrow(reg)) - 1L) %/% fac) + 1L
        c(tapply(reg[[cc[1]]], g, mean, na.rm = TRUE),
          tapply(reg[[cc[2]]], g, mean, na.rm = TRUE))
      }))
      ymax <- max(grp_vals, na.rm = TRUE)
      if (!is.finite(ymax) || ymax <= 0) ymax <- 1
      # Robust ceiling (see gam12_mark_clip_* in the config block). Applied only
      # when it gains real height, so regions without an outlier are drawn
      # exactly as before.
      ycap <- ymax
      if (is.finite(gam12_mark_clip_k)) {
        q <- stats::quantile(grp_vals, gam12_mark_clip_q, na.rm = TRUE)
        cand <- as.numeric(gam12_mark_clip_k * q)
        if (is.finite(cand) && cand > 0 && ymax / cand >= gam12_mark_clip_min) ycap <- cand
      }
      clipped <- ycap < ymax
      rng <- c(0, ycap * 1.04)
      # the largest round number that still fits inside the panel, so the top
      # break is always drawn
      ybr <- local({
        cand <- pretty(c(0, ycap), n = 3)
        cand <- cand[cand > 0 & cand <= ycap]
        if (!length(cand)) c(0, signif(ycap, 2)) else c(0, max(cand))
      })
      # Spans of the drawn curve that run off the ceiling, as rectangles half a
      # column wide so adjacent columns join into one continuous bar.
      clip_rect <- if (clipped) local({
        hi <- pmax(fine$wt, fine$ko)
        i  <- which(is.finite(hi) & hi > ycap)
        if (!length(i)) NULL else {
          hw <- (fine$x[2] - fine$x[1]) / 2
          data.table(xmin = fine$x[i] - hw, xmax = fine$x[i] + hw,
                     ymin = ycap * 0.955, ymax = ycap * 1.04)
        }
      }) else NULL

      ggplot(fine) +
        # replicate-range bands first, so the means read on top of them
        { if (gam12_mark_bands)
            geom_ribbon(aes(x = x, ymin = wt_lo, ymax = wt_hi), fill = geno_colors[["WT"]],
                        alpha = .20, na.rm = TRUE) } +
        { if (gam12_mark_bands)
            geom_ribbon(aes(x = x, ymin = ko_lo, ymax = ko_hi), fill = geno_colors[["HP1gKO"]],
                        alpha = .20, na.rm = TRUE) } +
        # Grey body from zero to the LOWER of the two means: the signal that is
        # present in both genotypes, drawn the way a browser draws a coverage
        # track. The coloured gap then sits on top of it, so the figure reads as
        # "this much signal, and this much of a difference" rather than as two
        # free-floating lines.
        geom_ribbon(aes(x = x, ymin = 0, ymax = pmin(wt, ko)),
                    fill = "grey72", colour = NA, na.rm = TRUE) +
        # the gap between the means, coloured by direction. Two ribbons each
        # clipped at the other mean, for the same reason the old panels used two:
        # one ribbon with a mapped fill joins its own points straight across
        # wherever the sign flips.
        geom_ribbon(aes(x = x, ymin = wt, ymax = pmax(ko, wt), fill = "up"), na.rm = TRUE) +
        geom_ribbon(aes(x = x, ymin = pmin(ko, wt), ymax = wt, fill = "down"), na.rm = TRUE) +
        # No lines on the coverage tracks. The grey body tops out at the lower of
        # the two means and the coloured band at the higher, so both means are
        # already the edges of a filled area - drawing them again added two more
        # strokes per panel, ten across the stack, for no information.
        # off-scale marker: a dark bar along the top wherever the curve is cut
        { if (!is.null(clip_rect))
            geom_rect(data = clip_rect,
                      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                      fill = "grey25", colour = NA, inherit.aes = FALSE) } +
        scale_fill_manual(values = c(up = "#b2182b", down = "#2166ac"),
                          breaks = c("up", "down"),
                          labels = c("more in cKO", "less in cKO"), name = NULL) +
        scale_x_continuous(expand = c(0, 0)) +
        # Two breaks, not three: a 0.525 in panel cannot seat three numbers
        # without them touching. They are placed explicitly rather than by
        # n.breaks, because with expand = c(0, 0) an automatic upper break
        # usually lands above the panel and is dropped, leaving panels labelled
        # only "0" - the reader then cannot tell how high the track goes.
        scale_y_continuous(breaks = ybr, expand = c(0, 0)) +
        coord_cartesian(xlim = c(from, to) / 1e6, ylim = rng, expand = FALSE) +
        labs(y = if (clipped)
                    sprintf("%s\ncoverage\n(peak %s)", mk,
                            formatC(ymax, format = "f",
                                    digits = if (ymax >= 10) 0 else 1))
                  else sprintf("%s\ncoverage", mk),
             x = if (last) sprintf("%s (Mb)", chrom) else NULL) +
        theme_m + theme(legend.position = if (show_key) "right" else "none",
                        legend.key.size = unit(9, "pt"), legend.text = element_text(size = 7),
                        # 2pt top and bottom instead of 1: at the reduced panel
                        # height the y-axis numbers of adjacent tracks were
                        # running into each other.
                        plot.margin = margin(2, 4, if (last) 3 else 2, 4),
                        axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5,
                                                    size = 7.5, lineheight = 0.95),
                        axis.text.y = element_text(size = 6),
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
    npmi_cap <- if (!show_raw) NULL else if (!is.null(npmi_cap)) npmi_cap else local({
      kk <- CJ(i = seq_len(nrow(W0)), j = seq_len(nrow(W0)))[j > i]
      kk <- kk[pos[j] - pos[i] <= maxdist]
      v  <- c(W0[cbind(kk$i, kk$j)], K0[cbind(kk$i, kk$j)])
      as.numeric(quantile(v[is.finite(v)], gam12_npmi_qcap, na.rm = TRUE))
    })
    tri <- if (show_raw) list(raw_tri(W0, pos, "WT", npmi_cap),
                              raw_tri(K0, pos, "HP1gKO", npmi_cap), p_map)
           else list(p_map)

    # --- insulation-score map strips --------------------------------------------
    # Block-averaged to the SAME display column count as the contact map above, so
    # a feature sits at the same x in both. Averaging rather than subsampling for
    # the reason given at the map: taking every k-th viewpoint aliases a smooth
    # field into stripes. The grid is then completed over position x square size
    # so that viewpoints with no full square (within one square size of a
    # chromosome end) read as grey missing data instead of as a white gap.
    ismap_panels <- local({
      ismap <- match.arg(ismap, gam12_ismap_modes)
      if (identical(ismap, "none")) return(NULL)
      isd <- load_meltron_is(chrom)
      if (is.null(isd) || !nrow(isd)) return(NULL)
      isd <- isd[viewpoint >= from & viewpoint < to]
      if (!nrow(isd)) {
        message("  insulation maps: no viewpoints in this region - strips omitted")
        return(NULL)
      }
      vps  <- sort(unique(isd$viewpoint))
      ifac <- max(1L, round(length(vps) / gam12_ismap_cols))
      idb  <- ifac * gam_binsize
      # Group by GENOMIC COORDINATE, not by rank among the viewpoints present, so
      # a gap in coverage does not shift everything after it.
      isd[, g := as.integer((viewpoint - from) %/% idb)]
      agg <- isd[, .(IS_WT = mean(IS_WT, na.rm = TRUE),
                     IS_KO = mean(IS_HP1gKO, na.rm = TRUE)),
                 by = .(g, is_size)]
      agg <- agg[CJ(g = 0L:as.integer((to - from - 1) %/% idb),
                    is_size = gam12_is_env$sizes), on = .(g, is_size)]
      agg[, `:=`(x = (from + (g + 0.5) * idb) / 1e6,
                 d = IS_KO - IS_WT)]
      agg[!is.finite(IS_WT), IS_WT := NA_real_]
      agg[!is.finite(IS_KO), IS_KO := NA_real_]
      agg[!is.finite(d),     d     := NA_real_]

      sz   <- gam12_is_env$sizes
      hgt  <- if (length(sz) > 1) min(diff(sz)) / 1000 else 80
      strip <- function(col, ylab, lim, leg) {
        ggplot(agg, aes(x, is_size / 1000,
                        fill = pmax(-lim, pmin(lim, get(col))))) +
          geom_tile(width = idb / 1e6, height = hgt) +
          # COMPACT COLOURBAR. The strips add two more collected legends to a
          # right-hand column that already held three, and at the default bar
          # height the stack is taller than a 16:9 canvas - the bottom legend
          # (LADs) was being cut off. Shortening only these two keeps all five.
          scale_fill_gradientn(colours = diff_cols, limits = c(-lim, lim),
                               na.value = "grey92", name = leg,
                               guide = guide_colourbar(
                                 barheight = grid::unit(0.85, "in"),
                                 barwidth  = grid::unit(0.15, "in"))) +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_continuous(expand = c(0, 0),
                             breaks = range(sz) / 1000) +
          coord_cartesian(xlim = c(from, to) / 1e6) +
          labs(y = ylab) +
          # 5 pt of top/bottom margin, not 1: with expand = c(0, 0) the two
          # square-size breaks sit exactly ON the panel edge, and at 1 pt the
          # top one was drawn half into the contact map above it.
          theme_m + theme(axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          axis.text.y = element_text(size = 5),
                          axis.ticks.y = element_line(linewidth = 0.2,
                                                      colour = "grey60"),
                          panel.grid = element_blank(),
                          plot.margin = margin(5, 4, 5, 4)) + ylab_h
      }
      # One shared limit and one legend title across whichever strips are drawn,
      # so patchwork's guides = "collect" gives them ONE key - which is also the
      # only way they can be read against each other. The difference getting its
      # own tighter limit is the bug fixed above; do not reintroduce it.
      #
      # The y tick labels are bare numbers at size 5, so the unit is named once,
      # in the TOP strip's left-hand label, rather than on every axis.
      lab <- function(main, top) if (top) paste0(main, "\n(square size, kb)") else main
      one <- identical(ismap, "diff")
      c(if (!one) list(strip("IS_WT", lab("insulation, WT", TRUE),  gam12_is_lim, "insulation\nscore"),
                       strip("IS_KO", lab("insulation, cKO", FALSE), gam12_is_lim, "insulation\nscore")),
        list(strip("d", lab("insulation\ncKO - WT", one), gam12_is_dlim, "insulation\nscore")))
    })
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
    ismap_h <- if (is.null(ismap_panels)) numeric(0)
               else rep(gam12_h_ismap, length(ismap_panels))
    rep_h   <- if (is.null(p_rep)) numeric(0) else gam12_h_repeat
    track_h <- c(ismap_h, hot_h, ctcf_h, lad_h, rep_h,
                 rep(gam12_h_mark, length(panels)), gam12_h_gene)
    heights <- grid::unit(c(rep(map_h_in, n_tri), track_h), "in")
    # guides = "collect" is essential here, not cosmetic. Without it every panel
    # keeps its own legend, and patchwork centres each one on ITS OWN ROW - so the
    # LAD legend (two keys) was being centred on a 0.22 in row and spilled over
    # the legends above and below it. Collected and top-justified, they stack.
    wrap_plots(c(tri, ismap_panels, hot_panels,
                 if (is.null(p_ctcf)) NULL else list(p_ctcf),
                 if (is.null(p_lad))  NULL else list(p_lad),
                 if (is.null(p_rep))  NULL else list(p_rep),
                 panels, list(p_gene)), ncol = 1, heights = heights) +
      plot_layout(guides = "collect") +
      plot_annotation(
        title = sprintf("%s:%.2f-%.2f Mb", chrom, from / 1e6, to / 1e6),
        theme = theme(plot.caption = element_text(hjust = 0, size = 7),
                      plot.title = element_text(face = "bold", size = 13),
                      legend.spacing.y = grid::unit(if (slide) 6 else 12, "pt"),
                      legend.box.spacing = grid::unit(10, "pt")),
        caption = if (slide) sprintf(
          # ONE LINE, AND IT MUST FIT. ggplot draws plot.caption as a single
          # line and never wraps it, so a caption longer than the canvas is
          # silently cut off mid-word - which is what the longer version with
          # the insulation sentence spelled out separately did. The colour
          # convention is the same for the map and the strips (red = more
          # contact in the cKO), so one clause covers both.
          # ONE LINE, AND IT MUST FIT. ggplot never wraps plot.caption, so every
          # addition here costs characters somewhere else - adding the repeat
          # sentence overflowed it and cut the measured fraction mid-number.
          # Keep this under ~210 characters at 13.33 in.
          "GAM %s, cKO vs WT, 20 kb. Map: delta-z (cKO - WT)%s; red = more contact in cKO%s. Repeats: LINE/SINE/LTR fraction per bin. Marks: mean coverage. NPMI measured: %.0f%%.",
          chrom,
          if (is.null(ismap_panels)) "" else
            sprintf("; %s: insulation by square size (240-1040 kb)",
                    if (length(ismap_panels) > 1) "strips" else "strip"),
          if (is.null(ismap_panels)) "" else " (condensing)",
          100 * npmi_cov) else paste0(
        if (npmi_cov < gam12_cov_warn)
          sprintf("DATA COVERAGE WARNING: only %.0f%% of NPMI cells in this region are measured. White areas of the contact maps are MISSING DATA, not absent contact.\n", 100 * npmi_cov)
        else sprintf("NPMI cells measured in this region: %.0f%%. ", 100 * npmi_cov),
        "GAM: CA1 pyramidal cells, HP1g WT/cKO, 20 kb, |delta-z|>2, max NPMI>0.1. ",
        "ChIP: Emx1-Cre cortex+hippocampus (pyramidal neurons), per-bin log2FC (KO/WT).\n",
        "Contact map block-averaged and 3x3 smoothed for display. ",
        "MARK TRACKS ARE COVERAGE, NOT A RATIO, drawn at a display bin chosen from the region width (~400 bins across it, floored at the ",
        "native 20 kb). Each panel overlays the two genotypes: the solid lines are the mean per-bin coverage ",
        "across that genotype's libraries (6 WT, 4 cKO - cKO rep2 excluded at QC). The shaded gap between them is the effect, coloured ",
        "red where the cKO is higher and blue where it is lower - the same convention as the contact maps above. Read it as: height is ",
        "how much signal is there, and the gap is the size of the change. Whether a gap beats the spread between replicates is answered ",
        "genome-wide in the ChIP section rather than per region (gap/spread: H4K20me3 0.84, H3K9me3 0.54, H3K9me2 0.30, H3K36me3 0.12, ",
        "H3K4me3 0.04).\n",
        "This replaced a log2 fold-change panel, which drew a ratio on a per-mark symmetric axis so every mark filled its panel equally ",
        "whatever the size of its effect - H3K9me2 (genome-wide median -0.18) read as dramatically as H4K20me3 (-0.50). The genome-wide ",
        "medians are in the ChIP section of the report; this figure is for seeing where signal actually is.\n",
        "REPEAT PANEL: LINE, SINE and LTR as the fraction of each bin they cover (RepeatMasker, on the same 20 kb grid), all three ",
        "on ONE shared y axis so that which class dominates where is readable - the hotspot result is that domains losing contacts are ",
        "LINE-rich (1.25x) and SINE-DEPLETED (0.80x), a comparison that disappears if each class gets its own auto-scaled axis. The ",
        "ceiling is twice the 99th percentile of the values drawn, capped at 1, so one repeat-pileup bin cannot flatten the region.\n",
        "THE Y SCALE IS PER MARK and starts at zero. H3K4me3 is punctate and high, H3K9me2 broad and low, so one shared scale would ",
        "flatten the punctate marks to a line (H3K9me2 and H3K9me3 are the exception and share one axis, so that the comparison ",
        "between them is a comparison of the drawing). Compare marks by the gap relative to their own height, not across panels.\n",
        "OFF-SCALE SPANS. The ceiling of a mark panel is twice the 99th percentile of the values drawn in it, not their maximum, so a ",
        "single repeat-pileup bin cannot flatten a whole region (chr11:3.04-3.16 Mb, the first assembled sequence on chr11 and against ",
        "the pericentromeric satellite, reaches H3K9me3 67.8 where the rest of chr11:1-25 Mb has a 99th percentile of 3.7). The ceiling ",
        "is only applied where it gains at least 25% of panel height, so most regions are drawn to their true maximum and nothing is ",
        "clipped. Where it is applied, the cut spans carry a dark bar along the top of the panel and the true peak is printed in the ",
        "axis title.\n",
        "Diff. contacts = count of significant differential contacts per 20 kb bin within 5 Mb (the dashed line on the map); ",
        "net = gained - lost, same sign convention as the difference map. Those are COUNTS, every bin has one. ",
        "HOTSPOT BINS is the separate CALL on that count: bins at or above the collaborators' proposed 50-75 per-window ",
        "cutoff, drawn at 50 (the permissive end; 75 would drop 73% of them). These are PER-BIN calls on the delivered count ",
        "and are NOT the hotspot DOMAINS the results tables use - those are q90 of a 5-bin smoothed track, gap-merged, 5-bin ",
        "minimum, and the two sets agree at a Jaccard of only 0.39 (lost) / 0.40 (gained). ",
        "LADs: constitutive cLAD/ciLAD, Peric-Hupkes 2010 (GSE17051), lifted mm9->mm39. ",
        "CTCF: ENCODE ENCSR877MSN conservative IDR peaks, BULK hippocampus (not sorted CA1), lifted mm10->mm39; one tick per peak.",
        if (is.null(ismap_panels)) "" else paste0(
        "\nINSULATION STRIPS (GAM_2c, the Pombo display). Each row is the mean NPMI in the square [v-s, v) x (v, v+s] crossing position v, ",
        "at square sizes s from 240 to 1040 kb, with the chromosome mean at that size subtracted. Red = the cKO has MORE contact density ",
        "across that position at that scale (condensing), blue = less (melting) - the same sign convention as the difference map above. ",
        "A vertical band of one colour is a change present at every scale; one confined to the bottom rows is local reorganisation only. ",
        "The strips are block-averaged to the same display columns as the map, and grey is missing data (no full square fits there)."
        ))) &
      theme(legend.justification = "top") -> out
    # Exact figure height: every row plus the title/caption/axis furniture.
    chrome <- if (slide) gam12_h_chrome_slide(length(ismap_panels))
              else gam12_h_chrome
    attr(out, "gam_map_in") <- list(n = n_tri, h = map_h_in,
                                    total = n_tri * map_h_in + sum(track_h) + chrome)
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
    list("chr2",  83.9e6, 91.5e6,  "chr2_Olfr"),        # 224 Olfr genes, 84.91-90.48
    # UNBIASED ADDITIONS, picked off gam6_large_blocks.tsv by width rather than
    # by candidate gene. Two gaps in the set above:
    #   - chr13:71-97 is the LARGEST block in the genome (25.8 Mb, lost) and had
    #     no figure at all. chr11:2.7-23.3 (20.6 Mb) was being treated as the
    #     headline lost block; this one is bigger.
    #   - every block region here was a LOST block. The report's own result is
    #     that lost domains are LAD/LINE/quiescent while GAINED domains are
    #     indistinguishable from rotated background - an asymmetry asserted in
    #     the tables and never once shown. chr3:143-159.7 is the widest gained
    #     block (15.8 Mb) and runs to the end of the chromosome.
    list("chr13", 70.0e6, 98.0e6,  "chr13_71-97Mb_lostblock"),
    list("chr3",  143.0e6, NULL,   "chr3_144-160Mb_gainedblock"),
    # The next three by block width, so the figure set covers the top of the
    # unbiased ranking rather than stopping at whatever was picked first.
    list("chr2",  136.0e6, 158.0e6, "chr2_137-157Mb_lostblock"),
    list("chr5",  46.0e6,  67.0e6,  "chr5_47-66Mb_lostblock"),
    list("chr10", 94.0e6,  111.5e6, "chr10_96-110Mb_gainedblock")
  )
  FIGW <- 11; SIDE <- 2.6
  # Spot-render one region without re-rendering all of them:
  #   GAM7_REGIONS=chr18_Pcdh Rscript GAM_7_Figures.R
  # Each region loads NPMI matrices, so rendering the whole list to check one
  # layout change is an hour of I/O. Comma-separated, matched against the region
  # name; unset means all of them.
  # EVERY VARIANT, IN ONE PASS. A region used to need three separate runs with
  # three different environment variables to produce its three renders, and the
  # masked set silently kept a stale palette for half a day because of exactly
  # that - a re-render of "all the figures" was not all of them. The variants are
  # now a list this loop walks, so one run leaves a consistent set on disk.
  #
  #   <name>            masked difference map, full triangle
  #   <name>_band       masked, clipped to a band (skipped when the full
  #                     triangle is already compact - see band_maxdist)
  #   <name>_alldz      every measured cell drawn, full triangle
  #   <name>_band_alldz every measured cell drawn, band
  #   <name>_slide      masked, 16:9, re-laid-out for presenting
  #
  # GAM7_VARIANTS restricts it (comma-separated, e.g. "slide" or "full,slide").
  # GAM7_SLIDE=1 and GAM7_DZ_SIG are kept as shorthands for the same thing so
  # existing invocations still work.
  all_variants <- list(
    list(id = "full",       mode = "full", dz = TRUE,  slide = FALSE),
    list(id = "band",       mode = "band", dz = TRUE,  slide = FALSE),
    list(id = "alldz",      mode = "full", dz = FALSE, slide = FALSE),
    list(id = "band_alldz", mode = "band", dz = FALSE, slide = FALSE),
    list(id = "slide",      mode = "full", dz = TRUE,  slide = TRUE))
  variants <- all_variants
  if (identical(Sys.getenv("GAM7_SLIDE", "0"), "1"))
    variants <- Filter(function(v) v$slide, all_variants)
  if (nzchar(Sys.getenv("GAM7_DZ_SIG", "")))
    variants <- Filter(function(v) identical(v$dz, !identical(Sys.getenv("GAM7_DZ_SIG"), "0")), variants)
  vwant <- Sys.getenv("GAM7_VARIANTS", "")
  if (nzchar(vwant)) {
    vv <- trimws(strsplit(vwant, ",")[[1]])
    variants <- Filter(function(v) v$id %in% vv, all_variants)
  }
  # Exported BEFORE the GAM7_REGIONS filter: section A2 ranks the whole region
  # set against the genome, and a spot-render of one region must not silently
  # redraw that figure with one bar in it.
  gam12_region_list <<- regions
  only <- Sys.getenv("GAM7_REGIONS", "")
  if (nzchar(only)) {
    want <- trimws(strsplit(only, ",")[[1]])
    regions <- Filter(function(r) r[[4]] %in% want, regions)
    message("GAM7_REGIONS set: rendering ", length(regions), " of the region list")
  }
  # Only some regions are worth a slide; the rest would be a 16:9 render nobody
  # asked for. A region named explicitly overrides that, because naming it is
  # the request.
  slide_regions <- c("chr11_1-25Mb_lostblock", "chr11_100-125Mb", "chr16_full",
                     "chr3_30-85Mb", "chr18_Pcdh",
                     "chr13_71-97Mb_lostblock", "chr3_144-160Mb_gainedblock",
                     "chr2_137-157Mb_lostblock", "chr5_47-66Mb_lostblock",
                     "chr10_96-110Mb_gainedblock")
  message("Rendering ", length(regions), " region(s) x ", length(variants),
          " variant(s): ", paste(vapply(variants, `[[`, "", "id"), collapse = ", "))
  for (r in regions) {
    for (v in variants) {
      if (v$slide && !nzchar(only) && !(r[[4]] %in% slide_regions)) next
      nm <- r[[4]]
      if (v$mode == "band") nm <- paste0(nm, "_band")
      if (!v$dz)            nm <- paste0(nm, "_alldz")
      if (v$slide)          nm <- paste0(nm, "_slide")
      md <- NULL
      if (v$slide) {
        # Solve the band distance FROM the 16:9 constraint instead of picking one.
        # A full triangle's map row is half the panel width tall, which alone
        # busts the budget. Total height is 1 map row + tracks + chrome, and we
        # want it to equal width * 9/16. So: map_h = target - tracks - chrome,
        # then invert band_maxdist()'s geometry, md = 2 * (map_h / panel_w) * span.
        lim   <- region_span(r[[1]], r[[2]], r[[3]])
        span  <- lim[2] - lim[1]
        # Every row that gam_region_plot() will actually draw has to be counted
        # here, or the solved band and the saved height disagree and the figure
        # overflows its canvas.
        trk   <- 0.75 + ismap_n(gam12_ismap_mode) * gam12_h_ismap +
                 gam12_h_hotcall + gam12_h_ctcf + gam12_h_lad + gam12_h_repeat +
                 2 * gam12_h_mark + gam12_h_gene
        pw    <- max(2, 13.33 - SIDE)
        map_h <- max(1.2, 13.33 * 9 / 16 - trk -
                          gam12_h_chrome_slide(ismap_n(gam12_ismap_mode)))
        md    <- min(2 * (map_h / pw) * span, span)
      } else if (v$mode == "band") {
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
      p <- tryCatch(
             if (v$slide)
               gam_region_plot(r[[1]], r[[2]], r[[3]], maxdist = md,
                               fig_width = 13.33, side_allow = SIDE,
                               show_raw = FALSE, dz_sig = v$dz,
                               marks = c("H3K9me3", "H4K20me3"),
                               slide = TRUE, ismap = gam12_ismap_mode)
             else
               gam_region_plot(r[[1]], r[[2]], r[[3]], maxdist = md,
                               fig_width = FIGW, side_allow = SIDE,
                               dz_sig = v$dz, npmi_cap = gam12_npmi_cap),
             error = function(e) { message("  skipped: ", conditionMessage(e)); NULL })
      if (!is.null(p)) {
        a <- attr(p, "gam_map_in")
        h <- if (is.null(a)) 14 else a$total
        save_fig(p, nm, subdir = "regions",
                 width = if (v$slide) 13.33 else FIGW, height = h)
      }
    }
  }
  cat("\nGAM_12 complete. gam_region_plot() is reusable for any region;",
      "pass maxdist to set the band height.\n")
})


# ==============================================================================
# SECTION A2 - the plotted regions, ranked, against a size-matched null
# ==============================================================================
# WHY THIS EXISTS. There are now fifteen region figures and no way to tell from
# the set which of them is actually extreme. A bare ranking of our own picks
# would be circular - we chose them, so of course they rank - so each region is
# ALSO placed against every same-width window in the genome. That percentile is
# the honest number: it is what makes "Pcdh decompacts" and "43.7% of 580 kb
# windows lose more short-range contact than Pcdh does" both true statements,
# and it is the one a collaborator will ask for.
#
# The statistic is mean NET differential contacts per bin (gained - lost), from
# the delivered per-bin counts. Signed, so direction is visible; per bin, so a
# 25 Mb block and a 1.8 Mb cluster are comparable.
local({
  fig_sub <- "regions"
  msg <- function(...) cat(sprintf(...), "\n", sep = "")
  # Self-contained: section A's X and region_span() live inside its own local(),
  # so this block loads what it needs rather than reaching into that scope.
  g1 <- file.path(rds_dir, "gam_bins_annotated.rds")
  if (!file.exists(g1) || !exists("gam12_region_list")) {
    msg("section A2 skipped - gam_bins_annotated.rds or the region list is missing")
    return(invisible())
  }
  A2 <- as.data.table(readRDS(g1)$annot); setorder(A2, chr, start)
  A2[, net := n_contacts_gained - n_contacts_lost]
  rspan <- function(chrom, from = NULL, to = NULL) {
    d <- A2[chr == chrom]
    c(if (is.null(from)) min(d$start) else from,
      if (is.null(to))   max(d$end)   else to)
  }
  R <- rbindlist(lapply(gam12_region_list, function(r) {
    lim <- rspan(r[[1]], r[[2]], r[[3]])
    d <- A2[chr == r[[1]] & start >= lim[1] & end <= lim[2]]
    if (!nrow(d)) return(NULL)
    nb <- nrow(d)
    # Size-matched null: every contiguous nb-bin window on every autosome,
    # stepped by a tenth of the window so the windows are not all independent
    # but do cover the genome. frollmean over the per-chromosome ordering gives
    # the same thing far more cheaply than an explicit loop.
    step <- max(1L, nb %/% 10L)
    nullv <- unlist(lapply(split(A2$net, A2$chr), function(v) {
      if (length(v) < nb) return(numeric(0))
      m <- frollmean(v, nb)
      m <- m[seq(nb, length(m), by = step)]
      m[is.finite(m)]
    }), use.names = FALSE)
    obs <- mean(d$net)
    data.table(name = r[[4]], chr = r[[1]],
               Mb = sprintf("%.0f-%.0f", lim[1] / 1e6, lim[2] / 1e6),
               width_Mb = (lim[2] - lim[1]) / 1e6, n_bins = nb,
               mean_net = obs,
               mean_lost = mean(d$n_contacts_lost),
               mean_gained = mean(d$n_contacts_gained),
               pct = 100 * mean(nullv < obs), n_null = length(nullv))
  }))
  if (!nrow(R)) { msg("section A2 skipped - no regions resolved"); return(invisible()) }
  fwrite(R, file.path(tables_dir, "gam12_region_ranking.tsv"), sep = "\t")
  R[, lab := sprintf("%s  (%s Mb)", name, Mb)]
  R[, dirn := fifelse(mean_net > 0, "gained in cKO", "lost in cKO")]
  setorder(R, mean_net)
  R[, lab := factor(lab, levels = lab)]
  # Two panels on one row of labels: the effect, and where it sits in the genome.
  # Separate panels rather than one axis because they are different units and a
  # dual axis would invite reading one off the other.
  p1 <- ggplot(R, aes(mean_net, lab, fill = dirn)) +
    geom_col(width = 0.68) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
    scale_fill_manual(values = c(`gained in cKO` = "#B2182B",
                                 `lost in cKO` = "#2166AC"), name = NULL) +
    labs(x = "mean net differential contacts per 20 kb bin (gained - lost)", y = NULL,
         subtitle = "Effect") +
    theme_minimal(base_size = 9) +
    theme(panel.grid.major.y = element_blank(), legend.position = "top")
  p2 <- ggplot(R, aes(pct, lab)) +
    geom_vline(xintercept = c(5, 50, 95), linetype = 2, linewidth = 0.25,
               colour = "grey70") +
    geom_segment(aes(x = 50, xend = pct, yend = lab), linewidth = 0.4,
                 colour = "grey55") +
    geom_point(aes(colour = dirn), size = 2.2) +
    scale_colour_manual(values = c(`gained in cKO` = "#B2182B",
                                   `lost in cKO` = "#2166AC"), guide = "none") +
    scale_x_continuous(limits = c(0, 100), breaks = c(0, 5, 25, 50, 75, 95, 100)) +
    labs(x = "percentile among same-width windows genome-wide", y = NULL,
         subtitle = "How extreme is that, really?") +
    theme_minimal(base_size = 9) +
    theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())
  pp <- (p1 | p2) + plot_layout(widths = c(1.35, 1)) +
    plot_annotation(
      title = "The plotted regions, ranked by magnitude of contact change",
      subtitle = paste("Mean net differential contacts per 20 kb bin. Right panel places each region against every",
                       "same-width window in the genome - a region we picked is not automatically an extreme one."),
      caption = paste0("Counts are the delivered per-bin differential contacts within 5 Mb. Null windows are stepped by a tenth of the",
                       " region width, so neighbouring windows overlap; this is a reference distribution, not an independent sample.",
                       "\nTable: gam12_region_ranking.tsv"),
      theme = theme(plot.title = element_text(face = "bold", size = 12),
                    plot.subtitle = element_text(size = 8.5),
                    plot.caption = element_text(hjust = 0, size = 7, colour = "grey35")))
  save_fig(pp, "region_ranking", subdir = fig_sub, width = 11,
           height = max(4.2, 0.30 * nrow(R) + 2.2))
  msg("  region_ranking: %d regions, null n = %d", nrow(R), R$n_null[1])
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
