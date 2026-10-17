# ================================================================
# GAM_2c - MELTRON: chromatin melting and condensation
# ----------------------------------------------------------------
#   cd GAM && Rscript GAM_2c_Meltron.R
#   GAM2C_CHROMS=chr19 Rscript GAM_2c_Meltron.R      (spot run)
#
# Needs only the NPMI matrices. Independent of GAM_1 and GAM_2.
#
# WHAT THIS IS
# ------------
# The Pombo lab's own instrument for "does this region hold its contacts
# together or has it come apart", from Winick-Ng et al. Nature 2021 (MELTRON,
# github.com/pombo-lab/Meltron) and its genome-wide extension MELTRONIC. It asks
# a different question from our hotspot track and answers it from a different
# quantity, which is exactly why it is worth having: the hotspot track counts the
# collaborators' per-contact |delta-z| > 2 calls, while melting compares the
# DISTRIBUTION of local contact density across a range of spatial scales. Nothing
# here depends on the delivered differential-contact calls.
#
# THE METHOD, AS PUBLISHED
# ------------------------
# 1. Insulation score (IS) at a viewpoint v and square size s: the mean contact
#    value in the square [v-s, v) x (v, v+s] - the block of contacts that cross
#    v. A bin is dropped if fewer than half the expected (s/res)^2 cells are
#    present. The diagonal and immediately adjacent bins are excluded, which this
#    square already does by construction.
# 2. Repeat over a range of square sizes, so the score describes a region across
#    spatial scales rather than at one arbitrary window.
# 3. For a region, pool its IS values over all viewpoints and all square sizes in
#    each sample and compare the two empirical cumulative distributions with a
#    ONE-SIDED Kolmogorov-Smirnov test; Bonferroni-correct across regions;
#    melting score = -log10(p_adj). Lower contact density in the test sample =
#    melting; higher = condensation, reported as a negative score here so the two
#    live on one axis.
#
# TWO DELIBERATE DEVIATIONS, both forced by the data
# --------------------------------------------------
# (a) IS IS ADDITIVE HERE, NOT A LOG RATIO. The published code computes
#     IS = log2(mean contacts in square) - log2(chromosome mean at that size).
#     That needs a positive quantity. NPMI is already a log-scale quantity
#     (pointwise mutual information) and goes negative, so log2 of it is
#     undefined wherever the square mean is at or below zero. We subtract the
#     chromosome mean at that square size instead. This is the same decision, for
#     the same reason, that GAM_2 makes for observed-minus-expected and that
#     CLAUDE.md records as settled point 8 - and it is why a melting score from
#     this pipeline must not be quoted against a published one as though the
#     underlying scale were the same. The SIGN and the RANKING are comparable;
#     the magnitude is not.
# (b) 20 kb VIEWPOINTS, NOT 40 kb. The published work bins to 40 kb. Our whole
#     analysis grid is 20 kb, so viewpoints at 20 kb join to everything else in
#     the repo without a regrid, and a finer viewpoint spacing is strictly more
#     information. The square sizes are unchanged and are multiples of both.
#
# THE P-VALUE IS NOT A CALIBRATED P-VALUE, AND IS NOT TREATED AS ONE
# ------------------------------------------------------------------
# The KS test assumes the pooled IS values are independent. They are emphatically
# not: adjacent viewpoints share most of their square, and the 11 square sizes at
# one viewpoint are nested inside each other. The published pipeline does this
# anyway and so does this one, so that the scores stay comparable - but a melting
# score here is a REscaled effect size, not evidence at that significance. This
# matters more for us than for the original: we have 2 biological replicates per
# genotype already merged, with no within-genotype variance to test against
# (CLAUDE.md gotcha 5 forbids exactly this kind of claim).
#
# So two things are reported beside every score:
#   - the KS statistic D, the actual effect size, which has no independence
#     assumption in it at all;
#   - for the region sets, an EMPIRICAL p from size-matched random regions drawn
#     on the same chromosome, which is the number to quote.
#
# WHAT THE TRACK SUPPORTS, AND AT WHAT SCALE
# ------------------------------------------
# The MAGNITUDE of the change is large almost everywhere (median |D| 0.42) -
# what a genome-wide effect looks like. The DIRECTION is regionally coherent but
# not finely so, and this is measured on every run into
# gam2c_meltron_spatial_coherence.tsv: aggregated to one point per sliding
# window, the sign changes between about a quarter of adjacent 500 kb windows,
# the median run of one direction is 3 windows (1.5 Mb) and the longest is 28
# (14 Mb).
#
# That scale is not an accident - it is the scale of the one-directional blocks
# GAM_4 calls, and those blocks come out melting (lost) and condensing (gained)
# at the floor of the size-matched null here. So the track carries real
# multi-megabase structure, while a SINGLE 500 kb window's melting/condensing
# label is weak evidence. Read regions, not windows; the genome-wide call
# columns describe the tails of a distribution rather than delineating domains.
#
# Outputs
#   results/rds/meltron_is.rds                     IS, both genotypes, all sizes
#   results/rds/meltron_scores.rds                 per-bin melting track + regions
#   results/tables/gam2c_meltron_regions.tsv       region-level scores + empirical p
#   results/tables/gam2c_meltron_genome_summary.tsv
#   results/tables/gam2c_meltron_spatial_coherence.tsv
#   results/igv/gam_melting_score.bedgraph
#   results/figures/meltron/*.{png,pdf}
# ================================================================
source("config.R")

msg <- function(fmt, ...) message(sprintf(fmt, ...))

gam2c_params <- list(
  version       = 1L,
  # 240-1040 kb in 80 kb steps. The cocaine-paper methods say "10 equidistant
  # sizes ranging from 240-1040 kb in 80 kb increments", which is internally
  # inconsistent - 240 to 1040 by 80 is 11 values, not 10. We take the two
  # endpoints and the stated increment, since those pin the scale range, and
  # every size is a multiple of both 20 kb and 40 kb.
  is_sizes_bp   = seq(240000L, 1040000L, by = 80000L),
  min_cell_frac = 0.5,      # square dropped below half its expected cells
  win_bins      = 25L,      # MELTRONIC sliding window, 25 x 20 kb = 500 kb
  smooth_bins   = 3L,       # published: smooth the score over 3 adjacent bins
  call_cut      = 5,        # published |score| cutoff - reported, NOT used (see below)
  call_frac     = 0.05,     # calls = this tail fraction of |signed KS D| each way
  max_na_frac   = 0.5,      # region skipped above this fraction missing
  n_null        = 500L      # size-matched random regions per region set
)

fig_sub  <- "meltron"
is_rds   <- file.path(rds_dir, "meltron_is.rds")
sc_rds   <- file.path(rds_dir, "meltron_scores.rds")
igv_dir  <- file.path(results_dir, "igv")
dir.create(igv_dir, recursive = TRUE, showWarnings = FALSE)

CHROMS <- local({
  e <- Sys.getenv("GAM2C_CHROMS", "")
  if (!nzchar(e)) return(gam_autosomes)
  w <- trimws(strsplit(e, ",")[[1]])
  bad <- setdiff(w, gam_autosomes)
  if (length(bad)) stop("GAM2C_CHROMS: not autosomes: ", paste(bad, collapse = ", "))
  msg("GAM2C_CHROMS set: %d of %d chromosomes - PARTIAL RUN", length(w),
      length(gam_autosomes))
  w
})
PARTIAL <- !identical(CHROMS, gam_autosomes)

# ---------------------------------------------------------------
# PART A - insulation scores
# ---------------------------------------------------------------
# Squares are computed from a 2D prefix sum (summed-area table) so that every
# square is O(1) regardless of its size. The naive loop is O(n * s^2) per size
# and, with 11 sizes up to 52 bins on each side, does not finish on chr1.
#
# The prefix sum is built IN PLACE on the loaded matrix, and only the cell-count
# prefix needs a second array, so the peak is about two matrices rather than
# four. One genotype at a time, per gotcha 6.
# ---------------------------------------------------------------
prefix2d <- function(M) {
  n <- ncol(M)
  for (j in seq_len(n)) M[, j] <- cumsum(M[, j])      # down each column
  for (j in seq_len(n)[-1]) M[, j] <- M[, j] + M[, j - 1L]  # then across
  M
}

# Sum over rows a..b, cols c..d from a prefix-sum table, vectorised over squares.
# Index 0 means "nothing above/left of the block", hence the zero guards.
rect_sum <- function(P, a, b, c, d) {
  n <- ncol(P)
  g <- function(i, j) {
    out <- numeric(length(i))
    ok <- i >= 1L & j >= 1L
    out[ok] <- P[cbind(i[ok], j[ok])]
    out
  }
  g(b, d) - g(a - 1L, d) - g(b, c - 1L) + g(a - 1L, c - 1L)
}

is_one <- function(genotype, chr) {
  m <- load_npmi(genotype, chr)
  n <- nrow(m); ids <- rownames(m)
  start0 <- as.integer(sub("^.+:([0-9]+)-([0-9]+)$", "\\1", ids))
  stopifnot(!is.unsorted(start0))

  ok <- !is.na(m)
  m[!ok] <- 0
  P <- prefix2d(m);            rm(m);  invisible(gc())
  Q <- prefix2d(ok * 1);       rm(ok); invisible(gc())

  out <- vector("list", length(gam2c_params$is_sizes_bp))
  for (k in seq_along(gam2c_params$is_sizes_bp)) {
    s <- gam2c_params$is_sizes_bp[k]
    w <- as.integer(s / gam_binsize)
    if (n <= 2 * w + 1) next
    v  <- (w + 1L):(n - w)                       # viewpoints with a full square
    sm <- rect_sum(P, v - w, v - 1L, v + 1L, v + w)
    cn <- rect_sum(Q, v - w, v - 1L, v + 1L, v + w)
    mu <- ifelse(cn >= gam2c_params$min_cell_frac * w * w, sm / cn, NA_real_)
    # ADDITIVE centring on the chromosome mean at this square size - see the
    # header for why this is a subtraction and not the published log2 ratio.
    out[[k]] <- data.table(chr = chr, viewpoint = start0[v], is_size = s,
                           IS = mu - mean(mu, na.rm = TRUE))
  }
  rm(P, Q); invisible(gc())
  rbindlist(out)
}

# The cache is only reusable if it covers at least the chromosomes this run
# wants. Silently reusing a partial cache for a full run is the kind of thing
# that produces a confident genome-wide number from one chromosome.
cache_ok <- local({
  if (!file.exists(is_rds) || nzchar(Sys.getenv("GAM2C_FORCE", ""))) return(FALSE)
  o <- readRDS(is_rds)
  if (!all(CHROMS %in% o$chroms)) {
    msg("Cache covers %d chromosome(s), this run needs %d - recomputing",
        length(o$chroms), length(CHROMS))
    return(FALSE)
  }
  if (!identical(o$params$is_sizes_bp, gam2c_params$is_sizes_bp) ||
      !identical(o$params$min_cell_frac, gam2c_params$min_cell_frac)) {
    msg("Cache was built with different IS parameters - recomputing")
    return(FALSE)
  }
  TRUE
})
if (cache_ok) {
  msg("Reusing %s (GAM2C_FORCE=1 to recompute)", is_rds)
  IS <- readRDS(is_rds)$is[chr %in% CHROMS]
} else {
  msg("Insulation scores: %d square sizes, %d-%d kb, at %d kb viewpoints",
      length(gam2c_params$is_sizes_bp), min(gam2c_params$is_sizes_bp) / 1000,
      max(gam2c_params$is_sizes_bp) / 1000, gam_binsize / 1000)
  acc <- list()
  for (chr in CHROMS) {
    msg("  %s", chr)
    a <- is_one("WT", chr); b <- is_one("HP1gKO", chr)
    setnames(a, "IS", "IS_WT"); setnames(b, "IS", "IS_HP1gKO")
    acc[[chr]] <- merge(a, b, by = c("chr", "viewpoint", "is_size"), all = TRUE)
    rm(a, b); invisible(gc())
  }
  IS <- rbindlist(acc); rm(acc); invisible(gc())
  setkey(IS, chr, viewpoint, is_size)
  saveRDS(list(is = IS, params = gam2c_params, partial = PARTIAL,
               chroms = CHROMS, built = Sys.time()), is_rds)
  msg("  wrote %s (%d rows)", is_rds, nrow(IS))
}

msg("  IS finite in both genotypes: %.1f%% of %d rows",
    100 * IS[, mean(is.finite(IS_WT) & is.finite(IS_HP1gKO))], nrow(IS))

# How often would the published log2 form have been defined at all? Reported
# rather than assumed, because it is the justification for deviation (a).
msg("  square means above zero (where log2 would be defined): WT %.1f%%, cKO %.1f%%",
    100 * IS[, mean(IS_WT > 0, na.rm = TRUE)],
    100 * IS[, mean(IS_HP1gKO > 0, na.rm = TRUE)])

# ---------------------------------------------------------------
# Scoring: one-sided KS in BOTH directions, keep the supported one, sign it.
#
# The published tool is run twice, once per direction, producing a melting score
# and a compaction score as separate positive numbers. Taking the smaller of the
# two p-values and signing by direction puts them on one axis, which is what the
# genome-wide version needs, and is the same quantity.
# ---------------------------------------------------------------
meltron_score <- function(x_wt, x_ko) {
  ok <- is.finite(x_wt) & is.finite(x_ko)
  if (sum(ok) < 20L) return(list(p = NA_real_, D = NA_real_, D_signed = NA_real_,
                                 dir = NA_real_, n = sum(ok)))
  a <- x_wt[ok]; b <- x_ko[ok]
  n <- length(a); m <- length(b)

  # The TWO ONE-SIDED KS statistics, from a single ordering pass:
  #   d_plus  = sup(F_WT - F_cKO)  - the cKO's values sit HIGHER  = condensing
  #   d_minus = sup(F_cKO - F_WT)  - the cKO's values sit LOWER   = melting
  # Computed directly rather than through two ks.test() calls, because this runs
  # once per 20 kb bin genome-wide and ks.test() re-sorts and re-allocates each
  # time. Verified against ks.test() - see the equivalence check below.
  o   <- order(c(a, b))
  cs  <- cumsum(c(rep(1 / n, n), rep(-1 / m, m))[o])
  d_plus  <- max(cs, 0)
  d_minus <- max(-cs, 0)

  # THE SIGNED STATISTIC IS THE DIFFERENCE OF THE TWO ONE-SIDED STATISTICS, not
  # the two-sided maximum with a sign attached. That earlier version flipped sign
  # between adjacent windows whenever d_plus and d_minus were both large - which
  # they often are, because two distributions can differ in both tails at once -
  # and produced a track that alternated melting/condensing every few bins. The
  # difference is continuous, lies in [-1, 1], and is zero when the two
  # distributions coincide.
  d_signed <- d_minus - d_plus
  D        <- max(d_plus, d_minus)                 # = the two-sided statistic
  dir      <- if (d_signed >= 0) 1 else -1
  d_one    <- if (dir > 0) d_minus else d_plus

  # Asymptotic one-sided KS p-value: exp(-2 * ne * d^2), ne = nm/(n+m). This is
  # exactly what ks.test() computes once n*m >= 10000, which is where we always
  # are: the smallest thing scored here is a sliding window of 25 viewpoints x
  # 11 square sizes = 275 values per genotype, so n*m = 75625. Below that
  # threshold ks.test() switches to the exact distribution and the two diverge,
  # which is why the check below pins exact = FALSE.
  # as.numeric() is load-bearing: length() returns an integer, and a region can
  # pool millions of values, so n * m overflows 32-bit integer arithmetic and
  # returns NA with a warning. That silently turned every large region's score
  # into NA while the KS statistic beside it stayed correct.
  ne <- as.numeric(n) * as.numeric(m) / (as.numeric(n) + as.numeric(m))
  p  <- exp(-2 * ne * d_one^2)
  stopifnot(is.finite(ne))
  list(p = p, D = D, D_signed = d_signed, dir = dir, n = sum(ok))
}

# EQUIVALENCE CHECK against ks.test(), run once at load. The hand-rolled
# statistic is the load-bearing part of this stage, so it is verified rather
# than asserted - and it is cheap enough to verify on every run.
local({
  set.seed(11)
  for (i in 1:30) {
    # meltron_score() drops positions where either input is NA, so it always
    # compares two equal-length vectors - give the check the same shape.
    L  <- sample(40:400, 1)
    aa <- rnorm(L)
    bb <- rnorm(L, mean = runif(1, -1, 1), sd = runif(1, 0.5, 2))
    mine <- meltron_score(aa, bb)
    ref  <- suppressWarnings(stats::ks.test(aa, bb))
    stopifnot(abs(mine$D - unname(ref$statistic)) < 1e-9)
    alt  <- if (mine$dir > 0) "less" else "greater"
    refp <- suppressWarnings(stats::ks.test(aa, bb, alternative = alt,
                                            exact = FALSE)$p.value)
    stopifnot(abs(mine$p - refp) < 1e-6 * max(1, refp) + 1e-9)
  }
  msg("  KS statistic and one-sided p verified against ks.test() on 30 random cases")
})

# ---------------------------------------------------------------
# One wide representation per chromosome, built once and used by BOTH the
# sliding-window track and the region scoring. Rows are viewpoints in genomic
# order, columns are square sizes, so any interval is a contiguous row range -
# which is what makes the size-matched null affordable (it draws tens of
# thousands of intervals, and each one would otherwise be a keyed subset of a
# multi-million-row table).
# ---------------------------------------------------------------
WIDE <- local({
  out <- list()
  for (ch in unique(IS$chr)) {
    d <- IS[chr == ch]
    A <- dcast(d, viewpoint ~ is_size, value.var = "IS_WT")
    B <- dcast(d, viewpoint ~ is_size, value.var = "IS_HP1gKO")
    stopifnot(identical(A$viewpoint, B$viewpoint))
    out[[ch]] <- list(vp = A$viewpoint,
                      wt = as.matrix(A[, -1, with = FALSE]),
                      ko = as.matrix(B[, -1, with = FALSE]))
  }
  out
})

# Score the rows r of chromosome ch.
score_rows <- function(ch, r) {
  w <- WIDE[[ch]]
  nullres <- list(p = NA_real_, D = NA_real_, D_signed = NA_real_,
                  dir = NA_real_, n = 0L)
  if (is.null(w) || !length(r)) return(nullres)
  r <- r[r >= 1L & r <= length(w$vp)]
  if (!length(r)) return(nullres)
  meltron_score(as.vector(w$wt[r, , drop = FALSE]),
                as.vector(w$ko[r, , drop = FALSE]))
}

# ---------------------------------------------------------------
# PART B - MELTRONIC: a genome-wide melting track
# ---------------------------------------------------------------
msg("MELTRONIC: sliding window of %d bins (%d kb), scored at every bin",
    gam2c_params$win_bins, gam2c_params$win_bins * gam_binsize / 1000)

W <- gam2c_params$win_bins
half <- W %/% 2L
track <- list()
for (ch in names(WIDE)) {
  w  <- WIDE[[ch]]
  nv <- length(w$vp)
  if (nv < W) next
  centres <- (half + 1L):(nv - half)
  # n*m for one window, against the threshold where ks.test() would switch to
  # the exact distribution and our asymptotic p would stop matching it.
  if (ch == names(WIDE)[1]) {
    npw <- W * ncol(w$wt)
    stopifnot(npw * npw >= 10000)
  }
  res <- vapply(centres, function(i) {
    s <- score_rows(ch, (i - half):(i + half))
    c(s$p, s$D, s$D_signed, s$dir, s$n)
  }, numeric(5))
  track[[ch]] <- data.table(chr = ch, viewpoint = w$vp[centres],
                            p_raw = res[1, ], D = res[2, ], D_signed = res[3, ],
                            dir = res[4, ], n_val = res[5, ])
}
TR <- rbindlist(track); rm(track); invisible(gc())
stopifnot(nrow(TR) > 0)

# Bonferroni across every window scored, as the published pipeline does across
# every region tested.
TR[, p_adj := p.adjust(p_raw, method = "bonferroni")]
TR[, score_raw := dir * pmin(-log10(p_adj), 300)]   # cap: p_adj underflows to 0
TR[!is.finite(score_raw), score_raw := NA_real_]

# Smooth over 3 adjacent bins, as published - applied to both quantities so the
# score and the statistic stay directly comparable.
setorder(TR, chr, viewpoint)
smooth3 <- function(v, k = gam2c_params$smooth_bins) {
  g   <- is.finite(v)
  num <- as.numeric(stats::filter(ifelse(g, v, 0), rep(1, k), sides = 2))
  den <- as.numeric(stats::filter(g * 1, rep(1, k), sides = 2))
  ifelse(is.finite(den) & den > 0, num / den, NA_real_)
}
TR[, melting_score := smooth3(score_raw), by = chr]
TR[, D_signed := smooth3(D_signed), by = chr]

# D_signed (computed in meltron_score) is the primary quantity of this track:
# bounded [-1, 1], the signed distance between the two distributions, carrying
# no independence assumption. The melting score is kept beside it for
# comparability with the published work, but it is NOT what the calls are made
# from, for a reason the data settled rather than theory: see below.

# THE PUBLISHED |score| > 5 CUTOFF DOES NOT TRANSFER TO THIS DATA, AND IS NOT USED.
# It calls about a third of the genome in the published setting. Here it calls
# roughly two thirds, because our windows pool many more values (11 square sizes
# x 25 viewpoints at 20 kb) into a test whose p-value assumes they are
# independent, so -log10(p_adj) saturates. The number is printed below so the
# discrepancy is visible rather than buried.
#
# Calls are therefore the tails of the genome's OWN signed-D distribution. This
# is a relative statement - "the most melted 5% of the genome" - and is labelled
# as such, which is the honest form given that no absolute calibration is
# available without within-genotype replicates (gotcha 5).
dcut <- as.numeric(quantile(abs(TR$D_signed), 1 - gam2c_params$call_frac,
                            na.rm = TRUE))
TR[, call := fcase(D_signed >=  dcut, "melting",
                   D_signed <= -dcut, "condensing",
                   default = "unchanged")]
TR[!is.finite(D_signed), call := NA_character_]
pub_frac <- TR[, mean(abs(melting_score) > gam2c_params$call_cut, na.rm = TRUE)]
TR[, start := viewpoint][, end := viewpoint + gam_binsize]
TR[, bin_id := sprintf("%s:%d-%d", chr, start, end)]

gsumm <- TR[, .(n = .N,
                melting = sum(call == "melting", na.rm = TRUE),
                condensing = sum(call == "condensing", na.rm = TRUE),
                pct_melting = 100 * mean(call == "melting", na.rm = TRUE),
                pct_condensing = 100 * mean(call == "condensing", na.rm = TRUE),
                median_D = median(D, na.rm = TRUE)), by = chr]
gsumm <- rbind(gsumm, TR[, .(chr = "ALL", n = .N,
                melting = sum(call == "melting", na.rm = TRUE),
                condensing = sum(call == "condensing", na.rm = TRUE),
                pct_melting = 100 * mean(call == "melting", na.rm = TRUE),
                pct_condensing = 100 * mean(call == "condensing", na.rm = TRUE),
                median_D = median(D, na.rm = TRUE))])
# IS THE DIRECTION SPATIALLY COHERENT? Measured, because the track invites being
# read as a domain caller and the answer decides whether that is legitimate.
# Aggregated to one point per sliding-window width (so overlapping windows are
# not counted repeatedly), how often does the sign change between neighbours,
# and how long are the runs of one sign?
coh <- local({
  dbin <- as.integer(gam2c_params$win_bins * gam_binsize)
  a <- TR[is.finite(D_signed)]
  a[, grp := (viewpoint %/% dbin)]
  a <- a[, .(m = mean(D_signed, na.rm = TRUE)), by = .(chr, grp)]
  setorder(a, chr, grp)
  flips <- a[, mean(diff(sign(m)) != 0), by = chr]
  # as.numeric on both: median() of an odd-length integer run returns an
  # integer and of an even-length one a double, and data.table rejects a column
  # whose type changes between groups.
  runs  <- a[, {r <- rle(sign(m))
                .(med = as.numeric(median(r$lengths)),
                  mx  = as.numeric(max(r$lengths)))}, by = chr]
  data.table(display_bin_kb = dbin / 1000,
             pct_sign_change = 100 * mean(flips$V1, na.rm = TRUE),
             median_run_bins = median(runs$med),
             longest_run_bins = max(runs$mx))
})
fwrite(coh, file.path(tables_dir, "gam2c_meltron_spatial_coherence.tsv"), sep = "\t")
msg("  direction: sign changes between %.0f%% of adjacent %d kb windows; runs median %.0f, max %.0f",
    coh$pct_sign_change, coh$display_bin_kb, coh$median_run_bins, coh$longest_run_bins)
msg("  => regionally coherent over megabases; a single window is weak evidence")

fwrite(gsumm, file.path(tables_dir, "gam2c_meltron_genome_summary.tsv"), sep = "\t")
msg("  calls from the signed-D tails (|D| >= %.3f): %.1f%% melting, %.1f%% condensing",
    dcut, gsumm[chr == "ALL", pct_melting], gsumm[chr == "ALL", pct_condensing])
msg("  median KS D = %.3f", gsumm[chr == "ALL", median_D])
msg("  for reference, the published |score| > %g cutoff would call %.1f%% of bins here",
    gam2c_params$call_cut, 100 * pub_frac)

# ---------------------------------------------------------------
# PART C - region-level MELTRON, with a size-matched empirical null
#
# The KS p-value is not calibrated (see the header), so every region set is also
# scored against random regions of the SAME WIDTH on the SAME CHROMOSOME. The
# empirical p is the fraction of those whose |score| reaches the observed one.
# That is the number to quote.
# ---------------------------------------------------------------
region_sets <- local({
  L <- list()
  add <- function(name, dt) {
    if (is.null(dt) || !nrow(dt)) return(invisible(NULL))
    dt <- as.data.table(dt)[chr %in% names(WIDE) & end > start]
    if (nrow(dt)) L[[name]] <<- dt[, .(chr, start = as.numeric(start),
                                       end = as.numeric(end))]
  }
  tt <- function(f) { p <- file.path(tables_dir, f)
                      if (file.exists(p)) fread(p) else NULL }
  d <- tt("hotspot_domains_lost_adj.tsv");   if (!is.null(d)) add("hotspot_lost", d)
  d <- tt("hotspot_domains_gained_adj.tsv"); if (!is.null(d)) add("hotspot_gained", d)
  d <- tt("gam6_large_blocks.tsv")
  if (is.null(d) && file.exists(file.path(results_dir, "bed", "gam6_large_blocks.bed")))
    d <- fread(file.path(results_dir, "bed", "gam6_large_blocks.bed"))
  if (!is.null(d) && "dirn" %in% names(d)) {
    add("large_blocks_lost",   d[dirn == "lost"])
    add("large_blocks_gained", d[dirn == "gained"])
  }
  # the named target regions, same coordinates GAM_7 plots
  add("Pcdh_cluster",  data.table(chr = "chr18", start = 36.6e6, end = 38.4e6))
  add("KRABZFP_chr13", data.table(chr = "chr13", start = 66.0e6, end = 68.3e6))
  add("KRABZFP_chr17", data.table(chr = "chr17", start = 21.0e6, end = 24.1e6))
  add("Vmn_chr7",      data.table(chr = "chr7",  start = 19.6e6, end = 24.0e6))
  add("Olfr_chr2",     data.table(chr = "chr2",  start = 84.9e6, end = 90.5e6))
  L
})
msg("Region-level MELTRON on %d region set(s)", length(region_sets))

# Rows of WIDE[[ch]] covered by [s, e)
rows_for <- function(ch, s_, e_) {
  w <- WIDE[[ch]]
  if (is.null(w)) return(integer(0))
  which(w$vp >= s_ & w$vp < e_)
}

# Pool the rows of a whole interval set, then score once: a set is the unit of
# interest, and a per-interval score would be dominated by interval width.
pool_rows <- function(R) {
  by_chr <- split(seq_len(nrow(R)), R$chr)
  lapply(by_chr, function(ii)
    sort(unique(unlist(lapply(ii, function(i) rows_for(R$chr[i], R$start[i], R$end[i]))))))
}
score_pool <- function(pl) {
  a <- unlist(lapply(names(pl), function(ch) as.vector(WIDE[[ch]]$wt[pl[[ch]], , drop = FALSE])))
  b <- unlist(lapply(names(pl), function(ch) as.vector(WIDE[[ch]]$ko[pl[[ch]], , drop = FALSE])))
  meltron_score(a, b)
}
to_score <- function(s) if (!is.finite(s$p)) NA_real_ else
  s$dir * min(-log10(max(s$p, 1e-300)), 300)

set.seed(42)
reg_rows <- list()
for (nm in names(region_sets)) {
  R   <- region_sets[[nm]]
  obs <- score_pool(pool_rows(R))
  obs_score <- to_score(obs)

  # SIZE-MATCHED NULL: every interval is redrawn at a random position on its OWN
  # chromosome, keeping each width and the number of intervals exactly. Widths
  # matter here for the same reason they do everywhere else in this repo - a
  # wide interval pools more values and so reaches a larger -log10(p) for the
  # same effect size, which is precisely the artefact the empirical p removes.
  nwid <- vapply(seq_len(nrow(R)), function(i) R$end[i] - R$start[i], 0)
  nullm <- vapply(seq_len(gam2c_params$n_null), function(b) {
    pl <- list()
    for (i in seq_len(nrow(R))) {
      ch <- R$chr[i]; w <- WIDE[[ch]]
      if (is.null(w)) next
      span <- max(w$vp) - min(w$vp)
      if (span <= nwid[i]) next
      st <- min(w$vp) + floor(runif(1) * (span - nwid[i]) / gam_binsize) * gam_binsize
      pl[[ch]] <- c(pl[[ch]], rows_for(ch, st, st + nwid[i]))
    }
    if (!length(pl)) return(c(NA_real_, NA_real_))
    pl <- lapply(pl, function(v) sort(unique(v)))
    sc <- score_pool(pl)
    c(to_score(sc), sc$D)
  }, numeric(2))
  null  <- nullm[1, ]
  nullD <- nullm[2, ]

  nn <- sum(is.finite(null))
  reg_rows[[nm]] <- data.table(
    set = nm, n_intervals = nrow(R), total_kb = sum(nwid) / 1000,
    n_values = obs$n, ks_D = obs$D,
    direction = fifelse(is.na(obs$dir), NA_character_,
                        fifelse(obs$dir > 0, "melting", "condensing")),
    meltron_score = obs_score,
    null_median = median(null, na.rm = TRUE),
    null_q95 = as.numeric(quantile(abs(null), 0.95, na.rm = TRUE)),
    null_D_median = median(nullD, na.rm = TRUE),
    null_n = nn,
    # Two empirical p-values, because the two quantities can disagree and the
    # disagreement is informative: the score rewards pooling more values, the KS
    # statistic does not. p_empirical_D is the one to quote.
    p_empirical = if (nn > 0 && is.finite(obs_score))
      (1 + sum(abs(null) >= abs(obs_score), na.rm = TRUE)) / (nn + 1) else NA_real_,
    p_empirical_D = if (sum(is.finite(nullD)) > 0 && is.finite(obs$D))
      (1 + sum(nullD >= obs$D, na.rm = TRUE)) / (sum(is.finite(nullD)) + 1) else NA_real_)
  msg("  %-20s %-11s D %.3f (null %.3f)  emp p(D) %-7s | score %6.1f  emp p %s",
      nm, reg_rows[[nm]]$direction, obs$D, reg_rows[[nm]]$null_D_median,
      format.pval(reg_rows[[nm]]$p_empirical_D, digits = 2), obs_score,
      format.pval(reg_rows[[nm]]$p_empirical, digits = 2))
}
REG <- rbindlist(reg_rows, fill = TRUE)
if (nrow(REG)) {
  setorder(REG, -meltron_score)
  fwrite(REG, file.path(tables_dir, "gam2c_meltron_regions.tsv"), sep = "\t")
}

saveRDS(list(track = TR, regions = REG, genome_summary = gsumm,
             params = gam2c_params, partial = PARTIAL, chroms = CHROMS,
             built = Sys.time()), sc_rds)
msg("  wrote %s", sc_rds)

# ---------------------------------------------------------------
# IGV track
# ---------------------------------------------------------------
local({
  d <- TR[is.finite(D_signed), .(chr, start, end, round(D_signed, 4))]
  p <- file.path(igv_dir, "gam_melting_score.bedgraph")
  lim <- as.numeric(quantile(d$V4, c(0.01, 0.99), na.rm = TRUE))
  writeLines(sprintf(paste0('track type=bedGraph name="MELTRON signed KS D" ',
                            'description="signed KS D; + = melting (density lost in cKO), ',
                            '- = condensing; %d kb window, %d square sizes" ',
                            'visibility=full color=178,24,43 altColor=33,102,172 ',
                            'viewLimits=%.2f:%.2f windowingFunction=mean'),
                     gam2c_params$win_bins * gam_binsize / 1000,
                     length(gam2c_params$is_sizes_bp), lim[1], lim[2]), p)
  fwrite(d, p, sep = "\t", col.names = FALSE, append = TRUE, quote = FALSE,
         scipen = 999)
  msg("  gam_melting_score.bedgraph (%d rows)", nrow(d))
})

# ---------------------------------------------------------------
# Figures
# ---------------------------------------------------------------
local({
  suppressPackageStartupMessages(library(ggplot2))

  p1 <- ggplot(TR[is.finite(D_signed)], aes(D_signed)) +
    geom_histogram(bins = 120, fill = "grey60") +
    geom_vline(xintercept = c(-1, 1) * dcut, linetype = 2, colour = "#b2182b") +
    labs(x = "signed KS statistic D  (+ melting in cKO, - condensing)",
         y = "20 kb bins",
         title = "Genome-wide melting, measured as an effect size",
         subtitle = sprintf(paste0("Dashed lines are the %.0f%% tails that define a call ",
                                   "(|D| >= %.3f). Median |D| = %.3f.\nThe published ",
                                   "|score| > %g cutoff would call %.0f%% of bins here,\n",
                                   "so it is not used - see the GAM_2c header."),
                            100 * gam2c_params$call_frac, dcut,
                            gsumm[chr == "ALL", median_D], gam2c_params$call_cut,
                            100 * pub_frac)) +
    theme_m
  save_fig(p1, "melting_score_distribution", subdir = fig_sub, width = 8, height = 5)

  # D against the score: how much of the score is effect size and how much is
  # just the number of values pooled.
  p2 <- ggplot(TR[is.finite(melting_score) & is.finite(D)], aes(D, abs(melting_score))) +
    geom_point(alpha = 0.05, size = 0.4, colour = "#2166ac") +
    labs(x = "KS statistic D (effect size, no independence assumption)",
         y = "|melting score|",
         title = "Why the calls come from D and not from the score",
         subtitle = paste0("Where the cloud is vertical, windows with the SAME distance ",
                           "between distributions get wildly different\nscores - the score is ",
                           "responding to how many values were pooled, not to the effect.")) +
    theme_m
  save_fig(p2, "melting_score_vs_D", subdir = fig_sub, width = 6.5, height = 5.5)

  show <- intersect(c("chr11", "chr16"), unique(TR$chr))
  if (length(show)) {
    # AGGREGATE FOR DISPLAY, to the window size itself. The track is smooth -
    # lag-1 autocorrelation 0.9989 between adjacent 20 kb bins - but a
    # chromosome holds several thousand of them, consecutive windows overlap by
    # 24 of 25 bins, and nothing finer than the window is resolved anyway.
    # Drawing every bin therefore renders one measurement thousands of times and
    # aliases into what looks like noise. One point per window width is the
    # honest display resolution.
    dbin <- as.integer(gam2c_params$win_bins * gam_binsize)
    d3 <- TR[chr %in% show & is.finite(D_signed)]
    d3[, grp := (viewpoint %/% dbin) * dbin]
    d3 <- d3[, .(D_signed = mean(D_signed, na.rm = TRUE)), by = .(chr, grp)]
    p3 <- ggplot(d3, aes(grp / 1e6, D_signed)) +
      geom_hline(yintercept = c(-1, 1) * dcut, linetype = 2,
                 colour = "grey60", linewidth = 0.3) +
      geom_area(aes(fill = D_signed > 0), colour = NA) +
      scale_fill_manual(values = c("TRUE" = "#b2182b", "FALSE" = "#2166ac"),
                        labels = c("TRUE" = "melting", "FALSE" = "condensing"),
                        name = NULL) +
      facet_wrap(~ chr, ncol = 1, scales = "free_x") +
      labs(x = "position (Mb)", y = "signed KS D",
           title = "Melting along the chromosome",
           subtitle = paste0(dbin / 1000, " kb display bins (one per sliding window); ",
                             "dashed lines are the genome-wide ",
                             round(100 * gam2c_params$call_frac), "% tails")) +
      theme_m + theme(legend.position = "top")
    save_fig(p3, "melting_track", subdir = fig_sub, width = 11, height = 5)
  }

  if (nrow(REG)) {
    d4 <- REG[is.finite(ks_D)]
    # ORDERED BY EXCESS OVER ITS OWN NULL, not by raw D. A single multi-Mb
    # interval reaches a sizeable D by itself, so ranking on D alone puts the
    # region sets that are indistinguishable from background at the top of the
    # chart - which is the opposite of what the chart is for.
    d4[, excess := ks_D - null_D_median]
    d4[, set_f := reorder(set, excess)]
    p4 <- ggplot(d4, aes(set_f, ks_D, fill = direction)) +
      geom_col() +
      geom_point(aes(y = null_D_median), shape = 4, size = 2.4, colour = "grey20") +
      geom_text(aes(y = pmax(ks_D, null_D_median) + 0.012,
                    label = sprintf("p=%.3f", p_empirical_D)),
                hjust = 0, size = 2.9, colour = "grey25") +
      scale_fill_manual(values = c("melting" = "#b2182b", "condensing" = "#2166ac"),
                        name = NULL, na.value = "grey70") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
      coord_flip() +
      labs(x = NULL, y = "KS statistic D (effect size)",
           title = "Region sets against size-matched random regions",
           subtitle = paste0("x marks the median D of ", gam2c_params$n_null,
                             " size-matched random region sets,\nsame widths and ",
                             "same chromosomes. A bar at its x is background.")) +
      theme_m + theme(legend.position = "top")
    save_fig(p4, "melting_region_sets", subdir = fig_sub, width = 8, height = 5.5)
  }
})

msg("GAM_2c complete.")
