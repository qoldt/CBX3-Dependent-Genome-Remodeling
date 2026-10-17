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
# TWO SCORES PER INTERVAL, NOT ONE SIGNED NUMBER
# ----------------------------------------------
# MELTRONIC emits melting_score AND condensation_score per interval, from two
# one-sided KS tests Bonferroni-corrected separately, and this follows it. The
# two tests are NOT complementary: where the two insulation distributions cross
# - one shifted in its lower tail, the other in its upper - both reject, and the
# interval genuinely carries both scores. Collapsing them into one signed value
# hides that case instead of resolving it. `net_score` is derived for display
# and is explicitly secondary.
#
# On this data the two never fire together at the published cutoff (0.0% of
# intervals), so the single-axis view happens to lose nothing here - but that is
# a result, printed on every run, not an assumption.
#
# THE WINDOW IS THE PUBLISHED ONE, AND GETTING IT WRONG BREAKS THE SCORE
# ---------------------------------------------------------------------
# MELTRONIC's intervals are 120 kb wide stepping one matrix bin
# (data/120kb_genomic_intervals.bed in its repo: 120 kb at a 40 kb step, i.e. 3
# bins at their 40 kb resolution). We match the PHYSICAL width, 6 bins at 20 kb.
#
# An earlier version of this stage used a 500 kb / 25-bin window, and that one
# choice broke everything downstream: it pooled 275 nested, overlapping values
# into a test that assumes independence, so -log10(p_adj) saturated at the cap
# and the published |score| > 5 cutoff called 69% of the genome. The conclusion
# drawn at the time - "the published cutoff does not transfer to this data" -
# was wrong, and was a statement about the window rather than about the data.
# With the published window the cutoff calls 26% melting and 26% condensing and
# is perfectly usable.
#
# THE P-VALUE STILL IS NOT A CALIBRATED P-VALUE
# ---------------------------------------------
# Even at the right window the KS test assumes the pooled IS values are
# independent, and they are not: the 11 square sizes at one viewpoint are nested
# inside each other. The published pipeline does this anyway and so does this
# one, so the scores stay comparable - but a score here is a rescaled effect
# size, not evidence at that significance. That matters more for us than for the
# original: we have 2 biological replicates per genotype already merged, with no
# within-genotype variance to test against (CLAUDE.md gotcha 5).
#
# So two things are reported beside every score:
#   - the one-sided KS statistics d_melt and d_cond, the actual effect sizes,
#     which carry no independence assumption at all;
#   - for the region sets, an EMPIRICAL p from size-matched random regions drawn
#     on the same chromosome, which is the number to quote.
#
# WHAT THE TRACK SUPPORTS, AND AT WHAT SCALE
# ------------------------------------------
# Measured into gam2c_meltron_spatial_coherence.tsv on every run: aggregated to
# one point per 120 kb window, the net score changes sign between ~15% of
# adjacent windows, the median run of one direction is 4 windows (480 kb) and
# the longest reaches 69 (8 Mb). So melting and condensing regions are discrete
# multi-megabase blocks rather than a fine-grained alternation - which is the
# scale of the one-directional blocks GAM_4 calls, and those come out melting
# (lost) and condensing (gained) at the floor of the size-matched null here.
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
  # MELTRONIC's own sliding window is 120 kb wide stepping one matrix bin
  # (data/120kb_genomic_intervals.bed in the MELTRONIC repo: 120 kb intervals at
  # a 40 kb step, i.e. 3 bins wide at their 40 kb resolution). We match the
  # PHYSICAL width, which at 20 kb is 6 bins. An earlier version used 25 bins
  # (500 kb) and that was the cause of the saturated scores: it pooled 275
  # values per window into a test assuming independence, where the published
  # window pools 30.
  win_bp        = 120000L,
  smooth_bins   = 3L,       # published: smooth the score over 3 adjacent bins
  call_cut      = 5,        # published |score| cutoff; usable here, see the header
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
  na <- list(p_melt = NA_real_, p_cond = NA_real_, d_melt = NA_real_,
             d_cond = NA_real_, D = NA_real_, D_signed = NA_real_,
             dir = NA_real_, n = sum(ok))
  if (sum(ok) < 10L) return(na)
  a <- x_wt[ok]; b <- x_ko[ok]
  n <- length(a); m <- length(b)

  # The TWO ONE-SIDED KS statistics, from a single ordering pass:
  #   d_cond = sup(F_WT - F_cKO)  - the cKO's values sit HIGHER = condensation
  #   d_melt = sup(F_cKO - F_WT)  - the cKO's values sit LOWER  = melting
  # Computed directly rather than through two ks.test() calls, because this runs
  # once per bin genome-wide. Verified against ks.test() below.
  o   <- order(c(a, b))
  cs  <- cumsum(c(rep(1 / n, n), rep(-1 / m, m))[o])
  d_cond <- max(cs, 0)
  d_melt <- max(-cs, 0)

  # BOTH SCORES ARE REPORTED, NOT ONE SIGNED NUMBER. This follows MELTRONIC,
  # which emits melting_score AND condensation_score per interval, and it is not
  # a formality: the two one-sided tests are NOT complementary. Where the two
  # empirical distributions cross - one shifted in its lower tail and the other
  # in its upper - both tests reject, and the interval genuinely has a melting
  # score and a condensation score at once. Collapsing that to a single signed
  # value hides the case rather than resolving it.
  ne <- as.numeric(n) * as.numeric(m) / (as.numeric(n) + as.numeric(m))
  stopifnot(is.finite(ne))
  pm <- exp(-2 * ne * d_melt^2)
  pc <- exp(-2 * ne * d_cond^2)
  list(p_melt = pm, p_cond = pc, d_melt = d_melt, d_cond = d_cond,
       D = max(d_melt, d_cond), D_signed = d_melt - d_cond,
       dir = if (d_melt >= d_cond) 1 else -1, n = sum(ok))
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
    for (alt in c("less", "greater")) {
      refk <- suppressWarnings(stats::ks.test(aa, bb, alternative = alt,
                                              exact = FALSE))
      got_d <- if (alt == "less") mine$d_melt else mine$d_cond
      got_p <- if (alt == "less") mine$p_melt else mine$p_cond
      stopifnot(abs(got_d - unname(refk$statistic)) < 1e-9)
      stopifnot(abs(got_p - refk$p.value) < 1e-6 * max(1, refk$p.value) + 1e-9)
    }
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
W <- as.integer(gam2c_params$win_bp / gam_binsize)
stopifnot(W >= 3L)
msg("MELTRONIC: %d kb sliding window (%d bins), stepping one bin, as published",
    gam2c_params$win_bp / 1000, W)

half <- W %/% 2L
track <- list()
for (ch in names(WIDE)) {
  w  <- WIDE[[ch]]
  nv <- length(w$vp)
  if (nv < W) next
  centres <- (half + 1L):(nv - half)
  res <- vapply(centres, function(i) {
    s <- score_rows(ch, (i - half):(i + half))
    c(s$p_melt, s$p_cond, s$d_melt, s$d_cond, s$n)
  }, numeric(5))
  track[[ch]] <- data.table(chr = ch, viewpoint = w$vp[centres],
                            p_melt = res[1, ], p_cond = res[2, ],
                            d_melt = res[3, ], d_cond = res[4, ], n_val = res[5, ])
}
TR <- rbindlist(track); rm(track); invisible(gc())
stopifnot(nrow(TR) > 0)

# Bonferroni across every interval scored, SEPARATELY per direction - exactly
# what MELTRONIC does with its two p-value columns.
TR[, `:=`(p_melt_adj = p.adjust(p_melt, method = "bonferroni"),
          p_cond_adj = p.adjust(p_cond, method = "bonferroni"))]
cap300 <- function(p) pmin(-log10(pmax(p, 1e-300)), 300)
TR[, `:=`(melting_score = cap300(p_melt_adj), condensation_score = cap300(p_cond_adj))]

# Smooth over 3 adjacent bins, as published.
setorder(TR, chr, viewpoint)
smooth3 <- function(v, k = gam2c_params$smooth_bins) {
  g   <- is.finite(v)
  num <- as.numeric(stats::filter(ifelse(g, v, 0), rep(1, k), sides = 2))
  den <- as.numeric(stats::filter(g * 1, rep(1, k), sides = 2))
  ifelse(is.finite(den) & den > 0, num / den, NA_real_)
}
TR[, melting_score      := smooth3(melting_score), by = chr]
TR[, condensation_score := smooth3(condensation_score), by = chr]
# A single display axis, derived rather than primary. Both components are kept.
TR[, net_score := melting_score - condensation_score]

TR[, `:=`(melts = melting_score > gam2c_params$call_cut,
          condenses = condensation_score > gam2c_params$call_cut)]
TR[, call := fcase(melts & condenses, "both",
                   melts, "melting",
                   condenses, "condensing",
                   default = "unchanged")]
TR[, start := viewpoint][, end := viewpoint + gam_binsize]
TR[, bin_id := sprintf("%s:%d-%d", chr, start, end)]

# CAN A REGION BE BOTH? Yes, and this is the number that says how often. The two
# one-sided tests are not complementary: where the two empirical distributions
# cross, both reject. Reported rather than hidden behind a signed value.
pct_both <- 100 * TR[, mean(call == "both", na.rm = TRUE)]

gcols <- function(d) d[, .(n = .N,
    melting    = sum(call == "melting",    na.rm = TRUE),
    condensing = sum(call == "condensing", na.rm = TRUE),
    both       = sum(call == "both",       na.rm = TRUE),
    pct_melting    = 100 * mean(call == "melting",    na.rm = TRUE),
    pct_condensing = 100 * mean(call == "condensing", na.rm = TRUE),
    pct_both       = 100 * mean(call == "both",       na.rm = TRUE),
    median_d_melt = median(d_melt, na.rm = TRUE),
    median_d_cond = median(d_cond, na.rm = TRUE))]
gsumm <- rbind(TR[, gcols(.SD), by = chr], cbind(chr = "ALL", gcols(TR)))

# Spatial coherence of the NET score, at one point per window width.
coh <- local({
  dbin <- as.integer(gam2c_params$win_bp)
  a <- TR[is.finite(net_score)]
  a[, grp := (viewpoint %/% dbin)]
  a <- a[, .(m = mean(net_score, na.rm = TRUE)), by = .(chr, grp)]
  setorder(a, chr, grp)
  flips <- a[, mean(diff(sign(m)) != 0), by = chr]
  runs  <- a[, {r <- rle(sign(m))
                .(med = as.numeric(median(r$lengths)),
                  mx  = as.numeric(max(r$lengths)))}, by = chr]
  data.table(display_bin_kb = dbin / 1000,
             pct_sign_change = 100 * mean(flips$V1, na.rm = TRUE),
             median_run_bins = median(runs$med),
             longest_run_bins = max(runs$mx))
})
fwrite(coh, file.path(tables_dir, "gam2c_meltron_spatial_coherence.tsv"), sep = "\t")
fwrite(gsumm, file.path(tables_dir, "gam2c_meltron_genome_summary.tsv"), sep = "\t")

msg("  at the published |score| > %g cutoff: %.1f%% melting, %.1f%% condensing, %.1f%% BOTH",
    gam2c_params$call_cut, gsumm[chr == "ALL", pct_melting],
    gsumm[chr == "ALL", pct_condensing], pct_both)
msg("  median KS D: melting %.3f, condensation %.3f",
    gsumm[chr == "ALL", median_d_melt], gsumm[chr == "ALL", median_d_cond])
msg("  net score: sign changes between %.0f%% of adjacent %d kb windows; runs median %.0f, max %.0f",
    coh$pct_sign_change, coh$display_bin_kb, coh$median_run_bins, coh$longest_run_bins)

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
to_score <- function(s) {
  if (!is.finite(s$p_melt) || !is.finite(s$p_cond)) return(NA_real_)
  min(-log10(max(s$p_melt, 1e-300)), 300) - min(-log10(max(s$p_cond, 1e-300)), 300)
}

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
    d_melt = obs$d_melt, d_cond = obs$d_cond,
    melting_score = min(-log10(max(obs$p_melt, 1e-300)), 300),
    condensation_score = min(-log10(max(obs$p_cond, 1e-300)), 300),
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
  d <- TR[is.finite(net_score), .(chr, start, end, round(net_score, 3))]
  p <- file.path(igv_dir, "gam_melting_score.bedgraph")
  lim <- as.numeric(quantile(abs(d$V4), 0.99, na.rm = TRUE))
  writeLines(sprintf(paste0('track type=bedGraph name="MELTRON net score" ',
                            'description="melting_score minus condensation_score; ',
                            '+ = melting, - = condensing; %d kb window, %d square sizes" ',
                            'visibility=full color=178,24,43 altColor=33,102,172 ',
                            'viewLimits=%.2f:%.2f windowingFunction=mean'),
                     gam2c_params$win_bp / 1000,
                     length(gam2c_params$is_sizes_bp), -lim, lim), p)
  fwrite(d, p, sep = "\t", col.names = FALSE, append = TRUE, quote = FALSE,
         scipen = 999)
  msg("  gam_melting_score.bedgraph (%d rows)", nrow(d))
})

# ---------------------------------------------------------------
# Figures
# ---------------------------------------------------------------
local({
  suppressPackageStartupMessages(library(ggplot2))

  # 1. BOTH SCORES, MIRRORED. This is the honest rendering and it answers the
  #    obvious objection directly: an interval can carry a melting score and a
  #    condensation score at the same time, because the two one-sided tests are
  #    not complementary. Drawing one signed value would hide exactly that.
  dbin <- as.integer(gam2c_params$win_bp)
  show <- intersect(c("chr11", "chr16"), unique(TR$chr))
  if (length(show)) {
    d3 <- TR[chr %in% show & is.finite(melting_score)]
    d3[, grp := (viewpoint %/% dbin) * dbin]
    d3 <- d3[, .(melting = mean(melting_score, na.rm = TRUE),
                 condensation = -mean(condensation_score, na.rm = TRUE)),
             by = .(chr, grp)]
    d3 <- melt(d3, id.vars = c("chr", "grp"), variable.name = "which",
               value.name = "score")
    p3 <- ggplot(d3, aes(grp / 1e6, score, fill = which)) +
      geom_hline(yintercept = c(-1, 1) * gam2c_params$call_cut, linetype = 2,
                 colour = "grey55", linewidth = 0.3) +
      geom_area(colour = NA) +
      scale_fill_manual(values = c(melting = "#b2182b", condensation = "#2166ac"),
                        name = NULL) +
      facet_wrap(~ chr, ncol = 1, scales = "free_x") +
      labs(x = "position (Mb)", y = "-log10 p (Bonferroni)",
           title = "Melting and condensation are two scores, not one signed one",
           subtitle = paste0(dbin / 1000, " kb window stepping one bin, as MELTRONIC. ",
                             "Dashed lines are the published |score| > ",
                             gam2c_params$call_cut, " cutoff.\nAn interval with ",
                             "both means its two insulation distributions CROSS - ",
                             "one tail up, the other down.")) +
      theme_m + theme(legend.position = "top")
    save_fig(p3, "melting_track", subdir = fig_sub, width = 11, height = 5.5)
  }

  # 2. How often is an interval called in BOTH directions at once?
  dd <- TR[is.finite(melting_score) & is.finite(condensation_score)]
  p1 <- ggplot(dd, aes(melting_score, condensation_score)) +
    geom_hline(yintercept = gam2c_params$call_cut, linetype = 2, colour = "grey55") +
    geom_vline(xintercept = gam2c_params$call_cut, linetype = 2, colour = "grey55") +
    geom_point(alpha = 0.05, size = 0.4, colour = "#2166ac") +
    labs(x = "melting score", y = "condensation score",
         title = "The two one-sided tests are not complementary",
         subtitle = sprintf(paste0("Top-right quadrant = called BOTH ways: %.1f%% of ",
                                   "intervals. Those are places where the two\ninsulation ",
                                   "distributions cross, so a single signed track would ",
                                   "misdescribe them."), pct_both)) +
    theme_m
  save_fig(p1, "melting_score_distribution", subdir = fig_sub, width = 6.5, height = 5.5)

  # 3. Score against effect size, per direction.
  d2 <- rbind(dd[, .(D = d_melt, score = melting_score, which = "melting")],
              dd[, .(D = d_cond, score = condensation_score, which = "condensation")])
  p2 <- ggplot(d2, aes(D, score)) +
    geom_point(alpha = 0.05, size = 0.4, colour = "#2166ac") +
    facet_wrap(~ which) +
    labs(x = "one-sided KS statistic D (effect size)", y = "-log10 p (Bonferroni)",
         title = "Score against effect size",
         subtitle = paste0("With the published ", gam2c_params$win_bp / 1000,
                           " kb window the score no longer saturates, so the ",
                           "published cutoff is usable.")) +
    theme_m
  save_fig(p2, "melting_score_vs_D", subdir = fig_sub, width = 8, height = 4.5)

  if (nrow(REG)) {
    d4 <- REG[is.finite(ks_D)]
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

# ============================================================================
# PART D - FIGURES THAT CONNECT MELTING TO THE REST OF THE ANALYSIS
# ----------------------------------------------------------------------------
# The three figures above describe the melting score itself. These relate it to
# everything else the project has measured, which is the only way it earns its
# place: a melting track nobody can connect to the hotspots, the marks or the
# domains is a number without a claim attached.
#
# They follow the published displays rather than inventing new ones. The Pombo
# plotting scripts (pombo-lab/Meltron code/) use (a) per-region ECDF panels -
# the direct picture of what the KS test compares - and (b) ComplexHeatmaps of
# regions x covariates with melting score as one column among compartment, RNA
# and ATAC. The insulation-by-square-size heatmap is the third: it is the only
# view that shows at WHICH SPATIAL SCALE a region changes, which a single
# collapsed score cannot.
# ============================================================================
local({
  suppressPackageStartupMessages(library(ggplot2))

  # ------------------------------------------------------------------
  # D1. INSULATION DIFFERENCE BY SQUARE SIZE - the heatmap-style track.
  # x = position, y = insulation square size, fill = IS(cKO) - IS(WT).
  # A region that melts at every scale is a block of red across the whole
  # panel; one that reorganises only locally is red at 240 kb and flat at
  # 1 Mb. The collapsed score cannot tell those apart.
  # ------------------------------------------------------------------
  for (rg in list(list(chr = "chr18", from = 36.0e6, to = 39.0e6, nm = "chr18_Pcdh"),
                  list(chr = "chr11", from = 3.0e6,  to = 25.0e6, nm = "chr11_lostblock"))) {
    d <- IS[chr == rg$chr & viewpoint >= rg$from & viewpoint < rg$to]
    if (!nrow(d)) next
    d[, diff := IS_HP1gKO - IS_WT]
    lim <- as.numeric(quantile(abs(d$diff), 0.98, na.rm = TRUE))
    p <- ggplot(d, aes(viewpoint / 1e6, is_size / 1000,
                       fill = pmax(-lim, pmin(lim, diff)))) +
      geom_raster() +
      scale_fill_gradientn(colours = c("#2166ac", "#f7f7f7", "#b2182b"),
                           limits = c(-lim, lim), na.value = "grey92",
                           name = "IS(cKO)\n- IS(WT)") +
      scale_y_continuous(breaks = gam2c_params$is_sizes_bp / 1000, expand = c(0, 0)) +
      scale_x_continuous(expand = c(0, 0)) +
      labs(x = sprintf("%s (Mb)", rg$chr), y = "insulation square size (kb)",
           title = sprintf("Where, and at what scale, contact density changes - %s", rg$nm),
           subtitle = paste0("Red = the cKO has MORE local contact density at that scale ",
                             "(condensing); blue = less (melting).\nA vertical band of one ",
                             "colour is a change at every scale; a change confined to the ",
                             "bottom rows is local reorganisation only.")) +
      theme_m + theme(panel.grid = element_blank())
    save_fig(p, paste0("is_by_scale_", rg$nm), subdir = fig_sub, width = 10, height = 4.2)
  }

  # ------------------------------------------------------------------
  # D2. MELTING AGAINST EVERYTHING ELSE WE HAVE MEASURED.
  # Per 20 kb bin, the net melting score against the project's own
  # annotations. This is the figure that says whether melting is telling us
  # anything the other tracks were not.
  # ------------------------------------------------------------------
  g1 <- tryCatch(readRDS(file.path(rds_dir, "gam_bins_annotated.rds")),
                 error = function(e) NULL)
  if (!is.null(g1)) {
    A <- as.data.table(g1$annot)
    J <- merge(TR[, .(bin_id, chr, net_score, melting_score, condensation_score)],
               A, by = "bin_id", suffixes = c("", ".a"))
    ov <- nrow(J) / nrow(TR)
    msg("  joined melting to annotation on %.1f%% of scored bins", 100 * ov)
    stopifnot(ov > 0.5)

    J[, hotspot := fcase(dom_lost_adj > 0 & dom_gained_adj > 0, "both",
                         dom_lost_adj > 0, "lost domain",
                         dom_gained_adj > 0, "gained domain",
                         default = "neither")]
    J[, net_contacts := n_contacts_gained - n_contacts_lost]
    # H4K20me3 per-bin log2FC, the headline ChIP effect
    if ("H4K20me3" %in% colnames(g1$l2fc))
      J[, k20 := g1$l2fc[match(bin_id, rownames(g1$l2fc)), "H4K20me3"]]

    cats <- list()
    cats$hotspot <- J[hotspot != "both", .(grp = "hotspot domain", lvl = hotspot,
                                           v = net_score)]
    if ("compartment_centred" %in% names(J))
      cats$chip <- J[!is.na(compartment_centred),
                     .(grp = "ChIP class", lvl = as.character(compartment_centred),
                       v = net_score)]
    hm <- grep("^hmm_", names(J), value = TRUE)
    if (length(hm)) {
      J[, top_state := sub("^hmm_", "", hm[max.col(as.matrix(.SD), ties.method = "first")]),
        .SDcols = hm]
      keep <- J[, .N, by = top_state][N > 2000, top_state]
      cats$state <- J[top_state %in% keep, .(grp = "ChromHMM state", lvl = top_state,
                                             v = net_score)]
    }
    D2 <- rbindlist(cats)
    if (nrow(D2)) {
      D2[, lvl := reorder(lvl, v, median, na.rm = TRUE)]
      p2 <- ggplot(D2[is.finite(v)], aes(lvl, v)) +
        geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
        geom_boxplot(outlier.shape = NA, fill = "grey88", linewidth = 0.3) +
        facet_wrap(~ grp, scales = "free_x", nrow = 1) +
        coord_cartesian(ylim = as.numeric(quantile(D2$v, c(0.01, 0.99), na.rm = TRUE))) +
        labs(x = NULL, y = "net melting score  (+ melting, - condensing)",
             title = "Melting against the project's own annotations",
             subtitle = paste0("Per 20 kb bin. If melting carried no information beyond ",
                               "the existing tracks these would all sit on zero.")) +
        theme_m + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7))
      save_fig(p2, "melting_vs_annotation", subdir = fig_sub, width = 11, height = 4.6)
    }

    # Continuous covariates: melting against contact change and against the
    # mark loss, binned so ~119k points do not become a black rectangle.
    cv <- list()
    cv$contacts <- J[is.finite(net_contacts),
                     .(x = net_contacts, y = net_score, w = "net differential contacts")]
    if ("k20" %in% names(J))
      cv$k20 <- J[is.finite(k20), .(x = k20, y = net_score, w = "H4K20me3 log2FC")]
    CV <- rbindlist(cv)
    if (nrow(CV)) {
      CV[, xb := cut(x, breaks = quantile(x, seq(0, 1, 0.05), na.rm = TRUE),
                     include.lowest = TRUE, labels = FALSE), by = w]
      SM <- CV[, .(x = median(x, na.rm = TRUE), y = median(y, na.rm = TRUE),
                   lo = quantile(y, 0.25, na.rm = TRUE),
                   hi = quantile(y, 0.75, na.rm = TRUE), n = .N), by = .(w, xb)]
      rho <- CV[, .(r = stats::cor(x, y, method = "spearman",
                                   use = "pairwise.complete.obs")), by = w]
      SM <- merge(SM, rho, by = "w")
      SM[, lab := sprintf("%s  (rho = %.3f)", w, r)]
      p3 <- ggplot(SM, aes(x, y)) +
        geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
        geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey85") +
        geom_line(colour = "#b2182b", linewidth = 0.7) +
        geom_point(size = 1.2, colour = "#b2182b") +
        facet_wrap(~ lab, scales = "free_x") +
        labs(x = NULL, y = "net melting score",
             title = "Melting against the two quantities the project is built on",
             subtitle = paste0("Each point is a 5% bin of the x variable; ribbon is the ",
                               "interquartile range of the melting score in that bin.")) +
        theme_m
      save_fig(p3, "melting_vs_covariates", subdir = fig_sub, width = 9, height = 4.4)
    }
  }

  # ------------------------------------------------------------------
  # D3. ECDF PANELS - the published diagnostic, and the direct picture of
  # what the KS test actually compares. The vertical gap between the two
  # curves IS the statistic.
  # ------------------------------------------------------------------
  if (exists("region_sets") && length(region_sets)) {
    ec <- rbindlist(lapply(names(region_sets), function(nm) {
      R <- region_sets[[nm]]
      pl <- pool_rows(R)
      a <- unlist(lapply(names(pl), function(ch) as.vector(WIDE[[ch]]$wt[pl[[ch]], , drop = FALSE])))
      b <- unlist(lapply(names(pl), function(ch) as.vector(WIDE[[ch]]$ko[pl[[ch]], , drop = FALSE])))
      ok <- is.finite(a) & is.finite(b)
      if (sum(ok) < 50) return(NULL)
      rbind(data.table(set = nm, genotype = "WT", v = a[ok]),
            data.table(set = nm, genotype = "cKO", v = b[ok]))
    }))
    if (nrow(ec)) {
      lab <- REG[, .(set, lab = sprintf("%s\nD = %.3f, emp p = %.3f", set, ks_D,
                                        p_empirical_D))]
      ec <- merge(ec, lab, by = "set")
      p4 <- ggplot(ec, aes(v, colour = genotype)) +
        stat_ecdf(linewidth = 0.8) +
        scale_colour_manual(values = c(WT = "#2166ac", cKO = "#b2182b"), name = NULL) +
        facet_wrap(~ lab, scales = "free_x") +
        labs(x = "insulation score (all square sizes pooled)",
             y = "cumulative probability",
             title = "What the Kolmogorov-Smirnov test actually compares",
             subtitle = paste0("The largest vertical gap between the two curves IS the ",
                               "statistic D. cKO curve to the LEFT = less local\ncontact ",
                               "density = melting. This is the published diagnostic ",
                               "(Meltron code/plot_ecdf.R).")) +
        theme_m + theme(legend.position = "top", strip.text = element_text(size = 7))
      save_fig(p4, "melting_ecdf_region_sets", subdir = fig_sub, width = 10, height = 7)
    }
  }
})

msg("GAM_2c complete.")
