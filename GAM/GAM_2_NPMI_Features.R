# ================================================================
# GAM_2 - NPMI per-bin features: insulation, distance-banded contacts, decay
# ----------------------------------------------------------------
# The HEAVY, run-once stage of the GAM pipeline. Parses every dense NPMI matrix
# (42 files, ~7.7 GB gzipped) ONE CHROMOSOME x ONE GENOTYPE at a time and reduces
# each to a small per-bin feature table, so later stages never touch the matrices:
#   * insulation score at several window sizes (one primary, rest recorded)
#   * per-bin short-range / long-range mean NPMI + log-log contact-decay slope
#   * the genome-wide contact-decay curve (mean NPMI per diagonal)
# Persists results/rds/npmi_bin_features.rds (the ONLY consumer interface) and
# per-chromosome reductions under results/rds/gam2_cache/ (resumability).
#
# NO A/B COMPARTMENTS. This stage used to call them (O-E, correlation,
# eigenvector selection, saddles) and that was removed on 2026-11-20: NPMI is
# already a normalised quantity, so an eigendecomposition of it does not
# recover compartments - its leading component was per-bin detectability. The
# collaborating lab confirmed this is the wrong input. Do not re-add it on these
# matrices; see CLAUDE.md "A/B compartments are not called".
#
# Needs only the NPMI matrices. Run from the GAM/ directory.
# ================================================================
source("config.R")

# ---------------------------------------------------------------
# Parameters local to this stage (config.R is not edited by this stage; if any
# of these become shared they should move there).
# ---------------------------------------------------------------
gam2_genos       <- c("WT", "HP1gKO")      # canonical; config maps HP1gKO -> cKO on disk
gam2_ins_kb      <- c("100kb" = 100, "200kb" = 200, "500kb" = 500)
gam2_ins_primary <- "200kb"
# Distance bands. SHORT = 1-50 bins (20 kb - 1 Mb): TAD / sub-TAD / loop scale,
# where contact is dominated by polymer proximity and boundary structure.
# LONG = >= 250 bins (>= 5 Mb): beyond any TAD, so contacts here report
# compartment- and territory-level mixing. 5 Mb is also gam_hotspot_maxdist, the
# window the differential hotspot counts were built on, so the two read together.
# The 1-5 Mb gap is deliberately left out so the bands do not share a transition zone.
gam2_short_max_bins <- 50L
gam2_long_min_bins  <- as.integer(5e6 / gam_binsize)
# Decay slope is fit over log-spaced distance shells from 40 kb to 10 Mb. Below
# 40 kb the first bins are dominated by the GAM window-size (sigma) smoothing;
# beyond 10 Mb few bins per chromosome have any signal and shells go noisy.
gam2_shell_edges <- local({
  e <- unique(round(2 * 1.5^(0:20)))
  c(e[e < 500], 500)
})
gam2_min_frac    <- 0.5   # min fraction of non-NA cells for an insulation window
gam2_cache_dir   <- file.path(rds_dir, "gam2_cache")
dir.create(gam2_cache_dir, showWarnings = FALSE, recursive = TRUE)
gam2_out         <- file.path(rds_dir, "npmi_bin_features.rds")

# Chromosome override, for cheap single-chromosome validation only:
#   GAM_2_CHROMS=chr19 Rscript GAM_2_NPMI_Features.R
gam2_chroms <- local({
  e <- Sys.getenv("GAM_2_CHROMS", unset = "")
  if (nzchar(e)) strsplit(e, ",")[[1]] else gam_npmi_chroms
})
stopifnot(all(gam2_chroms %in% gam_npmi_chroms))

# Anything that changes the maths invalidates cached reductions. v9
# (2026-11-20): compartment calling removed, so the cached reductions no longer
# carry pc1/saddle - the bump makes that visible rather than leaving old caches.
gam2_params <- list(version = 9L, binsize = gam_binsize, sigma = gam_sigma,
                    ins_kb = gam2_ins_kb, short_max = gam2_short_max_bins,
                    long_min = gam2_long_min_bins, shells = gam2_shell_edges,
                    min_frac = gam2_min_frac)

# ---------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------
# Sum/count of a vector over distance ranges [a, b] (in bins) around position i,
# both flanks, via cumulative sums: O(n) per column instead of O(n * shells).
# cs/cn are cumsum(c(0, x_or_0)) / cumsum(c(0, !is.na(x))).
flank_sums <- function(cs, cn, i, a, b, n) {
  dlo <- pmax(i + a, 1L); dhi <- pmin(i + b, n)          # downstream flank
  ulo <- pmax(i - b, 1L); uhi <- pmin(i - a, n)          # upstream flank
  dv <- dhi >= dlo; uv <- uhi >= ulo
  s <- ifelse(dv, cs[pmax(dhi, 0L) + 1L] - cs[dlo], 0) +
       ifelse(uv, cs[pmax(uhi, 0L) + 1L] - cs[ulo], 0)
  k <- ifelse(dv, cn[pmax(dhi, 0L) + 1L] - cn[dlo], 0) +
       ifelse(uv, cn[pmax(uhi, 0L) + 1L] - cn[ulo], 0)
  list(sum = s, n = k)
}

# ---------------------------------------------------------------
# Reduce ONE matrix (one genotype, one chromosome) to per-bin features.
# The matrix is loaded here and dies here: nothing n x n escapes this function,
# so at most one matrix is ever alive. chr1 is ~9.7k^2 ~ 760 MB as doubles and
# this function holds one copy, and WT and HP1gKO are
# deliberately NOT held together - the caller loops genotypes sequentially and
# gc()s between them.
# ---------------------------------------------------------------
reduce_npmi <- function(genotype, chr) {
  m  <- load_npmi(genotype, chr)
  n  <- nrow(m); ids <- rownames(m)
  stopifnot(n == ncol(m), identical(ids, colnames(m)))
  start0 <- as.integer(sub("^.+:([0-9]+)-([0-9]+)$", "\\1", ids))
  stopifnot(!is.unsorted(start0))

  # Unmappable bins are all-NA rows/cols in the source. They must STAY NA all the
  # way through: a 0 would read as "never co-segregates".
  nvalid <- colSums(!is.na(m))
  V <- nvalid > 0
  idx <- seq_len(n)

  # --- distance-banded NPMI + decay slope (raw NPMI) -------------------------
  # Per column i, cumulative sums give mean NPMI over any distance range cheaply.
  ne <- length(gam2_shell_edges)
  sh_a <- gam2_shell_edges[-ne]; sh_b <- gam2_shell_edges[-1] - 1L
  sh_d <- vapply(seq_along(sh_a), function(k) mean(sh_a[k]:sh_b[k]), 0) * gam_binsize
  a_all <- c(1L, gam2_long_min_bins, sh_a)
  b_all <- c(gam2_short_max_bins, n, sh_b)
  short <- long <- slope <- rep(NA_real_, n)
  lx <- log10(sh_d)
  for (i in which(V)) {
    v <- m[, i]; ok <- !is.na(v)
    cs <- c(0, cumsum(ifelse(ok, v, 0))); cn <- c(0, cumsum(ok))
    fs <- flank_sums(cs, cn, i, a_all, b_all, n)
    mu <- ifelse(fs$n > 0, fs$sum / fs$n, NA_real_)
    short[i] <- mu[1]; long[i] <- mu[2]
    sm <- mu[-(1:2)]; sn <- fs$n[-(1:2)]
    # log-log needs a positive mean: shells with NPMI <= 0 (noise floor) are
    # dropped rather than clipped; <5 usable shells -> NA, not a poor fit.
    use <- !is.na(sm) & sm > 0 & sn >= 3
    if (sum(use) >= 5) {
      x <- lx[use]; y <- log10(sm[use])
      slope[i] <- sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2)
    }
  }
  short[!is.finite(short)] <- NA; long[!is.finite(long)] <- NA

  # --- insulation (raw NPMI): mean NPMI across the boundary at bin i ----------
  # Square just off the diagonal: rows i-w..i-1 x cols i+1..i+w (Crane et al.).
  # LOW = few contacts across i = insulating. Centred on the chromosome median:
  # the absolute NPMI level depends on cell number / detection and is not
  # comparable between genotypes, the local contrast is. Edge bins without a
  # full window are NA rather than a truncated window.
  ins <- lapply(round(gam2_ins_kb * 1000 / gam_binsize), function(w) {
    out <- rep(NA_real_, n)
    if (n > 2 * w) for (i in (w + 1L):(n - w)) {
      if (!V[i]) next
      blk <- m[(i - w):(i - 1L), (i + 1L):(i + w)]
      ok <- !is.na(blk)
      if (mean(ok) >= gam2_min_frac) out[i] <- mean(blk[ok])
    }
    out - stats::median(out, na.rm = TRUE)
  })

  # --- contact-decay curve: mean NPMI on each diagonal ------------------------
  expd <- rep(NA_real_, n - 1L)
  for (d in seq_len(n - 1L)) {
    r <- seq_len(n - d)
    v <- m[(r + d - 1) * n + r]                # linear index of (row r, col r+d)
    if (any(!is.na(v))) expd[d] <- mean(v, na.rm = TRUE)
  }
  dd <- unique(pmin(round(10^seq(0, log10(n - 1), length.out = 80)), n - 1))
  curve <- data.table(chr = chr, genotype = genotype, dist_bp = dd * gam_binsize,
                      mean_npmi = expd[dd])
  rm(m); invisible(gc())

  feat <- data.table(bin_id = ids, chr = chr,
                     npmi_short = short, npmi_long = long, decay_slope = slope,
                     n_valid = as.integer(nvalid))
  for (w in names(ins)) set(feat, j = paste0("ins_", w), value = ins[[w]])
  list(feat = feat, curve = curve)
}

cache_file <- function(chr, g) file.path(gam2_cache_dir, sprintf("%s_%s.rds", chr, g))
read_cache <- function(chr, g) {
  f <- cache_file(chr, g)
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f)
  if (identical(x$params, gam2_params)) x else NULL
}

# ---------------------------------------------------------------
# Per-chromosome loop (resumable).
# Cache holds ONLY the small per-bin reduction (+ a thinned decay curve)
# : KB-MB per chromosome and genotype. The raw NPMI matrices
# are NEVER written to disk: summed over chromosomes and
# both genotypes they come to 13+ GB, which this machine's disk budget cannot
# take, and they are cheap to rebuild next to their size - only the one-matrix-
# at-a-time parse is slow, and that is exactly what this cache skips on a rerun.
# Do not "improve" this by caching the matrices.
# ---------------------------------------------------------------
for (chr in gam2_chroms) {
  for (g in gam2_genos) {
    if (!is.null(read_cache(chr, g))) { message(chr, " ", g, ": cached"); next }
    message(chr, " ", g, ": parsing ", format(Sys.time(), "%H:%M:%S"))
    res <- reduce_npmi(g, chr)
    res$params <- gam2_params
    saveRDS(res, cache_file(chr, g))
    # Release this chromosome's matrix before the next parse.
    rm(res); invisible(gc())
  }
}

# ---------------------------------------------------------------
# Assemble $features from the per-chromosome reductions
# ---------------------------------------------------------------
caches <- lapply(setNames(gam2_genos, gam2_genos), function(g)
  setNames(lapply(gam2_chroms, read_cache, g = g), gam2_chroms))
done <- gam2_chroms[vapply(gam2_chroms, function(ch) all(vapply(gam2_genos, function(g) !is.null(caches[[g]][[ch]]), TRUE)), TRUE)]
if (!length(done)) stop("no chromosome has both genotypes cached")

prim_col <- paste0("ins_", gam2_ins_primary)
per_geno <- lapply(setNames(gam2_genos, gam2_genos), function(g) {
  rbindlist(lapply(caches[[g]][done], `[[`, "feat"))
})
bins <- per_geno$WT[, .(bin_id, chr)]
stopifnot(identical(bins$bin_id, per_geno$HP1gKO$bin_id))   # same grid in both genotypes
bins[, `:=`(start = as.integer(sub("^.+:([0-9]+)-([0-9]+)$", "\\1", bin_id)),
            end   = as.integer(sub("^.+:([0-9]+)-([0-9]+)$", "\\2", bin_id)))]

features <- copy(bins)
for (g in gam2_genos) {
  d <- per_geno[[g]]
  features[, (paste0(c("insulation", "npmi_short", "npmi_long",
                       "decay_slope", "n_valid"), "_", g)) :=
             d[, .(get(prim_col), npmi_short, npmi_long, decay_slope, n_valid)]]
}
features[, insulation_delta := insulation_HP1gKO - insulation_WT]
setkey(features, bin_id)

other_w <- setdiff(names(gam2_ins_kb), gam2_ins_primary)
insulation_windows <- setNames(lapply(other_w, function(w) {
  dt <- data.table(bin_id = per_geno$WT$bin_id)
  for (g in gam2_genos) dt[, (paste0("insulation_", g)) := per_geno[[g]][[paste0("ins_", w)]]]
  attr(dt, "window_bins") <- as.integer(round(gam2_ins_kb[[w]] * 1000 / gam_binsize))
  setkey(dt, bin_id)
}), other_w)

curves <- rbindlist(lapply(gam2_genos, function(g) rbindlist(lapply(caches[[g]][done], `[[`, "curve"))))
genome_curve <- curves[, .(mean_npmi = mean(mean_npmi, na.rm = TRUE)), by = .(genotype, dist_bp)]

meta <- list(window_sizes_kb = gam2_ins_kb, primary_window = gam2_ins_primary,
             band_short_bins = c(1L, gam2_short_max_bins), band_long_min_bins = gam2_long_min_bins,
             decay_shell_edges_bins = gam2_shell_edges,
             chroms_done = done, chroms_requested = gam2_chroms,
             decay_curve = genome_curve,
             binsize = gam_binsize, sigma = gam_sigma, date = Sys.time(),
             # The METHOD version, carried on the combined output and not only
             # on the per-chromosome caches. Without it the output records no
             # trace of how it was computed, which is how a cache survived a
             # method change once and produced a byte-identical table in 60 s.
             # verify_GAM.R checks this against gam2_params$version.
             params = gam2_params)

# Only a full-genome run may overwrite the consumer interface; a single-chromosome
# validation (GAM_2_CHROMS) writes alongside it so it can never masquerade as the real thing.
out_path <- if (setequal(gam2_chroms, gam_npmi_chroms)) gam2_out else
  sub("\\.rds$", paste0("_partial.rds"), gam2_out)
saveRDS(list(features = features, insulation_windows = insulation_windows,
             meta = meta), out_path)
message("Wrote ", out_path, " (", nrow(features), " bins, ", length(done), " chromosomes)")

# ---------------------------------------------------------------
# Figures
# ---------------------------------------------------------------
if (nrow(curves)) {
  curves[, genotype := factor(genotype, levels = gam2_genos)]
  genome_curve[, genotype := factor(genotype, levels = gam2_genos)]
  p_dec <- ggplot() +
    geom_line(data = curves, aes(dist_bp / 1e6, mean_npmi, group = chr), colour = "grey75", linewidth = 0.3) +
    geom_line(data = genome_curve, aes(dist_bp / 1e6, mean_npmi, colour = genotype), linewidth = 1) +
    facet_wrap(~genotype) + scale_x_log10() +
    scale_colour_manual(values = geno_colors, guide = "none") +
    labs(x = "Genomic distance (Mb)", y = "Mean NPMI (diagonal mean)") + theme_m
  save_fig(p_dec, "contact_decay_by_genotype", subdir = "npmi", width = 10, height = 4.5)
}
