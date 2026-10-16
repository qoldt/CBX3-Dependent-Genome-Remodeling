# ================================================================
# GAM_2 - NPMI compartments: A/B calling, insulation, distance-banded contacts
# ----------------------------------------------------------------
# The HEAVY, run-once stage of the GAM pipeline. Parses every dense NPMI matrix
# (42 files, ~7.7 GB gzipped) ONE CHROMOSOME x ONE GENOTYPE at a time and reduces
# each to a small per-bin feature table, so GAM_3-5 never touch the matrices:
#   * observed - expected(distance) normalisation, Pearson correlation, PC1
#   * PC1 oriented per chromosome so positive = A = H3K4me3-high (WT reference)
#   * insulation score at several window sizes (one primary, rest recorded)
#   * per-bin short-range / long-range mean NPMI + log-log contact-decay slope
#   * saddle matrix per genotype (strata of PC1 quantile)
# Persists results/rds/npmi_bin_features.rds (the ONLY consumer interface) and
# per-chromosome reductions under results/rds/gam2_cache/ (resumability).
#
# Needs GAM_1's results/rds/gam_bins_annotated.rds for the H3K4me3 orientation.
# Runs after GAM_1; run from the GAM/ directory.
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
gam2_K           <- 20L   # saddle strata
gam2_min_frac    <- 0.5   # min fraction of non-NA cells for an insulation window
gam2_cache_dir   <- file.path(rds_dir, "gam2_cache")
dir.create(gam2_cache_dir, showWarnings = FALSE, recursive = TRUE)
gam2_out         <- file.path(rds_dir, "npmi_bin_features.rds")

# Chromosome override, for cheap single-chromosome validation only:
#   GAM_2_CHROMS=chr19 Rscript GAM_2_NPMI_Compartments.R
gam2_chroms <- local({
  e <- Sys.getenv("GAM_2_CHROMS", unset = "")
  if (nzchar(e)) strsplit(e, ",")[[1]] else gam_npmi_chroms
})
stopifnot(all(gam2_chroms %in% gam_npmi_chroms))

# Anything that changes the maths invalidates cached reductions.
gam2_params <- list(version = 1L, binsize = gam_binsize, sigma = gam_sigma,
                    ins_kb = gam2_ins_kb, short_max = gam2_short_max_bins,
                    long_min = gam2_long_min_bins, shells = gam2_shell_edges,
                    K = gam2_K, min_frac = gam2_min_frac, oe = "subtract")

# ---------------------------------------------------------------
# GAM_1 handoff: WT H3K4me3 per-bin signal, the orientation reference.
# ---------------------------------------------------------------
gam1_rds <- file.path(rds_dir, "gam_bins_annotated.rds")
if (!file.exists(gam1_rds))
  stop("GAM_1 output not found:\n  ", gam1_rds,
       "\nPC1 sign is arbitrary per chromosome and is oriented against WT H3K4me3 from",
       " GAM_1,\nso GAM_2 cannot run without it (no gene-density fallback, on purpose).",
       "\nRun it first:   cd GAM && Rscript GAM_1_Bins_and_Signal.R", call. = FALSE)
gam1 <- readRDS(gam1_rds)
stopifnot(all(c("signal", "stats") %in% names(gam1)))
h3k4_wt <- local({
  cols <- grep("^H3K4me3_WT", colnames(gam1$signal))   # sample_id = <mark>_<geno>_rep<n>
  if (length(cols)) {
    message("Orienting PC1 on ", length(cols), " WT H3K4me3 sample(s) from $signal")
    rowMeans(gam1$signal[, cols, drop = FALSE], na.rm = TRUE)
  } else {
    # Same mark, but pooled over genotypes: still an H3K4me3 reference, unlike
    # gene density. Loud, because WT-only is what was asked for.
    warning("no H3K4me3_WT columns in $signal; using $stats baseMean for H3K4me3",
            " (pooled over genotypes)", call. = FALSE)
    st <- gam1$stats[gam1$stats$mark == "H3K4me3", ]
    setNames(st$baseMean, st$bin_id)
  }
})
rm(gam1); invisible(gc())

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

# Leading eigenvector of the Pearson correlation matrix of the columns of `x`
# WITHOUT forming the n x n correlation matrix: C = Z'Z/(n-1) with Z the
# column-standardised x, so C v needs only products with x (x is symmetric, so
# x %*% u == t(x) %*% u). Mathematically identical to eigen(cor(x)); the
# explicit route would add another ~760 MB on chr1 and an O(n^3) eigen call.
pc1_corr <- function(x) {
  nv <- ncol(x)
  mu <- colMeans(x)
  s  <- vapply(seq_len(nv), function(j) stats::sd(x[, j]), 0)
  s[!is.finite(s) | s < 1e-12] <- Inf          # constant bin contributes nothing
  mv <- function(v, args = NULL) {
    w <- v / s
    u <- drop(x %*% w) - sum(mu * w)           # Z w
    drop(x %*% u - mu * sum(u)) / s / (nv - 1) # Z'(Z w) / (n-1)
  }
  if (requireNamespace("RSpectra", quietly = TRUE)) {
    e <- RSpectra::eigs_sym(mv, k = 2, n = nv, which = "LA")
    list(vec = e$vectors[, 1], val = e$values[1], gap = e$values[1] / e$values[2])
  } else {
    set.seed(1); v <- rnorm(nv); v <- v / sqrt(sum(v^2)); lam <- 0
    for (it in 1:2000) {
      w <- mv(v); lam_new <- sum(w * v); w <- w / sqrt(sum(w^2))
      if (abs(lam_new - lam) < 1e-10 * abs(lam_new) && it > 5) { v <- w; lam <- lam_new; break }
      v <- w; lam <- lam_new
    }
    list(vec = v, val = lam, gap = NA_real_)
  }
}

# PC1 sign: Spearman against a reference over bins finite in both. Rank-based so
# the few extreme H3K4me3 bins (promoter peaks) do not set the sign.
orient_sign <- function(pc1, ref) {
  ok <- is.finite(pc1) & is.finite(ref[names(pc1)])
  if (sum(ok) < 50) return(list(sign = NA_real_, rho = NA_real_, n = sum(ok)))
  rho <- suppressWarnings(stats::cor(pc1[ok], ref[names(pc1)][ok], method = "spearman"))
  list(sign = if (is.finite(rho) && rho < 0) -1 else 1, rho = rho, n = sum(ok))
}

# ---------------------------------------------------------------
# Reduce ONE matrix (one genotype, one chromosome) to per-bin features.
# The matrix is loaded here and dies here: nothing n x n escapes this function,
# so at most one matrix is ever alive. chr1 is ~9.7k^2 ~ 760 MB as doubles and
# this function peaks at ~2 copies (raw + O-E subset), so WT and HP1gKO are
# deliberately NOT held together - the caller loops genotypes sequentially and
# gc()s between them.
# ---------------------------------------------------------------
reduce_npmi <- function(genotype, chr, orient_ref, check_ref, can_orient) {
  m  <- load_npmi(genotype, chr)
  n  <- nrow(m); ids <- rownames(m)
  stopifnot(n == ncol(m), identical(ids, colnames(m)))
  start0 <- as.integer(sub("^.+:([0-9]+)-([0-9]+)$", "\\1", ids))
  stopifnot(!is.unsorted(start0))

  # Unmappable bins are all-NA rows/cols in the source. They must STAY NA all the
  # way through: a 0 would read as "never co-segregates" and drag PC1.
  nvalid <- colSums(!is.na(m))
  V <- nvalid > 0
  idx <- seq_len(n)

  # --- distance-banded NPMI + decay slope (raw NPMI, before O-E) -------------
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

  # --- observed - expected by genomic distance --------------------------------
  # Expected = mean over the d-th diagonal. SUBTRACTED, not divided: NPMI is
  # already a log-scale quantity (pointwise MI) and goes negative, so a ratio
  # would blow up or flip sign wherever a diagonal mean is near 0, and "O/E" of
  # a log-ratio is a difference anyway. Skipped alternative: LOESS on distance;
  # the diagonal mean is what Lieberman-Aiden used and needs no span choice.
  expd <- rep(NA_real_, n - 1L); expn <- integer(n - 1L)
  for (d in seq_len(n - 1L)) {
    r <- seq_len(n - d)
    up <- (r + d - 1) * n + r                  # linear index of (row r, col r+d)
    v <- m[up]; ok <- !is.na(v)
    expn[d] <- sum(ok)
    if (expn[d] > 0) { expd[d] <- mean(v[ok]); v <- v - expd[d] } else v[] <- NA
    m[up] <- v
    m[(r - 1) * n + r + d] <- v                # mirror to lower triangle
  }
  dd <- unique(pmin(round(10^seq(0, log10(n - 1), length.out = 80)), n - 1))
  curve <- data.table(chr = chr, genotype = genotype, dist_bp = dd * gam_binsize,
                      mean_npmi = expd[dd])

  # --- residual NA cells among MAPPABLE bins -----------------------------------
  # After dropping all-NA bins the only NA left is the diagonal (plus a stray
  # cell). Set them to 0 in O-E space = "exactly as expected for this distance",
  # i.e. neutral. That is NOT the unmappable-as-0 error: those bins are removed.
  # The stray cells are remembered so the saddle can exclude them from its means.
  na <- which(is.na(m))
  r <- (na - 1) %% n + 1; cc <- (na - 1) %/% n + 1
  stray <- r != cc & V[r] & V[cc]
  vmap <- cumsum(V)
  stray_r <- vmap[r[stray]]; stray_c <- vmap[cc[stray]]
  m[na] <- 0
  rm(na, r, cc, stray)
  oe <- if (all(V)) m else m[V, V, drop = FALSE]
  rm(m); invisible(gc())

  feat <- data.table(bin_id = ids, chr = chr, pc1 = NA_real_,
                     npmi_short = short, npmi_long = long, decay_slope = slope,
                     n_valid = as.integer(nvalid))
  for (w in names(ins)) set(feat, j = paste0("ins_", w), value = ins[[w]])

  # --- PC1 + orientation + saddle ----------------------------------------------
  sad <- NULL
  summ <- data.table(chr = chr, genotype = genotype, n_bins = n, n_mappable = sum(V),
                     oriented = FALSE, flip = NA, rho_orient = NA_real_,
                     rho_h3k4me3 = NA_real_, eig1 = NA_real_, eig_gap = NA_real_,
                     frac_A = NA_real_)
  if (can_orient && sum(V) > 2 * gam2_K) {
    pc <- pc1_corr(oe)
    p <- setNames(pc$vec, ids[V])
    # Median-centre BEFORE scaling. The eigenvector of a correlation matrix is not
    # mean-zero and on chr19 it sat ~2 SD above 0 (97% of bins "A" at a zero
    # cut), so sign(PC1) would not be an A/B call. Centring puts 0 at the middle
    # of the chromosome's own distribution, the same convention for both
    # genotypes; the cost is that a chromosome-wide A/B imbalance is not
    # represented (frac_A ~ 0.5 by construction) - fine for the WT-vs-KO question.
    p <- p - stats::median(p)
    p <- p / sqrt(mean(p^2))                        # RMS 1: comparable across chr/genotype
    # Sign by the reference for THIS genotype (WT: H3K4me3; HP1gKO: the WT PC1
    # of the same chromosome, so the delta is a real change and not a flip).
    o <- orient_sign(p, orient_ref)
    if (!is.finite(o$sign)) stop("cannot orient ", chr, " ", genotype,
                                 ": <50 bins overlap the reference")
    p <- p * o$sign
    if (!is.na(o$rho) && abs(o$rho) < 0.1)
      warning(chr, " ", genotype, ": |rho| with the orientation reference is only ",
              round(o$rho, 3), " - sign is weakly determined", call. = FALSE)
    ck <- orient_sign(p, check_ref)
    if (!is.na(ck$rho) && ck$rho < 0)
      warning(chr, " ", genotype, ": oriented PC1 anti-correlates with WT H3K4me3 (rho=",
              round(ck$rho, 3), ")", call. = FALSE)
    feat[match(names(p), bin_id), pc1 := p]
    summ[, `:=`(oriented = TRUE, flip = o$sign < 0, rho_orient = o$rho,
                rho_h3k4me3 = ck$rho, eig1 = pc$val, eig_gap = pc$gap,
                frac_A = mean(p > 0))]

    # Saddle: equal-count PC1 strata WITHIN the chromosome (rank-based) so it is a
    # single pass over the matrix. Genome-wide edges would need a second parse of
    # every matrix. S/C are kept as sums + counts so chromosomes pool exactly.
    nvv <- sum(V); K <- gam2_K
    s <- ceiling(rank(p, ties.method = "first") / nvv * K)
    Z <- matrix(0, nvv, K); Z[cbind(seq_len(nvv), s)] <- 1
    S <- crossprod(Z, oe %*% Z)
    nk <- tabulate(s, K)
    C <- outer(nk, nk) - diag(nk, K)                # remove the diagonal cells
    if (length(stray_r)) C <- C - matrix(tabulate((s[stray_r] - 1L) * K + s[stray_c], K * K), K)
    sad <- list(S = S, C = C)
  }
  list(feat = feat, curve = curve, saddle = sad, summary = summ)
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
# Cache holds ONLY the small per-bin reduction (+ a thinned decay curve and the
# KxK saddle sums): KB-MB per chromosome and genotype. The raw NPMI, O-E and
# correlation matrices are NEVER written to disk: summed over chromosomes and
# both genotypes they come to 13+ GB, which this machine's disk budget cannot
# take, and they are cheap to rebuild next to their size - only the one-matrix-
# at-a-time parse is slow, and that is exactly what this cache skips on a rerun.
# Do not "improve" this by caching the matrices.
# ---------------------------------------------------------------
for (chr in gam2_chroms) {
  # X/Y are not in GAM_1's autosome grid, so there is no H3K4me3 to orient them
  # by. An unoriented PC1 is meaningless for WT-vs-KO, so they get insulation and
  # contact bands but PC1 / compartment / saddle are left NA (not guessed).
  can_orient <- chr %in% gam_autosomes
  for (g in gam2_genos) {
    if (!is.null(read_cache(chr, g))) { message(chr, " ", g, ": cached"); next }
    message(chr, " ", g, ": parsing ", format(Sys.time(), "%H:%M:%S"))
    if (g == "WT") {
      ref <- h3k4_wt
    } else {
      wt <- read_cache(chr, "WT")
      ref <- if (!is.null(wt)) setNames(wt$feat$pc1, wt$feat$bin_id) else NULL
      stopifnot(!can_orient || !is.null(ref))
    }
    res <- reduce_npmi(g, chr, orient_ref = ref, check_ref = h3k4_wt, can_orient = can_orient)
    res$params <- gam2_params
    saveRDS(res, cache_file(chr, g))
    # Release this chromosome's matrices before the next parse.
    rm(res, ref); invisible(gc())
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
  d <- rbindlist(lapply(caches[[g]][done], `[[`, "feat"))
  d[, compartment := factor(fifelse(is.na(pc1), NA_character_, fifelse(pc1 > 0, "A", "B")),
                            levels = c("A", "B"))]
  d
})
bins <- per_geno$WT[, .(bin_id, chr)]
stopifnot(identical(bins$bin_id, per_geno$HP1gKO$bin_id))   # same grid in both genotypes
bins[, `:=`(start = as.integer(sub("^.+:([0-9]+)-([0-9]+)$", "\\1", bin_id)),
            end   = as.integer(sub("^.+:([0-9]+)-([0-9]+)$", "\\2", bin_id)))]

features <- copy(bins)
for (g in gam2_genos) {
  d <- per_geno[[g]]
  features[, (paste0(c("pc1", "compartment", "insulation", "npmi_short", "npmi_long",
                       "decay_slope", "n_valid"), "_", g)) :=
             d[, .(pc1, compartment, get(prim_col), npmi_short, npmi_long, decay_slope, n_valid)]]
}
features[, pc1_delta        := pc1_HP1gKO - pc1_WT]
features[, insulation_delta := insulation_HP1gKO - insulation_WT]
# Labelled WT_then_KO: A_B = A in WT, B in HP1gKO (a loss of A identity).
features[, compartment_switch := factor(
  fifelse(is.na(compartment_WT) | is.na(compartment_HP1gKO), NA_character_,
          paste0(compartment_WT, "_", compartment_HP1gKO)),
  levels = c("A_A", "A_B", "B_A", "B_B"))]
setkey(features, bin_id)

other_w <- setdiff(names(gam2_ins_kb), gam2_ins_primary)
insulation_windows <- setNames(lapply(other_w, function(w) {
  dt <- data.table(bin_id = per_geno$WT$bin_id)
  for (g in gam2_genos) dt[, (paste0("insulation_", g)) := per_geno[[g]][[paste0("ins_", w)]]]
  attr(dt, "window_bins") <- as.integer(round(gam2_ins_kb[[w]] * 1000 / gam_binsize))
  setkey(dt, bin_id)
}), other_w)

saddle <- lapply(setNames(gam2_genos, gam2_genos), function(g) {
  sd <- lapply(caches[[g]][done], `[[`, "saddle"); sd <- sd[!vapply(sd, is.null, TRUE)]
  if (!length(sd)) return(NULL)
  S <- Reduce(`+`, lapply(sd, `[[`, "S")); C <- Reduce(`+`, lapply(sd, `[[`, "C"))
  out <- S / C
  dimnames(out) <- list(paste0("Q", seq_len(gam2_K)), paste0("Q", seq_len(gam2_K)))
  out   # Q1 = lowest PC1 (B) ... QK = highest (A); mean O-E NPMI between strata
})

summ <- rbindlist(lapply(gam2_genos, function(g) rbindlist(lapply(caches[[g]][done], `[[`, "summary"))))
curves <- rbindlist(lapply(gam2_genos, function(g) rbindlist(lapply(caches[[g]][done], `[[`, "curve"))))
genome_curve <- curves[, .(mean_npmi = mean(mean_npmi, na.rm = TRUE)), by = .(genotype, dist_bp)]

meta <- list(n_pc1_bins = features[, sum(!is.na(pc1_WT) & !is.na(pc1_HP1gKO))],
             n_pc1_bins_by_genotype = c(WT = features[, sum(!is.na(pc1_WT))],
                                        HP1gKO = features[, sum(!is.na(pc1_HP1gKO))]),
             window_sizes_kb = gam2_ins_kb, primary_window = gam2_ins_primary,
             oe_method = "subtract diagonal mean (NPMI is log-scale); Pearson; PC1 via correlation-operator eigs",
             pc1_scale = "median-centred, RMS = 1 per chromosome x genotype; A = positive = H3K4me3-high",
             orientation = "WT: Spearman vs WT H3K4me3 ($signal); HP1gKO: vs WT PC1 of same chr",
             band_short_bins = c(1L, gam2_short_max_bins), band_long_min_bins = gam2_long_min_bins,
             decay_shell_edges_bins = gam2_shell_edges, saddle_K = gam2_K,
             saddle_strata = "per-chromosome equal-count PC1 quantiles, pooled by sums/counts",
             chroms_done = done, chroms_requested = gam2_chroms,
             unoriented_chroms = setdiff(done, gam_autosomes),
             decay_curve = genome_curve, orientation_summary = summ,
             binsize = gam_binsize, sigma = gam_sigma, date = Sys.time())

# Only a full-genome run may overwrite the consumer interface; a single-chromosome
# validation (GAM_2_CHROMS) writes alongside it so it can never masquerade as the real thing.
out_path <- if (setequal(gam2_chroms, gam_npmi_chroms)) gam2_out else
  sub("\\.rds$", paste0("_partial.rds"), gam2_out)
saveRDS(list(features = features, insulation_windows = insulation_windows,
             saddle = saddle, meta = meta), out_path)
message("Wrote ", out_path, " (", nrow(features), " bins, ", length(done), " chromosomes)")

fwrite(summ, file.path(tables_dir, "npmi_compartment_orientation.tsv"), sep = "\t")
fwrite(features[, .N, by = compartment_switch][order(compartment_switch)],
       file.path(tables_dir, "npmi_compartment_switch_counts.tsv"), sep = "\t")

# ---------------------------------------------------------------
# Figures
# ---------------------------------------------------------------
chr_levels <- gam_npmi_chroms[gam_npmi_chroms %in% done]
auto <- features[chr %in% gam_autosomes]
auto[, chr := factor(chr, levels = chr_levels)]

if (nrow(auto) && auto[, any(!is.na(pc1_WT))]) {
  # PC1 tracks, WT vs HP1gKO: filled area so the A (up) / B (down) read at a glance.
  trk <- melt(auto, id.vars = c("chr", "start"),
              measure.vars = c("pc1_WT", "pc1_HP1gKO"),
              variable.name = "genotype", value.name = "pc1")
  trk[, genotype := factor(sub("pc1_", "", genotype), levels = gam2_genos)]
  p_trk <- ggplot(trk[!is.na(pc1)], aes(start / 1e6, pc1, colour = genotype, fill = genotype)) +
    geom_area(alpha = 0.35, linewidth = 0.2, position = "identity") +
    facet_wrap(~chr, scales = "free_x", ncol = 3) +
    scale_colour_manual(values = geno_colors) + scale_fill_manual(values = geno_colors) +
    labs(x = "Position (Mb)", y = "PC1 (A > 0, H3K4me3-oriented)", colour = NULL, fill = NULL) +
    theme_m + theme(legend.position = "top")
  save_fig(p_trk, "pc1_tracks_by_chr", subdir = "compartments",
           width = 14, height = 2 + 1.6 * ceiling(length(unique(trk$chr)) / 3))

  # WT vs KO PC1: off-diagonal points that change sign are the switching bins.
  sc <- auto[!is.na(pc1_WT) & !is.na(pc1_HP1gKO)]
  r <- cor(sc$pc1_WT, sc$pc1_HP1gKO)
  p_sc <- ggplot(sc, aes(pc1_WT, pc1_HP1gKO, colour = compartment_switch)) +
    geom_hline(yintercept = 0, colour = "grey75") + geom_vline(xintercept = 0, colour = "grey75") +
    geom_point(size = 0.3, alpha = 0.35) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    scale_colour_disc(name = "WT -> HP1gKO") +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, label = sprintf("r = %.3f", r)) +
    labs(x = "PC1 WT", y = "PC1 HP1gKO") + theme_m + coord_equal()
  save_fig(p_sc, "pc1_scatter_WT_vs_HP1gKO", subdir = "compartments", width = 7, height = 6.5)
}

sad_long <- rbindlist(lapply(gam2_genos, function(g) {
  s <- saddle[[g]]; if (is.null(s)) return(NULL)
  data.table(genotype = factor(g, levels = gam2_genos),
             row = rep(seq_len(nrow(s)), ncol(s)), col = rep(seq_len(ncol(s)), each = nrow(s)),
             value = as.vector(s))
}))
if (nrow(sad_long)) {
  lim <- quantile(abs(sad_long$value), 0.98, na.rm = TRUE)    # clip so one hot corner does not wash the scale
  p_sad <- ggplot(sad_long, aes(col, row, fill = value)) + geom_raster() +
    facet_wrap(~genotype) + coord_equal() +
    scale_fill_heat0_div(limits = c(-lim, lim), oob = scales::squish, name = "mean O-E\nNPMI") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    labs(x = "PC1 stratum (B -> A)", y = "PC1 stratum (B -> A)") + theme_m
  save_fig(p_sad, "saddle_by_genotype", subdir = "compartments", width = 10, height = 5)
}

if (nrow(curves)) {
  curves[, genotype := factor(genotype, levels = gam2_genos)]
  genome_curve[, genotype := factor(genotype, levels = gam2_genos)]
  p_dec <- ggplot() +
    geom_line(data = curves, aes(dist_bp / 1e6, mean_npmi, group = chr), colour = "grey75", linewidth = 0.3) +
    geom_line(data = genome_curve, aes(dist_bp / 1e6, mean_npmi, colour = genotype), linewidth = 1) +
    facet_wrap(~genotype) + scale_x_log10() +
    scale_colour_manual(values = geno_colors, guide = "none") +
    labs(x = "Genomic distance (Mb)", y = "Mean NPMI (diagonal mean)") + theme_m
  save_fig(p_dec, "contact_decay_by_genotype", subdir = "compartments", width = 10, height = 4.5)
}
