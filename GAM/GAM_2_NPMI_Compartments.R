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
# version is bumped whenever the COMPARTMENT method changes, not just its
# parameters - a cached reduction carries the PC it was computed with, so a
# change to eigenvector selection must invalidate it too. v2: compartments
# called at gam_compartment_binsize with balance-first selection over 4
# eigenvectors (see select_compartment_pc).
gam2_params <- list(version = 2L, binsize = gam_binsize, sigma = gam_sigma,
                    ins_kb = gam2_ins_kb, short_max = gam2_short_max_bins,
                    long_min = gam2_long_min_bins, shells = gam2_shell_edges,
                    K = gam2_K, min_frac = gam2_min_frac, oe = "subtract",
                    comp_binsize = gam_compartment_binsize,
                    comp_select = "balance-first, k=4, frac_A in [0.2,0.8]")

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
# ACTIVE-CHROMATIN REFERENCE: mean of the per-bin z-scores of WT H3K4me3 and WT
# H3K36me3. This both SELECTS which eigenvector is the compartment axis and
# orients its sign (see select_compartment_pc() for the full rationale).
#
# Two marks, not one: H3K4me3 is punctate at promoters and thin across a 20 kb
# bin; H3K36me3 is broad over transcribed gene bodies and better matched to this
# resolution. Averaging their z-scores gives a steadier A-compartment anchor.
#
# Both are UNCHANGED in the KO (per-bin median log2FC +0.000 and -0.052, against
# -0.497 for H4K20me3), so the A/B axis is the same ruler in both genotypes.
# Anchoring on H3K9me3/H4K20me3 would be circular - the axis would move with the
# effect under study.
h3k4_wt <- local({
  zs <- function(v) { s <- stats::sd(v, na.rm = TRUE)
                      if (!is.finite(s) || s < 1e-12) return(rep(0, length(v)))
                      (v - mean(v, na.rm = TRUE)) / s }
  parts <- list(); used <- character(0)
  for (mk in c("H3K4me3", "H3K36me3")) {
    cols <- grep(paste0("^", mk, "_WT"), colnames(gam1$signal))
    if (length(cols)) {
      parts[[mk]] <- zs(rowMeans(gam1$signal[, cols, drop = FALSE], na.rm = TRUE))
      used <- c(used, paste0(mk, " (", length(cols), " WT libs)"))
    }
  }
  if (!length(parts)) {
    warning("no WT H3K4me3/H3K36me3 columns in $signal; falling back to pooled",
            " H3K4me3 baseMean - the A/B call will be weaker", call. = FALSE)
    st <- gam1$stats[gam1$stats$mark == "H3K4me3", ]
    return(setNames(st$baseMean, st$bin_id))
  }
  message("Active-chromatin reference from: ", paste(used, collapse = " + "))
  setNames(Reduce(`+`, parts) / length(parts), rownames(gam1$signal))
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
    e <- RSpectra::eigs_sym(mv, k = 4, n = nv, which = "LA")
    # BOTH eigenvectors are returned. Which one carries the A/B compartment is
    # decided downstream by correlation with an active-chromatin reference, NOT
    # by eigenvalue rank: in Hi-C-like data the leading eigenvector frequently
    # captures overall contact structure rather than compartments (on chr19 it
    # put 97% of bins on one side of zero, which is not an A/B split), and the
    # compartment signal then sits in the second. Selecting by marker is what
    # cooltools-style pipelines do with GC/gene density; here we have the real
    # marks, which is a stronger anchor.
    list(vecs = e$vectors[, 1:4, drop = FALSE], vals = e$values[1:4],
         vec = e$vectors[, 1], val = e$values[1], gap = e$values[1] / e$values[2])
  } else {
    set.seed(1); v <- rnorm(nv); v <- v / sqrt(sum(v^2)); lam <- 0
    for (it in 1:2000) {
      w <- mv(v); lam_new <- sum(w * v); w <- w / sqrt(sum(w^2))
      if (abs(lam_new - lam) < 1e-10 * abs(lam_new) && it > 5) { v <- w; lam <- lam_new; break }
      v <- w; lam <- lam_new
    }
    list(vecs = matrix(v, ncol = 1), vals = lam, vec = v, val = lam, gap = NA_real_)
  }
}

# Choose WHICH eigenvector is the compartment axis, and orient it.
#
# Reference = active chromatin, built from WT H3K4me3 AND H3K36me3 (mean of their
# per-bin z-scores). Two marks rather than one because they mark active chromatin
# differently and complementarily: H3K4me3 is punctate at promoters, which is a
# thin signal across a 20 kb bin, while H3K36me3 is broad over transcribed gene
# bodies and therefore better matched to this bin size. Their average is a more
# stable A-compartment anchor than either alone.
#
# Critically, these are the two marks that DO NOT change in the KO (per-bin median
# log2FC +0.000 and -0.052, against -0.497 for H4K20me3). Anchoring on them means
# the A/B axis is the same ruler in both genotypes. Anchoring on H3K9me3/H4K20me3
# would be circular: the axis would move with the very effect being measured.
#
# Returns the selected vector already sign-oriented so positive = A = active-high,
# plus the diagnostics needed to judge the call per chromosome.
select_compartment_pc <- function(pc, ids, ref, balance = c(0.2, 0.8)) {
  V <- pc$vecs
  rhos <- apply(V, 2, function(v) {
    names(v) <- ids
    ok <- is.finite(v) & is.finite(ref[ids])
    if (sum(ok) < 50) return(NA_real_)
    suppressWarnings(stats::cor(v[ok], ref[ids][ok], method = "spearman"))
  })
  # Fraction on the A side once oriented - the balance of each candidate.
  fracs <- vapply(seq_len(ncol(V)), function(k) {
    v <- V[, k]; sgn <- if (is.finite(rhos[k]) && rhos[k] < 0) -1 else 1
    mean(v * sgn > 0, na.rm = TRUE)
  }, 0)
  # SELECT ON BALANCE FIRST, then correlation. Ranking on |rho| alone picks the
  # wrong vector: measured across chromosomes, the LEADING eigenvector is
  # systematically degenerate (frac_A 0.04-0.13 - almost every bin on one side,
  # i.e. the overall contact-structure component) yet often correlates slightly
  # better with the active reference than the true compartment vector does. On
  # chr19 WT, PC1 had |rho| 0.32 at frac_A 0.038 while PC2 had |rho| 0.23 at
  # frac_A 0.64; on chr11 WT, PC2 reached |rho| 0.53 at frac_A 0.68. A vector
  # that assigns 96% of a chromosome to one compartment is not an A/B split
  # whatever it correlates with, so balance is a precondition, not a tiebreak.
  cand <- which(is.finite(rhos) & fracs > balance[1] & fracs < balance[2])
  if (!length(cand)) {
    # No balanced candidate: the axis is not determined for this chromosome.
    # Returning NA is the honest outcome - better than shipping a degenerate
    # vector that downstream code would treat as a compartment call.
    return(list(vec = setNames(rep(NA_real_, length(ids)), ids), which_pc = NA_integer_,
                rho = NA_real_, rho_all = rhos, frac_all = fracs, sign = NA_real_,
                frac_A = NA_real_, determined = FALSE))
  }
  k <- cand[which.max(abs(rhos[cand]))]
  v <- V[, k]; names(v) <- ids
  sgn <- if (is.finite(rhos[k]) && rhos[k] < 0) -1 else 1
  list(vec = v * sgn, which_pc = k, rho = abs(rhos[k]),
       rho_all = rhos, frac_all = fracs, sign = sgn, determined = TRUE,
       # A/B from the SIGN of the selected eigenvector, with no median-centring.
       # Centring was previously needed only because the wrong (all-positive)
       # eigenvector was being used; selecting the right one removes the need and
       # frees frac_A to differ from 0.5, which is what makes a net A<->B shift
       # detectable at all.
       frac_A = mean(v * sgn > 0, na.rm = TRUE))
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
                     frac_A = NA_real_, which_pc = NA_integer_,
                     rho_pc1 = NA_real_, rho_pc2 = NA_real_)
  if (can_orient && sum(V) > 2 * gam2_K) {
    # --- aggregate to the compartment resolution (see gam_compartment_binsize) --
    # Block-mean the O-E matrix over f x f blocks of native bins. Compartments
    # are a long-range property and are not resolvable at 20 kb; this asks the
    # eigendecomposition for what the data can actually support. NA-aware, so a
    # coarse bin built from partially unmappable native bins still contributes.
    f <- gam_compartment_factor
    nv_idx <- which(V)                      # native indices behind oe's rows
    grp    <- ((seq_along(nv_idx) - 1L) %/% f) + 1L   # coarse bin per oe row
    ngrp   <- max(grp)
    if (ngrp < 10L) stop(chr, " ", genotype, ": only ", ngrp,
                         " coarse bins - too few to call compartments")
    A <- matrix(0, nrow = ngrp, ncol = ncol(oe))
    cnt <- tabulate(grp, ngrp)
    oe0 <- oe; nna <- !is.na(oe0); oe0[!nna] <- 0
    A  <- rowsum(oe0, grp, reorder = FALSE)            # sums down rows
    An <- rowsum(nna * 1,  grp, reorder = FALSE)
    B  <- t(rowsum(t(A),  grp, reorder = FALSE))       # then across columns
    Bn <- t(rowsum(t(An), grp, reorder = FALSE))
    oe_c <- B / pmax(Bn, 1); oe_c[Bn == 0] <- NA
    rm(A, An, B, Bn, oe0, nna); invisible(gc())
    # keep only coarse bins with enough support
    keep_c <- rowSums(is.finite(oe_c)) > 0.5 * ncol(oe_c)
    if (sum(keep_c) < 10L) stop(chr, " ", genotype, ": too few usable coarse bins")
    oe_c <- oe_c[keep_c, keep_c, drop = FALSE]
    oe_c[!is.finite(oe_c)] <- 0
    # ids for the coarse bins = first native bin of each block; and a map from
    # each native bin back to its coarse bin (NA where the block was dropped)
    first_of <- nv_idx[match(seq_len(ngrp), grp)]
    agg_ids  <- ids[first_of][keep_c]
    cmap <- rep(NA_integer_, ngrp); cmap[keep_c] <- seq_len(sum(keep_c))
    agg_of_bin <- cmap[grp]

    pc <- pc1_corr(oe_c)
    sel <- select_compartment_pc(pc, agg_ids, orient_ref)
    # Expand the coarse-bin PC back onto the native 20 kb bins: every 20 kb bin
    # inherits the value of the 100 kb bin containing it. The compartment is a
    # property of the larger region, so this is an expansion of a coarser call,
    # not an interpolation - deliberately a step function.
    # An undetermined axis is a legitimate outcome, not an error: keep the
    # chromosome's insulation, distance bands and decay slope (none of which
    # depend on the PC) and leave pc1/compartment NA rather than shipping a
    # degenerate vector that downstream code would read as a compartment call.
    if (!isTRUE(sel$determined)) {
      warning(chr, " ", genotype, ": no balanced compartment eigenvector among the",
              " top 4 (frac_A ", paste(sprintf("%.2f", sel$frac_all), collapse = "/"),
              ") - pc1/compartment left NA", call. = FALSE)
      summ[, `:=`(which_pc = NA_integer_, rho_pc1 = sel$rho_all[1],
                  rho_pc2 = if (length(sel$rho_all) > 1) sel$rho_all[2] else NA_real_)]
    } else {
      p <- sel$vec[agg_of_bin]
      names(p) <- ids[V]
      p <- p[is.finite(p)]
      # NO median-centring. It was only ever needed because the all-positive
      # leading eigenvector was being used unconditionally; selecting the
      # compartment eigenvector by the active-mark reference removes the need, and
      # removing it is what frees frac_A to differ from 0.5 - without which a net
      # A<->B shift cannot be detected at all (the switch table would be forced
      # near-symmetric whatever the biology did).
      p <- p / sqrt(mean(p^2))                        # RMS 1: comparable across chr/genotype
      # Sign is already set by the selector against the active reference. For the
      # KO, re-orient against the WT PC of the same chromosome so the delta is a
      # real change rather than a flip.
      o <- orient_sign(p, orient_ref)   # p is already back on native bins
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
                  frac_A = mean(p > 0),
                  # which eigenvector carried the compartment signal, and how
                  # clearly - if which_pc is often 2, the leading eigenvector is
                  # structure rather than compartments and that belongs in the
                  # methods, not buried.
                  which_pc = sel$which_pc, rho_pc1 = sel$rho_all[1],
                  rho_pc2 = if (length(sel$rho_all) > 1) sel$rho_all[2] else NA_real_)]

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
