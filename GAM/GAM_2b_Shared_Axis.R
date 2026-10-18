# ================================================================
# GAM_2b - ONE SHARED A/B AXIS, with both genotypes scored on it
# ----------------------------------------------------------------
#   cd GAM && Rscript GAM_2b_Shared_Axis.R
#
# Needs: GAM_1 (the active-chromatin reference) and the NPMI matrices.
# Does NOT need GAM_2; it recomputes O-E from the matrices itself.
#
# WHY THIS EXISTS
# ---------------
# GAM_2 selects a compartment eigenvector INDEPENDENTLY IN EACH GENOTYPE and
# reports pc1_WT, pc1_HP1gKO and their difference. That difference is not
# interpretable, and the reason is recorded in CLAUDE.md (settled point 2): the
# two genotypes do not select the same eigenvector. Mean squared cosine between
# their leading 4-eigenvector subspaces is 0.27-0.57, the 4th principal angle is
# near-orthogonal, and the selected PC differs between genotypes on 4 of the 5
# chromosomes tested. pc1_delta therefore subtracts two different objects.
#
# The fix, also already established (settled point 3), is to stop asking each
# genotype for its own axis. Build ONE axis from the averaged correlation
# matrix, phase it once, and score both genotypes against that single ruler.
# Held-out validation showed this generalises rather than overfitting: weights
# fitted on odd coarse bins and evaluated on even gave chr16 0.555 -> 0.565,
# chr11 0.519 -> 0.553, chr1 0.435 -> 0.437, and k = 4 and k = 8 agreed.
#
# WHAT THE SCORE IS
# -----------------
# With `a` the shared axis and C_g a genotype's correlation matrix of the
# distance-corrected contact map, the per-bin score is
#
#     s_g = C_g %*% a
#
# i.e. how strongly bin i's correlation profile in genotype g aligns with the
# shared compartment pattern. This is one power-iteration step, so it degenerates
# correctly: if `a` happened to be C_g's own eigenvector, C_g a = lambda a and the
# score is that eigenvector back again. Unlike a per-genotype eigenvector it is
# defined by the SAME vector in both genotypes, so the difference of two scores is
# a difference of one measured quantity rather than of two different objects.
#
# AMPLITUDE IS DELIBERATELY REMOVED FROM THE TRACK
# ------------------------------------------------
# Each genotype's score is median-centred and scaled to RMS 1 per chromosome, so
# the track carries SHAPE only. This is not a cosmetic choice. lambda1 falls
# about 20% in the cKO (ratio 0.58-1.02, median 0.81), but that is a global
# property of the correlation matrix affecting every direction equally, not a
# compartment-specific effect - settled point 6, and normalising by each
# genotype's own lambda1 leaves a strength ratio of 0.953, i.e. nothing. If the
# amplitude were left in the track, that global drop would paint the cKO as
# uniformly weaker everywhere and reproduce exactly the claim that was retracted.
# The amplitude diagnostics are kept, per chromosome, in the summary table where
# they can be read as what they are.
#
# This stage does NOT re-answer the strength question. That was settled on axes
# external to both contact matrices (active chromatin p = 0.62, TSS density
# p = 0.74, both null across all 19 autosomes) and is not revisited here.
#
# Outputs
#   results/rds/shared_axis_scores.rds        per 20 kb bin, both genotypes
#   results/tables/gam2b_shared_axis_summary.tsv   per-chromosome diagnostics
#   results/igv/gam_shared_axis.bedgraph           the axis itself
#   results/igv/gam_shared_score_{WT,HP1gKO,delta}.bedgraph
#   results/figures/compartments/shared_axis_*.{png,pdf}
# ================================================================
source("config.R")

msg <- function(fmt, ...) message(sprintf(fmt, ...))

gam2b_params <- list(
  version      = 1L,
  binsize      = gam_compartment_binsize,   # 100 kb, see config.R
  k            = 4L,                        # eigenvectors offered to the selector
  min_support  = 0.5,                       # coarse bin kept if >50% of its row is finite
  chunk        = 1000L                      # columns per block-mean chunk
)

fig_sub    <- "compartments"
out_rds    <- file.path(rds_dir, "shared_axis_scores.rds")
igv_dir    <- file.path(results_dir, "igv")
dir.create(igv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# The ruler. Same function GAM_2 phases its per-genotype axes with, so the shared
# axis and the per-genotype axes point the same way and can be compared.
ref20 <- active_chromatin_ref()

# ---------------------------------------------------------------
# Block-mean a matrix down to coarse bins, NA-aware, in column chunks.
#
# Chunked because the obvious one-liner is what makes this stage run the machine
# out of memory: rowsum() needs a numeric matrix of counts the same shape as the
# input, and on chr1 the input is already ~762 MB, so materialising that second
# full matrix doubles the peak. Chunking bounds the extra to one chunk's worth.
# ---------------------------------------------------------------
block_mean <- function(m, grp, ngrp, chunk = gam2b_params$chunk) {
  n <- ncol(m)
  S <- matrix(0, ngrp, n); K <- matrix(0, ngrp, n)
  for (j0 in seq(1L, n, by = chunk)) {
    j  <- j0:min(j0 + chunk - 1L, n)
    b  <- m[, j, drop = FALSE]
    ok <- !is.na(b); b[!ok] <- 0
    S[, j] <- rowsum(b,      grp, reorder = FALSE)
    K[, j] <- rowsum(ok * 1, grp, reorder = FALSE)
  }
  S2 <- t(rowsum(t(S), grp, reorder = FALSE))
  K2 <- t(rowsum(t(K), grp, reorder = FALSE))
  rm(S, K); invisible(gc())
  out <- S2 / pmax(K2, 1)
  out[K2 == 0] <- NA_real_
  out
}

# ---------------------------------------------------------------
# One genotype, one chromosome -> the coarse observed-minus-expected matrix.
# The 20 kb matrix is loaded here and dies here.
# ---------------------------------------------------------------
coarse_oe <- function(genotype, chr) {
  m <- load_npmi(genotype, chr)
  n <- nrow(m); ids <- rownames(m)
  stopifnot(n == ncol(m), identical(ids, colnames(m)))
  start0 <- as.integer(sub("^.+:([0-9]+)-([0-9]+)$", "\\1", ids))
  stopifnot(!is.unsorted(start0))

  # O-E by SUBTRACTION of each diagonal's mean, not division. NPMI is already a
  # log-scale quantity (pointwise mutual information) and goes negative, so a
  # ratio blows up or flips sign wherever a diagonal mean is near zero. Same
  # convention as GAM_2; see CLAUDE.md settled point 8 for why this also means
  # our saddle statistic is additive and not comparable to published values.
  for (d in seq_len(n - 1L)) {
    r  <- seq_len(n - d)
    up <- (r + d - 1) * n + r
    v  <- m[up]; ok <- !is.na(v)
    if (any(ok)) v <- v - mean(v[ok]) else v[] <- NA_real_
    m[up] <- v
    m[(r - 1) * n + r + d] <- v
  }

  # GROUP BY GENOMIC COORDINATE, NOT BY RANK WITHIN THE MAPPABLE SET.
  # This is the single most important line in the file. Grouping by rank - which
  # is what you get from seq_along() over the mappable bins - gives WT and cKO
  # DIFFERENT coarse grids whenever their unmappable bins differ, so "the same"
  # coarse bin is a different stretch of chromosome in each genotype. An earlier
  # version of this analysis did exactly that and 2,888 coarse bins collapsed to
  # 67 shared ones instead of 578. Coordinates are the only thing both genotypes
  # agree on by construction.
  grp  <- as.integer(start0 %/% gam2b_params$binsize) + 1L
  ugrp <- sort(unique(grp))
  gidx <- match(grp, ugrp)
  ngrp <- length(ugrp)
  if (ngrp < 10L) stop(chr, " ", genotype, ": only ", ngrp, " coarse bins")

  oe_c <- block_mean(m, gidx, ngrp)
  rm(m); invisible(gc())

  keep <- rowSums(is.finite(oe_c)) > gam2b_params$min_support * ncol(oe_c)
  cstart <- (ugrp - 1L) * gam2b_params$binsize
  list(oe     = oe_c,
       ids    = sprintf("%s:%d-%d", chr, cstart, cstart + gam2b_params$binsize),
       cstart = cstart,
       keep   = keep)
}

# Correlation matrix of the columns, NA-free. Explicit here (unlike GAM_2, which
# avoids forming it) because at 100 kb it is a few thousand square and the whole
# point of this stage is to AVERAGE the two genotypes' correlation matrices,
# which needs both of them in hand.
cor_mat <- function(oe) {
  oe[!is.finite(oe)] <- 0
  C <- suppressWarnings(stats::cor(oe))
  C[!is.finite(C)] <- 0
  C
}

top_eigs <- function(C, k) {
  if (requireNamespace("RSpectra", quietly = TRUE)) {
    e <- RSpectra::eigs_sym(C, k = k, which = "LA")
    list(vecs = e$vectors, vals = e$values)
  } else {
    e <- eigen(C, symmetric = TRUE)
    list(vecs = e$vectors[, seq_len(k), drop = FALSE], vals = e$values[seq_len(k)])
  }
}

# ---------------------------------------------------------------
# Per chromosome
# ---------------------------------------------------------------
msg("Shared A/B axis at %d kb, k = %d eigenvectors",
    gam2b_params$binsize / 1000, gam2b_params$k)

scores <- list(); summ <- list()

# Spot-run one or a few chromosomes without the whole genome:
#   GAM2B_CHROMS=chr19,chr16 Rscript GAM_2b_Shared_Axis.R
# Each chromosome loads two NPMI matrices, so a full pass is heavy; this is for
# checking a change to the method, not for producing the outputs.
CHROMS <- local({
  e <- Sys.getenv("GAM2B_CHROMS", "")
  if (!nzchar(e)) return(gam_autosomes)
  w <- trimws(strsplit(e, ",")[[1]])
  bad <- setdiff(w, gam_autosomes)
  if (length(bad)) stop("GAM2B_CHROMS: not autosomes: ", paste(bad, collapse = ", "))
  msg("GAM2B_CHROMS set: %d of %d chromosomes - PARTIAL RUN, outputs are incomplete",
      length(w), length(gam_autosomes))
  w
})
PARTIAL <- !identical(CHROMS, gam_autosomes)

for (chr in CHROMS) {
  msg("  %s", chr)
  W <- coarse_oe("WT", chr)
  K <- coarse_oe("HP1gKO", chr)
  stopifnot(identical(W$ids, K$ids))      # guaranteed by coordinate grouping

  common <- which(W$keep & K$keep)
  if (length(common) < 10L) {
    warning(chr, ": only ", length(common), " coarse bins usable in both genotypes",
            " - skipped", call. = FALSE)
    next
  }
  ids <- W$ids[common]

  C_wt <- cor_mat(W$oe[common, common, drop = FALSE])
  C_ko <- cor_mat(K$oe[common, common, drop = FALSE])
  rm(W, K); invisible(gc())

  # THE SHARED AXIS: eigenvectors of the AVERAGED correlation matrix. Averaging
  # the correlation matrices rather than the contact maps keeps each genotype's
  # own distance correction intact.
  C_av <- (C_wt + C_ko) / 2
  e_av <- top_eigs(C_av, gam2b_params$k)

  # The active-chromatin reference, coarsened to the same grid by averaging the
  # 20 kb bins inside each 100 kb bin.
  ref_c <- local({
    nb <- as.integer(gam2b_params$binsize / gam_binsize)
    cs <- as.integer(sub("^.+:([0-9]+)-.*$", "\\1", ids))
    vapply(cs, function(s) {
      b <- sprintf("%s:%d-%d", chr, s + (0:(nb - 1L)) * gam_binsize,
                   s + (1:nb) * gam_binsize)
      v <- ref20[b]
      if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
    }, 0)
  })
  names(ref_c) <- ids

  sel <- select_compartment_pc(list(vecs = e_av$vecs), ids, ref_c)

  e_wt <- top_eigs(C_wt, 2L); e_ko <- top_eigs(C_ko, 2L)
  row <- data.table(
    chr = chr, n_coarse = length(ids),
    determined = isTRUE(sel$determined),
    which_pc = sel$which_pc, rho_ref = sel$rho, frac_A = sel$frac_A,
    lambda1_WT = e_wt$vals[1], lambda1_HP1gKO = e_ko$vals[1],
    lambda1_ratio = e_ko$vals[1] / e_wt$vals[1])

  if (!isTRUE(sel$determined)) {
    warning(chr, ": no balanced eigenvector on the averaged matrix among the top ",
            gam2b_params$k, " - no shared axis for this chromosome", call. = FALSE)
    row[, `:=`(rho_score_WT = NA_real_, rho_score_HP1gKO = NA_real_,
               rho_WT_vs_HP1gKO = NA_real_,
               rayleigh_WT = NA_real_, rayleigh_HP1gKO = NA_real_)]
    summ[[chr]] <- row
    rm(C_wt, C_ko, C_av); invisible(gc())
    next
  }

  a <- sel$vec                                   # oriented: positive = A = active
  s_wt <- drop(C_wt %*% a)
  s_ko <- drop(C_ko %*% a)

  # Amplitude along the shared axis, BEFORE normalisation, kept as a diagnostic
  # only. a'Ca / (a'a) is the Rayleigh quotient: how much of each genotype's
  # correlation structure lies along this direction.
  rq <- function(C) drop(crossprod(a, C %*% a)) / drop(crossprod(a))
  row[, `:=`(rayleigh_WT = rq(C_wt), rayleigh_HP1gKO = rq(C_ko))]

  # Shape only, per chromosome per genotype - see the header.
  nrm <- function(v) { v <- v - stats::median(v, na.rm = TRUE)
                       s <- sqrt(mean(v^2, na.rm = TRUE))
                       if (!is.finite(s) || s < 1e-12) v else v / s }
  s_wt_n <- nrm(s_wt); s_ko_n <- nrm(s_ko)

  row[, `:=`(rho_score_WT      = stats::cor(s_wt_n, ref_c, method = "spearman",
                                            use = "pairwise.complete.obs"),
             rho_score_HP1gKO  = stats::cor(s_ko_n, ref_c, method = "spearman",
                                            use = "pairwise.complete.obs"),
             rho_WT_vs_HP1gKO  = stats::cor(s_wt_n, s_ko_n, method = "spearman"))]
  summ[[chr]] <- row

  scores[[chr]] <- data.table(
    chr = chr, cstart = as.integer(sub("^.+:([0-9]+)-.*$", "\\1", ids)),
    shared_axis = as.numeric(a / sqrt(mean(a^2))),
    score_WT = s_wt_n, score_HP1gKO = s_ko_n,
    score_delta = s_ko_n - s_wt_n)

  rm(C_wt, C_ko, C_av); invisible(gc())
}

SUMM <- rbindlist(summ, fill = TRUE)
SC   <- rbindlist(scores, fill = TRUE)
stopifnot(nrow(SC) > 0)

# DOES IT ACTUALLY BEAT THE PER-GENOTYPE AXES? That is the whole claim of this
# stage, so it is measured here rather than asserted. GAM_2's own per-chromosome
# diagnostics carry |rho| against the same active-chromatin reference, so the two
# are directly comparable. Absent (GAM_2 not run) the columns are simply NA.
og <- file.path(tables_dir, "npmi_compartment_orientation.tsv")
if (file.exists(og)) {
  O <- fread(og)[chr %in% gam_autosomes]
  # SIGNED, not abs(). The sign IS the diagnostic: a NEGATIVE rho means that
  # genotype's A and B are swapped on that chromosome. abs() here erased exactly
  # two such cases - chr15 (-0.061) and chr19 (-0.134), both cKO - so the table
  # that exists to catch the inversion was hiding it.
  W <- dcast(O[, .(chr, genotype, rho = rho_h3k4me3)], chr ~ genotype,
             value.var = "rho")
  setnames(W, c("WT", "HP1gKO"), c("rho_pergeno_WT", "rho_pergeno_HP1gKO"),
           skip_absent = TRUE)
  P <- dcast(O[, .(chr, genotype, which_pc)], chr ~ genotype, value.var = "which_pc")
  setnames(P, c("WT", "HP1gKO"), c("pc_pergeno_WT", "pc_pergeno_HP1gKO"),
           skip_absent = TRUE)
  SUMM <- merge(SUMM, merge(W, P, by = "chr"), by = "chr", all.x = TRUE)
  SUMM[, beats_both := is.finite(rho_ref) &
         rho_ref > pmax(rho_pergeno_WT, rho_pergeno_HP1gKO, na.rm = TRUE)]
  SUMM[, same_pc_both_genotypes := pc_pergeno_WT == pc_pergeno_HP1gKO]
  setorder(SUMM, chr)
}

fwrite(SUMM, file.path(tables_dir, "gam2b_shared_axis_summary.tsv"), sep = "\t")
msg("  axis determined on %d of %d autosomes; median |rho| vs active chromatin %.3f",
    SUMM[, sum(determined)], nrow(SUMM), SUMM[determined == TRUE, median(rho_ref)])
msg("  usable (|rho| > 0.4): %d of %d", SUMM[determined == TRUE & rho_ref > 0.4, .N],
    nrow(SUMM))
if ("beats_both" %in% names(SUMM)) {
  msg("  beats BOTH per-genotype axes on %d of %d chromosomes",
      SUMM[beats_both == TRUE, .N], SUMM[is.finite(rho_pergeno_WT), .N])
  msg("  per-genotype median |rho| %.3f (WT) / %.3f (cKO) vs shared %.3f",
      SUMM[, median(rho_pergeno_WT, na.rm = TRUE)],
      SUMM[, median(rho_pergeno_HP1gKO, na.rm = TRUE)],
      SUMM[determined == TRUE, median(rho_ref)])
  msg("  the two genotypes picked the SAME eigenvector on %d of %d chromosomes",
      SUMM[same_pc_both_genotypes == TRUE, .N],
      SUMM[is.finite(pc_pergeno_WT), .N])
}

# ---------------------------------------------------------------
# Expand the 100 kb scores onto the 20 kb grid, the way GAM_2 expands its PC.
# A step function on purpose: the compartment is a property of the larger
# region, so this is the expansion of a coarser call, not an interpolation.
# ---------------------------------------------------------------
nb <- as.integer(gam2b_params$binsize / gam_binsize)
FINE <- SC[rep(seq_len(.N), each = nb)]
FINE[, start := cstart + (seq_len(.N) - 1L) %% nb * gam_binsize, by = .(chr, cstart)]
FINE[, end := start + gam_binsize]
FINE[, bin_id := sprintf("%s:%d-%d", chr, start, end)]
FINE[, cstart := NULL]
setcolorder(FINE, c("bin_id", "chr", "start", "end"))

# Only keep bins that exist on the 20 kb analysis grid, and say how many of the
# grid's bins got a score. A silent drop here would look like a clean result.
grid <- gam_bin_id(gam_bins())
FINE <- FINE[bin_id %in% grid]
cov_frac <- nrow(FINE) / length(grid)
msg("  scored %d of %d 20 kb autosomal bins (%.1f%%)", nrow(FINE), length(grid),
    100 * cov_frac)
# Only a full pass can be checked against the whole grid. On a partial run this
# would fail by construction, and silently lowering the bar would retire a real
# guard, so it is skipped explicitly and the run is labelled.
if (!PARTIAL) stopifnot(cov_frac > 0.80)

saveRDS(list(fine = FINE, coarse = SC, summary = SUMM, partial = PARTIAL,
             chroms = CHROMS, params = gam2b_params, built = Sys.time()), out_rds)
msg("  wrote %s", out_rds)

# ---------------------------------------------------------------
# IGV tracks
# ---------------------------------------------------------------
wr_bg <- function(dt, file, name, desc, lo, hi, col = "178,24,43",
                  alt = "33,102,172") {
  p <- file.path(igv_dir, file)
  writeLines(sprintf(paste0('track type=bedGraph name="%s" description="%s" ',
                            'visibility=full color=%s altColor=%s ',
                            'viewLimits=%.3f:%.3f windowingFunction=mean'),
                     name, desc, col, alt, lo, hi), p)
  fwrite(dt, p, sep = "\t", col.names = FALSE, append = TRUE, quote = FALSE,
         scipen = 999)
  msg("  %s (%d rows)", file, nrow(dt))
}

bg <- function(col) FINE[is.finite(get(col)), .(chr, start, end, round(get(col), 4))]
lim <- function(col) as.numeric(quantile(FINE[[col]], c(0.01, 0.99), na.rm = TRUE))

wr_bg(bg("shared_axis"), "gam_shared_axis.bedgraph", "GAM shared A/B axis",
      "averaged correlation matrix, phased on WT H3K4me3+H3K36me3; + = A",
      lim("shared_axis")[1], lim("shared_axis")[2])
wr_bg(bg("score_WT"), "gam_shared_score_WT.bedgraph", "GAM A/B score (WT)",
      "WT correlation structure projected on the shared axis, RMS 1 per chr",
      lim("score_WT")[1], lim("score_WT")[2])
wr_bg(bg("score_HP1gKO"), "gam_shared_score_HP1gKO.bedgraph", "GAM A/B score (cKO)",
      "cKO correlation structure projected on the shared axis, RMS 1 per chr",
      lim("score_HP1gKO")[1], lim("score_HP1gKO")[2])
wr_bg(bg("score_delta"), "gam_shared_score_delta.bedgraph", "GAM A/B score delta",
      "cKO minus WT on the shared axis; shape only, amplitude removed",
      lim("score_delta")[1], lim("score_delta")[2])

# ---------------------------------------------------------------
# Figures
# ---------------------------------------------------------------
local({
  suppressPackageStartupMessages(library(ggplot2))

  # 1. The axis and both scores along two representative chromosomes.
  show <- intersect(c("chr16", "chr11"), SUMM[determined == TRUE, chr])
  if (length(show)) {
    d <- melt(SC[chr %in% show],
              id.vars = c("chr", "cstart"),
              measure.vars = c("shared_axis", "score_WT", "score_HP1gKO"),
              variable.name = "track", value.name = "v")
    d[, track := factor(track, levels = c("shared_axis", "score_WT", "score_HP1gKO"),
                        labels = c("shared axis", "WT score", "cKO score"))]
    p <- ggplot(d, aes(cstart / 1e6, v)) +
      geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
      geom_area(aes(fill = v > 0), colour = NA) +
      scale_fill_manual(values = c("TRUE" = "#b2182b", "FALSE" = "#2166ac"),
                        labels = c("TRUE" = "A side", "FALSE" = "B side"),
                        name = NULL) +
      facet_grid(track ~ chr, scales = "free_x", switch = "y") +
      labs(x = "position (Mb)", y = NULL,
           title = "One shared A/B axis, both genotypes scored on it",
           subtitle = paste0("Axis from the averaged correlation matrix at ",
                             gam2b_params$binsize / 1000,
                             " kb, phased on WT H3K4me3 + H3K36me3. Scores are shape only",
                             " (RMS 1 per chromosome):\namplitude is removed because the ",
                             "~20% fall in lambda1 is global, not compartment-specific.")) +
      theme_m + theme(legend.position = "top",
                      axis.text.y = element_text(size = 7),
                      strip.placement = "outside")
    save_fig(p, "shared_axis_tracks", subdir = fig_sub, width = 11, height = 6)
  }

  # 2. The comparison that justifies the stage: shared axis against each
  #    genotype's own axis, per chromosome, on the same reference.
  if ("rho_pergeno_WT" %in% names(SUMM)) {
    d2 <- melt(SUMM[determined == TRUE,
                    .(chr, `shared axis` = rho_ref, `WT own axis` = rho_pergeno_WT,
                      `cKO own axis` = rho_pergeno_HP1gKO)],
               id.vars = "chr", variable.name = "axis", value.name = "rho")
    ord <- SUMM[determined == TRUE][order(rho_ref), chr]
    d2[, chr := factor(chr, levels = ord)]
    p2 <- ggplot(d2, aes(chr, rho, colour = axis, group = axis)) +
      geom_hline(yintercept = 0.4, linetype = 2, colour = "grey50") +
      geom_line(aes(group = chr), colour = "grey80", linewidth = 0.4) +
      geom_point(size = 2.2) +
      scale_colour_manual(values = c("shared axis" = "#b2182b",
                                     "WT own axis" = "#2166ac",
                                     "cKO own axis" = "#67a9cf"), name = NULL) +
      coord_flip() +
      labs(x = NULL, y = "|rho| vs active chromatin (higher = better determined axis)",
           title = "One shared axis is better determined\nthan either genotype's own",
           subtitle = paste0("Same reference for all three. Dashed line is the |rho| > 0.4 bar.\n",
                             "Per-genotype axes clear it on ",
                             SUMM[is.finite(rho_pergeno_WT) &
                                  (rho_pergeno_WT > 0.4 | rho_pergeno_HP1gKO > 0.4), .N],
                             " of ", SUMM[is.finite(rho_pergeno_WT), .N],
                             " autosomes in either genotype;\nthe shared axis clears it on ",
                             SUMM[determined == TRUE & rho_ref > 0.4, .N],
                             " and beats both on ", SUMM[beats_both == TRUE, .N], ".")) +
      theme_m + theme(legend.position = "top")
    save_fig(p2, "shared_axis_vs_pergenotype", subdir = fig_sub, width = 7.5, height = 7)
  }

  # 3. WT vs cKO score. If the two genotypes express the shared compartment
  #    pattern the same way, this is a tight diagonal; departures from it are
  #    what score_delta picks up.
  p3 <- ggplot(SC, aes(score_WT, score_HP1gKO)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = 2) +
    geom_point(alpha = 0.08, size = 0.5, colour = "#2166ac") +
    labs(x = "WT score on the shared axis", y = "cKO score on the shared axis",
         title = "Both genotypes on one ruler",
         subtitle = sprintf("Spearman rho = %.3f over %d coarse bins; median |delta| = %.3f",
                            stats::cor(SC$score_WT, SC$score_HP1gKO, method = "spearman",
                                       use = "pairwise.complete.obs"),
                            nrow(SC), median(abs(SC$score_delta), na.rm = TRUE))) +
    theme_m
  save_fig(p3, "shared_axis_wt_vs_ko", subdir = fig_sub, width = 6, height = 6)
})

msg("GAM_2b complete.")
