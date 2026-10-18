# ================================================================
# GAM_3 - Does contact remodelling coincide with the mark loss?
# ----------------------------------------------------------------
#     cd GAM && Rscript GAM_3_Coupling.R             (after GAM_1, GAM_2)
#
# Merges the three stages that all answer one question and previously
# contradicted each other ACROSS FILES, where nothing forced a reconciliation.
# Previously GAM_3, GAM_7 and GAM_4.
#
#   A (was GAM_3)  Per-bin Spearman of ChIP log2FC against differential-contact
#                  count. A clean NULL: rho -0.015 raw, -0.046 partial on
#                  coverage, flat under every variant. Bins losing the most
#                  H4K20me3 (median -0.764) and the least (-0.222) have identical
#                  contact densities. n = 119,952 - a powered negative.
#   B (was GAM_7)  Retracts A's premise. The contact variable is integrated over
#                  5 Mb while the ChIP variable varies at 20 kb, a ~250x scale
#                  mismatch, so A cannot distinguish "no relationship" from "this
#                  contact variable cannot resolve that scale". B recounts
#                  locally from the PAIRS; the null survives at every distance
#                  (-0.042 to +0.032).
#   C (was GAM_4)  Contact density per ChIP compartment. Corroborates A
#                  categorically - the three compartments are indistinguishable.
#                  Its compartment-SWITCH half is not usable: the A/B axis is
#                  weakly determined (CLAUDE.md).
#
# WHERE THE POSITIVE RESULT IS: H4K20me3 loss tracks BASELINE architecture, not
# contact CHANGE - pc1 +0.196, npmi_short -0.166 after coverage control.
# Insulation does NOT survive the control (-0.180 to -0.066); do not report it.
#
# STILL OWED, and why this is a structural merge and not a rewrite: the figures.
# They pose a question in the title and answer it in the subtitle with a
# machine-generated verdict string, carry raw column names on the axes
# (n_contacts_lost, partial_baseMean, frac_LINE: Q0_none), need two sentences to
# define the headline y-axis, and magnify a ~0.01 effect to full panel height.
# Having the three framings in one file is the precondition for the fix - one
# synthesis figure plus plain-language axes - not the fix itself.
#
# Output file names are unchanged (coupling_*, gam7_*, gam4_*).
# ================================================================
source("config.R")
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(GenomicRanges); library(IRanges)
  library(patchwork)   # section D composes the synthesis panel; `+` needs it attached
})



# ==============================================================================
# SECTION A - per-bin coupling, all variants               (was GAM_3)
# ==============================================================================

local({
  gam1_rds <- file.path(rds_dir, "gam_bins_annotated.rds")
  gam2_rds <- file.path(rds_dir, "npmi_bin_features.rds")
  if (!file.exists(gam1_rds))
    stop("GAM_1 handoff not found:\n  ", gam1_rds,
         "\nRun it first:   cd GAM && Rscript GAM_1_Bins_and_Signal.R", call. = FALSE)

  # --- Analysis constants (local to this stage; not pipeline parameters) --------
  N_BOOT         <- 500L   # chromosome-bootstrap replicates for rho CIs
  FLAT_RHO       <- 0.05   # |rho| below this counts as "flat" - see section 3
  MIN_STRATUM_N  <- 300L   # drop strata thinner than this (rho unstable below it)
  MIN_N          <- 30L
  N_DECILES      <- 10L
  set.seed(42)

  # ----------------------------------------------------------------
  # Load + validate
  # ----------------------------------------------------------------
  g1 <- readRDS(gam1_rds)
  stopifnot(is.list(g1), all(c("stats", "annot") %in% names(g1)))
  A <- as.data.table(g1$annot)
  S <- as.data.table(g1$stats)
  need_a <- c("bin_id", "chr", "n_contacts_lost", "n_contacts_gained",
              "frac_LINE", "frac_SINE", "frac_LTR")
  need_s <- c("bin_id", "mark", "baseMean", "log2FoldChange")
  stopifnot(all(need_a %in% names(A)), all(need_s %in% names(S)),
            !anyDuplicated(A$bin_id))
  hmm_cols <- grep("^hmm_", names(A), value = TRUE)
  if (!length(hmm_cols)) warning("no hmm_* columns in $annot: ChromHMM stratification skipped")

  # A silent all-zero contact table would read as a clean null for every mark, so
  # fail loudly instead (same reasoning as the zero-overlap asserts elsewhere).
  stopifnot(sum(A$n_contacts_lost) > 0, sum(A$n_contacts_gained) > 0)

  # KO-ORIENTED, matching config.R:468 and GAM_4. This file was last touched
  # 2026-10-17 and was MISSED when the convention was flipped on 2026-10-18, so
  # it shipped the collaborators' WT-minus-KO sign while everything else used
  # ours. The symptom was two opposite numbers for one quantity on one page of
  # the report: coupling_spearman.tsv gave rho(H4K20me3, net) = -0.0104 and the
  # GAM_7 table gave +0.012.
  A[, net := log2((n_contacts_gained + 1) / (n_contacts_lost + 1))]
  # Dominant ChromHMM state per bin = the state with the largest coverage fraction.
  # Chosen over "any overlap" membership because every 20 kb bin overlaps several
  # states, which would make the strata non-exclusive and double-count bins.
  if (length(hmm_cols)) {
    hm <- as.matrix(A[, ..hmm_cols]); hm[is.na(hm)] <- 0
    A[, hmm_dom := ifelse(rowSums(hm) > 0,
                          sub("^hmm_", "", hmm_cols[max.col(hm, ties.method = "first")]),
                          "none")]
  } else A[, hmm_dom := "none"]
  A[, hmm_dom := factor(hmm_dom)]

  S[, mark := factor(as.character(mark), levels = mark_levels)]
  S <- S[!is.na(mark)]
  D <- merge(S, A, by = "bin_id", all.x = TRUE, sort = FALSE)
  stopifnot(nrow(D) == nrow(S), !anyNA(D$chr))   # every stats row found its annot row
  marks <- intersect(mark_levels, as.character(unique(D$mark)))
  by_mark <- setNames(lapply(marks, function(m) D[mark == m]), marks)
  cat(sprintf("\nGAM_3: %d bins x %d marks (%s)\n", nrow(A), length(marks),
              paste(marks, collapse = ", ")))

  # GAM_2 is OPTIONAL: this stage must be useful before the heavy NPMI stage ran.
  F2 <- NULL
  if (file.exists(gam2_rds)) {
    g2 <- readRDS(gam2_rds)
    stopifnot("features" %in% names(g2))
    F2 <- as.data.table(g2$features)
    mm <- mean(A$bin_id %in% F2$bin_id)
    message(sprintf("GAM_2 features found: %.1f%% of autosomal bins matched", 100 * mm))
    stopifnot(mm > 0.95)   # a botched id format gives ~0% - a "null" by accident
  } else {
    message("GAM_2 NPMI features not found (", gam2_rds, ") - skipping section 6. ",
            "Run GAM_2 to enable the pc1/insulation/NPMI/decay couplings.")
  }

  # ----------------------------------------------------------------
  # Statistical helpers
  # ----------------------------------------------------------------
  # Resample CHROMOSOMES, not bins. Neighbouring 20 kb bins are strongly
  # autocorrelated (both log2FC and contact density are smooth along the genome),
  # so the effective n is far below 123k and the textbook CI/p for rho is wildly
  # over-confident. Rejected alternative: per-bin bootstrap (same flaw, treats bins
  # as independent). The chromosome bootstrap with 19 units is coarse but honest.
  # Ranks are computed ONCE globally and the resample takes Pearson on those ranks;
  # this is an approximation to re-ranking each replicate but ~100x cheaper and
  # the difference is negligible at this n.
  chr_boot <- function(chr, f, B = N_BOOT) {
    idx <- split(seq_along(chr), chr); nc <- length(idx)
    vapply(seq_len(B), function(b)
      f(unlist(idx[sample.int(nc, nc, replace = TRUE)], use.names = FALSE)), 0)
  }

  # Spearman rho + analytic p (reported but NOT to be trusted at n~1e5: see above)
  # + chromosome-bootstrap CI. `ranked = TRUE` means x/y are already rank-residuals
  # from a partial correlation (then no re-ranking; k covariates cost df).
  # log10_p is returned because p underflows to 0 for the strong couplings.
  rho_row <- function(x, y, chr, ranked = FALSE, k = 0L, boot = TRUE) {
    ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; chr <- chr[ok]
    n <- length(x)
    empty <- data.table(rho = NA_real_, p = NA_real_, log10_p = NA_real_, n = n,
                        ci_lo = NA_real_, ci_hi = NA_real_)
    if (n < MIN_N || sd(x) == 0 || sd(y) == 0) return(empty)
    if (!ranked) { x <- rank(x); y <- rank(y) }
    r  <- cor(x, y); df <- n - 2 - k
    lp <- pt(-abs(r * sqrt(df / (1 - r^2))), df, log.p = TRUE) + log(2)
    ci <- c(NA_real_, NA_real_)
    if (boot) ci <- quantile(chr_boot(chr, function(i) cor(x[i], y[i])),
                             c(0.025, 0.975), names = FALSE)
    data.table(rho = r, p = min(1, exp(lp)), log10_p = lp / log(10), n = n,
               ci_lo = ci[1], ci_hi = ci[2])
  }

  # Difference in rho between two predictors of the SAME outcome (joint bootstrap,
  # so the shared outcome noise cancels): is lost spatially distinguishable from
  # gained in how it tracks this mark?
  rho_diff <- function(xa, xb, y, chr, ranked = FALSE) {
    ok <- is.finite(xa) & is.finite(xb) & is.finite(y)
    xa <- xa[ok]; xb <- xb[ok]; y <- y[ok]; chr <- chr[ok]
    if (length(y) < MIN_N) return(data.table(diff = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_))
    if (!ranked) { xa <- rank(xa); xb <- rank(xb); y <- rank(y) }
    d  <- cor(xa, y) - cor(xb, y)
    ci <- quantile(chr_boot(chr, function(i) cor(xa[i], y[i]) - cor(xb[i], y[i])),
                   c(0.025, 0.975), names = FALSE)
    data.table(diff = d, ci_lo = ci[1], ci_hi = ci[2])
  }

  # Rank-based partial correlation: regress the RANKS of both variables on the
  # covariates and correlate the residuals. Alternative (ppcor::pcor spearman) does
  # the same arithmetic but would add a dependency for no gain.
  rres <- function(v, X) lm.fit(X, v)$residuals

  # Density bins: zeros get their own level, positives split into n quantile groups
  # (equal-count). Contact counts are integers with a mass at 0 and heavy ties, so
  # plain deciles would split tied values arbitrarily or collapse; this keeps tied
  # bins together and lets "no remodelled contacts" be read off as its own point.
  # Breaks are de-duplicated, so a heavily tied variable yields fewer than n groups.
  dens_bins <- function(x, n = N_DECILES, zero_sep = TRUE) {
    out <- rep(NA_integer_, length(x)); ok <- is.finite(x)
    pos <- ok & (if (zero_sep) x > 0 else TRUE)
    if (zero_sep) out[ok & x <= 0] <- 0L
    if (sum(pos) > 0) {
      br <- unique(quantile(x[pos], seq(0, 1, length.out = n + 1), names = FALSE))
      out[pos] <- if (length(br) < 2) 1L else as.integer(cut(x[pos], br, include.lowest = TRUE))
    }
    out
  }

  # Median CI from order statistics (distribution-free, no resampling needed)
  med_ci <- function(v) {
    v <- sort(v[is.finite(v)]); n <- length(v)
    if (n < 5) return(c(NA_real_, NA_real_))
    c(v[max(1, qbinom(0.025, n, 0.5))], v[min(n, qbinom(0.975, n, 0.5) + 1)])
  }

  # Binned-quantile curve of y (mark log2FC) against a contact variable.
  # WHY CURVES: a 123k-point scatter is an unreadable blob (and mostly ties on x);
  # the per-bin mean/median with a CI shows the dose-response shape directly. The
  # 2D-density plot below is the honest "show me the data" companion.
  qcurve <- function(x, y, zero_sep = TRUE) {
    b <- dens_bins(x, zero_sep = zero_sep); pct <- 100 * ecdf(x[is.finite(x)])(x)
    dt <- data.table(b, x, y, pct)[!is.na(b) & is.finite(y)]
    dt[, {
      m <- mean(y); se <- sd(y) / sqrt(.N); mc <- med_ci(y)
      .(n = .N, mean_contacts = mean(x), pct_mid = mean(pct), mean_l2fc = m,
        mean_lo = m - 1.96 * se, mean_hi = m + 1.96 * se,
        median_l2fc = median(y), median_lo = mc[1], median_hi = mc[2])
    }, by = .(bin = b)][order(bin)]
  }

  fmt_rho <- function(r, lo, hi) sprintf("%+.3f [%+.3f, %+.3f]", r, lo, hi)

  # ================================================================
  # 1 + 2.  COUPLING  (raw, mappable-only, partial) for lost / gained / net
  # ================================================================
  # Variants per mark and direction:
  #   raw                  all bins with a finite log2FC
  #   mappable_only        bins with >=1 differential contact in EITHER direction.
  #                        Zero in both is the signature of an unmappable bin, where
  #                        the contact count is 0 by construction and the ChIP
  #                        log2FC is pure noise - so it is the cheapest mappability
  #                        sensitivity check.
  #   partial_baseMean     controlling for the mark's own coverage (depth/mappability
  #                        proxy; log2FC variance and contact detectability both
  #                        scale with it)
  #   partial_covariates   baseMean + frac_LINE/SINE/LTR + dominant ChromHMM state
  #   given_other          lost | gained and gained | lost. Lost and gained
  #                        are near-equal in TOTAL size (25.6M vs 25.0M pairs), so
  #                        the interesting question is not which is bigger but
  #                        whether they are spatially DISTINGUISHABLE. If lost and
  #                        gained are both just "number of remodelled contacts in
  #                        this bin" they will track each other and any mark, and
  #                        only the conditional (given_other) rho can say whether
  #                        either carries independent information.
  message("\n== 1/2: coupling (Spearman, chromosome-bootstrap CI) ==")
  cov_X <- function(d, which) {
    X <- cbind(1, rank(d$baseMean))
    if (which == "baseMean") return(X)
    cbind(X, rank(d$frac_LINE), rank(d$frac_SINE), rank(d$frac_LTR),
          model.matrix(~ hmm_dom, droplevels(d[, .(hmm_dom)]))[, -1, drop = FALSE])
  }

  cs_list <- list(); diff_list <- list()
  for (m in marks) {
    d <- by_mark[[m]]; y <- d$log2FoldChange; ch <- d$chr
    add <- function(vr, dir, r) cs_list[[length(cs_list) + 1]] <<-
      cbind(data.table(mark = m, variant = vr, direction = dir), r)
    xs <- list(lost = d$n_contacts_lost, gained = d$n_contacts_gained, net = d$net)

    for (dir in names(xs)) add("raw", dir, rho_row(xs[[dir]], y, ch))

    mp <- (d$n_contacts_lost + d$n_contacts_gained) > 0
    for (dir in names(xs)) add("mappable_only", dir, rho_row(xs[[dir]][mp], y[mp], ch[mp]))

    for (wh in c("baseMean", "covariates")) {
      cc <- is.finite(y) & is.finite(d$baseMean) &
        is.finite(d$frac_LINE) & is.finite(d$frac_SINE) & is.finite(d$frac_LTR)
      dd <- d[cc]; X <- cov_X(dd, wh); k <- ncol(X) - 1L
      ry <- rres(rank(dd$log2FoldChange), X)
      for (dir in c("lost", "gained")) {
        xr <- rres(rank(dd[[if (dir == "lost") "n_contacts_lost" else "n_contacts_gained"]]), X)
        add(paste0("partial_", wh), dir, rho_row(xr, ry, dd$chr, ranked = TRUE, k = k))
      }
    }

    # lost | gained and gained | lost (rank-residualised on the other direction)
    cc <- is.finite(y); rl <- rank(d$n_contacts_lost[cc]); rg <- rank(d$n_contacts_gained[cc])
    ry <- rank(y[cc]); chc <- ch[cc]
    Xg <- cbind(1, rg); Xl <- cbind(1, rl)
    add("given_other", "lost",   rho_row(rres(rl, Xg), rres(ry, Xg), chc, ranked = TRUE, k = 1L))
    add("given_other", "gained", rho_row(rres(rg, Xl), rres(ry, Xl), chc, ranked = TRUE, k = 1L))

    # lost vs gained: is one a better spatial predictor of this mark than the other?
    # (raw only: the given_other residuals use different covariate sets per arm, so
    # a paired difference there has no single shared outcome to cancel noise on.)
    diff_list[[length(diff_list) + 1]] <- cbind(
      data.table(mark = m, variant = "raw"),
      rho_diff(d$n_contacts_lost, d$n_contacts_gained, y, ch))
    message("  ", m, " done")
  }
  cs <- rbindlist(cs_list)
  cs[, mark := factor(mark, levels = marks)]
  fwrite(cs, file.path(tables_dir, "coupling_spearman.tsv"), sep = "\t")

  cat("\n=== Spearman rho: per-bin log2FC vs differential-contact density (raw) ===\n")
  cat("    [95% CI = chromosome bootstrap; analytic p at n~1e5 is over-confident, see log10_p]\n")
  print(cs[variant == "raw", .(mark, direction, rho = round(rho, 3),
                               ci = fmt_rho(rho, ci_lo, ci_hi), log10_p = round(log10_p, 1), n)],
        row.names = FALSE)

  # ---- Direction asymmetry: lost vs gained per bin -----------------------------
  message("\n== 2: direction asymmetry ==")
  dd <- rbindlist(diff_list)
  dd[, mark := factor(mark, levels = marks)]
  dd[, lost_minus_gained_distinguishable := is.finite(ci_lo) & (ci_lo > 0 | ci_hi < 0)]
  setnames(dd, "diff", "rho_lost_minus_rho_gained")
  rl_lg <- rho_row(A$n_contacts_lost, A$n_contacts_gained, A$chr)
  asym <- data.table(
    total_lost = sum(A$n_contacts_lost), total_gained = sum(A$n_contacts_gained),
    frac_bins_lost_gt_gained = mean(A$n_contacts_lost > A$n_contacts_gained),
    frac_bins_gained_gt_lost = mean(A$n_contacts_gained > A$n_contacts_lost),
    frac_bins_equal = mean(A$n_contacts_lost == A$n_contacts_gained),
    rho_lost_vs_gained = rl_lg$rho, rho_lost_vs_gained_lo = rl_lg$ci_lo,
    rho_lost_vs_gained_hi = rl_lg$ci_hi,
    net_median = median(A$net), net_iqr = IQR(A$net),
    wilcox_paired_p = tryCatch(suppressWarnings(
      wilcox.test(A$n_contacts_lost, A$n_contacts_gained, paired = TRUE)$p.value),
      error = function(e) NA_real_))
  fwrite(asym, file.path(tables_dir, "coupling_direction_asymmetry.tsv"), sep = "\t")
  fwrite(dd,   file.path(tables_dir, "coupling_direction_difference.tsv"), sep = "\t")
  cat(sprintf(paste0("\nPer-bin lost vs gained: rho = %s; lost>gained in %.1f%% of bins, ",
                     "gained>lost in %.1f%%, equal in %.1f%%\n"),
              fmt_rho(rl_lg$rho, rl_lg$ci_lo, rl_lg$ci_hi), 100 * asym$frac_bins_lost_gt_gained,
              100 * asym$frac_bins_gained_gt_lost, 100 * asym$frac_bins_equal))
  cat("(totals are near-equal, so the question is spatial distinguishability, not size)\n")
  cat("rho(lost) - rho(gained) per mark (paired chromosome bootstrap):\n")
  print(dd[, .(mark, diff = round(rho_lost_minus_rho_gained, 3),
               ci = sprintf("[%+.3f, %+.3f]", ci_lo, ci_hi),
               distinguishable = lost_minus_gained_distinguishable)], row.names = FALSE)

  # ---- Quantile curves (centred on each mark's own median) ----------------------
  # Centred, not raw, log2FC: H4K20me3 sits near -0.7 genome-wide while the active
  # marks sit near 0, so raw curves would be separated by the GLOBAL shift and a
  # flat-but-offset line would look like a difference. Subtracting each mark's own
  # median leaves only the dose-response SHAPE, which is what the specificity
  # control is about (rank correlation is shift-invariant, so this matches rho).
  curves <- rbindlist(lapply(marks, function(m) {
    d <- by_mark[[m]]; med <- median(d$log2FoldChange, na.rm = TRUE)
    rbindlist(lapply(c("lost", "gained", "net"), function(dir) {
      x <- switch(dir, lost = d$n_contacts_lost, gained = d$n_contacts_gained, net = d$net)
      cu <- qcurve(x, d$log2FoldChange, zero_sep = dir != "net")
      cu[, `:=`(mark = m, direction = dir, mark_median = med)]
    }))
  }))
  curves[, mark := factor(mark, levels = marks)]
  for (cc in c("mean_l2fc", "mean_lo", "mean_hi", "median_l2fc", "median_lo", "median_hi"))
    curves[, paste0(cc, "_centred") := get(cc) - mark_median]
  fwrite(curves, file.path(tables_dir, "coupling_quantile_curves.tsv"), sep = "\t")

  curve_plot <- function(dirs, title, subtitle, xlab) {
    ggplot(curves[direction %in% dirs],
           aes(pct_mid, mean_l2fc_centred, colour = mark, fill = mark)) +
      geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
      geom_ribbon(aes(ymin = mean_lo_centred, ymax = mean_hi_centred), alpha = 0.15, colour = NA) +
      geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
      facet_wrap(~ direction, nrow = 1) +
      scale_colour_disc() + scale_fill_disc() +
      labs(title = title, subtitle = subtitle, x = xlab,
           y = "mean log2FC, centred on mark's genome-wide median", colour = NULL, fill = NULL) +
      theme_m + theme(legend.position = "bottom")
  }
  sub_curve <- paste0("bins grouped by contact density (zero-contact bins first, then deciles of the ",
                      "positive counts); ribbon = 95% CI of the mean.\nCentred per mark so the global ",
                      "loss does not offset the lines; a flat line = no coupling.")
  save_fig(curve_plot(c("lost", "gained"), "Mark change vs differential-contact density",
                      sub_curve, "percentile of contact density (bin mean)"),
           "coupling_quantile_curves", subdir = "coupling", width = 10, height = 5)
  save_fig(curve_plot("net", "Mark change vs net contact remodelling",
                      "net = log2((gained+1)/(lost+1)); positive = more contacts gained than lost (more in the cKO)",
                      "percentile of net measure (bin mean)"),
           "coupling_quantile_curves_net", subdir = "coupling", width = 6, height = 5)

  # ---- 2D density (hexbin) ------------------------------------------------------
  # The honest companion to the curves: shows the full joint distribution, including
  # the dense core and the heavy tie-column at x = 0, without 123k overplotted points.
  # hexbin is an optional dependency of ggplot2::geom_hex; fall back to rectangular
  # bins rather than failing the stage over a cosmetic choice.
  have_hex <- requireNamespace("hexbin", quietly = TRUE)
  dens_layer <- function() if (have_hex) geom_hex(bins = 45) else geom_bin2d(bins = 45)
  hexdat <- rbindlist(lapply(marks, function(m) {
    d <- by_mark[[m]][is.finite(log2FoldChange)]
    lim <- quantile(d$log2FoldChange, c(0.005, 0.995))   # clip tails so the core is visible
    d[log2FoldChange >= lim[1] & log2FoldChange <= lim[2]]
  }))
  for (dir in c("lost", "gained")) {
    col <- paste0("n_contacts_", dir)
    hexdat[, xx := log1p(get(col))]
    g <- ggplot(hexdat, aes(xx, log2FoldChange)) + dens_layer() +
      geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 6), se = FALSE,
                  colour = "grey20", linewidth = 0.6) +
      facet_wrap(~ mark, scales = "free_y", nrow = 1) +
      scale_fill_heat0(trans = "log10", name = "bins") +
      labs(title = sprintf("Mark log2FC vs %s-contact density (2D density)", dir),
           subtitle = "log1p(count) on x; y clipped to the 0.5-99.5 percentile per mark; line = GAM smooth",
           x = sprintf("log(1 + n_contacts_%s)", dir), y = "log2FC (HP1gKO / WT)") +
      theme_m
    save_fig(g, paste0("coupling_density2d_", dir), subdir = "coupling", width = 15, height = 4)
  }
  g_lg <- ggplot(A, aes(log1p(n_contacts_lost), log1p(n_contacts_gained))) + dens_layer() +
    geom_abline(slope = 1, intercept = 0, colour = "grey30", linetype = "dashed") +
    scale_fill_heat0(trans = "log10", name = "bins") +
    labs(title = "Per-bin contacts lost vs gained",
         subtitle = sprintf("rho = %.2f; totals are near-equal, the question is whether they are spatially distinct",
                            rl_lg$rho),
         x = "log(1 + n_contacts_lost)", y = "log(1 + n_contacts_gained)") + theme_m
  save_fig(g_lg, "lost_vs_gained_per_bin", subdir = "coupling", width = 6.5, height = 5.5)

  # ================================================================
  # 3.  SPECIFICITY CONTROL  - THE FALSIFICATION TEST
  # ================================================================
  # H4K20me3 (and more mildly H3K9me3) should track contact-loss density; the
  # active/euchromatic marks H3K4me3, H3K36me3, H3K9me2 should be FLAT against it
  # (MINUTE_4 shows them at |log2FC| < 0.1 over every repeat class). If the active
  # marks ALSO track contact loss, the coupling is not about heterochromatin - it
  # is a shared bin-level technical axis (mappability / coverage / depth) that moves
  # every mark's log2FC and the contact counts together. That is the failure mode
  # this test exists to catch, so it is reported whichever way it falls.
  #
  # "Flat" is an EFFECT-SIZE criterion (|rho| < FLAT_RHO), not "p > 0.05": at
  # n ~ 1.2e5 the analytic SE of rho is ~0.003 and spatial autocorrelation does not
  # rescue it, so a p-value rule would call a rho of 0.02 "coupled". The bootstrap CI
  # is reported next to it so the reader can judge for themselves. FLAT_RHO = 0.05
  # (explaining <0.3% of rank variance) is a judgement call and is stated in the output.
  #
  # The verdict is computed on BOTH the raw and the baseMean-partial correlations,
  # because a control that is flat only after removing coverage would mean the
  # raw coupling is partly a coverage effect, and the reader should see that.
  message("\n== 3: SPECIFICITY CONTROL (falsification test) ==")
  role <- c(H3K4me3 = "control", H3K36me3 = "control", H3K9me2 = "control",
            H3K9me3 = "mild_target", H4K20me3 = "target")

  judge <- function(vr) {
    r <- cs[variant == vr & direction == "lost"]
    g <- cs[variant == vr & direction == "gained", .(mark, rho_gained = rho,
                                                      gained_lo = ci_lo, gained_hi = ci_hi)]
    r <- merge(r, g, by = "mark")[order(mark)]
    r[, role := role[as.character(mark)]]
    r[, status := ifelse(!is.finite(rho), "NA",
                  ifelse(abs(rho) < FLAT_RHO, "flat",
                         ifelse(ci_lo > 0 | ci_hi < 0,
                                ifelse(rho < 0, "coupled_negative", "coupled_positive"),
                                "nonflat_CI_includes_0")))]
    tgt <- r[mark == "H4K20me3"]
    ctl <- r[role == "control"]
    verdict <- if (!nrow(tgt) || !nrow(ctl)) {
      "UNDETERMINED: H4K20me3 or the active-mark controls are missing from the input"
    } else {
      tgt_coupled <- abs(tgt$rho) >= FLAT_RHO && (tgt$ci_lo > 0 || tgt$ci_hi < 0)
      bad <- ctl[abs(rho) >= FLAT_RHO]
      weaker <- nrow(bad) > 0 && all(abs(bad$rho) < 0.5 * abs(tgt$rho))
      nm <- function(x) paste(sprintf("%s (rho %+.3f)", x$mark, x$rho), collapse = ", ")
      if (tgt_coupled && !nrow(bad))
        sprintf("SPECIFIC: all active-mark controls flat (|rho| < %.2f); H4K20me3 coupled (rho %+.3f).",
                FLAT_RHO, tgt$rho)
      else if (tgt_coupled && weaker)
        sprintf("AMBIGUOUS: control(s) %s are not flat but < half the H4K20me3 effect (rho %+.3f); a partial artifact component cannot be excluded.",
                nm(bad), tgt$rho)
      else if (tgt_coupled)
        sprintf("NOT SPECIFIC: control(s) %s are not flat (H4K20me3 rho %+.3f); coupling is not distinguishable from a mappability/coverage artifact.",
                nm(bad), tgt$rho)
      else if (nrow(bad))
        sprintf("NO H4K20me3 COUPLING and controls %s not flat: bin-level technical structure is present and H4K20me3 does not exceed it.", nm(bad))
      else
        sprintf("NO COUPLING: H4K20me3 (rho %+.3f) and all active-mark controls are flat.", tgt$rho)
    }
    r[, ratio_to_H4K20me3 := abs(rho) / abs(tgt$rho)]
    r[, same_sign_as_H4K20me3 := sign(rho) == sign(tgt$rho)]
    r[, overall_verdict := verdict]
    r[, FLAT_RHO := FLAT_RHO]
    r[, variant := vr]
    r[, .(variant, mark, role, rho_lost = rho, ci_lo, ci_hi, log10_p, n,
          rho_gained, gained_lo, gained_hi, status, ratio_to_H4K20me3,
          same_sign_as_H4K20me3, FLAT_RHO, overall_verdict)]
  }
  spec <- rbindlist(lapply(c("raw", "mappable_only", "partial_baseMean", "partial_covariates"), judge))

  # Per-chromosome replication of the sign: a coupling that is real should point the
  # same way on most of the 19 autosomes, not be carried by one or two.
  chr_rho <- rbindlist(lapply(marks, function(m) {
    d <- by_mark[[m]]
    rbindlist(lapply(c("lost", "gained"), function(dir) {
      d[, rho_row(get(paste0("n_contacts_", dir)), log2FoldChange, chr, boot = FALSE), by = chr
      ][, `:=`(mark = m, direction = dir)]
    }))
  }))
  fwrite(chr_rho, file.path(tables_dir, "coupling_per_chromosome.tsv"), sep = "\t")
  chr_sum <- chr_rho[direction == "lost", .(n_chr = sum(is.finite(rho)),
                                             n_chr_negative = sum(rho < 0, na.rm = TRUE),
                                             chr_rho_min = min(rho, na.rm = TRUE),
                                             chr_rho_max = max(rho, na.rm = TRUE)), by = mark]
  spec <- merge(spec, chr_sum, by = "mark", all.x = TRUE, sort = FALSE)
  spec[, mark := factor(mark, levels = marks)]
  spec[, variant := factor(variant, levels = c("raw", "mappable_only", "partial_baseMean", "partial_covariates"))]
  setorder(spec, variant, mark)
  spec[, variant := as.character(variant)]
  fwrite(spec, file.path(tables_dir, "coupling_specificity.tsv"), sep = "\t")

  cat("\n", strrep("=", 74), "\n SPECIFICITY CONTROL: contact-LOSS density vs per-bin log2FC\n",
      strrep("=", 74), "\n", sep = "")
  cat(sprintf(" flat := |rho| < %.2f (effect-size rule; p-values are uninformative at n~1e5)\n\n", FLAT_RHO))
  for (vr in unique(spec$variant)) {
    s <- spec[variant == vr]
    cat(sprintf("--- %s ---\n", vr))
    print(s[, .(mark, role, rho_lost = fmt_rho(rho_lost, ci_lo, ci_hi), status,
                rho_gained = round(rho_gained, 3),
                `x H4K20me3` = round(ratio_to_H4K20me3, 2),
                chr_neg = sprintf("%d/%d", n_chr_negative, n_chr))], row.names = FALSE)
    cat(" VERDICT [", vr, "]: ", s$overall_verdict[1], "\n\n", sep = "")
  }
  cat(strrep("=", 74), "\n")

  verdict_raw <- spec[variant == "raw", overall_verdict][1]
  verdict_adj <- spec[variant == "partial_baseMean", overall_verdict][1]
  if (!identical(substr(verdict_raw, 1, 8), substr(verdict_adj, 1, 8)))
    message("NOTE: the verdict CHANGES after controlling for baseMean (raw vs partial) - ",
            "the raw coupling is partly a coverage effect. Read both.")

  spec_plot <- cs[variant %in% c("raw", "partial_baseMean", "partial_covariates") &
                    direction %in% c("lost", "gained")]
  spec_plot[, variant := factor(variant, levels = c("raw", "partial_baseMean", "partial_covariates"))]
  g_spec <- ggplot(spec_plot, aes(mark, rho, colour = mark, shape = variant)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -FLAT_RHO, ymax = FLAT_RHO,
             fill = "grey85", alpha = 0.6) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey45") +
    geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), position = position_dodge(width = 0.6),
                    linewidth = 0.5) +
    facet_wrap(~ direction, nrow = 1) +
    scale_colour_disc(guide = "none") +
    labs(title = "Specificity control: does only H4K20me3 (+/- H3K9me3) track contact remodelling?",
         subtitle = paste0(paste(strwrap(paste("raw:", verdict_raw), 120), collapse = "\n"), "\n",
                           "grey band = 'flat' (|rho| < ", FLAT_RHO, "); bars = chromosome-bootstrap 95% CI"),
         x = NULL, y = "Spearman rho (log2FC vs contact density)", shape = NULL) +
    theme_m + theme(legend.position = "bottom", plot.subtitle = element_text(size = 8))
  save_fig(g_spec, "specificity_control_rho", subdir = "coupling", width = 11, height = 5.5)

  # All five marks on ONE panel for direct comparison, lost only (the control as
  # specified). Same data as the curves figure, restricted to the direction at issue.
  g_spec_curve <- curve_plot("lost", "Specificity control: all five marks vs contact-LOSS density",
                             paste(strwrap(paste("raw:", verdict_raw), 110), collapse = "\n"),
                             "percentile of n_contacts_lost (bin mean)") +
    theme(plot.subtitle = element_text(size = 8))
  save_fig(g_spec_curve, "specificity_control_curves", subdir = "coupling", width = 7, height = 5.5)

  # ================================================================
  # 4.  CONFOUNDER STRATIFICATION  (raw AND stratified)
  # ================================================================
  # A coupling that survives stratification is a much stronger claim than a raw
  # one: ChromHMM state, repeat content and coverage each plausibly drive BOTH the
  # contact counts (mappability, heterochromatin compaction) and the ChIP log2FC
  # (H4K20me3 loss is deepest over repeats/quiescent chromatin per MINUTE_3/4), so
  # the raw rho could be pure confounding. Per stratum we recompute rho; the
  # within-stratum estimates are then pooled by inverse-variance on Fisher-z
  # (weights n-3), giving ONE stratified rho comparable to the raw one. The pooled
  # SE ignores spatial autocorrelation, so only its estimate is meaningful, not its p.
  # Strata are defined on the CONFOUNDER only, never on the outcome: stratifying on
  # compartment_*/log2FC itself would truncate the outcome range and bias rho to 0.
  message("\n== 4: confounder stratification ==")
  strat_of <- function(d) {
    s <- list()
    if (length(hmm_cols)) s[["ChromHMM_dominant_state"]] <- as.character(d$hmm_dom)
    for (rp in c("LINE", "SINE", "LTR")) {
      b <- dens_bins(d[[paste0("frac_", rp)]], n = 4)
      s[[paste0("frac_", rp)]] <- ifelse(is.na(b), NA_character_, ifelse(b == 0, "Q0_none", paste0("Q", b)))
    }
    b <- dens_bins(d$baseMean, n = 4, zero_sep = FALSE)
    s[["baseMean"]] <- ifelse(is.na(b), NA_character_, paste0("Q", b))
    s
  }
  strat_rows <- list()
  for (m in marks) {
    d <- by_mark[[m]]; ss <- strat_of(d)
    for (dir in c("lost", "gained")) {
      x <- d[[paste0("n_contacts_", dir)]]; y <- d$log2FoldChange
      r0 <- rho_row(x, y, d$chr, boot = FALSE)
      strat_rows[[length(strat_rows) + 1]] <- cbind(
        data.table(mark = m, direction = dir, confounder = "none", stratum = "<raw>"), r0)
      for (cf in names(ss)) {
        st <- ss[[cf]]
        rows <- rbindlist(lapply(sort(unique(na.omit(st))), function(lv) {
          i <- which(st == lv)
          if (length(i) < MIN_STRATUM_N) return(NULL)
          cbind(data.table(stratum = lv), rho_row(x[i], y[i], d$chr[i], boot = FALSE))
        }))
        if (!nrow(rows)) next
        rows[, `:=`(mark = m, direction = dir, confounder = cf)]
        strat_rows[[length(strat_rows) + 1]] <- rows[, names(strat_rows[[1]]), with = FALSE]
        ok <- rows[is.finite(rho) & abs(rho) < 1]
        if (nrow(ok)) {
          w <- ok$n - 3; z <- sum(w * atanh(ok$rho)) / sum(w); se <- 1 / sqrt(sum(w))
          strat_rows[[length(strat_rows) + 1]] <- data.table(
            mark = m, direction = dir, confounder = cf, stratum = "<pooled within>",
            rho = tanh(z), p = 2 * pnorm(-abs(z / se)), log10_p = 2 * pnorm(-abs(z / se), log.p = TRUE) / log(10),
            n = sum(ok$n), ci_lo = tanh(z - 1.96 * se), ci_hi = tanh(z + 1.96 * se))
        }
      }
    }
  }
  strat <- rbindlist(strat_rows, use.names = TRUE)
  strat[, mark := factor(mark, levels = marks)]
  fwrite(strat, file.path(tables_dir, "coupling_stratified.tsv"), sep = "\t")

  cat("\n=== Raw vs pooled-within-stratum rho (lost); a coupling that survives is stronger ===\n")
  cmp <- dcast(strat[direction == "lost" & stratum %in% c("<raw>", "<pooled within>")],
               mark ~ confounder, value.var = "rho")
  cmp[, (setdiff(names(cmp), "mark")) := lapply(.SD, round, 3), .SDcols = setdiff(names(cmp), "mark")]
  setnames(cmp, "none", "raw")
  print(cmp, row.names = FALSE)
  cat("(partial correlations controlling baseMean / baseMean+repeats+ChromHMM are in coupling_spearman.tsv)\n")
  cat("\nPartial rho, lost vs log2FC:\n")
  print(dcast(cs[direction == "lost" & variant %in% c("raw", "mappable_only", "partial_baseMean",
                                                      "partial_covariates", "given_other")],
              mark ~ variant, value.var = "rho")[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)],
        row.names = FALSE)

  # Strata-level tile figure: rho per (confounder:stratum) x mark, lost | gained
  sp <- strat[stratum != "<pooled within>"][, lab := ifelse(stratum == "<raw>", "raw",
                                                          paste0(confounder, ": ", stratum))]
  sp[, lab := factor(lab, levels = rev(unique(lab)))]
  # FIXED LIMIT, not max(abs(rho)). The data-driven limit was being set by the
  # SMALLEST strata - ChromHMM Het is 666 bins, 0.54% of the genome - which
  # painted every other tile near-white and made the panel read as "nothing
  # couples anywhere".
  lim <- 0.10
  if (nrow(sp) && "n" %in% names(sp))
    sp[, lab := sprintf("%s  (n=%s)", lab, format(n, big.mark = ","))]
  sp[, lab := factor(lab, levels = rev(unique(lab)))]
  g_strat <- ggplot(sp, aes(mark, lab, fill = pmax(-lim, pmin(lim, rho)))) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%+.2f", rho)), size = 2.4) +
    facet_wrap(~ direction, nrow = 1) +
    scale_fill_heat0_div(limits = c(-lim, lim), name = "Spearman\nrho") +
    labs(title = "Does contact change predict mark change WITHIN matched chromatin?",
         x = NULL, y = NULL,
         subtitle = paste("NOT a map of where contacts change. Each tile is the rho between per-bin ChIP log2FC and",
                          "per-bin contact count\nINSIDE one stratum, so any association that defines the stratum",
                          "(LINE content, ChromHMM state) is removed by construction.\nContact loss IS strongly",
                          "LINE-rich - mean lost counts rise 18.3 to 20.8 from LINE Q1 to Q4 - and that contrast is",
                          "between rows, not within them.\nFill clipped at +-0.10; values printed in full.")) +
    theme_m + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))
  save_fig(g_strat, "coupling_stratified", subdir = "coupling", width = 10, height = 9)

  pp <- strat[stratum %in% c("<raw>", "<pooled within>")]
  pp[, what := ifelse(stratum == "<raw>", "raw", paste0("pooled: ", confounder))]
  g_pool <- ggplot(pp, aes(mark, rho, colour = what)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -FLAT_RHO, ymax = FLAT_RHO, fill = "grey85", alpha = 0.5) +
    geom_point(position = position_dodge(width = 0.7), size = 2) +
    facet_wrap(~ direction, nrow = 1) + scale_colour_disc() +
    labs(title = "Raw vs stratified (pooled within strata) coupling", x = NULL,
         y = "Spearman rho", colour = NULL) + theme_m
  save_fig(g_pool, "coupling_raw_vs_pooled", subdir = "coupling", width = 10, height = 5)

  # ================================================================
  # 5.  CONTACT DENSITY BY COMPARTMENT  (outcome side, descriptive)
  # ================================================================
  # compartment_centred is the PRIMARY labelling (the absolute cutoffs are not
  # shift-invariant against a genome-wide loss; see ChIP/results/COMPARTMENT_CUTOFFS.md).
  # Here the compartment is the GROUPING and contact density the response, so it is
  # legitimate - unlike stratifying the rho by compartment, which conditions on the
  # outcome. Effect size is the AUC (P[lost in group > lost in stable]), 0.5 = none.
  message("\n== 5: contact density by compartment ==")
  comp_rows <- list()
  for (cc in intersect(c("compartment_centred", "compartment_abs"), names(A))) {
    a <- A[!is.na(get(cc))]
    if (!nrow(a) || uniqueN(a[[cc]]) < 2) next
    for (dir in c("lost", "gained", "net")) {
      v <- switch(dir, lost = a$n_contacts_lost, gained = a$n_contacts_gained, net = a$net)
      g <- factor(a[[cc]]); kw <- tryCatch(kruskal.test(v ~ g)$p.value, error = function(e) NA_real_)
      for (lv in levels(g)) {
        i <- g == lv; ref <- g == "stable"
        auc <- if (lv != "stable" && any(ref) && any(i)) {
          w <- suppressWarnings(wilcox.test(v[i], v[ref]))
          unname(w$statistic) / (sum(i) * sum(ref))
        } else NA_real_
        comp_rows[[length(comp_rows) + 1]] <- data.table(
          scheme = cc, direction = dir, group = lv, n_bins = sum(i),
          mean = mean(v[i]), median = median(v[i]), auc_vs_stable = auc, kruskal_p = kw)
      }
    }
  }
  if (length(comp_rows)) {
    comp <- rbindlist(comp_rows)
    fwrite(comp, file.path(tables_dir, "coupling_by_compartment.tsv"), sep = "\t")
    cat("\nContact density by compartment_centred (AUC vs stable; 0.5 = no difference):\n")
    print(comp[scheme == "compartment_centred" & direction %in% c("lost", "gained"),
               .(direction, group, n_bins, mean = round(mean, 2), median, auc_vs_stable = round(auc_vs_stable, 3))],
          row.names = FALSE)
    pc <- A[!is.na(compartment_centred)][, .(compartment_centred, lost = n_contacts_lost,
                                             gained = n_contacts_gained)]
    pc <- melt(pc, id.vars = "compartment_centred", variable.name = "direction", value.name = "count")
    g_comp <- ggplot(pc, aes(compartment_centred, log1p(count), fill = compartment_centred)) +
      geom_violin(colour = NA, alpha = 0.6, scale = "width") +
      geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.3) +
      facet_wrap(~ direction, nrow = 1) + scale_fill_disc(guide = "none") +
      labs(title = "Differential-contact density by ChIP compartment (median-centred)",
           x = NULL, y = "log(1 + contacts per bin)") + theme_m
    save_fig(g_comp, "contacts_by_compartment", subdir = "coupling", width = 8, height = 4.5)
  } else message("  no compartment labels in $annot - skipped")

  # ================================================================
  # 6.  NPMI-DERIVED FEATURES  (only if GAM_2 has run)
  # ================================================================
  # CAVEAT carried into the output: the contact counts and these features are both
  # derived from the SAME NPMI matrices, so contact density vs insulation/decay is
  # internal consistency, not independent evidence. The independent test is
  # ChIP log2FC vs the feature, which is what the table leads with.
  if (!is.null(F2)) {
    message("\n== 6: NPMI-derived features ==")
    Fj <- copy(F2[match(A$bin_id, bin_id)])
    Fj[, bin_id := A$bin_id]
    for (b in c("npmi_short", "npmi_long", "decay_slope", "pc1", "insulation")) {
      cw <- paste0(b, "_WT"); ck <- paste0(b, "_HP1gKO"); cd <- paste0(b, "_delta")
      if (all(c(cw, ck) %in% names(Fj)) && !(cd %in% names(Fj))) Fj[, (cd) := get(ck) - get(cw)]
    }
    feat_all <- c("pc1_delta", "insulation_delta", "npmi_short_delta", "npmi_long_delta",
                  "decay_slope_delta",
                  unlist(lapply(c("WT", "HP1gKO"), function(g)
                    paste0(c("pc1", "insulation", "npmi_short", "npmi_long", "decay_slope"), "_", g))))
    feats <- intersect(feat_all, names(Fj))
    if (length(setdiff(feat_all, feats)))
      message("  GAM_2 lacks: ", paste(setdiff(feat_all, feats), collapse = ", "))

    fr <- rbindlist(lapply(marks, function(m) {
      d <- by_mark[[m]]; i <- match(d$bin_id, Fj$bin_id)
      rbindlist(lapply(feats, function(f)
        cbind(data.table(mark = m, feature = f), rho_row(Fj[[f]][i], d$log2FoldChange, d$chr))))
    }))
    fr[, mark := factor(mark, levels = marks)]
    fwrite(fr, file.path(tables_dir, "coupling_npmi_features.tsv"), sep = "\t")

    # contact density vs feature: the internal-consistency check (see caveat)
    fc <- rbindlist(lapply(c("lost", "gained"), function(dir)
      rbindlist(lapply(feats, function(f)
        cbind(data.table(direction = dir, feature = f),
              rho_row(A[[paste0("n_contacts_", dir)]], Fj[[f]], A$chr, boot = FALSE))))))
    fwrite(fc, file.path(tables_dir, "coupling_contacts_vs_npmi_features.tsv"), sep = "\t")

    cat("\nlog2FC vs NPMI feature, Spearman rho (independent of the contact calls):\n")
    print(dcast(fr[feature %in% grep("_delta$", feats, value = TRUE)], feature ~ mark,
                value.var = "rho")[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)],
          row.names = FALSE)

    fr[, feature := factor(feature, levels = rev(feats))]
    g_feat <- ggplot(fr, aes(mark, feature, fill = rho)) + geom_tile(colour = "white") +
      geom_text(aes(label = sprintf("%+.2f", rho)), size = 2.6) +
      scale_fill_heat0_div(limits = c(-1, 1) * max(abs(fr$rho), na.rm = TRUE), name = "Spearman\nrho") +
      labs(title = "Mark log2FC vs NPMI-derived 3D features", x = NULL, y = NULL,
           subtitle = "deltas are HP1gKO - WT; per-genotype rows are the baseline state") +
      theme_m + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))
    save_fig(g_feat, "coupling_npmi_features", subdir = "coupling", width = 7, height = 6)

    if ("compartment_switch" %in% names(Fj)) {
      sw <- rbindlist(lapply(marks, function(m) {
        d <- by_mark[[m]]; gm <- median(d$log2FoldChange, na.rm = TRUE)
        d <- data.table(l2fc = d$log2FoldChange,
                        sw = Fj$compartment_switch[match(d$bin_id, Fj$bin_id)])[!is.na(sw)]
        d[, .(mark = m, n = .N, median_l2fc = median(l2fc, na.rm = TRUE),
              centred = median(l2fc, na.rm = TRUE) - gm), by = .(switch = sw)]
      }))
      fwrite(sw, file.path(tables_dir, "coupling_compartment_switch.tsv"), sep = "\t")
      sw[, mark := factor(mark, levels = marks)]
      g_sw <- ggplot(sw, aes(switch, centred, fill = mark)) + geom_col(position = "dodge") +
        geom_hline(yintercept = 0, linewidth = 0.3) + scale_fill_disc() +
        labs(title = "Median mark change by A/B compartment switch (WT -> HP1gKO)",
             x = "compartment switch", y = "median log2FC, centred on mark median", fill = NULL) + theme_m
      save_fig(g_sw, "coupling_compartment_switch", subdir = "coupling", width = 8, height = 4.5)
    }
  }

  # ---- Session summary ---------------------------------------------------------
  cat("\nGAM_3 complete. Headline (raw): ", verdict_raw, "\n", sep = "")
  cat("             (baseMean-partial): ", verdict_adj, "\n", sep = "")
  cat("Tables  -> ", tables_dir, "\nFigures -> ", file.path(fig_dir, "coupling"), "\n", sep = "")
})


# ==============================================================================
# SECTION B - local recount from the pairs, by distance    (was GAM_7)
# ==============================================================================

local({
  g1 <- file.path(rds_dir, "gam_bins_annotated.rds")
  if (!file.exists(g1)) stop("GAM_1 output not found:\n  ", g1,
    "\nRun it first:  cd GAM && Rscript GAM_1_Bins_and_Signal.R", call. = FALSE)
  x <- readRDS(g1)
  a <- as.data.table(x$annot); setorder(a, chr, start)
  for (mk in colnames(x$l2fc)) a[[paste0("l2fc_", mk)]] <- x$l2fc[match(a$bin_id, rownames(x$l2fc)), mk]

  # Distance bands. The Hotspots track collapses everything under 5 Mb into one
  # number and discards everything above it; the pairs run to 57 Mb. These bands
  # separate the scales a compaction change and a compartment change would act on.
  BANDS <- list("0.02-1Mb" = c(2e4, 1e6), "1-5Mb" = c(1e6, 5e6),
                "5-10Mb" = c(5e6, 1e7), ">10Mb" = c(1e7, Inf))

  bin_index <- function(chrom, pos) {
    # pos is a 0-based start from the bg2; bins are gam_binsize wide from 0
    paste0(chrom, ":", (pos %/% gam_binsize) * gam_binsize)
  }
  a[, key := paste0(chr, ":", start)]
  setkey(a, key)

  count_one <- function(chrom, dir) {
    f <- diff_file(chrom, dir)
    if (!file.exists(f)) return(NULL)
    d <- fread(f, select = c("chrom_x","start_x","chrom_y","start_y","bindist"),
               showProgress = FALSE)
    if (!nrow(d)) return(NULL)
    d[, band := cut(bindist, breaks = c(0, 1e6, 5e6, 1e7, Inf),
                    labels = names(BANDS), right = FALSE)]
    # A pair contributes to BOTH of its endpoint bins: the question is "how many
    # differential contacts does this bin participate in", so each end counts.
    both <- rbindlist(list(
      d[, .(key = bin_index(chrom_x, start_x), band)],
      d[, .(key = bin_index(chrom_y, start_y), band)]))
    out <- both[, .N, by = .(key, band)]
    rm(d, both); gc(FALSE)
    out
  }

  cat("Counting pairs per bin, per distance band (", length(gam_autosomes), "chromosomes )\n")
  res <- list()
  for (ch in gam_autosomes) {
    for (dir in c("lost", "gained")) {
      o <- count_one(ch, dir)
      if (!is.null(o)) { o[, dirn := dir]; res[[paste(ch, dir)]] <- o }
    }
    cat("  ", ch, "done\n")
  }
  cnt <- rbindlist(res); rm(res); gc(FALSE)

  # wide: one column per direction x band
  w <- dcast(cnt, key ~ dirn + band, value.var = "N", fill = 0)
  a <- merge(a, w, by = "key", all.x = TRUE)
  for (j in setdiff(names(w), "key")) set(a, which(is.na(a[[j]])), j, 0L)

  a[, local_lost   := rowSums(.SD), .SDcols = grep("^lost_",   names(a), value = TRUE)]
  a[, local_gained := rowSums(.SD), .SDcols = grep("^gained_", names(a), value = TRUE)]

  cat("\n=== local (pairs) vs delivered (0-5Mb window) contact counts ===\n")
  cat(sprintf("  delivered lost : median %4d  IQR %5.0f\n", median(a$n_contacts_lost), IQR(a$n_contacts_lost)))
  cat(sprintf("  local     lost : median %4d  IQR %5.0f\n", median(a$local_lost), IQR(a$local_lost)))
  cat(sprintf("  rho(delivered, local) = %+.3f  <- how much the 5 Mb window smears\n",
              cor(a$n_contacts_lost, a$local_lost, method = "spearman")))

  # --- the re-test: does a LOCAL contact measure couple to the ChIP change? -----
  marks <- colnames(x$l2fc)
  cp <- rbindlist(lapply(marks, function(mk) {
    v <- a[[paste0("l2fc_", mk)]]
    rbindlist(lapply(c(list(all = NULL), BANDS), function(bd) NULL))  # placeholder
  }))
  rows <- list()
  for (mk in marks) {
    v <- a[[paste0("l2fc_", mk)]]
    for (dir in c("lost","gained")) {
      rows[[length(rows)+1]] <- data.table(mark = mk, direction = dir, band = "all (local)",
        rho = cor(v, a[[paste0("local_", dir)]], method = "spearman", use = "complete.obs"))
      for (bn in names(BANDS)) {
        col <- paste0(dir, "_", bn)
        if (col %in% names(a))
          rows[[length(rows)+1]] <- data.table(mark = mk, direction = dir, band = bn,
            rho = cor(v, a[[col]], method = "spearman", use = "complete.obs"))
      }
    }
    rows[[length(rows)+1]] <- data.table(mark = mk, direction = "delivered", band = "0-5Mb window",
      rho = cor(v, a$n_contacts_lost, method = "spearman", use = "complete.obs"))
  }
  cp <- rbindlist(rows)
  fwrite(cp, file.path(tables_dir, "gam7_coupling_by_distance.tsv"), sep = "\t")

  cat("\n=== coupling of per-bin ChIP log2FC to LOCAL contact counts, by distance ===\n")
  pr <- dcast(cp[direction != "delivered"], mark + direction ~ band, value.var = "rho")
  print(pr, row.names = FALSE, digits = 3)
  cat("\n  (GAM_3's delivered-track value, for comparison)\n")
  print(cp[direction == "delivered", .(mark, rho = round(rho, 3))], row.names = FALSE)

  saveRDS(list(bins = a, coupling = cp, bands = BANDS),
          file.path(rds_dir, "gam7_local_contacts.rds"))
  fwrite(a[, c("bin_id","chr","start","end","local_lost","local_gained",
               grep("^(lost|gained)_", names(a), value = TRUE)), with = FALSE],
         file.path(tables_dir, "gam7_local_contact_counts.tsv.gz"), sep = "\t")
  cat("\nGAM_7 complete.\n")
})


# ==============================================================================
# SECTION C - compartments vs contacts                     (was GAM_4)
# ==============================================================================

local({
  # --- Inputs ------------------------------------------------------
  gam1_rds <- file.path(rds_dir, "gam_bins_annotated.rds")
  gam2_rds <- file.path(rds_dir, "npmi_bin_features.rds")
  if (!file.exists(gam1_rds))
    stop("GAM_1 output not found:\n  ", gam1_rds,
         "\nRun it first:   cd GAM && Rscript GAM_1_Bins_and_Signal.R", call. = FALSE)
  g1 <- readRDS(gam1_rds)
  d  <- copy(as.data.table(g1$annot))
  stopifnot(all(c("bin_id", "n_contacts_lost", "n_contacts_gained", "tad_boundary",
                  "compartment_abs", "compartment_centred") %in% names(d)))
  # H4K20me3 per-bin log2FC from GAM_1's wide matrix, joined on bin_id (not row
  # order: the matrix and annot are separate objects)
  d[, h4_l2fc := g1$l2fc[bin_id, "H4K20me3"]]

  have_feat <- file.exists(gam2_rds)
  if (have_feat) {
    g2 <- readRDS(gam2_rds)
    want <- intersect(c("bin_id", "pc1_WT", "pc1_HP1gKO", "pc1_delta", "compartment_WT",
                        "compartment_HP1gKO", "compartment_switch", "insulation_delta"),
                      names(g2$features))
    d <- merge(d, g2$features[, ..want], by = "bin_id", all.x = TRUE, sort = FALSE)
    # features also cover chrX/Y; the left join keeps only the autosomal GAM_1 grid
  } else {
    message("GAM_4: ", gam2_rds, " not found - skipping the A/B-switch and insulation",
            " sections (run GAM_2 first: cd GAM && Rscript GAM_2_NPMI_Compartments.R).")
  }
  grp_levels <- c("stable", "H4K20me3_only", "co_loss")
  grp_cols   <- setNames(disc_pal(3), grp_levels)

  # --- Scheme choice: READ THIS before changing the grouping --------
  # PRIMARY = compartment_centred; compartment_abs is a SENSITIVITY only.
  # Both marks lose genome-wide (H4K20me3 median log2FC ~ -0.68, H3K9me3 ~ -0.42), so
  # the absolute cutoffs (-0.5 / -0.3 / 0.15) mostly ask "did it lose?", to which most
  # of the genome says yes. Group sizes and even their ORDERING then depend on where
  # the global shift sits (co_loss 46,921 vs H4K20me3_only 8,481 under abs; the
  # reverse once each mark is centred on its own median) - see
  # ChIP/results/COMPARTMENT_CUTOFFS.md. A 3D result that held only under the
  # absolute labels would be a statement about the global shift, not about the
  # compartments, so we report both and only trust what survives the centred scheme.
  schemes <- c(centred = "compartment_centred", abs = "compartment_abs")
  PRIMARY <- "centred"

  # --- Matching covariate for size_perm() --------------------------
  # size_perm() was written for REGIONS OF VARYING WIDTH and matches on width
  # deciles. Every GAM bin is 20 kb, so a width match is degenerate (one unique
  # break -> cut() errors). The same machinery works for any numeric nuisance
  # covariate passed in the "width" slots, so we match on per-bin REPEAT CONTENT
  # instead: repeats drive both the H4K20me3/H3K9me3 labels (IAP/L1 are the losers)
  # and mappability/contact counts, so a naive group-vs-all comparison would partly
  # be a repeat-density comparison. A tiny jitter breaks the heavy ties at 0 so the
  # deciles are unique; seeded so the null is reproducible.
  d[, rep_cov := rowSums(.SD, na.rm = TRUE), .SDcols = c("frac_LINE", "frac_SINE", "frac_LTR")]
  matched_perm <- function(idx, flag, cov, nperm = 2000, seed = 42) {
    ok <- is.finite(cov) & !is.na(flag)
    idx <- idx[ok[idx]]
    if (length(idx) < 10) return(NULL)
    set.seed(seed); cj <- cov[ok] + runif(sum(ok), 0, 1e-6)
    map <- cumsum(ok)                    # row -> position within the ok-universe
    qi  <- map[idx]
    size_perm(mean(flag[ok][qi]), cj[qi], cj, as.numeric(flag[ok]), nperm = nperm, seed = seed)
  }

  # Continuous version: is the group's MEAN of x different from a repeat-content-matched
  # draw from `universe` (default all bins)? size_perm() only returns an upper-tail p,
  # so run it on x and -x (same seed -> identical draws) for a two-sided p.
  # WHY NOT a Wilcoxon of group vs stable: the delivered GAM data merge 2 replicates
  # per genotype, so there is no within-genotype variance, and 123k bins are not
  # independent observations (neighbouring bins share contacts). A bin-level p-value
  # would be wildly anti-conservative; the matched null asks the right question
  # ("more than bins of the same repeat content") and we report the EFFECT SIZE.
  matched_effect <- function(idx, x, cov, universe = rep(TRUE, length(x)), nperm = 2000, seed = 42) {
    ok <- is.finite(cov) & is.finite(x) & universe
    idx <- idx[ok[idx]]
    if (length(idx) < 10 || sum(ok) < 50) return(NULL)
    set.seed(seed); cj <- cov[ok] + runif(sum(ok), 0, 1e-6)
    qi <- cumsum(ok)[idx]; xx <- x[ok]
    hi <- size_perm(mean(xx[qi]),  cj[qi], cj,  xx, nperm = nperm, seed = seed)
    lo <- size_perm(mean(-xx[qi]), cj[qi], cj, -xx, nperm = nperm, seed = seed)
    data.table(n_bins = length(idx), obs_mean = hi$obs_frac, matched_mean = hi$exp_frac,
               diff = hi$obs_frac - hi$exp_frac, z = hi$z,
               perm_p = min(1, 2 * min(hi$perm_p, lo$perm_p)))
  }

  # ================================================================
  # 1. Contact LOSS / GAIN density per group
  # ================================================================
  # Direction (config.R): selWT = stronger in WT = LOST in HP1gKO -> n_contacts_lost;
  # selcKO = GAINED in HP1gKO -> n_contacts_gained. Columns are KO-centric.
  # The hotspot counts are DENSE (almost every mappable bin has a value, 1-162), so the
  # signal is count MAGNITUDE; compare distributions, not "has any".
  metrics <- intersect(c("n_contacts_lost", "n_contacts_gained"), names(d))
  summ <- rbindlist(lapply(names(schemes), function(sc) {
    x <- copy(d); x[, grp := factor(get(schemes[[sc]]), levels = grp_levels)]
    x <- x[!is.na(grp)]
    out <- x[, c(list(n_bins = .N,
                      median_lost = as.numeric(median(n_contacts_lost, na.rm = TRUE)),
                      mean_lost   = mean(n_contacts_lost,   na.rm = TRUE),
                      median_gained = as.numeric(median(n_contacts_gained, na.rm = TRUE)),
                      mean_gained   = mean(n_contacts_gained,   na.rm = TRUE),
                      median_H4K20me3_l2fc = median(h4_l2fc, na.rm = TRUE),
                      frac_tad_boundary = mean(tad_boundary, na.rm = TRUE),
                      mean_repeat_content = mean(rep_cov, na.rm = TRUE)),
                 if (have_feat) list(median_pc1_delta = median(pc1_delta, na.rm = TRUE),
                                     median_insulation_delta = median(insulation_delta, na.rm = TRUE))),
             by = grp]
    out[, scheme := sc][]
  }), fill = TRUE)

  # Each group against a repeat-content-matched draw of ALL labelled-or-not bins (not
  # against `stable`: `stable` is itself a small, repeat-poor set and a stable-vs-group
  # contrast would be confounded by that). Effect sizes + permutation p, see
  # matched_effect().
  vs_matched <- rbindlist(lapply(names(schemes), function(sc) {
    g <- factor(d[[schemes[[sc]]]], levels = grp_levels)
    mets <- c(metrics, "h4_l2fc", if (have_feat) c("pc1_delta", "insulation_delta"))
    rbindlist(lapply(grp_levels, function(gp) {
      rbindlist(lapply(mets, function(m) {
        r <- matched_effect(which(g == gp), d[[m]], d$rep_cov)
        if (is.null(r)) NULL else cbind(data.table(scheme = sc, group = gp, metric = m), r)
      }))
    }))
  }))
  vs_matched[, p_adj_BH := p.adjust(perm_p, "BH"), by = scheme]

  # ================================================================
  # 2. size_perm(): overlap-shaped questions, matched on repeat content
  # ================================================================
  # "Is the group enriched for the top decile of contact-loss bins?" is an overlap
  # question. A naive hypergeometric/Fisher test would compare against ALL bins, which
  # differ systematically in repeat density (see matched_perm above). The null here
  # draws bins with the same repeat-content distribution as the group.
  d[, top_lost   := n_contacts_lost   >= quantile(n_contacts_lost,   0.9, na.rm = TRUE)]
  d[, top_gained := n_contacts_gained >= quantile(n_contacts_gained, 0.9, na.rm = TRUE)]
  flags <- c("top_lost", "top_gained", "tad_boundary")
  if (have_feat) {
    d[, switch_AB := compartment_switch == "A_B"]   # WT A -> KO B
    d[, switch_BA := compartment_switch == "B_A"]   # WT B -> KO A
    flags <- c(flags, "switch_AB", "switch_BA")
  }
  perm_res <- rbindlist(lapply(names(schemes), function(sc) {
    g <- factor(d[[schemes[[sc]]]], levels = grp_levels)
    rbindlist(lapply(grp_levels, function(gp) {
      idx <- which(g == gp)
      rbindlist(lapply(flags, function(fl) {
        r <- matched_perm(idx, d[[fl]], d$rep_cov)
        if (is.null(r)) return(NULL)
        cbind(data.table(scheme = sc, group = gp, flag = fl, n_bins = length(idx)), r)
      }))
    }))
  }))
  perm_res[, p_adj_BH := p.adjust(perm_p, "BH"), by = scheme]

  # ================================================================
  # 3. A/B switching: does the mark loss move chromatin between compartments?
  # ================================================================
  switch_tab <- NULL; cor_tab <- NULL
  if (have_feat) {
    # Report BOTH counts and the within-group fraction. Counts alone are dominated by
    # group size (stable >> co_loss); fractions are what answer "does this group switch
    # more than the others". The genome-wide row is the baseline to compare against.
    switch_tab <- rbindlist(lapply(names(schemes), function(sc) {
      x <- copy(d); x[, grp := factor(get(schemes[[sc]]), levels = grp_levels)]
      x <- x[!is.na(compartment_switch)]
      one <- function(xx, lab) {
        ct <- xx[, .N, by = compartment_switch]
        ct <- merge(data.table(compartment_switch = levels(d$compartment_switch)), ct,
                    by = "compartment_switch", all.x = TRUE)
        ct[is.na(N), N := 0L]
        ct[, `:=`(group = lab, group_total = sum(N), fraction = N / sum(N))][]
      }
      out <- rbind(rbindlist(lapply(grp_levels, function(gp) one(x[grp == gp], gp))),
                   one(x, "all_bins"))
      out[, scheme := sc][]
    }))
    setcolorder(switch_tab, c("scheme", "group", "compartment_switch", "N", "group_total", "fraction"))

    # PC1 shift vs H4K20me3 log2FC, genome-wide, per bin. Spearman: both are
    # heavy-tailed and we want monotonic association, not linearity. If mark loss
    # moves chromatin A->B, deeper H4K20me3 loss should go with a more NEGATIVE
    # pc1_delta (A positive), i.e. a POSITIVE rho.
    ok <- is.finite(d$pc1_delta) & is.finite(d$h4_l2fc)
    ct <- suppressWarnings(cor.test(d$pc1_delta[ok], d$h4_l2fc[ok], method = "spearman"))
    # NO P-VALUE. cor.test treats the bins as independent draws; they are not -
    # the differential-contact track has autocorrelation r = 0.59 at 20 kb, and
    # GAM/CLAUDE.md forbids exactly this. The reported p was 8.4e-27 on 117,274
    # bins for rho = -0.031, which is a statement about n, not about an effect.
    # The effect size is the number; the chromosome bootstrap elsewhere in this
    # file is what carries the uncertainty.
    cor_tab <- data.table(n_bins = sum(ok), spearman_rho = unname(ct$estimate))
    # also within WT compartment: a genome-wide rho can be carried by A-vs-B baseline
    # differences in both variables, so stratify
    if ("compartment_WT" %in% names(d))
      cor_tab <- rbind(cor_tab, rbindlist(lapply(c("A", "B"), function(cp) {
        k <- ok & d$compartment_WT %in% cp
        cc <- suppressWarnings(cor.test(d$pc1_delta[k], d$h4_l2fc[k], method = "spearman"))
        data.table(n_bins = sum(k), spearman_rho = unname(cc$estimate),
                   stratum = paste0("WT_", cp))   # no p: see the genome_wide row
      })), fill = TRUE)
    cor_tab[is.na(stratum), stratum := "genome_wide"]
  }

  # ================================================================
  # 4. Insulation change per group, and at TAD boundaries
  # ================================================================
  # insulation_delta = HP1gKO - WT (GAM_2; interpret the sign per GAM_2's definition of
  # the insulation score). Boundaries are where a change matters most, so we split each
  # group by tad_boundary; the all-bins row is the baseline a boundary should be read
  # against, because boundaries differ from non-boundary bins in insulation anyway.
  ins_tab <- NULL
  if (have_feat && "insulation_delta" %in% names(d)) {
    ins_tab <- rbindlist(lapply(names(schemes), function(sc) {
      x <- copy(d); x[, grp := factor(get(schemes[[sc]]), levels = grp_levels)]
      x <- x[is.finite(insulation_delta) & !is.na(tad_boundary)]
      rbind(x[!is.na(grp), .(n = .N, median_insulation_delta = median(insulation_delta),
                             mean_insulation_delta = mean(insulation_delta)),
              by = .(group = as.character(grp), tad_boundary)],
            x[, .(group = "all_bins", n = .N, median_insulation_delta = median(insulation_delta),
                  mean_insulation_delta = mean(insulation_delta)), by = tad_boundary])[
              , scheme := sc][]
    }))
    # group's boundary bins vs repeat-content-matched BOUNDARY bins (universe = all
    # boundary bins), and likewise for non-boundary bins: does the group change its
    # insulation more than boundaries/non-boundaries of the same repeat content do?
    ins_bd <- rbindlist(lapply(names(schemes), function(sc) {
      g <- factor(d[[schemes[[sc]]]], levels = grp_levels)
      rbindlist(lapply(grp_levels, function(gp) {
        rbindlist(lapply(c(TRUE, FALSE), function(bd) {
          r <- matched_effect(which(g == gp & d$tad_boundary %in% bd), d$insulation_delta,
                              d$rep_cov, universe = d$tad_boundary %in% bd)
          if (is.null(r)) NULL else cbind(data.table(scheme = sc, group = gp, tad_boundary = bd), r)
        }))
      }))
    }))
    ins_bd[, p_adj_BH := p.adjust(perm_p, "BH"), by = scheme]
    fwrite(ins_bd, file.path(tables_dir, "gam4_insulation_boundary_vs_matched.tsv"), sep = "\t")
  }

  # ================================================================
  # 5. Tables
  # ================================================================
  fwrite(summ,      file.path(tables_dir, "gam4_group_summary.tsv"),        sep = "\t")
  fwrite(vs_matched, file.path(tables_dir, "gam4_group_vs_matched_background.tsv"), sep = "\t")
  fwrite(perm_res,  file.path(tables_dir, "gam4_size_perm_overlap.tsv"),    sep = "\t")
  if (!is.null(switch_tab)) fwrite(switch_tab, file.path(tables_dir, "gam4_compartment_switch.tsv"), sep = "\t")
  if (!is.null(cor_tab))    fwrite(cor_tab,    file.path(tables_dir, "gam4_pc1_vs_H4K20me3.tsv"),   sep = "\t")
  if (!is.null(ins_tab))    fwrite(ins_tab,    file.path(tables_dir, "gam4_insulation_by_group.tsv"), sep = "\t")

  cat("\n=== GAM_4 per-group summary (", PRIMARY, " = primary; abs = sensitivity) ===\n", sep = "")
  print(summ, row.names = FALSE)
  cat("\n=== size-matched overlap (matched on repeat content) ===\n"); print(perm_res, row.names = FALSE)
  if (!is.null(switch_tab)) { cat("\n=== A/B switch table ===\n"); print(switch_tab, row.names = FALSE) }
  if (!is.null(cor_tab))    { cat("\n=== PC1 shift vs H4K20me3 log2FC ===\n"); print(cor_tab, row.names = FALSE) }

  # ================================================================
  # Figures -> figures/compartments_vs_contacts
  # ================================================================
  sub_dir <- "compartments_vs_contacts"
  long_c <- rbindlist(lapply(names(schemes), function(sc) {
    x <- d[, .(bin_id, grp = factor(get(schemes[[sc]]), levels = grp_levels),
               lost = n_contacts_lost, gained = n_contacts_gained)]
    x <- melt(x[!is.na(grp)], id.vars = c("bin_id", "grp"), variable.name = "direction",
              value.name = "n")
    x[, scheme := factor(sc, levels = names(schemes))][]
  }))
  long_c[, direction := factor(direction, levels = c("lost", "gained"),
                               labels = c("lost in HP1gKO (selWT)", "gained in HP1gKO (selcKO)"))]
  g_contacts <- ggplot(long_c, aes(grp, log1p(n), fill = grp)) +
    geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.15, linewidth = 0.3) +
    facet_grid(scheme ~ direction, labeller = labeller(scheme = c(
      centred = "median-centred (primary)", abs = "absolute cutoffs (sensitivity)"))) +
    scale_fill_manual(values = grp_cols, guide = "none") +
    labs(title = "Differential contact density per chromatin-mark-loss group",
         subtitle = "per 20 kb bin; the absolute-cutoff row is a sensitivity check, not the primary result",
         x = NULL, y = "log1p(number of differential contacts)") +
    theme_m
  save_fig(g_contacts, "contact_density_by_group", sub_dir, width = 9, height = 6)

  pp <- copy(perm_res)[, `:=`(group = factor(group, levels = grp_levels),
                              scheme = factor(scheme, levels = names(schemes)),
                              sig = ifelse(p_adj_BH < 0.05, "*", ""))]
  g_perm <- ggplot(pp, aes(flag, log2(fold), fill = group)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_text(aes(label = sig), position = position_dodge(0.8), vjust = -0.3) +
    facet_wrap(~scheme, ncol = 1) +
    scale_fill_manual(values = grp_cols, name = NULL) +
    labs(title = "Overlap vs a repeat-content-matched null (size_perm)",
         subtitle = "log2(observed / expected); * = BH < 0.05; top_* = top decile of contact counts",
         x = NULL, y = "log2 fold over matched null") +
    theme_m + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_fig(g_perm, "size_perm_overlap", sub_dir, width = 8, height = 6)

  if (!is.null(switch_tab)) {
    st <- copy(switch_tab)[, `:=`(group = factor(group, levels = c(grp_levels, "all_bins")),
                                  scheme = factor(scheme, levels = names(schemes)))]
    g_sw <- ggplot(st, aes(group, fraction, fill = compartment_switch)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = ifelse(fraction > 0.04, N, "")),
                position = position_stack(vjust = 0.5), size = 3) +
      facet_wrap(~scheme, nrow = 1) +
      scale_fill_disc(name = "WT -> HP1gKO") +
      labs(title = "A/B compartment switching per group",
           subtitle = "bar = fraction of the group's bins; label = bin count (A positive = H3K4me3-high)",
           x = NULL, y = "fraction of group") +
      theme_m + theme(axis.text.x = element_text(angle = 30, hjust = 1))
    save_fig(g_sw, "compartment_switch_fraction", sub_dir, width = 10, height = 4.5)

    dd <- d[is.finite(pc1_delta) & is.finite(h4_l2fc)]
    g_sc <- ggplot(dd, aes(h4_l2fc, pc1_delta)) +
      geom_bin2d(bins = 80) + scale_fill_heat0(trans = "log10", name = "bins") +
      geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
      geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55") +
      labs(title = "PC1 shift vs H4K20me3 change, per bin",
           subtitle = sprintf("Spearman rho = %.3f (n = %s); a positive rho = deeper loss goes with a shift toward B",
                              cor_tab$spearman_rho[cor_tab$stratum == "genome_wide"],
                              format(cor_tab$n_bins[cor_tab$stratum == "genome_wide"], big.mark = ",")),
           x = "H4K20me3 log2FC (HP1gKO / WT)", y = "PC1 delta (HP1gKO - WT)") +
      theme_m
    save_fig(g_sc, "pc1_delta_vs_H4K20me3", sub_dir, width = 7, height = 5.5)

    g_pc <- ggplot(d[!is.na(get(schemes[[PRIMARY]])) & is.finite(pc1_delta)],
                   aes(factor(get(schemes[[PRIMARY]]), levels = grp_levels), pc1_delta,
                       fill = factor(get(schemes[[PRIMARY]]), levels = grp_levels))) +
      geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
      geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.15, linewidth = 0.3) +
      scale_fill_manual(values = grp_cols, guide = "none") +
      labs(title = "PC1 shift per group (median-centred groups)",
           subtitle = "delta = HP1gKO - WT; positive = toward A", x = NULL, y = "PC1 delta") +
      theme_m
    save_fig(g_pc, "pc1_delta_by_group", sub_dir, width = 5.5, height = 4.5)
  }

  if (!is.null(ins_tab)) {
    it <- d[!is.na(get(schemes[[PRIMARY]])) & is.finite(insulation_delta) & !is.na(tad_boundary)]
    it[, `:=`(grp = factor(get(schemes[[PRIMARY]]), levels = grp_levels),
              where = factor(ifelse(tad_boundary, "TAD boundary", "other bins"),
                             levels = c("other bins", "TAD boundary")))]
    g_ins <- ggplot(it, aes(grp, insulation_delta, fill = where)) +
      geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
      geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.15, linewidth = 0.3,
                   position = position_dodge(0.8)) +
      scale_fill_disc(name = NULL) +
      labs(title = "Insulation change per group, split by TAD boundary",
           subtitle = "delta = HP1gKO - WT (median-centred groups)", x = NULL,
           y = "insulation delta") +
      theme_m
    save_fig(g_ins, "insulation_delta_by_group", sub_dir, width = 7, height = 4.5)
  }

  message("GAM_4 done: tables in ", tables_dir, ", figures in ", file.path(fig_dir, sub_dir))
})


# ==============================================================================
# SECTION D - the synthesis figure
# ==============================================================================
# The three sections above answer one question three ways and, until they were
# merged, disagreed only in presentation. This panel says what each one actually
# found, in plain language, on one page - so a reader does not have to
# reconstruct it from three sets of axis labels and a machine-generated verdict
# string.
#
# Deliberately NOT a new analysis. It reads the tables the sections just wrote.
local({
  f <- function(x) file.path(tables_dir, x)
  need <- c("coupling_specificity.tsv", "gam7_coupling_by_distance.tsv",
            "gam10_class_change.tsv")
  if (!all(file.exists(f(need)))) {
    message("synthesis figure skipped - needs: ", paste(need, collapse = ", "))
    return(invisible(NULL))
  }
  FLAT <- 0.05

  # --- panel 1: per-bin correlation, all five marks -------------------------
  sp <- fread(f("coupling_specificity.tsv"))[variant == "raw"]
  sp[, mark := factor(mark, levels = mark_levels)]
  p1 <- ggplot(sp, aes(rho_lost, mark)) +
    annotate("rect", xmin = -FLAT, xmax = FLAT, ymin = -Inf, ymax = Inf,
             fill = "grey90") +
    geom_vline(xintercept = 0, linewidth = .3, colour = "grey50") +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y",
                  width = .18, colour = "grey35") +
    geom_point(size = 2, colour = "#b2182b") +
    # The x range must be wider than FLAT or the "too small to matter" band
    # fills the panel and stops reading as a band.
    coord_cartesian(xlim = c(-0.13, 0.13)) +
    labs(subtitle = "1. Does a bin that loses more mark\n    also change its contacts more?",
         x = "correlation across 20 kb bins", y = NULL,
         caption = "grey band = too small to matter") +
    theme_m + theme(plot.caption = element_text(hjust = 0, size = 7, colour = "grey40"))

  # --- panel 2: the same question recounted locally, by distance ------------
  g7 <- fread(f("gam7_coupling_by_distance.tsv"))[mark == "H4K20me3" & direction == "lost"]
  g7 <- g7[band != "0-5Mb window"]
  g7[, band := factor(band, levels = c("all (local)", "0.02-1Mb", "1-5Mb", "5-10Mb", ">10Mb"))]
  p2 <- ggplot(g7, aes(rho, band)) +
    annotate("rect", xmin = -FLAT, xmax = FLAT, ymin = -Inf, ymax = Inf, fill = "grey90") +
    geom_vline(xintercept = 0, linewidth = .3, colour = "grey50") +
    geom_point(size = 2, colour = "#b2182b") +
    coord_cartesian(xlim = c(-0.13, 0.13)) +
    labs(subtitle = "2. Same question, but counting contacts\n    locally instead of over 5 Mb",
         x = "correlation across 20 kb bins", y = NULL,
         caption = "H4K20me3; the 5 Mb window was a ~250x scale mismatch") +
    theme_m + theme(plot.caption = element_text(hjust = 0, size = 7, colour = "grey40"))

  # --- panel 3: where the effect actually is --------------------------------
  cc <- fread(f("gam10_class_change.tsv"))[cls %in% c("inter_het", "eu_eu")]
  cc[, cls := factor(fifelse(cls == "eu_eu", "euchromatin", "heterochromatin"),
                     levels = c("heterochromatin", "euchromatin"))]
  ord <- c("<0.2Mb", "0.2-1Mb", "1-5Mb", "5-10Mb", ">10Mb")
  cc <- cc[band %in% ord][, band := factor(band, levels = ord)]
  p3 <- ggplot(cc, aes(band, pct, colour = cls, group = cls)) +
    geom_hline(yintercept = 0, linewidth = .3, colour = "grey50") +
    geom_line() + geom_point(size = 1.8) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = .15) +
    scale_colour_manual(values = c(heterochromatin = "#b2182b", euchromatin = "#2166ac"),
                        name = NULL) +
    labs(subtitle = "3. Forget the marks' magnitude - just ask which\n    kind of chromatin changed",
         x = "distance between the two bins", y = "change in contact strength (%)",
         caption = "bars = 95% CI from a bootstrap over the 19 autosomes") +
    theme_m + theme(legend.position = "bottom",
                    plot.caption = element_text(hjust = 0, size = 7, colour = "grey40"))

  pp <- p1 + p2 + p3 + patchwork::plot_layout(widths = c(1, 1, 1.25)) +
    patchwork::plot_annotation(
      title = "Does the H4K20me3 loss have a contact counterpart? Three ways of asking",
      subtitle = paste0(
        "Panels 1 and 2 find nothing, and that is a real negative, not a weak one: n = 119,952 bins, and the extremes make it concrete - ",
        "bins losing the most\nH4K20me3 (median -0.76) and the least (-0.22) have identical contact densities. Panel 2 exists because panel 1's contact measure is ",
        "integrated over\n5 Mb while the mark varies at 20 kb, so its null could have been a resolution limit; recounted locally, the null holds at every distance.\n\n",
        "Panel 3 is where the signal is, and it is a different question: not 'how much mark did this bin lose' but 'what kind of chromatin is this'. ",
        "Heterochromatin\ngains contacts while euchromatin loses them. So the architecture does change - it just does not track the size of the mark loss bin by bin."),
      theme = theme(plot.title = element_text(face = "bold", size = 13),
                    plot.subtitle = element_text(size = 8.5, colour = "grey25")))
  save_fig(pp, "coupling_synthesis", subdir = "coupling", width = 13.5, height = 6.4)
  cat("\nSaved the synthesis figure: figures/coupling/coupling_synthesis\n")
})
