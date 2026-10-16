# ================================================================
# GAM_3 - Contact coupling: does 3D contact remodelling coincide with the
#         H4K20me3 loss?   (the PRIMARY question of the GAM integration)
# ----------------------------------------------------------------
# Joins, on bin_id, the per-bin ChIP log2FC (GAM_1) to the per-bin differential
# contact density (GAM_1 $annot) and, when present, the NPMI-derived features
# (GAM_2). Sections:
#   1. Coupling            Spearman + binned-quantile curves + 2D density
#   2. Direction asymmetry lost and gained treated SEPARATELY, plus a net measure
#   3. SPECIFICITY CONTROL the falsification test (active marks must be flat)
#   4. Confounder strata   ChromHMM state / repeat content / coverage, + partials
#   5. Compartment view    contact density per compartment_centred group
#   6. NPMI features       pc1/insulation/npmi/decay (only if GAM_2 has run)
#
# DIRECTION SEMANTICS (config.R): selWT = stronger in WT = LOST in the KO;
# selcKO = GAINED. n_contacts_lost / n_contacts_gained are already named from the
# KO's point of view, so "lost" below always means "present in WT, weaker in KO".
#
# Runs after GAM_1 (required) and GAM_2 (optional). Run from GAM/ .
# ================================================================
source("config.R")

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

A[, net := log2((n_contacts_lost + 1) / (n_contacts_gained + 1))]
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
                    "net = log2((lost+1)/(gained+1)); positive = more contacts lost than gained",
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
lim <- max(abs(sp$rho), na.rm = TRUE)
g_strat <- ggplot(sp, aes(mark, lab, fill = rho)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%+.2f", rho)), size = 2.4) +
  facet_wrap(~ direction, nrow = 1) +
  scale_fill_heat0_div(limits = c(-lim, lim), name = "Spearman\nrho") +
  labs(title = "Coupling within confounder strata", x = NULL, y = NULL,
       subtitle = "rho of per-bin log2FC vs contact density; strata dropped below the minimum size") +
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
