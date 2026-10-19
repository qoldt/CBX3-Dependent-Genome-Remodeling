# ================================================================
# QC — hotspot domain calling: why these parameters
# ----------------------------------------------------------------
# Standalone QC. Reads GAM_1's rds, writes one multi-panel figure plus a
# parameter-sweep table. Run after GAM_1:
#     cd GAM && Rscript QC_Hotspot_Domains.R
#
# WHAT QUESTION THIS ANSWERS. The delivered Hotspots files give, per 20 kb bin,
# the number of significant differential contacts within 0-5 Mb of that bin. The
# track is DENSE (93% of autosomal bins carry at least one) and SMOOTH, so
# neighbouring bins are far from independent. That creates a choice: analyse the
# per-bin count as a continuous covariate, or threshold it and merge adjacent
# bins into domains and analyse those?
#
# Both are legitimate and they answer different questions, so the pipeline keeps
# both. This script exists to show a reader - including a collaborating lab that
# did not choose these parameters - what the alternatives looked like and why the
# defaults in GAM/config.R (gam_hotspot_params) were chosen rather than asserted.
#
# The three panels, and the argument each one carries:
#   A  Autocorrelation vs genomic lag. Sets how wide a smoothing window is
#      defensible: averaging is justified while neighbouring bins still share
#      signal, and becomes fusion of distinct domains once they do not.
#   B  Gap-tolerant merging sweep. Shows that absorbing sub-threshold gaps does
#      consolidate fragments, but chains into implausibly long domains if pushed.
#   C  Smooth-then-threshold sweep. Shows the same consolidation achieved without
#      the chaining pathology, which is why this is the default.
#   D  Run-length distributions: raw vs the chosen parameters.
#   E  A worked genomic example - raw counts, the smoothed track, the threshold
#      and the called domains on real coordinates.
# ================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(patchwork) })

rds <- file.path(rds_dir, "gam_bins_annotated.rds")
if (!file.exists(rds))
  stop("GAM_1 output not found:\n  ", rds,
       "\nRun it first:  cd GAM && Rscript GAM_1_Bins_and_Signal.R", call. = FALSE)
x <- readRDS(rds)
a <- as.data.table(x$annot); setorder(a, chr, start)
cat("Loaded", nrow(a), "bins\n")

# Direction under study. "lost" (selWT = stronger in WT = lost in the KO) is the
# one the headline question is about; the same sweep for "gained" is written to
# the table so the choice is not silently direction-specific.
DIRS <- c(lost = "n_contacts_lost", gained = "n_contacts_gained")

# ---- helper: run lengths with optional gap absorption -----------------------
run_lengths <- function(hot, chr, gap = 0L) {
  data.table(chr = chr, hot = hot)[, {
    h <- hot
    if (gap > 0L && length(h) > 2L) {
      r <- rle(h); L <- r$lengths; V <- r$values; k <- 2L
      while (k < length(V)) {
        if (!V[k] && L[k] <= gap && V[k - 1L] && V[k + 1L]) {
          L[k - 1L] <- L[k - 1L] + L[k] + L[k + 1L]
          V <- V[-c(k, k + 1L)]; L <- L[-c(k, k + 1L)]
        } else k <- k + 1L
      }
      h <- rep(V, L)
    }
    r2 <- rle(h); .(len = r2$lengths[r2$values])
  }, by = chr]$len
}
smooth_within_chr <- function(v, chr, w) {
  if (w <= 1L) return(as.numeric(v))
  data.table(chr = chr, v = as.numeric(v))[, {
    s <- as.numeric(stats::filter(v, rep(1 / w, w), sides = 2))
    .(s = ifelse(is.na(s), v, s))
  }, by = chr]$s
}
summarise_runs <- function(L, label, strategy) {
  if (!length(L)) L <- 0L
  data.table(strategy = strategy, setting = label, n_domains = length(L),
             median_bins = as.numeric(median(L)), mean_bins = mean(L),
             max_bins = max(L), max_Mb = max(L) * gam_binsize / 1e6,
             frac_bins_in_ge5 = sum(L[L >= 5]) / max(1, sum(L)))
}

# ---- A: autocorrelation, per chromosome then averaged ----------------------
# Averaging the per-chromosome ACFs rather than computing one genome-wide ACF
# avoids treating a chromosome boundary as a 20 kb step.
LAGMAX <- 30L
acf_dt <- rbindlist(lapply(names(DIRS), function(d) {
  m <- a[, .(ac = list(stats::acf(get(DIRS[[d]]), lag.max = LAGMAX, plot = FALSE)$acf[, , 1])), by = chr]
  data.table(direction = d, lag = 0:LAGMAX, r = colMeans(do.call(rbind, m$ac)))
}))
acf_dt[, kb := lag * gam_binsize / 1000]
chosen_kb <- gam_hotspot_params$smooth_bins * gam_binsize / 1000
pA <- ggplot(acf_dt[lag > 0], aes(kb, r, colour = direction)) +
  geom_hline(yintercept = 0, linewidth = .3, colour = "grey70") +
  geom_vline(xintercept = chosen_kb, linetype = 2, colour = "grey40") +
  geom_line(linewidth = .7) + geom_point(size = 1.1) +
  annotate("text", x = chosen_kb, y = max(acf_dt$r[acf_dt$lag > 0]),
           label = paste0("  chosen window ", chosen_kb, " kb"), hjust = 0, size = 3, colour = "grey30") +
  scale_colour_disc(name = "contacts") +
  labs(title = "A. Neighbouring bins share signal",
       subtitle = "Per-chromosome autocorrelation of the hotspot count track, averaged over chromosomes.\nSmoothing is defensible while r is substantial; beyond that it fuses distinct domains.",
       x = "genomic lag (kb)", y = "autocorrelation r") + theme_m

# ---- B/C: the two consolidation strategies ---------------------------------
q <- gam_hotspot_params$q
sweep <- rbindlist(lapply(names(DIRS), function(d) {
  v <- a[[DIRS[[d]]]]
  thr <- quantile(v, q, na.rm = TRUE)
  gaps <- rbindlist(lapply(0:4, function(g)
    cbind(direction = d, summarise_runs(run_lengths(v >= thr, a$chr, g),
          paste0("gap<=", g), "gap-merge raw track"))))
  smos <- rbindlist(lapply(c(1L, 3L, 5L, 7L, 11L, 15L), function(w) {
    s <- smooth_within_chr(v, a$chr, w)
    cbind(direction = d, summarise_runs(run_lengths(s >= quantile(s, q, na.rm = TRUE), a$chr, 0L),
          paste0(w * gam_binsize / 1000, " kb"), "smooth then threshold"))
  }))
  rbind(gaps, smos)
}))
fwrite(sweep, file.path(tables_dir, "qc_hotspot_domain_parameter_sweep.tsv"), sep = "\t")
cat("Wrote", file.path(tables_dir, "qc_hotspot_domain_parameter_sweep.tsv"), "\n")

mk_sweep_panel <- function(strat, title, sub, chosen) {
  d <- sweep[strategy == strat]
  d[, setting := factor(setting, levels = unique(setting))]
  # Facet by direction rather than dodge: the med/max annotations are the point of
  # this panel and they collide illegibly between dodged bars.
  ggplot(d, aes(setting, n_domains, fill = direction)) +
    geom_col(width = .7, show.legend = FALSE) +
    geom_text(aes(label = sprintf("med %.0f bins\nmax %.1f Mb", median_bins, max_Mb)),
              vjust = -0.25, size = 2.6, lineheight = .95, colour = "grey25") +
    facet_wrap(~ direction, ncol = 1, scales = "free_y") +
    scale_y_continuous(expand = expansion(mult = c(0, .3))) +
    scale_fill_disc() +
    labs(title = title, subtitle = sub, x = NULL, y = "domains called") +
    theme_m + theme(axis.text.x = element_text(
      face = ifelse(levels(d$setting) == chosen, "bold", "plain")))
}
pB <- mk_sweep_panel("gap-merge raw track", "B. Gap-merging the raw track",
  paste0("Threshold at q", q * 100, ", then absorb sub-threshold gaps. Consolidates fragments,\n",
         "but chains into multi-Mb domains if pushed, and the gap size is arbitrary."), "gap<=1")
pC <- mk_sweep_panel("smooth then threshold", "C. Smoothing first (the default)",
  paste0("Running mean within chromosome, then threshold at q", q * 100, ".\n",
         "Same consolidation, no chaining; window set by panel A."),
  paste0(chosen_kb, " kb"))

# ---- D: run-length distributions, raw vs chosen -----------------------------
v <- a[[DIRS[["lost"]]]]
raw_L <- run_lengths(v >= quantile(v, q, na.rm = TRUE), a$chr, 0L)
dom   <- gam_hotspot_domains(v, a$chr, a$start, a$end)
cho_L <- dom$domains$n_bins
pD <- ggplot(rbind(data.table(set = "raw track, strict adjacency", len = raw_L),
                   data.table(set = "chosen parameters", len = cho_L)),
             aes(len, fill = set)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = .65) +
  geom_vline(xintercept = gam_hotspot_params$min_bins, linetype = 2, colour = "grey40") +
  # log y: the raw track's singleton spike (~4,800) is 2 orders of magnitude above
  # the chosen-parameter distribution, which would otherwise be invisible here.
  scale_y_log10() + scale_x_continuous(limits = c(0, 40)) + scale_fill_disc(name = NULL) +
  labs(title = "D. Domain size distribution (contacts lost)",
       subtitle = paste0("Dashed = the ", gam_hotspot_params$min_bins, "-bin (",
                         gam_hotspot_params$min_bins * gam_binsize / 1000, " kb) size floor.\n",
                         "Raw thresholding is dominated by singletons."),
       x = "domain width (20 kb bins)", y = "domains (log scale)") + theme_m

# ---- E: a worked example on real coordinates -------------------------------
# The region with the most called domain territory, so the panel shows the
# caller working rather than an empty stretch.
best <- dom$domains[, .(bp = sum(end - start)), by = chr][order(-bp)][1, chr]
reg  <- a[chr == best][order(start)]
ctr  <- dom$domains[chr == best][which.max(n_bins)]
w0   <- max(min(reg$start), ctr$start - 2.0e6); w1 <- min(max(reg$end), ctr$end + 2.0e6)
sub  <- reg[start >= w0 & end <= w1]
sub[, smoothed := smooth_within_chr(get(DIRS[["lost"]]), chr, gam_hotspot_params$smooth_bins)]
dsub <- dom$domains[chr == best & end > w0 & start < w1]
pE <- ggplot(sub, aes(start / 1e6)) +
  { if (nrow(dsub)) geom_rect(data = dsub, inherit.aes = FALSE,
      aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf),
      fill = gaby_cols[3], alpha = .45) } +
  geom_col(aes(y = get(DIRS[["lost"]])), fill = "grey70", width = gam_binsize / 1e6) +
  geom_line(aes(y = smoothed), colour = hat_cols[2], linewidth = .7) +
  geom_hline(yintercept = dom$params$cutoff, linetype = 2, colour = "grey30") +
  labs(title = paste0("E. Worked example: ", best, " ",
                      round(w0 / 1e6, 1), "-", round(w1 / 1e6, 1), " Mb"),
       subtitle = paste0("Grey = raw per-bin count; red = ", chosen_kb,
                         " kb running mean; dashed = q", q * 100, " threshold (",
                         round(dom$params$cutoff, 1), "); shaded = called domains."),
       x = paste0("position on ", best, " (Mb)"), y = "differential contacts lost") + theme_m

fig <- (pA | pB) / (pC | pD) / pE +
  patchwork::plot_annotation(
    title = "Hotspot domain calling QC - GAM/ChIP integration",
    # What the sweep is for, beyond the parameters themselves: the per-bin count
    # track is retained UNMODIFIED as the primary continuous covariate, and
    # domains are the independent unit for enrichment tests, because adjacent
    # bins are not independent (panel A). Both an adjusted and an UNADJUSTED
    # domain set are carried through the pipeline.
    subtitle = paste0("Defaults in GAM/config.R$gam_hotspot_params: q=", q,
      ", smooth=", gam_hotspot_params$smooth_bins, " bins (", chosen_kb, " kb), gap<=",
      gam_hotspot_params$gap, ", min=", gam_hotspot_params$min_bins, " bins."),
    theme = theme(plot.title = element_text(face = "bold", size = 14)))
save_fig(fig, "hotspot_domain_calling_qc", subdir = "qc", width = 15, height = 15)

cat("\n=== chosen parameters, contacts lost ===\n")
cat(sprintf("  threshold (q%.0f of smoothed track): %.1f\n", q * 100, dom$params$cutoff))
cat(sprintf("  domains kept (>=%d bins): %d\n", gam_hotspot_params$min_bins, nrow(dom$domains)))
cat(sprintf("  median %.0f bins (%.0f kb); max %.1f Mb\n", median(cho_L),
            median(cho_L) * gam_binsize / 1000, max(cho_L) * gam_binsize / 1e6))
cat(sprintf("  genome in called domains: %.1f%%\n",
            100 * sum(dom$domains$end - dom$domains$start) / sum(a$end - a$start)))
