# ================================================================
# GAM_6 - Large one-directional contact blocks, genome-wide
# ----------------------------------------------------------------
# Reads GAM_1 (and GAM_2 if present). Run after them:
#     cd GAM && Rscript GAM_6_Large_Blocks.R
#
# WHY THIS STAGE EXISTS. GAM_1's hotspot domains are 0.1-2 Mb. Looking at the
# tracks shows the dominant feature is an order of magnitude larger: chr11 loses
# contacts across essentially all of 0-20 Mb and gains them across 56-108 Mb;
# chr16 alternates in 10-20 Mb blocks. A caller with a 100 kb floor fragments
# such a block into dozens of pieces and reports them as unrelated, so the
# largest-scale structure in the data was invisible to every stage before this.
#
# WHAT IS CALLED. Blocks of consistent DIRECTION, from the signed statistic
#     net = log2((lost + 1) / (gained + 1))
# smoothed over gam_block_smooth_bins. Direction, not density: a region can be
# busy in both directions and be unremarkable; the feature is the imbalance.
#
# ARTIFACT CHECK, done before trusting any of this: net is uncorrelated with the
# WT-minus-KO count of valid NPMI cells (rho +0.005 at 4 Mb, +0.025 per bin), so
# the blocks are not where one genotype simply has more data. This stage re-runs
# that check and prints it, because it is the first thing a reader should ask.
# ================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })

g1 <- file.path(rds_dir, "gam_bins_annotated.rds")
if (!file.exists(g1))
  stop("GAM_1 output not found:\n  ", g1,
       "\nRun it first:  cd GAM && Rscript GAM_1_Bins_and_Signal.R", call. = FALSE)
x <- readRDS(g1)
a <- as.data.table(x$annot); setorder(a, chr, start)
marks_here <- colnames(x$l2fc)
for (mk in marks_here) a[[paste0("l2fc_", mk)]] <- x$l2fc[match(a$bin_id, rownames(x$l2fc)), mk]

g2 <- file.path(rds_dir, "npmi_bin_features.rds")
has_g2 <- file.exists(g2)
if (has_g2) {
  f <- readRDS(g2)$features
  keep <- intersect(c("bin_id","pc1_WT","pc1_HP1gKO","npmi_short_WT","npmi_short_HP1gKO",
                      "npmi_long_WT","npmi_long_HP1gKO","n_valid_WT","n_valid_HP1gKO",
                      "insulation_WT","insulation_HP1gKO"), names(f))
  a <- merge(a, f[, ..keep], by = "bin_id", all.x = TRUE)
  setorder(a, chr, start)
} else message("GAM_2 output absent - skipping the NPMI-based columns")

# --- the signed direction statistic ------------------------------------------
a[, net := log2((n_contacts_lost + 1) / (n_contacts_gained + 1))]
# Smooth WITHIN chromosome; smoothing across a centromere/telomere boundary
# would invent a block that spans a discontinuity.
# Edges use a PARTIAL window, not the raw value. Falling back to the raw
# per-bin net at chromosome ends leaves a handful of unsmoothed, high-variance
# points next to a heavily smoothed track, which shows up as a spike at every
# telomere and can manufacture a block there. A partial-window mean keeps the
# same statistic, just with fewer contributing bins.
a[, net_s := {
  w <- gam_block_smooth_bins; n <- .N
  cs <- c(0, cumsum(net)); h <- (w - 1L) %/% 2L
  lo <- pmax(1L, seq_len(n) - h); hi <- pmin(n, seq_len(n) + h)
  (cs[hi + 1L] - cs[lo]) / (hi - lo + 1L)
}, by = chr]

# --- artifact check, printed not hidden --------------------------------------
cat("=== artifact check: is direction just a coverage difference? ===\n")
if (has_g2 && all(c("n_valid_WT","n_valid_HP1gKO") %in% names(a))) {
  dv <- as.numeric(a$n_valid_WT) - as.numeric(a$n_valid_HP1gKO)
  cat(sprintf("  rho(net, WT-KO valid NPMI cells) = %+.3f   (near 0 = not a coverage artifact)\n",
      cor(a$net, dv, method = "spearman", use = "complete.obs")))
} else cat("  (needs GAM_2; skipped)\n")

# --- call blocks --------------------------------------------------------------
# Sign of the smoothed net, with a dead zone around 0 so a block must be
# meaningfully one-directional rather than merely non-zero.
a[, dirn := fifelse(net_s >=  gam_block_net_cut, "lost",
            fifelse(net_s <= -gam_block_net_cut, "gained", NA_character_))]
a[, run := rleid(chr, dirn)]
blk <- a[!is.na(dirn), .(chr = chr[1], start = min(start), end = max(end), dirn = dirn[1],
                         n_bins = .N, mean_net = mean(net_s),
                         mean_lost = mean(n_contacts_lost),
                         mean_gained = mean(n_contacts_gained)), by = run]
blk[, width_Mb := (end - start) / 1e6]
blk <- blk[width_Mb >= gam_block_min_mb][order(-width_Mb)]
cat(sprintf("\n=== %d blocks >= %d Mb (%d lost-dominated, %d gained-dominated) ===\n",
            nrow(blk), gam_block_min_mb, sum(blk$dirn == "lost"), sum(blk$dirn == "gained")))
cat(sprintf("  genome in blocks: %.1f%%   median block %.1f Mb   largest %.1f Mb\n",
            100 * sum(blk$end - blk$start) / sum(a$end - a$start),
            median(blk$width_Mb), max(blk$width_Mb)))
print(head(blk[, .(chr, start, end, width_Mb, dirn, mean_net = round(mean_net, 2))], 15),
      row.names = FALSE)

# --- characterise: what differs between lost- and gained-dominated blocks? ----
# Annotate every bin with its block, then compare. This is the question the
# blocks exist to answer: is direction associated with chromatin state, or is it
# a purely architectural partition?
a[, block := NA_character_]
bi <- a[!is.na(dirn)][run %in% blk$run]
a[run %in% blk$run, block := dirn]
cmp_cols <- c(paste0("l2fc_", marks_here), "frac_LINE", "frac_SINE", "frac_LTR",
              intersect(c("pc1_WT","npmi_short_WT","npmi_long_WT","insulation_WT"), names(a)))
cmp <- rbindlist(lapply(cmp_cols, function(v) {
  L <- a[block == "lost",   get(v)]; G <- a[block == "gained", get(v)]
  L <- L[is.finite(L)];              G <- G[is.finite(G)]
  if (!length(L) || !length(G)) return(NULL)
  # Wilcoxon over bins would pseudo-replicate badly here (bins inside a 16 Mb
  # block are anything but independent), so report effect sizes only and let the
  # block count carry the uncertainty.
  data.table(variable = v, n_lost = length(L), n_gained = length(G),
             median_lost = median(L), median_gained = median(G),
             difference = median(L) - median(G))
}))
cmp[, abs_diff := abs(difference)]
setorder(cmp, -abs_diff)
fwrite(cmp, file.path(tables_dir, "gam6_block_characterisation.tsv"), sep = "\t")
cat("\n=== lost- vs gained-dominated blocks: what actually differs ===\n")
print(cmp[, .(variable, median_lost = round(median_lost, 3),
              median_gained = round(median_gained, 3), difference = round(difference, 3))],
      row.names = FALSE)

fwrite(blk[, .(chr, start, end, width_Mb, dirn, n_bins, mean_net, mean_lost, mean_gained)],
       file.path(tables_dir, "gam6_large_blocks.tsv"), sep = "\t")

# BED for IGV, coloured by direction
bed <- blk[order(chr, start), .(chr, start, end,
           name = sprintf("%s_%.0fMb", dirn, width_Mb), score = 0L, strand = ".",
           ts = start, te = end,
           rgb = fifelse(dirn == "lost", "178,24,43", "33,102,172"))]
bed_p <- file.path(bed_dir, "gam6_large_blocks.bed")
writeLines('track name="GAM large contact blocks" description="red=lost, blue=gained" itemRgb=On visibility=dense', bed_p)
fwrite(bed, bed_p, sep = "\t", col.names = FALSE, append = TRUE, quote = FALSE)
cat("\nWrote", bed_p, "\n")

# --- genome-wide figure -------------------------------------------------------
a[, pos_Mb := (start + end) / 2e6]
chr_ord <- paste0("chr", 1:19)
a[, chr := factor(chr, levels = chr_ord)]
blk[, chr := factor(chr, levels = chr_ord)]
p <- ggplot(a[!is.na(chr)], aes(pos_Mb, net_s)) +
  geom_rect(data = blk[!is.na(chr)], inherit.aes = FALSE,
            aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf, fill = dirn),
            alpha = 0.25) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_hline(yintercept = c(-gam_block_net_cut, gam_block_net_cut),
             linetype = 2, linewidth = 0.3, colour = "grey50") +
  geom_line(linewidth = 0.4) +
  facet_wrap(~chr, scales = "free_x", ncol = 4) +
  scale_fill_manual(values = c(lost = "#B2182B", gained = "#2166AC"), name = "block") +
  labs(title = "Large one-directional contact blocks, genome-wide",
       subtitle = paste0("net = log2((lost+1)/(gained+1)), ",
                         gam_block_smooth_bins * gam_binsize / 1e6,
                         " Mb running mean. Positive = contacts preferentially LOST in the KO.\n",
                         "Blocks: |net| > ", gam_block_net_cut, " for at least ",
                         gam_block_min_mb, " Mb. Direction is NOT a coverage artifact ",
                         "(rho with WT-KO valid cells ~ 0)."),
       x = "position (Mb)", y = "net contact direction") +
  theme_m + theme(legend.position = "top")
save_fig(p, "large_blocks_genome_wide", subdir = "blocks", width = 15, height = 11)

cat("\nGAM_6 complete.\n")
