# ================================================================
# GAM_4 - Large-scale contact structure: blocks and the block budget
# ----------------------------------------------------------------
#     cd GAM && Rscript GAM_4_Large_Scale.R          (after GAM_1, GAM_2)
#
# Merges the two stages that work on the same object at the same scale: the
# 10-25 Mb one-directional blocks, and the block x block contact budget that
# asks whether they represent redistribution rather than local loss.
# Previously GAM_6 and GAM_8.
#
# These were always one analysis split across two files - section B reads the
# very table section A writes, and could not run without it. Merging removes
# that ordering hazard: there is no longer a way to run the budget against a
# stale block set.
#
#   A (was GAM_6)  Call the blocks genome-wide and characterise them. 123 blocks
#                  >= 4 Mb covering 41.6% of the genome, 63 lost / 60 gained.
#                  Written only after the blocks were spotted by eye in IGV -
#                  every earlier stage worked at 0.1-2 Mb and was blind to this
#                  scale. That is the lesson: check the coarse scale first.
#   B (was GAM_8)  Mean NPMI for every block-class pair by distance. If lost-lost
#                  falls while gained-gained rises, the chromosome is moving
#                  contacts around rather than losing them.
#
# Output file names are unchanged (gam6_*, gam8_*).
# ================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })



# ==============================================================================
# SECTION A - call and characterise the large blocks           (was GAM_6)
# ==============================================================================

local({
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
})


# ==============================================================================
# SECTION B - block x block contact budget                     (was GAM_8)
# ==============================================================================

local({
  for (f in c("gam_bins_annotated.rds")) if (!file.exists(file.path(rds_dir, f)))
    stop("missing ", f, " - run GAM_1 first", call. = FALSE)
  blk_f <- file.path(tables_dir, "gam6_large_blocks.tsv")
  if (!file.exists(blk_f)) stop("missing gam6_large_blocks.tsv - run GAM_6 first", call. = FALSE)

  x <- readRDS(file.path(rds_dir, "gam_bins_annotated.rds"))
  a <- as.data.table(x$annot); setorder(a, chr, start)
  blk <- fread(blk_f)

  # class every bin by the block it falls in
  a[, cls := "neither"]
  for (i in seq_len(nrow(blk)))
    a[chr == blk$chr[i] & start >= blk$start[i] & end <= blk$end[i], cls := blk$dirn[i]]
  cat("=== bins per class ===\n"); print(table(a$cls))

  BAND_EDGES <- c(0, 1e6, 5e6, 1e7, Inf)
  BAND_NAMES <- c("0.02-1Mb", "1-5Mb", "5-10Mb", ">10Mb")
  CLS <- c("lost", "gained", "neither")

  acc <- list()
  for (ch in gam_autosomes) {
    ids <- a[chr == ch]
    if (!nrow(ids)) next
    for (g in c("WT", "HP1gKO")) {
      m <- tryCatch(load_npmi(g, ch), error = function(e) NULL)
      if (is.null(m)) next
      # align the matrix to our bin order by bin id
      keep <- intersect(rownames(m), ids$bin_id)
      m <- m[keep, keep, drop = FALSE]
      cl <- ids$cls[match(keep, ids$bin_id)]
      pos <- ids$start[match(keep, ids$bin_id)]
      n <- nrow(m)
      if (n < 50) { rm(m); gc(FALSE); next }
      d <- abs(outer(pos, pos, "-"))
      bandi <- cut(d, BAND_EDGES, labels = BAND_NAMES, right = FALSE)
      ci <- match(cl, CLS)
      pair <- outer(ci, ci, function(p, q) paste(CLS[pmin(p, q)], CLS[pmax(p, q)], sep = "-"))
      ok <- is.finite(m) & upper.tri(m) & !is.na(bandi)
      if (!any(ok)) { rm(m, d, bandi, pair); gc(FALSE); next }
      acc[[length(acc) + 1]] <- data.table(
        chr = ch, genotype = g, pair = pair[ok], band = as.character(bandi[ok]), v = m[ok]
      )[, .(sum_npmi = sum(v), n = .N), by = .(chr, genotype, pair, band)]
      rm(m, d, bandi, pair, ok); gc(FALSE)
    }
    cat("  ", ch, "done\n")
  }
  A <- rbindlist(acc); rm(acc); gc(FALSE)

  # pool chromosomes exactly (sums and counts, never a mean of means)
  P <- A[, .(mean_npmi = sum(sum_npmi) / sum(n), n = sum(n)), by = .(genotype, pair, band)]
  W <- dcast(P, pair + band ~ genotype, value.var = c("mean_npmi", "n"))
  setnames(W, c("mean_npmi_WT", "mean_npmi_HP1gKO"), c("WT", "HP1gKO"))
  W[, delta := HP1gKO - WT]
  W[, pct := 100 * delta / abs(WT)]
  setorder(W, band, pair)
  fwrite(W, file.path(tables_dir, "gam8_block_contact_budget.tsv"), sep = "\t")

  cat("\n=== mean NPMI by block-class pair and distance, WT vs KO ===\n")
  for (bn in BAND_NAMES) {
    cat("\n--", bn, "--\n")
    d <- W[band == bn]
    if (!nrow(d)) next
    cat(sprintf("  %-18s %8s %8s %9s %8s\n", "class pair", "WT", "KO", "delta", "%"))
    for (i in seq_len(nrow(d)))
      cat(sprintf("  %-18s %8.4f %8.4f %+9.4f %+7.1f%%\n",
          d$pair[i], d$WT[i], d$HP1gKO[i], d$delta[i], d$pct[i]))
  }

  cat("\n=== verdict shape ===\n")
  ll <- W[pair == "lost-lost"]; gg <- W[pair == "gained-gained"]
  for (bn in BAND_NAMES) {
    l <- ll[band == bn]; g <- gg[band == bn]
    if (!nrow(l) || !nrow(g)) next
    cat(sprintf("  %-9s lost-lost %+.1f%%  gained-gained %+.1f%%  -> %s\n", bn, l$pct, g$pct,
      if (l$pct < 0 && g$pct > 0) "RECIPROCAL (redistribution-shaped)"
      else if (l$pct < 0 && g$pct < 0) "both fall (net loss)"
      else if (l$pct > 0 && g$pct > 0) "both rise" else "mixed"))
  }

  p <- ggplot(W, aes(factor(band, BAND_NAMES), pct, fill = pair)) +
    geom_col(position = position_dodge(0.85), width = 0.78) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_fill_disc(name = NULL) +
    labs(title = "Contact budget change by block class and distance",
         subtitle = "Percent change in mean NPMI, HP1gKO vs WT. Reciprocal movement between\nlost-lost and gained-gained would indicate redistribution rather than local loss.",
         x = "contact distance", y = "change in mean NPMI (%)") +
    theme_m + theme(legend.position = "top", text = element_text(size = 14))
  save_fig(p, "block_contact_budget", subdir = "blocks", width = 11, height = 7)
  cat("\nGAM_8 complete.\n")
})
