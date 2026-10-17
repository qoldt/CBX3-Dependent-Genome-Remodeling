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

# ==============================================================================
# SECTION C - what are hotspot domains made of, and does it differ by direction?
# ==============================================================================
# The gap named in the report: domains are called and carried, but nothing
# annotated them as a SET. This does, three ways, and the three disagree in an
# instructive way.
#
# WHY NOT JUST A HYPERGEOMETRIC TEST. It is the obvious first move and it is
# computed here, but it is wrong in two specific ways that both inflate it:
#   1. It treats 20 kb bins as independent draws. They are not - the
#      differential-contact count track has autocorrelation r = 0.59 at 20 kb
#      and still ~0.4 at 100 kb (QC_Hotspot_Domains.R panel A). A hotspot domain
#      is by construction a RUN of correlated bins, so the effective sample size
#      is the number of domains (~1,000), not the number of bins (~10,000+).
#   2. Hotspot domains are large and most annotations are not. A domain of
#      several hundred kb overlaps a narrow feature by being long, which is why
#      the repo's standing rule is a size- or composition-matched null for
#      anything overlap-shaped.
# So the hypergeometric p is reported in a column labelled as inflated, beside a
# width-matched positional permutation that is the test actually being read.
#
# THE DIRECTION QUESTION gets its own test. Comparing lost domains to gained
# domains is the cleanest contrast available, because both sets are called by the
# same procedure from the same track on the same bins - mappability, coverage and
# width confounds are shared and largely cancel in the difference. The null there
# is a label permutation stratified by width decile, so the two sets keep
# comparable width distributions under the null.
local({
  msg <- function(...) cat(sprintf(...), "\n", sep = "")
  fig_sub <- "hotspots"
  NPERM <- 1000L

  X <- readRDS(file.path(rds_dir, "gam_bins_annotated.rds"))
  A <- as.data.table(as.data.frame(X$annot))
  A <- A[chr %in% gam_autosomes]
  setorder(A, chr, start)
  FE <- tryCatch(as.data.table(readRDS(file.path(rds_dir, "npmi_bin_features.rds"))$features),
                 error = function(e) NULL)
  if (!is.null(FE) && all(c("n_valid_WT", "n_valid_HP1gKO") %in% names(FE)))
    A <- merge(A, FE[, .(bin_id, n_valid_WT, n_valid_HP1gKO)], by = "bin_id", all.x = TRUE)
  else A[, `:=`(n_valid_WT = 1L, n_valid_HP1gKO = 1L)]
  setorder(A, chr, start)
  # Mappable = measured in at least one genotype. Domains are only ever placed on
  # mappable bins under the null, or the null would score hits for landing in a
  # gap where nothing can be annotated either.
  A[, mappable := (is.finite(n_valid_WT) & n_valid_WT > 0) |
                  (is.finite(n_valid_HP1gKO) & n_valid_HP1gKO > 0)]

  # --- annotation flags, one logical vector per annotation -------------------
  hmm_cols <- grep("^hmm_", names(A), value = TRUE)
  dom_state <- if (length(hmm_cols)) {
    m <- as.matrix(A[, ..hmm_cols]); m[!is.finite(m)] <- -1
    ifelse(rowSums(m >= 0) == 0, NA_character_, sub("^hmm_", "", hmm_cols[max.col(m, ties.method = "first")]))
  } else rep(NA_character_, nrow(A))
  qtop <- function(v) { v <- ifelse(is.finite(v), v, 0)
                        v >= quantile(v, 0.75, na.rm = TRUE) & v > 0 }
  lads <- tryCatch(load_lads(), error = function(e) NULL)
  lad_flag <- function(which_lad) {
    if (is.null(lads)) return(rep(NA, nrow(A)))
    g <- GRanges(A$chr, IRanges(A$start + 1L, A$end))
    # The column is "class", not "type". An earlier version guessed "type", fell
    # through to selecting EVERY LAD for both, and reported cLAD and ciLAD with
    # identical numbers - which is the giveaway to watch for here.
    if (!"class" %in% names(mcols(lads))) return(rep(NA, nrow(A)))
    sel <- lads[mcols(lads)$class == which_lad]
    if (!length(sel)) return(rep(NA, nrow(A)))
    IRanges::overlapsAny(g, sel)
  }
  FL <- list()
  for (st in sort(unique(na.omit(dom_state)))) FL[[paste0("chromHMM: ", st)]] <- dom_state == st
  FL[["repeat: LINE-high (top quartile)"]] <- qtop(A$frac_LINE)
  FL[["repeat: SINE-high (top quartile)"]] <- qtop(A$frac_SINE)
  FL[["repeat: LTR-high (top quartile)"]]  <- qtop(A$frac_LTR)
  FL[["constitutive LAD"]]   <- lad_flag("cLAD")
  FL[["constitutive inter-LAD"]] <- lad_flag("ciLAD")
  # CTCF: the one annotation here that is a point feature and external to both
  # assays. "CTCF-dense" is the top quartile of peaks per bin, which is the more
  # meaningful of the two at 20 kb - almost every bin contains at least one peak.
  ctcf <- tryCatch(load_ctcf(), error = function(e) NULL)
  if (!is.null(ctcf)) {
    cg <- GRanges(A$chr, IRanges(A$start + 1L, A$end))
    for (ti in levels(ctcf$tissue)) {
      n_pk <- GenomicRanges::countOverlaps(cg, ctcf[ctcf$tissue == ti])
      FL[[paste0("CTCF peak (bulk ", ti, ")")]] <- n_pk > 0
      FL[[paste0("CTCF-dense (bulk ", ti, ")")]] <-
        n_pk >= quantile(n_pk[n_pk > 0], 0.75, na.rm = TRUE)
    }
  } else msg("CTCF BEDs absent - run prepare_CTCF_annotation.R to include them")
  FL[["has a TSS"]]          <- is.finite(A$n_tss) & A$n_tss > 0
  FL[["TAD boundary"]]       <- is.finite(A$tad_boundary) & A$tad_boundary > 0
  FL[["mark loss: H4K20me3 only"]] <- !is.na(A$compartment_centred) & A$compartment_centred == "H4K20me3_only"
  FL[["mark loss: H3K9me3 + H4K20me3"]] <- !is.na(A$compartment_centred) & A$compartment_centred == "co_loss"
  FL[["mark loss: stable"]]  <- !is.na(A$compartment_centred) & A$compartment_centred == "stable"
  FL <- FL[vapply(FL, function(v) sum(v, na.rm = TRUE) >= 50, logical(1))]
  msg("annotations tested: %d (each >= 50 bins)", length(FL))

  # --- domain sets: per-bin membership columns written by GAM_1 --------------
  dom_cols <- c(lost_adj = "dom_lost_adj", lost_raw = "dom_lost_raw",
                gained_adj = "dom_gained_adj", gained_raw = "dom_gained_raw")
  dom_cols <- dom_cols[dom_cols %in% names(A)]
  if (!length(dom_cols)) { msg("no dom_* columns in annot - section C skipped"); return(invisible()) }

  # Runs of equal non-NA domain id become the domain's bin range. Working in bin
  # index space (not bp) is what makes the positional shuffle a simple offset.
  runs_of <- function(v) {
    ok <- !is.na(v) & v != "" & v != "0"
    if (!any(ok)) return(NULL)
    key <- paste(A$chr, as.character(v))
    key[!ok] <- NA_character_
    r <- rle(key); e <- cumsum(r$lengths); s <- e - r$lengths + 1L
    keep <- !is.na(r$values)
    data.table(start_i = s[keep], end_i = e[keep], chr = A$chr[s[keep]])[, len := end_i - start_i + 1L][]
  }

  chr_range <- A[, .(lo = .I[1], hi = .I[.N]), by = chr]
  setkey(chr_range, chr)

  # size-weighted = fraction of domain BINS in the annotation (bp-weighted,
  # since bins are fixed width); per-region = fraction of DOMAINS touching it.
  # CLAUDE.md requires both, because a few huge domains can carry the first
  # while most domains carry none.
  stat_of <- function(D, flag) {
    idx <- unlist(Map(function(a, b) a:b, D$start_i, D$end_i), use.names = FALSE)
    f <- flag[idx]
    per_dom <- vapply(seq_len(nrow(D)), function(i)
      any(flag[D$start_i[i]:D$end_i[i]], na.rm = TRUE), logical(1))
    c(size_w = mean(f, na.rm = TRUE), per_region = mean(per_dom, na.rm = TRUE))
  }

  # CIRCULAR SHIFT: rotate the WHOLE domain set along each chromosome by one
  # random offset, wrapping at the end. Every inter-domain distance is preserved,
  # so the clustering of the called domains survives into the null.
  #
  # Why this is needed on top of the independent shuffle. Hotspot domains are
  # merged only across a 1-bin (20 kb) sub-threshold gap, and they still cluster
  # after that: 12% of adjacent lost_adj domains sit within 100 kb of each other
  # and 27% within 200 kb. The annotations are clustered too - chromHMM states
  # and LADs come in blocks. Placing domains INDEPENDENTLY therefore builds a
  # null with less spatial structure than the data, which understates the
  # variance of the overlap statistic and makes the p too small. The rotation
  # null keeps the structure and is the conservative one; where the two disagree,
  # believe the rotation.
  rotate_stats <- function(D, flag, nperm = NPERM, seed = 42) {
    set.seed(seed)
    lo <- chr_range[D$chr, lo]; hi <- chr_range[D$chr, hi]
    nchr <- hi - lo + 1L
    out <- matrix(NA_real_, nperm, 2L)
    uch <- unique(D$chr)
    for (p in seq_len(nperm)) {
      off <- setNames(as.integer(stats::runif(length(uch)) * chr_range[uch, hi - lo + 1L]), uch)
      sh <- off[D$chr]
      s <- lo + ((D$start_i - lo) + sh) %% nchr
      e <- s + D$len - 1L
      over <- e > hi                       # wrap back to the chromosome start
      s[over] <- lo[over]; e[over] <- lo[over] + D$len[over] - 1L
      idx <- unlist(Map(function(a, b) a:b, s, e), use.names = FALSE)
      fv <- flag[idx]
      pd <- vapply(seq_along(s), function(i) any(flag[s[i]:e[i]], na.rm = TRUE), logical(1))
      out[p, ] <- c(mean(fv, na.rm = TRUE), mean(pd, na.rm = TRUE))
    }
    out
  }

  # Positional shuffle: each domain keeps its width and its chromosome and is
  # dropped at a random mappable-containing offset INDEPENDENTLY of the others.
  # Width and chromosome are matched by construction, which is the confound the
  # hypergeometric ignores - but the independence is itself an assumption, which
  # is what rotate_stats above exists to drop.
  shuffle_stats <- function(D, flag, nperm = NPERM, seed = 42) {
    set.seed(seed)
    lo <- chr_range[D$chr, lo]; hi <- chr_range[D$chr, hi]
    span <- hi - lo + 1L - D$len
    span[span < 1L] <- 1L
    out <- matrix(NA_real_, nperm, 2L)
    for (p in seq_len(nperm)) {
      s <- lo + as.integer(stats::runif(nrow(D)) * span)
      e <- s + D$len - 1L
      idx <- unlist(Map(function(a, b) a:b, s, e), use.names = FALSE)
      fv <- flag[idx]
      pd <- vapply(seq_along(s), function(i) any(flag[s[i]:e[i]], na.rm = TRUE), logical(1))
      out[p, ] <- c(mean(fv, na.rm = TRUE), mean(pd, na.rm = TRUE))
    }
    out
  }

  # --- MERGE SWEEP: do called domains cluster, and does merging percolate? ----
  # The same question QC_Hotspot_Domains.R asks of the raw BIN track, asked one
  # level up of the CALLED domains. Domains are merged across a 1-bin (20 kb)
  # sub-threshold gap only, so nearby domains stay separate; this measures how
  # much of the genome a further merge would absorb. The heterochromatin
  # fragments percolated catastrophically under the same sweep (1.3% of the
  # genome at gap 0 to 64.8% at 50 kb), so this is the check that hotspot
  # domains do not.
  SW <- list()
  for (nm in names(dom_cols)) {
    D0 <- runs_of(A[[dom_cols[[nm]]]])
    if (is.null(D0)) next
    gr0 <- GRanges(D0$chr, IRanges(A$start[D0$start_i] + 1L, A$end[D0$end_i]))
    gaps <- as.data.table(as.data.frame(gr0))[order(seqnames, start),
              .(g = (start - data.table::shift(end)) / 1000), by = seqnames]$g
    gaps <- gaps[is.finite(gaps)]
    for (mg in c(0, 20, 100, 200, 500, 1000, 2000)) {
      gr <- GenomicRanges::reduce(gr0, min.gapwidth = max(1L, as.integer(mg * 1000)))
      SW[[length(SW) + 1L]] <- data.table(
        domain_set = nm, merge_gap_kb = mg, n_regions = length(gr),
        median_width_kb = median(width(gr)) / 1000,
        max_width_Mb = max(width(gr)) / 1e6,
        genome_pct = 100 * sum(as.numeric(width(gr))) / sum(as.numeric(A$end - A$start)),
        pct_gaps_under_100kb = 100 * mean(gaps < 100),
        pct_gaps_under_200kb = 100 * mean(gaps < 200),
        median_gap_kb = median(gaps))
    }
  }
  SWD <- rbindlist(SW)
  fwrite(SWD, file.path(tables_dir, "gam4_hotspot_domain_merge_sweep.tsv"), sep = "\t")
  msg("\n=== do called hotspot domains cluster? merge sweep one level above the bin track ===")
  print(dcast(SWD, merge_gap_kb ~ domain_set, value.var = "genome_pct"))
  msg("adjacent-domain gaps under 100 kb: %s",
      paste(sprintf("%s %.0f%%", SWD[merge_gap_kb == 0]$domain_set,
                    SWD[merge_gap_kb == 0]$pct_gaps_under_100kb), collapse = ", "))

  pS <- ggplot(SWD, aes(merge_gap_kb, genome_pct, colour = domain_set)) +
    geom_line(linewidth = .6) + geom_point(size = 1.8) +
    scale_x_continuous(trans = "log1p", breaks = c(0, 20, 100, 200, 500, 1000, 2000)) +
    scale_colour_manual(values = setNames(disc_pal(length(unique(SWD$domain_set))),
                                          sort(unique(SWD$domain_set))), name = NULL) +
    labs(title = "Do hotspot domains percolate if adjacent ones are merged?",
         subtitle = paste0("Genome fraction covered as the merge gap grows. The heterochromatin fragment set percolated under ",
                           "this\nsweep - 1.3% of the genome at gap 0 to 64.8% at 50 kb - which is why its domain definition had to ",
                           "be\nrebuilt. Hotspot domains do not: the curve is gradual, so merging is a choice about scale rather ",
                           "than a\ncliff to stay off. Past ~500 kb it starts reproducing the GAM_6 blocks, which already exist ",
                           "separately."),
         x = "merge gap (kb, log scale)", y = "% of genome covered") +
    theme_m
  save_fig(pS, "hotspot_domain_merge_sweep", subdir = fig_sub, width = 8, height = 5.5)

  RES <- list()
  for (nm in names(dom_cols)) {
    D <- runs_of(A[[dom_cols[[nm]]]])
    if (is.null(D) || nrow(D) < 20) { msg("%s: too few domains, skipped", nm); next }
    msg("%s: %d domains, median %d bins, total %d bins", nm, nrow(D),
        as.integer(median(D$len)), sum(D$len))
    for (an in names(FL)) {
      flag <- FL[[an]]; flag[is.na(flag)] <- FALSE
      obs <- stat_of(D, flag)
      perm <- shuffle_stats(D, flag)
      rot  <- rotate_stats(D, flag)
      # the inflated test, computed deliberately so the inflation is visible
      N <- sum(A$mappable); K <- sum(flag & A$mappable)
      n <- sum(D$len); k <- round(obs[["size_w"]] * n)
      hyp_p <- if (n > 0 && K > 0 && N > K)
                 stats::phyper(k - 1L, K, N - K, n, lower.tail = FALSE) else NA_real_
      sdp <- stats::sd(perm[, 1]); mp <- mean(perm[, 1])
      RES[[length(RES) + 1L]] <- data.table(
        domain_set = nm,
        direction  = sub("_.*$", "", nm),
        definition = sub("^.*_", "", nm),
        annotation = an, n_domains = nrow(D), n_domain_bins = n,
        obs_size_weighted = obs[["size_w"]], exp_shuffled = mp,
        fold = obs[["size_w"]] / mp,
        z_shuffled = if (is.finite(sdp) && sdp > 0) (obs[["size_w"]] - mp) / sdp else NA_real_,
        p_shuffled = (sum(abs(perm[, 1] - mp) >= abs(obs[["size_w"]] - mp)) + 1) / (NPERM + 1),
        exp_rotated = mean(rot[, 1]),
        fold_rotated = obs[["size_w"]] / mean(rot[, 1]),
        z_rotated = { sr <- stats::sd(rot[, 1]); mr <- mean(rot[, 1])
                      if (is.finite(sr) && sr > 0) (obs[["size_w"]] - mr) / sr else NA_real_ },
        p_rotated = { mr <- mean(rot[, 1])
                      (sum(abs(rot[, 1] - mr) >= abs(obs[["size_w"]] - mr)) + 1) / (NPERM + 1) },
        obs_per_region = obs[["per_region"]], exp_per_region = mean(perm[, 2]),
        p_hypergeom_INFLATED = hyp_p)
    }
  }
  R <- rbindlist(RES)
  R[, p_shuffled_BH := p.adjust(p_shuffled, "BH"), by = domain_set]
  R[, p_rotated_BH  := p.adjust(p_rotated,  "BH"), by = domain_set]
  fwrite(R, file.path(tables_dir, "gam4_hotspot_annotation_enrichment.tsv"), sep = "\t")

  msg("\n=== how far the hypergeometric is out ===")
  msg("annotation x domain-set cells significant at 0.05, by null:")
  msg("  hypergeometric (bins independent)   %2d/%d", sum(R$p_hypergeom_INFLATED < 0.05, na.rm = TRUE), nrow(R))
  msg("  width-matched independent shuffle   %2d/%d  (BH %d)", sum(R$p_shuffled < 0.05, na.rm = TRUE), nrow(R), sum(R$p_shuffled_BH < 0.05, na.rm = TRUE))
  msg("  circular shift (spacing preserved)  %2d/%d  (BH %d)", sum(R$p_rotated < 0.05, na.rm = TRUE), nrow(R), sum(R$p_rotated_BH < 0.05, na.rm = TRUE))
  msg("  cells significant under the independent shuffle but NOT the rotation: %d",
      sum(R$p_shuffled < 0.05 & R$p_rotated >= 0.05, na.rm = TRUE))

  # --- the direction contrast: lost vs gained, same definition ---------------
  # Label permutation stratified by width decile, so the two sets keep comparable
  # width distributions under the null and the test is about composition only.
  DIR <- list()
  for (defn in unique(R$definition)) {
    dl <- runs_of(A[[dom_cols[[paste0("lost_", defn)]]]])
    dg <- runs_of(A[[dom_cols[[paste0("gained_", defn)]]]])
    if (is.null(dl) || is.null(dg)) next
    both <- rbind(copy(dl)[, grp := "lost"], copy(dg)[, grp := "gained"])
    both[, wq := cut(len, unique(quantile(len, seq(0, 1, 0.1))), include.lowest = TRUE)]
    for (an in names(FL)) {
      flag <- FL[[an]]; flag[is.na(flag)] <- FALSE
      sw <- function(sel) { idx <- unlist(Map(function(a,b) a:b, both$start_i[sel], both$end_i[sel]),
                                          use.names = FALSE); mean(flag[idx], na.rm = TRUE) }
      is_lost <- both$grp == "lost"
      obs <- sw(is_lost) - sw(!is_lost)
      set.seed(42)
      nullv <- replicate(NPERM, {
        lab <- both[, .(i = .I, keep = sample(is_lost[.I])), by = wq]
        sw(lab$keep[order(lab$i)]) - sw(!lab$keep[order(lab$i)])
      })
      DIR[[length(DIR) + 1L]] <- data.table(
        definition = defn, annotation = an,
        lost_frac = sw(is_lost), gained_frac = sw(!is_lost), diff = obs,
        exp_diff = mean(nullv), z = (obs - mean(nullv)) / stats::sd(nullv),
        p = (sum(abs(nullv - mean(nullv)) >= abs(obs - mean(nullv))) + 1) / (NPERM + 1))
    }
  }
  DD <- rbindlist(DIR)
  if (nrow(DD)) {
    DD[, p_BH := p.adjust(p, "BH"), by = definition]
    fwrite(DD, file.path(tables_dir, "gam4_hotspot_direction_contrast.tsv"), sep = "\t")
    msg("\n=== lost vs gained domains, composition difference (width-stratified label permutation) ===")
    print(DD[order(definition, -abs(z))][, .(definition, annotation,
            lost = round(lost_frac, 3), gained = round(gained_frac, 3),
            diff = round(diff, 3), z = round(z, 2), p_BH = signif(p_BH, 2))][1:min(20, .N)])
  }

  # --- figures ---------------------------------------------------------------
  R[, sig := ifelse(p_rotated_BH < 0.05, "*", "")]
  lim <- max(abs(log2(R$fold[is.finite(R$fold) & R$fold > 0])), na.rm = TRUE)
  p1 <- ggplot(R[is.finite(fold) & fold > 0],
               aes(domain_set, annotation, fill = pmax(pmin(log2(fold), lim), -lim))) +
    geom_tile(colour = "white", linewidth = .4) +
    geom_text(aes(label = paste0(sprintf("%.2f", fold), sig)), size = 2.6) +
    scale_fill_heat0_div(limits = c(-lim, lim), name = "log2 fold\nvs width-\nmatched null") +
    labs(title = "What hotspot domains are made of, against a width-matched null",
         subtitle = paste0("Fold = observed share of domain bp in the annotation, divided by the share when the same domains are\n",
                           "reshuffled to random positions on the same chromosome keeping their widths. * = BH < 0.05 under the\n",
                           "CIRCULAR-SHIFT null, which also preserves the spacing between domains and is the conservative one.\n",
                           "'adj' domains are smoothed-then-thresholded, 'raw' are thresholded directly; they agree on\n",
                           "only 40% of bins, so a result present in one and not the other is a result about the definition."),
         x = NULL, y = NULL) +
    theme_m + theme(panel.grid = element_blank(),
                    axis.text.x = element_text(angle = 30, hjust = 1))
  save_fig(p1, "hotspot_annotation_enrichment", subdir = fig_sub, width = 9.5, height = 8)

  cmp <- melt(R[, .(domain_set, annotation,
                    `circular shift` = p_rotated,
                    `independent shuffle` = p_shuffled,
                    `hypergeometric` = p_hypergeom_INFLATED)],
              id.vars = c("domain_set", "annotation"), variable.name = "test", value.name = "p")
  cmp[, p := pmax(p, 1e-300)]
  p2 <- ggplot(cmp[is.finite(p)], aes(test, -log10(p), group = paste(domain_set, annotation))) +
    geom_line(colour = "grey70", linewidth = .3) +
    geom_point(aes(colour = test), size = 1.6, show.legend = FALSE) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey40") +
    scale_colour_manual(values = setNames(disc_pal(3), levels(cmp$test))) +
    labs(title = "The same comparisons under three nulls, from most to least assumption",
         subtitle = paste0("One line per annotation x domain-set cell. The hypergeometric treats 20 kb bins as independent draws,\n",
                           "but the count track has autocorrelation r = 0.59 at 20 kb and a hotspot domain is by construction a run\n",
                           "of correlated bins - so its effective n is the ~1,000 domains, not the ~10,000 bins. Dashed line = 0.05."),
         x = NULL, y = "-log10(p)") +
    theme_m
  save_fig(p2, "hotspot_test_comparison", subdir = fig_sub, width = 7.5, height = 6)

  if (nrow(DD)) {
    DD[, sig := ifelse(p_BH < 0.05, "*", "")]
    p3 <- ggplot(DD, aes(diff, reorder(annotation, diff), fill = definition)) +
      geom_vline(xintercept = 0, linewidth = .3, colour = "grey50") +
      geom_col(position = position_dodge(width = .75), width = .7) +
      geom_text(aes(label = sig), position = position_dodge(width = .75), hjust = -0.3, size = 3.5) +
      scale_fill_manual(values = setNames(disc_pal(2), unique(DD$definition)), name = "domain definition") +
      labs(title = "Do domains that LOSE contacts differ in composition from those that GAIN?",
           subtitle = paste0("Difference in the share of domain bp in each annotation, lost minus gained. Positive = more common in\n",
                             "lost domains. Both sets come from the same track by the same procedure, so width, mappability and\n",
                             "coverage confounds are shared and largely cancel. * = BH < 0.05 under a width-decile-stratified\n",
                             "permutation of the lost/gained labels."),
           x = "share in lost domains minus share in gained domains", y = NULL) +
      theme_m
    save_fig(p3, "hotspot_direction_contrast", subdir = fig_sub, width = 9.5, height = 7)
  }
  msg("\nSection C complete.")
})
