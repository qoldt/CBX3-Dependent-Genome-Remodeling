# ================================================================
# GAM_1 - the shared 20 kb coordinate system: ChIP signal + annotation per bin
# ----------------------------------------------------------------
# WHY THIS STAGE EXISTS. GAM contacts live on a fixed 20 kb grid; the ChIP result
# lives on peaks. Every cross-assay question ("do bins that lose contacts also
# lose H4K20me3?") needs both on ONE unit. Projecting peak-level results onto
# bins would inherit the peak set's selection (only called regions, 0 elsewhere);
# quantifying the bigWigs DIRECTLY over the grid gives a value for every bin and
# is the same move MINUTE_7 makes for gene bodies. This stage writes the table
# that GAM_2..GAM_5 join on (key: bin_id, e.g. "chr19:0-20000").
#
# WHAT IT DOES: bw_loci() of all five marks over the 123,133 autosomal bins ->
# DESeq2 per mark (sizeFactors = 1, HP1gKO vs WT) -> per-bin annotation
# (ChromHMM, repeats, TAD boundary, TSS, hotspot counts) -> MINUTE_3-style
# compartment labels per bin. Output: results/rds/gam_bins_annotated.rds
#
# Run from GAM/:   cd GAM && Rscript GAM_1_Bins_and_Signal.R
# SUSTAINED ALL-CORE JOB (11 GB of bigWigs). Raw per-bin coverage is cached in
# results/counts/, so a re-run (or a run that only drops samples) skips it.
# ================================================================
source("config.R")

# --- ChIP inputs, resolved WITHOUT sourcing ../ChIP/config.R ------------------
# Sample sheet, ChIP data roots, COUNT_SCALE, core count and the MINUTE_3
# compartment cutoffs all come from config.R - the single source of truth. In
# particular gam_load_samples() applies the SAME default HP1gKO-rep2 exclusion
# as the ChIP side, which is what keeps these per-bin log2FC comparable to the
# published peak-level numbers. config.R also documents why we must NOT
# source("../ChIP/config.R") here (it would redirect every GAM output into
# ChIP/results).
samples      <- gam_load_samples()
excluded_ids <- attr(samples, "excluded_ids")
COUNT_SCALE  <- gam_count_scale
n_cores      <- gam_cores
stopifnot(setequal(unique(samples$mark), mark_levels))

# --- 1. The master grid -------------------------------------------------------
bins <- gam_bins()
bins$bin_id <- gam_bin_id(bins)
bin_ids <- bins$bin_id
stopifnot(length(bins) == 123133L, !anyDuplicated(bin_ids))
cat(sprintf("Master grid: %d autosomal %d bp bins\n", length(bins), gam_binsize))

# --- 2. Quantify every mark over the bins -------------------------------------
# bigWigs are NCBI-style ("1"), the grid is UCSC ("chr1"): quantify on a flipped
# copy. bw_loci() de-duplicates and may reorder, so results are mapped back by
# coordinate key (as MINUTE_7's bw_loci_aligned) instead of trusted positionally.
bins_ncbi <- bins; mcols(bins_ncbi) <- NULL
bins_ncbi <- to_ncbi(bins_ncbi)
stopifnot(length(bins_ncbi) == length(bins))
bin_key <- paste0(seqnames(bins_ncbi), ":", start(bins_ncbi), "-", end(bins_ncbi))

bw_bins <- function(bwfile) {
  o <- wigglescout::bw_loci(bwfile, bins_ncbi)
  v <- mcols(o)[[1]]
  names(v) <- paste0(seqnames(o), ":", start(o), "-", end(o))
  unname(v[bin_key])
}

count_tables <- list()
for (mark in mark_levels) {
  cat("\n=== Binned coverage:", mark, "===\n")
  s  <- samples[as.character(samples$mark) == mark, , drop = FALSE]
  bw <- setNames(s$bigwig_path, s$sample_id)

  # Reuse the cache only if it covers every wanted sample on the SAME grid;
  # columns are samples so a sensitivity run that drops samples needs no re-read.
  cached <- file.path(counts_dir, paste0(mark, "_bin_counts.tsv.gz"))
  if (file.exists(cached)) {
    cdt <- fread(cached)
    if (all(names(bw) %in% names(cdt)) && nrow(cdt) == length(bins) &&
        identical(as.character(cdt$bin_id), bin_ids)) {
      cat(" - Reusing cached counts:", cached, "\n")
      count_tables[[mark]] <- as.matrix(cdt[, names(bw), with = FALSE])
      rownames(count_tables[[mark]]) <- bin_ids
      next
    }
    cat(" - Cached counts unusable (samples or grid differ); re-quantifying\n")
  }

  miss <- !file.exists(bw)
  if (any(miss)) stop("missing bigWig(s) - download per the ChIP README or set MINUTE_DATA:\n  ",
                      paste(bw[miss], collapse = "\n  "), call. = FALSE)
  res <- mclapply(seq_along(bw), function(i) {
    message(sprintf("[%d/%d] %s", i, length(bw), names(bw)[i]))
    tryCatch(bw_bins(bw[[i]]), error = function(e) {
      message("Error in ", names(bw)[i], ": ", e$message); rep(NA_real_, length(bins)) })
  }, mc.cores = n_cores, mc.preschedule = FALSE)
  m <- do.call(cbind, lapply(res, function(x)
    if (is.numeric(x) && length(x) == length(bins)) x else rep(NA_real_, length(bins))))
  colnames(m) <- names(bw); rownames(m) <- bin_ids
  # A wholly-NA column is a failed read, not biology; stop rather than let
  # rowMeans(na.rm) quietly average over fewer samples.
  bad <- colSums(!is.na(m)) == 0
  if (any(bad)) stop("bw_loci failed for: ", paste(names(bw)[bad], collapse = ", "), call. = FALSE)
  fwrite(data.table(bin_id = bin_ids, as.data.frame(m)), cached, sep = "\t", quote = FALSE)
  cat("Saved:", cached, "\n")
  count_tables[[mark]] <- m
}

# --- 3. UN-ROUNDED per-bin signal (plotting / PC1 orientation / saddles) ------
# Same logic as mean_coverage_by_peak() in ChIP/config.R: rounding to integers
# collapses coverage onto a coarse lattice, so $signal is taken from the raw
# matrices. Columns follow the sample sheet; one column per sample, all marks.
signal <- do.call(cbind, unname(count_tables))
signal <- signal[, samples$sample_id, drop = FALSE]
rownames(signal) <- bin_ids

# --- 4. DESeq2 per mark --------------------------------------------------------
# sizeFactors = 1: the bigWigs are ALREADY input-scaled. DESeq2's median-of-ratios
# assumes most bins are unchanged and would divide the genome-wide H4K20me3 loss
# straight out - erasing the result. Contrast HP1gKO vs WT, WT as reference.
#
# COUNT_SCALE: mean coverage per 20 kb bin is ~0.2-3, and round() on that (as
# MINUTE_1 does per peak) leaves a handful of integers, a degenerate dispersion
# trend and a fake Poisson floor (the MINUTE_7 problem). Scaling by a constant
# before rounding leaves every log2 ratio unchanged; baseMean is divided back so
# it stays on the coverage scale. The only row filter is all-zero, which DESeq2
# cannot fit; MINUTE_1's rowSums > 10 would here discard most low-coverage bins.
stats_list <- list()
for (mark in mark_levels) {
  cat("\n=== DESeq2:", mark, "===\n")
  raw <- count_tables[[mark]]
  # NA = bin absent from the bigWig (no signal), so 0 for the model. Only
  # partial NAs reach here (wholly-NA samples already stopped).
  cat(sprintf(" - %.2f%% of cells NA -> 0\n", 100 * mean(is.na(raw))))
  raw[is.na(raw)] <- 0
  cnt  <- round(raw * COUNT_SCALE)
  keep <- which(rowSums(cnt) > 0)
  geno <- samples$genotype[match(colnames(cnt), samples$sample_id)]
  meta <- data.frame(row.names = colnames(cnt),
                     condition = factor(as.character(geno), levels = c("WT", "HP1gKO")))
  dds <- DESeqDataSetFromMatrix(countData = cnt[keep, , drop = FALSE], colData = meta,
                                design = ~ condition)
  sizeFactors(dds) <- rep(1, ncol(cnt))
  # parametric dispersion fit can fail on this many near-flat rows; fall back
  # to the local then mean trend rather than aborting a multi-hour stage.
  dds <- tryCatch(estimateDispersions(dds, quiet = TRUE), error = function(e) {
    message("parametric fit failed (", conditionMessage(e), "); trying local")
    tryCatch(estimateDispersions(dds, fitType = "local", quiet = TRUE), error = function(e2) {
      message("local fit failed; using mean"); estimateDispersions(dds, fitType = "mean", quiet = TRUE) })
  })
  dds <- nbinomWaldTest(dds, quiet = TRUE)
  r <- as.data.frame(results(dds, contrast = c("condition", "HP1gKO", "WT")))
  r <- r[match(seq_along(bin_ids), keep), , drop = FALSE]   # NA rows for dropped bins
  stats_list[[mark]] <- data.table(
    bin_id = bin_ids, mark = factor(mark, levels = mark_levels),
    baseMean = r$baseMean / COUNT_SCALE, log2FoldChange = r$log2FoldChange,
    pvalue = r$pvalue, padj = r$padj)
  cat(sprintf(" - %d of %d bins tested\n", length(keep), length(bin_ids)))
}
stats <- rbindlist(stats_list)

l2fc <- sapply(mark_levels, function(mk) stats_list[[mk]]$log2FoldChange)
dimnames(l2fc) <- list(bin_ids, mark_levels)

# --- 5. Per-bin annotation ----------------------------------------------------
# Annotation files are the same external ChIP annotation set (not committed).
need <- file.path(annot_dir, c("LINE.mm39.bed", "SINE.mm39.bed", "LTR.mm39.bed",
                               "TAD_boundaries_mm39.bed", "ChromHMM_18state_CortexP0_mm39.bed"))
if (any(!file.exists(need)))
  stop("annotation file(s) missing (download per the ChIP README or set MINUTE_DATA):\n  ",
       paste(need[!file.exists(need)], collapse = "\n  "), call. = FALSE)

# Fraction of each bin covered by a set of intervals. reduce() first so
# overlapping repeat copies are not double-counted (a fraction must be <= 1).
frac_covered <- function(bins, gr) {
  gr <- reduce(gr, ignore.strand = TRUE)
  f  <- numeric(length(bins))
  h  <- findOverlaps(bins, gr, ignore.strand = TRUE)
  if (length(h)) {
    w <- width(pintersect(bins[queryHits(h)], gr[subjectHits(h)], ignore.strand = TRUE))
    a <- tapply(w, queryHits(h), sum)
    f[as.integer(names(a))] <- as.numeric(a)
  }
  pmin(f / width(bins), 1)
}

# RepeatMasker UCSC tables: genoName/genoStart/genoEnd are columns 6-8. fread of
# those three is far cheaper than ChIP's read.table over all 17 columns, and the
# non-numeric header row (if any) drops out in the as.numeric filter, as there.
repeat_gr <- function(f) {
  d <- fread(f, header = FALSE, select = 6:8, colClasses = "character", sep = "\t")
  setnames(d, c("chr", "s", "e"))
  d[, `:=`(s = suppressWarnings(as.numeric(s)), e = suppressWarnings(as.numeric(e)))]
  d <- d[!is.na(s) & !is.na(e) & chr %in% seqlevels(bins)]
  GRanges(d$chr, IRanges(d$s + 1, d$e))   # BED 0-based -> 1-based
}
bg <- bins; mcols(bg) <- NULL
frac_rep <- sapply(c(LINE = "LINE", SINE = "SINE", LTR = "LTR"), function(cl)
  frac_covered(bg, repeat_gr(file.path(annot_dir, paste0(cl, ".mm39.bed")))))
colnames(frac_rep) <- paste0("frac_", colnames(frac_rep))

# ChromHMM: merge overlapping segments within a state first (the file has a few
# overlapping same-state segments that would otherwise double-count), as
# ChIP/config.R load_chromHMM(); then coverage fraction per state per bin.
hmm <- import(file.path(annot_dir, "ChromHMM_18state_CortexP0_mm39.bed"), format = "bed")
hmm$state <- as.character(hmm$name)
suppressWarnings(seqlevelsStyle(hmm) <- "UCSC")
hmm_states <- sort(unique(hmm$state))
frac_hmm <- sapply(hmm_states, function(s) frac_covered(bg, hmm[hmm$state == s]))
colnames(frac_hmm) <- paste0("hmm_", hmm_states)

tad <- import(file.path(annot_dir, "TAD_boundaries_mm39.bed"), format = "bed")
suppressWarnings(seqlevelsStyle(tad) <- "UCSC")
tad_boundary <- overlapsAny(bg, tad, ignore.strand = TRUE)

# TSS count per bin: UNIQUE transcript start positions, so isoforms sharing a
# TSS count once (otherwise gene-rich bins are inflated by isoform number).
txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
tss  <- unique(resize(transcripts(txdb), 1, fix = "start"))
mcols(tss) <- NULL; tss <- unique(granges(tss, use.names = FALSE))
n_tss <- countOverlaps(bg, tss, ignore.strand = TRUE)

# Overlap-rate assertions: a seqname-style mismatch yields ZERO overlaps, which
# reads downstream as a clean null, so each annotation must visibly hit.
stopifnot(mean(frac_rep[, "frac_LINE"] > 0) > 0.5, mean(frac_rep[, "frac_SINE"] > 0) > 0.5,
          mean(frac_rep[, "frac_LTR"] > 0) > 0.5,
          mean(rowSums(frac_hmm) > 0) > 0.8,
          any(tad_boundary), sum(n_tss) > 0)

# --- 6. Differential contact counts (hotspots) --------------------------------
# lost = selWT (stronger in WT, so LOST in HP1gKO); gained = selcKO. Columns are
# named from the KO's point of view; load_hotspots() already applies this.
hot <- load_hotspots()
hot[, bin_id := paste0(chr, ":", start, "-", end)]
in_hot <- bin_ids %in% hot$bin_id
hi <- match(bin_ids, hot$bin_id)
n_lost   <- ifelse(is.na(hi), 0L, as.integer(hot$n_contacts_lost[hi]))
n_gained <- ifelse(is.na(hi), 0L, as.integer(hot$n_contacts_gained[hi]))

# --- 7. Compartments: MAPPED from the MINUTE_3 peak labels --------------------
# These labels are assigned on the PEAKS and transferred to bins by overlap, NOT
# recomputed from per-bin log2FC. That distinction is the whole correctness of
# this block.
#
# WHY. MINUTE_3's cutoffs (-0.5 / -0.3 / 0.15) were tuned on the PEAK log2FC
# distribution, whose H4K20me3 IQR is 0.459. Averaging over 20 kb bins mixes
# peak signal with flanking background and compresses that IQR to 0.242 - 1.9x
# narrower. Re-applying an absolute -0.5 cutoff to the centred bin values
# therefore demands h4 < median - 0.5 = -0.997, which sits BELOW the 5th
# percentile of the bin distribution: it selected 147 bins (0.12%), giving 11
# H4K20me3_only and 114 co_loss where the peak analysis has thousands. A cutoff
# is not scale-free, and bins are a different scale.
#
# Mapping instead keeps the GAM compartments identical in definition to the
# published ChIP ones, which is the point - GAM_4 asks whether THOSE compartments
# remodel, so they must be the same objects.
#
# A bin gets the label of the labelled peak covering most of it (ties -> first).
# Bins overlapping no labelled peak stay NA, as most will: the labels are defined
# only on the shared H3K9me3/H4K20me3 peak set.
ar3 <- readRDS(require_chip_rds())
p20 <- ar3[["H4K20me3"]]; p9 <- ar3[["H3K9me3"]]
p20$k9_l2fc <- p9$log2FoldChange[match(p20$peak_id, p9$peak_id)]
p20 <- p20[is.finite(p20$k9_l2fc) & is.finite(p20$log2FoldChange), ]

assign_comp <- function(h4, k9) {
  out <- rep(NA_character_, length(h4))
  ok  <- !is.na(h4) & !is.na(k9)
  out[ok & h4 < H4_LOST_CUT & k9 < K9_LOSS_CUT]          <- "co_loss"
  out[ok & h4 < H4_LOST_CUT & abs(k9) < UNCH_CUT]        <- "H4K20me3_only"
  out[ok & abs(h4) < UNCH_CUT & abs(k9) < UNCH_CUT]      <- "stable"
  factor(out, levels = c("stable", "H4K20me3_only", "co_loss"))
}
# Peak-level medians, for the centred scheme - computed on peaks, where the
# cutoffs belong.
pmed_h4 <- median(p20$log2FoldChange, na.rm = TRUE)
pmed_k9 <- median(p20$k9_l2fc,        na.rm = TRUE)
pk_abs <- assign_comp(p20$log2FoldChange, p20$k9_l2fc)
pk_cen <- assign_comp(p20$log2FoldChange - pmed_h4, p20$k9_l2fc - pmed_k9)

pk_gr <- to_ucsc(GRanges(as.character(p20$chr),
                         IRanges(as.numeric(p20$start), as.numeric(p20$end))))
# Transfer by widest overlap. A 20 kb bin can span several labelled peaks, so
# "widest" is a deliberate choice over "first" - it makes the label reflect the
# territory the bin mostly is.
transfer <- function(lab) {
  h <- findOverlaps(bins, pk_gr, ignore.strand = TRUE)
  if (!length(h)) return(factor(rep(NA_character_, length(bins)), levels = levels(lab)))
  w  <- width(pintersect(bins[queryHits(h)], pk_gr[subjectHits(h)], ignore.strand = TRUE))
  lv <- as.character(lab)[subjectHits(h)]
  keep <- !is.na(lv)
  if (!any(keep)) return(factor(rep(NA_character_, length(bins)), levels = levels(lab)))
  d  <- data.table(q = queryHits(h)[keep], l = lv[keep], w = w[keep])
  agg <- d[, .(w = sum(w)), by = .(q, l)][order(-w)][, .SD[1], by = q]
  out <- rep(NA_character_, length(bins)); out[agg$q] <- agg$l
  factor(out, levels = levels(lab))
}
comp_abs <- transfer(pk_abs)
comp_cen <- transfer(pk_cen)
cat("Compartment labels mapped from", length(pk_gr), "labelled peaks\n")

# Bin-NATIVE sensitivity label, quantile-based so it IS scale-appropriate: the
# bottom quintile of each mark's own bin distribution counts as "lost", the
# middle three fifths as "unchanged". Reported alongside so a reader can see
# whether the peak-mapped grouping and a bin-intrinsic one agree; not primary,
# because it is not the published compartment definition.
h4 <- l2fc[, "H4K20me3"]; k9 <- l2fc[, "H3K9me3"]
med_h4 <- median(h4, na.rm = TRUE); med_k9 <- median(k9, na.rm = TRUE)
q4 <- quantile(h4, c(0.2, 0.4, 0.6), na.rm = TRUE)
q9 <- quantile(k9, c(0.2, 0.4, 0.6), na.rm = TRUE)
comp_binq <- {
  out <- rep(NA_character_, length(h4)); ok <- !is.na(h4) & !is.na(k9)
  out[ok & h4 < q4[1] & k9 < q9[1]]                                  <- "co_loss"
  out[ok & h4 < q4[1] & k9 >= q9[2] & k9 <= q9[3]]                   <- "H4K20me3_only"
  out[ok & h4 >= q4[2] & h4 <= q4[3] & k9 >= q9[2] & k9 <= q9[3]]    <- "stable"
  factor(out, levels = c("stable", "H4K20me3_only", "co_loss"))
}

# --- 6b. Hotspot DOMAINS, both adjusted and unadjusted -----------------------
# The per-bin counts above stay the primary continuous covariate. These add the
# second, independent-unit view. Two parameterisations are emitted so every
# domain-level claim downstream can be checked against the no-adjustment reading
# (see gam_hotspot_params / gam_hotspot_params_raw in config.R for the rationale).
# gam_hotspot_domains() groups by chromosome internally, so the property it needs
# is that starts increase WITHIN each chromosome - not a global sort. Asserting a
# global order would be wrong here anyway: gam_bins() returns natural seqlevel
# order (chr1, chr2, ... chr19) while order() on the character seqnames gives
# lexicographic order (chr1, chr10, chr11, ...), so the two legitimately differ.
stopifnot(!any(data.table(c = as.character(seqnames(bins)), s = start(bins))[
  , any(diff(s) <= 0), by = c]$V1))
dom_out <- list()
for (dir in c("lost", "gained")) {
  cnt <- if (dir == "lost") n_lost else n_gained
  for (variant in c("adj", "raw")) {
    pr <- if (variant == "adj") gam_hotspot_params else gam_hotspot_params_raw
    d  <- gam_hotspot_domains(cnt, as.character(seqnames(bins)),
                              start(bins) - 1L, end(bins),
                              q = pr$q, smooth_bins = pr$smooth_bins,
                              gap = pr$gap, min_bins = pr$min_bins)
    dom_out[[paste0(dir, "_", variant)]] <- d
    cat(sprintf("Hotspot domains %-6s %-3s: %5d domains, median %.0f bins, %.1f%% of genome\n",
        dir, variant, nrow(d$domains), median(c(d$domains$n_bins, 0)),
        100 * sum(d$domains$end - d$domains$start) / sum(end(bins) - (start(bins) - 1L))))
  }
}
for (k in names(dom_out))
  fwrite(dom_out[[k]]$domains, file.path(tables_dir, paste0("hotspot_domains_", k, ".tsv")), sep = "\t")

annot <- data.table(
  bin_id = bin_ids, chr = as.character(seqnames(bins)),
  start = start(bins) - 1L, end = end(bins),   # 0-based half-open, as the GAM files
  n_contacts_lost = n_lost, n_contacts_gained = n_gained,
  as.data.table(frac_rep), as.data.table(frac_hmm),
  tad_boundary = tad_boundary, n_tss = as.integer(n_tss),
  compartment_abs = comp_abs, compartment_centred = comp_cen,
  compartment_bin_quintile = comp_binq,
  # domain membership: NA where the bin is not in a called domain
  dom_lost_adj    = dom_out$lost_adj$bin$domain_id,
  dom_lost_raw    = dom_out$lost_raw$bin$domain_id,
  dom_gained_adj  = dom_out$gained_adj$bin$domain_id,
  dom_gained_raw  = dom_out$gained_raw$bin$domain_id,
  smoothed_lost   = dom_out$lost_adj$bin$smoothed,
  smoothed_gained = dom_out$gained_adj$bin$smoothed)
setkey(annot, bin_id)

# --- 8. Sanity assertions -----------------------------------------------------
stopifnot(
  length(unique(annot$bin_id)) == length(gam_bins()),
  nrow(annot) == length(gam_bins()),
  identical(rownames(l2fc), bin_ids),
  mean(in_hot) > 0.9)                       # hotspot files cover ~every mappable bin
for (mk in mark_levels) {
  nz <- mean(rowMeans(signal[, samples$sample_id[samples$mark == mk], drop = FALSE],
                      na.rm = TRUE) > 0, na.rm = TRUE)
  cat(sprintf("%-9s bins with non-zero signal: %.1f%%\n", mk, 100 * nz))
  stopifnot(nz > 0.9)
}
cat(sprintf("Bins with a hotspot row: %.1f%%\n", 100 * mean(in_hot)))

# Printed, not asserted: the 20 kb grid will not reproduce the peak-level number,
# but a sign flip or ~0 means the sizeFactors / exclusion wiring is wrong.
cat(sprintf("\nH4K20me3 per-bin median log2FC = %.3f  (peak-level reference: -0.678 with HP1gKO rep2 excluded)\n",
            med_h4))
cat(sprintf("H3K9me3  per-bin median log2FC = %.3f\n", med_k9))
cat("Compartments (absolute / centred):\n")
print(rbind(abs = table(comp_abs, useNA = "ifany"), centred = table(comp_cen, useNA = "ifany")))

# --- 9. Outputs ----------------------------------------------------------------
summ <- stats[, {
  q <- quantile(log2FoldChange, c(.05, .25, .5, .75, .95), na.rm = TRUE)
  .(n_bins = .N, n_tested = sum(!is.na(log2FoldChange)),
    median_log2FC = q[[3]], q05 = q[[1]], q25 = q[[2]], q75 = q[[4]], q95 = q[[5]],
    frac_p05 = mean(pvalue < 0.05, na.rm = TRUE))
}, by = mark]
fwrite(summ, file.path(tables_dir, "gam_bins_l2fc_summary.tsv"), sep = "\t")
print(summ)

out <- list(
  bins   = {b <- bins; mcols(b) <- DataFrame(bin_id = bin_ids); b},
  stats  = stats,
  l2fc   = l2fc,
  signal = signal,
  annot  = annot,
  meta   = list(
    binsize = gam_binsize,
    n_samples = table(samples$genotype[!duplicated(samples$sample_id)]),
    n_samples_per_mark_genotype = table(samples$mark, samples$genotype),
    excluded_samples = excluded_ids,
    count_scale = COUNT_SCALE,
    cutoffs = c(H4_LOST_CUT = H4_LOST_CUT, K9_LOSS_CUT = K9_LOSS_CUT, UNCH_CUT = UNCH_CUT),
    genome_median_l2fc = c(H4K20me3 = med_h4, H3K9me3 = med_k9),
    date = Sys.Date()))
saveRDS(out, file.path(rds_dir, "gam_bins_annotated.rds"))
cat("Saved:", file.path(rds_dir, "gam_bins_annotated.rds"), "\n")
