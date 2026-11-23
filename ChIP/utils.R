# ================================================================
# MINUTE pipeline - libraries, sample sheet loading, loaders and helpers
# ----------------------------------------------------------------
# Sourced by config.R (last line), so a stage only ever needs
#     source("config.R")
# This file holds CODE ONLY. Every value it uses (paths, peak files, cutoffs,
# gene-family definitions) is defined in config.R; nothing here should be edited
# to change an analysis. Cross-assay helpers (palettes + scales, save_fig /
# save_base_fig, ht_raster_device, hyper_row, size_perm, seqname flips,
# mm39_bins, GAM bin-id parsing) live in ../shared/common.R.
# ================================================================

# --- Libraries (union used across the pipeline) ---
# setup.R holds the dependency MANIFEST (CRAN / Bioconductor / GitHub) and the
# installer. Sourcing it only defines functions; the call below installs whatever
# is missing, so a fresh clone works on the first `Rscript run_MINUTE.R`.
# Set MINUTE_AUTO_INSTALL=0 to turn the auto-install off (missing packages then
# raise an error listing them, and `Rscript setup.R` installs them on demand).
source("setup.R")        # shim -> ../shared/setup.R (manifest + installer)
minute_ensure_deps()

# Cross-assay helpers: core library attachments, colour palettes + scales,
# save_fig()/save_base_fig(), ht_raster_device, hyper_row(), size_perm(),
# seqname flips (to_ucsc/to_ncbi), mm39_bins(), GAM bin-id parsing.
# SINGLE SOURCE OF TRUTH for figure colours - stage scripts must not write hex.
source("../shared/common.R")

# wigglescout <= 0.20 was sensitive to its future/furrr stack and needed pins
# (furrr 0.2.3 / future 1.23.0 / globals 0.14.0), restarting R afterwards so newer
# loaded versions unload. From 0.21 it no longer depends on future/furrr at all,
# and setup.R installs the current release. On an OLD wigglescout, if bw_loci()
# fails with "values must be length 1, but FUN(X[[i]]) result is length N", the
# pins have drifted (or a stray future::plan() is set) - reinstall those versions
# / run future::plan("sequential"); do NOT edit the R scripts.
# ChIP-SPECIFIC libraries. The core set (GenomicRanges, GenomeInfoDb,
# rtracklayer, data.table, ggplot2, circlize, ltc) is attached by
# shared/common.R above - do not re-attach it here.
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(AnnotationDbi)
  library(TxDb.Mmusculus.UCSC.mm39.knownGene)
  library(org.Mm.eg.db)
  library(wigglescout)
  library(parallel)
  library(DESeq2)
  library(dplyr)
  library(ggrepel)
  library(ComplexHeatmap)
})

# --- Output tree -------------------------------------------------
for (d in c(counts_dir, rds_dir, tables_dir, fig_dir, bed_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# --- Sample sheet, minus the excluded samples --------------------
# Builds `samples` (sample_id, mark, genotype, replicate, bigwig, bigwig_path)
# and `marks` from samples_tsv and exclude_samples (both in config.R). Genotype
# and replicate come from the sheet - NOT hard-coded in the stages.
# exclude_samples: sample_id values and/or "<genotype>:<replicate>" pairs, which
# expand across every mark.
samples <- local({
  s <- as.data.frame(data.table::fread(samples_tsv))
  if (length(exclude_samples)) {
    pairs <- grepl(":", exclude_samples, fixed = TRUE)
    drop  <- s$sample_id %in% exclude_samples[!pairs]
    for (p in exclude_samples[pairs]) {
      gr <- strsplit(p, ":", fixed = TRUE)[[1]]
      drop <- drop | (as.character(s$genotype) == gr[1] &
                      as.character(s$replicate) == gr[2])
    }
    message("Excluding ", sum(drop), " sample(s): ",
            paste(s$sample_id[drop], collapse = ", "))
    s <- s[!drop, , drop = FALSE]
  }
  s$genotype    <- factor(s$genotype, levels = c("WT", "HP1gKO"))
  s$bigwig_path <- file.path(bigwig_dir, s$bigwig)
  s
})
marks <- unique(as.character(samples$mark))   # preserves sample-sheet order

# --- bigWig paths for one mark, named by sample_id ---
# Order follows the sample sheet (WT rep1-6, then HP1gKO rep2-6).
get_bigwig_files <- function(mark) {
  s <- samples[as.character(samples$mark) == mark, , drop = FALSE]
  setNames(s$bigwig_path, s$sample_id)
}

# --- Significance: the one definition every stage uses ---
# Vectorised significance flag for a data.frame that has columns
# ChIP, log2FoldChange, pvalue, against sig_thresholds (config.R).
# Marks not listed there are never significant.
# NA log2FC/pvalue -> FALSE (matches the previous dplyr::filter behaviour).
is_significant <- function(df) {
  lfc_cut <- vapply(as.character(df$ChIP),
                    function(m) if (!is.null(sig_thresholds[[m]])) sig_thresholds[[m]]$lfc else Inf,
                    numeric(1))
  p_cut   <- vapply(as.character(df$ChIP),
                    function(m) if (!is.null(sig_thresholds[[m]])) sig_thresholds[[m]]$p else -Inf,
                    numeric(1))
  ok <- !is.na(df$log2FoldChange) & !is.na(df$pvalue) &
        abs(df$log2FoldChange) > lfc_cut & df$pvalue < p_cut
  unname(ok)
}

# --- ChromHMM ----------------------------------------------------
# Load the segmentation as a GRanges (state in $state). The file has a few
# overlapping same-state segments, which would double-count coverage; merge
# overlapping segments within each state so per-state coverage is well defined
# (fraction of a peak in any one state cannot exceed 1).
load_chromHMM <- function() {
  gr <- rtracklayer::import(chromhmm_bed, format = "bed")   # 4th BED col -> $name
  mcols(gr)$state <- as.character(mcols(gr)$name)
  red <- unlist(reduce(split(gr, mcols(gr)$state)))
  mcols(red)$state <- names(red)
  names(red) <- NULL
  red
}

# Per-peak coverage FRACTION of each ChromHMM state. Domains span many states,
# so we keep the full breakdown rather than a single size-biased label.
# gr: peaks GRanges (any seqnames style). Returns a matrix [peaks x states]; each
# value = fraction of that peak's width covered by the state (rows sum to <= 1).
# Downstream: dominant label + coverage-based enrichment.
chromHMM_coverage <- function(gr, chromhmm_gr) {
  suppressWarnings(seqlevelsStyle(chromhmm_gr) <- seqlevelsStyle(gr)[1])
  states <- sort(unique(as.character(mcols(chromhmm_gr)$state)))
  m <- matrix(0, nrow = length(gr), ncol = length(states),
              dimnames = list(NULL, states))
  hits <- findOverlaps(gr, chromhmm_gr, ignore.strand = TRUE)
  if (length(hits) > 0) {
    ov  <- pintersect(gr[queryHits(hits)], chromhmm_gr[subjectHits(hits)],
                      ignore.strand = TRUE)
    agg <- data.table(q = queryHits(hits),
                      s = as.character(mcols(chromhmm_gr)$state[subjectHits(hits)]),
                      w = width(ov))[, .(w = sum(w)), by = .(q, s)]
    m[cbind(agg$q, match(agg$s, states))] <- agg$w
    m <- pmin(m / width(gr), 1)   # fraction; cap at 1 (rare cross-state overlaps)
  }
  m
}

# SIZE-WEIGHTED per-state coverage enrichment + permutation test. Complements
# the per-region (unweighted) Wilcoxon tests by weighting each region by its
# length, so the result reflects the changed genomic TERRITORY rather than
# treating a 2 Mb domain and a 2 kb peak equally.
#   cov_mat : regions x states coverage fractions (the hmm_<state> columns)
#   size    : region length (bp or kb - only relative weights matter)
#   group   : logical, TRUE = foreground (e.g. significant / cluster / H3K9me3-loss)
# Returns per state: weighted mean coverage in fg vs bg (= fraction of that
# group's territory in the state), log2 ratio, and a permutation p-value that
# shuffles region LABELS (domain-level; a bp-level test would pseudo-replicate).
chromHMM_size_weighted <- function(cov_mat, size, group, nperm = 1000, seed = 42) {
  cov_mat <- as.matrix(cov_mat)
  keep    <- is.finite(size) & !is.na(group)
  cov_mat <- cov_mat[keep, , drop = FALSE]; size <- size[keep]; group <- as.logical(group)[keep]
  covW <- cov_mat * size                       # region x state, bp-in-state
  totW <- sum(size)
  colW <- colSums(covW, na.rm = TRUE)          # total state-bp over all regions
  gv   <- as.numeric(group)
  fgbp <- as.numeric(crossprod(gv, covW))      # state-bp in foreground
  sf   <- sum(gv * size); sb <- totW - sf
  fg   <- fgbp / sf; bg <- (colW - fgbp) / sb
  obsd <- fg - bg
  set.seed(seed)
  n <- length(group); k <- sum(group); cnt <- numeric(ncol(cov_mat))
  for (i in seq_len(nperm)) {
    g <- numeric(n); g[sample.int(n, k)] <- 1
    fb <- as.numeric(crossprod(g, covW)); sff <- sum(g * size)
    cnt <- cnt + (abs(fb / sff - (colW - fb) / (totW - sff)) >= abs(obsd))
  }
  data.frame(state = sub("^hmm_", "", colnames(cov_mat)),
             w_cov_fg = fg, w_cov_bg = bg,
             w_log2_ratio = log2((fg + 1e-4) / (bg + 1e-4)),
             perm_p = (cnt + 1) / (nperm + 1),
             stringsAsFactors = FALSE)
}

# Dominant state + its purity from a coverage matrix. Convenience single label
# = the state covering the largest fraction of the peak, plus that fraction
# (purity). NA where a peak overlaps no ChromHMM segment.
chromHMM_dominant <- function(cov) {
  rs <- rowSums(cov)
  mc <- max.col(cov, ties.method = "first")
  data.frame(
    chromHMM_state  = ifelse(rs > 0, colnames(cov)[mc], NA_character_),
    chromHMM_purity = ifelse(rs > 0, cov[cbind(seq_len(nrow(cov)), mc)], NA_real_),
    stringsAsFactors = FALSE
  )
}

# --- UN-ROUNDED mean coverage per peak (a PLOTTING covariate) ----------------
# MINUTE_1 rounds the per-peak mean coverage to integers to satisfy DESeq2, which
# collapses `baseMean` onto a coarse k/n_samples lattice - for H4K20me3 that is
# only ~214 distinct values across 116k peaks, with one value shared by >9,000
# peaks. That is far too quantised to encode as a point size.
#
# The counts TSV is written BEFORE that rounding, so it keeps full precision.
# Because sizeFactors = 1, the row mean of that file is EXACTLY what `baseMean`
# would have been had we never rounded - the same quantity, un-degraded.
#
# Use this for AESTHETICS ONLY (point size, colour). All statistics must keep
# using DESeq2's baseMean/log2FoldChange/pvalue, which are computed from the
# rounded integer counts.
# Returns a named numeric vector (peak_id -> mean scaled coverage), or NULL if
# the counts file is missing.
mean_coverage_by_peak <- function(mark) {
  f <- file.path(counts_dir, paste0(mark, "_bigwig_counts.tsv"))
  if (!file.exists(f)) {
    warning("mean_coverage_by_peak(): no counts file for ", mark, " - ", f)
    return(NULL)
  }
  dt  <- data.table::fread(f)
  ids <- as.character(dt[[1]])                       # peak_id column
  m   <- as.matrix(dt[, -1, with = FALSE])
  stats::setNames(rowMeans(m, na.rm = TRUE), ids)
}

# --- Repeat + TAD annotation -------------------------------------
# Load a UCSC-style repeat BED (RepeatMasker table) into a GRanges.
loadRepeatBED <- function(filepath) {
  expected_colnames <- c("bin", "swScore", "milliDiv", "milliDel", "milliIns",
                         "genoName", "genoStart", "genoEnd", "genoLeft", "strand",
                         "repName", "repClass", "repFamily", "repStart", "repEnd",
                         "repLeft", "id")
  df <- read.table(filepath, sep = "\t", header = FALSE, comment.char = "#", quote = "",
                   fill = TRUE, stringsAsFactors = FALSE,
                   col.names = expected_colnames, colClasses = rep("character", 17))
  df$genoStart <- suppressWarnings(as.numeric(df$genoStart))
  df$genoEnd   <- suppressWarnings(as.numeric(df$genoEnd))
  df <- df[!is.na(df$genoStart) & !is.na(df$genoEnd), ]
  GRanges(
    seqnames = df$genoName,
    ranges   = IRanges(start = df$genoStart + 1, end = df$genoEnd),
    strand   = df$strand,
    repName  = df$repName,
    repClass = df$repClass,
    repFamily = df$repFamily
  )
}

# Repeat + TAD annotation used by MINUTE_1 and MINUTE_2.
# Returns a list; the stage scripts assign to line_gr/sine_gr/ltr_gr/tad_gr.
load_annotation <- function() {
  ann <- list(
    line = loadRepeatBED(repeat_bed_files$LINE),
    sine = loadRepeatBED(repeat_bed_files$SINE),
    ltr  = loadRepeatBED(repeat_bed_files$LTR),
    tad  = import(tad_bed_file, format = "bed")
  )
  seqlevelsStyle(ann$tad) <- "UCSC"
  ann
}

# --- One silencing-relevant exon PER GENE for each family (mm39 TxDb) ---
# Returns a GRanges with one interval per family member gene (the exon selected by
# family_exon_rule, strand-aware), with `family` and `gene` metadata. Seqnames are
# UCSC; set the style at query time. Families and rules: config.R.
family_exons <- function() {
  txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
  exg  <- exonsBy(txdb, by = "gene")                       # by Entrez ID
  sym  <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, keys = names(exg),
            column = "SYMBOL", keytype = "ENTREZID"))
  sym  <- ifelse(is.na(sym), "", sym)
  pick <- function(ex, rule) {                             # ex = one gene's exons
    if (rule == "largest") return(ex[which.max(width(ex))])
    ex <- ex[if (as.character(strand(ex))[1] == "-") order(-start(ex)) else order(start(ex))]
    if (rule == "first") ex[1] else ex[length(ex)]         # "first" = 5', else "last" = 3'
  }
  parts <- lapply(names(gene_families), function(fn) {
    idx  <- which(grepl(gene_families[[fn]], sym))
    if (!length(idx)) return(NULL)
    rule <- family_exon_rule[[fn]]
    sel  <- lapply(idx, function(i) {
      g <- pick(exg[[i]], rule); mcols(g) <- NULL
      mcols(g)$family <- fn; mcols(g)$gene <- sym[i]; g
    })
    do.call(c, sel)
  })
  do.call(c, Filter(Negate(is.null), parts))
}
