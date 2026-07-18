# ================================================================
# MINUTE pipeline — shared configuration & parameters
# ----------------------------------------------------------------
# Sourced by MINUTE_1, MINUTE_2, MINUTE_3 and by run_MINUTE.R.
# SINGLE SOURCE OF TRUTH for parameters — edit values here only;
# the stage scripts must not redefine them.
#
# Run the pipeline FROM the ChIP/ directory (the folder holding
# this file); all repo-relative paths below assume that:
#     cd ChIP && Rscript run_MINUTE.R
#
# I/O layout
#   Inputs (read-only):
#     samples.tsv                sample sheet — SINGLE SOURCE OF TRUTH for samples
#     data/peaks/                consensus / master peak BEDs      (committed, small)
#     $MINUTE_DATA/bigwig/       scaled mm39 bigWigs               (external, large)
#     $MINUTE_DATA/annotation/   LINE / SINE / LTR / TAD BEDs      (external)
#   Outputs (generated, git-ignored except results/rds/):
#     results/counts/  results/rds/  results/tables/  results/figures/  results/bed/
# ================================================================

# --- Roots -------------------------------------------------------
# Large raw inputs (bigWigs + repeat/TAD annotation) are NOT committed. Download
# them from the Google Drive links in the README into the repo-relative folders
# ChIP/data/bigwig and ChIP/data/annotation (the defaults below; scripts run from
# ChIP/). To use a copy staged elsewhere (e.g. a Synology/external mount), set the
# MINUTE_DATA env var:  MINUTE_DATA=/path/to/store Rscript run_MINUTE.R
data_ext   <- path.expand(Sys.getenv("MINUTE_DATA", unset = "data"))
bigwig_dir <- file.path(data_ext, "bigwig")
annot_dir  <- file.path(data_ext, "annotation")

# Repo-relative input / output dirs (scripts are run from ChIP/)
peaks_dir   <- "data/peaks"
results_dir <- "results"
counts_dir  <- file.path(results_dir, "counts")
rds_dir     <- file.path(results_dir, "rds")
tables_dir  <- file.path(results_dir, "tables")
fig_dir     <- file.path(results_dir, "figures")
bed_dir     <- file.path(results_dir, "bed")
for (d in c(counts_dir, rds_dir, tables_dir, fig_dir, bed_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# File written by MINUTE_1, read by MINUTE_2 and MINUTE_3
annotated_rds <- file.path(rds_dir, "annotated_results_H3K9me3_H4K20me3_2000bp_merged.rds")

# --- Libraries (union used across the pipeline) ---
# wigglescout is sensitive to its future/furrr stack. This is the ONLY combination
# confirmed working here; restart R after installing so newer loaded versions unload:
#   remotes::install_version("furrr",   version = "0.2.3")
#   remotes::install_version("future",  version = "1.23.0")
#   remotes::install_version("globals", version = "0.14.0")
#   remotes::install_github("cnluzon/wigglescout")
# If bw_loci() fails with "values must be length 1, but FUN(X[[i]]) result is
# length N", the pins have drifted (or a stray future::plan() is set) — reinstall
# the versions above / run future::plan("sequential"); do NOT edit the R scripts.
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(AnnotationDbi)
  library(TxDb.Mmusculus.UCSC.mm39.knownGene)
  library(org.Mm.eg.db)
  library(GenomicRanges)
  library(rtracklayer)
  library(wigglescout)
  library(parallel)
  library(data.table)
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
})

# --- Sample sheet: SINGLE SOURCE OF TRUTH for samples ------------
# Columns: sample_id, mark, genotype (WT/HP1gKO), replicate, bigwig (filename).
# Genotype and replicate are read from here — NOT hard-coded in the stages.
samples <- as.data.frame(data.table::fread("samples.tsv"))
samples$genotype    <- factor(samples$genotype, levels = c("WT", "HP1gKO"))
samples$bigwig_path <- file.path(bigwig_dir, samples$bigwig)
marks <- unique(as.character(samples$mark))   # preserves sample-sheet order

# --- Associated peak files per mark (under data/peaks) ---
regions <- list(
  "H3K4me3"  = file.path(peaks_dir, "master_peaks_consensus/H3K4me3_consensus_masterPeak.bed"),
  "H3K9me2"  = file.path(peaks_dir, "merged_H3K9me2_masterPeak.bed"),
  "H3K9me3"  = file.path(peaks_dir, "2000bp_merge_H3K9me3_H4K20me3.bed"),
  "H4K20me3" = file.path(peaks_dir, "2000bp_merge_H3K9me3_H4K20me3.bed"),
  "H3K36me3" = file.path(peaks_dir, "merged_H3K36me3_masterPeak.bed")
)

# --- Repeat / TAD annotation files (external, under $MINUTE_DATA/annotation) ---
repeat_bed_files <- list(
  LINE = file.path(annot_dir, "LINE.mm39.bed"),
  SINE = file.path(annot_dir, "SINE.mm39.bed"),
  LTR  = file.path(annot_dir, "LTR.mm39.bed")
)
tad_bed_file <- file.path(annot_dir, "TAD_boundaries_mm39.bed")

# --- Helper: bigWig paths for one mark, named by sample_id ---
# Order follows the sample sheet (WT rep1-6, then HP1gKO rep2-6).
get_bigwig_files <- function(mark) {
  s <- samples[as.character(samples$mark) == mark, , drop = FALSE]
  setNames(s$bigwig_path, s$sample_id)
}

# --- Helper: load a UCSC-style repeat BED into a GRanges ---
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

# --- Helper: load repeat + TAD annotation used by MINUTE_1 and MINUTE_2 ---
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

# --- Significance thresholds: SINGLE SOURCE OF TRUTH ---
# A peak is significant if |log2FoldChange| > lfc AND pvalue < p, per mark.
sig_thresholds <- list(
  "H3K4me3"  = list(lfc = 0.5, p = 0.05),
  "H3K9me2"  = list(lfc = 0.5, p = 0.05),
  "H3K9me3"  = list(lfc = 0.5, p = 0.10),
  "H4K20me3" = list(lfc = 0.5, p = 0.20),
  "H3K36me3" = list(lfc = 0.5, p = 0.20)
)

# Vectorised significance flag for a data.frame that has columns
# ChIP, log2FoldChange, pvalue. Marks not listed above are never significant.
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
