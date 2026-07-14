# ================================================================
# MINUTE pipeline — shared configuration & parameters
# ----------------------------------------------------------------
# Sourced by MINUTE_1, MINUTE_2, MINUTE_3 and by run_MINUTE.R.
# This is the SINGLE SOURCE OF TRUTH for parameters — edit values
# here only; the stage scripts must not redefine them.
# ================================================================

# --- Working directory & paths ---
setwd("~/SynologyDrive/MINUTE/")
bigwigDir <- "~/SynologyDrive/MINUTE/bigwig/"
outputDir <- "counts"
if (!dir.exists(outputDir)) dir.create(outputDir)

# File written by MINUTE_1, read by MINUTE_2 and MINUTE_3
annotated_rds <- "annotated_results_H3K9me3_H4K20me3_2000bp_merged.rds"

# --- Libraries (union used across the pipeline) ---
# For wigglescout you may need these pinned versions of furrr & future:
#   remotes::install_version("furrr",   version = "0.2.3")
#   remotes::install_version("future",  version = "1.23.0")
#   remotes::install_version("globals", version = "0.14.0")
#   remotes::install_github("cnluzon/wigglescout")
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

# --- ChIP marks and corresponding BigWig file prefixes ---
chips <- list(
  "H3K4me3"  = "A4435_2_H3K4me3_S2_L008",
  "H3K9me2"  = "A4435_3_H3K9me2_S3_L008",
  "H3K9me3"  = "A4435_4_H3K9me3_S4_L008",
  "H4K20me3" = "A4435_5_H4K20me3_S5_L008",
  "H3K36me3" = "A4435_6_H3K36me3_S6_L008"
)

# --- Associated peak files per mark ---
regions <- list(
  "H3K4me3"  = "peaks/master_peaks_consensus/H3K4me3_consensus_masterPeak.bed",
  "H3K9me2"  = "peaks/merged_H3K9me2_masterPeak.bed",
  "H3K9me3"  = "peaks/2000bp_merge_H3K9me3_H4K20me3.bed",
  "H4K20me3" = "peaks/2000bp_merge_H3K9me3_H4K20me3.bed",
  "H3K36me3" = "peaks/merged_H3K36me3_masterPeak.bed"
)

# --- Repeat / TAD annotation files ---
repeat_bed_files <- list(LINE = "LINE.mm39.bed", SINE = "SINE.mm39.bed", LTR = "LTR.mm39.bed")
tad_bed_file     <- "TAD_boundaries_mm39.bed"

# --- Helper: build the 11 BigWig paths for a given mark prefix ---
get_bigwig_files <- function(prefix, dir = bigwigDir) {
  names <- c(
    "WT_1", "WT_2", "WT_3", "WT_4", "WT_5", "WT_6",
    "HP1gKO_2", "HP1gKO_3", "HP1gKO_4", "HP1gKO_5", "HP1gKO_6"
  )
  suffixes <- c(
    "_WT_rep1.mm39.scaled.bw", "_WT_rep2.mm39.scaled.bw", "_WT_rep3.mm39.scaled.bw",
    "_WT_rep4.mm39.scaled.bw", "_WT_rep5.mm39.scaled.bw", "_WT_rep6.mm39.scaled.bw",
    "_HP1gKO_rep2.mm39.scaled.bw", "_HP1gKO_rep3.mm39.scaled.bw", "_HP1gKO_rep4.mm39.scaled.bw",
    "_HP1gKO_rep5.mm39.scaled.bw", "_HP1gKO_rep6.mm39.scaled.bw"
  )
  setNames(paste0(dir, prefix, suffixes), names)
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
