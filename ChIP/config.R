# ================================================================
# MINUTE pipeline - configuration: input files, output tree, parameters
# ----------------------------------------------------------------
# Sourced by every MINUTE stage and by run_MINUTE.R. SINGLE SOURCE OF TRUTH for
# values - edit them here only; stage scripts must not redefine them.
#
#   config.R        values only (this file)
#   samples.tsv     the sample sheet (single source of truth for samples)
#   utils.R         libraries, sample loading, loaders, helper functions
#                   (sourced at the bottom, so stages need only config.R)
#   PARAMETERS.md   why each non-obvious value is set the way it is
#
# Run FROM the ChIP/ directory:   cd ChIP && Rscript run_MINUTE.R
# ================================================================

# ---- Large inputs (NOT committed - download per the README) --------
# MINUTE_DATA=/path/to/store uses a copy staged elsewhere.
data_ext   <- path.expand(Sys.getenv("MINUTE_DATA", unset = "data"))
bigwig_dir <- file.path(data_ext, "bigwig")        # scaled mm39 bigWigs
annot_dir  <- file.path(data_ext, "annotation")    # repeat / TAD / ChromHMM / KAP1 BEDs

# ---- Committed inputs ----------------------------------------------
peaks_dir   <- "data/peaks"                        # consensus / master peak BEDs
samples_tsv <- "samples.tsv"

# ---- Output tree (git-ignored, regenerable) ------------------------
# MINUTE_RESULTS=results_alt writes a sensitivity run to a separate tree.
results_dir <- path.expand(Sys.getenv("MINUTE_RESULTS", unset = "results"))
counts_dir  <- file.path(results_dir, "counts")
rds_dir     <- file.path(results_dir, "rds")
tables_dir  <- file.path(results_dir, "tables")
fig_dir     <- file.path(results_dir, "figures")   # save_fig() reads this by name
bed_dir     <- file.path(results_dir, "bed")
# The handoff written by MINUTE_1 and read by every later stage
annotated_rds <- file.path(rds_dir, "annotated_results_H3K9me3_H4K20me3_2000bp_merged.rds")

# ---- Sample exclusion ------------------------------------------------
# Default: HP1gKO rep2 is excluded (technical failure; PARAMETERS.md).
# MINUTE_EXCLUDE: comma-separated sample_ids and/or <genotype>:<replicate>
# pairs; "none" keeps every sample.
exclude_samples <- trimws(strsplit(Sys.getenv("MINUTE_EXCLUDE", unset = "HP1gKO:2"), ",")[[1]])
if (identical(tolower(exclude_samples), "none")) exclude_samples <- character(0)
exclude_samples <- exclude_samples[nzchar(exclude_samples)]

# ---- Regions quantified per mark (under data/peaks) ----------------
# H3K9me2 uses the shared H3K9me3/H4K20me3 regions, not its own peaks, since
# 2026-11-20 (MACS3 models it poorly). Provenance of every file:
# results/METHODS.md "Peak / region sets".
regions <- list(
  "H3K4me3"  = file.path(peaks_dir, "master_peaks_consensus/H3K4me3_consensus_masterPeak.bed"),
  "H3K9me2"  = file.path(peaks_dir, "2000bp_merge_H3K9me3_H4K20me3.bed"),
  "H3K9me3"  = file.path(peaks_dir, "2000bp_merge_H3K9me3_H4K20me3.bed"),
  "H4K20me3" = file.path(peaks_dir, "2000bp_merge_H3K9me3_H4K20me3.bed"),
  "H3K36me3" = file.path(peaks_dir, "merged_H3K36me3_masterPeak.bed")
)

# ---- Annotation (under $MINUTE_DATA/annotation) ----------------------
repeat_bed_files <- list(                          # UCSC RepeatMasker mm39
  LINE = file.path(annot_dir, "LINE.mm39.bed"),
  SINE = file.path(annot_dir, "SINE.mm39.bed"),
  LTR  = file.path(annot_dir, "LTR.mm39.bed")
)
tad_bed_file <- file.path(annot_dir, "TAD_boundaries_mm39.bed")
chromhmm_bed <- file.path(annot_dir, "ChromHMM_18state_CortexP0_mm39.bed")  # 18-state cortex P0
# KAP1/TRIM28 peaks, lifted from ChIP-Atlas mm10 (README). MINUTE_6 skips if absent.
kap1_neural_bed  <- file.path(annot_dir, "KAP1_Neural_mm39.bed")   # Neuro-2a (GSE110032)
kap1_allcell_bed <- file.path(annot_dir, "KAP1_allCell_mm39.bed")  # ChIP-Atlas all-cell aggregate

# ---- Significance (is_significant) -----------------------------------
# A peak is significant if |log2FoldChange| > lfc AND raw pvalue < p, per mark.
sig_thresholds <- list(
  "H3K4me3"  = list(lfc = 0.5, p = 0.05),
  "H3K9me2"  = list(lfc = 0.5, p = 0.05),
  "H3K9me3"  = list(lfc = 0.5, p = 0.10),
  "H4K20me3" = list(lfc = 0.5, p = 0.20),
  "H3K36me3" = list(lfc = 0.5, p = 0.20)
)

# ---- Gene families silenced by H3K9me3 / HUSH (family_exons) ---------
# Symbol-prefix definitions, quantified over their exons in MINUTE_5.
gene_families <- list(
  Protocadherin = "^Pcdh(a|b|g)",     # clustered protocadherins (chr18 a/b/g clusters)
  `KRAB-ZFP`    = "^Zfp",             # zinc-finger proxy for KRAB-ZFPs (large tandem clusters)
  Vomeronasal   = "^Vmn",            # Vmn1r / Vmn2r vomeronasal receptors
  Olfactory     = "^(Olfr|Or[0-9])"  # olfactory receptors: current (Or*) + legacy (Olfr*) names
)
# Which exon carries the silencing, per family (strand-aware).
family_exon_rule <- list(
  Protocadherin = "first",            # 5'-most: the large variable exon
  `KRAB-ZFP`    = "last",             # 3'-most: the ZNF array
  Vomeronasal   = "largest",          # the single coding exon
  Olfactory     = "largest"
)

# ---- Code ------------------------------------------------------------
source("utils.R")
