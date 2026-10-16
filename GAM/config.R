# ================================================================
# GAM pipeline - shared configuration & parameters
# ----------------------------------------------------------------
# Sourced by every GAM_* stage and by run_GAM.R.
# SINGLE SOURCE OF TRUTH for GAM parameters - edit values here only;
# the stage scripts must not redefine them.
#
# Cross-assay pieces (colour palettes, figure savers, hyper_row, size_perm,
# seqname flips, mm39_bins, GAM bin-id parsing) live in ../shared/common.R.
#
# Run the pipeline FROM the GAM/ directory (the folder holding this file):
#     cd GAM && Rscript run_GAM.R
#
# I/O layout
#   Inputs (read-only, NOT committed - ~7.7 GB, see README):
#     $GAM_DATA/NPMI_matrices/{WT,cKO}_at20Kb/   dense square NPMI matrices
#     $GAM_DATA/Diff_contacts/                   .bg2 differential contact pairs
#     $GAM_DATA/Hotspots/                        per-bin differential contact counts
#   Input from the ChIP assay (produced by ChIP/run_MINUTE.R):
#     ../ChIP/results/rds/annotated_results_*.rds
#   Outputs (generated, git-ignored in full - all regenerable by run_GAM.R):
#     results/counts/  results/rds/  results/tables/  results/figures/  results/bed/
# ================================================================

# --- Dependencies + cross-assay helpers --------------------------
source("../shared/setup.R")
# GAM-only dependency: every GAM input is gzipped, and data.table::fread()
# delegates .gz decompression to R.utils. Without it every loader below fails
# with "fread() requires 'R.utils'". The ChIP assay reads plain BEDs/bigWigs, so
# this belongs to GAM rather than the shared manifest.
# RSpectra: GAM_2 gets PC1 from the correlation OPERATOR via eigs_sym() without
# ever materialising the n x n correlation matrix (~760 MB on chr1) or paying an
# O(n^3) eigen(). It has a power-iteration fallback, so without this entry a
# fresh clone still runs but silently takes the slower, less accurate path -
# which is exactly the kind of difference that should not depend on what happens
# to be installed.
minute_add_deps(cran = c("R.utils", "RSpectra", "hexbin"))
minute_ensure_deps()
source("../shared/common.R")

# GAM-SPECIFIC libraries. The core set (GenomicRanges, GenomeInfoDb,
# rtracklayer, data.table, ggplot2, circlize, ltc) comes from shared/common.R.
suppressPackageStartupMessages({
  library(TxDb.Mmusculus.UCSC.mm39.knownGene)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(wigglescout)   # bw_loci(): ChIP signal over the GAM bin grid
  library(DESeq2)
  library(dplyr)
  library(parallel)
  library(ComplexHeatmap)
})

# --- Roots -------------------------------------------------------
# The GAM data ships inside this folder, so GAM_DATA defaults to "." . Set it to
# use a copy staged elsewhere (e.g. an external mount):
#     GAM_DATA=/Volumes/store/GAM Rscript run_GAM.R
gam_data   <- path.expand(Sys.getenv("GAM_DATA", unset = "."))
npmi_dir   <- file.path(gam_data, "NPMI_matrices")
diff_dir   <- file.path(gam_data, "Diff_contacts")
hotspot_dir<- file.path(gam_data, "Hotspots")

# results_dir is overridable so a sensitivity run can be written to a SEPARATE
# tree without disturbing the committed run:
#     GAM_RESULTS=results_alt Rscript run_GAM.R
results_dir <- path.expand(Sys.getenv("GAM_RESULTS", unset = "results"))
counts_dir  <- file.path(results_dir, "counts")
rds_dir     <- file.path(results_dir, "rds")
tables_dir  <- file.path(results_dir, "tables")
fig_dir     <- file.path(results_dir, "figures")   # save_fig() reads this by name
bed_dir     <- file.path(results_dir, "bed")
for (d in c(counts_dir, rds_dir, tables_dir, fig_dir, bed_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# --- GAM parameters: these MATCH THE FILENAMES on disk -----------
# Changing one here changes which files the path builders below resolve to, so a
# re-run at a different sigma / z / NPMI cutoff is a config edit, not a code edit.
gam_binsize        <- 20000L     # segregation_at20000
gam_sigma          <- "2.6"      # sigma.2.6
gam_seg_date       <- "2026-02-23"
gam_seg_prefix     <- "HP1g_mm39"
gam_z              <- 2          # z2       - differential contact z-score cutoff
gam_npmi           <- 0.1        # npmi0.1  - minimum NPMI in the selected genotype
gam_maxdist_label  <- "inf"      # maxdistinf
gam_hotspot_maxdist<- 5e6        # dist0-5000000 - window for the hotspot counts

# --- Genotype labels: THE ONE PLACE THE BRIDGE LIVES -------------
# Canonical repo-wide labels are WT / HP1gKO (see ../CLAUDE.md). The GAM files
# say "cKO" for HP1gKO. Everything downstream of this file uses the canonical
# name; only the path builders below ever see "cKO".
#
# GOTCHA: the NPMI filenames are IDENTICAL between the two genotype folders -
# genotype is encoded ONLY by the directory name (WT_at20Kb / cKO_at20Kb). A
# matrix file copied out of its folder has lost its genotype, so never resolve
# one by basename.
gam_genotypes <- c(WT = "WT", HP1gKO = "cKO")

# --- Chromosomes -------------------------------------------------
# The differential products (Diff_contacts, Hotspots) are AUTOSOME-ONLY.
# The NPMI matrices additionally carry chrX and chrY. Keep cross-assay
# integration on the autosomes; X/Y are available for NPMI-only analyses.
gam_autosomes   <- paste0("chr", 1:19)
gam_npmi_chroms <- c(gam_autosomes, "chrX", "chrY")

# --- Direction semantics: READ THIS BEFORE USING THE DIFF FILES --
# The upstream caller labels each significant differential contact by the
# genotype in which it is SELECTED (i.e. stronger):
#     selWT  -> contact stronger in WT   -> LOST   in HP1gKO
#     selcKO -> contact stronger in cKO  -> GAINED in HP1gKO
# Because the whole study is framed as "what does losing HP1g do?", every
# downstream column is named from the KO's point of view: n_contacts_lost /
# n_contacts_gained. This mapping is the single most reversible mistake in the
# GAM integration - getting it backwards inverts the headline result.
gam_sel <- c(lost = "selWT", gained = "selcKO")

# --- Path builders -----------------------------------------------
npmi_file <- function(genotype, chr) {
  genotype <- match.arg(genotype, names(gam_genotypes))
  file.path(npmi_dir, paste0(gam_genotypes[[genotype]], "_at", gam_binsize / 1000, "Kb"),
            sprintf("%s.segregation_at%d.%s.passed_qc.curated.sigma.%s.NPMI.%s.squared.txt.gz",
                    gam_seg_prefix, gam_binsize, gam_seg_date, gam_sigma, chr))
}
diff_file <- function(chr, sel) {
  sel <- if (sel %in% names(gam_sel)) gam_sel[[sel]] else sel
  file.path(diff_dir, sprintf("significant_WT_vs_cKO_%s_%s_maxdist%s_z%s_npmi%s.bg2.gz",
                              chr, sel, gam_maxdist_label, gam_z, gam_npmi))
}
hotspot_file <- function(chr, sel) {
  sel <- if (sel %in% names(gam_sel)) gam_sel[[sel]] else sel
  file.path(hotspot_dir,
            sprintf("significant_WT_vs_cKO_%s_%s_maxdist%s_z%s_npmi%s_dist0-%d.bg.gz",
                    chr, sel, gam_maxdist_label, gam_z, gam_npmi, as.integer(gam_hotspot_maxdist)))
}

# --- Handoff from the ChIP assay ---------------------------------
# Produced by ChIP/run_MINUTE.R (stage 1). Override the ChIP results tree with
# CHIP_RESULTS if the ChIP run was written elsewhere (MINUTE_RESULTS=...).
chip_results_dir <- path.expand(Sys.getenv("CHIP_RESULTS", unset = "../ChIP/results"))
chip_rds <- file.path(chip_results_dir, "rds",
                      "annotated_results_H3K9me3_H4K20me3_2000bp_merged.rds")

require_chip_rds <- function() {
  if (!file.exists(chip_rds))
    stop("ChIP handoff not found:\n  ", chip_rds,
         "\nThe GAM stages need the annotated peak-level results from the ChIP",
         " pipeline.\nRun it first:   cd ../ChIP && Rscript run_MINUTE.R",
         "\nOr point CHIP_RESULTS at an existing ChIP results tree.", call. = FALSE)
  chip_rds
}

# --- ChIP inputs that GAM re-quantifies --------------------------
# GAM_1 puts ChIP signal on the GAM bin grid, so it needs the ChIP bigWigs, the
# ChIP annotation BEDs and the ChIP sample sheet.
#
# WHY NOT source("../ChIP/config.R"): both configs define results_dir /
# counts_dir / rds_dir / fig_dir / tables_dir by the SAME names, so sourcing the
# ChIP one here would silently redirect every GAM output into ChIP/results (and
# re-attach the whole ChIP library set). The few things GAM needs are therefore
# re-derived below, and they MUST track ChIP/config.R: same MINUTE_DATA root,
# same sample sheet, same default exclusion.
minute_data <- path.expand(Sys.getenv("MINUTE_DATA", unset = "../ChIP/data"))
bigwig_dir  <- file.path(minute_data, "bigwig")
annot_dir   <- file.path(minute_data, "annotation")
chip_samples_tsv <- "../ChIP/samples.tsv"

# --- Sample sheet + exclusion (mirrors ChIP/config.R) ------------
# DEFAULT: HP1gKO rep2 is EXCLUDED, exactly as on the ChIP side - it is the
# highest-signal library in all five marks and discordant with rep3, its own
# technical twin, so it is treated as a technical failure. Keeping the same
# default is what makes the per-bin log2FC comparable to the published
# peak-level numbers; MINUTE_EXCLUDE=none reproduces the all-library run.
# Syntax: comma-separated sample_ids and/or <genotype>:<replicate> pairs.
gam_load_samples <- function() {
  sm <- as.data.frame(data.table::fread(chip_samples_tsv))
  ex <- trimws(strsplit(Sys.getenv("MINUTE_EXCLUDE", unset = "HP1gKO:2"), ",")[[1]])
  if (identical(tolower(ex), "none")) ex <- character(0)
  ex <- ex[nzchar(ex)]
  excluded <- character(0)
  if (length(ex)) {
    pairs <- grepl(":", ex, fixed = TRUE)
    drop  <- sm$sample_id %in% ex[!pairs]
    for (q in ex[pairs]) {
      gr <- strsplit(q, ":", fixed = TRUE)[[1]]
      drop <- drop | (as.character(sm$genotype) == gr[1] &
                      as.character(sm$replicate) == gr[2])
    }
    excluded <- sm$sample_id[drop]
    message("Excluding ", sum(drop), " sample(s): ", paste(excluded, collapse = ", "))
    sm <- sm[!drop, , drop = FALSE]
  }
  sm$genotype    <- factor(sm$genotype, levels = c("WT", "HP1gKO"))
  sm$bigwig_path <- file.path(bigwig_dir, sm$bigwig)
  attr(sm, "excluded_ids") <- excluded
  sm
}

# --- Quantification parameters -----------------------------------
# COUNT_SCALE: DESeq2 needs integers, but mean coverage over a 20 kb bin is
# often below 1, so plain rounding would collapse most bins onto 0/1. Scaling by
# 10 before rounding preserves the gradation; the factor cancels in log2FC, and
# baseMean is divided back out so it stays on the coverage scale. Same value and
# reasoning as ChIP/MINUTE_7_H3K36me3_Genes.R.
gam_count_scale <- 10

# Cores for the per-sample bigWig reads. Leaves one core free; GAM_1 is a
# sustained all-core job on a laptop, so lower this (GAM_CORES=4) when running
# it alongside anything else.
gam_cores <- as.integer(Sys.getenv("GAM_CORES", unset = max(1, parallel::detectCores() - 1)))

# --- Compartment cutoffs: MINUTE_3's, VERBATIM -------------------
# Do not retune these here. They define the stable / H4K20me3_only / co_loss
# groups, and the whole point of reusing them is that the GAM bins are labelled
# by the SAME rule as the ChIP peaks. Note they are ABSOLUTE and therefore not
# shift-invariant against a genome-wide loss, which is why every GAM stage
# treats the median-CENTRED variant as primary (ChIP/results/COMPARTMENT_CUTOFFS.md).
H4_LOST_CUT <- -0.5    # H4K20me3 log2FC below this = "H4K20me3-lost"
K9_LOSS_CUT <- -0.3    # H3K9me3  log2FC below this = "H3K9me3-loss"
UNCH_CUT    <-  0.15   # |log2FC| below this = "unchanged"

# --- The shared coordinate system --------------------------------
# The 20 kb grid that BOTH assays are projected onto. mm39_bins() comes from
# shared/common.R; UCSC seqnames, autosomes by default (matching the GAM
# differential products). This is the join key for the whole integration.
gam_bins <- function(chroms = gam_autosomes) mm39_bins(gam_binsize, chroms)

# --- Loader: per-bin differential contact counts (BOTH directions) ----------
# Reads all Hotspots files and returns ONE table on the GAM bin grid:
#   chr, start, end, n_contacts_lost, n_contacts_gained
# Coordinates are the file's own (0-based half-open, BED-like) - converted to
# 1-based GRanges only via parse_gam_bin()/makeGRangesFromDataFrame downstream.
#
# IMPORTANT: these files are DENSE, not a significant-vs-not hit list - ~115k
# bins per direction is essentially every mappable 20 kb autosomal bin. The
# signal is the COUNT MAGNITUDE (observed range 1-162), so treat it as a
# continuous per-bin contact-remodelling density. Bins absent from a file have
# no significant differential contacts, so they are filled with 0, not NA.
load_hotspots <- function(chroms = gam_autosomes) {
  one <- function(sel) {
    rbindlist(lapply(chroms, function(ch) {
      f <- hotspot_file(ch, sel)
      if (!file.exists(f)) stop("missing hotspot file: ", f)
      dt <- fread(f, col.names = c("chr", "start", "end", "count"))
      dt
    }))
  }
  lost   <- one("lost")[,   .(chr, start, end, n_contacts_lost   = count)]
  gained <- one("gained")[, .(chr, start, end, n_contacts_gained = count)]
  out <- merge(lost, gained, by = c("chr", "start", "end"), all = TRUE)
  out[is.na(n_contacts_lost),   n_contacts_lost   := 0L]
  out[is.na(n_contacts_gained), n_contacts_gained := 0L]
  setorder(out, chr, start)
  out[]
}

# --- Loader: differential contact PAIRS for one chromosome/direction --------
# ~25.6M pairs genome-wide per direction, so this is deliberately per-chromosome
# and never loads everything at once. `sel` is "lost"/"gained" (canonical) or a
# raw "selWT"/"selcKO". Columns are the file's own; see the README for the full
# list (set1_npmi/set2_npmi are WT/cKO NPMI, deltaz is the z-score difference).
load_diff_contacts <- function(chr, sel, select = NULL) {
  f <- diff_file(chr, sel)
  if (!file.exists(f)) stop("missing differential contact file: ", f)
  fread(f, select = select)
}

# --- Loader: one dense NPMI matrix -------------------------------------------
# Returns a numeric matrix with GAM bin-id dimnames ("chr19:0-20000"), which
# parse_gam_bin() turns into a GRanges. ONE CHROMOSOME AT A TIME, on purpose:
# chr1 is ~9.7k x 9.7k, about 760 MB as doubles, so holding both genotypes for a
# large chromosome at once is how this stage runs the machine out of memory.
# Diagonal/unmappable cells are NA in the source and stay NA here.
load_npmi <- function(genotype, chr) {
  f <- npmi_file(genotype, chr)
  if (!file.exists(f)) stop("missing NPMI matrix: ", f)
  dt <- fread(f, header = TRUE)
  ids <- as.character(dt[[1]])
  m <- as.matrix(dt[, -1, with = FALSE])
  rownames(m) <- ids
  if (!identical(colnames(m), ids))
    warning("NPMI matrix ", basename(f), ": row and column ids differ in order")
  m
}
