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
minute_add_deps(cran = c("R.utils", "RSpectra", "hexbin", "patchwork"))
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
# The consensus peak BEDs are committed in the ChIP tree (small), unlike the
# bigWigs and annotation. GAM_10 reads the merged H3K9me3/H4K20me3 set from here
# to define heterochromatin domains.
chip_peaks_dir <- "../ChIP/data/peaks"
het_domain_bed <- file.path(chip_peaks_dir, "2000bp_merge_H3K9me3_H4K20me3.bed")

# Lamina-associated domains, Peric-Hupkes 2010 (GSE17051), lifted mm9 -> mm39 by
# prepare_LAD_annotation.R. CONSTITUTIVE sets: cLAD is LAD in all four profiled
# cell types, ciLAD is LAD in none. They are deliberately NOT neuronal - they are
# the cell-type-invariant core, which is what makes them an axis independent of
# our own marks and of either genotype. Per-cell-type (NPC, astrocyte) calls
# exist in the same GEO entry if cell-type specificity ever matters.
lad_bed   <- file.path(annot_dir, "cLAD_mm39.bed")
cilad_bed <- file.path(annot_dir, "ciLAD_mm39.bed")
load_lads <- function() {
  if (!file.exists(lad_bed) || !file.exists(cilad_bed))
    stop("LAD BEDs missing - run: Rscript prepare_LAD_annotation.R", call. = FALSE)
  # One GRanges carrying a `class` column, not a list: every consumer so far
  # wants to overlap both sets against the same bins and keep the label.
  a <- rtracklayer::import(lad_bed);   a$class <- "cLAD"
  b <- rtracklayer::import(cilad_bed); b$class <- "ciLAD"
  g <- c(a, b)
  g$class <- factor(g$class, levels = c("cLAD", "ciLAD"))
  g
}

# --- TAD partitions -----------------------------------------------------------
# Two independent partitions, kept side by side on purpose. Any "this locus talks
# to distal TADs differently" claim depends on where the boundaries are, so a
# result that only holds under one partition is a result about the partition.
#
#   EXTERNAL  published boundaries, fixed across genotypes. Non-circular, but a
#             different cell type and blind to anything the KO does.
#   OWN       insulation minima from our own GAM, per genotype. Self-consistent
#             with the contact data and free to move between genotypes, but it is
#             our own call, so it cannot also be evidence about boundary strength.
tad_boundary_bed <- file.path(annot_dir, "TAD_boundaries_mm39.bed")

# Intervals are the spans BETWEEN consecutive boundaries. The deposited
# boundaries are ~150 kb wide, so a TAD runs from the end of one to the start of
# the next; the boundary bodies themselves belong to neither TAD.
load_tads_external <- function(chroms = gam_autosomes) {
  if (!file.exists(tad_boundary_bed))
    stop("missing TAD boundaries: ", tad_boundary_bed, call. = FALSE)
  b <- data.table::fread(tad_boundary_bed, header = FALSE, select = 1:3,
                         col.names = c("chr", "start", "end"))
  b[, chr := ifelse(grepl("^chr", chr), chr, paste0("chr", chr))]
  b <- b[chr %in% chroms]
  data.table::setorder(b, chr, start)
  out <- b[, {
    s <- end[-.N]; e <- start[-1]
    .(start = s, end = e)
  }, by = chr]
  out <- out[end > start]
  out[, tad_id := sprintf("%s:%d-%d", chr, start, end)]
  out[]
}

# Insulation here is mean NPMI ACROSS the boundary at bin i (GAM_2), so a LOW
# value means few cross-boundary contacts, i.e. a strong boundary - boundaries
# are local MINIMA. Prominence is measured in local sd so one threshold works on
# chromosomes with different contact density.
gam_tad_insul <- list(smooth_bins = 5L, win = 10L, prominence = 0.5, min_bins = 10L)

call_tads_insulation <- function(genotype = c("WT", "HP1gKO"),
                                 chroms = gam_autosomes,
                                 pars = gam_tad_insul) {
  genotype <- match.arg(genotype)
  f <- file.path(results_dir, "rds", "npmi_bin_features.rds")
  if (!file.exists(f)) stop("run GAM_2 first: ", f, call. = FALSE)
  ft <- data.table::as.data.table(readRDS(f)$features)
  col <- paste0("insulation_", genotype)
  ft <- ft[chr %in% chroms, .(chr, start, end, v = get(col))]
  data.table::setorder(ft, chr, start)
  bd <- ft[, {
    y <- v
    g <- is.finite(y)
    num <- as.numeric(stats::filter(ifelse(g, y, 0), rep(1, pars$smooth_bins), sides = 2))
    den <- as.numeric(stats::filter(g * 1, rep(1, pars$smooth_bins), sides = 2))
    sm  <- ifelse(is.finite(den) & den > 0, num / den, NA_real_)
    n <- length(sm); idx <- seq_len(n)
    isMin <- vapply(idx, function(i) {
      w <- seq(max(1L, i - pars$win), min(n, i + pars$win))
      is.finite(sm[i]) && sum(is.finite(sm[w])) > pars$win &&
        sm[i] <= min(sm[w], na.rm = TRUE)
    }, logical(1))
    sdv <- stats::sd(sm, na.rm = TRUE)
    prom <- vapply(idx, function(i) {
      w <- seq(max(1L, i - pars$win), min(n, i + pars$win))
      if (!is.finite(sm[i]) || !is.finite(sdv) || sdv <= 0) return(NA_real_)
      (max(sm[w], na.rm = TRUE) - sm[i]) / sdv
    }, numeric(1))
    keep <- which(isMin & is.finite(prom) & prom >= pars$prominence)
    # enforce a minimum separation, keeping the deeper minimum of a close pair
    if (length(keep) > 1) {
      ord <- keep[order(sm[keep])]
      taken <- integer(0)
      for (i in ord) if (!any(abs(i - taken) < pars$min_bins)) taken <- c(taken, i)
      keep <- sort(taken)
    }
    .(bidx = keep, bstart = start[keep], bend = end[keep])
  }, by = chr]
  # boundaries -> intervals between them
  out <- bd[, { s <- bend[-.N]; e <- bstart[-1]; .(start = s, end = e) }, by = chr]
  out <- out[end > start]
  out[, tad_id := sprintf("%s:%d-%d", chr, start, end)]
  out[]
}

# --- Heterochromatin domains, merged -------------------------------------------
# het_domain_bed is 126,197 fragments of MEAN WIDTH 3.5 kb (444 Mb total). It is a
# merged PEAK set, not a domain set: GAM_10 reduce()s it and keeps whatever single
# fragment happens to exceed a width cut, which is not the same thing as a domain.
# Gap-merging first is what turns it into one. Parameters are swept in GAM_13 -
# these are the chosen values, and every domain-level result is reported against
# the unmerged set too, as the hotspot domains already are.
# CHOSEN FROM THE SWEEP IN GAM_13. The fragments percolate: genome coverage runs
# 1.3% (gap 0) -> 4.5% (10 kb) -> 17.5% (20 kb) -> 64.8% (50 kb) at min_width
# 100 kb, because once the gap exceeds the typical spacing between peaks
# everything fuses into chromosome-length blocks. 10 kb sits below that
# transition and yields 331 domains over 110 Mb, median 160 kb - genuinely large
# domains rather than the whole genome.
#
# Worth stating plainly: this BED contains very little large-domain structure.
# At gap 0 only 33 Mb sits in pieces >= 100 kb, so the 444 Mb it covers is almost
# entirely sub-50 kb peaks. GAM_10's "domains >= 100 kb" were whichever single
# fragments happened to be long, not merged domains.
gam_het_merge <- list(gap = 10000L, min_width = 100000L)

merge_het_domains <- function(gap = gam_het_merge$gap,
                              min_width = gam_het_merge$min_width,
                              chroms = gam_autosomes) {
  if (!file.exists(het_domain_bed)) stop("missing ", het_domain_bed, call. = FALSE)
  dm <- data.table::fread(het_domain_bed, header = FALSE, select = 1:3,
                          col.names = c("chr", "start", "end"))
  dm[, chr := ifelse(grepl("^chr", chr), chr, paste0("chr", chr))]
  dm <- dm[chr %in% chroms]
  g <- GenomicRanges::reduce(
         GenomicRanges::GRanges(dm$chr, IRanges::IRanges(dm$start + 1, dm$end)),
         min.gapwidth = gap, ignore.strand = TRUE)
  g <- g[GenomicRanges::width(g) >= min_width]
  g
}

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

# --- Hotspot DOMAINS: the second analysis level -----------------------------
# The hotspot counts are a dense, SMOOTH track, not a sparse hit list, so there
# are two legitimate units and they answer different questions:
#
#   per-bin count  - the continuous covariate. Highest resolution, no threshold
#                    chosen, right for correlation/coupling (GAM_3). Neighbouring
#                    bins are NOT independent (measured autocorrelation r = 0.59
#                    at 20 kb, 0.40 at 100 kb, 0.26 at 500 kb), so any CI over
#                    bins must resample chromosomes, not bins.
#   domain         - the independent unit. Right for enrichment/overlap tests and
#                    for "how much territory remodels", because counting adjacent
#                    bins as independent observations pseudo-replicates - the same
#                    reason chromHMM_size_weighted() permutes region labels rather
#                    than base pairs.
#
# WHY SMOOTH BEFORE THRESHOLDING. Thresholding the raw track fragments it: at the
# 90th percentile the median run is 1 bin and only 22% of hotspot bins fall in
# domains >= 100 kb, because a single bin dipping below the cutoff splits a real
# domain in two. Absorbing sub-threshold gaps fixes some of that (gap<=2 -> 62%)
# but chains through long quiet stretches if pushed (gap<=3 produces a 3.6 Mb
# "domain"), and the gap size is arbitrary. A 100 kb (5-bin) running mean first
# gives median 4 bins and 80% of hotspot bins in >= 100 kb domains without that
# pathology. 100 kb is chosen to match the autocorrelation structure: still
# r ~ 0.40, so averaging over it is justified, whereas a 500 kb window (r ~ 0.26)
# would fuse genuinely separate domains.
#
# Defaults are deliberately explicit so a sensitivity run is a parameter change:
#   q           quantile of the smoothed track that counts as "hot"
#   smooth_bins running-mean width in bins (5 = 100 kb)
#   gap         sub-threshold bins absorbed between hot runs (1 = 20 kb)
#   min_bins    minimum domain size kept (5 = 100 kb)
# --- Large contact-direction blocks -------------------------------------------
# A SECOND, much coarser scale than gam_hotspot_params. The hotspot domains are
# 0.1-2 Mb; eyeballing the tracks shows one-directional territory at 10-20 Mb -
# e.g. chr11 is almost entirely contact-LOSING over 0-20 Mb and almost entirely
# contact-GAINING over 56-108 Mb. A caller with a 100 kb floor and a q90 cutoff
# cannot see that: it fragments a 16 Mb block into dozens of small domains and
# reports them as unrelated.
#
# These blocks are called on the signed NET direction,
#     net = log2((lost + 1) / (gained + 1)),
# rather than on either density alone, because the feature is the IMBALANCE: a
# block can be busy in both directions and still be unremarkable. net > 0 means
# contacts are preferentially lost in the KO, net < 0 gained.
#
# Verified NOT to be a coverage artifact before building on it: the correlation
# between net and the WT-minus-KO count of valid NPMI cells is +0.005 at 4 Mb
# and +0.025 per bin, i.e. nothing. Blocks are not where data is missing.
gam_block_smooth_bins <- 101L    # ~2 Mb running mean; blocks are a Mb-scale feature
gam_block_net_cut     <- 0.5     # |net| above this counts as one-directional (~1.4x)
gam_block_min_mb      <- 4       # minimum block size to report

# --- Compartment-calling resolution -------------------------------------------
# A/B compartments are called at a COARSER resolution than everything else.
#
# WHY. The first real GAM_2 run called them at the native 20 kb and the axis was
# not usable: only 11 of 38 chromosome x genotype combinations reached |rho| > 0.4
# (25 reach 0.3, 3 reach 0.5; an oriented axis IS found on all 38, but it is PC1
# only on chr7 and PC2-PC4 elsewhere - see npmi_compartment_orientation.tsv)
# against the active-chromatin reference, 13 of 38 gave an implausible frac_A
# (<0.1 or >0.9, e.g. chr3 WT at 0.989), and the PC1-vs-PC2 choice split exactly
# 19/19 - a coin flip. Axis quality did not track chromosome size (rho = -0.04),
# so it was not a power problem, and the O-E step was verified correct (raw NPMI
# does decay with distance, 0.69 at 20 kb to 0.12 at 20 Mb, so it is not being
# double-corrected).
#
# The remaining explanation is resolution. A/B compartments are conventionally
# called at 100 kb - 1 Mb; at 20 kb the correlation matrix is dominated by
# TAD-scale local structure rather than the compartment signal, which is a
# long-range property. Aggregating to 100 kb (5 x 5 blocks of the native matrix)
# before the eigendecomposition asks the method for what it can actually resolve.
#
# This affects ONLY the compartment axis. Insulation, the distance bands and the
# decay slope stay at 20 kb, where they belong. The resulting PC and A/B call are
# mapped back onto the 20 kb grid so every downstream join still works on bin_id.
gam_compartment_binsize <- 100000L
stopifnot(gam_compartment_binsize %% gam_binsize == 0)
gam_compartment_factor  <- as.integer(gam_compartment_binsize / gam_binsize)   # 5

gam_hotspot_params <- list(q = 0.90, smooth_bins = 5L, gap = 1L, min_bins = 5L)

# UNADJUSTED variant, carried alongside the smoothed one so that every
# domain-level result can be reported both ways. This is threshold-and-nothing-
# else: no smoothing, strict adjacency, no size floor. It is the conservative
# reading - it makes no assumption that a sub-threshold bin inside a run belongs
# to the domain - at the cost of fragmenting real domains into singletons (at
# q90, median run length 1 bin, 78% of hotspot bins in runs under 100 kb).
#
# Keeping both is the same instinct as reporting enrichments per-region AND
# size-weighted on the ChIP side: where the two agree the conclusion does not
# rest on the smoothing choice, and where they disagree that disagreement is
# itself the finding and belongs in the write-up.
#
# Note the truly unadjusted view is the per-bin count itself, which the pipeline
# always retains as n_contacts_lost / n_contacts_gained and uses as the primary
# continuous covariate. This variant is about DOMAIN construction only.
gam_hotspot_params_raw <- list(q = 0.90, smooth_bins = 1L, gap = 0L, min_bins = 1L)

# Returns a list: $bin (per-bin domain_id + smoothed value + hot flag, in input
# order) and $domains (one row per kept domain). Input must be sorted by chr,start.
gam_hotspot_domains <- function(counts, chr, start, end,
                                q           = gam_hotspot_params$q,
                                smooth_bins = gam_hotspot_params$smooth_bins,
                                gap         = gam_hotspot_params$gap,
                                min_bins    = gam_hotspot_params$min_bins) {
  dt <- data.table::data.table(i = seq_along(counts), chr = as.character(chr),
                               start = start, end = end, n = as.numeric(counts))
  # Running mean WITHIN each chromosome - smoothing across a chromosome boundary
  # would invent signal where there is a discontinuity.
  dt[, sm := {
    if (smooth_bins <= 1L) n else {
      v <- as.numeric(stats::filter(n, rep(1 / smooth_bins, smooth_bins), sides = 2))
      ifelse(is.na(v), n, v)      # edges keep their raw value
    }
  }, by = chr]
  cut_at <- stats::quantile(dt$sm, q, na.rm = TRUE)
  dt[, hot := sm >= cut_at]

  # Absorb short sub-threshold gaps, then label contiguous runs per chromosome.
  dt[, dom := {
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
    rr <- rle(h)
    id <- rep(NA_integer_, length(h))
    ends <- cumsum(rr$lengths); begs <- ends - rr$lengths + 1L
    w <- which(rr$values)
    if (length(w)) for (j in seq_along(w)) id[begs[w[j]]:ends[w[j]]] <- j
    id
  }, by = chr]
  dt[!is.na(dom), dom_key := paste0(chr, ":", dom)]

  # NOTE the name: this must NOT be called `dom`, because dt already has a `dom`
  # COLUMN and data.table evaluates the filter below inside dt's scope - the
  # column would shadow this table and `dom$dom_key` would hit an integer vector.
  domtab <- dt[!is.na(dom), .(chr = chr[1], start = min(start), end = max(end),
                              n_bins = .N, mean_count = mean(n), max_count = max(n),
                              sum_count = sum(n)), by = dom_key]
  domtab <- domtab[n_bins >= min_bins]
  domtab[, width_kb := (end - start) / 1000]
  # Drop membership for domains that failed the size filter, so the per-bin
  # domain_id and the $domains table can never disagree.
  keep_keys <- domtab$dom_key
  dt[!dom_key %in% keep_keys, dom_key := NA_character_]
  data.table::setorder(dt, i)
  list(bin = dt[, .(domain_id = dom_key, smoothed = sm, hot)],
       domains = domtab[order(chr, start)],
       params = list(q = q, cutoff = cut_at, smooth_bins = smooth_bins,
                     gap = gap, min_bins = min_bins))
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
