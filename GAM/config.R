# ================================================================
# GAM pipeline - configuration: input files, output tree, parameters
# ----------------------------------------------------------------
# Sourced by every GAM_* stage and by run_GAM.R. SINGLE SOURCE OF TRUTH for
# values - edit them here only; stage scripts must not redefine them.
#
#   config.R        values only (this file)
#   utils.R         libraries, path builders, loaders, helper functions
#                   (sourced at the bottom, so stages need only config.R)
#   PARAMETERS.md   why each non-obvious value is set the way it is
#
# Run FROM the GAM/ directory:   cd GAM && Rscript run_GAM.R
# ================================================================

# ---- Input data (read-only, NOT committed - ~7.7 GB, see README) ----
# GAM_DATA defaults to this folder; set it to use a copy elsewhere.
gam_data    <- path.expand(Sys.getenv("GAM_DATA", unset = "."))
npmi_dir    <- file.path(gam_data, "NPMI_matrices")   # {WT,cKO}_at20Kb/ dense NPMI
diff_dir    <- file.path(gam_data, "Diff_contacts")   # .bg2 differential contact pairs
hotspot_dir <- file.path(gam_data, "Hotspots")        # per-bin differential counts

# ---- Output tree (git-ignored, regenerable) -----------------------
# GAM_RESULTS=results_alt writes a sensitivity run to a separate tree.
results_dir <- path.expand(Sys.getenv("GAM_RESULTS", unset = "results"))
counts_dir  <- file.path(results_dir, "counts")
rds_dir     <- file.path(results_dir, "rds")
tables_dir  <- file.path(results_dir, "tables")
fig_dir     <- file.path(results_dir, "figures")   # save_fig() reads this by name
bed_dir     <- file.path(results_dir, "bed")

# ---- GAM file parameters: these MATCH THE FILENAMES on disk -------
# Changing one changes which files the path builders resolve to; it recomputes
# nothing. The significance cutoffs are baked in upstream, not chosen here.
gam_binsize         <- 20000L        # segregation_at20000
gam_sigma           <- "2.6"         # sigma.2.6
gam_seg_date        <- "2026-02-23"
gam_seg_prefix      <- "HP1g_mm39"
gam_z               <- 2             # z2      - differential contact z-score cutoff
gam_npmi            <- 0.1           # npmi0.1 - minimum NPMI in the selected genotype
gam_maxdist_label   <- "inf"         # maxdistinf
gam_hotspot_maxdist <- 5e6           # dist0-5000000 - window for the hotspot counts

# ---- Genotype labels: the one place "cKO" appears ------------------
# Canonical repo-wide is WT / HP1gKO; the GAM files say cKO. Genotype is encoded
# ONLY by the directory name - the NPMI file names are identical in both.
gam_genotypes <- c(WT = "WT", HP1gKO = "cKO")

# ---- Direction of the differential calls: READ BEFORE USING -------
# selWT  = contact stronger in WT  = LOST   in HP1gKO
# selcKO = contact stronger in cKO = GAINED in HP1gKO
# Getting this backwards inverts the headline result.
gam_sel <- c(lost = "selWT", gained = "selcKO")

# ---- Chromosomes --------------------------------------------------
# Diff_contacts and Hotspots are autosome-only; NPMI also has chrX and chrY.
gam_autosomes   <- paste0("chr", 1:19)
gam_npmi_chroms <- c(gam_autosomes, "chrX", "chrY")

# ---- Inputs from the ChIP assay -----------------------------------
# Not via source("../ChIP/config.R"): it defines results_dir / fig_dir / ... by
# the same names and would redirect every GAM output into ChIP/results. These
# MUST track ChIP/config.R (same MINUTE_DATA root, sample sheet, exclusion).
chip_results_dir <- path.expand(Sys.getenv("CHIP_RESULTS", unset = "../ChIP/results"))
chip_rds         <- file.path(chip_results_dir, "rds",
                              "annotated_results_H3K9me3_H4K20me3_2000bp_merged.rds")
minute_data      <- path.expand(Sys.getenv("MINUTE_DATA", unset = "../ChIP/data"))
bigwig_dir       <- file.path(minute_data, "bigwig")
annot_dir        <- file.path(minute_data, "annotation")
chip_samples_tsv <- "../ChIP/samples.tsv"
# Samples excluded by default: HP1gKO rep2, as on the ChIP side
# (MINUTE_EXCLUDE=none keeps all; read by gam_load_samples()).
chip_peaks_dir   <- "../ChIP/data/peaks"
het_domain_bed   <- file.path(chip_peaks_dir, "2000bp_merge_H3K9me3_H4K20me3.bed")

# ---- External annotation (in annot_dir) ---------------------------
# Provenance and why these sets: PARAMETERS.md "Annotation".
lad_bed          <- file.path(annot_dir, "cLAD_mm39.bed")    # constitutive LAD, Peric-Hupkes 2010
cilad_bed        <- file.path(annot_dir, "ciLAD_mm39.bed")   # constitutive inter-LAD
ctcf_hip_bed     <- file.path(annot_dir, "CTCF_hippocampus_mm39.bed")  # ENCODE, BULK tissue
ctcf_ctx_bed     <- file.path(annot_dir, "CTCF_cortex_mm39.bed")       # ENCODE, BULK tissue
tad_boundary_bed <- file.path(annot_dir, "TAD_boundaries_mm39.bed")    # published boundaries

# ---- ChIP quantification on the GAM grid (GAM_1) -------------------
gam_count_scale <- 10      # x coverage before rounding for DESeq2; cancels in log2FC
gam_cores <- as.integer(Sys.getenv("GAM_CORES", unset = max(1, parallel::detectCores() - 1)))

# ChIP mark-loss groups (stable / H4K20me3_only / co_loss): MINUTE_3's cutoffs
# VERBATIM - do not retune here. Absolute, so not shift-invariant; the
# median-centred variant is primary downstream.
H4_LOST_CUT <- -0.5        # H4K20me3 log2FC below this = "H4K20me3-lost"
K9_LOSS_CUT <- -0.3        # H3K9me3  log2FC below this = "H3K9me3-loss"
UNCH_CUT    <-  0.15       # |log2FC| below this = "unchanged"

# ---- Own TAD calls from insulation minima (call_tads_insulation) ---
gam_tad_insul <- list(smooth_bins = 5L, win = 10L, prominence = 0.5, min_bins = 10L)

# ---- Heterochromatin domains (merge_het_domains) ------------------
# Chosen from the GAM_13 sweep: 10 kb sits below the gap at which fragments
# fuse into chromosome-length blocks. PARAMETERS.md "Heterochromatin domains".
gam_het_merge <- list(gap = 10000L, min_width = 100000L)

# ---- Chromatin classes (gam_chromatin_classes) --------------------
# LEVEL 1 euchromatin | heterochromatin; LEVEL 2 nests inside it. Three states
# are assigned against their ChromHMM name because our own WT ChIP contradicts
# it (Quies2/Quies3 -> H3K9me3, QuiesG -> Tx): PARAMETERS.md "Chromatin classes".
gam_chromatin_state_map <- c(
  # --- euchromatin ---
  Tss = "TSS", TssFlnk = "TSS",
  Enh = "Enhancer", EnhG = "Enhancer", EnhLo = "Enhancer",
  EnhPr = "Enhancer", EnhPois = "Enhancer",
  Tx = "Tx", TxWk = "Tx", QuiesG = "Tx",
  # --- heterochromatin ---
  TssBiv = "Polycomb", ReprPC = "Polycomb", ReprPCWk = "Polycomb",
  Het = "H3K9me3", Quies2 = "H3K9me3", Quies3 = "H3K9me3",
  Quies = "Quiescent", Quies4 = "Quiescent")
gam_chromatin_level1 <- c(TSS = "euchromatin", Enhancer = "euchromatin", Tx = "euchromatin",
                          Polycomb = "heterochromatin", H3K9me3 = "heterochromatin",
                          Quiescent = "heterochromatin")
gam_class2_levels <- c("TSS", "Enhancer", "Tx", "Polycomb", "H3K9me3", "Quiescent")
gam_class1_levels <- c("euchromatin", "heterochromatin")
gam_het_peak_cut  <- 0.5   # bin peak coverage >= this -> H3K9me3, whatever ChromHMM says

# ---- Distance bands for every contact-by-separation analysis ------
# The >50Mb band is open-ended, so it is not comparable across chromosomes.
gam_band_edges <- c(0, 2e5, 1e6, 5e6, 1e7, 5e7, Inf)
gam_band_names <- c("<0.2Mb", "0.2-1Mb", "1-5Mb", "5-10Mb", "10-50Mb", ">50Mb")

# ---- Hotspot domains (gam_hotspot_domains) ------------------------
#   q quantile counted as hot | smooth_bins running mean (5 = 100 kb)
#   gap sub-threshold bins absorbed | min_bins minimum domain size
gam_hotspot_params     <- list(q = 0.90, smooth_bins = 5L, gap = 1L, min_bins = 5L)
gam_hotspot_params_raw <- list(q = 0.90, smooth_bins = 1L, gap = 0L, min_bins = 1L)  # unsmoothed variant

# ---- Large contact-direction blocks (GAM_4) -----------------------
# Called on net = log2((gained + 1) / (lost + 1)); net > 0 = gained in the cKO.
gam_block_smooth_bins <- 101L    # ~2 Mb running mean
gam_block_net_cut     <- 0.5     # |net| above this counts as one-directional (~1.4x)
gam_block_min_mb      <- 4       # minimum block size to report

# ---- Code ----------------------------------------------------------
source("utils.R")
