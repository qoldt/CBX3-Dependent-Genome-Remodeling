# ================================================================
# prepare_CTCF_annotation.R - CTCF peaks for mm39
# ----------------------------------------------------------------
#     cd GAM && Rscript prepare_CTCF_annotation.R
#
# WHY. The repo has two independent TAD partitions - published boundaries
# (load_tads_external) and insulation minima called from our own NPMI
# (call_tads_insulation) - and no third, orthogonal line of evidence for where a
# boundary is. CTCF is that evidence: it is the architectural protein that marks
# boundaries, it is measured by a different assay in a different lab, and it
# cannot move with HP1g loss in our data because it is external to it. It also
# gives the hotspot-domain composition test (GAM_4 section C) an annotation that
# is a point feature rather than a broad domain.
#
# SOURCE. ENCODE, conservative IDR thresholded peaks, mm10:
#   ENCSR877MSN  ENCFF093ZRU  layer of hippocampus, B6CASTF1/J male adult 2 mo
#   ENCSR644VYX  ENCFF395PXZ  left cerebral cortex,  B6CASTF1/J female adult 2 mo
#
# CELL TYPE - READ THIS BEFORE QUOTING ANYTHING FROM THIS TRACK. Both are BULK
# TISSUE, not sorted CA1 pyramidal neurons. No CA1-pyramidal-specific CTCF
# ChIP-seq is public as far as could be found. CTCF occupancy is largely
# cell-type-invariant, which is what makes bulk a usable proxy, but this track
# must be labelled "bulk hippocampus"/"bulk cortex" wherever it appears and must
# never be described as CA1. Both tissues are taken because they match the two
# assays between them: GAM is CA1 pyramidal cells (hippocampus) and the ChIP is
# Emx1-Cre cortex + hippocampus.
#
# ASSEMBLY. ENCODE ships mm10; everything here is mm39. Lifted with rtracklayer
# and UCSC's mm10ToMm39 chain - NOT CrossMap, which is not installed on this
# machine (the KAP1 BEDs in the same directory were done with CrossMap elsewhere,
# per ChIP/README.md). The chain is already cached at
# ChIP/data/annotation/mm10ToMm39.over.chain.gz; it is downloaded if absent.
#
# NOT reduce()d, and NO width floor. prepare_LAD_annotation.R collapses its
# lifted ranges with reduce() and drops anything under 10 kb as liftOver
# shrapnel, which is right for LADs (hundreds of kb) and would destroy this set:
# CTCF peaks are ~200-500 bp point features, so reduce() would fuse neighbouring
# sites into one and a 10 kb floor would delete every peak in the file. Peaks are
# kept individual here; only zero-width artefacts are dropped.
#
# Seqnames stay UCSC (chr1), matching the chain, the GAM bins and the other
# annotation BEDs. Output goes to annot_dir, which is git-ignored - the script is
# the committed artifact, not the BED.
# ================================================================
source("config.R")
suppressPackageStartupMessages({ library(rtracklayer); library(GenomeInfoDb) })

work <- file.path(tempdir(), "ctcf"); dir.create(work, showWarnings = FALSE)

# Prefer the cached chain next to the other annotation; fall back to UCSC.
chain_f <- file.path(work, "mm10ToMm39.over.chain")
cached  <- file.path(annot_dir, "mm10ToMm39.over.chain.gz")
if (!file.exists(chain_f)) {
  src_gz <- file.path(work, "mm10ToMm39.over.chain.gz")
  if (file.exists(cached)) {
    file.copy(cached, src_gz, overwrite = TRUE)
    message("chain: using the cached copy in ", annot_dir)
  } else {
    utils::download.file(
      "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz",
      src_gz, quiet = TRUE)
    message("chain: downloaded from UCSC")
  }
  R.utils::gunzip(src_gz, destname = chain_f, overwrite = TRUE, remove = FALSE)
}
ch <- import.chain(chain_f)

enc <- "https://www.encodeproject.org/files/%s/@@download/%s.bed.gz"
sets <- list(hippocampus = c(file = "ENCFF093ZRU", expt = "ENCSR877MSN"),
             cortex      = c(file = "ENCFF395PXZ", expt = "ENCSR644VYX"))

for (nm in names(sets)) {
  acc <- sets[[nm]][["file"]]
  src <- file.path(work, paste0(acc, ".bed.gz"))
  if (!file.exists(src)) utils::download.file(sprintf(enc, acc, acc), src, quiet = TRUE)

  # ENCODE narrowPeak: 10 columns. import(format="narrowPeak") keeps signalValue
  # and qValue, which are worth carrying so a peak set can be thresholded later.
  g10 <- tryCatch(import(src, format = "narrowPeak"),
                  error = function(e) import(src, format = "bed"))
  g10 <- GenomeInfoDb::keepStandardChromosomes(g10, pruning.mode = "coarse")
  n_in <- length(g10)

  # unlist() without reduce(): one lifted range per source peak wherever the lift
  # is clean. A peak that fragments contributes its pieces; sort for tidiness.
  lo <- sort(unlist(liftOver(g10, ch)))
  lo <- lo[width(lo) > 0]
  n_out <- length(lo)

  out <- file.path(annot_dir, sprintf("CTCF_%s_mm39.bed", nm))
  export(lo, out, format = "bed")
  message(sprintf("%-12s %s/%s  mm10 %6d peaks -> mm39 %6d (%.1f%% kept), median width %3.0f bp  -> %s",
                  nm, sets[[nm]][["expt"]], acc, n_in, n_out,
                  100 * n_out / n_in, median(width(lo)), basename(out)))
  if (n_out < 0.9 * n_in)
    warning(nm, ": only ", round(100 * n_out / n_in), "% of peaks survived the lift",
            " - expected >95% for mm10->mm39", call. = FALSE)
}

message("\nBulk tissue, NOT sorted CA1 - see the header before quoting this track.",
        "\nLoad with load_ctcf(); overlay on the 20 kb grid with to_ucsc() + findOverlaps.")
