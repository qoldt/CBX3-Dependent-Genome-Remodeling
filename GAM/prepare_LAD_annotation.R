# ================================================================
# Prepare lamina-associated domain (LAD) annotation in mm39
# ----------------------------------------------------------------
#     cd GAM && Rscript prepare_LAD_annotation.R
#
# WHY. Our own A/B compartment axis is weakly determined (only 11 of 38
# chromosome x genotype axes reach |rho| > 0.4), and it is selected against our
# own ChIP marks, which makes any mark-vs-compartment statement partly circular.
# A published LAD annotation is external to both genotypes and to all five
# marks, so it is a cleaner axis for "is this peripheral heterochromatin".
#
# SOURCE. Peric-Hupkes et al. 2010, GEO GSE17051 - LaminB1 DamID across mouse
# ESC, neural precursor, astrocyte and MEF.
#   cLAD_regions  - LAD in ALL FOUR cell types (constitutive)
#   ciLAD_regions - LAD in NONE of them (constitutive inter-LAD)
#
# CELL TYPE, STATED PLAINLY: these two sets are CONSTITUTIVE, i.e. the
# cell-type-invariant core. They are NOT a neuronal annotation. The same GEO
# entry carries per-cell-type calls (including neural precursor and astrocyte)
# in GSE17051_HMM_state_calls_per_probe.txt.gz, and GSE36132 adds Lamin A DamID
# in NPC and astrocyte. Neither is a mature CA1 pyramidal neuron, which is what
# the GAM is. Use the constitutive sets for robustness; switch to the
# brain-lineage calls when cell-type specificity matters, and say which was used.
#
# ASSEMBLY. The deposited BEDs are mm9 (chr1 runs to 197.2 Mb, the mm9 length).
# Lifted mm9 -> mm39 directly with rtracklayer and UCSC's mm9ToMm39 chain;
# fragments under 10 kb after lifting are dropped as liftOver shrapnel. 99% of
# base pairs survive in both sets, so this is not a lossy conversion.
# ================================================================
source("config.R")
suppressPackageStartupMessages({ library(rtracklayer); library(GenomeInfoDb) })

work <- file.path(tempdir(), "lad"); dir.create(work, showWarnings = FALSE)
geo  <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE17nnn/GSE17051/suppl/GSE17051_%s_regions.bed.gz"
chain_url <- "https://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm39.over.chain.gz"

chain_gz <- file.path(work, "mm9ToMm39.over.chain.gz")
chain_f  <- sub("[.]gz$", "", chain_gz)
if (!file.exists(chain_f)) {
  utils::download.file(chain_url, chain_gz, quiet = TRUE)
  R.utils::gunzip(chain_gz, destname = chain_f, overwrite = TRUE, remove = FALSE)
}
ch <- import.chain(chain_f)

for (nm in c("cLAD", "ciLAD")) {
  src <- file.path(work, paste0(nm, ".bed.gz"))
  if (!file.exists(src)) utils::download.file(sprintf(geo, nm), src, quiet = TRUE)
  g9 <- GenomeInfoDb::keepStandardChromosomes(import(src, format = "bed"),
                                              pruning.mode = "coarse")
  lo <- reduce(sort(unlist(liftOver(g9, ch))))
  lo <- lo[width(lo) >= 10000]
  out <- file.path(annot_dir, paste0(nm, "_mm39.bed"))
  export(lo, out, format = "bed")
  message(sprintf("%-6s mm9 %5d regions %6.1f Mb -> mm39 %5d regions %6.1f Mb (%.0f%% bp kept)  -> %s",
          nm, length(g9), sum(as.numeric(width(g9))) / 1e6, length(lo),
          sum(as.numeric(width(lo))) / 1e6,
          100 * sum(as.numeric(width(lo))) / sum(as.numeric(width(g9))), out))
}
message("\nOverlay onto the 20 kb grid with to_ucsc() + findOverlaps; see GAM_13 for the tests.")
