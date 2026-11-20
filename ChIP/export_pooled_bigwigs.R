# ================================================================
# export_pooled_bigwigs.R - one pooled bigWig per mark x genotype, for IGV
# ----------------------------------------------------------------
#   cd ChIP && Rscript export_pooled_bigwigs.R
#   POOL_MARKS=H3K9me2 Rscript export_pooled_bigwigs.R     (one mark only)
#
# Utility, not a pipeline stage: nothing downstream reads these files. They exist
# so the GAM IGV session (GAM/export_igv.R) can show ChIP signal next to the
# peak BEDs.
#
# WHAT "POOLED" MEANS HERE: the MEAN of the scaled bigWigs of every
# sample the analysis uses, i.e. the sample sheet AFTER the exclusion in
# config.R (HP1gKO rep1 has no usable library and rep2 is excluded at QC). The
# MINUTE bigWigs are already scaled to be comparable across samples, so a mean
# keeps that scale; a sum would not. WT pools six libraries and HP1gKO four, so
# the mean (not the sum) is also what keeps the two genotypes comparable.
#
# RESOLUTION: 50 bp bins (the mean over each bin). At single-base resolution the
# mean of 4-6 libraries fragments into so many runs that one pooled file was
# 2.2 GB against 0.2-0.3 GB per input; 50 bp is far below anything readable
# for a histone mark in a browser and keeps each file small.
#
# Output: results/igv/pooled/<mark>_<genotype>_pooled.bw, UCSC seqnames
# (chr1...), autosomes + chrX + chrY. results/ is git-ignored; this script is the
# committed artifact.
# ================================================================
source("config.R")
suppressPackageStartupMessages(library(rtracklayer))

out_dir <- file.path(results_dir, "igv", "pooled")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

marks <- unique(samples$mark)
e <- Sys.getenv("POOL_MARKS", "")
if (nzchar(e)) marks <- intersect(marks, trimws(strsplit(e, ",")[[1]]))
stopifnot(length(marks) > 0)

chroms <- c(as.character(1:19), "X", "Y")    # NCBI, as stored in the bigWigs
bin_bp <- 50L

for (mk in marks) for (g in c("WT", "HP1gKO")) {
  s <- samples[samples$mark == mk & samples$genotype == g, ]
  files <- s$bigwig_path
  stopifnot(length(files) > 0, all(file.exists(files)))
  out <- file.path(out_dir, sprintf("%s_%s_pooled.bw", mk, g))
  message(sprintf("%s %s: mean of %d bigWigs (%s)", mk, g, length(files),
                  paste(s$sample_id, collapse = ", ")))
  si  <- seqinfo(BigWigFile(files[1]))
  chr <- intersect(chroms, seqnames(si))
  # One chromosome at a time, so only one chromosome's coverage per library is
  # in memory at once.
  cov <- lapply(setNames(chr, chr), function(ch) {
    rng <- GRanges(ch, IRanges(1, seqlengths(si)[[ch]]))
    acc <- NULL
    for (f in files) {
      x <- import(f, which = rng, as = "RleList")[[ch]]
      acc <- if (is.null(acc)) x else acc + x
    }
    acc <- acc / length(files)
    # mean over fixed bins, then merge adjacent equal bins back into runs
    n  <- length(acc)
    st <- seq.int(1L, n, by = bin_bp)
    m  <- viewMeans(Views(acc, start = st, end = pmin(st + bin_bp - 1L, n)))
    r  <- Rle(round(m, 4))
    en <- pmin(st[cumsum(runLength(r))] + bin_bp - 1L, n)
    ss <- st[c(1L, head(cumsum(runLength(r)), -1L) + 1L)]
    GRanges(paste0("chr", ch), IRanges(ss, en), score = runValue(r))
  })
  lens <- setNames(seqlengths(si)[chr], paste0("chr", chr))
  gr <- unlist(GRangesList(unname(cov)))
  gr <- gr[score(gr) != 0]
  seqlevels(gr) <- names(lens); seqlengths(gr) <- lens
  export(gr, BigWigFile(out))
  stopifnot(file.exists(out), file.size(out) > 1e4)
  message(sprintf("  wrote %s (%.0f MB)", out, file.size(out) / 1e6))
}
