# ================================================================
# Minimal IGV export - ONE track per data type
# ----------------------------------------------------------------
#   cd GAM && Rscript export_igv.R
#
# Writes results/igv/ : the GAM tracks, pooled ChIP signal and peak sets, plus a session file that loads them in a
# sensible order. Deliberately minimal - one representative of each KIND of
# data, not every mark, genotype, direction and variant. The point is a view
# you can open and read, not a complete dump; the rest is in the .rds and the
# tables.
#
#   1 chip_H4K20me3_log2FC.bedgraph  the ChIP effect   (the mark that moves)
#   2 gam_contacts_lost.bedgraph     the GAM effect    (differential contacts)
#   3 gam_hotspot_domains.bed        called domains    (the adjusted variant)
#   4 gam_melting_score.bedgraph     MELTRON           (from GAM_2c)
#   5 chip_compartments.bed          the ChIP mark-loss groups (NOT A/B)
#
# Everything is mm39 with UCSC seqnames, which is what IGV expects. bedGraph
# coordinates are 0-based half-open, matching $annot$start/$end as stored, so no
# conversion is applied here - see gam_bin_id() in shared/common.R for why that
# convention is load-bearing.
# ================================================================
source("config.R")

igv_dir <- file.path(results_dir, "igv")
dir.create(igv_dir, recursive = TRUE, showWarnings = FALSE)

g1 <- file.path(rds_dir, "gam_bins_annotated.rds")
if (!file.exists(g1))
  stop("GAM_1 output not found:\n  ", g1,
       "\nRun it first:  cd GAM && Rscript GAM_1_Bins_and_Signal.R", call. = FALSE)
x <- readRDS(g1)
a <- as.data.table(x$annot); setorder(a, chr, start)

wr <- function(dt, file, track_line) {
  p <- file.path(igv_dir, file)
  writeLines(track_line, p)
  fwrite(dt, p, sep = "\t", col.names = FALSE, append = TRUE, quote = FALSE,
         scipen = 999)
  message("  ", file, "  (", nrow(dt), " rows)")
  invisible(p)
}

message("Writing IGV tracks to ", igv_dir)

# --- 1. ChIP effect: per-bin H4K20me3 log2FC --------------------------------
# One mark, the one carrying the result. Diverging around 0; the viewRange is
# set from the data's own 1st-99th percentile so the track is readable without
# fiddling, rather than squashed by a handful of outlier bins.
l2 <- x$l2fc[, "H4K20me3"][match(a$bin_id, rownames(x$l2fc))]
d1 <- data.table(a$chr, a$start, a$end, round(l2, 4))[is.finite(V4)]
rng <- quantile(d1$V4, c(0.01, 0.99), na.rm = TRUE)
wr(d1, "chip_H4K20me3_log2FC.bedgraph", sprintf(
  paste0('track type=bedGraph name="H4K20me3 log2FC (KO/WT)" ',
         'description="per 20kb bin, DESeq2 sizeFactors=1" visibility=full ',
         'color=178,24,43 altColor=33,102,172 viewLimits=%.2f:%.2f windowingFunction=mean'),
  rng[1], rng[2]))

# --- 2. GAM effect: differential contacts lost ------------------------------
# selWT = stronger in WT = LOST in the KO. One direction only; "gained" is the
# mirror and lives in the tables.
d2 <- data.table(a$chr, a$start, a$end, a$n_contacts_lost)
wr(d2, "gam_contacts_lost.bedgraph", sprintf(
  paste0('track type=bedGraph name="GAM contacts lost (selWT)" ',
         'description="significant differential contacts within 0-5Mb, z>2 npmi>0.1" ',
         'visibility=full color=70,70,70 viewLimits=0:%d windowingFunction=mean'),
  as.integer(quantile(d2$V4, 0.99))))

# --- 3. Called hotspot domains ----------------------------------------------
# The adjusted variant (smoothed, q90, gap<=1, >=100kb). The unadjusted variant
# is deliberately NOT exported: the two agree on only ~40% of bins, and showing
# both here invites reading one as confirming the other.
dom_f <- file.path(tables_dir, "hotspot_domains_lost_adj.tsv")
if (file.exists(dom_f)) {
  dm <- fread(dom_f)
  d3 <- data.table(dm$chr, dm$start, dm$end,
                   sprintf("dom%04d_n%d", seq_len(nrow(dm)), dm$n_bins),
                   pmin(1000, round(dm$mean_count * 10)), ".")
  wr(d3, "gam_hotspot_domains.bed",
     paste0('track name="Hotspot domains (lost, adjusted)" ',
            'description="score = mean contacts x10; see QC_Hotspot_Domains.R" ',
            'visibility=dense color=200,90,40 useScore=1'))
} else message("  (no hotspot domain table - run GAM_1)")

# --- 4. MELTRON melting track (written by GAM_2c) -----------------------------
# No A/B compartment track: compartment calling from NPMI was removed on
# 2026-11-20 (CLAUDE.md "A/B compartments are not called").
ms <- file.path(igv_dir, "gam_melting_score.bedgraph")
message(if (file.exists(ms)) "  gam_melting_score.bedgraph  (from GAM_2c)"
        else "  (no melting track - run GAM_2c_Meltron.R)")

# --- 5. ChIP differential-loss compartments ---------------------------------
# The median-centred labels (the primary ones). itemRgb so the three groups are
# distinguishable at a glance; colours match the project palette.
hex2rgb <- function(h) paste(as.vector(col2rgb(h)), collapse = ",")
cols <- c(stable = hex2rgb("grey65"),
          H4K20me3_only = hex2rgb(gaby_cols[3]),
          co_loss = hex2rgb(gaby_cols[5]))
cp <- a[!is.na(compartment_centred),
        .(chr, start, end, name = as.character(compartment_centred))]
if (nrow(cp)) {
  d5 <- data.table(cp$chr, cp$start, cp$end, cp$name, 0L, ".",
                   cp$start, cp$end, unname(cols[cp$name]))
  wr(d5, "chip_compartments.bed",
     paste0('track name="ChIP compartments (median-centred)" ',
            'description="grey=stable, blue=H4K20me3_only, plum=co_loss" ',
            'visibility=dense itemRgb=On'))
}

# --- 6. Large one-directional contact blocks (GAM_6) -------------------------
# The dominant large-scale feature: 10-25 Mb territory that is consistently
# losing or gaining contacts. Copied from GAM_6's bed so IGV and the tables
# cannot drift apart. Red = GAINED in the cKO, blue = lost - the repo convention.
blk_src <- file.path(bed_dir, "gam6_large_blocks.bed")
if (file.exists(blk_src)) {
  file.copy(blk_src, file.path(igv_dir, "gam_large_blocks.bed"), overwrite = TRUE)
  message("  gam_large_blocks.bed  (copied from GAM_6)")
} else message("  (no large-block bed - run GAM_4_Large_Scale.R)")

# --- 7. ChIP: pooled signal per mark x genotype, and the peak sets ------------
# The pooled bigWigs are written by ../ChIP/export_pooled_bigwigs.R (mean of the
# analysed libraries, 50 bp bins) and are referenced in place rather than
# copied: they are ~0.15-0.4 GB each. The session uses a relative path, so the
# repo layout must stay as it is for IGV to find them.
#
# Peak sets are copied with UCSC seqnames, because the ChIP BEDs are NCBI ("1")
# and the rest of the session is UCSC. Two kinds, named so they cannot be
# confused: `regions_*` are the region sets the ChIP analysis QUANTIFIES (from
# ChIP/config.R $regions - H3K9me2, H3K9me3 and H4K20me3 all share the 2 kb
# merge since 2026-11-20); `consensus_*` are MergePeaks_to_masterPeak.R's
# >=2-replicate MACS3 sets, shown for reference and not used by the analysis.
chip_marks <- c("H3K4me3", "H3K9me2", "H3K9me3", "H4K20me3", "H3K36me3")
pooled_rel <- file.path("..", "..", "..", "ChIP", "results", "igv", "pooled")
pooled <- unlist(lapply(chip_marks, function(mk)
  file.path(pooled_rel, sprintf("%s_%s_pooled.bw", mk, c("WT", "HP1gKO")))))
names(pooled) <- unlist(lapply(chip_marks, function(mk) paste(mk, c("WT", "HP1gKO"))))
pooled <- pooled[file.exists(file.path(igv_dir, pooled))]
message(sprintf("  %d pooled ChIP bigWigs referenced%s", length(pooled),
                if (!length(pooled)) " (none - run ../ChIP/export_pooled_bigwigs.R)" else ""))

peak_dir <- file.path(igv_dir, "peaks")
dir.create(peak_dir, showWarnings = FALSE)
peak_sets <- c(
  regions_H3K9me2_H3K9me3_H4K20me3 = file.path(chip_peaks_dir, "2000bp_merge_H3K9me3_H4K20me3.bed"),
  regions_H3K4me3  = file.path(chip_peaks_dir, "master_peaks_consensus/H3K4me3_consensus_masterPeak.bed"),
  regions_H3K36me3 = file.path(chip_peaks_dir, "merged_H3K36me3_masterPeak.bed"),
  consensus_H3K9me2  = file.path(chip_peaks_dir, "master_peaks_consensus/H3K9me2_consensus_masterPeak.bed"),
  consensus_H3K9me3  = file.path(chip_peaks_dir, "master_peaks_consensus/H3K9me3_consensus_masterPeak.bed"),
  consensus_H4K20me3 = file.path(chip_peaks_dir, "master_peaks_consensus/H4K20me3_consensus_masterPeak.bed"),
  consensus_H3K36me3 = file.path(chip_peaks_dir, "master_peaks_consensus/H3K36me3_consensus_masterPeak.bed"),
  # the H3K9me2 set used until 2026-11-20 - no committed recipe; kept so the
  # change can be looked at
  retired_H3K9me2_merged = file.path(chip_peaks_dir, "merged_H3K9me2_masterPeak.bed"))
peak_files <- character()
for (nm in names(peak_sets)) {
  f <- peak_sets[[nm]]
  if (!file.exists(f)) { message("  (missing peak set: ", f, ")"); next }
  pk <- fread(f, header = FALSE, select = 1:3)
  pk[, V1 := ifelse(grepl("^chr", V1), V1, paste0("chr", V1))]
  pk[V1 == "chrMT", V1 := "chrM"]
  out <- file.path("peaks", paste0(nm, ".bed"))
  wr(pk, out, sprintf('track name="%s" description="%s" visibility=dense', nm,
                      if (startsWith(nm, "regions_")) "quantified by the ChIP analysis"
                      else if (startsWith(nm, "consensus_")) "MACS3 >=2-replicate consensus, reference only"
                      else "former H3K9me2 regions, not used since 2026-11-20"))
  peak_files <- c(peak_files, out)
}

# --- session file ------------------------------------------------------------
# Relative paths, so the folder can be moved or sent as a unit.
files <- c("gam_large_blocks.bed", "chip_H4K20me3_log2FC.bedgraph",
           "gam_contacts_lost.bedgraph", "gam_hotspot_domains.bed",
           "gam_melting_score.bedgraph", "chip_compartments.bed")
files <- files[file.exists(file.path(igv_dir, files))]
files <- c(files, unname(pooled), peak_files)
writeLines(c(
  '<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
  '<Session genome="mm39" locus="chr11:1-110000000" version="8">',
  '  <Resources>',
  sprintf('    <Resource path="%s"/>', files),
  '  </Resources>',
  '</Session>'), file.path(igv_dir, "gam_chip_minimal.xml"))

writeLines(c(
  "# Minimal IGV view - GAM/ChIP integration (mm39)",
  "",
  "Open `gam_chip_minimal.xml` in IGV (File > Open Session), or drag the",
  "individual files onto an mm39 genome. One track per data type:",
  "",
  "| Track | What it is |",
  "|---|---|",
  "| H4K20me3 log2FC | the ChIP effect, per 20 kb bin, KO vs WT |",
  "| GAM contacts lost | differential contacts stronger in WT (lost in the KO) |",
  "| Hotspot domains | called remodelled domains, adjusted variant |",
  "| MELTRON signed KS D | + = contact density lost (melting), - = condensing |",
  "| ChIP compartments | mark-loss groups stable / H4K20me3_only / co_loss (not A/B) |",
  "| GAM large blocks | 10-25 Mb one-directional territory (red lost, blue gained) |",
  "| `<mark>_<genotype>_pooled.bw` | mean ChIP signal of the analysed libraries, 50 bp bins (from `ChIP/export_pooled_bigwigs.R`, referenced in place) |",
  "| `peaks/regions_*` | the region sets the ChIP analysis quantifies; H3K9me2/H3K9me3/H4K20me3 share one |",
  "| `peaks/consensus_*` | MACS3 >=2-replicate consensus peaks, reference only |",
  "| `peaks/retired_H3K9me2_merged` | the H3K9me2 regions used until 2026-11-20 |",
  "",
  "Default locus is the whole of chr11, which shows the dominant feature most",
  "clearly: contacts are lost across essentially all of 0-25 Mb and gained",
  "across 55-110 Mb. That large-block structure covers 41.6% of the genome.",
  "",
  "There is no A/B compartment track: NPMI is a normalised quantity and an",
  "eigendecomposition of it does not recover compartments.",
  "",
  "The melting track is the signed Kolmogorov-Smirnov statistic D, NOT the",
  "published -log10(p) melting score. On this data that score saturates and",
  "its |score| > 5 cutoff calls most of the genome; GAM_2c's header explains why.",
  "",
  "Read together: H4K20me3 falls almost everywhere, while contact direction is",
  "organised into very large blocks - the two are uncoupled at bin level",
  "(GAM_3: rho -0.015). Zoom out to see the blocks; zoom in and they dissolve.",
  "Regenerate with `cd GAM && Rscript export_igv.R`."),
  file.path(igv_dir, "README.md"))

message("Session: ", file.path(igv_dir, "gam_chip_minimal.xml"))
