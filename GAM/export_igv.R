# ================================================================
# Minimal IGV export - ONE track per data type
# ----------------------------------------------------------------
#   cd GAM && Rscript export_igv.R
#
# Writes results/igv/ : five tracks plus a session file that loads them in a
# sensible order. Deliberately minimal - one representative of each KIND of
# data, not every mark, genotype, direction and variant. The point is a view
# you can open and read, not a complete dump; the rest is in the .rds and the
# tables.
#
#   1 chip_H4K20me3_log2FC.bedgraph  the ChIP effect   (the mark that moves)
#   2 gam_contacts_lost.bedgraph     the GAM effect    (differential contacts)
#   3 gam_hotspot_domains.bed        called domains    (the adjusted variant)
#   4 gam_pc1_WT.bedgraph            architecture      (compartment axis, WT)
#   5 chip_compartments.bed          the ChIP groups   (colour-coded)
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

# --- 4. Architecture: the SHARED A/B axis ------------------------------------
# GAM_2b writes gam_shared_axis.bedgraph and the two per-genotype scores
# directly, so nothing is rebuilt here - this only reports whether they exist.
#
# The per-genotype PC1 track that used to live here has been RETIRED from the
# session. It showed one genotype's own selected eigenvector, and the two
# genotypes select the same eigenvector on only 8 of 19 autosomes, so putting
# WT's axis on screen invited a WT-vs-cKO comparison that is not defined. The
# shared axis is one ruler both genotypes are scored against, and it is better
# determined on 16 of 19 autosomes. gam_pc1_WT.bedgraph is still written by
# GAM_2's own outputs if you want it; it is simply not in the session.
sa <- file.path(igv_dir, "gam_shared_axis.bedgraph")
message(if (file.exists(sa)) "  gam_shared_axis.bedgraph  (from GAM_2b)"
        else "  (no shared axis - run GAM_2b_Shared_Axis.R)")
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

# --- session file ------------------------------------------------------------
# Relative paths, so the folder can be moved or sent as a unit.
files <- c("gam_large_blocks.bed", "chip_H4K20me3_log2FC.bedgraph",
           "gam_contacts_lost.bedgraph", "gam_hotspot_domains.bed",
           "gam_shared_axis.bedgraph", "gam_shared_score_delta.bedgraph",
           "gam_melting_score.bedgraph", "chip_compartments.bed")
files <- files[file.exists(file.path(igv_dir, files))]
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
  "| GAM shared A/B axis | one compartment axis both genotypes are scored on |",
  "| GAM A/B score delta | cKO minus WT on that axis; shape only |",
  "| MELTRON signed KS D | + = contact density lost (melting), - = condensing |",
  "| ChIP compartments | stable / H4K20me3_only / co_loss |",
  "| GAM large blocks | 10-25 Mb one-directional territory (red lost, blue gained) |",
  "",
  "Default locus is the whole of chr11, which shows the dominant feature most",
  "clearly: contacts are lost across essentially all of 0-25 Mb and gained",
  "across 55-110 Mb. That large-block structure covers 41.6% of the genome.",
  "",
  "The compartment track is the SHARED axis from GAM_2b, not a per-genotype PC1.",
  "Each genotype selects its own eigenvector, and they pick the same one on only",
  "8 of 19 autosomes, so a WT-PC1-vs-cKO-PC1 comparison is not defined. The",
  "shared axis is better determined than either genotype's own on 16 of 19",
  "autosomes, but it still only reaches |rho| > 0.4 against active chromatin on",
  "10 of 19 - read it as weak evidence, and see `../../CLAUDE.md` first.",
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
