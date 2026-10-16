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

# --- 4. Architecture: PC1, WT ------------------------------------------------
# One genotype. NA on chromosomes where no balanced compartment eigenvector was
# found - those gaps are real and should stay visible rather than be filled.
g2 <- file.path(rds_dir, "npmi_bin_features.rds")
if (file.exists(g2)) {
  f <- readRDS(g2)$features
  m <- merge(a[, .(bin_id, chr, start, end)], f[, .(bin_id, pc1_WT)], by = "bin_id")
  setorder(m, chr, start)
  d4 <- m[is.finite(pc1_WT), .(chr, start, end, round(pc1_WT, 4))]
  wr(d4, "gam_pc1_WT.bedgraph",
     paste0('track type=bedGraph name="GAM PC1 (WT, A/B axis)" ',
            'description="100kb compartment call expanded to 20kb; axis weakly ',
            'determined on most chromosomes - see GAM/CLAUDE.md" visibility=full ',
            'color=33,102,172 altColor=178,24,43 viewLimits=-2:2 windowingFunction=mean'))
} else message("  (no NPMI features - run GAM_2)")

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
# cannot drift apart. Red = lost in the KO, blue = gained.
blk_src <- file.path(bed_dir, "gam6_large_blocks.bed")
if (file.exists(blk_src)) {
  file.copy(blk_src, file.path(igv_dir, "gam_large_blocks.bed"), overwrite = TRUE)
  message("  gam_large_blocks.bed  (copied from GAM_6)")
} else message("  (no large-block bed - run GAM_6_Large_Blocks.R)")

# --- session file ------------------------------------------------------------
# Relative paths, so the folder can be moved or sent as a unit.
files <- c("gam_large_blocks.bed", "chip_H4K20me3_log2FC.bedgraph",
           "gam_contacts_lost.bedgraph", "gam_hotspot_domains.bed",
           "gam_pc1_WT.bedgraph", "chip_compartments.bed")
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
  "| GAM PC1 (WT) | compartment axis; gaps = no usable axis on that chromosome |",
  "| ChIP compartments | stable / H4K20me3_only / co_loss |",
  "| GAM large blocks | 10-25 Mb one-directional territory (red lost, blue gained) |",
  "",
  "Default locus is the whole of chr11, which shows the dominant feature most",
  "clearly: contacts are lost across essentially all of 0-25 Mb and gained",
  "across 55-110 Mb. That large-block structure covers 41.6% of the genome.",
  "",
  "The PC1 track is NOT reliable genome-wide - only 11 of 38 chromosome x",
  "genotype axes reach |rho| > 0.4. See `../../CLAUDE.md` before using it.",
  "",
  "Read together: H4K20me3 falls almost everywhere, while contact direction is",
  "organised into very large blocks - the two are uncoupled at bin level",
  "(GAM_3: rho -0.015). Zoom out to see the blocks; zoom in and they dissolve.",
  "Regenerate with `cd GAM && Rscript export_igv.R`."),
  file.path(igv_dir, "README.md"))

message("Session: ", file.path(igv_dir, "gam_chip_minimal.xml"))
