# ================================================================
# GAM pipeline - libraries, path builders, loaders and helper functions
# ----------------------------------------------------------------
# Sourced by config.R (last line), so a stage only ever needs
#     source("config.R")
# This file holds CODE ONLY. Every value it uses (paths, cutoffs, parameter
# lists) is defined in config.R; nothing here should be edited to change an
# analysis. Cross-assay helpers (palettes, save_fig, hyper_row, size_perm,
# seqname flips, mm39_bins, GAM bin-id parsing) live in ../shared/common.R.
# ================================================================

# --- Dependencies + cross-assay helpers --------------------------
source("../shared/setup.R")
# GAM-only dependency: every GAM input is gzipped, and data.table::fread()
# delegates .gz decompression to R.utils. Without it every loader below fails
# with "fread() requires 'R.utils'". The ChIP assay reads plain BEDs/bigWigs, so
# this belongs to GAM rather than the shared manifest.
minute_add_deps(cran = c("R.utils", "hexbin", "patchwork"))
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

# --- Output tree -------------------------------------------------
for (d in c(counts_dir, rds_dir, tables_dir, fig_dir, bed_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# --- Derived values that need the palette helpers ----------------
# Colours for the six level-2 chromatin classes. From the shared discrete
# palette, not raw hex, per the repo rule - with the deliberate grey exception
# for Quiescent, which is the "unmarked / not-a-category" class the greys exist
# for. One definition so the region track and the exported BED cannot drift.
gam_class2_colours <- local({
  cats <- c("TSS", "Enhancer", "Tx", "Polycomb", "H3K9me3")
  p <- tryCatch(disc_pal(length(cats)), error = function(e) NULL)
  if (is.null(p)) p <- c("#B2182B", "#EF8A62", "#4D9221", "#762A83", "#2166AC")
  c(stats::setNames(p, cats), Quiescent = "grey80")
})

# --- Path builders -----------------------------------------------
# Genotype is encoded ONLY by the directory name (WT_at20Kb / cKO_at20Kb); the
# NPMI file names are identical between genotypes, so never resolve one by
# basename. `sel` is "lost"/"gained" (canonical) or a raw "selWT"/"selcKO".
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
require_chip_rds <- function() {
  if (!file.exists(chip_rds))
    stop("ChIP handoff not found:\n  ", chip_rds,
         "\nThe GAM stages need the annotated peak-level results from the ChIP",
         " pipeline.\nRun it first:   cd ../ChIP && Rscript run_MINUTE.R",
         "\nOr point CHIP_RESULTS at an existing ChIP results tree.", call. = FALSE)
  chip_rds
}

# --- Sample sheet + exclusion (mirrors ChIP/config.R) ------------
# Default exclusion is HP1gKO rep2, as on the ChIP side (see config.R).
# MINUTE_EXCLUDE: comma-separated sample_ids and/or <genotype>:<replicate>
# pairs; "none" keeps every library.
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

# --- The shared coordinate system --------------------------------
# The 20 kb grid that BOTH assays are projected onto. mm39_bins() comes from
# shared/common.R; UCSC seqnames, autosomes by default (matching the GAM
# differential products). This is the join key for the whole integration.
gam_bins <- function(chroms = gam_autosomes) mm39_bins(gam_binsize, chroms)

# --- Annotation loaders ------------------------------------------
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

load_ctcf <- function() {
  if (!file.exists(ctcf_hip_bed) || !file.exists(ctcf_ctx_bed))
    stop("CTCF BEDs missing - run: Rscript prepare_CTCF_annotation.R", call. = FALSE)
  # One GRanges carrying a `tissue` column, not a list - same shape as
  # load_lads(), because every consumer wants to overlap both against the same
  # bins and keep the label. Levels fixed so drop = FALSE scales keep both keys.
  a <- rtracklayer::import(ctcf_hip_bed); a$tissue <- "hippocampus"
  b <- rtracklayer::import(ctcf_ctx_bed); b$tissue <- "cortex"
  g <- c(a, b)
  g$tissue <- factor(g$tissue, levels = c("hippocampus", "cortex"))
  g
}

# --- TAD partitions ----------------------------------------------
# Two independent partitions, kept side by side on purpose (PARAMETERS.md):
# a published boundary set (external) and our own insulation minima (own).
#
# External: intervals are the spans BETWEEN consecutive boundaries. The
# deposited boundaries are ~150 kb wide, so a TAD runs from the end of one to
# the start of the next; the boundary bodies themselves belong to neither TAD.
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

# Own: insulation is mean NPMI ACROSS the boundary at bin i (GAM_2), so a LOW
# value means few cross-boundary contacts, i.e. a strong boundary - boundaries
# are local MINIMA. Prominence is measured in local sd so one threshold works on
# chromosomes with different contact density. Parameters: gam_tad_insul.
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

# --- Heterochromatin domains, merged -----------------------------
# Gap-merges the H3K9me3/H4K20me3 peak fragments into domains. Why a peak set
# needs merging before it is a domain set, and how gap / min_width were chosen:
# PARAMETERS.md "Heterochromatin domains".
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

# --- Chromatin classes (two levels, one definition) --------------
# Returns one row per input bin: bin_id, het_frac, class2, class1.
# `annot` must be the gam_bins_annotated.rds $annot table (hmm_* coverage
# columns plus chr/start/end/bin_id). Pure function of its input, so every
# stage that calls it gets the same partition - which is the whole point of it
# living here rather than being rebuilt in each stage, as is_het was.
# The state -> class mapping and the peak override are in config.R; why each
# state is assigned as it is: PARAMETERS.md "Chromatin classes".
gam_chromatin_classes <- function(annot, peak_cut = gam_het_peak_cut) {
  a <- data.table::as.data.table(annot)
  stopifnot(all(c("bin_id", "chr", "start", "end") %in% names(a)))
  hm <- grep("^hmm_", names(a), value = TRUE)
  st <- sub("^hmm_", "", hm)
  unmapped <- setdiff(st, names(gam_chromatin_state_map))
  if (length(unmapped))
    stop("ChromHMM states with no level-2 assignment: ", paste(unmapped, collapse = ", "),
         "\n  Add them to gam_chromatin_state_map in config.R.", call. = FALSE)

  # Coverage of each level-2 group = sum of its states' coverage fractions.
  grp <- gam_chromatin_state_map[st]
  M <- as.matrix(a[, ..hm])
  G <- vapply(gam_class2_levels, function(g) {
    k <- which(grp == g)
    if (!length(k)) rep(0, nrow(M)) else rowSums(M[, k, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(M)))

  # Peak-based H3K9me3 coverage, the same geometry GAM_5 used for is_het.
  dm <- data.table::fread(het_domain_bed, header = FALSE, select = 1:3,
                          col.names = c("chr", "start", "end"))
  dm[, chr := ifelse(grepl("^chr", chr), chr, paste0("chr", chr))]
  dm <- dm[chr %in% unique(a$chr)]
  red <- GenomicRanges::reduce(
           GenomicRanges::GRanges(dm$chr, IRanges::IRanges(dm$start + 1, dm$end)),
           ignore.strand = TRUE)
  b <- GenomicRanges::GRanges(a$chr, IRanges::IRanges(a$start + 1, a$end))
  h <- GenomicRanges::findOverlaps(b, red, ignore.strand = TRUE)
  w <- GenomicRanges::width(GenomicRanges::pintersect(
         b[S4Vectors::queryHits(h)], red[S4Vectors::subjectHits(h)], ignore.strand = TRUE))
  v <- numeric(length(b))
  agg <- tapply(w, S4Vectors::queryHits(h), sum)
  v[as.integer(names(agg))] <- as.numeric(agg)
  het_frac <- pmin(v / GenomicRanges::width(b), 1)

  cls2 <- gam_class2_levels[max.col(G, ties.method = "first")]
  # A bin with no ChromHMM coverage at all has no dominant group; max.col would
  # hand it column 1 (TSS) on a row of zeros, which would invent active
  # chromatin out of missing annotation.
  cls2[rowSums(G, na.rm = TRUE) <= 0] <- NA_character_
  cls2[het_frac >= peak_cut] <- "H3K9me3"

  out <- data.table::data.table(
    bin_id   = a$bin_id,
    het_frac = het_frac,
    class2   = factor(cls2, levels = gam_class2_levels),
    class1   = factor(unname(gam_chromatin_level1[cls2]), levels = gam_class1_levels))
  # A seqname-style mismatch zeroes the peak overlap and every other failure
  # mode here empties a class, so assert rather than warn (repo rule).
  stopifnot(sum(out$class2 == "H3K9me3", na.rm = TRUE) > 1000,
            sum(out$class1 == "euchromatin", na.rm = TRUE) > 1000)
  out[]
}

# Writes the level-2 chromatin partition as an itemRgb BED for IGV. Merges
# adjacent bins of the same class so the file is runs, not 123,133 single bins.
write_chromatin_state_bed <- function(annot, path) {
  cc <- gam_chromatin_classes(annot)
  a <- data.table::as.data.table(annot)[, .(bin_id, chr, start, end)]
  a <- merge(a, cc[, .(bin_id, class2)], by = "bin_id")
  a <- a[!is.na(class2)]
  data.table::setorder(a, chr, start)
  # rleid over (chr, class) collapses runs; a gap in bins does not merge because
  # the runs are then written with their own start/end.
  a[, run := data.table::rleid(chr, class2)]
  r <- a[, .(chr = chr[1], start = min(start), end = max(end),
             name = as.character(class2[1])), by = run]
  rgb <- grDevices::col2rgb(gam_class2_colours[r$name])
  r[, `:=`(score = 0L, strand = ".", ts = start, te = end,
           itemRgb = paste(rgb[1, ], rgb[2, ], rgb[3, ], sep = ","))]
  writeLines(paste0('track name="GAM chromatin states (level 2)" ',
                    'description="TSS/Enhancer/Tx = euchromatin; ',
                    'Polycomb/H3K9me3/Quiescent = heterochromatin" itemRgb=On'), path)
  data.table::fwrite(r[, .(chr, start, end, name, score, strand, ts, te, itemRgb)],
                     path, sep = "\t", col.names = FALSE, append = TRUE, quote = FALSE)
  invisible(r)
}

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

# --- Hotspot DOMAINS: the second analysis level ------------------
# Returns a list: $bin (per-bin domain_id + smoothed value + hot flag, in input
# order) and $domains (one row per kept domain). Input must be sorted by chr,start.
# Smooth (running mean, within chromosome) -> threshold at quantile q -> absorb
# sub-threshold gaps of <= `gap` bins -> keep runs of >= min_bins. Why smooth
# before thresholding, and why both a smoothed and a raw variant exist:
# PARAMETERS.md "Hotspot domains".
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
# examples/npmi_minimal_example.R walks through this and GAM_2's use of it.
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
