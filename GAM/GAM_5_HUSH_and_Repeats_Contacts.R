# ================================================================
# GAM_5 - Do the derepressed HUSH / CBX3 targets remodel their 3D contacts?
# ----------------------------------------------------------------
# MINUTE_4/5 show that the HUSH/CBX3-silenced gene families (clustered
# protocadherins, KRAB-ZFP, vomeronasal, olfactory) and young repeats (IAP, young
# L1) lose H3K9me3 + H4K20me3 in the HP1gKO. Here we map those targets onto the
# 20 kb GAM grid and ask whether their bins change contacts / compartment
# relative to MATCHED background bins.
#
# HEADLINE QUESTION: do the derepressed clusters - clustered protocadherins
# especially, which uniquely also GAIN H3K4me3 - DECOMPACT, i.e. lose short-range
# contacts (npmi_short_delta < 0) and shift toward A (pc1_delta > 0)?
#
# Reads results/rds/gam_bins_annotated.rds (GAM_1) and
# results/rds/npmi_bin_features.rds (GAM_2); joins on bin_id. Locus heatmaps read
# NPMI for the one chromosome needed (load_npmi). Runs after GAM_1/GAM_2.
# ================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(ComplexHeatmap) })

gam1_rds <- file.path(rds_dir, "gam_bins_annotated.rds")
gam2_rds <- file.path(rds_dir, "npmi_bin_features.rds")
if (!file.exists(gam1_rds))
  stop("GAM_1 output not found:\n  ", gam1_rds,
       "\nRun it first:   cd GAM && Rscript GAM_1_Bins_and_Signal.R", call. = FALSE)
if (!file.exists(gam2_rds))
  stop("GAM_2 output not found:\n  ", gam2_rds,
       "\nRun it first:   cd GAM && Rscript GAM_2_NPMI_Compartments.R", call. = FALSE)
g1 <- readRDS(gam1_rds); g2 <- readRDS(gam2_rds)

bins <- g1$bins
d <- copy(as.data.table(g1$annot))
stopifnot(identical(d$bin_id, mcols(bins)$bin_id) || setequal(d$bin_id, mcols(bins)$bin_id))
feat <- as.data.table(g2$features)
# 3D features: pc1_delta / insulation_delta already exist; short- and long-range NPMI
# deltas are derived here (HP1gKO - WT, canonical labels - never "cKO").
for (v in c("npmi_short", "npmi_long", "decay_slope"))
  if (all(paste0(v, c("_WT", "_HP1gKO")) %in% names(feat)))
    feat[, (paste0(v, "_delta")) := get(paste0(v, "_HP1gKO")) - get(paste0(v, "_WT"))]
keep <- intersect(c("bin_id", "pc1_delta", "insulation_delta", "npmi_short_delta",
                    "npmi_long_delta", "decay_slope_delta", "compartment_WT",
                    "compartment_switch", "n_valid_WT", "n_valid_HP1gKO"), names(feat))
d <- merge(d, feat[, ..keep], by = "bin_id", all.x = TRUE, sort = FALSE)
setkey(d, NULL)
d <- d[match(mcols(bins)$bin_id, bin_id)]     # restore GRanges order: row i == bins[i]
stopifnot(identical(d$bin_id, mcols(bins)$bin_id))

# ================================================================
# 1. Gene-family bins
# ================================================================
# The four families are defined in ChIP/config.R (gene_families, family_exon_rule,
# family_exons()). We CANNOT source ChIP/config.R: it redefines results_dir/fig_dir/
# tables_dir and would silently redirect every GAM output into ChIP/results. So we
# take the exon GRanges MINUTE_5 already persisted - the identical objects the ChIP
# result was computed on, which is the whole point (no second definition of "a
# family member" to drift out of sync). The local reimplementation below is only a
# fallback for a tree where MINUTE_5 has not been run, and MUST mirror
# ChIP/config.R if those definitions ever change.
fam_rds <- file.path(chip_results_dir, "rds", "family_exon_signal.rds")
if (file.exists(fam_rds)) {
  fam_gr <- readRDS(fam_rds)$fam_gr                 # NCBI seqnames ("18")
  message("Gene families: using MINUTE_5's persisted exons (", fam_rds, ")")
} else {
  message("Gene families: ", fam_rds, " absent - rebuilding with a LOCAL copy of ",
          "ChIP/config.R::family_exons(). Run ChIP MINUTE_5 to use the canonical set.")
  local_gene_families <- list(Protocadherin = "^Pcdh(a|b|g)", `KRAB-ZFP` = "^Zfp",
                              Vomeronasal = "^Vmn", Olfactory = "^(Olfr|Or[0-9])")
  local_exon_rule <- list(Protocadherin = "first", `KRAB-ZFP` = "last",
                          Vomeronasal = "largest", Olfactory = "largest")
  local_family_exons <- function() {
    txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
    exg  <- exonsBy(txdb, by = "gene")
    sym  <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, keys = names(exg),
              column = "SYMBOL", keytype = "ENTREZID"))
    sym  <- ifelse(is.na(sym), "", sym)
    pick <- function(ex, rule) {
      if (rule == "largest") return(ex[which.max(width(ex))])
      ex <- ex[if (as.character(strand(ex))[1] == "-") order(-start(ex)) else order(start(ex))]
      if (rule == "first") ex[1] else ex[length(ex)]
    }
    parts <- lapply(names(local_gene_families), function(fn) {
      idx <- which(grepl(local_gene_families[[fn]], sym))
      if (!length(idx)) return(NULL)
      do.call(c, lapply(idx, function(i) {
        g <- pick(exg[[i]], local_exon_rule[[fn]]); mcols(g) <- NULL
        mcols(g)$family <- fn; mcols(g)$gene <- sym[i]; g }))
    })
    do.call(c, Filter(Negate(is.null), parts))
  }
  fam_gr <- local_family_exons()                    # UCSC already
}
fam_gr <- to_ucsc(fam_gr)                           # persisted = NCBI, bins = UCSC
fam_names <- c("Protocadherin", "KRAB-ZFP", "Vomeronasal", "Olfactory")
fam_names <- intersect(fam_names, unique(as.character(mcols(fam_gr)$family)))

# exon -> containing 20 kb bin(s). An exon straddling a boundary lands in both bins;
# a target SET is the unique bins, so a 30-gene cluster inside 10 bins is 10 bins,
# not 30 (counting genes would pseudo-replicate neighbouring bins that share contacts).
ov <- suppressWarnings(findOverlaps(fam_gr, bins, ignore.strand = TRUE))
exon_rate <- length(unique(queryHits(ov))) / length(fam_gr)
# Not 1.0: the grid is autosomal and several Zfp genes sit on chrX. A botched
# NCBI/UCSC flip gives ~0, which would read downstream as a clean null.
stopifnot(length(ov) > 0, exon_rate > 0.7,
          all(vapply(fam_names, function(f)
            any(as.character(mcols(fam_gr)$family[queryHits(ov)]) == f), logical(1))))
message(sprintf("Family exons overlapping a GAM bin: %.1f%% (%d / %d)",
                100 * exon_rate, length(unique(queryHits(ov))), length(fam_gr)))
fam_bins <- lapply(setNames(fam_names, fam_names), function(f) {
  h <- ov[as.character(mcols(fam_gr)$family[queryHits(ov)]) == f]
  sort(unique(subjectHits(h)))                      # row indices into d / bins
})

# ================================================================
# 2. Repeat-subfamily bins
# ================================================================
# IAP / young L1 subfamilies - the ones carrying the ChIP result (MINUTE_4: deepest
# H4K20me3 loss over IAPEz-int, IAPLTR1a_Mm, L1MdT, L1MdGf, L1MdA). The rmsk-format
# LINE/LTR BEDs are external ChIP inputs (not committed); skip the repeat sets with a
# message if they are absent. Read with fread, keeping only the 5 columns we use.
# annot_dir comes from config.R (the single source of truth for paths).
read_rmsk <- function(f) {
  # UCSC rmsk dump: 1-2 '#filter' comment lines, then a '#bin ...' header. fread's
  # auto-detection trips on the comments, so anchor on the header text instead.
  x <- fread(f, sep = "\t", skip = "swScore", header = TRUE, quote = "",
             select = c("genoName", "genoStart", "genoEnd", "repName"), fill = TRUE)
  setnames(x, c("chr", "start", "end", "repName"))
  x[is.finite(start) & is.finite(end)]
}
rep_defs <- list(IAPEz_int = c("LTR", "IAPEz-int"), IAPLTR1a_Mm = c("LTR", "IAPLTR1a_Mm"),
                 L1MdT = c("LINE", "^L1MdT"), L1MdGf = c("LINE", "^L1MdGf"),
                 L1MdA = c("LINE", "^L1MdA"))
rep_files <- c(LINE = file.path(annot_dir, "LINE.mm39.bed"), LTR = file.path(annot_dir, "LTR.mm39.bed"))
rep_bins <- list()
for (cls in c("LINE", "LTR")) {
  if (!file.exists(rep_files[[cls]])) {
    message("GAM_5: ", rep_files[[cls]], " not found - skipping ", cls, " subfamilies",
            " (download the ChIP annotation BEDs, see ChIP/README.md).")
    next
  }
  rm_dt <- read_rmsk(rep_files[[cls]])
  for (nm in names(rep_defs)[vapply(rep_defs, `[`, "", 1) == cls]) {
    cp <- rm_dt[grepl(rep_defs[[nm]][2], repName)]
    if (!nrow(cp)) next
    # Assign each copy to ONE bin by its midpoint: IAPs are ~7 kb and L1s up to 6 kb,
    # so overlap-assignment would put one copy in 2-3 bins and inflate the set with
    # neighbours that are not really "carrying" the element. rmsk start is 0-based.
    mid <- GRanges(cp$chr, IRanges(start = floor((cp$start + cp$end) / 2) + 1, width = 1))
    h   <- suppressWarnings(findOverlaps(mid, bins, ignore.strand = TRUE))
    # same guard as for the families: a chr-style mismatch would give zero here
    stopifnot(length(h) / length(mid) > 0.8)
    rep_bins[[nm]] <- sort(unique(subjectHits(h)))
    message(sprintf("  %-12s %6d copies -> %5d bins", nm, nrow(cp), length(rep_bins[[nm]])))
  }
}
target_sets <- c(fam_bins, rep_bins)
target_type <- setNames(c(rep("gene family", length(fam_bins)), rep("repeat", length(rep_bins))),
                        names(target_sets))

# ================================================================
# 3. Matched background
# ================================================================
# A naive comparison of target bins vs ALL bins is confounded: olfactory/vomeronasal/
# KRAB-ZFP/Pcdh loci and IAP/L1 copies all sit in particular chromatin states, are
# repeat-rich, and differ in how well the bin is covered, and each of those alone
# shifts contact counts and NPMI. A raw difference would partly measure "this is
# quiescent/repeat-dense chromatin", not "this is a HUSH target". So each target bin
# is matched to background bins in the SAME stratum of
#   dominant ChromHMM state  x  repeat-content quantile  x  coverage quantile
# where repeat content = frac_LINE+frac_SINE+frac_LTR and coverage = valid NPMI cells
# per bin (GAM_2 n_valid_*; low-coverage bins have noisier contacts and delta values).
# The null redraws, per permutation, as many background bins per stratum as there are
# target bins in it (candidates exclude the target set itself), and recomputes the
# mean of each metric. We report the EFFECT SIZE (obs - matched mean, and z); the
# data are 2 pooled replicates per genotype with no within-genotype variance, so a
# bin-level test would treat 123k non-independent bins as replicates - only the
# matched-background permutation is defensible.
hmm_cols <- grep("^hmm_", names(d), value = TRUE)
state <- if (length(hmm_cols)) {
  m <- as.matrix(d[, ..hmm_cols]); m[!is.finite(m)] <- -1
  ifelse(rowSums(m >= 0) == 0, "none", sub("^hmm_", "", hmm_cols[max.col(m, ties.method = "first")]))
} else rep("none", nrow(d))
qcut <- function(x, k) {
  x[!is.finite(x)] <- median(x, na.rm = TRUE)
  b <- unique(quantile(x, seq(0, 1, length.out = k + 1), na.rm = TRUE))
  if (length(b) < 3) return(factor(rep("all", length(x))))
  cut(x, b, include.lowest = TRUE)
}
d[, rep_cov := rowSums(.SD, na.rm = TRUE), .SDcols = c("frac_LINE", "frac_SINE", "frac_LTR")]
cov_v <- if (all(c("n_valid_WT", "n_valid_HP1gKO") %in% names(d))) d$n_valid_WT + d$n_valid_HP1gKO else rep(1, nrow(d))
strata <- paste(state, as.integer(qcut(d$rep_cov, 4)), as.integer(qcut(cov_v, 3)), sep = "|")

# Metrics. Direction (config.R): selWT = stronger in WT = LOST in HP1gKO. All deltas
# are HP1gKO - WT. Counts are heavy-tailed, so compare log1p.
d[, `:=`(lost_log1p = log1p(n_contacts_lost), gained_log1p = log1p(n_contacts_gained))]
d[, net_remodel := gained_log1p - lost_log1p]        # >0 : more contacts gained than lost
d[, top_lost   := as.numeric(n_contacts_lost   >= quantile(n_contacts_lost,   0.9, na.rm = TRUE))]
d[, top_gained := as.numeric(n_contacts_gained >= quantile(n_contacts_gained, 0.9, na.rm = TRUE))]
if ("compartment_WT" %in% names(d)) d[, in_A_WT := as.numeric(compartment_WT == "A")]
if ("compartment_switch" %in% names(d)) {
  d[, switch_AB := as.numeric(compartment_switch == "A_B")]
  d[, switch_BA := as.numeric(compartment_switch == "B_A")]
}
metric_cols <- intersect(c("lost_log1p", "gained_log1p", "net_remodel", "top_lost", "top_gained",
                           "pc1_delta", "npmi_short_delta", "npmi_long_delta", "decay_slope_delta",
                           "insulation_delta", "in_A_WT", "switch_AB", "switch_BA"), names(d))
M <- as.matrix(d[, ..metric_cols])

NPERM <- 1000    # per set; stratified draws, p resolution 1e-3
match_perm <- function(target, seed = 42) {
  cand <- setdiff(seq_len(nrow(d)), target)
  by_s <- split(cand, strata[cand])
  need <- table(strata[target])
  need <- need[names(need) %in% names(by_s) & need > 0]
  n_matched <- sum(need)
  obs <- colMeans(M[target[strata[target] %in% names(need)], , drop = FALSE], na.rm = TRUE)
  set.seed(seed)
  perm <- t(vapply(seq_len(NPERM), function(i) {
    smp <- unlist(lapply(names(need), function(s) {
      pool <- by_s[[s]]; pool[sample.int(length(pool), need[[s]], replace = length(pool) < need[[s]])] }))
    colMeans(M[smp, , drop = FALSE], na.rm = TRUE)
  }, numeric(length(metric_cols))))
  expm <- colMeans(perm); sdp <- apply(perm, 2, sd)
  pv <- vapply(seq_along(obs), function(j)
    (sum(abs(perm[, j] - expm[j]) >= abs(obs[j] - expm[j])) + 1) / (NPERM + 1), numeric(1))
  data.table(metric = metric_cols, n_target_bins = length(target), n_matched = n_matched,
             observed = obs, matched_mean = expm, diff = obs - expm,
             z = (obs - expm) / sdp, perm_p = pv,
             naive_all_bins_mean = colMeans(M[cand, , drop = FALSE], na.rm = TRUE))
}
res <- rbindlist(lapply(names(target_sets), function(s) {
  cbind(data.table(set = s, type = target_type[[s]]), match_perm(target_sets[[s]]))
}))
res[, naive_diff := observed - naive_all_bins_mean]   # the confounded comparison, for contrast
res[, p_adj_BH := p.adjust(perm_p, "BH"), by = metric]
res[, set := factor(set, levels = names(target_sets))]
setorder(res, set, metric)
fwrite(res, file.path(tables_dir, "gam5_target_vs_matched_background.tsv"), sep = "\t")
fwrite(data.table(set = names(target_sets), type = target_type, n_bins = lengths(target_sets)),
       file.path(tables_dir, "gam5_target_sets.tsv"), sep = "\t")

# size_perm() cross-check for the two binary headline flags. size_perm() matches on ONE
# numeric covariate (it was written for width), so this matches on repeat content only;
# it is the coarser of the two nulls and is here to show the stratified result is not
# an artefact of the stratification. A tiny jitter breaks the ties at 0.
sp <- rbindlist(lapply(names(target_sets), function(s) {
  rbindlist(lapply(intersect(c("top_lost", "switch_AB"), names(d)), function(fl) {
    ok <- is.finite(d[[fl]]) & is.finite(d$rep_cov)
    qi <- cumsum(ok)[target_sets[[s]][ok[target_sets[[s]]]]]
    set.seed(42); cj <- d$rep_cov[ok] + runif(sum(ok), 0, 1e-6)
    if (length(qi) < 10) return(NULL)
    cbind(data.table(set = s, flag = fl, n_bins = length(qi)),
          size_perm(mean(d[[fl]][ok][qi]), cj[qi], cj, d[[fl]][ok], nperm = 1000))
  }))
}))
fwrite(sp, file.path(tables_dir, "gam5_size_perm_binary_flags.tsv"), sep = "\t")

# ================================================================
# 4. Headline: do the derepressed clusters DECOMPACT?
# ================================================================
# Decompaction = short-range contacts LOST (npmi_short_delta below matched) AND
# compartment shift toward A (pc1_delta above matched), each BH < 0.05. A mixed or
# null outcome is reported as such - this is a test, not a label.
hl <- dcast(res[metric %in% c("npmi_short_delta", "pc1_delta", "net_remodel")],
            set + type ~ metric, value.var = c("diff", "p_adj_BH"))
if (all(c("diff_npmi_short_delta", "diff_pc1_delta") %in% names(hl))) {
  hl[, loses_short_range := diff_npmi_short_delta < 0 & p_adj_BH_npmi_short_delta < 0.05]
  hl[, shifts_toward_A   := diff_pc1_delta        > 0 & p_adj_BH_pc1_delta        < 0.05]
  hl[, verdict := fifelse(loses_short_range & shifts_toward_A, "DECOMPACTS (short-range loss + toward A)",
                  fifelse(loses_short_range, "loses short-range contacts, no A shift",
                  fifelse(shifts_toward_A, "shifts toward A, no short-range loss", "no decompaction signal")))]
}
fwrite(hl, file.path(tables_dir, "gam5_decompaction_verdict.tsv"), sep = "\t")
cat("\n=== HEADLINE: do the derepressed clusters decompact? (vs matched background) ===\n")
print(hl, row.names = FALSE)
cat("\nPer-set effects (obs - matched mean; naive_diff = the confounded all-bin contrast):\n")
print(res[metric %in% c("lost_log1p", "gained_log1p", "pc1_delta", "npmi_short_delta"),
          .(set, metric, n_target_bins, diff = round(diff, 4), z = round(z, 2),
            perm_p, p_adj_BH = round(p_adj_BH, 4), naive_diff = round(naive_diff, 4))],
      row.names = FALSE)

# ================================================================
# 5. Figures -> figures/hush_contacts
# ================================================================
sub_dir <- "hush_contacts"
shown <- c(metric_cols[metric_cols %in% c("lost_log1p", "gained_log1p", "net_remodel", "pc1_delta",
                                          "npmi_short_delta", "npmi_long_delta", "insulation_delta",
                                          "switch_AB", "switch_BA")])
pd <- res[metric %in% shown]
pd[, `:=`(metric = factor(metric, levels = shown), star = ifelse(p_adj_BH < 0.05, "*", ""))]
g_z <- ggplot(pd, aes(metric, set, fill = pmax(pmin(z, 6), -6))) +
  geom_tile(colour = "white") +
  geom_text(aes(label = paste0(sprintf("%.1f", z), star)), size = 3) +
  scale_fill_heat0_div(name = "z vs\nmatched") +
  labs(title = "Contact remodelling at HUSH/CBX3 targets vs matched background bins",
       subtitle = "z of (target mean - matched-background mean); * = BH < 0.05; deltas are HP1gKO - WT; lost/gained from selWT/selcKO",
       x = NULL, y = NULL) +
  theme_m + theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_fig(g_z, "target_vs_matched_z", sub_dir, width = 9, height = 5)

if (all(c("diff_npmi_short_delta", "diff_pc1_delta") %in% names(hl))) {
  g_q <- ggplot(hl, aes(diff_pc1_delta, diff_npmi_short_delta, colour = type, label = set)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55") +
    geom_point(size = 3) + geom_text(vjust = -0.9, size = 3, show.legend = FALSE) +
    scale_colour_disc(name = NULL) +
    labs(title = "Decompaction test: A shift vs short-range contact change",
         subtitle = "decompaction = bottom-right (PC1 toward A, short-range NPMI down); both axes vs matched background",
         x = "PC1 delta, target - matched", y = "short-range NPMI delta, target - matched") +
    theme_m
  save_fig(g_q, "decompaction_quadrant", sub_dir, width = 7, height = 5)
}

# Raw distributions: targets vs one matched draw vs all bins - shows what matching removed.
set.seed(7)
dist_dt <- rbindlist(lapply(names(target_sets), function(s) {
  tg  <- target_sets[[s]]; cand <- setdiff(seq_len(nrow(d)), tg); by_s <- split(cand, strata[cand])
  need <- table(strata[tg]); need <- need[names(need) %in% names(by_s)]
  bg  <- unlist(lapply(names(need), function(k) { p <- by_s[[k]]; p[sample.int(length(p), need[[k]], replace = length(p) < need[[k]])] }))
  rbind(data.table(set = s, source = "target", v = d$lost_log1p[tg]),
        data.table(set = s, source = "matched", v = d$lost_log1p[bg]),
        data.table(set = s, source = "all bins", v = d$lost_log1p[sample(cand, min(5000, length(cand)))]))
}))
dist_dt[, `:=`(set = factor(set, levels = names(target_sets)),
               source = factor(source, levels = c("target", "matched", "all bins")))]
g_d <- ggplot(dist_dt, aes(set, v, fill = source)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.15, linewidth = 0.3, position = position_dodge(0.8)) +
  scale_fill_disc(name = NULL) +
  labs(title = "Contacts LOST in HP1gKO (selWT) at targets, matched background, all bins",
       subtitle = "a naive target-vs-all-bins gap that closes against 'matched' is composition, not HUSH biology",
       x = NULL, y = "log1p(contacts lost per bin)") +
  theme_m + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(g_d, "contacts_lost_target_matched_all", sub_dir, width = 9, height = 4.8)

# ---- Per-locus NPMI heatmaps: WT | HP1gKO | delta -------------------------------
# One chromosome at a time (load_npmi), subset to the locus, drop the big matrix
# before loading the next genotype. NPMI is >0 for co-segregating windows; a shared
# 0..q99 sequential scale across genotypes keeps WT and KO comparable (separate
# scales would make a uniform loss look like nothing). Delta is diverging, symmetric.
locus_heatmap <- function(locus, label, file) {
  ch <- as.character(seqnames(locus))
  sub <- list()
  for (g in names(gam_genotypes)) {
    m  <- load_npmi(g, ch)
    gr <- parse_gam_bin(rownames(m))
    k  <- which(overlapsAny(gr, locus))
    sub[[g]] <- m[k, k, drop = FALSE]
    rm(m); gc(verbose = FALSE)
  }
  stopifnot(nrow(sub$WT) > 3, identical(rownames(sub$WT), rownames(sub$HP1gKO)))
  loc_gr <- parse_gam_bin(rownames(sub$WT))
  has_gene <- overlapsAny(loc_gr, fam_gr, ignore.strand = TRUE)
  ann <- HeatmapAnnotation(`family gene` = has_gene, which = "column",
                           col = list(`family gene` = c(`TRUE` = disc_pal(5)[2], `FALSE` = "grey90")),
                           show_legend = FALSE, show_annotation_name = FALSE, simple_anno_size = unit(2, "mm"))
  top <- unname(quantile(unlist(sub), 0.99, na.rm = TRUE))
  cf  <- heat_col_fun(seq(0, top, length.out = length(heat_cols)))
  mk  <- function(m, nm, ttl, col, show_ann, leg = TRUE) Heatmap(
    m, name = nm, col = col, na_col = "grey95", cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE, show_column_names = FALSE, column_title = ttl,
    top_annotation = if (show_ann) ann else NULL, use_raster = TRUE, raster_device = ht_raster_device,
    show_heatmap_legend = leg, heatmap_legend_param = list(title = nm))
  dl  <- sub$HP1gKO - sub$WT; dm <- max(abs(dl), na.rm = TRUE) * 0.6
  hl_ <- mk(sub$WT, "NPMI", "WT", cf, TRUE) + mk(sub$HP1gKO, "NPMI ", "HP1gKO", cf, TRUE, leg = FALSE) +
    mk(dl, "delta", "HP1gKO - WT", heat_col_fun(seq(-dm, dm, length.out = length(heat_cols))), FALSE)
  save_base_fig(function() draw(hl_, column_title = sprintf("%s  %s:%d-%d", label, ch, start(locus), end(locus)),
                                merge_legends = TRUE),
                file, sub_dir, width = 13, height = 5, dpi = 200)
}
fam_u <- fam_gr[as.character(seqnames(fam_gr)) %in% gam_autosomes]
pad <- 10L * gam_binsize                                  # 200 kb context either side
locus_of <- function(gr) {
  r <- range(gr); r <- resize(r, width(r) + 2 * pad, fix = "center")
  GRanges(seqnames(r), IRanges(max(1L, start(r)), end(r)))
}
pcdh <- fam_u[mcols(fam_u)$family == "Protocadherin" & as.character(seqnames(fam_u)) == "chr18"]
loci <- list()
if (length(pcdh)) loci[["Pcdh_chr18"]] <- list(locus_of(pcdh), "Clustered protocadherins")
kz <- fam_u[mcols(fam_u)$family == "KRAB-ZFP"]
if (length(kz)) {
  # densest KRAB-ZFP cluster = the merged window (gaps < 250 kb) holding the most genes
  rr <- reduce(kz, min.gapwidth = 250e3, ignore.strand = TRUE)
  rr$n <- countOverlaps(rr, kz, ignore.strand = TRUE)
  best <- rr[which.max(rr$n)]
  loci[["KRABZFP_densest"]] <- list(locus_of(best), sprintf("KRAB-ZFP cluster (%d genes)", best$n))
}
for (nm in names(loci)) {
  tryCatch(locus_heatmap(loci[[nm]][[1]], loci[[nm]][[2]], paste0("npmi_locus_", nm)),
           error = function(e) message("GAM_5: locus heatmap ", nm, " skipped: ", conditionMessage(e)))
}

message("GAM_5 done: tables in ", tables_dir, ", figures in ", file.path(fig_dir, sub_dir))
