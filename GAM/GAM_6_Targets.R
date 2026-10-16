# ================================================================
# GAM_6 - The derepressed targets: HUSH families, repeats, protocadherin
# ----------------------------------------------------------------
#     cd GAM && Rscript GAM_6_Targets.R              (after GAM_1, GAM_2, GAM_4)
#
# Merges the three stages that ask what happens to the loci the ChIP side found
# derepressed. Previously GAM_5, GAM_9 and GAM_14.
#
#   A (was GAM_5)   All four clustered families and the young repeats, against
#                   composition-matched background bins. Protocadherin was the
#                   only set with a decompaction signature.
#   B (was GAM_9)   Protocadherin at native 20 kb from the contact PAIRS, because
#                   A tested 29 bins with genome-wide parameters and half its
#                   verdict leaned on chr18's compartment axis, which is unusable.
#   C (was GAM_14)  The cluster is approximately one TAD, so: how does that TAD
#                   talk to OTHER TADs? Separation-matched permutation over every
#                   other chr18 TAD, under two independent TAD partitions.
#
# READ A AND C TOGETHER. A found protocadherin uniquely decompacted; C finds the
# Pcdh TAD indistinguishable from other chr18 TADs, intra or distal, under both
# partitions (p 0.67-0.88). Those are compatible - A compared against matched
# background BINS genome-wide, C against neighbouring TADs on the same
# chromosome at matched separation - but any claim about this locus should state
# which comparison it rests on.
#
# Output file names are unchanged (gam5_*, gam9_*, gam14_*).
# ================================================================
source("config.R")
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(GenomicRanges); library(IRanges)
  library(TxDb.Mmusculus.UCSC.mm39.knownGene); library(org.Mm.eg.db)
})





# ==============================================================================
# SECTION A - families and repeats vs matched background   (was GAM_5)
# ==============================================================================

local({
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
      # Restrict to the chromosomes the bin grid actually covers BEFORE the guard
      # below. The complete UCSC rmsk table spans scaffolds and chrX/chrY, none of
      # which exist in the autosome-only GAM grid, so ~10-15% of copies cannot map
      # for entirely legitimate reasons. Counting those as failures would make the
      # seqname guard fire on correct data and hide the thing it is for.
      n_all <- nrow(cp)
      cp <- cp[chr %in% as.character(seqlevels(bins))]
      if (!nrow(cp)) {
        warning(nm, ": no copies on the gridded chromosomes - check seqname style",
                call. = FALSE); next
      }
      # Assign each copy to ONE bin by its midpoint: IAPs are ~7 kb and L1s up to 6 kb,
      # so overlap-assignment would put one copy in 2-3 bins and inflate the set with
      # neighbours that are not really "carrying" the element. rmsk start is 0-based.
      mid <- GRanges(cp$chr, IRanges(start = floor((cp$start + cp$end) / 2) + 1, width = 1))
      h   <- suppressWarnings(findOverlaps(mid, bins, ignore.strand = TRUE))
      # Same guard as for the families, now applied only to copies on gridded
      # chromosomes: a chr-style mismatch would still drive this to zero.
      stopifnot(length(h) / length(mid) > 0.8)
      rep_bins[[nm]] <- sort(unique(subjectHits(h)))
      message(sprintf("  %-12s %6d of %6d copies on gridded chr -> %5d bins",
                      nm, nrow(cp), n_all, length(rep_bins[[nm]])))
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
})


# ==============================================================================
# SECTION B - protocadherin at native resolution           (was GAM_9)
# ==============================================================================

local({
  x <- readRDS(file.path(rds_dir, "gam_bins_annotated.rds"))
  a <- as.data.table(x$annot); setorder(a, chr, start)
  for (mk in colnames(x$l2fc)) a[[paste0("l2fc_", mk)]] <- x$l2fc[match(a$bin_id, rownames(x$l2fc)), mk]

  # --- locate the cluster from the gene models, not from a hard-coded interval --
  txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
  gn <- genes(txdb, single.strand.genes.only = TRUE)
  sym <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, keys = names(gn),
           column = "SYMBOL", keytype = "ENTREZID"))
  pc <- gn[which(grepl("^Pcdh(a|b|g)", ifelse(is.na(sym), "", sym)))]
  pc <- pc[as.character(seqnames(pc)) == "chr18"]
  if (!length(pc)) stop("no chr18 Pcdh genes found", call. = FALSE)
  lo <- min(start(pc)); hi <- max(end(pc))
  PAD <- 2e6
  w0 <- max(0, lo - PAD); w1 <- hi + PAD
  cat(sprintf("=== Pcdh cluster chr18 %.2f-%.2f Mb (%d genes); window %.2f-%.2f Mb ===\n",
              lo/1e6, hi/1e6, length(pc), w0/1e6, w1/1e6))

  inw  <- a[chr == "chr18" & start >= w0 & end <= w1]
  incl <- a[chr == "chr18" & end > lo & start < hi]
  cat(sprintf("  cluster bins: %d   window bins: %d\n", nrow(incl), nrow(inw)))

  # --- ChIP marks over the cluster vs the rest of chr18 and the genome ---------
  cat("\n=== ChIP log2FC: cluster vs chr18 vs genome ===\n")
  cat(sprintf("  %-10s %9s %9s %9s\n", "mark", "cluster", "chr18", "genome"))
  for (mk in colnames(x$l2fc)) {
    v <- paste0("l2fc_", mk)
    cat(sprintf("  %-10s %+9.3f %+9.3f %+9.3f\n", mk,
        median(incl[[v]], na.rm = TRUE),
        median(a[chr == "chr18"][[v]], na.rm = TRUE),
        median(a[[v]], na.rm = TRUE)))
  }

  # --- local contacts from the pairs, at native resolution, by distance --------
  g7 <- file.path(rds_dir, "gam7_local_contacts.rds")
  if (file.exists(g7)) {
    b7 <- as.data.table(readRDS(g7)$bins)
    cols <- grep("^(lost|gained)_", names(b7), value = TRUE)
    m <- merge(a[chr == "chr18", .(bin_id, start, end)],
               b7[, c("bin_id", cols, "local_lost", "local_gained"), with = FALSE], by = "bin_id")
    m[, inclust := start >= lo & end <= hi]
    cat("\n=== local differential contacts, cluster vs rest of chr18 (per bin) ===\n")
    cat(sprintf("  %-14s %10s %10s %8s\n", "band", "cluster", "chr18 rest", "ratio"))
    for (cc in cols) {
      a1 <- mean(m[inclust == TRUE][[cc]]); a2 <- mean(m[inclust == FALSE][[cc]])
      cat(sprintf("  %-14s %10.1f %10.1f %8.2f\n", cc, a1, a2, a1 / max(a2, 1e-9)))
    }
  } else message("GAM_7 output absent - skipping the pairs-based cluster counts")

  # --- native-resolution NPMI across the window, WT vs KO ----------------------
  mats <- list()
  for (g in c("WT", "HP1gKO")) {
    mm <- load_npmi(g, "chr18")
    ids <- intersect(inw$bin_id, rownames(mm))
    mats[[g]] <- mm[ids, ids, drop = FALSE]
    rm(mm); gc(FALSE)
  }
  ids <- intersect(rownames(mats$WT), rownames(mats$HP1gKO))
  W <- mats$WT[ids, ids]; K <- mats$HP1gKO[ids, ids]
  # bin ids are "chr18:36000000-36020000": take the START, not the whole range.
  # sub(".*:", "") alone leaves "36000000-36020000", which coerces to NA and makes
  # every downstream mean NaN.
  pos <- as.numeric(sub("-.*", "", sub(".*:", "", ids)))
  stopifnot(all(is.finite(pos)))
  D <- K - W
  cat(sprintf("\n=== native 20 kb NPMI over the window (%d x %d bins) ===\n", length(ids), length(ids)))
  inc <- pos >= lo & pos <= hi
  sub_mean <- function(M, r, c) mean(M[r, c][upper.tri(M[r, c])], na.rm = TRUE)
  cat(sprintf("  within-cluster   NPMI  WT %.4f  KO %.4f  delta %+.4f (%+.1f%%)\n",
      sub_mean(W, inc, inc), sub_mean(K, inc, inc),
      sub_mean(K, inc, inc) - sub_mean(W, inc, inc),
      100 * (sub_mean(K, inc, inc) - sub_mean(W, inc, inc)) / abs(sub_mean(W, inc, inc))))
  fl <- !inc
  cat(sprintf("  cluster-to-flank NPMI  WT %.4f  KO %.4f  delta %+.4f (%+.1f%%)\n",
      mean(W[inc, fl], na.rm = TRUE), mean(K[inc, fl], na.rm = TRUE),
      mean(K[inc, fl], na.rm = TRUE) - mean(W[inc, fl], na.rm = TRUE),
      100 * (mean(K[inc, fl], na.rm = TRUE) - mean(W[inc, fl], na.rm = TRUE)) / abs(mean(W[inc, fl], na.rm = TRUE))))
  cat(sprintf("  within-flank     NPMI  WT %.4f  KO %.4f  delta %+.4f (%+.1f%%)\n",
      sub_mean(W, fl, fl), sub_mean(K, fl, fl),
      sub_mean(K, fl, fl) - sub_mean(W, fl, fl),
      100 * (sub_mean(K, fl, fl) - sub_mean(W, fl, fl)) / abs(sub_mean(W, fl, fl))))

  dt <- CJ(i = seq_along(ids), j = seq_along(ids))[i < j]
  dt[, `:=`(x = pos[i] / 1e6, y = pos[j] / 1e6, d = D[cbind(i, j)])]
  lim <- quantile(abs(dt$d), 0.98, na.rm = TRUE)
  p <- ggplot(dt[is.finite(d)], aes(x, y, fill = pmax(-lim, pmin(lim, d)))) +
    geom_raster() +
    annotate("rect", xmin = lo/1e6, xmax = hi/1e6, ymin = lo/1e6, ymax = hi/1e6,
             fill = NA, colour = "black", linewidth = 0.5) +
    scale_fill_gradient2(low = heat_low, mid = heat_mid, high = heat_high,
                         midpoint = 0, name = "KO - WT\nNPMI") +
    coord_fixed() +
    labs(title = "Protocadherin cluster, chr18 - native 20 kb NPMI difference",
         subtitle = sprintf("KO minus WT. Box = Pcdha/b/g cluster (%.2f-%.2f Mb). Red = more contact in KO.",
                            lo/1e6, hi/1e6),
         x = "chr18 (Mb)", y = "chr18 (Mb)") +
    theme_m + theme(text = element_text(size = 14))
  save_fig(p, "pcdh_native_npmi_difference", subdir = "pcdh", width = 9, height = 8)

  fwrite(data.table(bin_id = ids, pos = pos, in_cluster = inc,
                    npmi_row_WT = rowMeans(W, na.rm = TRUE),
                    npmi_row_KO = rowMeans(K, na.rm = TRUE)),
         file.path(tables_dir, "gam9_pcdh_native_profile.tsv"), sep = "\t")
  cat("\nGAM_9 complete.\n")
})


# ==============================================================================
# SECTION C - the Pcdh TAD against other TADs              (was GAM_14)
# ==============================================================================

local({
  set.seed(1)
  N_PERM   <- 2000L
  SEP_BRK  <- c(0, 2e6, 5e6, 1e7, 2e7, Inf)
  SEP_NM   <- c("<2Mb", "2-5Mb", "5-10Mb", "10-20Mb", ">20Mb")
  CHR      <- "chr18"
  fig_sub  <- "pcdh_tad"
  msg <- function(...) cat(sprintf(...), "\n", sep = "")

  # --- where is the Pcdh cluster ------------------------------------------------
  pcdh <- local({
    g <- genes(TxDb.Mmusculus.UCSC.mm39.knownGene, single.strand.genes.only = TRUE)
    sym <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, keys = names(g),
             column = "SYMBOL", keytype = "ENTREZID"))
    g <- g[!is.na(sym) & grepl("^Pcdh(a|b|g)", sym)]
    g <- g[as.character(seqnames(g)) == CHR]
    stopifnot(length(g) > 10)
    GRanges(CHR, IRanges(min(start(g)), max(end(g))))
  })
  msg("Pcdh cluster: %s:%.2f-%.2f Mb (%.0f kb)", CHR,
      start(pcdh) / 1e6, end(pcdh) / 1e6, width(pcdh) / 1e3)

  # --- NPMI for the chromosome, once --------------------------------------------
  MW <- load_npmi("WT", CHR); MK <- load_npmi("HP1gKO", CHR)
  ids <- intersect(rownames(MW), rownames(MK))
  MW <- MW[ids, ids]; MK <- MK[ids, ids]
  bpos <- as.numeric(sub("-.*", "", sub(".*:", "", ids)))
  msg("NPMI loaded: %d bins", length(ids))

  # mean NPMI over the bins of one TAD pair, and over a TAD's own interior
  pair_mean <- function(M, i, j) {
    if (!length(i) || !length(j)) return(NA_real_)
    v <- M[i, j, drop = FALSE]
    if (identical(i, j)) v <- v[upper.tri(v)]
    v <- v[is.finite(v)]
    if (length(v) < 5) NA_real_ else mean(v)
  }

  analyse <- function(tads, label) {
    t <- as.data.table(tads)[chr == CHR]
    setorder(t, start)
    if (nrow(t) < 10) { msg("%s: only %d TADs on %s, skipped", label, nrow(t), CHR); return(NULL) }
    # bins of each TAD
    t[, idx := .(list(integer(0)))]
    blist <- lapply(seq_len(nrow(t)), function(k)
      which(bpos >= t$start[k] & bpos < t$end[k]))
    keep <- which(lengths(blist) >= 3)      # a TAD needs interior to have a mean
    t <- t[keep]; blist <- blist[keep]
    t[, mid := (start + end) / 2]
    # which TAD holds the Pcdh cluster: the one with the largest overlap
    ovl <- pmax(0, pmin(t$end, end(pcdh)) - pmax(t$start, start(pcdh)))
    q0  <- which.max(ovl)
    msg("\n=== %s: %d usable TADs on %s; Pcdh sits in %s:%.2f-%.2f Mb (%.0f%% of the cluster) ===",
        label, nrow(t), CHR, CHR, t$start[q0] / 1e6, t$end[q0] / 1e6,
        100 * ovl[q0] / width(pcdh))

    # --- all query x partner pair means ---
    n <- nrow(t)
    grid <- CJ(q = seq_len(n), p = seq_len(n))[q != p]
    grid[, sep := abs(t$mid[q] - t$mid[p])]
    grid[, `:=`(WT = mapply(function(a, b) pair_mean(MW, blist[[a]], blist[[b]]), q, p),
                KO = mapply(function(a, b) pair_mean(MK, blist[[a]], blist[[b]]), q, p))]
    grid <- grid[is.finite(WT) & is.finite(KO) & WT > 0]
    grid[, pct := (KO - WT) / WT * 100]
    grid[, sep_band := cut(sep, SEP_BRK, labels = SEP_NM, right = FALSE)]

    # --- intra-TAD, every TAD ---
    intra <- data.table(
      q = seq_len(n), width = t$end - t$start,
      WT = vapply(seq_len(n), function(k) pair_mean(MW, blist[[k]], blist[[k]]), 0),
      KO = vapply(seq_len(n), function(k) pair_mean(MK, blist[[k]], blist[[k]]), 0))
    intra <- intra[is.finite(WT) & is.finite(KO) & WT > 0][, pct := (KO - WT) / WT * 100]

    # --- the test: separation-stratified mean per query TAD ---
    # Average within separation band first, then across bands, so a TAD is not
    # scored by how many near partners it happens to have.
    per_q <- grid[, .(m = mean(pct)), by = .(q, sep_band)][
               , .(stratified = mean(m), n_bands = .N), by = q][n_bands >= 3]
    obs <- per_q[q == q0, stratified]
    nullv <- per_q[q != q0, stratified]
    if (!length(obs) || length(nullv) < 20) { msg("  too few comparable TADs"); return(NULL) }
    p_emp <- (1 + sum(abs(nullv - mean(nullv)) >= abs(obs - mean(nullv)))) / (1 + length(nullv))
    msg("  Pcdh TAD -> distal, separation-matched: %+.2f%%   (other chr18 TADs: median %+.2f%%, IQR %+.2f to %+.2f)",
        obs, median(nullv), quantile(nullv, .25), quantile(nullv, .75))
    msg("  two-sided empirical p vs %d other TADs = %.3f", length(nullv), p_emp)

    io <- intra[q == q0, pct]; inull <- intra[q != q0, pct]
    p_intra <- if (length(io) && length(inull) >= 20)
      (1 + sum(abs(inull - mean(inull)) >= abs(io - mean(inull)))) / (1 + length(inull)) else NA_real_
    msg("  Pcdh TAD intra: %+.2f%%   (other TADs median %+.2f%%), empirical p = %.3f",
        ifelse(length(io), io, NA), median(inull), p_intra)

    list(
      label = label,
      summary = data.table(
        partition = label, n_tads = nrow(t),
        pcdh_tad = sprintf("%s:%d-%d", CHR, t$start[q0], t$end[q0]),
        pcdh_tad_Mb = (t$end[q0] - t$start[q0]) / 1e6,
        frac_cluster_in_tad = ovl[q0] / width(pcdh),
        intra_pct = ifelse(length(io), io, NA_real_),
        intra_other_median = median(inull), intra_p = p_intra,
        distal_pct = obs, distal_other_median = median(nullv), distal_p = p_emp,
        n_null_tads = length(nullv)),
      by_sep = grid[q == q0, .(pct = mean(pct), n = .N), by = sep_band][
                 , partition := label][],
      by_sep_null = grid[q != q0, .(pct = mean(pct)), by = .(q, sep_band)][
                 , partition := label][],
      profile = grid[q == q0, .(partition = label, sep, pct,
                                partner_mid = t$mid[p])],
      nullv = data.table(partition = label, stratified = nullv),
      obs = obs)
  }

  parts <- list(external = load_tads_external())
  for (g in c("WT", "HP1gKO")) {
    tt <- tryCatch(call_tads_insulation(g), error = function(e) { msg("insulation TADs (%s) failed: %s", g, conditionMessage(e)); NULL })
    if (!is.null(tt)) parts[[paste0("own_", g)]] <- tt
  }
  for (nm in names(parts)) msg("%s: %d TADs genome-wide", nm, nrow(parts[[nm]]))

  R <- Filter(Negate(is.null), lapply(names(parts), function(nm) analyse(parts[[nm]], nm)))
  stopifnot(length(R) > 0)

  SUM <- rbindlist(lapply(R, `[[`, "summary"))
  fwrite(SUM, file.path(tables_dir, "gam14_pcdh_tad_summary.tsv"), sep = "\t")
  msg("\n=== summary across partitions ==="); print(SUM)

  BY <- rbindlist(lapply(R, `[[`, "by_sep"))
  fwrite(BY, file.path(tables_dir, "gam14_pcdh_by_separation.tsv"), sep = "\t")
  PROF <- rbindlist(lapply(R, `[[`, "profile"))
  fwrite(PROF, file.path(tables_dir, "gam14_pcdh_partner_profile.tsv"), sep = "\t")

  # --- figures ------------------------------------------------------------------
  NULLD <- rbindlist(lapply(R, `[[`, "nullv"))
  OBS <- SUM[, .(partition, distal_pct)]
  p <- ggplot(NULLD, aes(stratified)) +
    geom_histogram(bins = 40, fill = "grey80", colour = NA) +
    geom_vline(data = OBS, aes(xintercept = distal_pct), colour = "#b2182b", linewidth = .8) +
    facet_wrap(~ partition, scales = "free_y", nrow = 1) +
    labs(title = "Pcdh TAD vs distal TADs, against every other chr18 TAD at matched separation",
         subtitle = paste0("Grey = the separation-stratified contact change of each other chr18 TAD against its own ",
                           "partners.\nRed = the Pcdh TAD. Averaging within separation band first stops a TAD being ",
                           "scored by how many near\npartners it happens to have."),
         x = "separation-matched change in mean NPMI, HP1gKO vs WT (%)", y = "TADs") +
    theme_m
  save_fig(p, "pcdh_vs_null", subdir = fig_sub, width = 11, height = 4.6)

  BYN <- rbindlist(lapply(R, `[[`, "by_sep_null"))
  bq <- BYN[, .(lo = quantile(pct, .25, na.rm = TRUE), md = median(pct, na.rm = TRUE),
                hi = quantile(pct, .75, na.rm = TRUE)), by = .(partition, sep_band)]
  p <- ggplot() +
    geom_hline(yintercept = 0, linewidth = .3, colour = "grey50") +
    geom_ribbon(data = bq, aes(as.integer(sep_band), ymin = lo, ymax = hi, group = partition),
                fill = "grey85") +
    geom_line(data = bq, aes(as.integer(sep_band), md, group = partition), colour = "grey45") +
    geom_line(data = BY, aes(as.integer(sep_band), pct, group = partition), colour = "#b2182b", linewidth = .8) +
    geom_point(data = BY, aes(as.integer(sep_band), pct), colour = "#b2182b", size = 1.8) +
    facet_wrap(~ partition, nrow = 1) +
    scale_x_continuous(breaks = seq_along(SEP_NM), labels = SEP_NM) +
    labs(title = "Where along the chromosome does the Pcdh TAD change its contacts?",
         subtitle = paste0("Red = the Pcdh TAD. Grey band = interquartile range over all other chr18 TADs at the same ",
                           "separation,\ngrey line their median - so the grey IS the distance control."),
         x = "separation from the Pcdh TAD", y = "% change in mean NPMI") +
    theme_m
  save_fig(p, "pcdh_by_separation", subdir = fig_sub, width = 11, height = 4.6)

  p <- ggplot(PROF, aes(partner_mid / 1e6, pct)) +
    geom_hline(yintercept = 0, linewidth = .3, colour = "grey50") +
    geom_vline(xintercept = (start(pcdh) + end(pcdh)) / 2e6, linetype = 2, colour = "grey40") +
    geom_point(aes(colour = pct > 0), size = 1.3) +
    scale_colour_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "#2166ac"), guide = "none") +
    facet_wrap(~ partition, nrow = 1) +
    labs(title = "Pcdh TAD contact change with each partner TAD along chr18",
         subtitle = "Dashed line = the Pcdh cluster. Each point is one partner TAD; red = more contact in the KO.",
         x = "partner TAD midpoint (Mb)", y = "% change in mean NPMI") +
    theme_m
  save_fig(p, "pcdh_partner_profile", subdir = fig_sub, width = 11, height = 4.6)

  msg("\nGAM_14 complete.")
})
