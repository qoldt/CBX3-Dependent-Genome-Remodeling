# ================================================================
# MINUTE_6 — KAP1 / TRIM28 intersection with the differential-loss domains
# ----------------------------------------------------------------
# Does the H4K20me3/H3K9me3 loss occur where KAP1 (which recruits HUSH/SETDB1)
# binds? Tests overlap of each differential-loss group (and the HUSH/CBX3 gene-
# family exons) with KAP1 peaks against a SIZE-MATCHED permutation null - because
# our domains span kb-Mb and KAP1 peaks are narrow, a naive overlap test would be
# confounded by domain length. Runs against both KAP1 tracks (Neuro-2a neural +
# ChIP-Atlas all-cell). Reads external mm39 KAP1 BEDs from data/annotation;
# skips gracefully if they are not present. No bigWig re-read.
# ================================================================
source("config.R")
suppressPackageStartupMessages({
  library(GenomicRanges); library(rtracklayer); library(ggplot2); library(data.table)
  library(TxDb.Mmusculus.UCSC.mm39.knownGene)
})

kap1_tracks <- c(Neural = kap1_neural_bed, allCell = kap1_allcell_bed)
kap1_tracks <- kap1_tracks[file.exists(kap1_tracks)]

if (length(kap1_tracks) == 0) {
  message("MINUTE_6: no KAP1 mm39 BEDs found in ", annot_dir,
          " - skipping (see README 'KAP1 intersection' to prepare them).")
} else {

  # --- Size-matched permutation: is `query` overlap higher than a random set of
  # regions drawn from `universe` with the SAME width distribution as `query`? ---
  size_perm <- function(obs, query_w, uni_w, uni_ov, nbins = 20, nperm = 2000, seed = 42) {
    brks <- unique(quantile(uni_w, seq(0, 1, length.out = nbins + 1)))
    ub   <- cut(uni_w,   brks, include.lowest = TRUE)
    qb   <- cut(query_w, brks, include.lowest = TRUE)
    need <- table(qb); bins <- split(seq_along(uni_w), ub)
    set.seed(seed); perm <- numeric(nperm)
    for (i in seq_len(nperm)) {
      smp <- unlist(lapply(names(need), function(b)
        if (!is.na(b) && need[[b]] > 0 && length(bins[[b]]) > 0)
          sample(bins[[b]], need[[b]], replace = TRUE)))
      perm[i] <- mean(uni_ov[smp])
    }
    data.frame(obs_frac = obs, exp_frac = mean(perm), fold = obs / mean(perm),
               z = (obs - mean(perm)) / sd(perm),
               perm_p = (sum(perm >= obs) + 1) / (nperm + 1))
  }

  ar  <- readRDS(annotated_rds)
  # (1) differential-loss groups; universe = all measured merged-set peaks
  k20 <- ar[["H4K20me3"]]; k9 <- ar[["H3K9me3"]]
  k20$k9 <- k9$log2FoldChange[match(k20$peak_id, k9$peak_id)]
  k20 <- k20[is.finite(k20$k9), ]
  k20$group <- with(k20, ifelse(log2FoldChange < -0.5 & k9 < -0.3,      "co_loss",
                        ifelse(log2FoldChange < -0.5 & abs(k9) < 0.15,  "H4K20me3_only",
                        ifelse(abs(log2FoldChange) < 0.15 & abs(k9) < 0.15, "stable", "other"))))
  uni_grp <- GRanges(k20$chr, IRanges(k20$start, k20$end)); mcols(uni_grp)$group <- k20$group
  wg <- width(uni_grp)

  # (2) HUSH/CBX3 gene-family exons; universe = all mm39 gene exons
  fam_gr <- if (file.exists(file.path(rds_dir, "family_exon_signal.rds")))
    readRDS(file.path(rds_dir, "family_exon_signal.rds"))$fam_gr else family_exons()
  seqlevelsStyle(fam_gr) <- "NCBI"; fam_gr <- keepStandardChromosomes(fam_gr, pruning.mode = "coarse")
  all_ex <- reduce(unlist(exonsBy(TxDb.Mmusculus.UCSC.mm39.knownGene, by = "gene"), use.names = FALSE))
  seqlevelsStyle(all_ex) <- "NCBI"; all_ex <- keepStandardChromosomes(all_ex, pruning.mode = "coarse")
  we <- width(all_ex)

  res <- do.call(rbind, lapply(names(kap1_tracks), function(tk) {
    kap    <- import(kap1_tracks[[tk]], format = "BED")
    ov_grp <- overlapsAny(uni_grp, kap, ignore.strand = TRUE)   # precompute once per track
    ov_ex  <- overlapsAny(all_ex,  kap, ignore.strand = TRUE)
    # differential-loss groups
    grp <- do.call(rbind, lapply(c("co_loss", "H4K20me3_only", "stable"), function(gp) {
      idx <- which(mcols(uni_grp)$group == gp)
      cbind(track = tk, set = gp, n = length(idx),
            size_perm(mean(ov_grp[idx]), wg[idx], wg, ov_grp))
    }))
    # gene families (query = family exons; universe = all exons)
    fam <- do.call(rbind, lapply(unique(as.character(mcols(fam_gr)$family)), function(fm) {
      q <- fam_gr[as.character(mcols(fam_gr)$family) == fm]
      cbind(track = tk, set = paste0("family:", fm), n = length(q),
            size_perm(mean(overlapsAny(q, kap, ignore.strand = TRUE)), width(q), we, ov_ex))
    }))
    rbind(grp, fam)
  }))
  fwrite(res, file.path(tables_dir, "kap1_intersection.tsv"), sep = "\t")
  cat("=== KAP1 intersection (observed vs size-matched null) ===\n"); print(res, row.names = FALSE)

  # Plot: observed vs expected KAP1 overlap fraction, fold + perm p annotated
  rd <- as.data.table(res)
  rd[, set := factor(set, levels = unique(set))]
  pd <- melt(rd[, .(track, set, observed = obs_frac, expected = exp_frac)],
             id.vars = c("track", "set"), variable.name = "type", value.name = "frac")
  g <- ggplot(pd, aes(set, frac, fill = type)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_text(data = rd, inherit.aes = FALSE,
              aes(set, obs_frac, label = sprintf("%.1fx\np=%.2g", fold, perm_p)),
              vjust = -0.3, size = 2.6, lineheight = 0.9) +
    facet_wrap(~track, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c(observed = "#b2182b", expected = "grey70"), name = NULL) +
    labs(title = "KAP1/TRIM28 overlap vs size-matched null (permutation, N=2000)",
         subtitle = "fold = observed / size-matched expected; groups: merged-peak universe, families: all-exon universe",
         x = NULL, y = "fraction overlapping a KAP1 peak") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))
  save_fig(g, "kap1_intersection", "differential_loss", width = 9, height = 7)
  message("Saved: kap1_intersection.tsv + differential_loss/kap1_intersection.{png,pdf}")

  # ================================================================
  # (B) KAP1 x GENOTYPE EFFECT — does KAP1 binding *predict* the CBX3-dependent
  #     loss? The overlap test above only asks "where is KAP1"; here we split
  #     regions/exons into KAP1-bound vs unbound and compare the KO log2FC of
  #     each mark (Wilcoxon). If bound regions lose MORE, KAP1 occupancy is
  #     coupled to the loss; if the distributions coincide, the loss is
  #     KAP1-independent. This is the direct integration of the static KAP1 ChIP
  #     with the genotype effect.
  # ================================================================
  eps <- 0.25
  # per-region mean log2FC for a mark from a signal matrix whose colnames are
  # <mark>_<genotype>_rep<k> (used for the family exons, which carry no DESeq2 fit)
  mark_l2fc <- function(mat, mk) {
    p <- strsplit(colnames(mat), "_"); mkc <- vapply(p, `[`, "", 1); gtc <- vapply(p, `[`, "", 2)
    wt <- rowMeans(mat[, mkc == mk & gtc == "WT",     drop = FALSE], na.rm = TRUE)
    ko <- rowMeans(mat[, mkc == mk & gtc == "HP1gKO", drop = FALSE], na.rm = TRUE)
    log2((ko + eps) / (wt + eps))
  }
  MARKS_B <- c("H4K20me3", "H3K9me3")

  # (B1) genome-wide merged peaks — use the DESeq2 log2FC already in the .rds
  gw <- do.call(rbind, lapply(names(kap1_tracks), function(tk) {
    kap <- import(kap1_tracks[[tk]], format = "BED")
    do.call(rbind, lapply(MARKS_B, function(mk) {
      g  <- ar[[mk]]
      bd <- overlapsAny(GRanges(g$chr, IRanges(g$start, g$end)), kap, ignore.strand = TRUE)
      ok <- is.finite(g$log2FoldChange)
      data.frame(track = tk, mark = mk, set = "all peaks",
                 bound = bd[ok], log2FC = g$log2FoldChange[ok])
    }))
  }))

  # (B2) HUSH/CBX3 family exons — log2FC from the bw_loci signal matrix
  famfc <- data.frame(family = as.character(mcols(fam_gr)$family),
                      H4K20me3 = mark_l2fc(readRDS(file.path(rds_dir, "family_exon_signal.rds"))$fam_mat, "H4K20me3"),
                      H3K9me3  = mark_l2fc(readRDS(file.path(rds_dir, "family_exon_signal.rds"))$fam_mat, "H3K9me3"))
  famB <- do.call(rbind, lapply(names(kap1_tracks), function(tk) {
    kap <- import(kap1_tracks[[tk]], format = "BED")
    b   <- overlapsAny(fam_gr, kap, ignore.strand = TRUE)
    do.call(rbind, lapply(MARKS_B, function(mk)
      data.frame(track = tk, mark = mk, set = famfc$family, bound = b, log2FC = famfc[[mk]])))
  }))

  dd <- as.data.table(rbind(gw, famB))
  dd[, set := factor(set, levels = c("all peaks", names(gene_families)))]
  # Wilcoxon bound vs unbound per (track, mark, set); require >=5 each side
  effB <- dd[, {
    nb <- sum(bound); nu <- sum(!bound)
    pv <- if (nb >= 5 && nu >= 5) suppressWarnings(wilcox.test(log2FC ~ bound))$p.value else NA_real_
    .(n_bound = nb, n_unbound = nu,
      med_bound   = if (nb > 0) median(log2FC[bound])  else NA_real_,
      med_unbound = if (nu > 0) median(log2FC[!bound]) else NA_real_,
      delta_median = (if (nb > 0) median(log2FC[bound]) else NA_real_) -
                     (if (nu > 0) median(log2FC[!bound]) else NA_real_),
      p = pv)
  }, by = .(track, mark, set)]
  effB[, padj := p.adjust(p, "BH")]
  fwrite(effB, file.path(tables_dir, "kap1_genotype_effect.tsv"), sep = "\t")
  cat("\n=== KAP1 x genotype effect (KO log2FC by KAP1-bound vs unbound) ===\n")
  print(effB[order(track, mark, set)], row.names = FALSE)

  # Plot: KO log2FC distribution split by KAP1-bound, mark x track facets
  lab <- effB[!is.na(delta_median),
              .(track, mark, set, y = 0.9,
                txt = sprintf("Δ%+.2f%s", delta_median,
                              ifelse(is.na(p), "", ifelse(p < 0.05, "*", ""))))]
  dd[, ngrp := .N, by = .(track, mark, set, bound)]   # drop tiny (n<5) boxes (single-point flat lines)
  gB <- ggplot(dd[ngrp >= 5], aes(set, log2FC, fill = bound)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
    geom_boxplot(outlier.shape = NA, position = position_dodge(0.8), width = 0.7, linewidth = 0.3) +
    geom_text(data = lab, inherit.aes = FALSE, aes(set, y, label = txt), size = 2.4) +
    facet_grid(mark ~ track) +
    scale_fill_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "grey75"),
                      labels = c(`TRUE` = "KAP1-bound", `FALSE` = "unbound"), name = NULL) +
    coord_cartesian(ylim = c(-2, 1)) +
    labs(title = "Does KAP1 binding predict the CBX3-dependent loss?",
         subtitle = "KO/WT log2FC split by KAP1 overlap; Δ = median(bound) − median(unbound), * Wilcoxon p<0.05",
         x = NULL, y = "log2 fold-change (HP1gKO / WT)") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))
  save_fig(gB, "kap1_genotype_effect", "differential_loss", width = 9, height = 7)
  message("Saved: kap1_genotype_effect.tsv + differential_loss/kap1_genotype_effect.{png,pdf}")
}
