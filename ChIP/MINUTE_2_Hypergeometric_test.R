# ================================================================
# MINUTE_2 - Hypergeometric enrichment of significant peaks vs annotations
# Runs standalone after MINUTE_1, or in sequence via run_MINUTE.R
# ================================================================
source("config.R")

# DESeq2 + annotation results produced by MINUTE_1
annotated_results <- readRDS(annotated_rds)

# Repeat + TAD annotation (shared loader from config.R)
ann     <- load_annotation()
line_gr <- ann$line
sine_gr <- ann$sine
ltr_gr  <- ann$ltr
tad_gr  <- ann$tad


# ================================
# Enrichment (Hypergeometric) per ChIP and Annotation
# ================================

combined_sig <- do.call(rbind, lapply(names(annotated_results), function(mark) {
  df <- annotated_results[[mark]]
  df$ChIP <- mark
  df$peak_id <- paste0(df$chr, ":", df$start, "-", df$end)
  df
}))



# Keep only significant peaks (thresholds defined in config.R)
combined_sig <- combined_sig[is_significant(combined_sig), ]


try({
  seqlevelsStyle(line_gr) <- "UCSC"
  seqlevelsStyle(sine_gr) <- "UCSC"
  seqlevelsStyle(ltr_gr)  <- "UCSC"
  seqlevelsStyle(tad_gr)  <- "UCSC"
}, silent = TRUE)

# Helper to compute one hypergeometric row + OR + direction
hyper_row <- function(chip, annotation, all_df, sig_df, flag_all) {
  # Universe counts
  N <- nrow(all_df)
  K <- sum(flag_all, na.rm = TRUE)
  n <- nrow(sig_df)
  
  # Map sig -> all via peak_id
  idx_sig_in_all <- match(sig_df$peak_id, all_df$peak_id)
  flag_sig <- flag_all[idx_sig_in_all]
  k <- sum(flag_sig, na.rm = TRUE)
  
  # Hypergeometric (one-sided for enrichment)
  if (is.na(N) || N == 0 || is.na(n) || is.na(K) || K < 0) {
    pval <- NA_real_
    expected <- NA_real_
    enrich <- NA_real_
  } else {
    expected <- if (N > 0) n * (K / N) else NA_real_
    pval <- if (n > 0 && K > 0 && N > K) phyper(k - 1, K, N - K, n, lower.tail = FALSE) else NA_real_
    enrich <- if (!is.na(expected) && expected > 0) k / expected else NA_real_
  }
  
  # 2x2 table for OR
  a <- k
  b <- n - k
  c <- K - k
  d <- N - K - b
  # continuity correction if any cell is zero or negative (safety)
  needs_cc <- any(c(a, b, c, d) <= 0)
  if (needs_cc) {
    a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d <- d + 0.5
  }
  or_val <- (a * d) / (b * c)
  if (!is.finite(or_val)) or_val <- NA_real_
  
  # Direction from median log2FC among significant peaks overlapping the annotation
  # (enriched if median > 0, depleted if median <= 0; NA if no overlapping significant peaks)
  if (!is.null(sig_df$log2FoldChange) && k > 0) {
    l2fc_med <- stats::median(sig_df$log2FoldChange[which(flag_sig)], na.rm = TRUE)
    direction <- if (is.finite(l2fc_med) && !is.na(l2fc_med) && l2fc_med > 0) "enriched" else "depleted"
  } else {
    l2fc_med <- NA_real_
    direction <- NA_character_
  }
  
  data.frame(
    ChIP = chip,
    annotation = annotation,
    N_total = N,
    K_in_annotation = K,
    n_signif = n,
    k_signif_in_annotation = k,
    expected_overlap = expected,
    enrichment = enrich,
    odds_ratio = or_val,
    l2fc_median_sig_in_annotation = l2fc_med,
    direction = direction,
    p_value = pval,
    stringsAsFactors = FALSE
  )
}

# Build per-ChIP enrichment with direction
enrichment_list <- lapply(names(annotated_results), function(mark) {
  all_df <- annotated_results[[mark]]
  all_df$ChIP <- mark
  all_df$peak_id <- paste0(all_df$chr, ":", all_df$start, "-", all_df$end)
  
  peak_gr <- GRanges(seqnames = all_df$chr, ranges = IRanges(all_df$start, all_df$end))
  seqlevelsStyle(peak_gr) <- "UCSC"
  
  # Flags over ALL measured peaks
  flag_TAD <- isTRUE(all_df$overlaps_with_Tad_boundary) | (all_df$overlaps_with_Tad_boundary %in% TRUE)
  
  flag_LINE <- rep(FALSE, length(peak_gr))
  flag_SINE <- rep(FALSE, length(peak_gr))
  flag_LTR  <- rep(FALSE, length(peak_gr))
  if (exists("line_gr")) flag_LINE[queryHits(findOverlaps(peak_gr, line_gr, ignore.strand = TRUE))] <- TRUE
  if (exists("sine_gr")) flag_SINE[queryHits(findOverlaps(peak_gr, sine_gr, ignore.strand = TRUE))] <- TRUE
  if (exists("ltr_gr"))  flag_LTR[ queryHits(findOverlaps(peak_gr, ltr_gr,  ignore.strand = TRUE))] <- TRUE
  
  fam <- if ("repeat_family" %in% names(all_df)) all_df$repeat_family else rep("", nrow(all_df))
  fam <- ifelse(is.na(fam), "", fam)
  flag_ERVK      <- grepl("\\bERVK\\b", fam)
  flag_ERVL_MaLR <- grepl("\\bERVL-MaLR\\b", fam)
  flag_ERV1      <- grepl("\\bERV1\\b", fam)
  
  nm <- if ("repeat_name" %in% names(all_df)) all_df$repeat_name else rep("", nrow(all_df))
  nm <- ifelse(is.na(nm), "", nm)
  flag_IAPEz_int <- grepl("\\bIAPEz-int\\b", nm)
  flag_IAPLTR1a  <- grepl("\\bIAPLTR1a_Mm\\b", nm)
  
  # Significant set for this ChIP
  sig_df <- combined_sig %>%
    dplyr::filter(ChIP == mark) %>%
    dplyr::mutate(peak_id = paste0(chr, ":", start, "-", end))
  
  rows <- list(
    hyper_row(mark, "TAD_boundary",        all_df, sig_df, flag_TAD),
    hyper_row(mark, "LINE",                all_df, sig_df, flag_LINE),
    hyper_row(mark, "SINE",                all_df, sig_df, flag_SINE),
    hyper_row(mark, "LTR",                 all_df, sig_df, flag_LTR),
    hyper_row(mark, "LTR:ERVK",            all_df, sig_df, flag_ERVK),
    hyper_row(mark, "LTR:ERVL-MaLR",       all_df, sig_df, flag_ERVL_MaLR),
    hyper_row(mark, "LTR:ERV1",            all_df, sig_df, flag_ERV1),
    hyper_row(mark, "LTR:ERVK:IAPEz-int",  all_df, sig_df, flag_IAPEz_int),
    hyper_row(mark, "LTR:ERVK:IAPLTR1a_Mm",   all_df, sig_df, flag_IAPLTR1a)
  )
  do.call(rbind, rows)
})

enrichment_df <- do.call(rbind, enrichment_list) %>%
  dplyr::arrange(ChIP, annotation) %>%
  dplyr::group_by(ChIP) %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup()

# Persist the enrichment table
fwrite(enrichment_df, file.path(tables_dir, "enrichment_hypergeometric.tsv"), sep = "\t")

# Inspect
enrichment_df

library(dplyr)
library(ggplot2)

# Prepare plotting data --------------------------------------------
plot_df <- enrichment_df %>%
  # keep rows we can actually plot
  filter(!is.na(odds_ratio), is.finite(odds_ratio),
         !is.na(p_adj_BH), !is.na(direction)) %>%
  mutate(
    neglog10FDR = -log10(pmax(p_adj_BH, .Machine$double.xmin)),  # avoid Inf
    direction = factor(direction, levels = c("depleted", "enriched")),
    # order annotations within each ChIP by OR (optional but helpful)
    annotation = as.factor(annotation)
  )

# If you’d like a consistent annotation order across facets, uncomment:
# plot_df$annotation <- factor(plot_df$annotation,
#                              levels = unique(plot_df$annotation))

# Plot --------------------------------------------------------------
p_enrich <- ggplot(plot_df,
                   aes(x = annotation,
                       y = odds_ratio,
                       size = neglog10FDR,
                       color = direction)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(alpha = 0.9) +
  facet_wrap(~ ChIP, nrow = 1) +
  scale_y_log10() +
  scale_color_manual(values = c(depleted = "steelblue3",
                                enriched = "firebrick2")) +
  labs(x = "Annotation",
       y = "Odds ratio (log10 scale)",
       size = expression(-log[10]("FDR")),
       color = NULL) +
  theme_bw(base_size = 12) +
  coord_flip()+
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )+
  theme_minimal()

ggsave(file.path(fig_dir, "enrichment_dotplot.pdf"), p_enrich, width = 12, height = 5)
p_enrich


# ================================================================
# ChromHMM chromatin-state enrichment (coverage-based)
# ----------------------------------------------------------------
# ChromHMM is annotated per peak in MINUTE_1 as per-state COVERAGE FRACTIONS
# (hmm_<state> columns) rather than a single size-biased label, because the
# H4K20me3/H3K9me3 domains span many states. For each mark x state we compare
# the mean coverage fraction in the significant set vs the background
# (non-significant peaks): log2 ratio of means = effect size, Wilcoxon = test.
# ================================================================

hmm_cols   <- grep("^hmm_", names(annotated_results[[1]]), value = TRUE)
hmm_states <- sub("^hmm_", "", hmm_cols)

chromhmm_enrich <- do.call(rbind, lapply(names(annotated_results), function(mark) {
  d   <- annotated_results[[mark]]
  d$ChIP <- mark
  sig <- is_significant(d)
  # per-region (unweighted) Wilcoxon
  pr <- do.call(rbind, lapply(hmm_states, function(st) {
    col <- paste0("hmm_", st)
    a <- d[[col]][sig]; b <- d[[col]][!sig]           # background = non-significant
    data.frame(ChIP = mark, state = st, n_sig = sum(sig),
               mean_cov_sig = mean(a, na.rm = TRUE), mean_cov_bg = mean(b, na.rm = TRUE),
               log2_ratio = log2((mean(a, na.rm = TRUE) + 1e-4) / (mean(b, na.rm = TRUE) + 1e-4)),
               p_value = tryCatch(suppressWarnings(wilcox.test(a, b)$p.value),
                                  error = function(e) NA_real_),
               stringsAsFactors = FALSE)
  }))
  # size-weighted (by domain length) with permutation test
  sw <- chromHMM_size_weighted(d[, hmm_cols], d$end - d$start, sig)
  merge(pr, sw, by = "state")
}))

chromhmm_enrich <- chromhmm_enrich %>%
  dplyr::group_by(ChIP) %>%
  dplyr::mutate(p_adj_BH   = p.adjust(p_value, method = "BH"),
                perm_p_adj = p.adjust(perm_p,  method = "BH")) %>%
  dplyr::ungroup()

fwrite(chromhmm_enrich, file.path(tables_dir, "enrichment_chromHMM.tsv"), sep = "\t")

# Two dot plots: PER-REGION (Wilcoxon) and SIZE-WEIGHTED (permutation) - report both
hmm_dotplot <- function(df, lr, fdr, ttl, sub) {
  df %>% dplyr::filter(is.finite(.data[[lr]]), !is.na(.data[[fdr]])) %>%
    dplyr::mutate(neglog10FDR = -log10(pmax(.data[[fdr]], .Machine$double.xmin)),
                  state = factor(state, levels = sort(unique(state)))) %>%
    ggplot(aes(x = ChIP, y = state, size = neglog10FDR, colour = .data[[lr]])) +
    geom_point() +
    scale_colour_gradient2(low = "steelblue3", mid = "grey90", high = "firebrick2",
                           midpoint = 0, name = expression(log[2]~"ratio")) +
    scale_size_continuous(name = expression(-log[10]("FDR"))) +
    labs(title = ttl, subtitle = sub, x = NULL, y = "ChromHMM state") +
    theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank())
}
ggsave(file.path(fig_dir, "enrichment_chromHMM_dotplot.pdf"),
       hmm_dotplot(chromhmm_enrich, "log2_ratio", "p_adj_BH",
                   "ChromHMM enrichment of significant peaks - PER-REGION",
                   "each region equal; colour = log2 mean-coverage ratio (sig vs bg); size = FDR"),
       width = 8, height = 7)
ggsave(file.path(fig_dir, "enrichment_chromHMM_sizeweighted_dotplot.pdf"),
       hmm_dotplot(chromhmm_enrich, "w_log2_ratio", "perm_p_adj",
                   "ChromHMM enrichment of significant peaks - SIZE-WEIGHTED",
                   "weighted by domain length; colour = log2 territory ratio (sig vs bg); size = perm FDR"),
       width = 8, height = 7)
message("Saved: enrichment_chromHMM.tsv + per-region & size-weighted dot plots")


# ================================================================
# ChromHMM: H3K9me3-loss vs H3K9me3-unchanged among H4K20me3-lost regions
# ----------------------------------------------------------------
# What chromatin STATE distinguishes regions where H3K9me3 co-falls with
# H4K20me3 from those losing H4K20me3 only? Restrict to H4K20me3-lost regions
# (shared 2 kb-merged peak set), split on the H3K9me3 change, then compare the
# per-state coverage fractions (mean + Wilcoxon). Thresholds are tunable.
# ================================================================
h4_lost_cut <- -0.5    # H4K20me3 log2FC below this = "H4K20me3-lost" population
k9_loss_cut <- -0.3    # H3K9me3 log2FC below this = "H3K9me3-loss"
k9_unch_cut <-  0.15   # |H3K9me3 log2FC| below this = "H3K9me3-unchanged"

k20 <- annotated_results[["H4K20me3"]]
k9  <- annotated_results[["H3K9me3"]]
k20$k9_l2fc <- k9$log2FoldChange[match(k20$peak_id, k9$peak_id)]
pop <- k20[is.finite(k20$k9_l2fc) & k20$log2FoldChange < h4_lost_cut, ]
pop$k9group <- ifelse(pop$k9_l2fc < k9_loss_cut, "H3K9me3_loss",
                ifelse(abs(pop$k9_l2fc) < k9_unch_cut, "H3K9me3_unchanged", NA))
pop <- pop[!is.na(pop$k9group), ]
message(sprintf("H4K20me3-lost: %d | H3K9me3-loss: %d | H3K9me3-unchanged: %d",
                nrow(pop), sum(pop$k9group == "H3K9me3_loss"),
                sum(pop$k9group == "H3K9me3_unchanged")))

hmm_cols <- grep("^hmm_", names(pop), value = TRUE)
k9split <- do.call(rbind, lapply(hmm_cols, function(col) {
  a <- pop[[col]][pop$k9group == "H3K9me3_loss"]
  b <- pop[[col]][pop$k9group == "H3K9me3_unchanged"]
  data.frame(state = sub("^hmm_", "", col),
             mean_cov_loss = mean(a, na.rm = TRUE),
             mean_cov_unchanged = mean(b, na.rm = TRUE),
             log2_ratio = log2((mean(a, na.rm = TRUE) + 1e-4) / (mean(b, na.rm = TRUE) + 1e-4)),
             p_value = tryCatch(suppressWarnings(wilcox.test(a, b)$p.value),
                                error = function(e) NA_real_),
             stringsAsFactors = FALSE)
}))
# size-weighted (by domain length); foreground = H3K9me3-loss, background = unchanged
sw9 <- chromHMM_size_weighted(pop[, hmm_cols], pop$end - pop$start,
                              pop$k9group == "H3K9me3_loss")
k9split <- merge(k9split, sw9, by = "state")
k9split$p_adj_BH   <- p.adjust(k9split$p_value, method = "BH")
k9split$perm_p_adj <- p.adjust(k9split$perm_p,  method = "BH")
fwrite(k9split, file.path(tables_dir, "enrichment_chromHMM_H3K9me3_loss_vs_unchanged.tsv"), sep = "\t")

k9_lollipop <- function(lr, fdr, ttl) {
  k9split %>%
    dplyr::filter(is.finite(.data[[lr]]), !is.na(.data[[fdr]])) %>%
    dplyr::mutate(neglog10FDR = -log10(pmax(.data[[fdr]], .Machine$double.xmin)),
                  state = reorder(state, .data[[lr]])) %>%
    ggplot(aes(.data[[lr]], state, size = neglog10FDR, colour = .data[[lr]])) +
    geom_vline(xintercept = 0, linetype = "dashed") + geom_point() +
    scale_colour_gradient2(low = "steelblue3", mid = "grey85", high = "firebrick2",
                           midpoint = 0, name = expression(log[2]~"ratio")) +
    scale_size_continuous(name = expression(-log[10]("FDR"))) +
    labs(title = ttl,
         subtitle = "red/right = higher where H3K9me3 also falls; blue/left = higher where H3K9me3 stays",
         x = "log2 coverage ratio (H3K9me3-loss / unchanged)", y = "ChromHMM state") +
    theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank())
}
ggsave(file.path(fig_dir, "enrichment_chromHMM_H3K9me3_loss_vs_unchanged.pdf"),
       k9_lollipop("log2_ratio", "p_adj_BH",
                   "ChromHMM: H3K9me3-loss vs unchanged (H4K20me3-lost) - PER-REGION"),
       width = 8, height = 7)
ggsave(file.path(fig_dir, "enrichment_chromHMM_H3K9me3_loss_vs_unchanged_sizeweighted.pdf"),
       k9_lollipop("w_log2_ratio", "perm_p_adj",
                   "ChromHMM: H3K9me3-loss vs unchanged (H4K20me3-lost) - SIZE-WEIGHTED"),
       width = 8, height = 7)
message("Saved: enrichment_chromHMM_H3K9me3_loss_vs_unchanged.{tsv, per-region pdf, size-weighted pdf}")


