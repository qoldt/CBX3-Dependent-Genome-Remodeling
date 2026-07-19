# ================================================================
# MINUTE_5 - Differential H3K9me3 vs H4K20me3 loss in the HP1gKO
# ----------------------------------------------------------------
# H3K9me3 and H4K20me3 share the 2 kb-merged peak set, so their DESeq2
# log2FoldChange align by peak_id. Split the shared regions into three groups
# by the two marks' change and characterise each by chromatin state (ChromHMM),
# repeats, and nearby genes - to resolve WHERE H4K20me3 is lost independently
# of H3K9me3.
#
# Operates purely on annotated_results columns (no bigWig re-read).
# Runs after MINUTE_1 (needs the annotated .rds).
#
# QC note: HP1gKO replicate 2 has globally elevated signal (a scaling residual,
# uniform across all marks) - with sizeFactors=1 this makes the measured losses
# CONSERVATIVE. Do NOT "correct" via DESeq2 median-of-ratios (erases the global
# H4K20me3 loss). See memory: gene-family-hush-derepression.
# ================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(data.table) })

theme_m <- theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank())

# --- Group thresholds (tunable) ---
H4_LOST_CUT <- -0.5    # H4K20me3 log2FC below this = "H4K20me3-lost"
K9_LOSS_CUT <- -0.3    # H3K9me3 log2FC below this = "H3K9me3-loss"
UNCH_CUT    <-  0.15   # |log2FC| below this = "unchanged"

ar  <- readRDS(annotated_rds)
k20 <- ar[["H4K20me3"]]
k9  <- ar[["H3K9me3"]]
k20$k9_l2fc <- k9$log2FoldChange[match(k20$peak_id, k9$peak_id)]
d <- k20[is.finite(k20$k9_l2fc) & is.finite(k20$log2FoldChange), ]
d$kb <- (d$end - d$start) / 1000

d$group <- with(d, ifelse(log2FoldChange < H4_LOST_CUT & k9_l2fc < K9_LOSS_CUT,      "co_loss",
                   ifelse(log2FoldChange < H4_LOST_CUT & abs(k9_l2fc) < UNCH_CUT,    "H4K20me3_only",
                   ifelse(abs(log2FoldChange) < UNCH_CUT & abs(k9_l2fc) < UNCH_CUT,  "stable", NA))))
g <- d[!is.na(d$group), ]
g$group <- factor(g$group, levels = c("stable", "H4K20me3_only", "co_loss"))
grp_cols <- c(stable = "grey65", H4K20me3_only = "#1f78b4", co_loss = "#e31a1c")

cat("=== group sizes ===\n"); print(table(g$group))

hmm_cols <- grep("^hmm_", names(g), value = TRUE)

# ---------------------------------------------------------------
# 1. Group overview: size, log2FC
# ---------------------------------------------------------------
overview <- g %>% group_by(group) %>%
  summarise(n = n(),
            median_kb = round(median(kb), 2),
            median_H4K20me3_l2fc = round(median(log2FoldChange), 2),
            median_H3K9me3_l2fc  = round(median(k9_l2fc), 2),
            median_baseMean = round(median(baseMean), 2), .groups = "drop")
fwrite(overview, file.path(tables_dir, "diffloss_group_overview.tsv"), sep = "\t")
cat("\n=== group overview ===\n"); print(as.data.frame(overview))

g_size <- ggplot(g, aes(group, kb, fill = group)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
  scale_y_log10() + scale_fill_manual(values = grp_cols, guide = "none") +
  labs(title = "Domain size per differential-loss group", x = NULL, y = "domain size (kb, log10)") +
  theme_m
save_fig(g_size, "diffloss_domain_size", "differential_loss", width = 6, height = 4)

# ---------------------------------------------------------------
# 2. ChromHMM: mean coverage + enrichment (per-region Wilcoxon AND size-weighted)
# ---------------------------------------------------------------
# mean coverage per group x state (heatmap)
hmm_cov <- g %>% group_by(group) %>%
  summarise(across(all_of(hmm_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
fwrite(hmm_cov, file.path(tables_dir, "diffloss_chromHMM_coverage.tsv"), sep = "\t")
cov_long <- melt(as.data.table(hmm_cov), id.vars = "group", variable.name = "state", value.name = "cov")
cov_long[, state := sub("^hmm_", "", state)]
g_cov <- ggplot(cov_long, aes(state, group, fill = cov)) +
  geom_tile(colour = "white") + scale_fill_viridis_c(name = "mean\ncoverage") +
  labs(title = "Mean ChromHMM state coverage per group", x = "ChromHMM state", y = NULL) +
  theme_m + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(g_cov, "diffloss_chromHMM_coverage", "differential_loss", width = 9, height = 3.2)

# enrichment: fg vs bg, per-region Wilcoxon + size-weighted permutation
enrich_chromhmm <- function(df, fg, contrast) {
  pr <- do.call(rbind, lapply(hmm_cols, function(cl) {
    a <- df[[cl]][fg]; b <- df[[cl]][!fg]
    data.frame(state = sub("^hmm_", "", cl),
               mean_fg = mean(a, na.rm = TRUE), mean_bg = mean(b, na.rm = TRUE),
               log2_ratio = log2((mean(a, na.rm = TRUE) + 1e-4) / (mean(b, na.rm = TRUE) + 1e-4)),
               p_value = tryCatch(suppressWarnings(wilcox.test(a, b)$p.value),
                                  error = function(e) NA_real_), stringsAsFactors = FALSE)
  }))
  sw <- chromHMM_size_weighted(df[, hmm_cols], df$kb, fg)
  m <- merge(pr, sw, by = "state")
  m$contrast <- contrast
  m$p_adj_BH <- p.adjust(m$p_value, "BH"); m$perm_p_adj <- p.adjust(m$perm_p, "BH")
  m
}
run_contrast <- function(fgname, bgname) {
  sub <- g[g$group %in% c(fgname, bgname), ]
  enrich_chromhmm(sub, sub$group == fgname, paste0(fgname, "_vs_", bgname))
}
hmm_enrich <- rbind(
  run_contrast("H4K20me3_only", "stable"),
  run_contrast("co_loss",       "stable"),
  run_contrast("H4K20me3_only", "co_loss"))
fwrite(hmm_enrich, file.path(tables_dir, "diffloss_chromHMM_enrichment.tsv"), sep = "\t")

hmm_dot <- function(lr, fdr, ttl) {
  hmm_enrich %>% filter(is.finite(.data[[lr]]), !is.na(.data[[fdr]])) %>%
    mutate(neglog10FDR = -log10(pmax(.data[[fdr]], .Machine$double.xmin)),
           state = factor(state, levels = sort(unique(state)))) %>%
    ggplot(aes(contrast, state, size = neglog10FDR, colour = .data[[lr]])) +
    geom_point() +
    scale_colour_gradient2(low = "steelblue3", mid = "grey90", high = "firebrick2",
                           midpoint = 0, name = expression(log[2]~"ratio")) +
    scale_size_continuous(name = expression(-log[10]("FDR"))) +
    labs(title = ttl, x = NULL, y = "ChromHMM state") +
    theme_m + theme(axis.text.x = element_text(angle = 20, hjust = 1))
}
save_fig(hmm_dot("log2_ratio", "p_adj_BH", "ChromHMM enrichment - PER-REGION"),
         "diffloss_chromHMM_dotplot", "differential_loss", width = 7, height = 6)
save_fig(hmm_dot("w_log2_ratio", "perm_p_adj", "ChromHMM enrichment - SIZE-WEIGHTED (by domain length)"),
         "diffloss_chromHMM_dotplot_sizeweighted", "differential_loss", width = 7, height = 6)

# ---------------------------------------------------------------
# 3. Repeats: composition + specific-family enrichment (Fisher, group vs stable)
# ---------------------------------------------------------------
rep_comp <- g %>% mutate(rc = ifelse(is.na(repeat_class), "NA", repeat_class)) %>%
  count(group, rc) %>% group_by(group) %>% mutate(frac = n / sum(n)) %>% ungroup()
fwrite(rep_comp, file.path(tables_dir, "diffloss_repeat_composition.tsv"), sep = "\t")
g_rep <- ggplot(rep_comp, aes(group, frac, fill = rc)) + geom_col() +
  labs(title = "Repeat-class composition per group", x = NULL, y = "fraction of regions", fill = "repeat_class") +
  theme_m
save_fig(g_rep, "diffloss_repeat_composition", "differential_loss", width = 6, height = 4)

fam <- ifelse(is.na(g$repeat_family), "", g$repeat_family)
nm  <- ifelse(is.na(g$repeat_name),   "", g$repeat_name)
rep_feats <- list(
  LINE = g$repeat_class == "LINE", SINE = g$repeat_class == "SINE",
  LTR = g$repeat_class == "LTR", none = g$repeat_class == "none",
  complex = g$repeat_class == "complex",
  `LTR:ERVK` = grepl("\\bERVK\\b", fam), `LTR:ERVL-MaLR` = grepl("\\bERVL-MaLR\\b", fam),
  `LTR:ERV1` = grepl("\\bERV1\\b", fam),
  `IAPEz-int` = grepl("IAPEz", nm), `IAPLTR` = grepl("IAPLTR", nm),
  # young/active mouse L1 (HUSH/TRIM28 substrates) - the mouse analogue of the
  # human L1HS/L1PA HUSH targets that can't be lifted over from hg38
  `LINE:L1MdA` = grepl("L1MdA", nm), `LINE:L1MdT` = grepl("L1MdT", nm),
  `LINE:L1MdGf` = grepl("L1MdGf", nm), `LINE:L1MdF` = grepl("L1MdF", nm),
  `LINE:young_L1(A/T/Gf)` = grepl("L1Md(A|T|Gf)", nm))
rep_enrich <- do.call(rbind, lapply(c("H4K20me3_only", "co_loss"), function(grp) {
  do.call(rbind, lapply(names(rep_feats), function(fn) {
    f <- rep_feats[[fn]]
    idx <- g$group %in% c(grp, "stable")           # restrict to grp vs stable

    tt  <- table(g$group[idx] == grp, f[idx])
    ftr <- suppressWarnings(fisher.test(tt))
    data.frame(group = grp, feature = fn,
               frac_group  = mean(f[g$group == grp]),
               frac_stable = mean(f[g$group == "stable"]),
               odds_ratio = unname(ftr$estimate), p_value = ftr$p.value, stringsAsFactors = FALSE)
  }))
}))
rep_enrich <- rep_enrich %>% group_by(group) %>% mutate(p_adj_BH = p.adjust(p_value, "BH")) %>% ungroup()
fwrite(rep_enrich, file.path(tables_dir, "diffloss_repeat_enrichment.tsv"), sep = "\t")
g_repe <- rep_enrich %>% filter(is.finite(odds_ratio), odds_ratio > 0) %>%
  ggplot(aes(log2(odds_ratio), feature, colour = group, size = -log10(pmax(p_adj_BH, 1e-300)))) +
  geom_vline(xintercept = 0, linetype = "dashed") + geom_point(alpha = 0.85) +
  scale_colour_manual(values = grp_cols[c("H4K20me3_only", "co_loss")]) +
  scale_size_continuous(name = expression(-log[10]("FDR"))) +
  labs(title = "Repeat-family enrichment vs stable regions",
       x = "log2 odds ratio (group vs stable)", y = NULL, colour = "group") + theme_m
save_fig(g_repe, "diffloss_repeat_enrichment", "differential_loss", width = 7, height = 4.5)

# ---------------------------------------------------------------
# 4. Nearby genes: region, distance-to-TSS, gene families, gene lists
# ---------------------------------------------------------------
reg_comp <- g %>% mutate(gr = ifelse(is.na(genomic_region), "NA", genomic_region)) %>%
  count(group, gr) %>% group_by(group) %>% mutate(frac = n / sum(n)) %>% ungroup()
fwrite(reg_comp, file.path(tables_dir, "diffloss_region_composition.tsv"), sep = "\t")
g_reg <- ggplot(reg_comp, aes(group, frac, fill = gr)) + geom_col() +
  labs(title = "Genomic-region composition per group", x = NULL, y = "fraction of regions", fill = "region") +
  theme_m
save_fig(g_reg, "diffloss_region_composition", "differential_loss", width = 6.5, height = 4)

g_tss <- ggplot(g, aes(group, abs(distance_to_tss) + 1, fill = group)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3) + scale_y_log10() +
  scale_fill_manual(values = grp_cols, guide = "none") +
  labs(title = "Distance to nearest TSS per group", x = NULL, y = "|distance to TSS| (bp, log10)") + theme_m
save_fig(g_tss, "diffloss_distance_to_tss", "differential_loss", width = 6, height = 4)

# gene-family overlap per group
sym <- ifelse(is.na(g$SYMBOL), "", as.character(g$SYMBOL))
fam_frac <- do.call(rbind, lapply(names(gene_families), function(fn) {
  flag <- grepl(gene_families[[fn]], sym)
  data.frame(family = fn,
             aggregate(flag, list(group = g$group), function(x) round(100 * mean(x), 2)) |>
               setNames(c("group", "pct")))
}))
fwrite(fam_frac, file.path(tables_dir, "diffloss_gene_family_pct.tsv"), sep = "\t")

# per-group gene lists (for downstream GO)
gene_summary <- data.frame(group = character(), n_genes = integer())
for (grp in levels(g$group)) {
  genes <- sort(unique(na.omit(g$SYMBOL[g$group == grp])))
  writeLines(genes, file.path(tables_dir, paste0("diffloss_genes_", grp, ".txt")))
  gene_summary <- rbind(gene_summary, data.frame(group = grp, n_genes = length(genes)))
}
cat("\n=== unique genes per group (exported for GO) ===\n"); print(gene_summary)
cat("\n=== gene-family % per group ===\n"); print(fam_frac)

message("MINUTE_5 done: diffloss_* tables + figures + per-group gene lists")
