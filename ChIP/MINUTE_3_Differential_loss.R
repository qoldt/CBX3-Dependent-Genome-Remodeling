# ================================================================
# MINUTE_3 - Differential H3K9me3 vs H4K20me3 loss in the HP1gKO
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
# QC note: HP1gKO replicate 2 has globally elevated signal, uniform across all
# marks. It is now EXCLUDED BY DEFAULT (config.R) - it is discordant with rep3,
# its own technical replicate, so it is a technical failure rather than a scaling
# residual. Losses were correspondingly CONSERVATIVE while it was included
# (H4K20me3 median log2FC -0.474 with it, -0.678 without).
# Do NOT "correct" via DESeq2 median-of-ratios (erases the global H4K20me3 loss).
#
# THRESHOLD CAVEAT: the cutoffs below are ABSOLUTE, so these groups are NOT
# shift-invariant - when a mark's global median moves, peaks migrate between
# groups with no change in biology. Excluding rep2 moved H3K9me3's median from
# -0.22 to -0.42 and restructured the compartments (stable 4,105 -> 931;
# H4K20me3_only 13,682 -> 8,481; co_loss 19,938 -> 46,921). With absolute cuts
# the question is "did it lose?" (against a genome-wide loss, nearly everything
# did); centring each mark on its own median instead asks "did it lose MORE than
# typical?", which is what a compartment claim actually needs.
# ================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(data.table) })
# geno_colors / theme_m / mark_levels come from config.R

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
grp_cols <- c(stable = "grey65", H4K20me3_only = gaby_cols[3], co_loss = gaby_cols[5])

cat("=== group sizes ===\n"); print(table(g$group))

# --- Absolute vs median-centred group definitions ---------------------------
# The cutoffs above are ABSOLUTE, so against a genome-wide loss they mostly ask
# "did it lose?" - to which nearly everything answers yes, and the answer moves
# whenever the global shift moves. Centring each mark on its own median instead
# asks "did it lose MORE than typical?", which is what a compartment claim
# needs and which is invariant to the shift.
#
# Measured across sample sets (all libraries vs HP1gKO rep2 excluded), group
# sizes move by a mean of 84% under absolute cutoffs but only 11% under
# median-centred ones - and the ordering flips: absolute cutoffs make co_loss
# dominate once rep2 is dropped, whereas median-centred keeps H4K20me3_only >
# co_loss in BOTH sample sets. The uncoupling result therefore survives, but
# only when stated relative to each mark's own median.
#
# This table reports BOTH schemes for the current run so the dependence is
# visible in the outputs, not just in the docs. The absolute scheme remains the
# one used downstream.
med_h4 <- median(d$log2FoldChange, na.rm = TRUE)
med_k9 <- median(d$k9_l2fc,        na.rm = TRUE)
grp_centred <- with(d, ifelse(log2FoldChange - med_h4 < H4_LOST_CUT & k9_l2fc - med_k9 < K9_LOSS_CUT, "co_loss",
                       ifelse(log2FoldChange - med_h4 < H4_LOST_CUT & abs(k9_l2fc - med_k9) < UNCH_CUT, "H4K20me3_only",
                       ifelse(abs(log2FoldChange - med_h4) < UNCH_CUT & abs(k9_l2fc - med_k9) < UNCH_CUT, "stable", NA))))
lv <- c("stable", "H4K20me3_only", "co_loss")
group_defs <- data.frame(
  group          = lv,
  n_absolute     = as.integer(table(factor(g$group,      levels = lv))),
  n_median_centred = as.integer(table(factor(grp_centred, levels = lv))),
  stringsAsFactors = FALSE)
group_defs$median_H4K20me3_used <- med_h4
group_defs$median_H3K9me3_used  <- med_k9
group_defs$n_peaks_considered   <- nrow(d)
fwrite(group_defs, file.path(tables_dir, "diffloss_group_definitions.tsv"), sep = "\t")
cat("\n=== group sizes: absolute vs median-centred cutoffs ===\n")
cat(sprintf("(marks centred on median H4K20me3 %+0.3f, H3K9me3 %+0.3f)\n", med_h4, med_k9))
print(group_defs[, c("group", "n_absolute", "n_median_centred")], row.names = FALSE)
cat("Saved: tables/diffloss_group_definitions.tsv\n")

# --- Figure: the two schemes side by side -----------------------------------
# The cutoff LINES are identical in both panels; what moves is the data. Left:
# raw log2FC, where the point cloud sits well below/left of the lines because
# both marks lost genome-wide - so the "H3K9me3 lost" line falls in the middle
# of the bulk and 63% of all peaks clear it. Right: each mark centred on its own
# median, so the same lines now cut the tails rather than the body.
scheme_df <- rbind(
  data.frame(scheme = "absolute (raw log2FC)",
             k9 = d$k9_l2fc, h4 = d$log2FoldChange,
             group = factor(g$group[match(d$peak_id, g$peak_id)], levels = lv)),
  data.frame(scheme = "median-centred (per mark)",
             k9 = d$k9_l2fc - med_k9, h4 = d$log2FoldChange - med_h4,
             group = factor(grp_centred, levels = lv)))
scheme_df$group <- addNA(scheme_df$group)
levels(scheme_df$group)[is.na(levels(scheme_df$group))] <- "unclassified"
set.seed(42)                       # subsample purely for render size
sub <- scheme_df[sample(nrow(scheme_df), min(60000, nrow(scheme_df))), ]

g_scheme <- ggplot(sub, aes(k9, h4, colour = group)) +
  geom_point(size = 0.25, alpha = 0.25) +
  geom_vline(xintercept = K9_LOSS_CUT, linewidth = 0.3, colour = "grey25") +
  geom_hline(yintercept = H4_LOST_CUT, linewidth = 0.3, colour = "grey25") +
  geom_vline(xintercept = c(-UNCH_CUT, UNCH_CUT), linetype = "dashed",
             linewidth = 0.3, colour = "grey25") +
  geom_hline(yintercept = c(-UNCH_CUT, UNCH_CUT), linetype = "dashed",
             linewidth = 0.3, colour = "grey25") +
  scale_colour_manual(values = c(grp_cols, unclassified = "grey85"), name = NULL,
                      guide = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  coord_cartesian(xlim = c(-2, 1.5), ylim = c(-2.5, 1.5)) +
  facet_wrap(~scheme) +
  labs(title = "Compartment cutoffs: absolute vs median-centred",
       subtitle = sprintf(paste("same cutoff lines in both panels - only the data moves.",
                                "medians: H4K20me3 %+0.2f, H3K9me3 %+0.2f"), med_h4, med_k9),
       x = "H3K9me3 log2FC", y = "H4K20me3 log2FC",
       caption = sprintf(paste0(
         "solid = 'lost' cutoffs (H4K20me3 < %s, H3K9me3 < %s); dashed = 'unchanged' band (+/-%s)\n",
         "LEFT: the bulk sits past the H3K9me3 line, so a typical peak is auto-labelled lost - ",
         "63%% of all peaks clear it, and 65%% of co_loss have H3K9me3 within +/-0.15 of the median.\n",
         "RIGHT: centring makes the cutoffs ask 'more than typical?' - only 12%% of peaks clear the same line."),
         H4_LOST_CUT, K9_LOSS_CUT, UNCH_CUT)) +
  theme_m + theme(plot.caption = element_text(hjust = 0, size = 8))
save_fig(g_scheme, "diffloss_cutoff_schemes", "differential_loss", width = 11, height = 5.5)

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
  geom_tile(colour = "white") + scale_fill_heat0(name = "mean\ncoverage") +
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
    scale_colour_heat0_div( name = expression(log[2]~"ratio")) +
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
  scale_fill_manual(values = repeat_palette, na.value = "grey80") +
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
  scale_fill_disc() +
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

message("MINUTE_3 done: diffloss_* tables + figures + per-group gene lists")
