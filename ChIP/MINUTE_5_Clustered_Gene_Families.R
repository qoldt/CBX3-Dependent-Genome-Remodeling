# ================================================================
# MINUTE_5 - Clustered gene families silenced by HUSH / CBX3
# ----------------------------------------------------------------
# The HUSH/CBX3-silenced families are genomic gene CLUSTERS - clustered
# protocadherins (chr18 Pcdha/b/g), KRAB-ZFP clusters, and the vomeronasal /
# olfactory receptor arrays. Quantify the input-scaled signal over each gene's
# silencing-relevant exon (bw_loci; KRAB-ZFP last exon, protocadherin first exon,
# Vmn/Olfr largest exon), test the KO change per mark, and ask whether H3K9me3
# co-loss (with H4K20me3) exceeds a matched genome-wide background.
#
# Persists results/rds/family_exon_signal.rds (consumed by MINUTE_6_Intersection)
# and per-family exon BEDs (deepTools inputs). Runs after MINUTE_1.
# geno_colors / theme_m / mark_levels come from config.R.
# ================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(data.table) })

fam_gr <- family_exons()                               # one relevant exon per gene
seqlevelsStyle(fam_gr) <- "NCBI"                       # match the bigWigs (as in MINUTE_1)
fam_gr <- keepStandardChromosomes(fam_gr, pruning.mode = "coarse")
cat(sprintf("\nGene-family target exons: %d genes across %d families\n",
            length(fam_gr), length(unique(mcols(fam_gr)$family))))
print(table(mcols(fam_gr)$family))

# quantify every sample over the family exons (bw_loci; parallel over samples)
fam_signal <- parallel::mclapply(seq_len(nrow(samples)), function(i) {
  bw <- samples$bigwig_path[i]
  if (!file.exists(bw)) return(rep(NA_real_, length(fam_gr)))
  tryCatch(mcols(bw_loci(bw, fam_gr))[[1]],
           error = function(e) rep(NA_real_, length(fam_gr)))
}, mc.cores = parallel::detectCores())
fam_mat <- do.call(cbind, fam_signal)
colnames(fam_mat) <- samples$sample_id
stopifnot(nrow(fam_mat) == length(fam_gr))
saveRDS(list(fam_gr = fam_gr, fam_mat = fam_mat),
        file.path(rds_dir, "family_exon_signal.rds"))

# Per-family BED of the target exons (deepTools inputs; NCBI seqnames matching
# the scaled bigWigs; strand-aware, name = gene symbol)
for (fam in unique(as.character(mcols(fam_gr)$family))) {
  gi  <- fam_gr[as.character(mcols(fam_gr)$family) == fam]
  bed <- data.frame(chrom = as.character(seqnames(gi)), start = start(gi) - 1L,
                    end = end(gi), name = as.character(mcols(gi)$gene),
                    score = 0L, strand = as.character(strand(gi)))
  bed <- bed[order(bed$chrom, bed$start), ]
  fwrite(bed, file.path(bed_dir, sprintf("family_%s_exons.bed",
         gsub("[^A-Za-z0-9]+", "_", fam))), sep = "\t", col.names = FALSE)
}
message("Saved per-family target-exon BEDs to ", bed_dir)

# per exon x mark: mean WT vs mean KO signal + log2 change
eps_f  <- 0.25
family <- as.character(mcols(fam_gr)$family)
gene   <- as.character(mcols(fam_gr)$gene)
long <- do.call(rbind, lapply(intersect(mark_levels, as.character(samples$mark)), function(mk) {
  wtc <- samples$sample_id[as.character(samples$mark) == mk & as.character(samples$genotype) == "WT"]
  koc <- samples$sample_id[as.character(samples$mark) == mk & as.character(samples$genotype) == "HP1gKO"]
  wt  <- rowMeans(fam_mat[, wtc, drop = FALSE], na.rm = TRUE)
  ko  <- rowMeans(fam_mat[, koc, drop = FALSE], na.rm = TRUE)
  data.frame(family = family, gene = gene,
             chr = as.character(seqnames(fam_gr)), start = start(fam_gr), end = end(fam_gr),
             mark = mk, wt_signal = wt, ko_signal = ko,
             log2FC = log2((ko + eps_f) / (wt + eps_f)), stringsAsFactors = FALSE)
}))
fwrite(long, file.path(tables_dir, "family_exon_signal.tsv"), sep = "\t")

# per family x mark summary: median change + paired Wilcoxon (KO vs WT over exons)
fam_summary <- long %>% group_by(family, mark) %>%
  summarise(n_exons = n(),
            median_WT = round(median(wt_signal, na.rm = TRUE), 2),
            median_KO = round(median(ko_signal, na.rm = TRUE), 2),
            median_log2FC = round(median(log2FC, na.rm = TRUE), 3),
            wilcox_p = tryCatch(suppressWarnings(
              wilcox.test(ko_signal, wt_signal, paired = TRUE)$p.value),
              error = function(e) NA_real_),
            .groups = "drop") %>%
  group_by(mark) %>% mutate(p_adj_BH = p.adjust(wilcox_p, "BH")) %>% ungroup()
fwrite(fam_summary, file.path(tables_dir, "family_exon_summary.tsv"), sep = "\t")
cat("\n=== Gene-family exon signal: median log2FC (KO/WT) per family x mark ===\n")
print(as.data.frame(fam_summary))

long$mark <- factor(long$mark, levels = mark_levels)
long$family <- factor(long$family, levels = names(gene_families))

# Plot A: per-exon change (log2 KO/WT) by mark, faceted by family
g_fam_chg <- ggplot(long, aes(mark, log2FC, fill = mark)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_boxplot(outlier.size = 0.25, outlier.alpha = 0.2, linewidth = 0.3) +
  facet_wrap(~family, nrow = 1) +
  scale_fill_viridis_d(option = "D", end = 0.9, guide = "none") +
  labs(title = "ChIP-signal change over gene-family exons (HP1gKO vs WT)",
       subtitle = "H3K9me3/H4K20me3 below 0 = loss of silencing marks over the family (CBX3/HUSH-dependent)",
       x = "Mark", y = "log2(KO / WT) per exon") +
  theme_m + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
save_fig(g_fam_chg, "family_exon_change", "gene_families", width = 14, height = 4)

# Plot B: WT vs KO signal DISTRIBUTION per gene, per family x mark (boxplot, so
# the across-gene spread is visible rather than a single median bar)
lvl <- melt(as.data.table(long[, c("family", "mark", "gene", "wt_signal", "ko_signal")]),
            id.vars = c("family", "mark", "gene"), variable.name = "Genotype", value.name = "signal")
lvl[, Genotype := factor(ifelse(Genotype == "wt_signal", "WT", "HP1gKO"), levels = c("WT", "HP1gKO"))]
g_fam_lvl <- ggplot(lvl[signal > 0], aes(mark, signal, fill = Genotype)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.15, linewidth = 0.3,
               position = position_dodge(0.8)) +
  facet_wrap(~family, nrow = 1) +          # log y is comparable across families
  scale_y_log10() +                        # signal spans orders of magnitude across marks
  scale_fill_manual(values = geno_colors) +
  labs(title = "ChIP signal over gene-family exons, per gene (WT vs HP1gKO)",
       subtitle = "boxes = across-gene distribution of per-gene mean signal (log10 y)",
       x = "Mark", y = "signal (log10)") +
  theme_m + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
save_fig(g_fam_lvl, "family_exon_signal_level", "gene_families", width = 14, height = 4)

# Plot C: summary dotplot (family x mark: colour = median log2FC, size = FDR)
g_fam_dot <- fam_summary %>%
  mutate(mark = factor(mark, levels = mark_levels),
         family = factor(family, levels = names(gene_families)),
         neglog10FDR = -log10(pmax(p_adj_BH, .Machine$double.xmin))) %>%
  ggplot(aes(mark, family, size = neglog10FDR, colour = median_log2FC)) +
  geom_point() +
  scale_colour_gradient2(low = "steelblue3", mid = "grey90", high = "firebrick2",
                         midpoint = 0, name = "median\nlog2FC") +
  scale_size_continuous(name = expression(-log[10]("FDR"))) +
  labs(title = "Gene-family exon change summary (HP1gKO vs WT)",
       subtitle = "colour = median log2(KO/WT); size = paired-Wilcoxon FDR", x = "Mark", y = NULL) +
  theme_m + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(g_fam_dot, "family_exon_dotplot", "gene_families", width = 7, height = 3.5)

# Plot D: MA plot (change vs abundance) per mark, coloured by family
g_fam_ma <- ggplot(long, aes((wt_signal + ko_signal) / 2, log2FC, colour = family)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_point(size = 0.4, alpha = 0.4) +
  scale_x_log10() +
  scale_colour_viridis_d(option = "D", end = 0.9, name = "Family") +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  facet_wrap(~mark, nrow = 1) +
  labs(title = "MA plot: gene-family exon change vs abundance",
       subtitle = "loss holds across the abundance range = not a signal-level artifact",
       x = "mean signal (WT+KO)/2 (log10)", y = "log2(KO / WT)") +
  theme_m
save_fig(g_fam_ma, "family_exon_MA", "gene_families", width = 14, height = 3.6)

# Plot E: per-gene REPLICATE heatmaps (rows = genes, cols = all 55 samples).
# Row z-scored raw signal; column-split by mark, top annotation = genotype;
# row-split by family. Shows per-gene, per-replicate WT->KO structure (not a mean).
suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize) })
col_meta <- data.frame(
  Sample   = colnames(fam_mat),
  Mark     = factor(as.character(samples$mark[match(colnames(fam_mat), samples$sample_id)]),
                    levels = mark_levels),
  Genotype = factor(as.character(samples$genotype[match(colnames(fam_mat), samples$sample_id)]),
                    levels = c("WT", "HP1gKO")))
row_fam <- factor(family, levels = names(gene_families))

draw_fam_hm <- function(mat, cols, ttl, name, w) {
  m <- t(scale(t(mat[, cols, drop = FALSE]))); m[!is.finite(m)] <- 0
  ha <- HeatmapAnnotation(Genotype = col_meta$Genotype[cols],
                          col = list(Genotype = geno_colors), show_annotation_name = TRUE)
  ht <- Heatmap(m, name = "row z", top_annotation = ha,
                column_split = droplevels(col_meta$Mark[cols]),
                row_split = row_fam, cluster_columns = FALSE, cluster_rows = TRUE,
                show_row_names = FALSE, show_column_names = FALSE,
                row_title_gp = gpar(fontsize = 9), column_title_gp = gpar(fontsize = 10),
                heatmap_legend_param = list(title = "row z-score"))
  save_base_fig(function() draw(ht, column_title = ttl), name, "gene_families",
                width = w, height = 11, dpi = 200)
}
draw_fam_hm(fam_mat, seq_len(ncol(fam_mat)),
            "Gene-family exon signal per replicate (all marks)",
            "family_exon_heatmap_allmarks", 11)
sil <- which(col_meta$Mark %in% c("H3K9me3", "H4K20me3"))
draw_fam_hm(fam_mat, sil,
            "Gene-family exon signal per replicate (silencing marks)",
            "family_exon_heatmap_silencing", 7)

message("Saved: family_exon_signal.{tsv,rds} + family_exon_summary.tsv")

# --- Family co-loss test: does H3K9me3 co-loss (with H4K20me3) differ by family? -
# Classify each gene by its H3K9me3 & H4K20me3 change (signal-based log2 KO/WT -
# noisier than DESeq2, but the family exons aren't all in the DESeq peak sets),
# then test whether the co-loss propensity differs across families.
loss_cut <- -0.3
cw <- dcast(as.data.table(long), family + gene ~ mark, value.var = "log2FC")
cw <- cw[is.finite(H3K9me3) & is.finite(H4K20me3)]
cw[, k20_lost := H4K20me3 < loss_cut]
cw[, k9_lost  := H3K9me3  < loss_cut]
cw[, class := fifelse(k20_lost & k9_lost,  "co-loss (both)",
              fifelse(k20_lost & !k9_lost, "H4K20me3-only",
              fifelse(!k20_lost & k9_lost, "H3K9me3-only", "neither")))]
cw[, family := factor(family, levels = names(gene_families))]
fwrite(cw, file.path(tables_dir, "family_coloss_classification.tsv"), sep = "\t")

# --- Matched genome-wide background: one exon from each of ~2000 random genes
#     (NOT overlapping any family locus), quantified the SAME way (bw_loci, same
#     eps, same loss_cut). This gives the GENOME-WIDE co-loss rate so each family
#     can be tested against BACKGROUND ("does this family co-lose more than the
#     genome at large?"), not only against the other families. -----------------
suppressPackageStartupMessages(library(TxDb.Mmusculus.UCSC.mm39.knownGene))
ex_by_gene <- exonsBy(TxDb.Mmusculus.UCSC.mm39.knownGene, by = "gene")
set.seed(1)
bg_gr <- unlist(endoapply(ex_by_gene[sample(names(ex_by_gene), min(2500, length(ex_by_gene)))],
                          function(e) e[which.max(width(e))]))            # largest exon, one/gene
seqlevelsStyle(bg_gr) <- "NCBI"
bg_gr <- keepStandardChromosomes(bg_gr, pruning.mode = "coarse")
bg_gr <- bg_gr[!overlapsAny(bg_gr, fam_gr, ignore.strand = TRUE)]         # exclude family loci
bg_gr <- head(bg_gr, 2000)
bg_signal <- parallel::mclapply(seq_len(nrow(samples)), function(i) {
  bw <- samples$bigwig_path[i]
  if (!file.exists(bw)) return(rep(NA_real_, length(bg_gr)))
  tryCatch(mcols(bw_loci(bw, bg_gr))[[1]], error = function(e) rep(NA_real_, length(bg_gr)))
}, mc.cores = parallel::detectCores())
bg_mat <- do.call(cbind, bg_signal); colnames(bg_mat) <- samples$sample_id
bg_l2 <- function(mk) {
  wtc <- samples$sample_id[as.character(samples$mark) == mk & as.character(samples$genotype) == "WT"]
  koc <- samples$sample_id[as.character(samples$mark) == mk & as.character(samples$genotype) == "HP1gKO"]
  log2((rowMeans(bg_mat[, koc, drop = FALSE], na.rm = TRUE) + eps_f) /
       (rowMeans(bg_mat[, wtc, drop = FALSE], na.rm = TRUE) + eps_f))
}
bg <- data.table(H3K9me3 = bg_l2("H3K9me3"), H4K20me3 = bg_l2("H4K20me3"))
bg <- bg[is.finite(H3K9me3) & is.finite(H4K20me3)]
bg[, k20_lost := H4K20me3 < loss_cut][, k9_lost := H3K9me3 < loss_cut]
bg_c <- bg[k20_lost == TRUE, sum(k9_lost)]; bg_n <- bg[k20_lost == TRUE, .N]
bg_rate <- bg_c / bg_n                                                    # background co-loss rate

# co-loss propensity = among H4K20me3-lost genes, fraction that ALSO lose H3K9me3
prop <- cw[k20_lost == TRUE, {
  k <- sum(k9_lost); bt <- binom.test(k, .N)
  .(n_H4K20me3_lost = .N, co_loss_frac = k / .N,
    lo = bt$conf.int[1], hi = bt$conf.int[2])
}, by = family][order(-co_loss_frac)]
# PRIMARY: each family vs genome-wide BACKGROUND — is co-loss (among H4K20me3-lost
# genes) more frequent in the family than in the matched random-gene background?
# OR = (family co-loss odds) / (background co-loss odds), oriented so OR>1 = family
# co-loses MORE than background. Fisher p from the 2x2; BH across families.
vsbg <- do.call(rbind, lapply(names(gene_families), function(fm) {
  a <- cw[k20_lost == TRUE & family == fm & k9_lost == TRUE,  .N]   # family: co-loss
  b <- cw[k20_lost == TRUE & family == fm & k9_lost == FALSE, .N]   # family: H4K20me3-only
  if (a + b == 0) return(NULL)
  or <- (a / max(b, 0.5)) / (max(bg_c, 0.5) / max(bg_n - bg_c, 0.5))
  p  <- fisher.test(matrix(c(a, b, bg_c, bg_n - bg_c), 2, byrow = TRUE))$p.value
  data.frame(family = fm, n_H4K20me3_lost = a + b, co_loss_frac = a / (a + b),
             background_frac = bg_rate, odds_ratio = or, p_value = p)
}))
vsbg$p_adj_BH <- p.adjust(vsbg$p_value, "BH")
fwrite(vsbg, file.path(tables_dir, "family_coloss_vs_background.tsv"), sep = "\t")

# SECONDARY: do the families differ from EACH OTHER? (chi-square + pairwise Fisher)
ctab <- cw[k20_lost == TRUE, table(family, k9_lost)]
chi  <- suppressWarnings(chisq.test(ctab))
fams <- as.character(prop$family)
pw <- if (length(fams) >= 2) do.call(rbind, combn(fams, 2, simplify = FALSE, FUN = function(pr) {
  sub <- cw[k20_lost == TRUE & family %in% pr]
  ft  <- fisher.test(table(factor(as.character(sub$family), levels = pr), sub$k9_lost))
  data.frame(family_a = pr[1], family_b = pr[2],
             odds_ratio = unname(ft$estimate), p_value = ft$p.value)
})) else NULL
if (!is.null(pw)) { pw$p_adj_BH <- p.adjust(pw$p_value, "BH")
  fwrite(pw, file.path(tables_dir, "family_coloss_pairwise_fisher.tsv"), sep = "\t") }
fwrite(prop, file.path(tables_dir, "family_coloss_propensity.tsv"), sep = "\t")
cat(sprintf("\n=== H3K9me3 co-loss among H4K20me3-lost genes ===\n(background rate = %.0f%% of %d random-gene loci; chi-square across families p = %.2g)\n",
            100 * bg_rate, bg_n, chi$p.value))
print(as.data.frame(vsbg))

# Co-locate a readable stats summary NEXT TO the figures (gene_families/)
gf_dir <- file.path(fig_dir, "gene_families")
if (!dir.exists(gf_dir)) dir.create(gf_dir, recursive = TRUE, showWarnings = FALSE)
stats_txt <- c(
  "Family co-loss test - H3K9me3 co-loss among H4K20me3-lost genes",
  sprintf("(loss = signal-based log2(KO/WT) < %.2f)", loss_cut), "",
  sprintf("PRIMARY - each family vs genome-wide background (%d random non-family gene loci):",
          bg_n),
  sprintf("  Background co-loss rate = %.1f%% (%d/%d H4K20me3-lost background genes also lose H3K9me3)",
          100 * bg_rate, bg_c, bg_n),
  "  OR>1 = family co-loses MORE than background; Fisher p, BH-adjusted:",
  capture.output(print(as.data.frame(vsbg), row.names = FALSE)), "",
  "SECONDARY - do families differ from each other?",
  sprintf("  Chi-square (family x also-lose-H3K9me3): X-squared = %.1f, df = %d, p = %.3g",
          chi$statistic, chi$parameter, chi$p.value),
  "  Co-loss propensity per family (fraction of H4K20me3-lost that also lose H3K9me3, 95% binom CI):",
  capture.output(print(as.data.frame(prop), row.names = FALSE)),
  "  Pairwise Fisher between families (BH-adjusted):",
  if (!is.null(pw)) capture.output(print(as.data.frame(pw), row.names = FALSE)) else "  (n/a)")
writeLines(stats_txt, file.path(gf_dir, "family_coloss_stats.txt"))

# Plot A: co-loss propensity per family vs the genome-wide background rate.
# Bars = fraction of a family's H4K20me3-lost genes that also lose H3K9me3;
# dashed line = background rate; * marks families significantly above background.
prop_sig <- merge(prop, as.data.table(vsbg)[, .(family, or = odds_ratio, padj = p_adj_BH)],
                  by = "family", all.x = TRUE)
g_prop <- ggplot(prop_sig, aes(reorder(family, -co_loss_frac), co_loss_frac, fill = family)) +
  geom_hline(yintercept = bg_rate, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25) +
  geom_text(aes(label = sprintf("%.0f%%%s\n(n=%d)", 100 * co_loss_frac,
                                ifelse(!is.na(padj) & padj < 0.05, "*", ""), n_H4K20me3_lost)),
            vjust = -0.3, size = 3, lineheight = 0.9) +
  scale_fill_viridis_d(option = "D", end = 0.9, guide = "none") +
  ylim(0, 1.15) +
  labs(title = "H3K9me3 co-loss propensity by family vs genome-wide background",
       subtitle = sprintf("bars = H4K20me3-lost genes also losing H3K9me3; dashed = background (%.0f%%); * = Fisher vs background BH<0.05 (family_coloss_stats.txt)",
                          100 * bg_rate),
       x = NULL, y = "fraction also losing H3K9me3") +
  theme_m
save_fig(g_prop, "family_coloss_propensity", "gene_families", width = 6.5, height = 4)

# Plot B: dH3K9me3 vs dH4K20me3 per gene, faceted by family (the quantitative heatmap)
cw[, class := factor(class, levels = c("co-loss (both)", "H4K20me3-only", "H3K9me3-only", "neither"))]
class_cols <- c("co-loss (both)" = "#b2182b", "H4K20me3-only" = "#2166ac",
                "H3K9me3-only" = "#f4a582", "neither" = "grey75")
g_cw <- ggplot(cw, aes(H4K20me3, H3K9me3, colour = class)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_point(size = 0.6, alpha = 0.5) +
  scale_colour_manual(values = class_cols, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  facet_wrap(~family, nrow = 1) +
  labs(title = "Per-gene H3K9me3 vs H4K20me3 change by family",
       subtitle = "lines at 0 (no change); colour = class at loss cutoff log2 KO/WT < -0.3 (bottom-left = co-loss)",
       x = "delta H4K20me3 (log2 KO/WT)", y = "delta H3K9me3 (log2 KO/WT)") +
  theme_m
save_fig(g_cw, "family_coloss_scatter", "gene_families", width = 14, height = 3.6)
