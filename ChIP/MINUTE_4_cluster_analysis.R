# ================================================================
# MINUTE_4 - Characterise the k-means clusters from MINUTE_3
# ----------------------------------------------------------------
# Runs after MINUTE_3 (which persists results/rds/cluster_analysis_inputs.rds).
# Answers: what IS each chromatin-signature cluster, and what genomic features
# distinguish them? Operates purely on the persisted cluster table + per-sample
# signal matrix - no bigWigs re-read.
#
# Inputs (from MINUTE_3):
#   sig_df      one row per clustered significant peak: Cluster + annotation
#               (genomic_region, repeat_class/family, SYMBOL, chromHMM_state,
#                hmm_<state> coverage fractions, log2FC, baseMean, peak_size_kb)
#   signal_mat  peaks x samples raw ChIP signal (rownames = peak_id)
#   col_annot   per-sample metadata (Sample, ChIP, Genotype, Replicate)
# ================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(data.table) })

inp        <- readRDS(file.path(rds_dir, "cluster_analysis_inputs.rds"))
sig_df     <- inp$sig_df
signal_mat <- inp$signal_mat
col_annot  <- inp$col_annot

sig_df <- sig_df[!is.na(sig_df$Cluster), ]
sig_df$Cluster <- factor(sig_df$Cluster, levels = sort(unique(sig_df$Cluster)))
clusters <- levels(sig_df$Cluster)
cat("Clusters:", paste(clusters, collapse = ", "),
    "| sizes:", paste(as.integer(table(sig_df$Cluster)), collapse = ", "), "\n")

geno_colors <- c(WT = "#1f78b4", HP1gKO = "#e31a1c")
theme_m     <- theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank())

# ================================================================
# 1. Signal profile - mean raw signal per cluster x mark x genotype.
#    Defines each cluster: which marks are present, and how WT->KO changes them.
# ================================================================
sm <- as.data.table(signal_mat, keep.rownames = "peak_id")
sm <- melt(sm, id.vars = "peak_id", variable.name = "Sample", value.name = "signal")
sm[, Sample   := as.character(Sample)]
sm[, Cluster  := sig_df$Cluster[match(peak_id, sig_df$peak_id)]]
sm[, ChIP     := as.character(col_annot$ChIP[match(Sample, col_annot$Sample)])]
sm[, Genotype := as.character(col_annot$Genotype[match(Sample, col_annot$Sample)])]

prof <- sm[!is.na(Cluster),
           .(mean_signal = mean(signal, na.rm = TRUE)),
           by = .(Cluster, ChIP, Genotype)]
prof[, Genotype := factor(Genotype, levels = c("WT", "HP1gKO"))]
prof[, ChIP := factor(ChIP, levels = c("H3K4me3", "H3K36me3", "H3K9me2", "H3K9me3", "H4K20me3"))]

g_prof <- ggplot(prof, aes(Cluster, mean_signal, fill = Genotype)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  facet_wrap(~ChIP, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = geno_colors) +
  labs(title = "Mean ChIP signal per cluster (WT vs HP1gKO)",
       subtitle = "which marks define each cluster, and which are lost in the knockout",
       x = "Cluster", y = "mean signal") +
  theme_m
ggsave(file.path(fig_dir, "cluster_signal_profile.png"), g_prof, width = 14, height = 4, dpi = 300)
message("Saved: cluster_signal_profile.png")

# ================================================================
# 2. Loss heatmap - log2(HP1gKO / WT) mean signal per cluster x mark.
# ================================================================
w <- dcast(prof, Cluster + ChIP ~ Genotype, value.var = "mean_signal")
eps <- 1e-3
w[, log2_KO_WT := log2((HP1gKO + eps) / (WT + eps))]

g_loss <- ggplot(w, aes(ChIP, Cluster, fill = log2_KO_WT)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", log2_KO_WT)), size = 3) +
  scale_fill_gradient2(low = "#2166ac", mid = "grey95", high = "#b2182b",
                       midpoint = 0, name = "log2 KO/WT") +
  labs(title = "Signal change per cluster (log2 HP1gKO / WT)",
       x = "Mark", y = "Cluster") +
  theme_m
ggsave(file.path(fig_dir, "cluster_loss_heatmap.png"), g_loss, width = 7, height = 4.5, dpi = 300)
message("Saved: cluster_loss_heatmap.png")

# ================================================================
# 3. ChromHMM profile - mean per-state coverage per cluster.
#    Shown two ways: absolute mean coverage, and z-scored per state (relative
#    enrichment across clusters, so genome-abundant Quies doesn't dominate).
# ================================================================
hmm_cols <- grep("^hmm_", names(sig_df), value = TRUE)
hmm_prof <- sig_df %>% group_by(Cluster) %>%
  summarise(across(all_of(hmm_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
hmm_mat <- as.matrix(hmm_prof[, hmm_cols]); rownames(hmm_mat) <- hmm_prof$Cluster
colnames(hmm_mat) <- sub("^hmm_", "", colnames(hmm_mat))

mk_hmm_long <- function(mat, value_name) {
  d <- as.data.table(mat, keep.rownames = "Cluster")
  d <- melt(d, id.vars = "Cluster", variable.name = "state", value.name = value_name)
  d$Cluster <- factor(d$Cluster, levels = clusters)
  d
}
hmm_abs <- mk_hmm_long(hmm_mat, "cov")
hmm_z   <- mk_hmm_long(scale(hmm_mat), "z")   # per-state z across clusters

g_hmm_abs <- ggplot(hmm_abs, aes(state, Cluster, fill = cov)) +
  geom_tile(colour = "white") +
  scale_fill_viridis_c(name = "mean\ncoverage") +
  labs(title = "Mean ChromHMM state coverage per cluster (absolute)",
       x = "ChromHMM state", y = "Cluster") +
  theme_m + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g_hmm_z <- ggplot(hmm_z, aes(state, Cluster, fill = z)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(low = "#2166ac", mid = "grey95", high = "#b2182b",
                       midpoint = 0, name = "z\n(per state)") +
  labs(title = "ChromHMM state enrichment per cluster (z-scored per state)",
       subtitle = "which cluster is relatively high/low in each state",
       x = "ChromHMM state", y = "Cluster") +
  theme_m + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "cluster_chromHMM_coverage_abs.png"), g_hmm_abs, width = 9, height = 4, dpi = 300)
ggsave(file.path(fig_dir, "cluster_chromHMM_coverage_zscore.png"), g_hmm_z, width = 9, height = 4, dpi = 300)
message("Saved: cluster_chromHMM_coverage_{abs,zscore}.png")

# ================================================================
# 4. Genomic composition per cluster: (a) repeat_class %  (b) genomic_region %
#    (per-cluster gene-family fractions removed - superseded by the dedicated
#     family exon analysis in part 7)
# ================================================================
comp_bar <- function(df, var, title) {
  d <- df %>% mutate(v = ifelse(is.na(.data[[var]]), "NA", as.character(.data[[var]]))) %>%
    count(Cluster, v) %>% group_by(Cluster) %>% mutate(frac = n / sum(n)) %>% ungroup()
  ggplot(d, aes(Cluster, frac, fill = v)) +
    geom_col() +
    labs(title = title, x = "Cluster", y = "fraction of peaks", fill = var) +
    theme_m
}
ggsave(file.path(fig_dir, "cluster_repeat_composition.png"),
       comp_bar(sig_df, "repeat_class", "Repeat-class composition per cluster"),
       width = 7, height = 4.5, dpi = 300)
ggsave(file.path(fig_dir, "cluster_region_composition.png"),
       comp_bar(sig_df, "genomic_region", "Genomic-region composition per cluster"),
       width = 8, height = 4.5, dpi = 300)
message("Saved: cluster_repeat_composition.png + cluster_region_composition.png")

# ================================================================
# 5. Cluster summary table
# ================================================================
summ <- sig_df %>% group_by(Cluster) %>%
  summarise(
    n_peaks       = n(),
    median_kb     = round(median(peak_size_kb, na.rm = TRUE), 1),
    median_l2fc   = round(median(log2FoldChange, na.rm = TRUE), 2),
    top_state     = names(sort(table(chromHMM_state), decreasing = TRUE))[1],
    top_repeat    = names(sort(table(repeat_class),  decreasing = TRUE))[1],
    pct_pericentro = round(100 * mean(start < 3e6), 1),
    pct_Pcdh      = round(100 * mean(grepl("^Pcdh", ifelse(is.na(SYMBOL), "", SYMBOL))), 1),
    pct_Zfp       = round(100 * mean(grepl("^Zfp",  ifelse(is.na(SYMBOL), "", SYMBOL))), 1),
    pct_Vmn       = round(100 * mean(grepl("^Vmn",  ifelse(is.na(SYMBOL), "", SYMBOL))), 1),
    .groups = "drop")
fwrite(summ, file.path(tables_dir, "cluster_summary.tsv"), sep = "\t")
cat("\n=== Cluster summary ===\n"); print(as.data.frame(summ))
message("Saved: cluster_summary.tsv")

# ================================================================
# 6. Per-cluster ChromHMM enrichment TEST
#    Formalises the part-3 coverage heatmaps: for each cluster x state, compare
#    the state's coverage in that cluster vs all OTHER clustered peaks (Wilcoxon
#    + BH FDR, log2 mean-coverage ratio). Tells which states each cluster is
#    distinctively enriched/depleted for, relative to the other clusters.
# ================================================================
clus_enrich <- do.call(rbind, lapply(clusters, function(cl) {
  inc <- sig_df$Cluster == cl
  # per-region (unweighted) Wilcoxon: this cluster vs all other clustered peaks
  pr <- do.call(rbind, lapply(hmm_cols, function(col) {
    a <- sig_df[[col]][inc]; b <- sig_df[[col]][!inc]
    data.frame(Cluster = cl, state = sub("^hmm_", "", col), n = sum(inc),
               mean_cov_in = mean(a, na.rm = TRUE), mean_cov_out = mean(b, na.rm = TRUE),
               log2_ratio = log2((mean(a, na.rm = TRUE) + 1e-4) / (mean(b, na.rm = TRUE) + 1e-4)),
               p_value = tryCatch(suppressWarnings(wilcox.test(a, b)$p.value),
                                  error = function(e) NA_real_),
               stringsAsFactors = FALSE)
  }))
  # size-weighted (by domain length) with permutation test
  sw <- chromHMM_size_weighted(sig_df[, hmm_cols], sig_df$peak_size_kb, inc)
  merge(pr, sw, by = "state")
}))
clus_enrich <- clus_enrich %>% group_by(Cluster) %>%
  mutate(p_adj_BH   = p.adjust(p_value, method = "BH"),
         perm_p_adj = p.adjust(perm_p,  method = "BH")) %>% ungroup()
fwrite(clus_enrich, file.path(tables_dir, "cluster_chromHMM_enrichment.tsv"), sep = "\t")

# Two dot plots: per-region (Wilcoxon) and size-weighted (permutation)
clus_dotplot <- function(lr, fdr, ttl, sub) {
  clus_enrich %>%
    filter(is.finite(.data[[lr]]), !is.na(.data[[fdr]])) %>%
    mutate(Cluster = factor(paste("Cluster", Cluster), levels = paste("Cluster", clusters)),
           neglog10FDR = -log10(pmax(.data[[fdr]], .Machine$double.xmin)),
           state = factor(state, levels = sort(unique(state)))) %>%
    ggplot(aes(Cluster, state, size = neglog10FDR, colour = .data[[lr]])) +
    geom_point() +
    scale_colour_gradient2(low = "steelblue3", mid = "grey90", high = "firebrick2",
                           midpoint = 0, name = expression(log[2]~"ratio")) +
    scale_size_continuous(name = expression(-log[10]("FDR"))) +
    labs(title = ttl, subtitle = sub, x = NULL, y = "ChromHMM state") +
    theme_m + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}
ggsave(file.path(fig_dir, "cluster_chromHMM_enrichment.png"),
       clus_dotplot("log2_ratio", "p_adj_BH",
                    "Per-cluster ChromHMM enrichment (cluster vs other clusters) - PER-REGION",
                    "colour = log2 coverage ratio; each region equal; size = FDR"),
       width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "cluster_chromHMM_enrichment_sizeweighted.png"),
       clus_dotplot("w_log2_ratio", "perm_p_adj",
                    "Per-cluster ChromHMM enrichment (cluster vs other clusters) - SIZE-WEIGHTED",
                    "colour = log2 territory ratio (weighted by domain length); size = perm FDR"),
       width = 8, height = 6, dpi = 300)
message("Saved: cluster_chromHMM_enrichment.{tsv, per-region png, size-weighted png}")

# ================================================================
# 7. Gene-family exon signal (HUSH / CBX3-silenced families)
#    KRAB-ZFP, clustered protocadherin (chr18), vomeronasal and olfactory-
#    receptor genes are normally silenced by H3K9me3 / HUSH. Rather than relying
#    on individual genes passing per-peak significance, quantify ChIP signal
#    DIRECTLY over the SILENCING-RELEVANT exon of each gene (config's
#    family_exon_rule: KRAB-ZFP = last exon, protocadherin = first exon,
#    Vmn/Olfr = single coding exon) and compare WT vs HP1gKO per family x mark.
#    One interval per gene.  << reads bigWigs, unlike parts 1-6 >>
# ================================================================
mark_levels <- c("H3K4me3", "H3K36me3", "H3K9me2", "H3K9me3", "H4K20me3")

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
ggsave(file.path(fig_dir, "family_exon_change.png"), g_fam_chg, width = 14, height = 4, dpi = 300)

# Plot B: WT vs KO signal DISTRIBUTION per gene, per family x mark (boxplot, so
# the across-gene spread is visible rather than a single median bar)
lvl <- melt(as.data.table(long[, c("family", "mark", "gene", "wt_signal", "ko_signal")]),
            id.vars = c("family", "mark", "gene"), variable.name = "Genotype", value.name = "signal")
lvl[, Genotype := factor(ifelse(Genotype == "wt_signal", "WT", "HP1gKO"), levels = c("WT", "HP1gKO"))]
g_fam_lvl <- ggplot(lvl, aes(mark, signal, fill = Genotype)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.15, linewidth = 0.3,
               position = position_dodge(0.8)) +
  facet_wrap(~family, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = geno_colors) +
  labs(title = "ChIP signal over gene-family exons, per gene (WT vs HP1gKO)",
       subtitle = "boxes = across-gene distribution of per-gene mean signal",
       x = "Mark", y = "signal") +
  theme_m + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
ggsave(file.path(fig_dir, "family_exon_signal_level.png"), g_fam_lvl, width = 14, height = 4, dpi = 300)

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
ggsave(file.path(fig_dir, "family_exon_dotplot.png"), g_fam_dot, width = 7, height = 3.5, dpi = 300)

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
ggsave(file.path(fig_dir, "family_exon_MA.png"), g_fam_ma, width = 14, height = 3.6, dpi = 300)

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

draw_fam_hm <- function(mat, cols, ttl, outfile, w) {
  m <- t(scale(t(mat[, cols, drop = FALSE]))); m[!is.finite(m)] <- 0
  ha <- HeatmapAnnotation(Genotype = col_meta$Genotype[cols],
                          col = list(Genotype = geno_colors), show_annotation_name = TRUE)
  ht <- Heatmap(m, name = "row z", top_annotation = ha,
                column_split = droplevels(col_meta$Mark[cols]),
                row_split = row_fam, cluster_columns = FALSE, cluster_rows = TRUE,
                show_row_names = FALSE, show_column_names = FALSE,
                row_title_gp = gpar(fontsize = 9), column_title_gp = gpar(fontsize = 10),
                heatmap_legend_param = list(title = "row z-score"))
  png(file.path(fig_dir, outfile), width = w, height = 11, units = "in", res = 200)
  draw(ht, column_title = ttl); dev.off()
}
draw_fam_hm(fam_mat, seq_len(ncol(fam_mat)),
            "Gene-family exon signal per replicate (all marks)",
            "family_exon_heatmap_allmarks.png", 11)
sil <- which(col_meta$Mark %in% c("H3K9me3", "H4K20me3"))
draw_fam_hm(fam_mat, sil,
            "Gene-family exon signal per replicate (silencing marks)",
            "family_exon_heatmap_silencing.png", 7)

message("Saved: family_exon_signal.{tsv,rds} + summary.tsv + 6 figures (box, level, dot, MA, 2 heatmaps)")
