# ================================================================
# MINUTE_7 — direct repeat-class signal change (bw_loci over repeat copies)
# ----------------------------------------------------------------
# Complements MINUTE_2's hypergeometric test. That test asks "are *significant
# peaks* over-represented on repeat class X" - it is peak-calling-dependent and
# magnitude-blind. Here we instead quantify the input-scaled signal DIRECTLY over
# the repeat copies themselves (bw_loci, exactly as for the gene-family exons in
# MINUTE_4) and measure how each mark CHANGES over each class, independent of peak
# calling. Because sizeFactors = 1, absolute signal is comparable across genotype.
#
# Two readouts per class x mark:
#   (1) AGGREGATE (primary) - mean signal across all copies per sample, then
#       log2(mean KO / mean WT). Averaging over copies cancels the per-copy
#       mappability noise of near-identical young repeats (the WT->KO change is
#       internally controlled, so this stays interpretable even where absolute
#       signal is unique-mapping-limited).
#   (2) PER-COPY (secondary) - per-copy log2FC distribution + paired Wilcoxon,
#       over copies with measurable signal (drops the unmappable bottom quartile).
# ================================================================
source("config.R")
suppressPackageStartupMessages({
  library(GenomicRanges); library(ggplot2); library(data.table); library(dplyr)
})

N_MAX <- 4000    # cap copies per class (subsample large classes for tractability)
eps   <- 0.25
set.seed(42)

line <- loadRepeatBED(repeat_bed_files$LINE)
ltr  <- loadRepeatBED(repeat_bed_files$LTR)
sine <- loadRepeatBED(repeat_bed_files$SINE)
pick <- function(gr, pat, field = "repName") gr[grepl(pat, mcols(gr)[[field]])]

# Classes: young L1 (HUSH substrates) + all L1, IAP/ERVK, older ERV, SINE controls
class_defs <- list(
  "L1MdA (young L1)" = pick(line, "^L1MdA"),
  "L1MdT (young L1)" = pick(line, "^L1MdT"),
  "L1MdGf (young L1)"= pick(line, "^L1MdGf"),
  "L1MdF"            = pick(line, "^L1MdF"),
  "L1 (all)"         = pick(line, "^L1$", "repFamily"),
  "IAPEz-int"        = pick(ltr,  "IAPEz-int"),
  "IAPLTR1a_Mm"      = pick(ltr,  "IAPLTR1a_Mm"),
  "ERVK (all)"       = pick(ltr,  "^ERVK$", "repFamily"),
  "ERVL-MaLR"        = pick(ltr,  "ERVL-MaLR", "repFamily"),
  "B1 SINE (ctrl)"   = pick(sine, "^B1"),
  "B2 SINE (ctrl)"   = pick(sine, "^B2"))
class_defs <- class_defs[vapply(class_defs, length, 0L) > 0]      # drop empty
cat("Repeat-class copy counts (pre-subsample):\n")
print(vapply(class_defs, length, 0L))

# subsample + combine into ONE GRanges (quantify in a single bw_loci pass/sample)
combined <- do.call(c, unname(lapply(names(class_defs), function(cl) {
  g <- class_defs[[cl]]
  if (length(g) > N_MAX) g <- g[sort(sample(length(g), N_MAX))]
  mcols(g) <- DataFrame(rep_class = rep(cl, length(g)))   # one label per copy
  g
})))
seqlevelsStyle(combined) <- "NCBI"
combined <- keepStandardChromosomes(combined, pruning.mode = "coarse")
combined <- combined[as.character(seqnames(combined)) %in% c(1:19, "X", "Y")]  # bigWig seqnames
cl_vec <- as.character(mcols(combined)$rep_class)

# The parent classes (L1 (all), ERVK (all)) share copies with their subfamilies, so
# `combined` has duplicate coordinates. Quantify the UNIQUE loci once, then expand
# back per-copy. Align each bw_loci result by coordinate key (robust to wigglescout
# dropping/reordering loci), so rep_mat always matches cl_vec length.
key    <- paste(seqnames(combined), start(combined), end(combined), sep = ":")
uni    <- !duplicated(key)
uni_gr <- combined[uni]; uni_key <- key[uni]
cat(sprintf("\nQuantifying %d unique loci (%d copies across %d classes) x %d samples ...\n",
            length(uni_gr), length(combined), length(class_defs), nrow(samples)))
uni_signal <- parallel::mclapply(seq_len(nrow(samples)), function(i) {
  bw <- samples$bigwig_path[i]
  if (!file.exists(bw)) return(rep(NA_real_, length(uni_gr)))
  r <- tryCatch(bw_loci(bw, uni_gr), error = function(e) NULL)
  if (is.null(r)) return(rep(NA_real_, length(uni_gr)))
  v  <- rep(NA_real_, length(uni_gr))
  v[match(paste(seqnames(r), start(r), end(r), sep = ":"), uni_key)] <- mcols(r)[[1]]
  v
}, mc.cores = parallel::detectCores())
uni_mat <- do.call(cbind, uni_signal); stopifnot(nrow(uni_mat) == length(uni_gr))
rep_mat <- uni_mat[match(key, uni_key), , drop = FALSE]   # expand to per-copy
colnames(rep_mat) <- samples$sample_id
stopifnot(nrow(rep_mat) == length(cl_vec))
saveRDS(list(gr = combined, mat = rep_mat, rep_class = cl_vec),
        file.path(rds_dir, "repeat_class_signal.rds"))

marks_here <- intersect(c("H3K4me3","H3K36me3","H3K9me2","H3K9me3","H4K20me3"),
                        as.character(samples$mark))

# ---- (1) AGGREGATE per class x mark: mean signal over copies -> log2(KO/WT) ----
agg <- do.call(rbind, lapply(marks_here, function(mk) {
  wtc <- samples$sample_id[as.character(samples$mark) == mk & as.character(samples$genotype) == "WT"]
  koc <- samples$sample_id[as.character(samples$mark) == mk & as.character(samples$genotype) == "HP1gKO"]
  # per-class mean signal per sample, then average within genotype
  do.call(rbind, lapply(unique(cl_vec), function(cl) {
    ix <- cl_vec == cl
    wt <- mean(rep_mat[ix, wtc], na.rm = TRUE); ko <- mean(rep_mat[ix, koc], na.rm = TRUE)
    data.frame(rep_class = cl, mark = mk, n_copies = sum(ix),
               mean_WT = wt, mean_KO = ko, log2FC = log2((ko + eps) / (wt + eps)))
  }))
}))
fwrite(agg, file.path(tables_dir, "repeat_signal_aggregate.tsv"), sep = "\t")
cat("\n=== Repeat-class AGGREGATE signal change (log2 KO/WT) ===\n")
print(dcast(as.data.table(agg), rep_class ~ mark, value.var = "log2FC"), row.names = FALSE)

# ---- (2) PER-COPY log2FC + paired Wilcoxon (copies with measurable signal) ----
percopy <- do.call(rbind, lapply(marks_here, function(mk) {
  wtc <- samples$sample_id[as.character(samples$mark) == mk & as.character(samples$genotype) == "WT"]
  koc <- samples$sample_id[as.character(samples$mark) == mk & as.character(samples$genotype) == "HP1gKO"]
  wt <- rowMeans(rep_mat[, wtc, drop = FALSE], na.rm = TRUE)
  ko <- rowMeans(rep_mat[, koc, drop = FALSE], na.rm = TRUE)
  data.frame(rep_class = cl_vec, mark = mk, wt = wt, ko = ko,
             abund = (wt + ko) / 2, log2FC = log2((ko + eps) / (wt + eps)))
}))
# drop the unmappable/empty bottom quartile per mark (near-zero in both genotypes)
percopy <- as.data.table(percopy)[, keep := abund > quantile(abund, 0.25, na.rm = TRUE), by = mark][keep == TRUE]
copy_summary <- percopy[, {
  p <- tryCatch(suppressWarnings(wilcox.test(ko, wt, paired = TRUE)$p.value), error = function(e) NA_real_)
  .(n = .N, median_log2FC = median(log2FC, na.rm = TRUE), wilcox_p = p)
}, by = .(rep_class, mark)][, p_adj_BH := p.adjust(wilcox_p, "BH"), by = mark]
fwrite(copy_summary, file.path(tables_dir, "repeat_signal_percopy_summary.tsv"), sep = "\t")

# ---- Plots -----------------------------------------------------------------
cls_order <- names(class_defs)
mk_order  <- marks_here
agg$rep_class <- factor(agg$rep_class, levels = cls_order); agg$mark <- factor(agg$mark, levels = mk_order)
percopy[, rep_class := factor(rep_class, levels = cls_order)][, mark := factor(mark, levels = mk_order)]

# Plot A (primary): aggregate log2FC tile heatmap, class x mark
gA <- ggplot(agg, aes(mark, rep_class, fill = log2FC)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%+.2f", log2FC)), size = 2.9) +
  scale_fill_gradient2(low = "#2166ac", mid = "grey95", high = "#b2182b", midpoint = 0,
                       name = "log2(KO/WT)") +
  scale_y_discrete(limits = rev(cls_order)) +
  labs(title = "Direct repeat-class signal change (aggregate over copies)",
       subtitle = "mean input-scaled signal over all copies, log2(KO/WT); negative = mark lost over the class",
       x = "Mark", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(gA, "repeat_signal_aggregate", "repeats", width = 7.5, height = 5)

# Plot B (secondary): per-copy log2FC distribution, faceted by mark
gB <- ggplot(percopy, aes(rep_class, log2FC)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_boxplot(aes(fill = mark), outlier.size = 0.2, outlier.alpha = 0.1, linewidth = 0.3) +
  facet_wrap(~mark, nrow = 1) +
  scale_fill_viridis_d(option = "D", end = 0.9, guide = "none") +
  coord_flip(ylim = c(-2, 1.5)) +
  labs(title = "Per-copy repeat signal change by class and mark",
       subtitle = "per-copy log2(KO/WT) over copies with measurable signal (bottom-quartile abundance dropped)",
       x = NULL, y = "log2(KO / WT) per copy") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank())
save_fig(gB, "repeat_signal_percopy", "repeats", width = 15, height = 4.5)

message("Saved: repeat_signal_{aggregate,percopy_summary}.tsv + rds/repeat_class_signal.rds")
message("Saved: repeats/repeat_signal_{aggregate,percopy}.{png,pdf}")
