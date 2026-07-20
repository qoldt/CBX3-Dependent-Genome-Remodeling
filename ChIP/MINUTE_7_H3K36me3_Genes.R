# ================================================================
# MINUTE_7 - H3K36me3 over GENE BODIES
# ----------------------------------------------------------------
# WHY THIS STAGE EXISTS. Per-peak testing is close to the worst case for
# H3K36me3: the mark is a broad co-transcriptional domain over gene bodies, but
# the peak set fragments it into pieces with a median width of ~0.3-1 kb, and
# each fragment is then tested on its own. The result is degenerate - only 15 of
# 77,770 peaks reach p < 0.05 (0.019%), where even a pure null predicts 5%. That
# is a ~260x deficit, i.e. the test is not merely underpowered but badly
# conservative, so "no significant peaks" carries almost no information.
#
# WHAT THIS DOES INSTEAD. Quantify H3K36me3 over whole gene bodies - the unit the
# mark actually occupies - and test per gene. Averaging over the whole domain
# rather than a 0.3 kb fragment reduces the per-unit noise, and there are ~10x
# fewer tests.
#
# WHY IT MATTERS BIOLOGICALLY. CBX3/HP1gamma is not only a heterochromatin
# protein: it is recruited to the bodies of actively transcribed genes and has
# reported roles in transcription elongation and co-transcriptional splicing. So
# a genic H3K36me3 effect is a plausible prior, not a fishing expedition. The
# per-peak data already hint at one - transcribed ChromHMM states lose more
# (Tx -0.093, TxWk -0.078) than distal intergenic (0.000) - but the effect is
# ~10x smaller than H4K20me3's and cannot be resolved peak-by-peak.
#
# WHAT THIS STAGE CANNOT DO. It cannot separate "H3K36me3 changed" from "the gene
# changed expression": H3K36me3 is deposited co-transcriptionally, so a drop can
# reflect less transcription rather than altered methylation. Distinguishing
# those needs RNA-seq. Every conclusion here is about the MARK, not the mechanism.
#
# Re-reads the H3K36me3 bigWigs (gene bodies are not in the peak-based counts).
# Independent of MINUTE_2-6; needs only config.R + the bigWigs.
# ================================================================
source("config.R")
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(data.table)
  library(TxDb.Mmusculus.UCSC.mm39.knownGene); library(org.Mm.eg.db)
})

MARK        <- "H3K36me3"
MIN_GENE_KB <- 1      # genes shorter than this carry too little domain to average
MIN_COV     <- 0.5    # WT mean coverage below this = not meaningfully marked
fig_sub     <- "h3k36me3_genes"

# --- 1. Gene bodies ---------------------------------------------------------
txdb  <- TxDb.Mmusculus.UCSC.mm39.knownGene
genes_gr <- genes(txdb, single.strand.genes.only = TRUE)
mcols(genes_gr)$symbol <- suppressMessages(AnnotationDbi::mapIds(
  org.Mm.eg.db, keys = names(genes_gr), column = "SYMBOL", keytype = "ENTREZID"))
genes_gr <- genes_gr[!is.na(mcols(genes_gr)$symbol) &
                     width(genes_gr) >= MIN_GENE_KB * 1000]
genes_gr <- keepStandardChromosomes(genes_gr, pruning.mode = "coarse")
seqlevelsStyle(genes_gr) <- "NCBI"
cat(sprintf("Gene bodies: %d (>= %g kb, standard chromosomes)\n", length(genes_gr), MIN_GENE_KB))

# --- 2. Quantify (bw_loci returns MEAN coverage, i.e. length-normalised) -----
bw <- get_bigwig_files(MARK)
cat(sprintf("Quantifying %s over gene bodies for %d samples...\n", MARK, length(bw)))
cov_list <- parallel::mclapply(seq_along(bw), function(i) {
  if (!file.exists(bw[[i]])) { message("Missing: ", bw[[i]]); return(rep(NA_real_, length(genes_gr))) }
  message(sprintf("[%d/%d] %s", i, length(bw), names(bw)[i]))
  tryCatch(mcols(wigglescout::bw_loci(bw[[i]], genes_gr))[[1]],
           error = function(e) { message("Error ", names(bw)[i], ": ", e$message)
                                 rep(NA_real_, length(genes_gr)) })
}, mc.cores = max(1, parallel::detectCores() - 1), mc.preschedule = FALSE)
cov_mat <- do.call(cbind, cov_list); colnames(cov_mat) <- names(bw)
rownames(cov_mat) <- mcols(genes_gr)$symbol
stopifnot(ncol(cov_mat) == length(bw))

gene_meta <- data.frame(
  symbol  = mcols(genes_gr)$symbol,
  entrez  = names(genes_gr),
  chr     = as.character(seqnames(genes_gr)),
  gene_kb = width(genes_gr) / 1000,
  stringsAsFactors = FALSE)

wt_ids  <- samples$sample_id[samples$mark == MARK & samples$genotype == "WT"]
ko_ids  <- samples$sample_id[samples$mark == MARK & samples$genotype == "HP1gKO"]
gene_meta$wt_cov <- rowMeans(cov_mat[, wt_ids, drop = FALSE], na.rm = TRUE)
gene_meta$ko_cov <- rowMeans(cov_mat[, ko_ids, drop = FALSE], na.rm = TRUE)

keep <- is.finite(gene_meta$wt_cov) & gene_meta$wt_cov >= MIN_COV & complete.cases(cov_mat)
cat(sprintf("Genes with WT mean coverage >= %g: %d of %d\n", MIN_COV, sum(keep), nrow(gene_meta)))
cov_mat   <- cov_mat[keep, , drop = FALSE]
gene_meta <- gene_meta[keep, ]

# --- 3. DESeq2 per gene -----------------------------------------------------
# Same conventions as MINUTE_1: sizeFactors = 1 (bigWigs are already
# input-scaled), contrast HP1gKO vs WT.
#
# COUNT_SCALE: mean coverage over a gene body is ~0.5-3 here, and rounding that
# straight to integers (as MINUTE_1 does per peak) is destructive twice over.
# First it collapses the values onto a handful of integers - with so few distinct
# values the mean-dispersion trend is degenerate and DESeq2's dispersion fit
# fails outright ("newsplit: out of vertex space" from locfit). Second, and more
# importantly, it imposes a FAKE POISSON NOISE FLOOR: a coverage of 2.0 averaged
# over a 30 kb gene body is a far more precise quantity than a raw count of 2,
# but the negative binomial treats it as the latter and inflates the variance
# accordingly. Scaling by a constant before rounding preserves every log2 ratio
# exactly (it cancels in the contrast) while restoring resolution and moving the
# Poisson floor closer to the real precision of the estimate.
COUNT_SCALE <- 10
counts <- round(cov_mat * COUNT_SCALE)
meta   <- data.frame(row.names = colnames(counts),
                     condition = factor(as.character(samples$genotype[
                       match(colnames(counts), samples$sample_id)]), levels = c("WT", "HP1gKO")))
dds <- DESeqDataSetFromMatrix(counts, meta, design = ~ condition)
sizeFactors(dds) <- rep(1, ncol(counts))
# Robust dispersion fitting: parametric -> local -> mean. These are not real
# counts, so the trend can still be badly behaved; fall back rather than fail,
# and say which fit was used so it is not silent.
fit_used <- "parametric"
dds <- tryCatch(DESeq(dds, quiet = TRUE), error = function(e) {
  message("  parametric/local dispersion fit failed (", conditionMessage(e),
          ") - falling back to fitType='mean'")
  fit_used <<- "mean"
  DESeq(dds, fitType = "mean", quiet = TRUE)
})
cat(sprintf("DESeq2 dispersion fitType used: %s (counts scaled x%d before rounding)\n",
            fit_used, COUNT_SCALE))
res <- as.data.frame(results(dds, contrast = c("condition", "HP1gKO", "WT")))
gene_res <- cbind(gene_meta, res[, c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")])
fwrite(gene_res, file.path(tables_dir, "h3k36me3_gene_body_results.tsv"), sep = "\t")

l <- gene_res$log2FoldChange
cat(sprintf("\n=== %s over gene bodies: n=%d  median log2FC %+0.3f  IQR [%+0.3f, %+0.3f]\n",
            MARK, nrow(gene_res), median(l, na.rm = TRUE),
            quantile(l, .25, na.rm = TRUE), quantile(l, .75, na.rm = TRUE)))
cat(sprintf("    %%down %.1f%%   |log2FC|>0.5 %.1f%%   p<0.05 %.2f%%   padj<0.1 %d\n",
            100 * mean(l < 0, na.rm = TRUE), 100 * mean(abs(l) > 0.5, na.rm = TRUE),
            100 * mean(gene_res$pvalue < 0.05, na.rm = TRUE), sum(gene_res$padj < 0.1, na.rm = TRUE)))

# --- 4. Did aggregation actually buy power? ---------------------------------
# The honest comparison against the per-peak run. pi0 is a Storey-style estimate
# of the NULL fraction from the flat right tail; pi0 well below 1 means signal is
# now detectable. If the gene-level p-value histogram is still depleted near 0,
# aggregation did NOT rescue this mark and the limit is the assay, not the unit.
pi0_of <- function(p) min(1, 2 * mean(p > 0.5, na.rm = TRUE))
peak_p <- readRDS(annotated_rds)[[MARK]]$pvalue
cmp <- data.frame(
  unit    = c("per peak (MINUTE_1)", "per gene body (this stage)"),
  n       = c(sum(is.finite(peak_p)), sum(is.finite(gene_res$pvalue))),
  frac_p05 = c(mean(peak_p < 0.05, na.rm = TRUE), mean(gene_res$pvalue < 0.05, na.rm = TRUE)),
  pi0     = c(pi0_of(peak_p), pi0_of(gene_res$pvalue)),
  n_padj10 = c(sum(readRDS(annotated_rds)[[MARK]]$padj < 0.1, na.rm = TRUE),
               sum(gene_res$padj < 0.1, na.rm = TRUE)))
fwrite(cmp, file.path(tables_dir, "h3k36me3_peak_vs_gene_power.tsv"), sep = "\t")
cat("\n=== does gene-body aggregation buy power? (5% of tests are expected below p<0.05 under a pure null) ===\n")
print(cmp, row.names = FALSE)

pdf_df <- rbind(data.frame(unit = "per peak (MINUTE_1)", p = peak_p),
                data.frame(unit = "per gene body (this stage)", p = gene_res$pvalue))
pdf_df <- pdf_df[is.finite(pdf_df$p), ]
g_p <- ggplot(pdf_df, aes(p, fill = unit)) +
  geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 1, 0.02), colour = NA) +
  scale_fill_disc(guide = "none") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey20", linewidth = 0.4) +
  facet_wrap(~unit, ncol = 1, scales = "free_y") +
  labs(title = paste(MARK, "p-value distribution: per peak vs per gene body"),
       subtitle = "dashed = uniform null. Depletion near 0 = conservative/miscalibrated, not just underpowered",
       x = "p-value", y = "density") +
  theme_m
save_fig(g_p, "h3k36me3_pvalue_peak_vs_gene", fig_sub, width = 7, height = 7)

# --- 5. Is the change transcription-coupled? --------------------------------
# H3K36me3 tracks transcription, so WT coverage is a rough expression proxy. If
# loss scales with WT level, the effect follows transcribed territory.
gene_res$wt_bin <- cut(gene_res$wt_cov, quantile(gene_res$wt_cov, seq(0, 1, .2), na.rm = TRUE),
                       include.lowest = TRUE, labels = paste0("Q", 1:5))
gene_res$len_bin <- cut(gene_res$gene_kb, c(0, 5, 20, 50, Inf),
                        labels = c("<5 kb", "5-20", "20-50", ">50"))
by_wt <- gene_res %>% filter(!is.na(wt_bin)) %>% group_by(wt_bin) %>%
  summarise(n = n(), median_log2FC = median(log2FoldChange, na.rm = TRUE),
            median_wt_cov = median(wt_cov), .groups = "drop")
by_len <- gene_res %>% filter(!is.na(len_bin)) %>% group_by(len_bin) %>%
  summarise(n = n(), median_log2FC = median(log2FoldChange, na.rm = TRUE), .groups = "drop")
fwrite(by_wt,  file.path(tables_dir, "h3k36me3_gene_by_wt_signal.tsv"), sep = "\t")
fwrite(by_len, file.path(tables_dir, "h3k36me3_gene_by_length.tsv"),    sep = "\t")
cat("\n=== median log2FC by WT signal quintile (expression proxy) ===\n"); print(as.data.frame(by_wt), row.names = FALSE)
cat("\n=== median log2FC by gene length ===\n"); print(as.data.frame(by_len), row.names = FALSE)

g_wt <- ggplot(gene_res[!is.na(gene_res$wt_bin), ], aes(wt_bin, log2FoldChange, fill = wt_bin)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_boxplot(outlier.size = 0.2, alpha = 0.85) +
  scale_fill_disc(guide = "none") +
  labs(title = paste(MARK, "gene-body change by WT signal level"),
       subtitle = "WT coverage is a proxy for transcription; H3K36me3 is deposited co-transcriptionally",
       x = "WT mean coverage quintile", y = "log2FoldChange (HP1gKO vs WT)") +
  theme_m
save_fig(g_wt, "h3k36me3_gene_change_by_wt_signal", fig_sub, width = 7, height = 5)

# --- 6. MA plot + labelled extremes ----------------------------------------
top <- gene_res %>% filter(is.finite(log2FoldChange), wt_cov > 0) %>%
  slice_max(abs(log2FoldChange), n = 25, with_ties = FALSE)
fwrite(gene_res %>% arrange(log2FoldChange) %>% head(100),
       file.path(tables_dir, "h3k36me3_gene_top100_down.tsv"), sep = "\t")
g_ma <- ggplot(gene_res[gene_res$wt_cov > 0, ], aes(wt_cov, log2FoldChange)) +
  geom_bin2d(bins = 80) +
  scale_fill_heat0(transform = "log10", name = "genes\nper bin") +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_hline(yintercept = median(l, na.rm = TRUE), linetype = "dashed", colour = "grey20") +
  ggrepel::geom_text_repel(data = top, aes(label = symbol), size = 2.4,
                           max.overlaps = Inf, segment.size = 0.2, colour = "black") +
  scale_x_log10() +
  labs(title = paste(MARK, "gene-body MA"),
       subtitle = sprintf("dashed = median log2FC (%+0.3f); n = %s genes",
                          median(l, na.rm = TRUE), format(nrow(gene_res), big.mark = ",")),
       x = "WT mean coverage over gene body (log10)", y = "log2FoldChange") +
  theme_m
save_fig(g_ma, "h3k36me3_gene_MA", fig_sub, width = 8, height = 5.5)

saveRDS(list(gene_res = gene_res, cov_mat = cov_mat),
        file.path(rds_dir, "h3k36me3_gene_body.rds"))
cat("\nSaved: tables/h3k36me3_*.tsv, figures/", fig_sub, "/, rds/h3k36me3_gene_body.rds\n", sep = "")
