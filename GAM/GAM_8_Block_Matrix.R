# ================================================================
# GAM_8 - Block x block contact budget: is this redistribution?
# ----------------------------------------------------------------
#     cd GAM && Rscript GAM_8_Block_Matrix.R      (after GAM_1, GAM_2, GAM_6)
#
# THE QUESTION. GAM_6 found 123 large blocks, 63 losing and 60 gaining contacts,
# covering 41.6% of the genome in near-reciprocal amounts (552 vs 462 Mb). NPMI
# is a normalised measure, so a global REDISTRIBUTION of contact budget - the KO
# simply moving contacts from one kind of territory to another - would produce
# that pattern without any local biology. With two merged biological replicates
# there is no variance to test against, so the discriminator has to be the shape
# of the change, not its significance.
#
# This stage measures the contact budget directly from the NPMI matrices,
# classifying every bin as belonging to a lost-dominated block, a gained-
# dominated block, or neither, then computing mean NPMI for each CLASS PAIR
# (lost-lost, lost-gained, gained-gained, ...) per distance band, per genotype.
#
# What the answers would look like:
#   - redistribution: lost-lost falls and gained-gained rises by similar amounts,
#     with the total roughly conserved.
#   - genuine local loss: lost-lost falls without a matching rise elsewhere.
#   - compartment switching: the OFF-DIAGONAL (lost-gained) changes most.
#
# One chromosome at a time; matrices are never held in pairs.
# ================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })

for (f in c("gam_bins_annotated.rds")) if (!file.exists(file.path(rds_dir, f)))
  stop("missing ", f, " - run GAM_1 first", call. = FALSE)
blk_f <- file.path(tables_dir, "gam6_large_blocks.tsv")
if (!file.exists(blk_f)) stop("missing gam6_large_blocks.tsv - run GAM_6 first", call. = FALSE)

x <- readRDS(file.path(rds_dir, "gam_bins_annotated.rds"))
a <- as.data.table(x$annot); setorder(a, chr, start)
blk <- fread(blk_f)

# class every bin by the block it falls in
a[, cls := "neither"]
for (i in seq_len(nrow(blk)))
  a[chr == blk$chr[i] & start >= blk$start[i] & end <= blk$end[i], cls := blk$dirn[i]]
cat("=== bins per class ===\n"); print(table(a$cls))

BAND_EDGES <- c(0, 1e6, 5e6, 1e7, Inf)
BAND_NAMES <- c("0.02-1Mb", "1-5Mb", "5-10Mb", ">10Mb")
CLS <- c("lost", "gained", "neither")

acc <- list()
for (ch in gam_autosomes) {
  ids <- a[chr == ch]
  if (!nrow(ids)) next
  for (g in c("WT", "HP1gKO")) {
    m <- tryCatch(load_npmi(g, ch), error = function(e) NULL)
    if (is.null(m)) next
    # align the matrix to our bin order by bin id
    keep <- intersect(rownames(m), ids$bin_id)
    m <- m[keep, keep, drop = FALSE]
    cl <- ids$cls[match(keep, ids$bin_id)]
    pos <- ids$start[match(keep, ids$bin_id)]
    n <- nrow(m)
    if (n < 50) { rm(m); gc(FALSE); next }
    d <- abs(outer(pos, pos, "-"))
    bandi <- cut(d, BAND_EDGES, labels = BAND_NAMES, right = FALSE)
    ci <- match(cl, CLS)
    pair <- outer(ci, ci, function(p, q) paste(CLS[pmin(p, q)], CLS[pmax(p, q)], sep = "-"))
    ok <- is.finite(m) & upper.tri(m) & !is.na(bandi)
    if (!any(ok)) { rm(m, d, bandi, pair); gc(FALSE); next }
    acc[[length(acc) + 1]] <- data.table(
      chr = ch, genotype = g, pair = pair[ok], band = as.character(bandi[ok]), v = m[ok]
    )[, .(sum_npmi = sum(v), n = .N), by = .(chr, genotype, pair, band)]
    rm(m, d, bandi, pair, ok); gc(FALSE)
  }
  cat("  ", ch, "done\n")
}
A <- rbindlist(acc); rm(acc); gc(FALSE)

# pool chromosomes exactly (sums and counts, never a mean of means)
P <- A[, .(mean_npmi = sum(sum_npmi) / sum(n), n = sum(n)), by = .(genotype, pair, band)]
W <- dcast(P, pair + band ~ genotype, value.var = c("mean_npmi", "n"))
setnames(W, c("mean_npmi_WT", "mean_npmi_HP1gKO"), c("WT", "HP1gKO"))
W[, delta := HP1gKO - WT]
W[, pct := 100 * delta / abs(WT)]
setorder(W, band, pair)
fwrite(W, file.path(tables_dir, "gam8_block_contact_budget.tsv"), sep = "\t")

cat("\n=== mean NPMI by block-class pair and distance, WT vs KO ===\n")
for (bn in BAND_NAMES) {
  cat("\n--", bn, "--\n")
  d <- W[band == bn]
  if (!nrow(d)) next
  cat(sprintf("  %-18s %8s %8s %9s %8s\n", "class pair", "WT", "KO", "delta", "%"))
  for (i in seq_len(nrow(d)))
    cat(sprintf("  %-18s %8.4f %8.4f %+9.4f %+7.1f%%\n",
        d$pair[i], d$WT[i], d$HP1gKO[i], d$delta[i], d$pct[i]))
}

cat("\n=== verdict shape ===\n")
ll <- W[pair == "lost-lost"]; gg <- W[pair == "gained-gained"]
for (bn in BAND_NAMES) {
  l <- ll[band == bn]; g <- gg[band == bn]
  if (!nrow(l) || !nrow(g)) next
  cat(sprintf("  %-9s lost-lost %+.1f%%  gained-gained %+.1f%%  -> %s\n", bn, l$pct, g$pct,
    if (l$pct < 0 && g$pct > 0) "RECIPROCAL (redistribution-shaped)"
    else if (l$pct < 0 && g$pct < 0) "both fall (net loss)"
    else if (l$pct > 0 && g$pct > 0) "both rise" else "mixed"))
}

p <- ggplot(W, aes(factor(band, BAND_NAMES), pct, fill = pair)) +
  geom_col(position = position_dodge(0.85), width = 0.78) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_disc(name = NULL) +
  labs(title = "Contact budget change by block class and distance",
       subtitle = "Percent change in mean NPMI, HP1gKO vs WT. Reciprocal movement between\nlost-lost and gained-gained would indicate redistribution rather than local loss.",
       x = "contact distance", y = "change in mean NPMI (%)") +
  theme_m + theme(legend.position = "top", text = element_text(size = 14))
save_fig(p, "block_contact_budget", subdir = "blocks", width = 11, height = 7)
cat("\nGAM_8 complete.\n")
