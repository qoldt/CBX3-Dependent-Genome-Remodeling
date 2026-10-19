# GAM_8_Block_Classes.R ------------------------------------------------------
# What chromatin are the contacts inside the large one-directional blocks?
#
# GAM_4 writes gam8_block_contact_budget.tsv, which shows that pairs with BOTH
# ends in a lost block fall much harder than background and both ends in a
# gained block rise. That table says nothing about WHICH chromatin those
# contacts are in. This stage adds the level-2 chromatin class of both ends to
# the same within-block pairs, then bootstraps over chromosomes.
#
# Needs GAM_1 (annotation) and GAM_4 (gam6_large_blocks.tsv). Heavy - it reads
# every autosomal NPMI matrix for both genotypes, one chromosome at a time, so
# run it under the safe-local-alignment skill.
#
# Writes gam16_block_contact_classes{,_perchrom}.tsv, gam16_block_class_ci.tsv
# and gam16_block_class_contrast_ci.tsv. Read with results/LARGE_BLOCKS.md.
#
# The gam16_ prefix follows the repo's legacy table numbering, which does NOT
# track the stage number - see GAM/CLAUDE.md.
source("config.R")
library(data.table)
setDTthreads(2)

ann <- readRDS(file.path(rds_dir, "gam_bins_annotated.rds"))
a   <- as.data.table(ann$annot)
cc  <- gam_chromatin_classes(a)
a   <- merge(a, cc, by = "bin_id", sort = FALSE)
setorder(a, chr, start)

blk <- fread(file.path(tables_dir, "gam6_large_blocks.tsv"))
cat(sprintf("blocks: %d lost / %d gained over %d chromosomes\n",
            sum(blk$dirn == "lost"), sum(blk$dirn == "gained"),
            uniqueN(blk$chr)))

BAND_EDGES <- gam_band_edges
BAND_NAMES <- gam_band_names
CL <- gam_class2_levels

acc <- list()
for (ch in sort(unique(blk$chr))) {
  B <- blk[chr == ch]
  A <- a[chr == ch]
  if (!nrow(A)) next
  # block id per bin: 0 = not in a block. Blocks do not overlap (they are runs).
  A[, blkid := 0L][, bdir := NA_character_]
  for (i in seq_len(nrow(B))) {
    A[start >= B$start[i] & start < B$end[i], `:=`(blkid = i, bdir = B$dirn[i])]
  }
  if (!any(A$blkid > 0L)) next

  k2  <- match(A$class2, CL)          # NA where unclassified
  pos <- A$start
  bid <- A$blkid
  bdr <- A$bdir

  for (g in c("WT", "HP1gKO")) {
    m <- load_npmi(g, ch)
    idx <- match(A$bin_id, rownames(m))
    keep <- !is.na(idx)
    m <- m[idx[keep], idx[keep], drop = FALSE]
    k2k <- k2[keep]; posk <- pos[keep]; bidk <- bid[keep]; bdrk <- bdr[keep]

    iu <- which(upper.tri(m), arr.ind = TRUE)
    i1 <- iu[, 1]; i2 <- iu[, 2]; rm(iu)
    val <- m[cbind(i1, i2)]
    rm(m); gc()
    # BOTH ends in the SAME block - the same definition the budget table uses
    # for lost-lost / gained-gained, so the two are directly comparable.
    same <- bidk[i1] > 0L & bidk[i1] == bidk[i2]
    ok <- same & is.finite(val) & !is.na(k2k[i1]) & !is.na(k2k[i2])
    if (!any(ok)) { rm(val, i1, i2); gc(); next }
    i1 <- i1[ok]; i2 <- i2[ok]; val <- val[ok]
    band <- cut(abs(posk[i2] - posk[i1]), BAND_EDGES, labels = BAND_NAMES, right = FALSE)
    lo <- pmin(k2k[i1], k2k[i2]); hi <- pmax(k2k[i1], k2k[i2])
    pair <- paste(CL[lo], CL[hi], sep = "-")
    d <- data.table(dirn = bdrk[i1], pair = pair, band = band, v = val)
    d <- d[!is.na(band)]
    acc[[length(acc) + 1]] <- d[, .(s = sum(v), n = .N), by = .(dirn, pair, band)][
      , genotype := g][, chr := ch]
    rm(d, val, i1, i2, band, pair, lo, hi); gc()
  }
  cat(sprintf("  %s done (%d block bins)\n", ch, sum(A$blkid > 0L)))
}

R <- rbindlist(acc)
fwrite(R, file.path(tables_dir, "gam16_block_contact_classes_perchrom.tsv"), sep = "\t")

P <- R[, .(s = sum(s), n = sum(n)), by = .(dirn, pair, band, genotype)]
W <- dcast(P, dirn + pair + band ~ genotype, value.var = c("s", "n"))
W[, `:=`(mean_WT = s_WT / n_WT, mean_HP1gKO = s_HP1gKO / n_HP1gKO)]
W[, pct := 100 * (mean_HP1gKO / mean_WT - 1)]
W[, band := factor(band, BAND_NAMES)]
setorder(W, dirn, band, -n_WT)
fwrite(W[, .(dirn, pair, band, mean_WT, mean_HP1gKO, n_WT, n_HP1gKO, pct)],
       file.path(tables_dir, "gam16_block_contact_classes.tsv"), sep = "\t")
cat("\nwrote gam16_block_contact_classes.tsv\n")

# --- bootstrap over chromosomes ---------------------------------------------
# The chromosome is the only resampling unit the design supports (2 vs 2 merged
# replicates, no per-contact variance). 2000 resamples, percentile interval.
R <- fread(file.path(tables_dir, "gam16_block_contact_classes_perchrom.tsv"))
set.seed(42)
CHR <- sort(unique(R$chr)); B <- 2000
pctf <- function(d) {
  w <- dcast(d[, .(s = sum(s), n = sum(n)), by = .(pair, band, genotype)],
             pair + band ~ genotype, value.var = c("s", "n"))
  w[, 100 * ((s_HP1gKO / n_HP1gKO) / (s_WT / n_WT) - 1)]
}
out <- list()
for (dd in c("lost", "gained")) {
  D <- R[dirn == dd]
  for (bd in gam_band_names) {
    for (pr in c("H3K9me3-H3K9me3", "Quiescent-Quiescent")) {
      sub <- D[pair == pr & band == bd]
      if (!nrow(sub) || sum(sub$n) < 2000) next
      bs <- vapply(seq_len(B), function(i) {
        cs <- sample(CHR, length(CHR), TRUE)
        x <- rbindlist(lapply(cs, function(c) sub[chr == c]))
        if (!nrow(x)) return(NA_real_)
        g <- x[, .(s = sum(s), n = sum(n)), by = genotype]
        if (nrow(g) < 2) return(NA_real_)
        100 * (g[genotype == "HP1gKO", s/n] / g[genotype == "WT", s/n] - 1)
      }, numeric(1))
      g <- sub[, .(s = sum(s), n = sum(n)), by = genotype]
      out[[length(out)+1]] <- data.table(
        dirn = dd, pair = pr, band = bd,
        pct = 100 * (g[genotype=="HP1gKO", s/n] / g[genotype=="WT", s/n] - 1),
        lo = quantile(bs, .025, na.rm = TRUE), hi = quantile(bs, .975, na.rm = TRUE),
        n_chr = uniqueN(sub$chr), n_pairs = sum(sub[genotype=="WT"]$n))
    }
  }
}
O <- rbindlist(out)
O[, band := factor(band, gam_band_names)]
setorder(O, dirn, band, pair)
fwrite(O, file.path(tables_dir, "gam16_block_class_ci.tsv"), sep = "\t")
print(O[, .(dirn, pair, band, pct = round(pct,2), lo = round(lo,2), hi = round(hi,2), n_chr)])

# --- the H3K9me3 vs quiescent contrast, on ONE shared chromosome draw -------
R <- fread(file.path(tables_dir, "gam16_block_contact_classes_perchrom.tsv"))
setkey(R, dirn, pair, band, chr, genotype)
set.seed(42); CHR <- sort(unique(R$chr)); B <- 2000
pc <- function(x) { g <- x[, .(s = sum(s), n = sum(n)), by = genotype]
  if (nrow(g) < 2) return(NA_real_)
  100 * (g[genotype == "HP1gKO", s/n] / g[genotype == "WT", s/n] - 1) }
out <- list()
for (dd in c("lost","gained")) for (bd in gam_band_names) {
  A <- R[dirn == dd & pair == "H3K9me3-H3K9me3" & band == bd]
  Q <- R[dirn == dd & pair == "Quiescent-Quiescent" & band == bd]
  if (!nrow(A) || !nrow(Q) || sum(A$n) < 5000) next
  obs <- pc(A) - pc(Q)
  # The SAME chromosome draw for both classes, so a chromosome-level offset
  # cancels inside each resample exactly as it does in the point estimate.
  bs <- vapply(seq_len(B), function(i) {
    cs <- sample(CHR, length(CHR), TRUE)
    pc(A[.(dd, "H3K9me3-H3K9me3", bd, cs)]) - pc(Q[.(dd, "Quiescent-Quiescent", bd, cs)])
  }, numeric(1))
  out[[length(out)+1]] <- data.table(dirn = dd, band = bd, contrast = obs,
    lo = quantile(bs,.025,na.rm=TRUE), hi = quantile(bs,.975,na.rm=TRUE),
    n_chr = uniqueN(A$chr))
}
O <- rbindlist(out); O[, band := factor(band, gam_band_names)]; setorder(O, dirn, band)
fwrite(O, file.path(tables_dir, "gam16_block_class_contrast_ci.tsv"), sep="\t")
print(O[, .(dirn, band, contrast = round(contrast,2), lo = round(lo,2), hi = round(hi,2), n_chr)])

cat("
GAM_8 complete.
")
