# ================================================================
# GAM_9 - The protocadherin cluster at native resolution
# ----------------------------------------------------------------
#     cd GAM && Rscript GAM_9_Pcdh_Fine.R        (after GAM_1, GAM_7)
#
# WHY SEPARATELY. GAM_5 found protocadherin to be the only target set with a
# decompaction signature, but it tested 29 bins with genome-wide parameters and
# half its verdict leaned on chr18's compartment axis, which is not usable
# (WT |rho| 0.107). The cluster deserves its own look at the native 20 kb
# resolution, using the contact PAIRS rather than the 0-5 Mb Hotspots summary,
# and without inheriting any genome-scale smoothing.
#
# It is also the one place where fine resolution is clearly warranted: the
# Pcdha/b/g clusters span ~1 Mb, which the 5 Mb delivered window cannot resolve
# at all.
# ================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(data.table)
  library(TxDb.Mmusculus.UCSC.mm39.knownGene); library(org.Mm.eg.db) })

x <- readRDS(file.path(rds_dir, "gam_bins_annotated.rds"))
a <- as.data.table(x$annot); setorder(a, chr, start)
for (mk in colnames(x$l2fc)) a[[paste0("l2fc_", mk)]] <- x$l2fc[match(a$bin_id, rownames(x$l2fc)), mk]

# --- locate the cluster from the gene models, not from a hard-coded interval --
txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
gn <- genes(txdb, single.strand.genes.only = TRUE)
sym <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, keys = names(gn),
         column = "SYMBOL", keytype = "ENTREZID"))
pc <- gn[which(grepl("^Pcdh(a|b|g)", ifelse(is.na(sym), "", sym)))]
pc <- pc[as.character(seqnames(pc)) == "chr18"]
if (!length(pc)) stop("no chr18 Pcdh genes found", call. = FALSE)
lo <- min(start(pc)); hi <- max(end(pc))
PAD <- 2e6
w0 <- max(0, lo - PAD); w1 <- hi + PAD
cat(sprintf("=== Pcdh cluster chr18 %.2f-%.2f Mb (%d genes); window %.2f-%.2f Mb ===\n",
            lo/1e6, hi/1e6, length(pc), w0/1e6, w1/1e6))

inw  <- a[chr == "chr18" & start >= w0 & end <= w1]
incl <- a[chr == "chr18" & end > lo & start < hi]
cat(sprintf("  cluster bins: %d   window bins: %d\n", nrow(incl), nrow(inw)))

# --- ChIP marks over the cluster vs the rest of chr18 and the genome ---------
cat("\n=== ChIP log2FC: cluster vs chr18 vs genome ===\n")
cat(sprintf("  %-10s %9s %9s %9s\n", "mark", "cluster", "chr18", "genome"))
for (mk in colnames(x$l2fc)) {
  v <- paste0("l2fc_", mk)
  cat(sprintf("  %-10s %+9.3f %+9.3f %+9.3f\n", mk,
      median(incl[[v]], na.rm = TRUE),
      median(a[chr == "chr18"][[v]], na.rm = TRUE),
      median(a[[v]], na.rm = TRUE)))
}

# --- local contacts from the pairs, at native resolution, by distance --------
g7 <- file.path(rds_dir, "gam7_local_contacts.rds")
if (file.exists(g7)) {
  b7 <- as.data.table(readRDS(g7)$bins)
  cols <- grep("^(lost|gained)_", names(b7), value = TRUE)
  m <- merge(a[chr == "chr18", .(bin_id, start, end)],
             b7[, c("bin_id", cols, "local_lost", "local_gained"), with = FALSE], by = "bin_id")
  m[, inclust := start >= lo & end <= hi]
  cat("\n=== local differential contacts, cluster vs rest of chr18 (per bin) ===\n")
  cat(sprintf("  %-14s %10s %10s %8s\n", "band", "cluster", "chr18 rest", "ratio"))
  for (cc in cols) {
    a1 <- mean(m[inclust == TRUE][[cc]]); a2 <- mean(m[inclust == FALSE][[cc]])
    cat(sprintf("  %-14s %10.1f %10.1f %8.2f\n", cc, a1, a2, a1 / max(a2, 1e-9)))
  }
} else message("GAM_7 output absent - skipping the pairs-based cluster counts")

# --- native-resolution NPMI across the window, WT vs KO ----------------------
mats <- list()
for (g in c("WT", "HP1gKO")) {
  mm <- load_npmi(g, "chr18")
  ids <- intersect(inw$bin_id, rownames(mm))
  mats[[g]] <- mm[ids, ids, drop = FALSE]
  rm(mm); gc(FALSE)
}
ids <- intersect(rownames(mats$WT), rownames(mats$HP1gKO))
W <- mats$WT[ids, ids]; K <- mats$HP1gKO[ids, ids]
# bin ids are "chr18:36000000-36020000": take the START, not the whole range.
# sub(".*:", "") alone leaves "36000000-36020000", which coerces to NA and makes
# every downstream mean NaN.
pos <- as.numeric(sub("-.*", "", sub(".*:", "", ids)))
stopifnot(all(is.finite(pos)))
D <- K - W
cat(sprintf("\n=== native 20 kb NPMI over the window (%d x %d bins) ===\n", length(ids), length(ids)))
inc <- pos >= lo & pos <= hi
sub_mean <- function(M, r, c) mean(M[r, c][upper.tri(M[r, c])], na.rm = TRUE)
cat(sprintf("  within-cluster   NPMI  WT %.4f  KO %.4f  delta %+.4f (%+.1f%%)\n",
    sub_mean(W, inc, inc), sub_mean(K, inc, inc),
    sub_mean(K, inc, inc) - sub_mean(W, inc, inc),
    100 * (sub_mean(K, inc, inc) - sub_mean(W, inc, inc)) / abs(sub_mean(W, inc, inc))))
fl <- !inc
cat(sprintf("  cluster-to-flank NPMI  WT %.4f  KO %.4f  delta %+.4f (%+.1f%%)\n",
    mean(W[inc, fl], na.rm = TRUE), mean(K[inc, fl], na.rm = TRUE),
    mean(K[inc, fl], na.rm = TRUE) - mean(W[inc, fl], na.rm = TRUE),
    100 * (mean(K[inc, fl], na.rm = TRUE) - mean(W[inc, fl], na.rm = TRUE)) / abs(mean(W[inc, fl], na.rm = TRUE))))
cat(sprintf("  within-flank     NPMI  WT %.4f  KO %.4f  delta %+.4f (%+.1f%%)\n",
    sub_mean(W, fl, fl), sub_mean(K, fl, fl),
    sub_mean(K, fl, fl) - sub_mean(W, fl, fl),
    100 * (sub_mean(K, fl, fl) - sub_mean(W, fl, fl)) / abs(sub_mean(W, fl, fl))))

dt <- CJ(i = seq_along(ids), j = seq_along(ids))[i < j]
dt[, `:=`(x = pos[i] / 1e6, y = pos[j] / 1e6, d = D[cbind(i, j)])]
lim <- quantile(abs(dt$d), 0.98, na.rm = TRUE)
p <- ggplot(dt[is.finite(d)], aes(x, y, fill = pmax(-lim, pmin(lim, d)))) +
  geom_raster() +
  annotate("rect", xmin = lo/1e6, xmax = hi/1e6, ymin = lo/1e6, ymax = hi/1e6,
           fill = NA, colour = "black", linewidth = 0.5) +
  scale_fill_gradient2(low = heat_low, mid = heat_mid, high = heat_high,
                       midpoint = 0, name = "KO - WT\nNPMI") +
  coord_fixed() +
  labs(title = "Protocadherin cluster, chr18 - native 20 kb NPMI difference",
       subtitle = sprintf("KO minus WT. Box = Pcdha/b/g cluster (%.2f-%.2f Mb). Red = more contact in KO.",
                          lo/1e6, hi/1e6),
       x = "chr18 (Mb)", y = "chr18 (Mb)") +
  theme_m + theme(text = element_text(size = 14))
save_fig(p, "pcdh_native_npmi_difference", subdir = "pcdh", width = 9, height = 8)

fwrite(data.table(bin_id = ids, pos = pos, in_cluster = inc,
                  npmi_row_WT = rowMeans(W, na.rm = TRUE),
                  npmi_row_KO = rowMeans(K, na.rm = TRUE)),
       file.path(tables_dir, "gam9_pcdh_native_profile.tsv"), sep = "\t")
cat("\nGAM_9 complete.\n")
