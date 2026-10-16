# ================================================================
# GAM_11 - Do domains that LOSE the mark decouple from each other?
# ----------------------------------------------------------------
#     cd GAM && Rscript GAM_11_Mark_Change_Contacts.R      (after GAM_1)
#
# THREE VARIABLES, which no earlier stage combined:
#   1. DOMAIN IDENTITY - which bins are heterochromatin, from the ChIP-defined
#      merged H3K9me3/H4K20me3 peak set.
#   2. CONTACTS - NPMI between a PAIR of bins, so a contact is classified by
#      BOTH of its endpoints rather than attributed to one bin.
#   3. WHETHER THE MARK CHANGED - the per-bin H4K20me3 log2FC.
#
# GAM_10 used (1) and (2) only. Pooling all heterochromatin mixes domains that
# lost the mark with domains that kept it, which would average away exactly the
# effect in question. This stage asks the causal version: do domains that LOSE
# H4K20me3 lose contacts WITH EACH OTHER, relative to domains that keep it?
#
# WHY TERTILES, NOT AN ABSOLUTE CUTOFF. H4K20me3 falls almost genome-wide
# (median per-bin log2FC -0.497, ~100% of bins down), so "lost" versus "kept" is
# only meaningful RELATIVE to that shift - the same reason the ChIP side reports
# median-centred compartments (see ../ChIP/results/COMPARTMENT_CUTOFFS.md). The
# split is therefore into tertiles computed WITHIN heterochromatin: het_lost is
# the third that falls most, het_kept the third that falls least. Both are
# heterochromatin; they differ only in how much mark they lost.
#
# STATISTICS - as GAM_10:
#   * every comparison is WITHIN a distance band (NPMI falls ~6x from 20 kb to
#     20 Mb, and the classes have different distance distributions);
#   * per-chromosome means first, then a bootstrap over the 19 autosomes - the
#     honest n is 19, not the number of pairs;
#   * difference-in-differences, so each class is its own control across
#     genotype and the baseline differences between strata cancel.
#
# THE CONFOUNDER TO WATCH, reported rather than assumed away: bins that lose
# more H4K20me3 are already more compact and more B-compartment in WT
# (GAM_3: pc1 +0.196, npmi_short -0.166 after coverage control). So het_lost and
# het_kept are not interchangeable at baseline. The DiD removes a constant
# baseline offset but not a difference in how the two respond; the WT baselines
# are printed so the reader can judge.
# ================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(GenomicRanges) })

gam11_het_cut <- 0.5
gam11_nboot   <- 2000L

g1 <- file.path(rds_dir, "gam_bins_annotated.rds")
if (!file.exists(g1)) stop("GAM_1 output not found:\n  ", g1,
  "\nRun it first:  cd GAM && Rscript GAM_1_Bins_and_Signal.R", call. = FALSE)
x <- readRDS(g1)
a <- as.data.table(x$annot); setorder(a, chr, start)
a[, h4 := x$l2fc[match(bin_id, rownames(x$l2fc)), "H4K20me3"]]
a[, k9 := x$l2fc[match(bin_id, rownames(x$l2fc)), "H3K9me3"]]

# --- variable 1: heterochromatin identity from the ChIP domains --------------
dm <- fread(het_domain_bed, header = FALSE, select = 1:3, col.names = c("chr","start","end"))
dm <- dm[grepl("^([0-9]+|chr[0-9]+)$", chr)]
dg <- to_ucsc(GRanges(ifelse(grepl("^chr", dm$chr), dm$chr, paste0("chr", dm$chr)),
                      IRanges(dm$start + 1, dm$end)))
bins <- GRanges(a$chr, IRanges(a$start + 1, a$end))
red <- reduce(dg, ignore.strand = TRUE)
h <- findOverlaps(bins, red, ignore.strand = TRUE)
w <- width(pintersect(bins[queryHits(h)], red[subjectHits(h)], ignore.strand = TRUE))
v <- numeric(length(bins)); agg <- tapply(w, queryHits(h), sum)
v[as.integer(names(agg))] <- as.numeric(agg)
a[, het_frac := pmin(v / width(bins), 1)]
a[, is_het := het_frac >= gam11_het_cut]
stopifnot(sum(a$is_het) > 1000)

# --- variable 3: did the mark change? tertiles WITHIN heterochromatin --------
ht <- a[is_het == TRUE & is.finite(h4)]
qs <- quantile(ht$h4, c(1/3, 2/3))
a[, stratum := "eu"]
a[is_het == TRUE & is.finite(h4) & h4 <= qs[1], stratum := "het_lost"]
a[is_het == TRUE & is.finite(h4) & h4 >  qs[1] & h4 <= qs[2], stratum := "het_mid"]
a[is_het == TRUE & is.finite(h4) & h4 >  qs[2], stratum := "het_kept"]

cat("=== strata (tertiles of H4K20me3 log2FC within heterochromatin) ===\n")
s <- a[stratum != "eu", .(bins = .N, median_H4K20me3 = round(median(h4, na.rm = TRUE), 3),
      median_H3K9me3 = round(median(k9, na.rm = TRUE), 3),
      median_het_frac = round(median(het_frac), 2)), by = stratum]
print(s[order(median_H4K20me3)], row.names = FALSE)
cat(sprintf("  tertile cuts at log2FC %.3f and %.3f; euchromatic bins: %d\n",
            qs[1], qs[2], sum(a$stratum == "eu")))

BAND_EDGES <- c(0, 2e5, 1e6, 5e6, 1e7, Inf)
BAND_NAMES <- c("<0.2Mb", "0.2-1Mb", "1-5Mb", "5-10Mb", ">10Mb")
KEEP <- c("het_lost-het_lost", "het_lost-het_kept", "het_kept-het_kept",
          "het_lost-eu", "het_kept-eu", "eu-eu")

acc <- list()
for (ch in gam_autosomes) {
  ids <- a[chr == ch]; if (nrow(ids) < 50) next
  for (g in c("WT", "HP1gKO")) {
    m <- tryCatch(load_npmi(g, ch), error = function(e) NULL)
    if (is.null(m)) next
    keep <- intersect(rownames(m), ids$bin_id)
    m <- m[keep, keep, drop = FALSE]
    k <- match(keep, ids$bin_id)
    st <- ids$stratum[k]; pos <- ids$start[k]; n <- length(keep)
    band <- cut(abs(outer(pos, pos, "-")), BAND_EDGES, labels = BAND_NAMES, right = FALSE)
    # unordered class label for the pair
    cls <- outer(st, st, function(p, q) paste(pmin(p, q), pmax(p, q), sep = "-"))
    ok <- is.finite(m) & upper.tri(m) & !is.na(band) & (cls %in% KEEP)
    if (!any(ok)) { rm(m, band, cls); gc(FALSE); next }
    acc[[length(acc) + 1]] <- data.table(chr = ch, genotype = g, cls = cls[ok],
      band = as.character(band[ok]), v = m[ok])[, .(s = sum(v), n = .N),
      by = .(chr, genotype, cls, band)]
    rm(m, band, cls, ok); gc(FALSE)
  }
  cat("  ", ch, "done\n")
}
A <- rbindlist(acc); rm(acc); gc(FALSE)

pool <- function(d) d[, .(mean_npmi = sum(s) / sum(n), n = sum(n)), by = .(genotype, cls, band)]
P <- pool(A)
W <- dcast(P, cls + band ~ genotype, value.var = c("mean_npmi", "n"))
setnames(W, c("mean_npmi_WT","mean_npmi_HP1gKO"), c("WT","KO"))
W[, pct := 100 * (KO - WT) / WT]

chrs <- unique(A$chr); set.seed(42)
boot <- rbindlist(lapply(seq_len(gam11_nboot), function(i) {
  smp <- sample(chrs, length(chrs), replace = TRUE)
  d <- rbindlist(lapply(seq_along(smp), function(j) { z <- A[chr == smp[j]]; z[, chr := j]; z }))
  wj <- dcast(pool(d), cls + band ~ genotype, value.var = "mean_npmi")
  wj[, .(cls, band, pct = 100 * (HP1gKO - WT) / WT, rep = i)]
}))
ci <- boot[, .(lo = quantile(pct, .025, na.rm = TRUE), hi = quantile(pct, .975, na.rm = TRUE)),
           by = .(cls, band)]
W <- merge(W, ci, by = c("cls","band"))
W[, band := factor(band, BAND_NAMES)]; setorder(W, band, cls)
fwrite(W, file.path(tables_dir, "gam11_stratum_pair_change.tsv"), sep = "\t")

cat("\n=== change in mean NPMI by stratum pair, within distance band ===\n")
cat("   95% CI from a bootstrap over the 19 autosomes\n")
for (bn in BAND_NAMES) {
  d <- W[band == bn]; if (!nrow(d)) next
  cat("\n--", bn, "--\n")
  cat(sprintf("  %-20s %8s %8s %9s %18s %13s\n","stratum pair","WT","KO","change","95% CI","pairs"))
  for (i in seq_len(nrow(d)))
    cat(sprintf("  %-20s %8.4f %8.4f %+8.1f%% [%+6.1f%%,%+6.1f%%] %13s\n",
        d$cls[i], d$WT[i], d$KO[i], d$pct[i], d$lo[i], d$hi[i],
        format(d$n_WT[i], big.mark = ",")))
}

# --- THE CONTRAST: lost-lost vs kept-kept -----------------------------------
dd <- dcast(boot[cls %in% c("het_lost-het_lost","het_kept-het_kept")],
            band + rep ~ cls, value.var = "pct")
setnames(dd, c("het_lost-het_lost","het_kept-het_kept"), c("LL","KK"))
dd[, diff := LL - KK]
D <- dd[, .(est = mean(diff, na.rm = TRUE),
            lo = quantile(diff, .025, na.rm = TRUE),
            hi = quantile(diff, .975, na.rm = TRUE),
            p = 2 * min(mean(diff <= 0, na.rm = TRUE), mean(diff >= 0, na.rm = TRUE))), by = band]
D <- merge(D, W[cls == "het_lost-het_lost", .(band, LL = pct)], by = "band")
D <- merge(D, W[cls == "het_kept-het_kept", .(band, KK = pct)], by = "band")
D[, band := factor(band, BAND_NAMES)]; setorder(D, band)
fwrite(D, file.path(tables_dir, "gam11_lost_vs_kept.tsv"), sep = "\t")

cat("\n=== DO MARK-LOSING DOMAINS DECOUPLE FROM EACH OTHER? ===\n")
cat("   (lost-lost change) minus (kept-kept change), same distance band\n")
cat(sprintf("  %-10s %11s %11s %12s %18s %7s\n","band","lost-lost","kept-kept","difference","95% CI","p"))
for (i in seq_len(nrow(D)))
  cat(sprintf("  %-10s %+10.1f%% %+10.1f%% %+11.1f%% [%+6.1f%%,%+6.1f%%] %7.3f\n",
      as.character(D$band[i]), D$LL[i], D$KK[i], D$est[i], D$lo[i], D$hi[i], D$p[i]))
cat("\n  Negative difference = mark-losing heterochromatin loses contact with other\n")
cat("  mark-losing heterochromatin MORE than mark-retaining heterochromatin does.\n")
cat("  A CI excluding 0 is the result; p is a bootstrap proportion, not parametric.\n")

cat("\n=== baseline check (the strata are NOT interchangeable in WT) ===\n")
b0 <- W[cls %in% c("het_lost-het_lost","het_kept-het_kept"), .(cls, band, WT = round(WT, 4))]
print(dcast(b0, band ~ cls, value.var = "WT"), row.names = FALSE)

pw <- W[cls %in% KEEP]
p <- ggplot(pw, aes(band, pct, fill = cls)) +
  geom_col(position = position_dodge(.85), width = .78) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(.85),
                width = .2, linewidth = .35) +
  geom_hline(yintercept = 0, linewidth = .3) +
  scale_fill_disc(name = NULL) +
  labs(title = "Contact change by what happened to the mark at both ends",
       subtitle = paste0("Heterochromatin split into thirds by H4K20me3 log2FC. Change in mean NPMI,\n",
                         "HP1gKO vs WT, within distance band. Bars are 95% CIs over the 19 autosomes."),
       x = "contact distance", y = "change in mean NPMI (%)") +
  theme_m + theme(legend.position = "top", text = element_text(size = 13))
save_fig(p, "mark_change_contact_classes", subdir = "heterochromatin", width = 12, height = 7)
cat("\nGAM_11 complete.\n")
