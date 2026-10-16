# ================================================================
# GAM_10 - Do heterochromatin domains lose contacts WITHIN themselves,
#          or WITH EACH OTHER?
# ----------------------------------------------------------------
#     cd GAM && Rscript GAM_10_Heterochromatin_Contacts.R     (after GAM_1)
#
# THE QUESTION. Every earlier test counted, per bin, the differential contacts
# that bin participates in, with no regard to WHAT IT CONTACTS. So a
# rearrangement in which heterochromatin decouples from other heterochromatin
# while each region keeps its total contact count is invisible to all of them -
# the counts are unchanged, only the partners move. This stage classifies a
# contact by BOTH of its endpoints.
#
# THE DOMAINS. data/peaks/2000bp_merge_H3K9me3_H4K20me3.bed, the ChIP-defined
# shared H3K9me3/H4K20me3 peak set: 126,197 domains, 444 Mb, 16.3% of the genome.
# They are mostly SMALL - median 1,210 bp - so "a contact within one domain" is
# the matrix diagonal for 99.5% of them and cannot be measured at 20 kb. The
# within-domain question is therefore asked only of domains >= gam10_min_domain,
# where an intra-domain contact spans at least two bins. The between-domain
# question uses every bin, scored by heterochromatin coverage.
#
# CLASSES (a contact is assigned by its two endpoint bins):
#   intra_dom  both bins inside the SAME large domain
#   inter_het  both bins heterochromatic, NOT the same large domain
#   het_eu     exactly one bin heterochromatic
#   eu_eu      neither
#
# STATISTICS - three things make or break this:
#  1. DISTANCE IS THE CONFOUNDER. Intra-domain contacts are short-range by
#     construction and inter-domain ones are long, while NPMI falls from ~0.69 at
#     20 kb to ~0.12 at 20 Mb. Comparing classes without matching distance
#     measures the distance distribution, not the biology. Everything below is
#     computed WITHIN a distance band.
#  2. THE HONEST n IS 19, NOT 25 MILLION. Each bin appears in thousands of pairs,
#     so contacts are nowhere near independent. Per-chromosome means are formed
#     first and the uncertainty comes from a bootstrap over the 19 autosomes.
#     A test over individual contacts would report absurd significance.
#  3. DIFFERENCE-IN-DIFFERENCES. Heterochromatic and euchromatic bins differ in
#     baseline NPMI, coverage and mappability, so a between-class comparison at
#     one genotype is confounded. Each class is instead its own control across
#     genotype: the statistic is (KO-WT)/WT within a class, and classes are
#     compared only through the DIFFERENCE of those changes.
# ================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(GenomicRanges) })

gam10_min_domain <- 100000L   # smallest domain in which an intra contact is measurable
gam10_het_cut    <- 0.5       # bin counts as heterochromatic above this coverage
gam10_nboot      <- 2000L

g1 <- file.path(rds_dir, "gam_bins_annotated.rds")
if (!file.exists(g1)) stop("GAM_1 output not found:\n  ", g1,
  "\nRun it first:  cd GAM && Rscript GAM_1_Bins_and_Signal.R", call. = FALSE)
x <- readRDS(g1)
a <- as.data.table(x$annot); setorder(a, chr, start)
a[, h4 := x$l2fc[match(bin_id, rownames(x$l2fc)), "H4K20me3"]]

# --- heterochromatin domains -> per-bin coverage and large-domain membership --
dom_bed <- het_domain_bed   # from config.R
if (!file.exists(dom_bed)) stop("missing ", dom_bed, call. = FALSE)
dm <- fread(dom_bed, header = FALSE, select = 1:3, col.names = c("chr","start","end"))
dm <- dm[grepl("^([0-9]+|chr[0-9]+)$", chr)]
dg <- to_ucsc(GRanges(ifelse(grepl("^chr", dm$chr), dm$chr, paste0("chr", dm$chr)),
                      IRanges(dm$start + 1, dm$end)))
bins <- GRanges(a$chr, IRanges(a$start + 1, a$end))
stopifnot(length(intersect(seqlevels(bins), seqlevels(dg))) > 0)

cov_frac <- {
  red <- reduce(dg, ignore.strand = TRUE)
  h <- findOverlaps(bins, red, ignore.strand = TRUE)
  w <- width(pintersect(bins[queryHits(h)], red[subjectHits(h)], ignore.strand = TRUE))
  v <- numeric(length(bins)); agg <- tapply(w, queryHits(h), sum)
  v[as.integer(names(agg))] <- as.numeric(agg); pmin(v / width(bins), 1)
}
a[, het_frac := cov_frac]
a[, is_het := het_frac >= gam10_het_cut]
stopifnot(sum(a$is_het) > 1000)   # a seqname-style failure would zero this

big <- dg[width(dg) >= gam10_min_domain]
cat(sprintf("=== domains: %d total, %d >= %d kb ===\n", length(dg), length(big),
            gam10_min_domain / 1000))
h2 <- findOverlaps(bins, big, ignore.strand = TRUE)
# a bin belongs to the large domain it overlaps most
wd <- width(pintersect(bins[queryHits(h2)], big[subjectHits(h2)], ignore.strand = TRUE))
dd <- data.table(q = queryHits(h2), s = subjectHits(h2), w = wd)[order(-w)][, .SD[1], by = q]
a[, dom_id := NA_integer_]; a[dd$q, dom_id := dd$s]
cat(sprintf("  bins heterochromatic (cov >= %.2f): %d of %d (%.1f%%)\n",
            gam10_het_cut, sum(a$is_het), nrow(a), 100 * mean(a$is_het)))
cat(sprintf("  bins inside a large domain: %d, across %d domains\n",
            sum(!is.na(a$dom_id)), uniqueN(a$dom_id, na.rm = TRUE)))

BAND_EDGES <- c(0, 2e5, 1e6, 5e6, 1e7, Inf)
BAND_NAMES <- c("<0.2Mb", "0.2-1Mb", "1-5Mb", "5-10Mb", ">10Mb")

# --- accumulate per chromosome x genotype x class x band ---------------------
acc <- list()
for (ch in gam_autosomes) {
  ids <- a[chr == ch]
  if (nrow(ids) < 50) next
  for (g in c("WT", "HP1gKO")) {
    m <- tryCatch(load_npmi(g, ch), error = function(e) NULL)
    if (is.null(m)) next
    keep <- intersect(rownames(m), ids$bin_id)
    m <- m[keep, keep, drop = FALSE]
    k <- match(keep, ids$bin_id)
    het <- ids$is_het[k]; dom <- ids$dom_id[k]; pos <- ids$start[k]
    n <- length(keep)
    d <- abs(outer(pos, pos, "-"))
    band <- cut(d, BAND_EDGES, labels = BAND_NAMES, right = FALSE)
    same_dom <- outer(dom, dom, function(p, q) !is.na(p) & !is.na(q) & p == q)
    both_het <- outer(het, het, "&")
    one_het  <- outer(het, het, "|") & !both_het
    cls <- matrix("eu_eu", n, n)
    cls[both_het] <- "inter_het"
    cls[same_dom] <- "intra_dom"
    cls[one_het]  <- "het_eu"
    ok <- is.finite(m) & upper.tri(m) & !is.na(band)
    if (!any(ok)) { rm(m, d, band, cls, same_dom, both_het, one_het); gc(FALSE); next }
    acc[[length(acc) + 1]] <- data.table(chr = ch, genotype = g, cls = cls[ok],
      band = as.character(band[ok]), v = m[ok])[, .(s = sum(v), n = .N),
      by = .(chr, genotype, cls, band)]
    rm(m, d, band, cls, same_dom, both_het, one_het, ok); gc(FALSE)
  }
  cat("  ", ch, "done\n")
}
A <- rbindlist(acc); rm(acc); gc(FALSE)
fwrite(A, file.path(tables_dir, "gam10_perchrom_sums.tsv"), sep = "\t")

# --- pooled means, percent change, and a chromosome bootstrap ----------------
pool <- function(d) d[, .(mean_npmi = sum(s) / sum(n), n = sum(n)), by = .(genotype, cls, band)]
P <- pool(A)
W <- dcast(P, cls + band ~ genotype, value.var = c("mean_npmi", "n"))
setnames(W, c("mean_npmi_WT","mean_npmi_HP1gKO"), c("WT","KO"))
W[, pct := 100 * (KO - WT) / WT]

chrs <- unique(A$chr)
set.seed(42)
boot <- rbindlist(lapply(seq_len(gam10_nboot), function(i) {
  smp <- sample(chrs, length(chrs), replace = TRUE)
  d <- rbindlist(lapply(seq_along(smp), function(j) { z <- A[chr == smp[j]]; z[, chr := j]; z }))
  p <- pool(d)
  w <- dcast(p, cls + band ~ genotype, value.var = "mean_npmi")
  w[, .(cls, band, pct = 100 * (HP1gKO - WT) / WT, rep = i)]
}))
ci <- boot[, .(lo = quantile(pct, .025, na.rm = TRUE), hi = quantile(pct, .975, na.rm = TRUE)),
           by = .(cls, band)]
W <- merge(W, ci, by = c("cls","band"))
setorder(W, band, cls)
fwrite(W, file.path(tables_dir, "gam10_class_change.tsv"), sep = "\t")

cat("\n=== mean NPMI change by contact class, within distance band ===\n")
cat("   (95% CI from a bootstrap over the 19 autosomes - n is 19, not the pair count)\n")
for (bn in BAND_NAMES) {
  d <- W[band == bn]; if (!nrow(d)) next
  cat("\n--", bn, "--\n")
  cat(sprintf("  %-10s %8s %8s %9s %18s %12s\n","class","WT","KO","change","95% CI","pairs"))
  for (i in seq_len(nrow(d)))
    cat(sprintf("  %-10s %8.4f %8.4f %+8.1f%% [%+6.1f%%,%+6.1f%%] %12s\n",
        d$cls[i], d$WT[i], d$KO[i], d$pct[i], d$lo[i], d$hi[i],
        format(d$n_WT[i], big.mark = ",")))
}

# --- the contrast the question asks for: intra vs inter, same band -----------
did <- boot[cls %in% c("intra_dom","inter_het")]
did <- dcast(did, band + rep ~ cls, value.var = "pct")
did[, diff := intra_dom - inter_het]
D <- did[, .(est = mean(diff, na.rm = TRUE),
             lo = quantile(diff, .025, na.rm = TRUE),
             hi = quantile(diff, .975, na.rm = TRUE),
             p_two_sided = 2 * min(mean(diff <= 0, na.rm = TRUE),
                                   mean(diff >= 0, na.rm = TRUE))), by = band]
D <- merge(D, W[cls == "intra_dom", .(band, intra = pct)], by = "band")
D <- merge(D, W[cls == "inter_het", .(band, inter = pct)], by = "band")
D[, band := factor(band, BAND_NAMES)]; setorder(D, band)
fwrite(D, file.path(tables_dir, "gam10_intra_vs_inter.tsv"), sep = "\t")

cat("\n=== WITHIN-domain vs BETWEEN-heterochromatin, difference-in-differences ===\n")
cat(sprintf("  %-10s %10s %10s %10s %18s %8s\n","band","intra","inter","intra-inter","95% CI","p"))
for (i in seq_len(nrow(D)))
  cat(sprintf("  %-10s %+9.1f%% %+9.1f%% %+10.1f%% [%+6.1f%%,%+6.1f%%] %8.3f\n",
      D$band[i], D$intra[i], D$inter[i], D$est[i], D$lo[i], D$hi[i], D$p_two_sided[i]))
cat("\n  A CI excluding 0 means the two classes changed by different amounts.\n")
cat("  p is the bootstrap two-sided proportion, not a parametric test.\n")

pw <- W[cls %in% c("intra_dom","inter_het","het_eu","eu_eu")]
pw[, band := factor(band, BAND_NAMES)]
pw[, cls := factor(cls, c("intra_dom","inter_het","het_eu","eu_eu"))]
p <- ggplot(pw, aes(band, pct, fill = cls)) +
  geom_col(position = position_dodge(.8), width = .72) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(.8),
                width = .22, linewidth = .4) +
  geom_hline(yintercept = 0, linewidth = .3) +
  scale_fill_disc(name = NULL) +
  labs(title = "Heterochromatin contacts: within a domain vs between domains",
       subtitle = paste0("Change in mean NPMI, HP1gKO vs WT, within each distance band.\n",
                         "Bars are 95% CIs from a bootstrap over the 19 autosomes."),
       x = "contact distance", y = "change in mean NPMI (%)") +
  theme_m + theme(legend.position = "top", text = element_text(size = 14))
save_fig(p, "heterochromatin_contact_classes", subdir = "heterochromatin", width = 11, height = 7)
cat("\nGAM_10 complete.\n")
