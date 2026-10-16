# ================================================================
# GAM_15 - Per-chromosome characterisation: chr3, 7, 11, 16, 18
# ----------------------------------------------------------------
#     cd GAM && Rscript GAM_15_Chromosome_Reports.R        (after GAM_1, GAM_2, GAM_6)
#
# WHY THESE FIVE. The data producers flagged chr11 and chr16 as showing the
# largest deviations; GAM_6 found 10-25 Mb one-directional blocks covering 41.6%
# of the genome, with chr11 split almost in half; chr3 and chr7 carry large
# blocks too; chr18 holds the protocadherin cluster. This stage answers "what is
# actually happening on this chromosome" with numbers rather than impressions,
# and writes one findings file per chromosome so the answer survives the session.
#
# A RETRACTION THIS STAGE EXISTS PARTLY TO PREVENT. An earlier agent reported
# chr11 as having "the strongest mark-contact coupling of any chromosome,
# rho +0.51 to +0.66, rank 1 of 19". That agent produced no output file and the
# numbers had no source. Directly computed, chr11 IS rank 1 of 19 for
# rho(H4K20me3 log2FC, net contact change) - but at rho = 0.128, with chr17 at
# -0.149, larger in the opposite direction, and a genome-wide value of 0.012.
# The spread across chromosomes is small and near-symmetric about zero. So the
# per-chromosome ranking is computed here for ALL 19, printed with its spread,
# and any claim about one chromosome must be read against that spread.
# ================================================================
source("config.R")
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(GenomicRanges); library(IRanges)
  library(patchwork)
})

set.seed(1)
FOCUS   <- c("chr3", "chr7", "chr11", "chr16", "chr18")
BANDS   <- list(c(2e4, 1e6), c(1e6, 5e6), c(5e6, 1e7), c(1e7, Inf))
BANDNM  <- c("0.02-1Mb", "1-5Mb", "5-10Mb", ">10Mb")
WIN     <- 2e6                      # smoothing window for the along-chromosome view
fig_sub <- "chromosomes"
msg <- function(...) cat(sprintf(...), "\n", sep = "")

X <- readRDS(file.path(rds_dir, "gam_bins_annotated.rds"))
A <- as.data.table(as.data.frame(X$annot)); setorder(A, chr, start)
for (mk in colnames(X$l2fc)) A[[mk]] <- X$l2fc[match(A$bin_id, rownames(X$l2fc)), mk]
A <- A[chr %in% gam_autosomes]
A[, net := as.numeric(n_contacts_gained) - as.numeric(n_contacts_lost)]

FE <- tryCatch(as.data.table(readRDS(file.path(rds_dir, "npmi_bin_features.rds"))$features),
               error = function(e) NULL)
if (!is.null(FE)) A <- merge(A, FE[, .(bin_id, insulation_WT, insulation_HP1gKO,
                                       insulation_delta, npmi_short_WT, npmi_short_HP1gKO)],
                             by = "bin_id", all.x = TRUE)
setorder(A, chr, start)

blocks <- tryCatch(fread(file.path(tables_dir, "gam6_large_blocks.tsv")),
                   error = function(e) NULL)
lads <- tryCatch(load_lads(), error = function(e) NULL)

# --- heterochromatin membership (same definition as GAM_13 bin-level) ---------
A[, het_frac := {
  dm <- fread(het_domain_bed, header = FALSE, select = 1:3,
              col.names = c("chr", "start", "end"))
  dm[, chr := ifelse(grepl("^chr", chr), chr, paste0("chr", chr))]
  dm <- dm[chr %in% gam_autosomes]
  red <- reduce(GRanges(dm$chr, IRanges(dm$start + 1, dm$end)), ignore.strand = TRUE)
  b <- GRanges(chr, IRanges(start + 1, end))
  h <- findOverlaps(b, red, ignore.strand = TRUE)
  w <- width(pintersect(b[queryHits(h)], red[subjectHits(h)], ignore.strand = TRUE))
  v <- numeric(length(b)); agg <- tapply(w, queryHits(h), sum)
  v[as.integer(names(agg))] <- as.numeric(agg)
  pmin(v / width(b), 1)
}]
A[, is_het := het_frac >= 0.5]
stopifnot(sum(A$is_het) > 1000)

if (!is.null(lads)) {
  for (cl in c("cLAD", "ciLAD")) {
    L <- reduce(lads[lads$class == cl], ignore.strand = TRUE)
    b <- GRanges(A$chr, IRanges(A$start + 1, A$end))
    A[[paste0("frac_", cl)]] <- {
      h <- findOverlaps(b, L, ignore.strand = TRUE)
      w <- width(pintersect(b[queryHits(h)], L[subjectHits(h)], ignore.strand = TRUE))
      v <- numeric(length(b)); agg <- tapply(w, queryHits(h), sum)
      v[as.integer(names(agg))] <- as.numeric(agg); pmin(v / width(b), 1)
    }
  }
}

# ==============================================================================
# 1. The genome-wide ranking these five must be read against
# ==============================================================================
rk <- A[, .(rho_H4K20me3 = as.numeric(suppressWarnings(
              cor(H4K20me3, net, method = "spearman", use = "complete.obs"))),
            rho_H3K9me3 = as.numeric(suppressWarnings(
              cor(H3K9me3, net, method = "spearman", use = "complete.obs"))),
            med_H4K20me3 = as.numeric(median(H4K20me3, na.rm = TRUE)),
            frac_het = mean(is_het), n_bins = .N), by = chr]
setorder(rk, -rho_H4K20me3); rk[, rank := seq_len(.N)]
fwrite(rk, file.path(tables_dir, "gam15_chromosome_ranking.tsv"), sep = "\t")
msg("=== per-chromosome rho(H4K20me3 log2FC, net contact change), all 19 ===")
print(rk[, .(chr, rho_H4K20me3 = round(rho_H4K20me3, 3),
             rho_H3K9me3 = round(rho_H3K9me3, 3), rank, focus = chr %in% FOCUS)])
msg("genome-wide rho = %.4f; spread across chromosomes %.3f to %.3f - small and near-symmetric,",
    as.numeric(cor(A$H4K20me3, A$net, method = "spearman", use = "complete.obs")),
    min(rk$rho_H4K20me3), max(rk$rho_H4K20me3))
msg("so a chromosome's RANK is weak evidence on its own.")

p <- ggplot(rk, aes(reorder(chr, rho_H4K20me3), rho_H4K20me3, fill = chr %in% FOCUS)) +
  geom_hline(yintercept = 0, linewidth = .3, colour = "grey50") +
  geom_col(width = .7) + coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "grey75"), guide = "none") +
  labs(title = "Mark-contact coupling by chromosome",
       subtitle = paste0("Spearman rho between per-bin H4K20me3 log2FC and net contact change. ",
                         "Red = the five examined here.\nThe whole spread is -0.15 to +0.13 around a ",
                         "genome-wide 0.012, so the ordering is weak evidence."),
       x = NULL, y = "Spearman rho") +
  theme_m
save_fig(p, "chromosome_ranking", subdir = fig_sub, width = 7, height = 6)

# ==============================================================================
# 2. Per-chromosome: contact budget by distance and chromatin class
# ==============================================================================
band_of <- function(d) {
  out <- rep(NA_integer_, length(d))
  for (i in seq_along(BANDS))
    out[is.na(out) & d >= BANDS[[i]][1] & d < BANDS[[i]][2]] <- i
  out
}

budget <- list()
for (ch in FOCUS) {
  sub <- A[chr == ch]
  M <- tryCatch(list(WT = load_npmi("WT", ch), KO = load_npmi("HP1gKO", ch)),
                error = function(e) NULL)
  if (is.null(M)) { msg("%s: no NPMI, skipped", ch); next }
  ids <- intersect(sub$bin_id, rownames(M$WT)); ids <- intersect(ids, rownames(M$KO))
  sub <- sub[match(ids, bin_id)]
  W <- M$WT[ids, ids]; K <- M$KO[ids, ids]
  iu <- which(upper.tri(W), arr.ind = TRUE)
  d  <- abs(sub$start[iu[, 2]] - sub$start[iu[, 1]])
  bd <- band_of(d)
  hi <- sub$is_het[iu[, 1]]; hj <- sub$is_het[iu[, 2]]
  cls <- fifelse(hi & hj, "het-het", fifelse(!hi & !hj, "eu-eu", "het-eu"))
  wv <- W[iu]; kv <- K[iu]
  ok <- is.finite(wv) & is.finite(kv) & !is.na(bd)
  dt <- data.table(band = BANDNM[bd[ok]], cls = cls[ok], wv = wv[ok], kv = kv[ok])
  budget[[ch]] <- dt[, .(WT = mean(wv), KO = mean(kv), n = .N,
                         pct = (mean(kv) - mean(wv)) / mean(wv) * 100),
                     by = .(band, cls)][, chr := ch][]
  rm(M, W, K, iu, wv, kv, dt); gc(FALSE)
  msg("%s: contact budget done (%d bins)", ch, length(ids))
}
BUD <- rbindlist(budget)
fwrite(BUD, file.path(tables_dir, "gam15_contact_budget.tsv"), sep = "\t")
msg("\n=== contact budget by chromosome, class and distance ===")
print(dcast(BUD, chr + cls ~ band, value.var = "pct"))

BUD[, band := factor(band, levels = BANDNM)]
p <- ggplot(BUD, aes(band, pct, colour = cls, group = cls)) +
  geom_hline(yintercept = 0, linewidth = .3, colour = "grey50") +
  geom_line() + geom_point(size = 1.6) +
  facet_wrap(~ chr, nrow = 1) +
  scale_colour_manual(values = c(`het-het` = "#b2182b", `eu-eu` = "#2166ac",
                                 `het-eu` = "grey60"), name = NULL) +
  labs(title = "Contact budget by chromatin class and separation",
       subtitle = paste0("Percent change in mean NPMI, HP1gKO vs WT. Genome-wide the pattern is het-het up ",
                         "and eu-eu down;\na chromosome that inverts that is the interesting case."),
       x = "genomic separation", y = "% change in mean NPMI") +
  theme_m
save_fig(p, "contact_budget_by_chromosome", subdir = fig_sub, width = 13, height = 4.4)

# ==============================================================================
# 3. Along-chromosome overview + a findings file per chromosome
# ==============================================================================
roll <- function(y, k) {
  g <- is.finite(y)
  num <- as.numeric(stats::filter(ifelse(g, y, 0), rep(1, k), sides = 2))
  den <- as.numeric(stats::filter(g * 1, rep(1, k), sides = 2))
  ifelse(is.finite(den) & den > 0, num / den, NA_real_)
}
kw <- max(3L, as.integer(WIN / gam_binsize))

dir.create(file.path(results_dir, "chromosomes"), showWarnings = FALSE, recursive = TRUE)

for (ch in FOCUS) {
  sub <- copy(A[chr == ch]); setorder(sub, start)
  sub[, `:=`(net_s = roll(net, kw), h4_s = roll(H4K20me3, kw),
             ins_s = if ("insulation_delta" %in% names(sub)) roll(insulation_delta, kw) else NA_real_)]
  sub[, mid := (start + end) / 2e6]

  tracks <- list(
    ggplot(sub, aes(mid, net_s)) +
      geom_hline(yintercept = 0, linewidth = .3, colour = "grey50") +
      geom_area(aes(fill = net_s > 0), na.rm = TRUE) +
      scale_fill_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "#2166ac"), guide = "none") +
      labs(y = "net diff.\ncontacts", x = NULL) + theme_m +
      theme(axis.text.x = element_blank(), axis.title.y = element_text(size = 7)),
    ggplot(sub, aes(mid, h4_s)) +
      geom_hline(yintercept = median(A$H4K20me3, na.rm = TRUE), linewidth = .3, colour = "grey40") +
      geom_line(colour = "#53426b", na.rm = TRUE) +
      labs(y = "H4K20me3\nlog2FC", x = NULL) + theme_m +
      theme(axis.text.x = element_blank(), axis.title.y = element_text(size = 7)),
    ggplot(sub, aes(mid, het_frac)) +
      geom_area(fill = "grey55", na.rm = TRUE) +
      labs(y = "het\nfraction", x = sprintf("%s (Mb)", ch)) + theme_m +
      theme(axis.title.y = element_text(size = 7)))
  if ("frac_cLAD" %in% names(sub)) {
    tracks <- append(tracks, list(
      ggplot(sub, aes(mid, frac_cLAD)) + geom_area(fill = "#53426b", na.rm = TRUE) +
        labs(y = "cLAD\nfraction", x = NULL) + theme_m +
        theme(axis.text.x = element_blank(), axis.title.y = element_text(size = 7))), after = 2)
  }
  pp <- wrap_plots(tracks, ncol = 1, heights = rep(1, length(tracks))) +
    plot_annotation(title = sprintf("%s overview", ch),
      subtitle = sprintf("All tracks smoothed over %.0f Mb. Net diff. contacts = gained - lost per 20 kb bin; red = more contact in the KO.", WIN / 1e6),
      theme = theme(plot.title = element_text(face = "bold", size = 13)))
  save_fig(pp, sprintf("%s_overview", ch), subdir = fig_sub, width = 11, height = 2 + 1.5 * length(tracks))

  # --- the findings file ---
  bl <- if (!is.null(blocks)) blocks[chr == ch] else NULL
  bu <- BUD[chr == ch]
  f <- file.path(results_dir, "chromosomes", sprintf("%s.md", ch))
  con <- file(f, "w")
  writeLines(c(
    sprintf("# %s", ch), "",
    sprintf("Generated %s by `GAM_15_Chromosome_Reports.R`. Numbers here are regenerable; do not hand-edit.", Sys.Date()), "",
    "## Mark-contact coupling, in context", "",
    sprintf("- rho(H4K20me3 log2FC, net contact change) = **%.3f**, rank %d of %d.",
            rk[chr == ch, rho_H4K20me3], rk[chr == ch, rank], nrow(rk)),
    sprintf("- The spread across all 19 autosomes is %.3f to %.3f around a genome-wide %.4f, so this rank is weak evidence on its own.",
            min(rk$rho_H4K20me3), max(rk$rho_H4K20me3),
            as.numeric(cor(A$H4K20me3, A$net, method = "spearman", use = "complete.obs"))),
    sprintf("- Median H4K20me3 log2FC here %.3f, against %.3f genome-wide.",
            rk[chr == ch, med_H4K20me3], median(A$H4K20me3, na.rm = TRUE)),
    sprintf("- %.1f%% of bins are heterochromatic (>=50%% coverage by the merged H3K9me3/H4K20me3 set), against %.1f%% genome-wide.",
            100 * rk[chr == ch, frac_het], 100 * mean(A$is_het)), "",
    "## Contact budget", "",
    "Percent change in mean NPMI, HP1gKO vs WT, by chromatin class and separation.", "",
    "| class | 0.02-1Mb | 1-5Mb | 5-10Mb | >10Mb |", "|---|---|---|---|---|"),
    con)
  for (cl in c("het-het", "het-eu", "eu-eu")) {
    r <- sapply(BANDNM, function(b) { v <- bu[cls == cl & band == b, pct]
                                      if (length(v)) sprintf("%+.1f%%", v) else "-" })
    writeLines(sprintf("| %s | %s |", cl, paste(r, collapse = " | ")), con)
  }
  writeLines(c("",
    sprintf("Genome-wide the pattern is het-het gaining and eu-eu losing; %s %s that.", ch,
            if (nrow(bu) && isTRUE(bu[cls == "het-het", mean(pct, na.rm = TRUE)] <
                                   bu[cls == "eu-eu", mean(pct, na.rm = TRUE)]))
              "**inverts**" else "follows"), ""), con)
  if (!is.null(bl) && nrow(bl)) {
    writeLines(c("## Large one-directional blocks (GAM_6)", "",
      sprintf("%d block(s) on this chromosome, %.1f Mb in total.", nrow(bl),
              sum(bl$end - bl$start) / 1e6), ""), con)
    writeLines("| start (Mb) | end (Mb) | width (Mb) | direction |", con)
    writeLines("|---|---|---|---|", con)
    for (i in seq_len(nrow(bl)))
      writeLines(sprintf("| %.2f | %.2f | %.2f | %s |", bl$start[i] / 1e6, bl$end[i] / 1e6,
                         (bl$end[i] - bl$start[i]) / 1e6,
                         if ("direction" %in% names(bl)) as.character(bl$direction[i])
                         else if ("class" %in% names(bl)) as.character(bl$class[i]) else "?"), con)
    writeLines("", con)
  }
  writeLines(c("## Caveats carried from the pipeline", "",
    "- GAM is 2 biological replicates per genotype, already merged, so there is no within-genotype variance; every claim here is an effect size against a matched or bootstrap null, never a per-contact p-value.",
    "- Per-bin ChIP log2FC at 20 kb is noise-dominated (per-bin sd 0.14-0.16 for the broad marks). Chromosome-level medians are sound; do not read local structure off a single bin.",
    "- The A/B compartment axis is weakly determined on most chromosomes and is not used here.",
    ""), con)
  close(con)
  msg("wrote %s", f)
}

msg("\nGAM_15 complete.")
