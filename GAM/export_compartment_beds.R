# ================================================================
# export_compartment_beds.R - A/B compartment BEDs for IGV
#   cd GAM && Rscript export_compartment_beds.R
# Utility, not a pipeline stage. Needs GAM_1 and GAM_2b.
# Writes results/bed/compartments_{WT,HP1gKO,diff}.bed, to be loaded next to
# results/igv/gam_chip_minimal.xml. Read the "A/B compartments" section of
# CLAUDE.md before quoting anything off these tracks.
# ================================================================
# Write A/B compartment BEDs for IGV from the SHARED axis (GAM_2b).
#
# WHY THE SHARED AXIS AND NOT THE PER-GENOTYPE PCs
#   GAM_2 selects an eigenvector independently per genotype. Even after the v8
#   balancing fix, measured on results/tables/npmi_compartment_orientation.tsv:
#   only 27 of 38 chromosome x genotype combinations orient at all, both
#   genotypes orient on just 10 of 19 autosomes, and of those 10 the two
#   genotypes pick the SAME eigenvector on only 6. A per-genotype BED would
#   therefore be blank on half the genome and would carry a "difference" that is
#   partly a change of eigenvector rather than a change of compartment.
#   (Inverted axes are no longer possible: GAM_2 blanks an anti-correlating
#   axis to NA and asserts the invariant before writing, so every oriented row
#   has rho_check > 0. That was not true when this file was first written.)
#   The shared axis is one vector scored in both genotypes and is determined on
#   18 of 19 autosomes at a median |rho| of 0.42 against active chromatin,
#   against 0.36 for the per-genotype axes.
#
# Coordinates: $start/$end in shared_axis_scores.rds are already 0-based
# half-open (gam_bin_id() convention). Written through unchanged - NO conversion.
source("config.R")
suppressPackageStartupMessages(library(data.table))
msg <- function(...) message(sprintf(...))

SA <- readRDS(file.path(rds_dir, "shared_axis_scores.rds"))
stopifnot(!isTRUE(SA$partial))
F <- copy(SA$fine)
stopifnot(all(c("chr","start","end","score_WT","score_HP1gKO") %in% names(F)))

# --- coverage assert: a silent drop here would look like a clean result -------
# Measured PER CHROMOSOME, against the chromosomes that actually HAVE a shared
# axis. A flat genome-wide threshold conflates two different things: a
# chromosome the selector legitimately refused (no eigenvector beats the
# pc_margin - GAM_2b writes determined = FALSE and scores nothing), and bins
# vanishing from a chromosome that does have an axis, which is the silent bug
# this guard exists to catch. Only the second is an error.
gbins   <- gam_bins()
grid    <- gam_bin_id(gbins)
grid_dt <- data.table(chr = as.character(GenomeInfoDb::seqnames(gbins)))[, .(N = .N), by = chr]
stopifnot(all(F$bin_id %in% grid))

have   <- F[, unique(chr)]
absent <- setdiff(grid_dt$chr, have)
if (length(absent))
  msg("NO SHARED AXIS on %d chromosome(s): %s - excluded from these tracks",
      length(absent), paste(absent, collapse = ", "))

cov_chr <- merge(F[, .(n_scored = .N), by = chr], grid_dt, by = "chr")
cov_chr[, cov := n_scored / N]
msg("shared-axis bins on the 20 kb autosomal grid: %d of %d (%.1f%%); per-chromosome coverage min %.3f on %s",
    nrow(F), sum(grid_dt$N), 100 * nrow(F) / sum(grid_dt$N),
    min(cov_chr$cov), cov_chr[which.min(cov), chr])
# Floor is 0.90, not 0.95. Measured 2026-10-18 the real per-chromosome range is
# 0.939 (chr7) to ~0.99; chr7/chr14/chr19 sit just under 0.95 from ordinary
# unmappable coarse bins, and chr7 is the repeat-richest autosome here
# (lambda1 235 vs ~130 typical). The old flat 0.95 was a genome-wide AVERAGE
# that passed only because the better chromosomes carried it. The message above
# prints the running minimum, so a drift shows up even when it does not trip.
bad <- cov_chr[cov <= 0.90]
if (nrow(bad))
  stop("bins missing from chromosomes that DO have a shared axis: ",
       paste(sprintf("%s %.1f%%", bad$chr, 100 * bad$cov), collapse = ", "))
F <- F[is.finite(score_WT) & is.finite(score_HP1gKO)]
setorder(F, chr, start)
# Gaps are expected (unmappable coarse bins are dropped). merge_runs() breaks a
# run at any gap, so no interval ever spans one - check that assumption holds.
ngap <- F[, sum(diff(start) != gam_binsize), by = chr][, sum(V1)]
msg("coverage gaps in the scored grid: %d (intervals will never span them)", ngap)

# --- calls --------------------------------------------------------------------
# A = positive side of the shared axis = the WT H3K4me3+H3K36me3-high side.
# NOTE: GAM_2b median-centres each genotype's score per chromosome, so frac(A) is
# 0.5 in BOTH genotypes BY CONSTRUCTION. These tracks can show WHICH bins swap;
# they cannot show a net A<->B shift, and must not be read as showing one.
F[, call_WT := fifelse(score_WT      > 0, "A", "B")]
F[, call_KO := fifelse(score_HP1gKO  > 0, "A", "B")]
F[, conf := pmin(abs(score_WT), abs(score_HP1gKO))]
msg("frac A  WT %.4f  cKO %.4f  (0.5 by construction)",
    mean(F$call_WT=="A"), mean(F$call_KO=="A"))

# --- orientation QC against a ruler with no ChIP in it (TSS density) ----------
g1 <- readRDS(file.path(rds_dir, "gam_bins_annotated.rds"))
A <- g1$annot[, .(bin_id, n_tss)]
Q <- merge(F[, .(bin_id, chr, call_WT, call_KO)], A, by = "bin_id")
ovl <- nrow(Q)/nrow(F)
msg("QC join (shared-axis bins -> GAM_1 annotation): %.3f", ovl)
stopifnot(ovl > 0.95)                        # cross-table join, asserted not warned
msg("TSS per bin: WT  A %.3f vs B %.3f   (ratio %.2f)",
    Q[call_WT=="A", mean(n_tss)], Q[call_WT=="B", mean(n_tss)],
    Q[call_WT=="A", mean(n_tss)]/Q[call_WT=="B", mean(n_tss)])
msg("TSS per bin: cKO A %.3f vs B %.3f   (ratio %.2f)",
    Q[call_KO=="A", mean(n_tss)], Q[call_KO=="B", mean(n_tss)],
    Q[call_KO=="A", mean(n_tss)]/Q[call_KO=="B", mean(n_tss)])
bad <- Q[, .(r = mean(n_tss[call_WT=="A"])/mean(n_tss[call_WT=="B"])), by=chr][r < 1]
msg("chromosomes where WT 'A' is NOT the TSS-richer side: %s",
    if (nrow(bad)) paste(bad$chr, collapse=", ") else "none")
bad2 <- Q[, .(r = mean(n_tss[call_KO=="A"])/mean(n_tss[call_KO=="B"])), by=chr][r < 1]
msg("chromosomes where cKO 'A' is NOT the TSS-richer side: %s",
    if (nrow(bad2)) paste(bad2$chr, collapse=", ") else "none")

# --- merge adjacent same-call bins into intervals ------------------------------
merge_runs <- function(d, lab, val) {
  d <- copy(d); setorder(d, chr, start)
  d[, brk := cumsum(c(TRUE, chr[-1] != chr[-.N] |
                        start[-1] != end[-.N] |
                        get(lab)[-1] != get(lab)[-.N]))]
  d[, .(start = min(start), end = max(end), lab = get(lab)[1],
        n_bins = .N, value = mean(get(val))), by = .(brk, chr)][, brk := NULL][]
}

RED  <- "178,24,43"    # repo convention: red = the more / active-in-cKO direction
BLUE <- "33,102,172"
GREY <- "140,140,140"

write_bed <- function(dt, file, name, desc, colmap) {
  p <- file.path(bed_dir, file)
  writeLines(sprintf('track name="%s" description="%s" itemRgb=On visibility=dense',
                     name, desc), p)
  out <- dt[, .(chr, start, end, nm, score, ".", start, end, colmap[lab])]
  fwrite(out, p, sep="\t", col.names=FALSE, append=TRUE, quote=FALSE, scipen=999)
  msg("  %-28s %6d intervals  (median %.0f kb)", file, nrow(dt),
      median(dt$end - dt$start)/1000)
  p
}
scale1000 <- function(v, cap) as.integer(pmin(1000, round(1000 * abs(v) / cap)))

# ---- WT ----------------------------------------------------------------------
W <- merge_runs(F, "call_WT", "score_WT")
cap <- quantile(abs(F$score_WT), 0.95)
W[, nm := sprintf("%s|%.2f", lab, value)][, score := scale1000(value, cap)]
f1 <- write_bed(W, "compartments_WT.bed", "GAM A/B compartments - WT",
  paste0("A/B from the SHARED axis (GAM_2b averaged-correlation axis, phased once on WT ",
         "H3K4me3+H3K36me3); A = positive = active-chromatin side. RED = A (active, ",
         "the direction that is 'more active' and so matches the repo's red=more-in-cKO ",
         "sense), BLUE = B. Score = |axis score| scaled to q95; faint intervals sit near ",
         "the zero crossing and are not reliable calls. WEAK AXIS: |rho| vs active ",
         "chromatin reaches 0.4 on only 10 of 19 autosomes (median 0.42); do not quote an ",
         "A/B assignment without that caveat."),
  c(A = RED, B = BLUE))

# ---- HP1gKO ------------------------------------------------------------------
K <- merge_runs(F, "call_KO", "score_HP1gKO")
cap <- quantile(abs(F$score_HP1gKO), 0.95)
K[, nm := sprintf("%s|%.2f", lab, value)][, score := scale1000(value, cap)]
f2 <- write_bed(K, "compartments_HP1gKO.bed", "GAM A/B compartments - HP1gKO (cKO)",
  paste0("Same SHARED axis as the WT track - one ruler, both genotypes, so the two are ",
         "comparable bin for bin; the per-genotype PC1s are NOT (the genotypes select the ",
         "same eigenvector on only 8 of 19 autosomes). RED = A, BLUE = B. Each genotype's ",
         "score is median-centred per chromosome, so exactly 50% of bins are A in BOTH ",
         "tracks BY CONSTRUCTION - a net A<->B shift cannot be read off these files."),
  c(A = RED, B = BLUE))

# ---- diff --------------------------------------------------------------------
D <- F[call_WT != call_KO]
D[, sw := paste0(call_WT, "_to_", call_KO)]
Dm <- merge_runs(D, "sw", "conf")
Dm[, nm := sprintf("%s|conf=%.2f", lab, value)]
Dm[, score := scale1000(value, quantile(F$conf, 0.95))]
msg("switching bins: %d of %d (%.1f%%)  A->B %d  B->A %d",
    nrow(D), nrow(F), 100*nrow(D)/nrow(F),
    sum(D$sw=="A_to_B"), sum(D$sw=="B_to_A"))
f3 <- write_bed(Dm, "compartments_diff.bed", "GAM A/B compartment CHANGE (cKO vs WT)",
  paste0("READ WITH EXTREME CAUTION - A COMPARTMENT 'CHANGE' IS NOT WELL DEFINED IN THIS ",
         "DATA. RED = B->A in the cKO (toward active), BLUE = A->B. 28.5% of bins change ",
         "call on the shared axis, but the two genotypes' scores on that one fixed axis ",
         "agree with each other at only Spearman rho 0.56 (median over autosomes) while ",
         "the axis itself reproduces at rho 0.98 between random halves of the bins. Among ",
         "bins where BOTH genotypes sit well away from the zero crossing (conf >= 0.77, ",
         "the q75 of min|score|) only 6.6% change call, i.e. 1,975 bins = 1.7% of the ",
         "genome - the other 94% of this file is sign noise at the A/B boundary, not a ",
         "compartment transition. There is also NO net direction: scores are median-centred ",
         "per chromosome, so A->B and B->A are forced equal (16,960 bins each). Filter on ",
         "the conf= value in the name. Use as a lead list only, never as a result."),
  c(A_to_B = BLUE, B_to_A = RED))

msg("\nconfidence filter on the diff track:")
for (q in c(0, .5, .75, .9)) {
  thr <- quantile(F$conf, q)
  msg("  conf >= q%02d (%.2f): %d switching bins (%.1f%% of all bins)",
      100*q, thr, D[conf>=thr, .N], 100*D[conf>=thr,.N]/nrow(F))
}
msg("\nwrote:\n  %s\n  %s\n  %s", f1, f2, f3)
