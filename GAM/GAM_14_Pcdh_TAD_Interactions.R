# ================================================================
# GAM_14 - The protocadherin cluster as a TAD: who does it talk to?
# ----------------------------------------------------------------
#     cd GAM && Rscript GAM_14_Pcdh_TAD_Interactions.R     (after GAM_1, GAM_2)
#
# THE QUESTION. GAM_5 and GAM_9 looked inside the Pcdh cluster. But the cluster
# is itself approximately one TAD, so the sharper question is how that TAD's
# contacts with OTHER, more distal TADs change in the KO - decompaction inside a
# domain and decoupling of that domain from the rest of the chromosome are
# different events, and only the second needs the inter-TAD view.
#
# THE CONFOUND, and the test. Contact frequency falls steeply with separation,
# so "Pcdh contacts distant TADs less in the KO" is uninformative unless it is
# compared against other TADs at the SAME separations. The test is therefore a
# permutation over QUERY TADs, stratified by separation: every other TAD on the
# chromosome plays the part of Pcdh in turn, and we ask where Pcdh's
# separation-matched mean change falls in that distribution. The null is "this
# TAD behaves like a typical chr18 TAD at matched distance", which is the null
# worth rejecting.
#
# TWO PARTITIONS, REPORTED SIDE BY SIDE. Where the boundaries are decides what
# counts as intra vs inter, so a result that holds under only one partition is a
# result about the partition, not the biology:
#   external - published mm39 boundaries, fixed across genotypes, non-circular
#              but from a different cell type and blind to anything the KO does
#   own      - insulation minima from our own GAM, per genotype, self-consistent
#              with the contact data and free to move, but it is our own call
# Both come from config.R (load_tads_external / call_tads_insulation).
# ================================================================
source("config.R")
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(GenomicRanges); library(IRanges)
  library(TxDb.Mmusculus.UCSC.mm39.knownGene); library(org.Mm.eg.db)
})

set.seed(1)
N_PERM   <- 2000L
SEP_BRK  <- c(0, 2e6, 5e6, 1e7, 2e7, Inf)
SEP_NM   <- c("<2Mb", "2-5Mb", "5-10Mb", "10-20Mb", ">20Mb")
CHR      <- "chr18"
fig_sub  <- "pcdh_tad"
msg <- function(...) cat(sprintf(...), "\n", sep = "")

# --- where is the Pcdh cluster ------------------------------------------------
pcdh <- local({
  g <- genes(TxDb.Mmusculus.UCSC.mm39.knownGene, single.strand.genes.only = TRUE)
  sym <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, keys = names(g),
           column = "SYMBOL", keytype = "ENTREZID"))
  g <- g[!is.na(sym) & grepl("^Pcdh(a|b|g)", sym)]
  g <- g[as.character(seqnames(g)) == CHR]
  stopifnot(length(g) > 10)
  GRanges(CHR, IRanges(min(start(g)), max(end(g))))
})
msg("Pcdh cluster: %s:%.2f-%.2f Mb (%.0f kb)", CHR,
    start(pcdh) / 1e6, end(pcdh) / 1e6, width(pcdh) / 1e3)

# --- NPMI for the chromosome, once --------------------------------------------
MW <- load_npmi("WT", CHR); MK <- load_npmi("HP1gKO", CHR)
ids <- intersect(rownames(MW), rownames(MK))
MW <- MW[ids, ids]; MK <- MK[ids, ids]
bpos <- as.numeric(sub("-.*", "", sub(".*:", "", ids)))
msg("NPMI loaded: %d bins", length(ids))

# mean NPMI over the bins of one TAD pair, and over a TAD's own interior
pair_mean <- function(M, i, j) {
  if (!length(i) || !length(j)) return(NA_real_)
  v <- M[i, j, drop = FALSE]
  if (identical(i, j)) v <- v[upper.tri(v)]
  v <- v[is.finite(v)]
  if (length(v) < 5) NA_real_ else mean(v)
}

analyse <- function(tads, label) {
  t <- as.data.table(tads)[chr == CHR]
  setorder(t, start)
  if (nrow(t) < 10) { msg("%s: only %d TADs on %s, skipped", label, nrow(t), CHR); return(NULL) }
  # bins of each TAD
  t[, idx := .(list(integer(0)))]
  blist <- lapply(seq_len(nrow(t)), function(k)
    which(bpos >= t$start[k] & bpos < t$end[k]))
  keep <- which(lengths(blist) >= 3)      # a TAD needs interior to have a mean
  t <- t[keep]; blist <- blist[keep]
  t[, mid := (start + end) / 2]
  # which TAD holds the Pcdh cluster: the one with the largest overlap
  ovl <- pmax(0, pmin(t$end, end(pcdh)) - pmax(t$start, start(pcdh)))
  q0  <- which.max(ovl)
  msg("\n=== %s: %d usable TADs on %s; Pcdh sits in %s:%.2f-%.2f Mb (%.0f%% of the cluster) ===",
      label, nrow(t), CHR, CHR, t$start[q0] / 1e6, t$end[q0] / 1e6,
      100 * ovl[q0] / width(pcdh))

  # --- all query x partner pair means ---
  n <- nrow(t)
  grid <- CJ(q = seq_len(n), p = seq_len(n))[q != p]
  grid[, sep := abs(t$mid[q] - t$mid[p])]
  grid[, `:=`(WT = mapply(function(a, b) pair_mean(MW, blist[[a]], blist[[b]]), q, p),
              KO = mapply(function(a, b) pair_mean(MK, blist[[a]], blist[[b]]), q, p))]
  grid <- grid[is.finite(WT) & is.finite(KO) & WT > 0]
  grid[, pct := (KO - WT) / WT * 100]
  grid[, sep_band := cut(sep, SEP_BRK, labels = SEP_NM, right = FALSE)]

  # --- intra-TAD, every TAD ---
  intra <- data.table(
    q = seq_len(n), width = t$end - t$start,
    WT = vapply(seq_len(n), function(k) pair_mean(MW, blist[[k]], blist[[k]]), 0),
    KO = vapply(seq_len(n), function(k) pair_mean(MK, blist[[k]], blist[[k]]), 0))
  intra <- intra[is.finite(WT) & is.finite(KO) & WT > 0][, pct := (KO - WT) / WT * 100]

  # --- the test: separation-stratified mean per query TAD ---
  # Average within separation band first, then across bands, so a TAD is not
  # scored by how many near partners it happens to have.
  per_q <- grid[, .(m = mean(pct)), by = .(q, sep_band)][
             , .(stratified = mean(m), n_bands = .N), by = q][n_bands >= 3]
  obs <- per_q[q == q0, stratified]
  nullv <- per_q[q != q0, stratified]
  if (!length(obs) || length(nullv) < 20) { msg("  too few comparable TADs"); return(NULL) }
  p_emp <- (1 + sum(abs(nullv - mean(nullv)) >= abs(obs - mean(nullv)))) / (1 + length(nullv))
  msg("  Pcdh TAD -> distal, separation-matched: %+.2f%%   (other chr18 TADs: median %+.2f%%, IQR %+.2f to %+.2f)",
      obs, median(nullv), quantile(nullv, .25), quantile(nullv, .75))
  msg("  two-sided empirical p vs %d other TADs = %.3f", length(nullv), p_emp)

  io <- intra[q == q0, pct]; inull <- intra[q != q0, pct]
  p_intra <- if (length(io) && length(inull) >= 20)
    (1 + sum(abs(inull - mean(inull)) >= abs(io - mean(inull)))) / (1 + length(inull)) else NA_real_
  msg("  Pcdh TAD intra: %+.2f%%   (other TADs median %+.2f%%), empirical p = %.3f",
      ifelse(length(io), io, NA), median(inull), p_intra)

  list(
    label = label,
    summary = data.table(
      partition = label, n_tads = nrow(t),
      pcdh_tad = sprintf("%s:%d-%d", CHR, t$start[q0], t$end[q0]),
      pcdh_tad_Mb = (t$end[q0] - t$start[q0]) / 1e6,
      frac_cluster_in_tad = ovl[q0] / width(pcdh),
      intra_pct = ifelse(length(io), io, NA_real_),
      intra_other_median = median(inull), intra_p = p_intra,
      distal_pct = obs, distal_other_median = median(nullv), distal_p = p_emp,
      n_null_tads = length(nullv)),
    by_sep = grid[q == q0, .(pct = mean(pct), n = .N), by = sep_band][
               , partition := label][],
    by_sep_null = grid[q != q0, .(pct = mean(pct)), by = .(q, sep_band)][
               , partition := label][],
    profile = grid[q == q0, .(partition = label, sep, pct,
                              partner_mid = t$mid[p])],
    nullv = data.table(partition = label, stratified = nullv),
    obs = obs)
}

parts <- list(external = load_tads_external())
for (g in c("WT", "HP1gKO")) {
  tt <- tryCatch(call_tads_insulation(g), error = function(e) { msg("insulation TADs (%s) failed: %s", g, conditionMessage(e)); NULL })
  if (!is.null(tt)) parts[[paste0("own_", g)]] <- tt
}
for (nm in names(parts)) msg("%s: %d TADs genome-wide", nm, nrow(parts[[nm]]))

R <- Filter(Negate(is.null), lapply(names(parts), function(nm) analyse(parts[[nm]], nm)))
stopifnot(length(R) > 0)

SUM <- rbindlist(lapply(R, `[[`, "summary"))
fwrite(SUM, file.path(tables_dir, "gam14_pcdh_tad_summary.tsv"), sep = "\t")
msg("\n=== summary across partitions ==="); print(SUM)

BY <- rbindlist(lapply(R, `[[`, "by_sep"))
fwrite(BY, file.path(tables_dir, "gam14_pcdh_by_separation.tsv"), sep = "\t")
PROF <- rbindlist(lapply(R, `[[`, "profile"))
fwrite(PROF, file.path(tables_dir, "gam14_pcdh_partner_profile.tsv"), sep = "\t")

# --- figures ------------------------------------------------------------------
NULLD <- rbindlist(lapply(R, `[[`, "nullv"))
OBS <- SUM[, .(partition, distal_pct)]
p <- ggplot(NULLD, aes(stratified)) +
  geom_histogram(bins = 40, fill = "grey80", colour = NA) +
  geom_vline(data = OBS, aes(xintercept = distal_pct), colour = "#b2182b", linewidth = .8) +
  facet_wrap(~ partition, scales = "free_y", nrow = 1) +
  labs(title = "Pcdh TAD vs distal TADs, against every other chr18 TAD at matched separation",
       subtitle = paste0("Grey = the separation-stratified contact change of each other chr18 TAD against its own ",
                         "partners.\nRed = the Pcdh TAD. Averaging within separation band first stops a TAD being ",
                         "scored by how many near\npartners it happens to have."),
       x = "separation-matched change in mean NPMI, HP1gKO vs WT (%)", y = "TADs") +
  theme_m
save_fig(p, "pcdh_vs_null", subdir = fig_sub, width = 11, height = 4.6)

BYN <- rbindlist(lapply(R, `[[`, "by_sep_null"))
bq <- BYN[, .(lo = quantile(pct, .25, na.rm = TRUE), md = median(pct, na.rm = TRUE),
              hi = quantile(pct, .75, na.rm = TRUE)), by = .(partition, sep_band)]
p <- ggplot() +
  geom_hline(yintercept = 0, linewidth = .3, colour = "grey50") +
  geom_ribbon(data = bq, aes(as.integer(sep_band), ymin = lo, ymax = hi, group = partition),
              fill = "grey85") +
  geom_line(data = bq, aes(as.integer(sep_band), md, group = partition), colour = "grey45") +
  geom_line(data = BY, aes(as.integer(sep_band), pct, group = partition), colour = "#b2182b", linewidth = .8) +
  geom_point(data = BY, aes(as.integer(sep_band), pct), colour = "#b2182b", size = 1.8) +
  facet_wrap(~ partition, nrow = 1) +
  scale_x_continuous(breaks = seq_along(SEP_NM), labels = SEP_NM) +
  labs(title = "Where along the chromosome does the Pcdh TAD change its contacts?",
       subtitle = paste0("Red = the Pcdh TAD. Grey band = interquartile range over all other chr18 TADs at the same ",
                         "separation,\ngrey line their median - so the grey IS the distance control."),
       x = "separation from the Pcdh TAD", y = "% change in mean NPMI") +
  theme_m
save_fig(p, "pcdh_by_separation", subdir = fig_sub, width = 11, height = 4.6)

p <- ggplot(PROF, aes(partner_mid / 1e6, pct)) +
  geom_hline(yintercept = 0, linewidth = .3, colour = "grey50") +
  geom_vline(xintercept = (start(pcdh) + end(pcdh)) / 2e6, linetype = 2, colour = "grey40") +
  geom_point(aes(colour = pct > 0), size = 1.3) +
  scale_colour_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "#2166ac"), guide = "none") +
  facet_wrap(~ partition, nrow = 1) +
  labs(title = "Pcdh TAD contact change with each partner TAD along chr18",
       subtitle = "Dashed line = the Pcdh cluster. Each point is one partner TAD; red = more contact in the KO.",
       x = "partner TAD midpoint (Mb)", y = "% change in mean NPMI") +
  theme_m
save_fig(p, "pcdh_partner_profile", subdir = fig_sub, width = 11, height = 4.6)

msg("\nGAM_14 complete.")
