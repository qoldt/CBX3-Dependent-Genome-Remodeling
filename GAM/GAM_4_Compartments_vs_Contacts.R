# ================================================================
# GAM_4 - Does the chromatin-mark loss compartment have an architectural
#         counterpart in 3D?
# ----------------------------------------------------------------
# MINUTE_3 sorts regions into stable / H4K20me3_only / co_loss by how H4K20me3
# and H3K9me3 change in the HP1gKO. GAM_1 projected those labels onto the 20 kb
# GAM grid (annot$compartment_centred / compartment_abs). Here we ask whether the
# same bins also remodel their 3D contacts and whether they move between A/B
# compartments:
#   1. differential contact LOSS / GAIN density per group
#   2. size-matched permutation tests for the overlap-shaped questions
#   3. A/B switching (pc1_delta, compartment_switch) per group, and PC1 shift vs
#      H4K20me3 log2FC genome-wide
#   4. insulation change per group and specifically at TAD boundaries
#   5. per-group summary table + figures
#
# Reads results/rds/gam_bins_annotated.rds (GAM_1, REQUIRED) and
# results/rds/npmi_bin_features.rds (GAM_2). Without GAM_2 the contact-density
# sections still run; only the A/B / insulation sections are skipped.
# Joins on bin_id. Runs after GAM_1 (and GAM_2). Run from GAM/.
# ================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })

# --- Inputs ------------------------------------------------------
gam1_rds <- file.path(rds_dir, "gam_bins_annotated.rds")
gam2_rds <- file.path(rds_dir, "npmi_bin_features.rds")
if (!file.exists(gam1_rds))
  stop("GAM_1 output not found:\n  ", gam1_rds,
       "\nRun it first:   cd GAM && Rscript GAM_1_Bins_and_Signal.R", call. = FALSE)
g1 <- readRDS(gam1_rds)
d  <- copy(as.data.table(g1$annot))
stopifnot(all(c("bin_id", "n_contacts_lost", "n_contacts_gained", "tad_boundary",
                "compartment_abs", "compartment_centred") %in% names(d)))
# H4K20me3 per-bin log2FC from GAM_1's wide matrix, joined on bin_id (not row
# order: the matrix and annot are separate objects)
d[, h4_l2fc := g1$l2fc[bin_id, "H4K20me3"]]

have_feat <- file.exists(gam2_rds)
if (have_feat) {
  g2 <- readRDS(gam2_rds)
  want <- intersect(c("bin_id", "pc1_WT", "pc1_HP1gKO", "pc1_delta", "compartment_WT",
                      "compartment_HP1gKO", "compartment_switch", "insulation_delta"),
                    names(g2$features))
  d <- merge(d, g2$features[, ..want], by = "bin_id", all.x = TRUE, sort = FALSE)
  # features also cover chrX/Y; the left join keeps only the autosomal GAM_1 grid
} else {
  message("GAM_4: ", gam2_rds, " not found - skipping the A/B-switch and insulation",
          " sections (run GAM_2 first: cd GAM && Rscript GAM_2_NPMI_Compartments.R).")
}
grp_levels <- c("stable", "H4K20me3_only", "co_loss")
grp_cols   <- setNames(disc_pal(3), grp_levels)

# --- Scheme choice: READ THIS before changing the grouping --------
# PRIMARY = compartment_centred; compartment_abs is a SENSITIVITY only.
# Both marks lose genome-wide (H4K20me3 median log2FC ~ -0.68, H3K9me3 ~ -0.42), so
# the absolute cutoffs (-0.5 / -0.3 / 0.15) mostly ask "did it lose?", to which most
# of the genome says yes. Group sizes and even their ORDERING then depend on where
# the global shift sits (co_loss 46,921 vs H4K20me3_only 8,481 under abs; the
# reverse once each mark is centred on its own median) - see
# ChIP/results/COMPARTMENT_CUTOFFS.md. A 3D result that held only under the
# absolute labels would be a statement about the global shift, not about the
# compartments, so we report both and only trust what survives the centred scheme.
schemes <- c(centred = "compartment_centred", abs = "compartment_abs")
PRIMARY <- "centred"

# --- Matching covariate for size_perm() --------------------------
# size_perm() was written for REGIONS OF VARYING WIDTH and matches on width
# deciles. Every GAM bin is 20 kb, so a width match is degenerate (one unique
# break -> cut() errors). The same machinery works for any numeric nuisance
# covariate passed in the "width" slots, so we match on per-bin REPEAT CONTENT
# instead: repeats drive both the H4K20me3/H3K9me3 labels (IAP/L1 are the losers)
# and mappability/contact counts, so a naive group-vs-all comparison would partly
# be a repeat-density comparison. A tiny jitter breaks the heavy ties at 0 so the
# deciles are unique; seeded so the null is reproducible.
d[, rep_cov := rowSums(.SD, na.rm = TRUE), .SDcols = c("frac_LINE", "frac_SINE", "frac_LTR")]
matched_perm <- function(idx, flag, cov, nperm = 2000, seed = 42) {
  ok <- is.finite(cov) & !is.na(flag)
  idx <- idx[ok[idx]]
  if (length(idx) < 10) return(NULL)
  set.seed(seed); cj <- cov[ok] + runif(sum(ok), 0, 1e-6)
  map <- cumsum(ok)                    # row -> position within the ok-universe
  qi  <- map[idx]
  size_perm(mean(flag[ok][qi]), cj[qi], cj, as.numeric(flag[ok]), nperm = nperm, seed = seed)
}

# Continuous version: is the group's MEAN of x different from a repeat-content-matched
# draw from `universe` (default all bins)? size_perm() only returns an upper-tail p,
# so run it on x and -x (same seed -> identical draws) for a two-sided p.
# WHY NOT a Wilcoxon of group vs stable: the delivered GAM data merge 2 replicates
# per genotype, so there is no within-genotype variance, and 123k bins are not
# independent observations (neighbouring bins share contacts). A bin-level p-value
# would be wildly anti-conservative; the matched null asks the right question
# ("more than bins of the same repeat content") and we report the EFFECT SIZE.
matched_effect <- function(idx, x, cov, universe = rep(TRUE, length(x)), nperm = 2000, seed = 42) {
  ok <- is.finite(cov) & is.finite(x) & universe
  idx <- idx[ok[idx]]
  if (length(idx) < 10 || sum(ok) < 50) return(NULL)
  set.seed(seed); cj <- cov[ok] + runif(sum(ok), 0, 1e-6)
  qi <- cumsum(ok)[idx]; xx <- x[ok]
  hi <- size_perm(mean(xx[qi]),  cj[qi], cj,  xx, nperm = nperm, seed = seed)
  lo <- size_perm(mean(-xx[qi]), cj[qi], cj, -xx, nperm = nperm, seed = seed)
  data.table(n_bins = length(idx), obs_mean = hi$obs_frac, matched_mean = hi$exp_frac,
             diff = hi$obs_frac - hi$exp_frac, z = hi$z,
             perm_p = min(1, 2 * min(hi$perm_p, lo$perm_p)))
}

# ================================================================
# 1. Contact LOSS / GAIN density per group
# ================================================================
# Direction (config.R): selWT = stronger in WT = LOST in HP1gKO -> n_contacts_lost;
# selcKO = GAINED in HP1gKO -> n_contacts_gained. Columns are KO-centric.
# The hotspot counts are DENSE (almost every mappable bin has a value, 1-162), so the
# signal is count MAGNITUDE; compare distributions, not "has any".
metrics <- intersect(c("n_contacts_lost", "n_contacts_gained"), names(d))
summ <- rbindlist(lapply(names(schemes), function(sc) {
  x <- copy(d); x[, grp := factor(get(schemes[[sc]]), levels = grp_levels)]
  x <- x[!is.na(grp)]
  out <- x[, c(list(n_bins = .N,
                    median_lost = as.numeric(median(n_contacts_lost, na.rm = TRUE)),
                    mean_lost   = mean(n_contacts_lost,   na.rm = TRUE),
                    median_gained = as.numeric(median(n_contacts_gained, na.rm = TRUE)),
                    mean_gained   = mean(n_contacts_gained,   na.rm = TRUE),
                    median_H4K20me3_l2fc = median(h4_l2fc, na.rm = TRUE),
                    frac_tad_boundary = mean(tad_boundary, na.rm = TRUE),
                    mean_repeat_content = mean(rep_cov, na.rm = TRUE)),
               if (have_feat) list(median_pc1_delta = median(pc1_delta, na.rm = TRUE),
                                   median_insulation_delta = median(insulation_delta, na.rm = TRUE))),
           by = grp]
  out[, scheme := sc][]
}), fill = TRUE)

# Each group against a repeat-content-matched draw of ALL labelled-or-not bins (not
# against `stable`: `stable` is itself a small, repeat-poor set and a stable-vs-group
# contrast would be confounded by that). Effect sizes + permutation p, see
# matched_effect().
vs_matched <- rbindlist(lapply(names(schemes), function(sc) {
  g <- factor(d[[schemes[[sc]]]], levels = grp_levels)
  mets <- c(metrics, "h4_l2fc", if (have_feat) c("pc1_delta", "insulation_delta"))
  rbindlist(lapply(grp_levels, function(gp) {
    rbindlist(lapply(mets, function(m) {
      r <- matched_effect(which(g == gp), d[[m]], d$rep_cov)
      if (is.null(r)) NULL else cbind(data.table(scheme = sc, group = gp, metric = m), r)
    }))
  }))
}))
vs_matched[, p_adj_BH := p.adjust(perm_p, "BH"), by = scheme]

# ================================================================
# 2. size_perm(): overlap-shaped questions, matched on repeat content
# ================================================================
# "Is the group enriched for the top decile of contact-loss bins?" is an overlap
# question. A naive hypergeometric/Fisher test would compare against ALL bins, which
# differ systematically in repeat density (see matched_perm above). The null here
# draws bins with the same repeat-content distribution as the group.
d[, top_lost   := n_contacts_lost   >= quantile(n_contacts_lost,   0.9, na.rm = TRUE)]
d[, top_gained := n_contacts_gained >= quantile(n_contacts_gained, 0.9, na.rm = TRUE)]
flags <- c("top_lost", "top_gained", "tad_boundary")
if (have_feat) {
  d[, switch_AB := compartment_switch == "A_B"]   # WT A -> KO B
  d[, switch_BA := compartment_switch == "B_A"]   # WT B -> KO A
  flags <- c(flags, "switch_AB", "switch_BA")
}
perm_res <- rbindlist(lapply(names(schemes), function(sc) {
  g <- factor(d[[schemes[[sc]]]], levels = grp_levels)
  rbindlist(lapply(grp_levels, function(gp) {
    idx <- which(g == gp)
    rbindlist(lapply(flags, function(fl) {
      r <- matched_perm(idx, d[[fl]], d$rep_cov)
      if (is.null(r)) return(NULL)
      cbind(data.table(scheme = sc, group = gp, flag = fl, n_bins = length(idx)), r)
    }))
  }))
}))
perm_res[, p_adj_BH := p.adjust(perm_p, "BH"), by = scheme]

# ================================================================
# 3. A/B switching: does the mark loss move chromatin between compartments?
# ================================================================
switch_tab <- NULL; cor_tab <- NULL
if (have_feat) {
  # Report BOTH counts and the within-group fraction. Counts alone are dominated by
  # group size (stable >> co_loss); fractions are what answer "does this group switch
  # more than the others". The genome-wide row is the baseline to compare against.
  switch_tab <- rbindlist(lapply(names(schemes), function(sc) {
    x <- copy(d); x[, grp := factor(get(schemes[[sc]]), levels = grp_levels)]
    x <- x[!is.na(compartment_switch)]
    one <- function(xx, lab) {
      ct <- xx[, .N, by = compartment_switch]
      ct <- merge(data.table(compartment_switch = levels(d$compartment_switch)), ct,
                  by = "compartment_switch", all.x = TRUE)
      ct[is.na(N), N := 0L]
      ct[, `:=`(group = lab, group_total = sum(N), fraction = N / sum(N))][]
    }
    out <- rbind(rbindlist(lapply(grp_levels, function(gp) one(x[grp == gp], gp))),
                 one(x, "all_bins"))
    out[, scheme := sc][]
  }))
  setcolorder(switch_tab, c("scheme", "group", "compartment_switch", "N", "group_total", "fraction"))

  # PC1 shift vs H4K20me3 log2FC, genome-wide, per bin. Spearman: both are
  # heavy-tailed and we want monotonic association, not linearity. If mark loss
  # moves chromatin A->B, deeper H4K20me3 loss should go with a more NEGATIVE
  # pc1_delta (A positive), i.e. a POSITIVE rho.
  ok <- is.finite(d$pc1_delta) & is.finite(d$h4_l2fc)
  ct <- suppressWarnings(cor.test(d$pc1_delta[ok], d$h4_l2fc[ok], method = "spearman"))
  cor_tab <- data.table(n_bins = sum(ok), spearman_rho = unname(ct$estimate), p = ct$p.value)
  # also within WT compartment: a genome-wide rho can be carried by A-vs-B baseline
  # differences in both variables, so stratify
  if ("compartment_WT" %in% names(d))
    cor_tab <- rbind(cor_tab, rbindlist(lapply(c("A", "B"), function(cp) {
      k <- ok & d$compartment_WT %in% cp
      cc <- suppressWarnings(cor.test(d$pc1_delta[k], d$h4_l2fc[k], method = "spearman"))
      data.table(n_bins = sum(k), spearman_rho = unname(cc$estimate), p = cc$p.value,
                 stratum = paste0("WT_", cp))
    })), fill = TRUE)
  cor_tab[is.na(stratum), stratum := "genome_wide"]
}

# ================================================================
# 4. Insulation change per group, and at TAD boundaries
# ================================================================
# insulation_delta = HP1gKO - WT (GAM_2; interpret the sign per GAM_2's definition of
# the insulation score). Boundaries are where a change matters most, so we split each
# group by tad_boundary; the all-bins row is the baseline a boundary should be read
# against, because boundaries differ from non-boundary bins in insulation anyway.
ins_tab <- NULL
if (have_feat && "insulation_delta" %in% names(d)) {
  ins_tab <- rbindlist(lapply(names(schemes), function(sc) {
    x <- copy(d); x[, grp := factor(get(schemes[[sc]]), levels = grp_levels)]
    x <- x[is.finite(insulation_delta) & !is.na(tad_boundary)]
    rbind(x[!is.na(grp), .(n = .N, median_insulation_delta = median(insulation_delta),
                           mean_insulation_delta = mean(insulation_delta)),
            by = .(group = as.character(grp), tad_boundary)],
          x[, .(group = "all_bins", n = .N, median_insulation_delta = median(insulation_delta),
                mean_insulation_delta = mean(insulation_delta)), by = tad_boundary])[
            , scheme := sc][]
  }))
  # group's boundary bins vs repeat-content-matched BOUNDARY bins (universe = all
  # boundary bins), and likewise for non-boundary bins: does the group change its
  # insulation more than boundaries/non-boundaries of the same repeat content do?
  ins_bd <- rbindlist(lapply(names(schemes), function(sc) {
    g <- factor(d[[schemes[[sc]]]], levels = grp_levels)
    rbindlist(lapply(grp_levels, function(gp) {
      rbindlist(lapply(c(TRUE, FALSE), function(bd) {
        r <- matched_effect(which(g == gp & d$tad_boundary %in% bd), d$insulation_delta,
                            d$rep_cov, universe = d$tad_boundary %in% bd)
        if (is.null(r)) NULL else cbind(data.table(scheme = sc, group = gp, tad_boundary = bd), r)
      }))
    }))
  }))
  ins_bd[, p_adj_BH := p.adjust(perm_p, "BH"), by = scheme]
  fwrite(ins_bd, file.path(tables_dir, "gam4_insulation_boundary_vs_matched.tsv"), sep = "\t")
}

# ================================================================
# 5. Tables
# ================================================================
fwrite(summ,      file.path(tables_dir, "gam4_group_summary.tsv"),        sep = "\t")
fwrite(vs_matched, file.path(tables_dir, "gam4_group_vs_matched_background.tsv"), sep = "\t")
fwrite(perm_res,  file.path(tables_dir, "gam4_size_perm_overlap.tsv"),    sep = "\t")
if (!is.null(switch_tab)) fwrite(switch_tab, file.path(tables_dir, "gam4_compartment_switch.tsv"), sep = "\t")
if (!is.null(cor_tab))    fwrite(cor_tab,    file.path(tables_dir, "gam4_pc1_vs_H4K20me3.tsv"),   sep = "\t")
if (!is.null(ins_tab))    fwrite(ins_tab,    file.path(tables_dir, "gam4_insulation_by_group.tsv"), sep = "\t")

cat("\n=== GAM_4 per-group summary (", PRIMARY, " = primary; abs = sensitivity) ===\n", sep = "")
print(summ, row.names = FALSE)
cat("\n=== size-matched overlap (matched on repeat content) ===\n"); print(perm_res, row.names = FALSE)
if (!is.null(switch_tab)) { cat("\n=== A/B switch table ===\n"); print(switch_tab, row.names = FALSE) }
if (!is.null(cor_tab))    { cat("\n=== PC1 shift vs H4K20me3 log2FC ===\n"); print(cor_tab, row.names = FALSE) }

# ================================================================
# Figures -> figures/compartments_vs_contacts
# ================================================================
sub_dir <- "compartments_vs_contacts"
long_c <- rbindlist(lapply(names(schemes), function(sc) {
  x <- d[, .(bin_id, grp = factor(get(schemes[[sc]]), levels = grp_levels),
             lost = n_contacts_lost, gained = n_contacts_gained)]
  x <- melt(x[!is.na(grp)], id.vars = c("bin_id", "grp"), variable.name = "direction",
            value.name = "n")
  x[, scheme := factor(sc, levels = names(schemes))][]
}))
long_c[, direction := factor(direction, levels = c("lost", "gained"),
                             labels = c("lost in HP1gKO (selWT)", "gained in HP1gKO (selcKO)"))]
g_contacts <- ggplot(long_c, aes(grp, log1p(n), fill = grp)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.15, linewidth = 0.3) +
  facet_grid(scheme ~ direction, labeller = labeller(scheme = c(
    centred = "median-centred (primary)", abs = "absolute cutoffs (sensitivity)"))) +
  scale_fill_manual(values = grp_cols, guide = "none") +
  labs(title = "Differential contact density per chromatin-mark-loss group",
       subtitle = "per 20 kb bin; the absolute-cutoff row is a sensitivity check, not the primary result",
       x = NULL, y = "log1p(number of differential contacts)") +
  theme_m
save_fig(g_contacts, "contact_density_by_group", sub_dir, width = 9, height = 6)

pp <- copy(perm_res)[, `:=`(group = factor(group, levels = grp_levels),
                            scheme = factor(scheme, levels = names(schemes)),
                            sig = ifelse(p_adj_BH < 0.05, "*", ""))]
g_perm <- ggplot(pp, aes(flag, log2(fold), fill = group)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = sig), position = position_dodge(0.8), vjust = -0.3) +
  facet_wrap(~scheme, ncol = 1) +
  scale_fill_manual(values = grp_cols, name = NULL) +
  labs(title = "Overlap vs a repeat-content-matched null (size_perm)",
       subtitle = "log2(observed / expected); * = BH < 0.05; top_* = top decile of contact counts",
       x = NULL, y = "log2 fold over matched null") +
  theme_m + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(g_perm, "size_perm_overlap", sub_dir, width = 8, height = 6)

if (!is.null(switch_tab)) {
  st <- copy(switch_tab)[, `:=`(group = factor(group, levels = c(grp_levels, "all_bins")),
                                scheme = factor(scheme, levels = names(schemes)))]
  g_sw <- ggplot(st, aes(group, fraction, fill = compartment_switch)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = ifelse(fraction > 0.04, N, "")),
              position = position_stack(vjust = 0.5), size = 3) +
    facet_wrap(~scheme, nrow = 1) +
    scale_fill_disc(name = "WT -> HP1gKO") +
    labs(title = "A/B compartment switching per group",
         subtitle = "bar = fraction of the group's bins; label = bin count (A positive = H3K4me3-high)",
         x = NULL, y = "fraction of group") +
    theme_m + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_fig(g_sw, "compartment_switch_fraction", sub_dir, width = 10, height = 4.5)

  dd <- d[is.finite(pc1_delta) & is.finite(h4_l2fc)]
  g_sc <- ggplot(dd, aes(h4_l2fc, pc1_delta)) +
    geom_bin2d(bins = 80) + scale_fill_heat0(trans = "log10", name = "bins") +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55") +
    labs(title = "PC1 shift vs H4K20me3 change, per bin",
         subtitle = sprintf("Spearman rho = %.3f (n = %s); a positive rho = deeper loss goes with a shift toward B",
                            cor_tab$spearman_rho[cor_tab$stratum == "genome_wide"],
                            format(cor_tab$n_bins[cor_tab$stratum == "genome_wide"], big.mark = ",")),
         x = "H4K20me3 log2FC (HP1gKO / WT)", y = "PC1 delta (HP1gKO - WT)") +
    theme_m
  save_fig(g_sc, "pc1_delta_vs_H4K20me3", sub_dir, width = 7, height = 5.5)

  g_pc <- ggplot(d[!is.na(get(schemes[[PRIMARY]])) & is.finite(pc1_delta)],
                 aes(factor(get(schemes[[PRIMARY]]), levels = grp_levels), pc1_delta,
                     fill = factor(get(schemes[[PRIMARY]]), levels = grp_levels))) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
    geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.15, linewidth = 0.3) +
    scale_fill_manual(values = grp_cols, guide = "none") +
    labs(title = "PC1 shift per group (median-centred groups)",
         subtitle = "delta = HP1gKO - WT; positive = toward A", x = NULL, y = "PC1 delta") +
    theme_m
  save_fig(g_pc, "pc1_delta_by_group", sub_dir, width = 5.5, height = 4.5)
}

if (!is.null(ins_tab)) {
  it <- d[!is.na(get(schemes[[PRIMARY]])) & is.finite(insulation_delta) & !is.na(tad_boundary)]
  it[, `:=`(grp = factor(get(schemes[[PRIMARY]]), levels = grp_levels),
            where = factor(ifelse(tad_boundary, "TAD boundary", "other bins"),
                           levels = c("other bins", "TAD boundary")))]
  g_ins <- ggplot(it, aes(grp, insulation_delta, fill = where)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
    geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.15, linewidth = 0.3,
                 position = position_dodge(0.8)) +
    scale_fill_disc(name = NULL) +
    labs(title = "Insulation change per group, split by TAD boundary",
         subtitle = "delta = HP1gKO - WT (median-centred groups)", x = NULL,
         y = "insulation delta") +
    theme_m
  save_fig(g_ins, "insulation_delta_by_group", sub_dir, width = 7, height = 4.5)
}

message("GAM_4 done: tables in ", tables_dir, ", figures in ", file.path(fig_dir, sub_dir))
