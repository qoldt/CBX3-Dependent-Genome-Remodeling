# ================================================================
# GAM_5 - Heterochromatin domains and their contacts
# ----------------------------------------------------------------
#     cd GAM && Rscript GAM_5_Heterochromatin.R      (after GAM_1, GAM_2)
#
# Merges the three stages that all ask "what happens to contacts involving
# heterochromatin", all define their classes from the same merged
# H3K9me3/H4K20me3 peak set, and all use the same statistics: per-chromosome
# means, bootstrap over the 19 autosomes, difference-in-differences. Previously
# GAM_10, GAM_11 and GAM_13.
#
# WHAT EACH SECTION ANSWERS, and how they relate:
#
#   A (was GAM_10)  Do heterochromatin domains lose contacts WITHIN themselves
#                   or BETWEEN each other? Classes intra_dom / inter_het / eu-eu.
#   B (was GAM_11)  Does it matter whether the domain LOST the mark in the KO?
#                   Heterochromatin split into tertiles by H4K20me3 log2FC.
#   C (was GAM_13)  Section A again, but on domains that are actually domains.
#                   The het_domain_bed is 126,197 fragments of mean width 3.5 kb
#                   - a peak set - so A's "large domains" were whichever single
#                   fragments happened to be long. C builds domains three ways
#                   and runs A's contrast on all three.
#
# SO SECTION C SUPERSEDES SECTION A's intra-vs-inter contrast, and A is kept
# because it also carries the eu-eu baseline and the per-chromosome sums that
# nothing else writes. Where they disagree, trust C. They currently agree: the
# intra/inter difference is indistinguishable from zero under every definition.
#
# NOT YET UNIFIED, deliberately. Sections A and B still build their own domains
# from the unmerged fragments, so this file is a structural merge with
# byte-identical outputs rather than a rewrite. Pointing A and B at C's merged
# domains would change their numbers - an improvement, but a separate change
# that should be made and reported on its own.
#
# Output file names are unchanged (gam10_*, gam11_*, gam13_*) so every existing
# consumer, the report generator included, keeps working.
# ================================================================
source("config.R")
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(GenomicRanges); library(IRanges)
})



# ==============================================================================
# SECTION A - intra-domain vs inter-domain contacts            (was GAM_10)
# ==============================================================================

local({
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
  # The two-level chromatin partition, from config.R so that every stage sees the
  # same one. It supersedes the old is_het/eu_eu binary, which called the 85% of
  # the genome that merely lacks an H3K9me3/H4K20me3 peak "euchromatin" and put
  # Polycomb in with it. het_frac below is the same peak coverage as before and
  # is kept because the H4K20me3 strata in section B still use it.
  cc <- gam_chromatin_classes(a)
  a <- merge(a, cc, by = "bin_id", sort = FALSE)
  setorder(a, chr, start)
  a[, is_het := class2 == "H3K9me3" & !is.na(class2)]
  stopifnot(sum(a$is_het) > 1000)   # a seqname-style failure would zero this
  cat(sprintf("=== chromatin classes ===\n  level 1: %s\n  level 2: %s\n",
      paste(sprintf("%s %d", levels(a$class1), table(a$class1)), collapse = "  "),
      paste(sprintf("%s %d", levels(a$class2), table(a$class2)), collapse = "  ")))
  cat(sprintf("  bins with no ChromHMM coverage, excluded from every pair: %d\n",
              sum(is.na(a$class2))))
  # The partition as an itemRgb BED, for loading beside the contact maps in IGV.
  local({
    bp <- file.path(bed_dir, "chromatin_states_level2.bed")
    r <- write_chromatin_state_bed(as.data.table(x$annot)[chr %in% gam_autosomes], bp)
    cat(sprintf("  wrote %s (%d runs from %d bins)\n", bp, nrow(r), nrow(a)))
  })

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

  BAND_EDGES <- gam_band_edges   # config.R - one definition for every stage
  BAND_NAMES <- gam_band_names

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
      k1 <- as.integer(ids$class1[k]); k2 <- as.integer(ids$class2[k])
      dom <- ids$dom_id[k]; pos <- ids$start[k]
      # Work on the upper-triangle INDEX rather than on n x n matrices. The old
      # version built four n x n logicals plus an n x n character matrix; on
      # chr1 that is five 9,758^2 objects at once, and the character one alone
      # is ~760 MB of pointers. Indexing costs two integer vectors instead.
      iu <- which(upper.tri(m), arr.ind = TRUE)
      i1 <- iu[, 1]; i2 <- iu[, 2]; rm(iu)
      val  <- m[cbind(i1, i2)]
      band <- cut(abs(pos[i2] - pos[i1]), BAND_EDGES, labels = BAND_NAMES, right = FALSE)
      # A bin with no ChromHMM coverage has no class, and a pair involving one
      # is dropped from BOTH levels rather than being given a default - there
      # are 3,218 such bins and defaulting them would invent a class.
      ok <- is.finite(val) & !is.na(band) & !is.na(k1[i1]) & !is.na(k1[i2])
      if (!any(ok)) { rm(m, i1, i2, val, band, ok); gc(FALSE); next }
      i1 <- i1[ok]; i2 <- i2[ok]; val <- val[ok]; band <- as.character(band[ok]); rm(ok)

      # Aggregate on INTEGER pair codes and attach the labels afterwards, on the
      # few hundred surviving rows. Pasting labels across ~47M pairs first would
      # allocate a 47M-element character vector per level for no gain.
      code <- function(z) pmin(z[i1], z[i2]) * 100L + pmax(z[i1], z[i2])
      agg <- function(cd) data.table(cd = cd, band = band, v = val)[
                            , .(s = sum(v), n = .N), by = .(cd, band)]
      lab <- function(d, levs) {
        d[, pair := paste(levs[cd %/% 100L], levs[cd %% 100L], sep = "-")][, cd := NULL][]
      }
      A1 <- lab(agg(code(k1)), gam_class1_levels)[, level := "L1"]
      A2 <- lab(agg(code(k2)), gam_class2_levels)[, level := "L2"]
      # The same-domain / different-domain split, which is only meaningful for a
      # pair of H3K9me3 bins and is the question section A exists to ask.
      kh <- which(gam_class2_levels == "H3K9me3")
      hh <- k2[i1] == kh & k2[i2] == kh
      A3 <- if (any(hh)) {
        sd <- !is.na(dom[i1][hh]) & !is.na(dom[i2][hh]) & dom[i1][hh] == dom[i2][hh]
        data.table(pair = fifelse(sd, "same domain", "different domains"),
                   band = band[hh], v = val[hh])[, .(s = sum(v), n = .N), by = .(pair, band)][
                   , level := "dom"][]
      } else NULL
      acc[[length(acc) + 1]] <- rbindlist(list(A1, A2, A3), use.names = TRUE)[
                                  , `:=`(chr = ch, genotype = g)][]
      rm(m, i1, i2, val, band, A1, A2, A3, hh); gc(FALSE)
    }
    cat("  ", ch, "done\n")
  }
  A <- rbindlist(acc); rm(acc); gc(FALSE)
  fwrite(A, file.path(tables_dir, "gam10_perchrom_sums.tsv"), sep = "\t")

  # --- pooled means, percent change, and a chromosome bootstrap ----------------
  pool <- function(d) d[, .(mean_npmi = sum(s) / sum(n), n = sum(n)),
                        by = .(genotype, level, pair, band)]
  P <- pool(A)
  W <- dcast(P, level + pair + band ~ genotype, value.var = c("mean_npmi", "n"))
  setnames(W, c("mean_npmi_WT","mean_npmi_HP1gKO"), c("WT","KO"))
  W[, pct := 100 * (KO - WT) / WT]

  chrs <- unique(A$chr)
  set.seed(42)
  boot <- rbindlist(lapply(seq_len(gam10_nboot), function(i) {
    smp <- sample(chrs, length(chrs), replace = TRUE)
    d <- rbindlist(lapply(seq_along(smp), function(j) { z <- A[chr == smp[j]]; z[, chr := j]; z }))
    p <- pool(d)
    w <- dcast(p, level + pair + band ~ genotype, value.var = "mean_npmi")
    w[, .(level, pair, band, pct = 100 * (HP1gKO - WT) / WT, rep = i)]
  }))
  ci <- boot[, .(lo = quantile(pct, .025, na.rm = TRUE), hi = quantile(pct, .975, na.rm = TRUE)),
             by = .(level, pair, band)]
  W <- merge(W, ci, by = c("level","pair","band"))
  setorder(W, level, band, pair)
  fwrite(W, file.path(tables_dir, "gam10_class_change.tsv"), sep = "\t")

  cat("\n=== mean NPMI change by contact class, within distance band ===\n")
  cat("   (95% CI from a bootstrap over the 19 autosomes - n is 19, not the pair count)\n")
  cat("   L1 = euchromatin/heterochromatin; L2 = the six sub-classes.\n")
  for (lv in c("L1", "L2")) for (bn in BAND_NAMES) {
    d <- W[level == lv & band == bn]; if (!nrow(d)) next
    setorder(d, -pct)
    cat("\n--", lv, bn, "--\n")
    cat(sprintf("  %-42s %8s %8s %9s %18s %14s\n","pair","WT","KO","change","95% CI","pairs"))
    for (i in seq_len(nrow(d)))
      cat(sprintf("  %-42s %8.4f %8.4f %+8.1f%% [%+6.1f%%,%+6.1f%%] %14s\n",
          d$pair[i], d$WT[i], d$KO[i], d$pct[i], d$lo[i], d$hi[i],
          format(d$n_WT[i], big.mark = ",")))
  }

  # --- the contrast the question asks for: intra vs inter, same band -----------
  # Now carried by the "dom" level, which splits H3K9me3-H3K9me3 pairs only.
  did <- boot[level == "dom"]
  did <- dcast(did, band + rep ~ pair, value.var = "pct")
  did[, diff := `same domain` - `different domains`]
  D <- did[, .(est = mean(diff, na.rm = TRUE),
               lo = quantile(diff, .025, na.rm = TRUE),
               hi = quantile(diff, .975, na.rm = TRUE),
               p_two_sided = 2 * min(mean(diff <= 0, na.rm = TRUE),
                                     mean(diff >= 0, na.rm = TRUE))), by = band]
  D <- merge(D, W[level == "dom" & pair == "same domain",        .(band, intra = pct)], by = "band")
  D <- merge(D, W[level == "dom" & pair == "different domains",   .(band, inter = pct)], by = "band")
  D[, band := factor(band, BAND_NAMES)]; setorder(D, band)
  fwrite(D, file.path(tables_dir, "gam10_intra_vs_inter.tsv"), sep = "\t")

  cat("\n=== WITHIN-domain vs BETWEEN-heterochromatin, difference-in-differences ===\n")
  cat(sprintf("  %-10s %10s %10s %10s %18s %8s\n","band","intra","inter","intra-inter","95% CI","p"))
  for (i in seq_len(nrow(D)))
    cat(sprintf("  %-10s %+9.1f%% %+9.1f%% %+10.1f%% [%+6.1f%%,%+6.1f%%] %8.3f\n",
        D$band[i], D$intra[i], D$inter[i], D$est[i], D$lo[i], D$hi[i], D$p_two_sided[i]))
  cat("\n  A CI excluding 0 means the two classes changed by different amounts.\n")
  cat("  p is the bootstrap two-sided proportion, not a parametric test.\n")

  # LEVEL 1 only. The six-way level-2 version is 21 pairs and belongs in the
  # table, not in a bar chart; GAM_7 draws the level-2 view per chromosome.
  pw <- W[level == "L1"]
  pw[, band := factor(band, BAND_NAMES)]
  pw[, pair := factor(pair, c("heterochromatin-heterochromatin",
                              "euchromatin-heterochromatin", "euchromatin-euchromatin"))]
  pw <- pw[!is.na(pair)]
  p <- ggplot(pw, aes(band, pct, fill = pair)) +
    geom_col(position = position_dodge(.8), width = .72) +
    geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(.8),
                  width = .22, linewidth = .4) +
    geom_hline(yintercept = 0, linewidth = .3) +
    scale_fill_disc(name = NULL) +
    labs(title = "Contact change by chromatin class and separation",
         subtitle = paste0("Change in mean NPMI, HP1gKO vs WT. Heterochromatin here is Polycomb + ",
                           "H3K9me3 + quiescent;\neuchromatin is TSS + enhancer + transcribed. ",
                           "Bars are 95% CIs over the 19 autosomes."),
         x = "contact distance", y = "change in mean NPMI (%)") +
    theme_m + theme(legend.position = "top", text = element_text(size = 14))
  save_fig(p, "heterochromatin_contact_classes", subdir = "heterochromatin", width = 11, height = 7)

  # --- LEVEL 2, the six self-pairs --------------------------------------------
  # The level-1 figure above shows heterochromatin beating euchromatin. It
  # cannot show the finding that matters most, which is that the three
  # heterochromatin sub-classes do NOT behave alike: H3K9me3 keeps gaining out
  # to 10 Mb while Polycomb and quiescent rise only below 0.2 Mb and then fall
  # with everything else. Only the self-pairs are drawn - the 15 mixed pairs are
  # in gam10_class_change.tsv and would make this unreadable.
  sp <- paste(gam_class2_levels, gam_class2_levels, sep = "-")
  pw2 <- W[level == "L2" & pair %in% sp]
  if (nrow(pw2)) {
    pw2[, band := factor(band, BAND_NAMES)]
    pw2[, cls := factor(sub("-.*$", "", pair), gam_class2_levels)]
    # DODGED ERROR BARS, NOT RIBBONS. Six translucent ribbons over six bands
    # overlap into a single wash that hides the lines they belong to, and the
    # widest one (TSS) stretched the y axis to -40 so every other class was
    # squeezed into the top third. Dodging separates the intervals and lets the
    # axis follow the estimates.
    dg <- position_dodge(width = 0.5)
    p2 <- ggplot(pw2, aes(band, pct, colour = cls, group = cls)) +
      geom_hline(yintercept = 0, linewidth = .3, colour = "grey50") +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = .45,
                    alpha = .45, position = dg) +
      geom_line(linewidth = .8, position = dg) +
      geom_point(size = 2.3, position = dg) +
      scale_colour_disc(name = NULL) +
      labs(title = "Only H3K9me3 heterochromatin keeps gaining contact past 1 Mb",
           subtitle = paste0("Self-pairs only (both bins in the same class), HP1gKO vs WT. ",
                             "Bars are 95% bootstrap CIs over the 19 autosomes.\n",
                             "Polycomb and quiescent chromatin rise below 0.2 Mb and then fall ",
                             "with euchromatin - the gain is not repression in general."),
           x = "distance between the two bins (Mb = megabases)",
           y = "change in mean NPMI (%)") +
      theme_m + theme(legend.position = "top", text = element_text(size = 13))
    save_fig(p2, "contact_classes_level2", subdir = "heterochromatin",
             width = 11, height = 7)
  }

  # --- LEVEL 2, the heterochromatin pairs INCLUDING the mixed ones ------------
  # The self-pair figure above cannot answer the question it raises. "H3K9me3
  # gains" has two readings - H3K9me3 bins find each OTHER, or H3K9me3 bins gain
  # contact with whatever is nearby - and they are told apart only by the MIXED
  # pair. H3K9me3-quiescent is the discriminator, and it was in the table and in
  # no figure.
  #
  # Six lines, not 21: the three heterochromatin classes crossed with themselves
  # and each other. That is the comparison the conclusion rests on and it stays
  # readable; the full 21 are in gam10_class_change.tsv.
  # PAIR NAMES ARE ORDERED BY LEVEL, NOT ALPHABETICALLY. The table builds a pair
  # as lo-hi on gam_class2_levels, so it says "Polycomb-H3K9me3"; sort() would
  # ask for "H3K9me3-Polycomb" and match nothing, and the figure would silently
  # lose a line rather than fail.
  hetc <- c("Polycomb", "H3K9me3", "Quiescent")
  hp   <- unique(unlist(lapply(hetc, function(x) vapply(hetc, function(y) {
            k <- sort(match(c(x, y), gam_class2_levels))
            paste(gam_class2_levels[k], collapse = "-") }, ""))))
  pw3  <- W[level == "L2" & pair %in% hp]
  if (nrow(pw3)) {
    pw3[, band := factor(band, BAND_NAMES)]
    # Self-pairs solid, mixed pairs dashed, so the two kinds are separable
    # without reading six legend entries against six colours.
    pw3[, mixed := sub("-.*$", "", pair) != sub("^.*-", "", pair)]
    pw3[, pair := factor(pair, hp)]
    dg3 <- position_dodge(width = 0.5)
    p3 <- ggplot(pw3, aes(band, pct, colour = pair, group = pair)) +
      geom_hline(yintercept = 0, linewidth = .3, colour = "grey50") +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = .45,
                    alpha = .45, position = dg3) +
      geom_line(aes(linetype = mixed), linewidth = .8, position = dg3) +
      geom_point(size = 2.3, position = dg3) +
      scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"),
                            guide = "none") +
      scale_colour_disc(name = NULL) +
      labs(title = "H3K9me3 gains contact with itself, not with the quiescent chromatin around it",
           subtitle = paste(strwrap(paste0(
             "The three heterochromatin classes, self-pairs (solid) and mixed pairs (dashed), ",
             "HP1gKO vs WT. Bars are 95% bootstrap CIs over the 19 autosomes. Past 0.2 Mb the ",
             "three separate: H3K9me3-H3K9me3 keeps rising, H3K9me3-quiescent is flat and ",
             "quiescent-quiescent falls - so the gain is specific to H3K9me3 self-association ",
             "and does not spread to its neighbours. Quiescent chromatin is 69% of the genome ",
             "and is what the LAD- and LINE-rich domains that LOSE contacts are made of."),
             width = 110), collapse = "\n"),
           x = "distance between the two bins (Mb = megabases)",
           y = "change in mean NPMI (%)") +
      theme_m + theme(legend.position = "top", text = element_text(size = 13))
    save_fig(p3, "contact_classes_level2_mixed", subdir = "heterochromatin",
             width = 11, height = 7.5)
  }
  cat("\nGAM_10 complete.\n")
})


# ==============================================================================
# SECTION B - does losing the mark decouple the domain?        (was GAM_11)
# ==============================================================================

local({
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

  BAND_EDGES <- gam_band_edges   # config.R - one definition for every stage
  BAND_NAMES <- gam_band_names
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
})


# ==============================================================================
# SECTION C - domains built properly, and A's contrast redone  (was GAM_13)
# ==============================================================================

local({
  set.seed(1)
  N_BOOT  <- 2000L
  # Keeps its own near-diagonal band (one 0.02-1Mb bin, deliberately excluding
  # the diagonal) but takes the same LONG-RANGE split as everything else: the
  # old open-ended >10Mb had nearly half its pairs beyond 50 Mb. See
  # gam_band_edges in config.R.
  BANDS   <- list(c(2e4, 1e6), c(1e6, 5e6), c(5e6, 1e7), c(1e7, 5e7), c(5e7, Inf))
  BANDNM  <- c("0.02-1Mb", "1-5Mb", "5-10Mb", "10-50Mb", ">50Mb")
  SZ_CUT  <- c(0, 2e5, 5e5, Inf)                 # domain size tertile-ish classes
  SZ_NM   <- c("100-200kb", "200-500kb", ">500kb")

  fig_sub <- "het_domains"
  msg <- function(...) cat(sprintf(...), "\n", sep = "")

  # --- inputs -------------------------------------------------------------------
  X <- readRDS(file.path(rds_dir, "gam_bins_annotated.rds"))
  A <- as.data.table(as.data.frame(X$annot)); setorder(A, chr, start)
  A[, h4 := X$l2fc[match(bin_id, rownames(X$l2fc)), "H4K20me3"]]
  A <- A[chr %in% gam_autosomes]
  bins <- GRanges(A$chr, IRanges(A$start + 1, A$end))

  # ==============================================================================
  # 1. The sweep - why the chosen gap, and what it costs
  # ==============================================================================
  msg("=== gap-merge sweep ===")
  gl <- sum(as.numeric(A[, .(w = max(end)), by = chr]$w))
  sweep <- rbindlist(lapply(c(0L, 2000L, 5000L, 10000L, 20000L, 50000L), function(gp)
    rbindlist(lapply(c(50000L, 100000L, 200000L), function(mw) {
      g <- merge_het_domains(gp, mw)
      data.table(gap = gp, min_width = mw, n = length(g),
                 total_Mb = sum(as.numeric(width(g))) / 1e6,
                 pct_genome = 100 * sum(as.numeric(width(g))) / gl,
                 median_kb = as.numeric(median(width(g))) / 1e3)
    }))))
  fwrite(sweep, file.path(tables_dir, "gam13_merge_sweep.tsv"), sep = "\t")
  print(sweep[min_width == 100000L])

  p <- ggplot(sweep, aes(gap / 1e3, pct_genome, colour = factor(min_width / 1e3))) +
    geom_line() + geom_point(size = 1.4) +
    geom_vline(xintercept = gam_het_merge$gap / 1e3, linetype = 2, colour = "grey40") +
    scale_colour_manual(values = disc_pal(3), name = "min width (kb)") +
    # The input is 126,197 peak fragments of mean width 3.5 kb, which is why the
    # percolation point is so early: past ~20 kb the gap exceeds the typical
    # spacing between peaks.
    labs(title = "Heterochromatin fragments percolate when gap-merged",
         subtitle = paste0("Dashed = the gap used downstream. Past ~20 kb the fragments fuse into ",
                           "chromosome-length blocks."),
         x = "merge gap (kb)", y = "% of autosomal genome in domains") +
    theme_m
  save_fig(p, "merge_sweep", subdir = fig_sub, width = 8, height = 5)

  # ==============================================================================
  # 2. The three domain definitions
  # ==============================================================================
  # Per-bin heterochromatin coverage from the raw fragments, used for the
  # bin-level definition below and as a covariate.
  het_raw <- {
    dm <- fread(het_domain_bed, header = FALSE, select = 1:3,
                col.names = c("chr", "start", "end"))
    dm[, chr := ifelse(grepl("^chr", chr), chr, paste0("chr", chr))]
    dm <- dm[chr %in% gam_autosomes]
    reduce(GRanges(dm$chr, IRanges(dm$start + 1, dm$end)), ignore.strand = TRUE)
  }
  A[, het_frac := {
    h <- findOverlaps(bins, het_raw, ignore.strand = TRUE)
    w <- width(pintersect(bins[queryHits(h)], het_raw[subjectHits(h)], ignore.strand = TRUE))
    v <- numeric(length(bins)); agg <- tapply(w, queryHits(h), sum)
    v[as.integer(names(agg))] <- as.numeric(agg)
    pmin(v / width(bins), 1)
  }]
  stopifnot(sum(A$het_frac > 0) > 10000)   # a seqname-style failure would zero this

  doms <- list(
    unmerged    = merge_het_domains(0L, gam_het_merge$min_width),
    peak_merged = merge_het_domains(gam_het_merge$gap, gam_het_merge$min_width),
    bin_level   = {
      # runs of heterochromatic bins, allowing a one-bin gap, then a width filter
      A[, is_het := het_frac >= 0.5]
      r <- A[, {
        v <- is_het
        # bridge single-bin holes so one dropout does not split a domain
        v2 <- v | (data.table::shift(v, 1L, FALSE) & data.table::shift(v, -1L, FALSE))
        rl <- rle(v2); e <- cumsum(rl$lengths); s <- e - rl$lengths + 1L
        k <- which(rl$values)
        .(start = start[s[k]], end = end[e[k]])
      }, by = chr]
      g <- GRanges(r$chr, IRanges(r$start + 1, r$end))
      g[width(g) >= gam_het_merge$min_width]
    }
  )
  dom_tab <- rbindlist(lapply(names(doms), function(n) data.table(
    definition = n, n = length(doms[[n]]),
    total_Mb = sum(as.numeric(width(doms[[n]]))) / 1e6,
    median_kb = as.numeric(median(width(doms[[n]]))) / 1e3,
    max_Mb = max(as.numeric(width(doms[[n]]))) / 1e6)))
  fwrite(dom_tab, file.path(tables_dir, "gam13_domain_definitions.tsv"), sep = "\t")
  msg("\n=== domain definitions ==="); print(dom_tab)

  # How much do they agree? Jaccard on bin membership, the same measure used for
  # the hotspot domains, so the two are directly comparable.
  memb <- lapply(doms, function(g) overlapsAny(bins, g, ignore.strand = TRUE))
  jac <- CJ(a = names(doms), b = names(doms))[a < b]
  jac[, jaccard := mapply(function(x, y) sum(memb[[x]] & memb[[y]]) / sum(memb[[x]] | memb[[y]]), a, b)]
  fwrite(jac, file.path(tables_dir, "gam13_definition_agreement.tsv"), sep = "\t")
  msg("\n=== agreement between definitions (Jaccard on bins) ==="); print(jac)

  # ==============================================================================
  # 3. Contact statistics per definition
  # ==============================================================================
  # Bootstrap over chromosomes: resample the 19 autosomes, recompute the
  # chromosome-mean statistic. The honest n is 19, not the pair count.
  boot_ci <- function(per_chr_vals, per_chr_n, B = N_BOOT) {
    ch <- names(per_chr_vals)
    est <- sum(per_chr_vals * per_chr_n, na.rm = TRUE) / sum(per_chr_n, na.rm = TRUE)
    bs <- replicate(B, {
      i <- sample(ch, length(ch), replace = TRUE)
      sum(per_chr_vals[i] * per_chr_n[i], na.rm = TRUE) / sum(per_chr_n[i], na.rm = TRUE)
    })
    c(est = est, lo = unname(quantile(bs, .025, na.rm = TRUE)),
      hi = unname(quantile(bs, .975, na.rm = TRUE)))
  }

  band_of <- function(d) {
    out <- rep(NA_integer_, length(d))
    for (i in seq_along(BANDS))
      out[is.na(out) & d >= BANDS[[i]][1] & d < BANDS[[i]][2]] <- i
    out
  }

  res_intra_inter <- list(); res_size <- list(); res_net <- list()

  for (dn in names(doms)) {
    g <- doms[[dn]]
    msg("\n=== %s: %d domains ===", dn, length(g))
    # map each bin to a domain id (widest overlap), and the domain's size class
    h <- findOverlaps(bins, g, ignore.strand = TRUE)
    ov <- width(pintersect(bins[queryHits(h)], g[subjectHits(h)], ignore.strand = TRUE))
    best <- data.table(q = queryHits(h), s = subjectHits(h), w = ov)
    setorder(best, q, -w); best <- best[!duplicated(q)]
    A[, dom := NA_integer_][best$q, dom := best$s]
    A[, dom_w := NA_real_][!is.na(dom), dom_w := width(g)[dom]]
    A[, dom_sz := cut(dom_w, SZ_CUT, labels = SZ_NM, right = FALSE)]

    per_chr <- list()
    for (ch in gam_autosomes) {
      idx <- A[chr == ch, which = TRUE]
      if (length(idx) < 50) next
      M <- tryCatch(list(WT = load_npmi("WT", ch), KO = load_npmi("HP1gKO", ch)),
                    error = function(e) NULL)
      if (is.null(M)) { msg("  %s: no matrix, skipped", ch); next }
      ids <- intersect(A$bin_id[idx], rownames(M$WT))
      ids <- intersect(ids, rownames(M$KO))
      if (length(ids) < 50) next
      sub <- A[match(ids, bin_id)]
      W <- M$WT[ids, ids]; K <- M$KO[ids, ids]
      pos <- sub$start
      n <- length(ids)
      iu <- which(upper.tri(W), arr.ind = TRUE)
      d  <- abs(pos[iu[, 2]] - pos[iu[, 1]])
      bd <- band_of(d)
      wv <- W[iu]; kv <- K[iu]
      di <- sub$dom[iu[, 1]]; dj <- sub$dom[iu[, 2]]
      cls <- fifelse(is.na(di) | is.na(dj), "other",
              fifelse(di == dj, "intra", "inter"))
      ok <- is.finite(wv) & is.finite(kv) & !is.na(bd)
      dt <- data.table(band = BANDNM[bd[ok]], cls = cls[ok], wv = wv[ok], kv = kv[ok],
                       szi = as.character(sub$dom_sz[iu[, 1]])[ok],
                       szj = as.character(sub$dom_sz[iu[, 2]])[ok])
      per_chr[[ch]] <- list(
        ii = dt[cls %in% c("intra", "inter"),
                .(WT = mean(wv), KO = mean(kv), n = .N), by = .(band, cls)],
        sz = dt[cls == "intra" & !is.na(szi),
                .(WT = mean(wv), KO = mean(kv), n = .N), by = .(band, sz = szi)])
      rm(M, W, K); gc(FALSE)
      msg("  %s done (%d bins)", ch, length(ids))
    }
    if (!length(per_chr)) next

    # --- intra vs inter, per band ---
    II <- rbindlist(lapply(names(per_chr), function(ch) {
      x <- per_chr[[ch]]$ii; x[, chr := ch][] }))
    out <- rbindlist(lapply(BANDNM, function(b) rbindlist(lapply(c("intra", "inter"), function(cl) {
      x <- II[band == b & cls == cl]
      if (!nrow(x)) return(NULL)
      v <- setNames((x$KO - x$WT) / x$WT * 100, x$chr); nn <- setNames(x$n, x$chr)
      ci <- boot_ci(v, nn)
      data.table(definition = dn, band = b, cls = cl, pct = ci["est"],
                 lo = ci["lo"], hi = ci["hi"], n_pairs = sum(x$n))
    }))))
    res_intra_inter[[dn]] <- out

    # --- difference-in-differences: intra minus inter ---
    dd <- rbindlist(lapply(BANDNM, function(b) {
      xi <- II[band == b & cls == "intra"]; xe <- II[band == b & cls == "inter"]
      ch <- intersect(xi$chr, xe$chr); if (length(ch) < 5) return(NULL)
      f <- function(i) {
        a <- xi[match(i, chr)]; e <- xe[match(i, chr)]
        wi <- sum(a$n); we <- sum(e$n)
        (sum((a$KO - a$WT) / a$WT * 100 * a$n) / wi) - (sum((e$KO - e$WT) / e$WT * 100 * e$n) / we)
      }
      est <- f(ch)
      bs <- replicate(N_BOOT, f(sample(ch, length(ch), replace = TRUE)))
      data.table(definition = dn, band = b, intra_minus_inter = est,
                 lo = quantile(bs, .025), hi = quantile(bs, .975),
                 p_two_sided = 2 * min(mean(bs <= 0), mean(bs >= 0)))
    }))
    res_net[[dn]] <- dd

    # --- by domain size ---
    SZ <- rbindlist(lapply(names(per_chr), function(ch) {
      x <- per_chr[[ch]]$sz; if (!nrow(x)) return(NULL); x[, chr := ch][] }))
    if (nrow(SZ)) {
      sz <- rbindlist(lapply(BANDNM, function(b) rbindlist(lapply(SZ_NM, function(s) {
        x <- SZ[band == b & sz == s]; if (nrow(x) < 5) return(NULL)
        v <- setNames((x$KO - x$WT) / x$WT * 100, x$chr); nn <- setNames(x$n, x$chr)
        ci <- boot_ci(v, nn)
        data.table(definition = dn, band = b, size = s, pct = ci["est"],
                   lo = ci["lo"], hi = ci["hi"], n_pairs = sum(x$n))
      }))))
      res_size[[dn]] <- sz
    }
  }

  II <- rbindlist(res_intra_inter); fwrite(II, file.path(tables_dir, "gam13_intra_vs_inter.tsv"), sep = "\t")
  DD <- rbindlist(res_net);        fwrite(DD, file.path(tables_dir, "gam13_did.tsv"), sep = "\t")
  SZ <- rbindlist(res_size);       fwrite(SZ, file.path(tables_dir, "gam13_by_domain_size.tsv"), sep = "\t")
  msg("\n=== intra vs inter (all definitions) ==="); print(II)
  msg("\n=== difference-in-differences ==="); print(DD)

  II[, band := factor(band, levels = BANDNM)]
  p <- ggplot(II, aes(band, pct, colour = cls, group = cls)) +
    geom_hline(yintercept = 0, linewidth = .3, colour = "grey50") +
    geom_point(position = position_dodge(.4), size = 1.8) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = .18, position = position_dodge(.4)) +
    facet_wrap(~ definition, nrow = 1) +
    scale_colour_manual(values = c(intra = "#b2182b", inter = "#2166ac"), name = NULL) +
    labs(title = "Heterochromatin contacts: within a domain vs between domains",
         subtitle = paste0("Percent change in mean NPMI, HP1gKO vs WT, by separation; bars are 95% CIs over the 19 ",
                           "autosomes.\nThe three panels are three definitions of a domain from the SAME peak set - ",
                           "where they disagree, the definition is doing the work."),
         x = "genomic separation", y = "% change in mean NPMI") +
    theme_m
  save_fig(p, "intra_vs_inter_by_definition", subdir = fig_sub, width = 11, height = 5)

  if (nrow(SZ)) {
    SZ[, band := factor(band, levels = BANDNM)][, size := factor(size, levels = SZ_NM)]
    p <- ggplot(SZ, aes(band, pct, colour = size, group = size)) +
      geom_hline(yintercept = 0, linewidth = .3, colour = "grey50") +
      geom_point(position = position_dodge(.5), size = 1.8) +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = .18, position = position_dodge(.5)) +
      facet_wrap(~ definition, nrow = 1) +
      scale_colour_manual(values = disc_pal(3), name = "domain size") +
      labs(title = "Within-domain contact change by domain size",
           subtitle = "Percent change in mean NPMI, HP1gKO vs WT. Bars are 95% CIs over the 19 autosomes.",
           x = "genomic separation", y = "% change in mean NPMI") +
      theme_m
    save_fig(p, "by_domain_size", subdir = fig_sub, width = 11, height = 5)
  }

  # ==============================================================================
  # 4. The external axis: how do these domains sit against constitutive LADs?
  # ==============================================================================
  lads <- tryCatch(load_lads(), error = function(e) NULL)
  if (!is.null(lads)) {
    lad_tab <- rbindlist(lapply(names(doms), function(dn) {
      g <- doms[[dn]]
      rbindlist(lapply(c("cLAD", "ciLAD"), function(cl) {
        L <- reduce(lads[lads$class == cl], ignore.strand = TRUE)
        h <- findOverlaps(g, L, ignore.strand = TRUE)
        w <- width(pintersect(g[queryHits(h)], L[subjectHits(h)], ignore.strand = TRUE))
        data.table(definition = dn, lad_class = cl,
                   pct_domain_bp = 100 * sum(as.numeric(w)) / sum(as.numeric(width(g))))
      }))
    }))
    fwrite(lad_tab, file.path(tables_dir, "gam13_domains_vs_lads.tsv"), sep = "\t")
    msg("\n=== domain bp inside constitutive LADs ==="); print(lad_tab)

    p <- ggplot(lad_tab, aes(definition, pct_domain_bp, fill = lad_class)) +
      geom_col(position = position_dodge(.8), width = .7) +
      scale_fill_manual(values = c(cLAD = "#53426b", ciLAD = "#cbbddc"), name = NULL) +
      labs(title = "Heterochromatin domains against constitutive LADs",
           # cLAD/ciLAD are external to both genotypes and to all five marks, so
           # unlike our own PC1 they cannot be circular - that is why they are here.
           subtitle = paste0("cLAD/ciLAD (Peric-Hupkes 2010) are the cell-type-invariant core, NOT a neuronal ",
                             "annotation."),
           x = NULL, y = "% of domain bp") +
      theme_m
    save_fig(p, "domains_vs_lads", subdir = fig_sub, width = 7.5, height = 5)
  }

  msg("\nGAM_13 complete.")
})
