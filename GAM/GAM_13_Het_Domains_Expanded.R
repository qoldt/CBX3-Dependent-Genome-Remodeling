# ================================================================
# GAM_13 - Heterochromatin domains, properly constructed
# ----------------------------------------------------------------
#     cd GAM && Rscript GAM_13_Het_Domains_Expanded.R      (after GAM_1, GAM_2)
#
# WHY. GAM_10 and GAM_11 define "heterochromatin domain" from het_domain_bed by
# reduce() plus a width filter. That BED is 126,197 fragments of MEAN WIDTH
# 3.5 kb - a merged PEAK set, not a domain set - so "domains >= 100 kb" meant
# whichever single peaks happened to be long, not merged domains. The sweep
# below shows there is very little large-domain structure in it at all: at gap 0
# only 33 Mb sits in pieces >= 100 kb, against the 444 Mb it covers in total.
#
# So this stage builds domains three ways and reports every contact result
# against all three. Where they disagree, that disagreement IS the result - the
# same discipline the hotspot domains already follow (adjusted vs unadjusted
# agree on only 40% of bins).
#
#   (a) UNMERGED   reduce() only, width >= gam10_min_domain. GAM_10's definition,
#                  carried forward as the baseline so this stage is comparable.
#   (b) PEAK-MERGED gap-merged fragments at the chosen gap. Sweep reported.
#   (c) BIN-LEVEL  runs of 20 kb bins whose heterochromatin coverage exceeds
#                  gam10_het_cut, gap-merged on the bin grid. This is the
#                  cleanest of the three because it lives on the same grid as
#                  the contacts and the ChIP signal, so nothing has to be
#                  re-overlapped to use it.
#
# The contact statistics mirror GAM_10/GAM_11 exactly - per-chromosome means,
# bootstrap over the 19 autosomes, difference-in-differences - because that is
# the framing that survived; a CI excluding zero means two classes changed by
# different amounts. Bins are strongly autocorrelated, so the chromosome, not
# the bin, is the resampling unit.
# ================================================================
source("config.R")
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(GenomicRanges); library(IRanges)
})

set.seed(1)
N_BOOT  <- 2000L
BANDS   <- list(c(2e4, 1e6), c(1e6, 5e6), c(5e6, 1e7), c(1e7, Inf))
BANDNM  <- c("0.02-1Mb", "1-5Mb", "5-10Mb", ">10Mb")
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
  labs(title = "Heterochromatin fragments percolate when gap-merged",
       subtitle = paste0("126,197 peak fragments, mean width 3.5 kb. Past ~20 kb the gap exceeds the ",
                         "typical spacing\nbetween peaks and they fuse into chromosome-length blocks. ",
                         "Dashed line = the gap used downstream."),
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
       subtitle = paste0("Percent change in mean NPMI, HP1gKO vs WT, by separation. Bars are 95% CIs from a ",
                         "bootstrap over the 19 autosomes.\nThe three panels are three ways of defining a ",
                         "domain from the same peak set - where they disagree, the definition is doing the work."),
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
         subtitle = paste0("cLAD/ciLAD (Peric-Hupkes 2010) are external to both genotypes and to all five ",
                           "marks, so unlike our own PC1\nthey cannot be circular. They are the ",
                           "cell-type-invariant core, NOT a neuronal annotation."),
         x = NULL, y = "% of domain bp") +
    theme_m
  save_fig(p, "domains_vs_lads", subdir = fig_sub, width = 7.5, height = 5)
}

msg("\nGAM_13 complete.")
