# ================================================================
# verify_GAM.R - is the analysis internally consistent?
#   cd GAM && Rscript verify_GAM.R
# Utility, not a pipeline stage. Reads only; writes nothing.
# Exit 0 = all checks pass, 1 = at least one FAIL.
# ================================================================
# WHY THIS EXISTS
#
# Every correction made on 2026-10-18 was the same defect wearing a different
# hat: a number written in PROSE disagreeing with the TABLE it cites. BED tracks
# hours older than the table they summarised, a methods note claiming "every
# statistic" was untouched two commits after that stopped being true, and
# make_report.R hard-coding "38 of 38" in prose so that re-running the report
# could not fix it.
#
# None of these are caught by any stopifnot() in the pipeline, because each
# stage was individually correct. The inconsistency only exists BETWEEN files.
#
# THE RULE THIS ENCODES: results/tables/*.tsv is the single source of truth.
# A statistic is computed here FROM the table, and every file that quotes it
# must contain the current value. Nothing expected is stored in this script or
# anywhere else, so the checker cannot itself go stale - if the analysis
# changes, the facts change with it and the prose is what fails.
#
# Note the polarity: a file passes by CONTAINING the current value, not by
# lacking old ones. A change-note that says "previously 38 of 38, now 27 of 38"
# passes, which is deliberate - that is good documentation, not staleness.
#
# WHAT THIS CANNOT COVER, and it is the important limitation (found 2026-10-19)
#
# The design only works for a statistic that HAS a table. A number computed
# once in a session and written straight into prose is invisible here: there is
# nothing to compute it from, so no fact exists and no check fails. Passing
# 39/39 means every statistic WITH a table agrees with it - not that every
# statistic in the prose is backed.
#
# ================================================================
source("config.R")
suppressPackageStartupMessages(library(data.table))

PASS <- 0L; FAIL <- 0L; SKIP <- 0L
ok   <- function(what, detail = "") { PASS <<- PASS + 1L
  message(sprintf("  PASS  %s%s", what, if (nzchar(detail)) paste0("  [", detail, "]") else "")) }
bad  <- function(what, detail = "") { FAIL <<- FAIL + 1L
  message(sprintf("  FAIL  %s%s", what, if (nzchar(detail)) paste0("\n          ", detail) else "")) }
skip <- function(what, why)         { SKIP <<- SKIP + 1L
  message(sprintf("  skip  %s  [%s]", what, why)) }

rd <- function(p) if (file.exists(p)) fread(p) else NULL

# ---------------------------------------------------------------
# A. FACTS: a statistic computed from a table, and the files that quote it.
#
# `render` turns the current value into the regexes that would legitimately
# appear in prose. A fact passes when EVERY file listed contains at least one.
# ---------------------------------------------------------------
facts <- list()
# `src` is the table the value was computed FROM. Recording it lets family B
# derive the report's dependencies instead of relying on a hand-kept list -
# which is what failed: a GAM_4 re-run rewrote gam4_hotspot_annotation_
# enrichment.tsv AFTER the report was built, and nothing noticed, because the
# report's dependency list named three tables and that was not one of them.
fact <- function(id, value, matcher, quoted_in, note = "", src = NULL) {
  facts[[id]] <<- list(id = id, value = value, matcher = matcher,
                       quoted_in = quoted_in, note = note, src = src)
}

# --- matchers: given a file's text and the current value, did it say it? ------
#
# Two kinds, because counts and measurements fail differently.

# COUNTS are exact, and the surrounding "of M" is what makes them checkable:
# "27" alone appears everywhere, "27 of 38" does not. Tolerates "N/M", an
# intervening markup tag, and the filler that reads naturally in prose ("6 of
# the 10 where both orient"). Deliberately loose about phrasing - a checker
# that fails on wording rather than content gets ignored, which is worse than
# not having one.
m_n_of_m <- function(m) function(txt, v) {
  rx <- sprintf("%d\\s*(of|/|out of)\\s*(<[^>]*>|the|only|just|\\*\\*|\\s)*%d", v, m)
  grepl(rx, txt, perl = TRUE)
}

# MEASUREMENTS cannot be string-matched, because how a number is rendered is a
# writing decision, not a fact. The report prints median H4K20me3 log2FC of
# -0.4975 as "-0.497" - truncated, not rounded - so sprintf("%.3f") produces
# "-0.498" and would fail on prose that is perfectly correct. So: pull every
# numeric token out of the file and ask whether ANY of them is within `tol` of
# the value. The risk is a coincidental match; the defence is that `tol` is
# small relative to the quantity, and that the failure actually being guarded
# against is a number wrong by a lot (38 against 27), not in its last digit.
m_near <- function(tol) function(txt, v) {
  nums <- suppressWarnings(as.numeric(unlist(regmatches(
    txt, gregexpr("-?[0-9]+\\.?[0-9]*", txt)))))
  nums <- nums[is.finite(nums)]
  length(nums) && min(abs(nums - v)) <= tol
}

# --- the ChIP effect sizes the whole study is conditioned on ------------------
L <- rd(file.path(tables_dir, "gam_bins_l2fc_summary.tsv"))
if (!is.null(L)) {
  fact("grid.n_bins", L$n_bins[1], function(txt, v)
         grepl(format(v, big.mark = ","), txt, fixed = TRUE) ||
         grepl(format(v), txt, fixed = TRUE),
       c("make_report.R"), "20 kb autosomal bins", src = "gam_bins_l2fc_summary.tsv")
  for (mk in L$mark) {
    v <- L[mark == mk, median_log2FC]
    fact(sprintf("chip.median_l2fc.%s", mk), v, m_near(0.0015),
         c("make_report.R", "../ABSTRACT.md"),
         sprintf("%s median log2FC on the GAM 20 kb grid", mk),
         src = "gam_bins_l2fc_summary.tsv")
  }
}

# --- the central contact result ----------------------------------------------
# These moved onto the two-level partition on 2026-10-19. "het-het" is now
# LEVEL 1 - Polycomb + H3K9me3 + quiescent - not the old "both bins carry an
# H3K9me3/H4K20me3 peak", and "eu-eu" is TSS + enhancer + transcribed rather
# than "neither bin carries a peak". The values therefore changed, and the
# documents quoting them were updated in the same commit.
CC <- rd(file.path(tables_dir, "gam10_class_change.tsv"))
if (!is.null(CC) && all(c("level", "pair", "band", "pct") %in% names(CC))) {
  L1 <- CC[level == "L1"]
  v <- L1[pair == "heterochromatin-heterochromatin" & band == "<0.2Mb", pct]
  if (length(v) == 1)
    fact("contacts.het_het_short", v, m_near(0.06),
         c("make_report.R", "CLAUDE.md", "../ABSTRACT.md"),
         "level-1 het-het below 0.2 Mb, % change in mean NPMI", src = "gam10_class_change.tsv")
  # Band edges moved on 2026-10-19: the old open-ended ">10Mb" pooled 10-195 Mb,
  # with nearly half its pairs beyond 50 Mb. This fact now names the top band
  # from config so it cannot silently key on a band that no longer exists.
  v2 <- L1[pair == "euchromatin-euchromatin" & band == tail(gam_band_names, 1), pct]
  if (length(v2) == 1)
    # Checked against the RENDERED report, not make_report.R: the class tables
    # are now generated from gam10_class_change.tsv, so the source legitimately
    # contains no literal to match. Same reason as the melting facts below.
    fact("contacts.eu_eu_long", v2, m_near(0.06),
         c("results/GAM_REPORT.html", "../ABSTRACT.md"),
         "level-1 eu-eu in the top distance band, % change in mean NPMI", src = "gam10_class_change.tsv")
  # The H3K9me3 sub-class on its own, which is where the effect concentrates and
  # which is the closest thing to the quantity the old binary measured.
  v3 <- CC[level == "L2" & pair == "H3K9me3-H3K9me3" & band == "<0.2Mb", pct]
  if (length(v3) == 1)
    fact("contacts.k9_k9_short", v3, m_near(0.06),
         c("results/GAM_REPORT.html", "../ABSTRACT.md"),
         "level-2 H3K9me3-H3K9me3 below 0.2 Mb, % change in mean NPMI",
         src = "gam10_class_change.tsv")
  # THE SEPARATION ITSELF, which is what the conclusion now rests on. Below
  # 0.2 Mb every heterochromatin pair rises together, so that band alone cannot
  # distinguish "H3K9me3 self-association" from "heterochromatin in general".
  # Past it they diverge - H3K9me3-H3K9me3 up, H3K9me3-quiescent flat,
  # quiescent-quiescent down - and ABSTRACT.md quotes all three. An earlier
  # conclusion attributed the gain to LAD/LINE-rich quiescent domains, i.e. to
  # the chromatin that is losing; these facts are here so that cannot recur
  # silently.
  for (f in list(
        list(id = "contacts.k9_k9_mid",     pr = "H3K9me3-H3K9me3",   bd = "0.2-1Mb"),
        list(id = "contacts.k9_k9_1_5",     pr = "H3K9me3-H3K9me3",   bd = "1-5Mb"),
        list(id = "contacts.k9_quies_short",pr = "H3K9me3-Quiescent", bd = "<0.2Mb"),
        list(id = "contacts.k9_quies_mid",  pr = "H3K9me3-Quiescent", bd = "0.2-1Mb"),
        list(id = "contacts.k9_quies_1_5",  pr = "H3K9me3-Quiescent", bd = "1-5Mb"),
        list(id = "contacts.quies_short",   pr = "Quiescent-Quiescent", bd = "<0.2Mb"),
        list(id = "contacts.quies_mid",     pr = "Quiescent-Quiescent", bd = "0.2-1Mb"),
        list(id = "contacts.quies_1_5",     pr = "Quiescent-Quiescent", bd = "1-5Mb"))) {
    vv <- CC[level == "L2" & pair == f$pr & band == f$bd, pct]
    if (length(vv) == 1)
      fact(f$id, vv, m_near(0.06), c("../ABSTRACT.md"),
           sprintf("level-2 %s at %s, %% change in mean NPMI", f$pr, f$bd),
           src = "gam10_class_change.tsv")
  }
}

# --- domain SIZE sets the response only in H3K9me3 domains (GAM_5 section D) --
# The claim the report and the abstract make is a CONTRAST: inside H3K9me3-
# dominant domains the short-range gain grows with domain size, and inside
# Quiescent-dominant ones it does not. Both arms are pinned, because quoting
# only the H3K9me3 arm would let the comparison it is half of go stale.
local({
  f <- file.path(tables_dir, "gam13_domain_size_change_ci.tsv")
  if (!file.exists(f)) return(invisible())
  DS <- fread(f)
  for (g in list(
      list(id = "domsize.k9_small",  dc = "H3K9me3-dominant",   sz = "<80kb",    bd = "<0.2Mb"),
      list(id = "domsize.k9_mid",    dc = "H3K9me3-dominant",   sz = "80-320kb", bd = "<0.2Mb"),
      list(id = "domsize.k9_large",  dc = "H3K9me3-dominant",   sz = ">320kb",   bd = "<0.2Mb"),
      list(id = "domsize.quies_small", dc = "Quiescent-dominant", sz = "<80kb",    bd = "<0.2Mb"),
      list(id = "domsize.quies_large", dc = "Quiescent-dominant", sz = ">320kb",   bd = "<0.2Mb"),
      list(id = "domsize.k9_large_mid_band", dc = "H3K9me3-dominant", sz = ">320kb", bd = "0.2-1Mb"))) {
    vv <- DS[kind == "intra" & dclass == g$dc & size == g$sz & band == g$bd, pct]
    if (length(vv) == 1)
      fact(g$id, vv, m_near(0.06), c("make_report.R"),
           sprintf("intra-domain %s %s at %s, %% change in mean NPMI",
                   g$dc, g$sz, g$bd),
           src = "gam13_domain_size_change_ci.tsv")
  }
})

# --- the per-chromosome contrast behind "17 of 19 autosomes" -----------------
# This was quoted in the plan stock-take and in ABSTRACT.md but computed nowhere
# that persists, so nothing could catch it drifting. The claim is that
# heterochromatin-heterochromatin contacts rise RELATIVE to euchromatin-
# euchromatin below 1 Mb on most autosomes. Recomputed here as the sign of
# (het-het - eu-eu) per chromosome, which is what makes it offset-invariant:
# the per-chromosome level of pct is dominated by a global offset and the
# contrast cancels it. Reads level 1 since 2026-10-19.
PB <- rd(file.path(tables_dir, "gam15_perchrom_class_band.tsv"))
if (!is.null(PB) && all(c("chr", "level", "pair", "band", "pct") %in% names(PB))) {
  w <- dcast(PB[level == "L1" & band %in% c("<0.2Mb", "0.2-1Mb")],
             chr + band ~ pair, value.var = "pct")
  if (all(c("euchromatin-euchromatin", "heterochromatin-heterochromatin") %in% names(w))) {
    w[, contrast := `heterochromatin-heterochromatin` - `euchromatin-euchromatin`]
    # Require it in BOTH sub-bands, not in their mean: a mean could be carried
    # by one band while the other disagrees, and the claim is about both.
    n_pos <- w[, .(all_pos = all(contrast > 0, na.rm = TRUE)), by = chr][, sum(all_pos)]
    fact("contacts.contrast_autosomes", n_pos, m_n_of_m(length(gam_autosomes)),
         c("../ABSTRACT.md"),
         "autosomes where het-het beats eu-eu in BOTH sub-bands below 1 Mb",
         src = "gam15_perchrom_class_band.tsv")
  }
}

# --- the large one-directional blocks and what is inside them -----------------
# results/LARGE_BLOCKS.md is a hand-written document full of numbers, which is
# exactly the shape of file this checker exists for.
local({
  B <- rd(file.path(tables_dir, "gam6_large_blocks.tsv"))
  if (!is.null(B) && nrow(B)) {
    cnt <- function(txt, v) grepl(sprintf("\\b%d\\b", v), txt, perl = TRUE)
    fact("blocks.n_lost", sum(B$dirn == "lost"), cnt,
         c("results/LARGE_BLOCKS.md"), "large blocks that LOSE contacts",
         src = "gam6_large_blocks.tsv")
    fact("blocks.n_gained", sum(B$dirn == "gained"), cnt,
         c("results/LARGE_BLOCKS.md"), "large blocks that GAIN contacts",
         src = "gam6_large_blocks.tsv")
  }
  G <- rd(file.path(tables_dir, "gam8_block_contact_budget.tsv"))
  if (!is.null(G) && nrow(G)) {
    v <- G[pair == "lost-lost" & band == "5-10Mb", pct]
    if (length(v) == 1)
      fact("blocks.lost_lost_5_10Mb", v, m_near(0.06), c("results/LARGE_BLOCKS.md"),
           "both ends in a lost block, 5-10 Mb, % change in mean NPMI",
           src = "gam8_block_contact_budget.tsv")
    v <- G[pair == "gained-gained" & band == "1-5Mb", pct]
    if (length(v) == 1)
      fact("blocks.gained_gained_1_5Mb", v, m_near(0.06), c("results/LARGE_BLOCKS.md"),
           "both ends in a gained block, 1-5 Mb, % change in mean NPMI",
           src = "gam8_block_contact_budget.tsv")
  }
  K <- rd(file.path(tables_dir, "gam16_block_class_contrast_ci.tsv"))
  if (!is.null(K) && nrow(K)) {
    v <- K[dirn == "lost" & band == "0.2-1Mb", contrast]
    if (length(v) == 1)
      fact("blocks.k9_vs_quies_lost", v, m_near(0.06), c("results/LARGE_BLOCKS.md"),
           "H3K9me3-H3K9me3 minus quiescent-quiescent INSIDE lost blocks, 0.2-1 Mb",
           src = "gam16_block_class_contrast_ci.tsv")
  }
})

# --- hotspot domains ----------------------------------------------------------
for (d in list(list(f = "hotspot_domains_lost_adj.tsv",   id = "hotspot.n_lost",
                    nm = "merged hotspot domains that LOSE contacts"),
               list(f = "hotspot_domains_gained_adj.tsv", id = "hotspot.n_gained",
                    nm = "merged hotspot domains that GAIN contacts"))) {
  H <- rd(file.path(tables_dir, d$f))
  if (!is.null(H))
    fact(d$id, nrow(H), function(txt, v)
           grepl(sprintf("\\b%d\\b", v), txt, perl = TRUE),
         c("make_report.R"), d$nm, src = d$f)
}

E <- rd(file.path(tables_dir, "gam4_hotspot_annotation_enrichment.tsv"))
if (!is.null(E) && "fold_rotated" %in% names(E)) {
  for (a in list(list(key = "constitutive LAD",                  id = "hotspot.fold_lad"),
                 list(key = "repeat: LINE-high (top quartile)",  id = "hotspot.fold_line"))) {
    v <- E[domain_set == "lost_adj" & annotation == a$key, fold_rotated]
    if (length(v) == 1)
      fact(a$id, v, m_near(0.006), c("make_report.R", "../ABSTRACT.md"),
           sprintf("lost domains vs rotation null: %s", a$key),
           src = "gam4_hotspot_annotation_enrichment.tsv")
  }
}
# THE THREE-NULL COMPARISON, SPLIT BY DIRECTION OF EFFECT. The report used to
# read the hypergeometric's 18-of-104 as it being conservative. It is computed
# one-sided on the upper tail, so it cannot call a depletion at all, and half
# the table is depletions - on the half it can see it calls MORE than the
# rotation, which is what the _INFLATED in the column name means. Both halves
# are checked so neither the headline nor the split can go stale.
if (!is.null(E) && all(c("fold_rotated", "p_hypergeom_INFLATED",
                         "p_shuffled", "p_rotated") %in% names(E))) local({
  cnt <- function(txt, v) grepl(sprintf("\\b%d\\b", v), txt, perl = TRUE)
  for (half in list(list(sel = quote(fold_rotated >  1), lab = "enriched"),
                    list(sel = quote(fold_rotated <= 1), lab = "depleted"))) {
    H <- E[eval(half$sel)]
    for (nl in list(list(col = "p_hypergeom_INFLATED", id = "hypergeom"),
                    list(col = "p_shuffled",           id = "shuffle"),
                    list(col = "p_rotated",            id = "rotation"))) {
      fact(sprintf("hotspot.nulls_%s_%s", nl$id, half$lab),
           sum(H[[nl$col]] < 0.05, na.rm = TRUE), cnt, c("make_report.R"),
           sprintf("%s cells significant at 0.05 under the %s null", half$lab, nl$id),
           src = "gam4_hotspot_annotation_enrichment.tsv")
    }
  }
})

# --- coupling -----------------------------------------------------------------
CS <- rd(file.path(tables_dir, "coupling_spearman.tsv"))
if (!is.null(CS)) {
  v <- CS[mark == "H4K20me3" & variant == "raw" & direction == "lost", rho]
  if (length(v) == 1)
    fact("coupling.h4k20me3_raw_lost", v, m_near(0.0006),
         c("make_report.R"), "the headline coupling rho", src = "coupling_spearman.tsv")
  n <- CS[mark == "H4K20me3" & variant == "raw" & direction == "lost", n]
  if (length(n) == 1)
    fact("coupling.n", n, function(txt, v)
           grepl(format(v, big.mark = ","), txt, fixed = TRUE),
         c("make_report.R"), "bins behind the coupling rho", src = "coupling_spearman.tsv")
  # The NET variant, added 2026-10-19. The melting section compares its own rho
  # against "the mark-vs-contact coupling", and the honest comparator there is
  # the net measure, because the melting rho is also against net. The lost-only
  # rho above is a different quantity and quoting it in that sentence invited
  # the reader to compare a signed number with one of the opposite sign.
  vn <- CS[mark == "H4K20me3" & variant == "raw" & direction == "net", rho]
  if (length(vn) == 1)
    fact("coupling.h4k20me3_raw_net", vn, m_near(0.0006),
         c("make_report.R"), "H4K20me3 vs NET contact change",
         src = "coupling_spearman.tsv")
}

# --- melting vs covariates ----------------------------------------------------
# These two rho values are drawn into the figure by GAM_2c and quoted in the
# report caption. Until 2026-10-19 the caption carried them as literals with the
# sign they had BEFORE `net` was flipped KO-oriented on 2026-10-18, so the report
# said -0.305 / +0.071 for values that are +0.305 / -0.071. The caption now reads
# the table, which makes that drift impossible; this check covers the table's
# continued existence and the sign, so a future flip cannot pass unnoticed.
MR <- rd(file.path(tables_dir, "gam2c_melting_covariate_rho.tsv"))
if (!is.null(MR) && all(c("w", "r") %in% names(MR))) {
  vnet <- MR[grepl("net differential contacts", w, fixed = TRUE), r]
  if (length(vnet) == 1)
    fact("melting.rho_net_contacts", vnet, m_near(0.002),
         c("results/GAM_REPORT.html"),
         "melting score vs net differential contacts",
         src = "gam2c_melting_covariate_rho.tsv")
  vk20 <- MR[grepl("H4K20me3", w, fixed = TRUE), r]
  if (length(vk20) == 1)
    fact("melting.rho_h4k20me3", vk20, m_near(0.002),
         c("results/GAM_REPORT.html"),
         "melting score vs H4K20me3 log2FC",
         src = "gam2c_melting_covariate_rho.tsv")
}

message("\n== A. quoted statistics match the tables they cite ==")
if (!length(facts)) skip("all fact checks", "no tables on disk")
for (f in facts) {
  miss <- character()
  for (fp in f$quoted_in) {
    if (!file.exists(fp)) { miss <- c(miss, sprintf("%s (absent)", fp)); next }
    txt <- paste(readLines(fp, warn = FALSE), collapse = "\n")
    # make_report.R writes minus signs as the HTML entity, so a raw read sees
    # "&minus;0.497" and a numeric scan pulls out +0.497 - the right magnitude
    # with the wrong sign, which is exactly the kind of near-miss that looks
    # like a finding. Normalise the MINUS entity and U+2212 only. Deliberately
    # NOT &ndash;, which is a range dash: turning "0.2&ndash;1 Mb" into
    # "0.2-1 Mb" would manufacture a -1 that was never written.
    txt <- gsub("&minus;|\u2212", "-", txt)
    if (!isTRUE(f$matcher(txt, f$value))) miss <- c(miss, fp)
  }
  lab <- sprintf("%s = %s", f$id,
                 if (f$value %% 1 == 0) format(f$value) else sprintf("%.3f", f$value))
  if (length(miss)) {
    bad(lab, sprintf("%s\n          not found in: %s", f$note, paste(miss, collapse = ", ")))
  } else {
    ok(lab, f$note)
  }
}

# ---------------------------------------------------------------
# B. FRESHNESS: a derived file must be newer than everything it is built from.
# This is what made the BEDs wrong - they were correct for an input that had
# since been replaced, so nothing in them was detectably false on its own.
#
# KNOWN FALSE POSITIVE, and why it is left in. This compares mtimes, so a stage
# re-run that produces IDENTICAL content still marks every consumer stale -
# GAM_2 re-run from cache to record its method version rewrote
# npmi_bin_features.rds with the same features, and that alone flagged GAM_3,
# GAM_4 and the report. Content hashing would avoid it, but it would need a
# record of what each output was built FROM, which nothing currently stores.
#
# Left as a FAIL rather than a warning because the resolution - re-run the
# consumer - is cheap, is the honest way to make the claim true, and exercises
# the repo rule that results/ must stay fully regenerable from committed code.
# A warning here would be ignored, which is how the BEDs went four hours stale.
# ---------------------------------------------------------------
fact_srcs <- function()
  c(unique(file.path(tables_dir,
      Filter(Negate(is.null), lapply(facts, `[[`, "src")))), "make_report.R")

message("\n== B. derived files are newer than their inputs ==")
mt <- function(p) if (file.exists(p)) file.mtime(p) else NA
deps <- list(
  # Derived, not hand-listed: the report must post-date every table any fact is
  # computed from, plus its own source.
  list(out = file.path(results_dir, "GAM_REPORT.html"),       ins = fact_srcs()),
  list(out = file.path(results_dir, "GAM_REPORT_artifact.html"), ins = fact_srcs()),
  # NOTE the name is misleading: gam6_block_characterisation.tsv is written by
  # GAM_4_Large_Scale.R, not GAM_6.
  list(out = file.path(tables_dir, "gam6_block_characterisation.tsv"),
       ins = c(file.path(rds_dir, "gam_bins_annotated.rds"),
               file.path(rds_dir, "npmi_bin_features.rds"))),
  list(out = file.path(tables_dir, "gam4_size_perm_overlap.tsv"),
       ins = c(file.path(rds_dir, "gam_bins_annotated.rds"),
               file.path(rds_dir, "npmi_bin_features.rds"))),
  list(out = file.path(tables_dir, "gam5_decompaction_verdict.tsv"),
       ins = file.path(rds_dir, "npmi_bin_features.rds")),
  list(out = file.path(tables_dir, "coupling_spearman.tsv"),
       ins = file.path(rds_dir, "npmi_bin_features.rds")))
for (d in deps) {
  to <- mt(d$out)
  if (is.na(to)) { skip(basename(d$out), "not built"); next }
  older <- d$ins[!is.na(vapply(d$ins, function(i) as.numeric(mt(i)), 0)) &
                 vapply(d$ins, function(i) isTRUE(mt(i) > to), TRUE)]
  if (length(older)) {
    bad(sprintf("%s is stale", basename(d$out)),
        sprintf("built %s, but these inputs are newer: %s",
                format(to, "%m-%d %H:%M"),
                paste(sprintf("%s (%s)", basename(older),
                              vapply(older, function(i) format(mt(i), "%m-%d %H:%M"), "")),
                      collapse = ", ")))
  } else {
    ok(sprintf("%s", basename(d$out)), format(to, "%m-%d %H:%M"))
  }
}

# ---------------------------------------------------------------
# C. INVARIANTS: things that must be true across files, not within one.
# ---------------------------------------------------------------
message("\n== C. cross-file invariants ==")

# C2. Cache versions in the stored rds match the script that would write them.
ver_in_script <- function(path, pat) {
  if (!file.exists(path)) return(NA_integer_)
  ln <- readLines(path, warn = FALSE)
  hit <- regmatches(ln, regexpr(pat, ln, perl = TRUE))
  hit <- unlist(hit)
  if (!length(hit)) return(NA_integer_)
  as.integer(gsub("\\D", "", hit[1]))
}
for (v in list(
  list(rds = file.path(rds_dir, "npmi_bin_features.rds"),
       src = "GAM_2_NPMI_Features.R", pat = "version\\s*=\\s*\\d+L", nm = "GAM_2"))) {
  if (!file.exists(v$rds)) { skip(sprintf("%s cache version", v$nm), "rds absent"); next }
  # GAM_2 stores params under meta.
  stored <- tryCatch({
    x <- readRDS(v$rds)
    if (!is.null(x$params$version)) x$params$version else x$meta$params$version
  }, error = function(e) NULL)
  want   <- ver_in_script(v$src, v$pat)
  if (is.null(stored) || is.na(want)) {
    skip(sprintf("%s method version", v$nm),
         "not recorded on this rds - populates on the next run of that stage")
    next
  }
  if (!identical(as.integer(stored), as.integer(want))) {
    bad(sprintf("%s cache version", v$nm),
        sprintf("rds holds v%s but %s declares v%s - the cache was reused across a method change",
                stored, v$src, want))
  } else {
    ok(sprintf("%s cache version", v$nm), sprintf("v%s", stored))
  }
}

# C3. Sign convention. Red = gained = positive = cKO-oriented, everywhere.
# GAM_4 and GAM_2c both violated this once, in opposite directions, and GAM_2c
# contradicted itself inside one file.
sign_files <- c("GAM_4_Large_Scale.R", "GAM_2c_Meltron.R", "GAM_7_Figures.R")
wrong <- character()
for (fp in sign_files) {
  if (!file.exists(fp)) next
  txt <- paste(readLines(fp, warn = FALSE), collapse = "\n")
  # the lost/blue and gained/red hex pair, as used repo-wide
  if (grepl("gained\\s*=\\s*\"#2166AC\"", txt) || grepl("lost\\s*=\\s*\"#B2182B\"", txt))
    wrong <- c(wrong, fp)
}
if (length(wrong)) {
  bad("gained=red / lost=blue in fill scales",
      paste("inverted in:", paste(wrong, collapse = ", ")))
} else {
  ok("gained=red / lost=blue in fill scales", paste(length(sign_files), "files"))
}

# C3b. The quantity, not just the palette. The sign bug that actually occurred
# was `net` computed WT-oriented while the legend said cKO - the colours were
# right and the number was backwards, which no palette check can see. Every
# definition of `net` must put GAINED in the numerator, so positive = gained =
# more in the cKO = red, consistently with everything else.
net_files <- c("GAM_3_Coupling.R", "GAM_4_Large_Scale.R")
netbad <- character()
for (fp in net_files) {
  if (!file.exists(fp)) next
  ln <- grep("net\\s*:=", readLines(fp, warn = FALSE), value = TRUE)
  ln <- grep("log2", ln, value = TRUE)
  if (!length(ln)) { netbad <- c(netbad, sprintf("%s (no net definition found)", fp)); next }
  # gained must appear before lost inside the ratio
  okf <- vapply(ln, function(l) {
    g <- regexpr("gained", l, fixed = TRUE); o <- regexpr("lost", l, fixed = TRUE)
    g > 0 && o > 0 && g < o }, TRUE)
  if (!all(okf)) netbad <- c(netbad, fp)
}
if (length(netbad)) {
  bad("net = log2((gained+1)/(lost+1)), i.e. positive = gained = cKO",
      paste("wrong or missing in:", paste(netbad, collapse = ", ")))
} else {
  ok("net = log2((gained+1)/(lost+1)), i.e. positive = gained = cKO",
     paste(length(net_files), "files"))
}

# C3c. Table prefix -> the stage that actually writes it.
#
# The gamN_ prefix on a table is a LEGACY SECTION number from a pre-consolidation
# numbering, not the current stage: GAM_4 writes gam4_* AND gam6_*, GAM_5 writes
# gam10_/gam11_/gam13_*, GAM_6 writes gam9_/gam14_*, GAM_7 writes gam12_/gam15_*.
# The stage files still carry the old labels in their completion messages for
# the same reason. This is a live trap - it made the first version of the
# freshness check above list gam6_block_characterisation.tsv as a GAM_6 output,
# and that produced a
# false alarm. Pinning the mapping here means the next person gets a failure
# instead of a wrong assumption.
writer_map <- list(
  "gam6_block_characterisation.tsv" = "GAM_4_Large_Scale.R",
  "gam4_hotspot_annotation_enrichment.tsv" = "GAM_4_Large_Scale.R",
  "gam10_class_change.tsv"          = "GAM_5_Heterochromatin.R",
  "gam12_block_ranking.tsv"         = "GAM_7_Figures.R",
  "gam12_k9_large_domains.tsv"      = "GAM_7_Figures.R",
  "gam13_domain_size_change_ci.tsv" = "GAM_5_Heterochromatin.R",
  "gam11_sig_vs_h4k20_rho.tsv"      = "GAM_5_Heterochromatin.R",
  "gam14_pcdh_tad_summary.tsv"      = "GAM_6_Targets.R",
  "coupling_spearman.tsv"           = "GAM_3_Coupling.R",
  "gam2c_melting_covariate_rho.tsv" = "GAM_2c_Meltron.R")
wrong_writer <- character()
for (tb in names(writer_map)) {
  src <- writer_map[[tb]]
  if (!file.exists(src)) next
  if (!any(grepl(tb, readLines(src, warn = FALSE), fixed = TRUE)))
    wrong_writer <- c(wrong_writer, sprintf("%s not written by %s", tb, src))
}
if (length(wrong_writer)) {
  bad("table prefixes map to the stage that writes them",
      paste(wrong_writer, collapse = "; "))
} else {
  ok("table prefixes map to the stage that writes them",
     sprintf("%d tables (prefix is a LEGACY section number, not the stage)",
             length(writer_map)))
}

message(sprintf("\n%s\n  %d passed, %d failed, %d skipped\n",
                strrep("=", 60), PASS, FAIL, SKIP))
if (FAIL > 0L) {
  message("Prose disagrees with the tables, or a derived file is stale.\n",
          "The TABLE is right by definition - fix the file that quotes it.\n")
  quit(status = 1L)
}
