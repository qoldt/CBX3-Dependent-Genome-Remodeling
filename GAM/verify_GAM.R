# ================================================================
# verify_GAM.R - is the analysis internally consistent?
#   cd GAM && Rscript verify_GAM.R
# Utility, not a pipeline stage. Reads only; writes nothing.
# Exit 0 = all checks pass, 1 = at least one FAIL.
# ================================================================
# WHY THIS EXISTS
#
# Every correction made on 2026-10-18 was the same defect wearing a different
# hat: a number written in PROSE disagreeing with the TABLE it cites.
#
#   - compartments_{WT,HP1gKO,diff}.bed were four hours older than the
#     orientation table they summarise, and their track description quoted
#     pre-fix statistics;
#   - export_compartment_beds.R's header described an anti-correlated axis as
#     something that ships, which GAM_2 now forbids with an assertion;
#   - BIN_BALANCING.md section 5 said "every statistic in results/tables/" was
#     untouched by balancing, which stopped being true two commits earlier;
#   - CLAUDE.md carried pre-v8 quality numbers and "currently version 2";
#   - make_report.R hard-coded "38 of 38" and "PC1 on only chr7" in prose, so
#     re-running the report could not fix them, and the published page
#     contradicted a table cited in its own paragraph.
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
fact <- function(id, value, render, quoted_in, note = "") {
  facts[[id]] <<- list(id = id, value = value, render = render,
                       quoted_in = quoted_in, note = note)
}

# "N of M", tolerating "N/M", an intervening markup tag or newline, and the
# short filler words that read naturally in prose ("6 of the 10 where both
# orient"). Kept deliberately loose: a checker that fails on phrasing rather
# than on content gets ignored, which is worse than not having it.
n_of_m <- function(n, m)
  sprintf("%d\\s*(of|/|out of)\\s*(<[^>]*>|the|only|just|\\*\\*|\\s)*%d", n, m)
# a number to 2 and 3 decimal places, e.g. 0.356 also matches "0.36"
dp <- function(x) sprintf("(%.3f|%.2f)", x, x)

O <- rd(file.path(tables_dir, "npmi_compartment_orientation.tsv"))
if (!is.null(O)) {
  Oa <- O[chr %in% gam_autosomes]
  ori <- Oa[oriented == TRUE]
  # n_all is AUTOSOMES x genotypes (38). nrow(O) would be 42 - it carries chrX
  # and chrY, which never orient and which no document counts.
  n_ori  <- nrow(ori); n_all <- nrow(Oa)
  n_pc1  <- sum(ori$which_pc == 1, na.rm = TRUE)
  both   <- ori[, .N, by = chr][N == 2, chr]
  same   <- sum(vapply(both, function(c) {
                  w <- ori[chr == c, which_pc]; length(unique(w)) == 1 }, TRUE))
  med    <- median(abs(ori$rho_orient), na.rm = TRUE)

  fact("compartment.oriented", n_ori, function(v) n_of_m(v, n_all),
       c("CLAUDE.md", "make_report.R", "results/COMPARTMENTS.md"),
       "oriented chromosome x genotype combinations")
  fact("compartment.pc1", n_pc1, function(v) n_of_m(v, n_ori),
       c("CLAUDE.md", "make_report.R", "results/COMPARTMENTS.md",
         "results/BIN_BALANCING.md"),
       "combinations whose compartment is on PC1")
  fact("compartment.both_genotypes", length(both), function(v) n_of_m(v, length(gam_autosomes)),
       c("CLAUDE.md", "make_report.R", "export_compartment_beds.R"),
       "autosomes where BOTH genotypes orient")
  fact("compartment.same_pc", same, function(v) n_of_m(v, length(both)),
       c("CLAUDE.md", "make_report.R", "export_compartment_beds.R"),
       "of those, where both pick the same eigenvector")
  fact("compartment.median_rho", med, dp,
       c("CLAUDE.md", "make_report.R"),
       "median |rho| of the per-genotype axes vs active chromatin")
}

S <- rd(file.path(tables_dir, "gam2b_shared_axis_summary.tsv"))
if (!is.null(S)) {
  det <- S[determined == TRUE]
  fact("shared.determined", nrow(det), function(v) n_of_m(v, nrow(S)),
       c("CLAUDE.md", "make_report.R", "results/COMPARTMENTS.md"),
       "autosomes with a shared axis")
  fact("shared.pc1", sum(det$which_pc == 1, na.rm = TRUE), function(v) n_of_m(v, nrow(det)),
       c("CLAUDE.md", "make_report.R", "results/COMPARTMENTS.md",
         "results/BIN_BALANCING.md"),
       "shared axes carried by PC1")
  fact("shared.median_rho", median(abs(det$rho_ref), na.rm = TRUE), dp,
       c("CLAUDE.md", "make_report.R", "results/COMPARTMENTS.md"),
       "median |rho| of the shared axis vs active chromatin")
}

message("\n== A. quoted statistics match the tables they cite ==")
if (!length(facts)) skip("all fact checks", "no compartment tables on disk")
for (f in facts) {
  rx <- f$render(f$value)
  miss <- character()
  for (fp in f$quoted_in) {
    if (!file.exists(fp)) { miss <- c(miss, sprintf("%s (absent)", fp)); next }
    txt <- paste(readLines(fp, warn = FALSE), collapse = "\n")
    if (!any(vapply(rx, function(r) grepl(r, txt, perl = TRUE), TRUE)))
      miss <- c(miss, fp)
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
# ---------------------------------------------------------------
message("\n== B. derived files are newer than their inputs ==")
mt <- function(p) if (file.exists(p)) file.mtime(p) else NA
deps <- list(
  list(out = file.path(bed_dir, "compartments_WT.bed"),
       ins = file.path(rds_dir, "shared_axis_scores.rds")),
  list(out = file.path(bed_dir, "compartments_HP1gKO.bed"),
       ins = file.path(rds_dir, "shared_axis_scores.rds")),
  list(out = file.path(bed_dir, "compartments_diff.bed"),
       ins = file.path(rds_dir, "shared_axis_scores.rds")),
  list(out = file.path(tables_dir, "gam2b_shared_axis_summary.tsv"),
       ins = file.path(rds_dir, "gam_bins_annotated.rds")),
  list(out = file.path(results_dir, "GAM_REPORT.html"),
       ins = c(file.path(tables_dir, "npmi_compartment_orientation.tsv"),
               file.path(tables_dir, "gam2b_shared_axis_summary.tsv"),
               "make_report.R")))
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

# C1. Every oriented axis agrees with its orientation check. GAM_2 asserts this
# before writing; re-checked here because a hand-edited table would slip past it.
if (!is.null(O) && "rho_check" %in% names(O)) {
  b <- O[oriented == TRUE & is.finite(rho_check) & rho_check < 0]
  if (nrow(b)) {
    bad("oriented axes all have rho_check > 0",
        sprintf("%d row(s) anti-correlate: %s", nrow(b),
                paste(b$chr, b$genotype, collapse = ", ")))
  } else {
    ok("oriented axes all have rho_check > 0",
       sprintf("%d rows", sum(O$oriented, na.rm = TRUE)))
  }
} else skip("rho_check invariant", "column absent")

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
  list(rds = file.path(rds_dir, "shared_axis_scores.rds"),
       src = "GAM_2b_Shared_Axis.R", pat = "version\\s*=\\s*\\d+L", nm = "GAM_2b"),
  list(rds = file.path(rds_dir, "npmi_bin_features.rds"),
       src = "GAM_2_NPMI_Compartments.R", pat = "version\\s*=\\s*\\d+L", nm = "GAM_2"))) {
  if (!file.exists(v$rds)) { skip(sprintf("%s cache version", v$nm), "rds absent"); next }
  # GAM_2b stores params at the top level; GAM_2 stores them under meta.
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

# C4. The BED tracks must not silently cover fewer chromosomes than they claim.
bedf <- file.path(bed_dir, "compartments_WT.bed")
if (file.exists(bedf) && !is.null(S)) {
  bchr <- unique(fread(bedf, skip = 1, header = FALSE)$V1)
  dchr <- S[determined == TRUE, chr]
  gap  <- setdiff(dchr, bchr)
  if (length(gap)) {
    bad("BED covers every chromosome with a shared axis",
        sprintf("missing: %s", paste(gap, collapse = ", ")))
  } else {
    ok("BED covers every chromosome with a shared axis",
       sprintf("%d chromosomes", length(bchr)))
  }
} else {
  skip("BED chromosome coverage", "BED or summary absent")
}

message(sprintf("\n%s\n  %d passed, %d failed, %d skipped\n",
                strrep("=", 60), PASS, FAIL, SKIP))
if (FAIL > 0L) {
  message("Prose disagrees with the tables, or a derived file is stale.\n",
          "The TABLE is right by definition - fix the file that quotes it.\n")
  quit(status = 1L)
}
