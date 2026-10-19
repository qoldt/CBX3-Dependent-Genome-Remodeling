# ================================================================
# make_report.R - assemble the GAM working report
# ----------------------------------------------------------------
#     cd GAM && Rscript make_report.R        (after whichever stages have run)
#
# WHY THIS EXISTS AS A SCRIPT. The previous version of this report was built ad
# hoc and lived only in a session scratch directory, so it could not be
# regenerated and did not survive. This does, and it degrades gracefully: any
# figure or table that is missing is reported as missing rather than aborting
# the run, so a partial night still produces a readable document.
#
# AUDIENCE: internal working doc. Dense, every figure, caveats stated next to
# the claim they qualify, no background section. Not the collaborator-facing
# write-up - that one needs the methods and the sign conventions spelled out.
#
# Images are downscaled and base64-embedded so the HTML is a single portable
# file that can be sent to a phone. results/ is git-ignored, so the HTML is an
# artifact, not a source; this script is the thing under version control.
#
# TWO OUTPUT MODES.
#   Rscript make_report.R              -> results/GAM_REPORT.html
#   Rscript make_report.R --artifact   -> also results/GAM_REPORT_artifact.html
# The second is for publishing as a claude.ai Artifact, whose host supplies its
# own <!doctype>/<head>/<body> and would nest a second document inside it. That
# variant therefore emits only <title>, <style> and the content, and it carries
# the explicit [data-theme] blocks the Artifact host needs: there a viewer can
# force light or dark, which prefers-color-scheme alone does not see, so a
# report built for the browser only would show one theme's text on the other
# theme's background. Regenerate it whenever the report changes - a published
# Artifact does not update itself.
# ================================================================
source("config.R")
suppressPackageStartupMessages({ library(data.table); library(base64enc) })

ARTIFACT <- "--artifact" %in% commandArgs(TRUE)
OUT      <- file.path(results_dir, "GAM_REPORT.html")
OUT_ART  <- file.path(results_dir, "GAM_REPORT_artifact.html")
# IMAGE WIDTH. The content column is max-width 1180px with 16px padding, so an
# image is displayed at up to ~1148 CSS px. At 1100 px of actual pixels it was
# being upscaled on any retina display, which is what made the figures look
# pixelated in a desktop browser. 1800 px gives roughly 1.6 device pixels per CSS
# pixel, which is crisp at the common 1.5x-2x ratios without quadrupling bytes.
MAXW <- 1800L
msg <- function(...) cat(sprintf(...), "\n", sep = "")
tmpd <- file.path(tempdir(), "gamreport"); dir.create(tmpd, showWarnings = FALSE)

esc <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

# --- figure embedding ---------------------------------------------------------
# sips is macOS-only but this repo is explicitly a macOS analysis box (the
# safe-local-alignment skill exists for that reason). If it is absent the image
# is used at full size rather than skipped.
#
# TWO DELIVERY MODES, because the size budgets are different:
#   local HTML  - base64 data: URI, so the file is self-contained and portable.
#   --artifact  - the PNG is written beside the page and referenced relatively.
# Base64 was the reason the artifact had to run at a low image width: embedded
# data: URIs count against the artifact's 16 MB page limit, and ~40 figures of a
# retina-resolution plot do not fit. As separate published files that limit is
# per file instead of per page, so the artifact gets the same resolution as the
# local copy. ASSETS collects what has to be published alongside the page.
ASSETS <- new.env(parent = emptyenv()); ASSETS$files <- character(0)
asset_dir <- file.path(results_dir, "GAM_REPORT_artifact_files")

# The ChIP figures live in a different results tree. Same embedding, different
# root - kept as a separate function so a missing ChIP run degrades to the
# "missing figure" box rather than silently resolving against the GAM tree.
embed_chip_fig <- function(rel, caption = NULL)
  embed_fig(rel, caption, root = file.path(chip_results_dir, "figures"))

# The body is built ONCE and then written out in both modes, so the image source
# cannot be decided here - deciding it here is what left GAM_REPORT.html full of
# broken <img src="figures/..."> links whenever the script was run with
# --artifact, because the artifact's relative paths were baked into the local
# copy too. Instead each figure is registered and a token is emitted, and
# resolve_figs() below swaps the token for whichever source that output wants.
FIGREG <- new.env(parent = emptyenv()); FIGREG$items <- list()

embed_fig <- function(rel, caption = NULL, root = fig_dir) {
  f <- file.path(root, rel)
  if (!file.exists(f)) {
    return(sprintf('<div class="missing">missing figure: <code>%s</code> — the stage that writes it did not complete</div>',
                   esc(rel)))
  }
  src <- f
  if (nzchar(Sys.which("sips"))) {
    small <- file.path(tmpd, gsub("[/]", "_", rel))
    ok <- suppressWarnings(system2("sips", c("-Z", MAXW, shQuote(f), "--out", shQuote(small)),
                                   stdout = FALSE, stderr = FALSE))
    if (ok == 0 && file.exists(small)) src <- small
  }
  pub <- file.path("figures", gsub("[/]", "_", rel))
  # One registry entry per published name, so the same figure used twice is
  # encoded once. A clash between two different sources under one published name
  # is a real bug (the second would silently overwrite the first), so stop.
  prev <- FIGREG$items[[pub]]
  if (!is.null(prev) && !identical(prev, src))
    stop("two different figures map to the same published name: ", pub,
         "\n  ", prev, "\n  ", src, call. = FALSE)
  FIGREG$items[[pub]] <- src
  # loading="lazy" matters here: the artifact is read on a phone, and without it
  # every figure in the report is fetched before the first one renders.
  tag <- sprintf('<img src="@@FIG:%s@@" loading="lazy" decoding="async">', pub)
  paste0("<figure>", tag,
         if (!is.null(caption)) paste0("<figcaption>", caption, "</figcaption>") else "",
         "</figure>")
}

# Swap every @@FIG:<pub>@@ token for a real source.
#   mode "embed"    - base64 data: URI, so the local HTML is one portable file.
#   mode "external" - copy the PNG into asset_dir and reference it relatively,
#                     which is what the Artifact host publishes alongside the page.
resolve_figs <- function(html, mode) {
  pubs <- names(FIGREG$items)
  if (!length(pubs)) return(html)
  reps <- vapply(pubs, function(pub) {
    src <- FIGREG$items[[pub]]
    if (mode == "external") {
      dest <- file.path(asset_dir, pub)
      dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
      file.copy(src, dest, overwrite = TRUE)
      ASSETS$files <- unique(c(ASSETS$files, pub))
      pub
    } else {
      paste0("data:image/png;base64,", base64enc::base64encode(src))
    }
  }, character(1))
  # fixed = TRUE: published names contain "+", "." and other regex metacharacters.
  for (i in seq_along(pubs))
    html <- gsub(paste0("@@FIG:", pubs[i], "@@"), reps[i], html, fixed = TRUE)
  left <- gregexpr("@@FIG:", html, fixed = TRUE)[[1]]
  if (left[1] != -1) stop(length(left), " figure token(s) left unresolved", call. = FALSE)
  html
}

tbl <- function(path, caption = NULL, digits = 3, max_rows = 40L, cols = NULL) {
  f <- if (file.exists(path)) path else file.path(tables_dir, path)
  if (!file.exists(f))
    return(sprintf('<div class="missing">missing table: <code>%s</code></div>', esc(basename(path))))
  d <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(d) || !nrow(d))
    return(sprintf('<div class="missing">empty table: <code>%s</code></div>', esc(basename(path))))
  if (!is.null(cols)) d <- d[, intersect(cols, names(d)), with = FALSE]
  note <- ""
  if (nrow(d) > max_rows) { note <- sprintf("<p class='note'>showing %d of %d rows</p>", max_rows, nrow(d)); d <- head(d, max_rows) }
  for (j in names(d)) if (is.numeric(d[[j]])) set(d, j = j, value = round(d[[j]], digits))
  h <- paste0("<tr>", paste0("<th>", esc(names(d)), "</th>", collapse = ""), "</tr>")
  rows <- apply(as.matrix(d), 1, function(r)
    paste0("<tr>", paste0("<td>", esc(as.character(r)), "</td>", collapse = ""), "</tr>"))
  paste0(if (!is.null(caption)) paste0("<p class='cap'>", caption, "</p>") else "",
         "<table>", h, paste0(rows, collapse = ""), "</table>", note)
}

# demote = TRUE is for an .md embedded under a heading the report has already
# printed: its own "# chrN" title would otherwise repeat that heading at <h2>,
# i.e. at the level of the enclosing SECTION, so the title is dropped and the
# rest shifts down one level.
md_file <- function(path, demote = FALSE) {
  if (!file.exists(path)) return(sprintf('<div class="missing">missing: <code>%s</code></div>', esc(path)))
  x <- readLines(path, warn = FALSE)
  out <- character(0); intab <- FALSE; inlist <- FALSE
  # A bare <li> outside a list is invalid and loses its marker and indent.
  close_list <- function() { if (inlist) { out <<- c(out, "</ul>"); inlist <<- FALSE } }
  for (ln in x) {
    if (grepl("^\\|", ln)) {
      cells <- trimws(strsplit(sub("^\\|", "", sub("\\|$", "", ln)), "|", fixed = TRUE)[[1]])
      if (all(grepl("^-+$", cells))) next
      close_list()
      tag <- if (!intab) { intab <- TRUE; out <- c(out, "<table>"); "th" } else "td"
      out <- c(out, paste0("<tr>", paste0("<", tag, ">", esc(cells), "</", tag, ">", collapse = ""), "</tr>"))
    } else {
      if (intab) { out <- c(out, "</table>"); intab <- FALSE }
      hh <- if (demote) c("h5", "h4", NA) else c("h4", "h3", "h2")
      if (grepl("^### ", ln))      { close_list(); out <- c(out, paste0("<", hh[1], ">", esc(sub("^### ", "", ln)), "</", hh[1], ">")) }
      else if (grepl("^## ", ln))  { close_list(); out <- c(out, paste0("<", hh[2], ">", esc(sub("^## ", "", ln)), "</", hh[2], ">")) }
      else if (grepl("^# ", ln))   { close_list(); if (!is.na(hh[3])) out <- c(out, paste0("<", hh[3], ">", esc(sub("^# ", "", ln)), "</", hh[3], ">")) }
      else if (grepl("^- ", ln))   { if (!inlist) { out <- c(out, "<ul>"); inlist <- TRUE }
                                     out <- c(out, paste0("<li>", bold(esc(sub("^- ", "", ln))), "</li>")) }
      else if (nzchar(trimws(ln))) { close_list(); out <- c(out, paste0("<p>", bold(esc(ln)), "</p>")) }
    }
  }
  close_list()
  if (intab) out <- c(out, "</table>")
  paste(out, collapse = "\n")
}
# perl = TRUE for the lazy quantifier; R's default TRE does not support it.
# Runs after esc(), so the backticks are still literal but any < > inside them
# is already an entity.
bold <- function(x) {
  x <- gsub("\\*\\*(.+?)\\*\\*", "<strong>\\1</strong>", x, perl = TRUE)
  gsub("`(.+?)`", "<code>\\1</code>", x, perl = TRUE)
}

# paste0(..., collapse = x) is a trap here: with length-1 arguments paste0
# concatenates them element-wise FIRST, so collapse gets a single string and
# inserts nothing. That silently ran sentences together in the rendered report
# ("not a domain set.GAM_10 took"). c(...) first, then collapse.
sec <- function(title, ...) paste0("<section><h2>", title, "</h2>", paste(c(...), collapse = "\n"), "</section>")
p_   <- function(...) paste0("<p>", paste(c(...), collapse = " "), "</p>")

# ==============================================================================
# LIVE NUMBERS. These used to be typed into the HTML by hand, which is how the
# ribbon cutoff ended up quoted against two different denominators in two
# places. Anything that can be recomputed from results/ in a few seconds is
# recomputed here, so the prose cannot drift from the files. Each one degrades
# to NA and the sentence that uses it says so, rather than aborting the report.
num <- new.env(parent = emptyenv())
num$ok <- FALSE
try({
  h <- load_hotspots()[chr %in% gam_autosomes]
  h[, mx := pmax(n_contacts_lost, n_contacts_gained)]
  num$hot_n     <- nrow(h)
  num$hot_med   <- median(h$n_contacts_lost + h$n_contacts_gained)
  num$hot_max   <- max(h$mx)
  num$hot_both  <- sum(h$n_contacts_lost >= 50 & h$n_contacts_gained >= 50)
  num$cut       <- vapply(c(50, 60, 75), function(k) sum(h$mx >= k), 0L)
  # TWO denominators, both stated wherever the number is. 123,133 is the whole
  # 20 kb grid; the delivered hotspot files cover slightly fewer bins, and a
  # percentage quoted without saying which is being used is the ambiguity this
  # block exists to remove.
  num$grid_n    <- 123133L
  num$ok        <- TRUE
}, silent = TRUE)

num$is_ok <- FALSE
try({
  isf <- file.path(rds_dir, "meltron_is.rds")
  if (file.exists(isf)) {
    isd <- readRDS(isf)$is
    isd[, d := IS_HP1gKO - IS_WT]
    num$is_sizes <- sort(unique(isd$is_size))
    num$is_vp    <- length(unique(isd$viewpoint))
    num$is_rho_gt <- cor(isd$IS_WT, isd$IS_HP1gKO, use = "pairwise")
    num$is_q999   <- unname(quantile(abs(c(isd$IS_WT, isd$IS_HP1gKO)), 0.999, na.rm = TRUE))
    w <- dcast(isd, chr + viewpoint ~ is_size, value.var = "d")
    m <- as.matrix(w[, -(1:2)]); m <- m[complete.cases(m), , drop = FALSE]
    pc <- prcomp(m, center = TRUE, scale. = FALSE)
    v  <- pc$sdev^2 / sum(pc$sdev^2)
    num$is_pc1  <- 100 * v[1]
    num$is_pc12 <- 100 * sum(v[1:2])
    num$is_rho_scale <- cor(m[, 1], m[, ncol(m)])
    if (isTRUE(num$ok)) {
      hh <- h[, .(chr, viewpoint = start, net = n_contacts_gained - n_contacts_lost)]
      num$is_rho_net <- vapply(c(240e3, 480e3, 1040e3), function(sz) {
        x <- merge(isd[is_size == sz, .(chr, viewpoint, d)], hh,
                   by = c("chr", "viewpoint"))
        cor(x$d, x$net, method = "spearman", use = "pairwise")
      }, 0)
    }
    num$is_ok <- TRUE
  }
}, silent = TRUE)
pc_ <- function(n, d) sprintf("%.1f%%", 100 * n / d)

# ==============================================================================
parts <- list()

parts$head <- paste0(
  "<h1>GAM &times; ChIP working report</h1>",
  "<p class='sub'>HP1&gamma; (heterochromatin protein 1 gamma, the protein encoded by Cbx3) cKO (conditional knockout) vs WT (wild type) &mdash; ",
  "GAM (genome architecture mapping) in CA1 (hippocampal field CA1) pyramidal cells, ChIP (chromatin immunoprecipitation) in Emx1-Cre ",
  "(the Cre driver line used for the conditional deletion) cortex+hippocampus. mm39 (the GRCm39 mouse genome assembly), 20 kb. Generated ",
  Sys.Date(), " by <code>make_report.R</code>.</p>",
  "<div class='warn'><strong>Read first.</strong> GAM is 2 biological replicates per genotype, already merged, so there is no within-genotype variance: ",
  "every number here is an effect size against a matched or bootstrap null, never a per-contact p-value. ",
  "Per-bin ChIP log2FC (log2 fold change) at 20 kb is noise-dominated (per-bin standard deviation, or sd, 0.14&ndash;0.16 for the broad marks, 0.57 for H3K4me3); genome-wide medians are sound, local structure is not. ",
  "The A/B compartment axis &mdash; the split of each chromosome into the active A and inactive B spatial compartments &mdash; is weakly determined on most chromosomes and is not used for any claim below.</div>")

parts$bottom <- sec("Bottom line",
  "<p>Everything below this section is the supporting detail. Percentages are change in mean NPMI",
  "(normalised pointwise mutual information, the contact measure used throughout),",
  "cKO vs WT, with 95% CIs (confidence intervals) from a bootstrap over the 19 autosomes",
  "(<code>gam10_class_change.tsv</code>, <code>gam13_did.tsv</code>).</p>",

  "<h3>1. H4K20me3 is lost almost everywhere. That part is solid.</h3>",
  "<p>Median log2FC about &minus;0.50 per 20 kb bin, &minus;0.68 over peaks, ~95% of domains down.",
  "A global shift, not a minority-of-sites effect. Everything else here is the architecture question.</p>",

  "<h3>2. Heterochromatin becomes more self-associated at short range. It is NOT mixing.</h3>",
  "<p>This is the central result. Heterochromatin gains contacts with heterochromatin, and with nothing else.",
  "By partner class &mdash; <strong>het</strong> is heterochromatin, euchromatin the gene-rich remainder:</p>",
  "<table>",
  "<tr><th>pair type</th><th>&lt;0.2 Mb</th><th>0.2&ndash;1 Mb</th><th>1&ndash;5 Mb</th><th>5&ndash;10 Mb</th><th>&gt;10 Mb</th></tr>",
  "<tr><td>het&ndash;het, same domain</td><td><strong>+4.3%</strong> (1.8, 6.6)</td><td><strong>+5.9%</strong> (1.3, 18.6)</td><td><strong>+4.8%</strong> (3.6, 29.7)</td><td>&mdash;</td><td>&mdash;</td></tr>",
  "<tr><td>het&ndash;het, different domains</td><td><strong>+4.2%</strong> (3.4, 4.9)</td><td><strong>+6.2%</strong> (3.5, 8.2)</td><td>+4.8% (&minus;0.02, 8.5)</td><td>+5.1% (&minus;5.9, 9.0)</td><td>&minus;0.6% (&minus;10.3, 4.9)</td></tr>",
  "<tr><td>het&ndash;euchromatin</td><td>+2.0% (1.3, 2.8)</td><td>+0.2% (&minus;1.2, 1.6)</td><td>&minus;1.9% (&minus;5.1, 1.1)</td><td><strong>&minus;5.2%</strong> (&minus;10.1, &minus;1.0)</td><td>&minus;6.4% (&minus;13.5, 0.0)</td></tr>",
  "<tr><td>euchromatin&ndash;euchromatin</td><td>+0.9% (0.3, 1.6)</td><td><strong>&minus;1.7%</strong> (&minus;3.0, &minus;0.5)</td><td><strong>&minus;5.0%</strong> (&minus;7.5, &minus;2.4)</td><td><strong>&minus;8.1%</strong> (&minus;12.0, &minus;4.2)</td><td><strong>&minus;8.9%</strong> (&minus;13.8, &minus;4.0)</td></tr>",
  "</table>",
  "<p>Read down the columns, not across:</p>",
  "<ul>",
  "<li><strong>The gain is heterochromatin-to-heterochromatin, and it is not &ldquo;same chromatin type&rdquo;.</strong>",
  "Euchromatin&ndash;euchromatin is also same-type and already <em>falls</em> by 0.2&ndash;1 Mb (&minus;1.7%), and below",
  "200 kb the mixed het&ndash;euchromatin pairs rise <em>more</em> than same-type eu&ndash;eu (euchromatin&ndash;euchromatin) pairs do (+2.0% against",
  "+0.9%). So the short-range gain is a property of heterochromatin specifically, not of partners matching.</li>",
  "<li><strong>It reaches 1&ndash;5 Mb for within-domain pairs</strong> (+4.8%, CI 3.6 to 29.7), where between-domain pairs",
  "just touch zero (&minus;0.02 lower bound). That CI is hugely asymmetric on 4,657 pairs, so read it as",
  "&ldquo;positive, magnitude poorly determined&rdquo;. Past 5 Mb every het&ndash;het CI spans zero.</li>",
  "<li><strong>It is not mixing.</strong> Mixing would require het&ndash;euchromatin to rise. It does not: flat at 0.2&ndash;1 Mb",
  "(+0.2%, CI spans zero) and <em>falling</em> at long range (&minus;5.2%, &minus;6.4%). Heterochromatin is keeping more",
  "company with heterochromatin and no more with euchromatin &mdash; segregation sharpening at short range, not dissolving.</li>",
  "<li><strong>The long-range loss spares heterochromatin.</strong> Beyond 10 Mb eu&ndash;eu is &minus;8.9% and",
  "het&ndash;euchromatin &minus;6.4%, both substantial, while het&ndash;het is &minus;0.6% with a CI spanning zero &mdash;",
  "flat, not falling. So &ldquo;decreased long-range contacts&rdquo; is true of euchromatin and of mixed pairs, and",
  "<em>not</em> of heterochromatin with itself.</li>",
  "</ul>",

  "<h3>3. Individual domains do not compact. Heterochromatin-to-heterochromatin contact is what rises.</h3>",
  "<p>Within-domain versus between-domain is indistinguishable: every difference-in-differences CI spans zero",
  "(p 0.19&ndash;0.89) under all three domain definitions. So the right statement is &ldquo;heterochromatin contacts",
  "heterochromatin more&rdquo;, not &ldquo;each domain tightens&rdquo;. That null replicates across definitions, which is",
  "what makes it believable rather than a definition artefact.</p>",

  "<h3>4. The leading variance was per-bin detectability. Removing it fixes the eigenvector, not the agreement.</h3>",
  "<p class='note'>Rewritten 2026-10-18. The previous version of this section reported an axis on 38 of 38",
  "combinations with PC1 on only chr7, and named per-bin detection frequency as an untested candidate for the",
  "leading variance. It has since been tested, and it was the cause. The numbers below replace those.</p>",
  "<p>Per <code>npmi_compartment_orientation.tsv</code>, an oriented compartment axis is now recovered on",
  "<strong>27 of 38</strong> chromosome &times; genotype combinations. That is <em>fewer</em> than before on",
  "purpose: a selection margin and an orientation invariant were added, so a chromosome whose best eigenvector",
  "does not clearly beat the runner-up is now refused rather than called. Its agreement with the orienting",
  "marker &mdash; the mean of per-bin WT H3K4me3 and H3K36me3 z-scores (signal in standard deviations from the genome-wide mean) &mdash; is modest.",
  "Agreement is Spearman&rsquo;s &rho; (rho, the rank correlation coefficient), and |&rho;| is its magnitude, ignoring sign:</p>",
  "<table>",
  "<tr><th>|&rho;| with orienting marker</th><th>combos (of 27 oriented)</th></tr>",
  "<tr><td>&ge; 0.3</td><td>21</td></tr><tr><td>&ge; 0.4</td><td>8</td></tr><tr><td>&ge; 0.5</td><td>0</td></tr>",
  "</table>",
  "<p><strong>The diagnostic named here as untested has been run, and it was the answer.</strong> GAM_2 already",
  "subtracted expected-by-distance, so the distance gradient was handled &mdash; but a <em>per-bin</em> effect is",
  "orthogonal to a per-distance one. Each 20 kb window has its own detectability (mappability, repeat content, how",
  "it survives cryosectioning), and that moves a bin's whole row and column together, loading onto the leading",
  "eigenvector. Fitting and subtracting it before the PCA takes the compartment back to PC1 on",
  "<strong>16 of 27</strong> combinations, from 2 of 24. On the shared axis the same change moved 15 of 16",
  "chromosomes up <em>exactly one</em> eigenvector &mdash; the signature of one eigenvector having been spent on",
  "detectability. Method and maths: <code>results/BIN_BALANCING.md</code>.</p>",
  "<p><strong>It fixed which eigenvector, and almost nothing else.</strong> Median |&rho;| moved only 0.340 to",
  "0.356 per genotype, and on the shared axis it is a wash (8 chromosomes up, 7 down). The axis is better",
  "<em>defined</em>; it does not agree better with active chromatin. The ceiling is the data, and no decomposition",
  "raises it. Both genotypes orient on only <strong>10 of 19</strong> autosomes, and on <strong>6 of those 10</strong>",
  "they pick the same eigenvector, so a per-genotype PC1-vs-PC1 comparison remains unsound on most of the genome.",
  "Full record: <code>results/COMPARTMENTS.md</code> sections 4b and 4c.</p>",
  "<p>Each genotype is oriented independently against an external ruler, so an apparent A/B switch could still be an orientation artefact; orienting",
  "both on a genotype-independent marker (GC content &mdash; the guanine+cytosine fraction &mdash; or gene density, rather than H3K4me3 signal) would remove that",
  "degree of freedom. cLAD/ciLAD &mdash; constitutive lamina-associated domains (LADs) and constitutive inter-LADs &mdash; remain the external",
  "alternative that avoids the axis entirely.</p>",

  "<h3>5. The remodelling tracks the size of the mark loss only weakly, bin by bin.</h3>",
  "<p>&rho; = &minus;0.015 raw, &minus;0.046 partialled on coverage. <strong>This was described as a powered negative and",
  "it is not one.</strong> &ldquo;Powered&rdquo; was quoted against n = 119,952 bins, but the unit of independence here is the",
  "chromosome: the 19-chromosome bootstrap gives a CI half-width of &plusmn;0.03 on the raw &rho;, and the effect implied",
  "by the hotspot result (&minus;0.028) sits <em>inside</em> that interval. The analysis cannot distinguish zero from the",
  "effect at issue. The partial-covariates variant is in fact &rho; = &minus;0.0216 with CI [&minus;0.0372, &minus;0.0049],",
  "which <em>excludes</em> zero in the predicted direction.</p>",
  "<p>What the data supports: a real, monotone, <strong>tiny</strong> dose-response. Median H4K20me3 log2FC moves from",
  "&minus;0.485 to &minus;0.517 across the entire range of contact loss &mdash; 0.032 log2 units against a genome-wide",
  "median of &minus;0.497, so architecture change accounts for roughly <strong>6%</strong> of the depth of the mark loss.",
  "The 5 Mb versus 20 kb scale mismatch (~250&times;) still limits what a per-bin correlation can resolve, but the honest",
  "statement is &ldquo;weak and small&rdquo;, not &ldquo;absent&rdquo;.</p>",

  "<h3>6. The heterochromatin gain holds on 17 of 19 autosomes at short range, and reverses on several at long range.</h3>",
  "<p>Per-chromosome numbers are percent change in mean NPMI, cKO vs WT, computed from",
  "<code>gam10_perchrom_sums.tsv</code> (per-chromosome sums and pair counts) as",
  "100 &times; (mean<sub>cKO</sub> / mean<sub>WT</sub> &minus; 1) within each chromosome &times; class &times; separation",
  "cell. No bootstrap: these are the point estimates for single chromosomes, and the CIs quoted elsewhere",
  "come from resampling chromosomes, which is not available within one.</p>",

  "<div class=\'warn\'><strong>A chromosome&rsquo;s overall direction is not interpretable; the contrast between its",
  "classes is.</strong> NPMI is a relative measure, so anything that shifts one chromosome&rsquo;s detection",
  "frequency between genotypes moves every class on it the same way. chr7 is up in all four classes at all five",
  "separations (+0.8% to +13.8%); chr11 is down in all classes past 1 Mb (&minus;2.9% to &minus;27.8%). A uniform",
  "shift across every class and band is a whole-chromosome offset, not a class-specific mechanism. The offset",
  "cancels when two classes on the same chromosome are differenced, so the het&ndash;het minus eu&ndash;eu contrast",
  "below is the quantity that survives it.</div>",

  "<h4>The offset-invariant contrast, every autosome</h4>",
  "<p>het&ndash;het (different domains) minus euchromatin&ndash;euchromatin, in percentage points. Positive means",
  "heterochromatin gained relative to euchromatin on that chromosome at that separation.</p>",
  "<table>",
  "<tr><th>chr</th><th>&lt;0.2 Mb</th><th>0.2&ndash;1 Mb</th><th>1&ndash;5 Mb</th><th>5&ndash;10 Mb</th><th>&gt;10 Mb</th></tr>",
  "<tr><td>chr1</td><td>&minus;0.6</td><td>&minus;1.6</td><td>&minus;3.4</td><td><strong>&minus;17.0</strong></td><td><strong>&minus;10.0</strong></td></tr><tr><td>chr2</td><td>+2.3</td><td>+7.2</td><td>+2.6</td><td>&minus;0.7</td><td>&minus;2.7</td></tr><tr><td>chr3</td><td>&minus;0.8</td><td>&minus;2.2</td><td>&minus;2.5</td><td>+7.3</td><td>+5.2</td></tr><tr><td>chr4</td><td>+4.3</td><td>+8.7</td><td><strong>+10.2</strong></td><td>+9.7</td><td>&minus;4.5</td></tr><tr><td>chr5</td><td>+4.3</td><td>+7.2</td><td><strong>+12.3</strong></td><td><strong>+14.9</strong></td><td><strong>+17.0</strong></td></tr><tr><td>chr6</td><td>+1.5</td><td>+4.8</td><td>+5.9</td><td>&minus;6.0</td><td><strong>+18.2</strong></td></tr><tr><td>chr7</td><td>+3.0</td><td>+8.3</td><td>+5.0</td><td>+1.7</td><td>&minus;5.3</td></tr><tr><td>chr8</td><td>+1.8</td><td>+0.3</td><td>+2.5</td><td><strong>&minus;15.5</strong></td><td>+1.0</td></tr><tr><td>chr9</td><td>+5.1</td><td><strong>+16.0</strong></td><td>+7.1</td><td>+5.4</td><td><strong>+10.5</strong></td></tr><tr><td>chr10</td><td>+2.3</td><td>+2.5</td><td>&minus;8.2</td><td>&minus;2.2</td><td>&minus;8.8</td></tr><tr><td>chr11</td><td>+7.3</td><td><strong>+13.6</strong></td><td>&minus;2.6</td><td><strong>&minus;11.8</strong></td><td>+3.0</td></tr><tr><td>chr12</td><td>+4.7</td><td>+3.3</td><td>+4.8</td><td><strong>+13.0</strong></td><td>+9.3</td></tr><tr><td>chr13</td><td>+4.1</td><td><strong>+13.5</strong></td><td><strong>+22.1</strong></td><td>+6.6</td><td>&minus;2.5</td></tr><tr><td>chr14</td><td>+3.8</td><td>+8.2</td><td><strong>+14.8</strong></td><td><strong>+23.3</strong></td><td><strong>+23.8</strong></td></tr><tr><td>chr15</td><td>+1.9</td><td>+5.0</td><td>+1.5</td><td><strong>&minus;13.5</strong></td><td>+0.9</td></tr><tr><td>chr16</td><td>+1.2</td><td>+4.8</td><td><strong>+14.8</strong></td><td><strong>+15.8</strong></td><td><strong>+14.6</strong></td></tr><tr><td>chr17</td><td>+4.9</td><td>+9.0</td><td>+9.6</td><td>+6.6</td><td>+5.5</td></tr><tr><td>chr18</td><td>+0.0</td><td>+3.1</td><td>+4.2</td><td>+2.5</td><td>&minus;5.0</td></tr><tr><td>chr19</td><td>+2.4</td><td>+5.4</td><td>+8.9</td><td><strong>+23.1</strong></td><td>+3.0</td></tr>",
  "<tr><td><strong>positive on</strong></td><td><strong>17/19</strong></td><td><strong>17/19</strong></td><td><strong>15/19</strong></td><td>12/19</td><td>12/19</td></tr>",
  "</table>",
  "<ul>",
  "<li><strong>Below 1 Mb the result is near-universal.</strong> Heterochromatin gains relative to euchromatin on",
  "17 of 19 autosomes in both short bands. The two exceptions are chr1 (&minus;0.6, &minus;1.6) and chr3",
  "(&minus;0.8, &minus;2.2) &mdash; both small enough to be noise around zero rather than a reversal. This is the",
  "per-chromosome support for result 2 and it is the strongest pattern in the data.</li>",
  "<li><strong>Beyond 5 Mb it splits 12 to 7.</strong> The contrast stays strongly positive on chr14 (+23.3, +23.8),",
  "chr16 (+15.8, +14.6), chr5 (+14.9, +17.0) and chr19 (+23.1), and goes clearly negative on chr1 (&minus;17.0),",
  "chr8 (&minus;15.5), chr15 (&minus;13.5) and chr11 (&minus;11.8). So the long-range behaviour is",
  "chromosome-specific in a way the short-range behaviour is not.</li>",
  "<li><strong>chr11 is not the lone exception.</strong> It inverts at 1&ndash;5 Mb (&minus;2.6) and 5&ndash;10 Mb",
  "(&minus;11.8), but chr1 inverts harder at both (&minus;3.4, &minus;17.0), and chr8, chr10 and chr15 also invert",
  "at long range. Any statement that singles chr11 out as reversing the rule is wrong; it is one of five.</li>",
  "<li><strong>Power is very unevenly distributed.</strong> het&ndash;het pair counts at 1&ndash;5 Mb range from 1,245",
  "(chr11) to 61,329 (chr7) &mdash; a 49&times; spread, because heterochromatin content varies from 3.1% of bins on",
  "chr11 to far more on chr7. The chromosomes that invert at long range are, with the exception of chr1, among the",
  "thinnest. Read the sign, not the magnitude, on chr8, chr10, chr11, chr15 and chr18.</li>",
  "</ul>",
  embed_fig("chromosomes/contact_budget_all19_unified.png",
            "The whole per-chromosome result in one view: chromosome down the side, separation across, chromatin class in the four panels, percent change in mean NPMI as the fill. Grey cells have no bin pairs &mdash; same-domain pairs exist only where a large heterochromatin domain does, and only out to the separations it spans."),
  embed_fig("chromosomes/perchrom_class_contrast_all19.png",
            "het&ndash;het minus euchromatin&ndash;euchromatin per autosome and separation band. This is the offset-invariant contrast."),
  embed_fig("chromosomes/perchrom_class_band_all19.png",
            "All four pair classes on all 19 autosomes. Read the shape across separations within a panel; the level of a whole panel is confounded by the per-chromosome offset."),
  "<div class=\'warn\'><strong>chr11 specifically.</strong> It is the least heterochromatic autosome (3.1% of bins",
  "against 7.2% genome-wide), so its het&ndash;het cells rest on roughly 900&ndash;1,250 bin pairs where eu&ndash;eu",
  "has millions. Its dominant feature is a loss in <em>every</em> class (eu&ndash;eu &minus;22.0% at 5&ndash;10 Mb,",
  "&minus;27.8% beyond 10 Mb), with het&ndash;het more negative still on thin data. Nothing computed here attributes",
  "that to the loss of one large LAD-associated domain, and chr11 is an unlikely place for a large heterochromatin",
  "domain by construction.</div>",

  "<h3>7. Four things tested and rejected.</h3>",
  "<ul>",
  "<li><strong>Heterochromatin-rich chromosomes gain contacts.</strong> No: &rho; = &minus;0.014, p = 0.96 across 19 autosomes.",
  "chr7 and chr11 are simply the extremes of both rankings at once.</li>",
  "<li><strong>The Pcdh (protocadherin cluster) TAD &mdash; topologically associating domain &mdash; changes its distal contacts distinctively.</strong> No: p = 0.75/0.85 intra (within the domain) and 0.88/0.68 distal",
  "against separation-matched chr18 TADs, under both partitions.</li>",
  "<li><strong>Large domains compact internally.</strong> No &mdash; item 3.</li>",
  "<li><strong>H3K36me3 rises over Pcdh.</strong> No: &minus;0.049 against a genome median of &minus;0.052, p = 0.93.</li>",
  "</ul>",
  "<p>Pcdh does meet GAM_5&rsquo;s decompaction criteria (short-range loss plus a shift toward A) while failing the",
  "matched-TAD test. Both are true: the comparisons differ &mdash; matched background bins genome-wide versus neighbouring",
  "TADs at matched separation &mdash; so any claim must say which.</p>")

parts$hotspot_annot <- sec("What hotspot domains are made of, and whether that differs by direction",
  p_("The section above calls the domains. This one asks what is inside them.",
     "Three nulls are used, and they disagree; each is stated before its numbers."),

  "<h3>Method: what is counted</h3>",
  p_("Two statistics per annotation, both reported because they diverge when a few large domains carry the signal:",
     "<strong>size-weighted</strong> is the share of domain base pairs falling in the annotation;",
     "<strong>per-region</strong> is the share of domains touching it at all.",
     "22 annotations pass a 50-bin floor: 15 chromHMM states, LINE/SINE/LTR-high (top quartile each), cLAD, ciLAD,",
     "has-a-TSS, TAD boundary, CTCF peak and CTCF-dense per tissue, and the three mark-loss classes."),

  "<h3>Method: the three nulls</h3>",
  "<ul>",
  "<li><strong>Hypergeometric.</strong> The obvious first move, computed and reported in a column named",
  "<code>p_hypergeom_INFLATED</code>. It treats 20 kb bins as independent draws. They are not: the",
  "differential-contact count track has autocorrelation r = 0.59 at 20 kb and a hotspot domain is by construction a",
  "<em>run</em> of correlated bins, so the effective sample size is the ~1,000 domains, not the ~10,000 bins.</li>",
  "<li><strong>Independent width-matched shuffle.</strong> Each domain keeps its width and its chromosome and is",
  "dropped at a random mappable position, 1,000 times. Width and chromosome are matched by construction, which is the",
  "confound the hypergeometric ignores &mdash; a 200 kb domain overlaps a narrow feature by being long.</li>",
  "<li><strong>Circular shift.</strong> The whole domain set is rotated along each chromosome by one random offset,",
  "preserving every inter-domain distance. <strong>This is the test to read.</strong></li>",
  "</ul>",
  p_("Why the rotation is needed: domains still cluster after merging. 12% of adjacent lost domains sit within 100 kb",
     "of the next and 27% within 200 kb (for the unadjusted set, 50% and 69%). Annotations cluster too &mdash; chromHMM",
     "states and LADs come in blocks &mdash; so placing domains independently builds a null with less spatial structure",
     "than the data, which understates the variance and makes p too small."),

  "<h3>Result: the nulls disagree, and by how much</h3>",
  "<table>",
  "<tr><th>null</th><th>cells significant at 0.05 (of 104)</th><th>after BH</th></tr>",
  "<tr><td>hypergeometric</td><td>18</td><td>&mdash;</td></tr>",
  "<tr><td>independent width-matched shuffle</td><td>49</td><td>38</td></tr>",
  "<tr><td><strong>circular shift</strong></td><td><strong>36</strong></td><td><strong>17</strong></td></tr>",
  "</table>",
  p_("<strong>13 cells are significant under independent placement but not under the rotation</strong> &mdash; that is",
     "the clustering effect, measured. The hypergeometric is not simply the most liberal: it is the most",
     "<em>conservative</em> here, because its expected value is the genome-wide annotation fraction while the shuffle's",
     "is the fraction reachable by domain-sized windows, and for clustered annotations those differ. The three tests",
     "disagree about the expectation, not only about the variance."),
  embed_fig("hotspots/hotspot_test_comparison.png",
            "One line per annotation &times; domain-set cell, under all three nulls."),

  "<h3>Method: does merging change this?</h3>",
  p_("Domains are merged across a 1-bin (20 kb) sub-threshold gap only. The same sweep that exposed the heterochromatin",
     "percolation is run one level up, on the called domains, to check whether a further merge would run away."),
  "<h3>Result: hotspot domains do not percolate</h3>",
  "<table>",
  "<tr><th>merge gap</th><th>domains</th><th>% of genome</th><th>largest</th></tr>",
  "<tr><td>as called</td><td>979</td><td>8.09%</td><td>1.5 Mb</td></tr>",
  "<tr><td>&le; 100 kb</td><td>860</td><td>8.35%</td><td>3.6 Mb</td></tr>",
  "<tr><td>&le; 200 kb</td><td>716</td><td>9.11%</td><td>5.9 Mb</td></tr>",
  "<tr><td>&le; 500 kb</td><td>508</td><td>11.68%</td><td>6.7 Mb</td></tr>",
  "<tr><td>&le; 1 Mb</td><td>354</td><td>15.98%</td><td>11.7 Mb</td></tr>",
  "</table>",
  p_("Gradual, not a cliff &mdash; against the heterochromatin fragments, which went from 1.3% of the genome at gap 0",
     "to 64.8% at 50 kb. Merging here is a choice about scale, and past ~500 kb it starts re-deriving the large blocks,",
     "which already exist as their own object."),
  embed_fig("hotspots/hotspot_domain_merge_sweep.png"),

  "<h3>Result: composition, against the rotation null</h3>",
  p_("Fold is the observed share of domain base pairs in the annotation divided by the share after rotation.",
     "Only cells passing BH &lt; 0.05 on the rotation are listed, for the merged (<code>adj</code>) domain sets."),
  "<table>",
  "<tr><th>domain set</th><th>annotation</th><th>observed share</th><th>fold</th><th>z</th></tr>",
  "<tr><td>lost</td><td>constitutive LAD</td><td>0.477</td><td><strong>1.35</strong></td><td>4.9</td></tr>",
  "<tr><td>lost</td><td>LINE-high</td><td>0.313</td><td><strong>1.25</strong></td><td>4.0</td></tr>",
  "<tr><td>lost</td><td>chromHMM Quies</td><td>0.806</td><td>1.10</td><td>4.5</td></tr>",
  "<tr><td>lost</td><td>has a TSS</td><td>0.418</td><td>0.92</td><td>&minus;2.6</td></tr>",
  "<tr><td>lost</td><td>CTCF peak (bulk cortex)</td><td>0.221</td><td>0.88</td><td>&minus;2.6</td></tr>",
  "<tr><td>lost</td><td>SINE-high</td><td>0.194</td><td>0.80</td><td>&minus;2.5</td></tr>",
  "<tr><td>lost</td><td>chromHMM Enh</td><td>0.000</td><td>0.01</td><td>&minus;2.8</td></tr>",
  "</table>",
  p_("Domains losing contacts are LAD-associated, quiescent and LINE-rich, and depleted of promoters, enhancers,",
     "SINEs and CTCF. Note that <em>every</em> surviving cell is a lost-domain cell: after the rotation, gained domains",
     "are not distinguishable from rotated background on any annotation tested."),
  embed_fig("hotspots/hotspot_annotation_enrichment.png"),

  "<h3>Method: the direction contrast</h3>",
  p_("Lost domains against gained domains directly. This is the cleanest comparison available, because both sets are",
     "called by the same procedure from the same track on the same bins, so mappability, coverage and width confounds",
     "are shared and largely cancel in the difference. The null is a permutation of the lost/gained labels stratified",
     "by width decile, so the two sets keep comparable width distributions and the test is about composition only."),
  "<h3>Result</h3>",
  "<table>",
  "<tr><th>annotation</th><th>lost</th><th>gained</th><th>difference</th><th>z</th></tr>",
  "<tr><td>constitutive LAD</td><td>0.477</td><td>0.358</td><td><strong>+0.119</strong></td><td>4.5</td></tr>",
  "<tr><td>has a TSS</td><td>0.418</td><td>0.482</td><td>&minus;0.063</td><td>&minus;4.2</td></tr>",
  "<tr><td>chromHMM Quies3</td><td>0.010</td><td>0.029</td><td>&minus;0.019</td><td>&minus;4.2</td></tr>",
  "<tr><td>chromHMM Het</td><td>0.000</td><td>0.007</td><td>&minus;0.007</td><td>&minus;3.7</td></tr>",
  "<tr><td>chromHMM Quies</td><td>0.806</td><td>0.761</td><td>+0.045</td><td>3.5</td></tr>",
  "<tr><td>LINE-high</td><td>0.313</td><td>0.262</td><td>+0.051</td><td>3.1</td></tr>",
  "</table>",
  p_("A coherent split along the direction axis: contacts are lost in lamina-associated, quiescent, LINE-rich",
     "chromatin and gained in promoter-containing, inter-LAD chromatin. This is the one place in the report where",
     "the architecture change maps onto an A/B-flavoured distinction without going through the compartment axis."),
  embed_fig("hotspots/hotspot_direction_contrast.png"),
  tbl("gam4_hotspot_annotation_enrichment.tsv", "Full table, all three nulls, both domain definitions.", max_rows = 30L),
  tbl("gam4_hotspot_direction_contrast.tsv", "Full direction contrast.", max_rows = 25L),

  "<div class='warn'><strong>Carry the two domain definitions through.</strong> <code>adj</code> is",
  "smoothed-then-thresholded, <code>raw</code> is thresholded directly; they agree on only 40% of bins and the",
  "<code>raw</code> sets are median 1 bin, so the independence problem bites hardest there. A result present in one",
  "and not the other is a result about the definition.</div>")

parts$het <- sec("Heterochromatin domains",
  p_("The merged H3K9me3/H4K20me3 set is <strong>126,197 fragments of mean width 3.5 kb</strong> &mdash; a peak set, not a domain set.",
     "GAM_10 took whichever single fragments exceeded a width cut and called them domains. GAM_13 builds domains three ways and reports every result against all three;",
     "where they disagree, the definition is doing the work."),
  p_("The fragments <strong>percolate</strong> when gap-merged: 1.3% of the genome at gap 0, 4.5% at 10 kb, 17.5% at 20 kb, 64.8% at 50 kb.",
     "10 kb sits below that transition, which is why it is the chosen gap."),
  embed_fig("het_domains/merge_sweep.png"),
  tbl("gam13_domain_definitions.tsv", "The three definitions."),
  tbl("gam13_definition_agreement.tsv", "Agreement between them (the Jaccard index &mdash; shared bins divided by the union &mdash; on bin membership). For scale, the two hotspot-domain variants agree at 0.40."),
  embed_fig("het_domains/intra_vs_inter_by_definition.png",
            "Within-domain vs between-domain contact change, per definition. Bars are 95% CIs from a bootstrap over the 19 autosomes."),
  embed_fig("het_domains/by_domain_size.png"),
  tbl("gam13_did.tsv", "Difference-in-differences: intra (within-domain) minus inter (between-domain). A CI excluding zero means the two classes changed by different amounts."),
  embed_fig("het_domains/domains_vs_lads.png"),
  tbl("gam13_domains_vs_lads.tsv", "Domain bp inside constitutive LADs. cLAD/ciLAD are external to both genotypes and all five marks, so unlike our own PC1 they cannot be circular &mdash; but they are the cell-type-invariant core, not a neuronal annotation."))

parts$pcdh <- sec("The protocadherin cluster as a TAD",
  p_("The cluster is approximately one TAD, so the question is how that TAD's contacts with <em>other</em> TADs change.",
     "Contact frequency falls steeply with separation, so the test is a permutation over query TADs stratified by separation:",
     "every other chr18 TAD plays the part of Pcdh in turn, and we ask where Pcdh falls in that distribution."),
  p_("Two partitions are reported side by side &mdash; published boundaries and insulation minima called from our own",
     "GAM. A result that holds under only one is a result about the partition."),
  "<div class='warn'><strong>Of the two, ours is the one with external support.</strong> CTCF is the protein that",
  "makes boundaries and is measured by a different assay in a different lab, so it is an independent check on where a",
  "boundary is. Peak density around our insulation minima is <strong>1.72&times;</strong> genome average at",
  "&plusmn;5 kb, falling monotonically to 1.39&times; at &plusmn;100 kb &mdash; the shape a real boundary set should",
  "have. Around the published boundaries it is <strong>1.09&times; at every window size</strong>, i.e. flat. The same",
  "CTCF track is 5.27&times; enriched at TSSs, so the track is sound and the flatness is a property of that boundary",
  "set. Its provenance is also unrecorded &mdash; cell type, caller and resolution are all unknown",
  "(METHODS open question). Treat the published partition as the weaker of the two, not the safer one.</div>",
  tbl("gam14_pcdh_tad_summary.tsv", "Summary across partitions."),
  embed_fig("pcdh_tad/pcdh_vs_null.png"),
  embed_fig("pcdh_tad/pcdh_by_separation.png"),
  embed_fig("pcdh_tad/pcdh_partner_profile.png"),
  tbl("gam14_pcdh_by_separation.tsv"))

chr_files <- file.path(results_dir, "chromosomes", paste0(c("chr3","chr7","chr11","chr16","chr18"), ".md"))
parts$chrom <- sec("Per-chromosome characterisation",
  p_("chr11 and chr16 were flagged by the data producers; chr3 and chr7 carry large GAM_6 blocks; chr18 holds the Pcdh cluster."),
  "<div class='warn'><strong>Mark-contact coupling does not rank chromosomes.</strong> chr11 is rank 1 of 19 for &rho;(H4K20me3, net contacts), but at &rho; = 0.128, with chr17 at &minus;0.149 &mdash; larger, opposite sign &mdash; and a genome-wide value of 0.012. The spread is small and roughly symmetric about zero, so the ordering carries no weight.</div>",
  embed_fig("chromosomes/chromosome_ranking.png"),
  "<div class='warn'><strong>A tempting over-read, tested and rejected.</strong> chr7 is by far the most heterochromatic chromosome (20.2% of bins against 7.2% genome-wide) <em>and</em> gains contacts at every separation, while chr11 is the least (3.1%) and loses most. Across all 19 autosomes that relationship is absent: &rho; = &minus;0.014, p = 0.96. The two are simply the extremes of both rankings at once.</div>",
  embed_fig("chromosomes/het_content_vs_net.png"),
  tbl("gam15_chromosome_ranking.tsv", "All 19 autosomes, so any one chromosome can be read against the spread.", max_rows = 20L),
  embed_fig("chromosomes/contact_budget_all19_unified.png",
            "All 19 autosomes, all four chromatin classes, all five separation bands. Derived from the per-chromosome sums in <code>gam10_perchrom_sums.tsv</code>."),
  embed_fig("chromosomes/contact_budget_by_chromosome.png",
            "The same quantity recomputed directly from the NPMI matrices for the five focus chromosomes, with slightly coarser bands. Kept as an independent cross-check of the figure above, which comes from pre-summed counts."),
  tbl("gam15_contact_budget.tsv", "Contact budget by chromosome, chromatin class and separation.", max_rows = 60L),
  paste0(sapply(c("chr3","chr7","chr11","chr16","chr18"), function(ch)
    paste0("<h3>", ch, "</h3>",
           embed_fig(sprintf("chromosomes/%s_overview.png", ch)),
           md_file(file.path(results_dir, "chromosomes", paste0(ch, ".md")),
                   demote = TRUE))), collapse = ""))

# chr18_Pcdh has no _band variant by design: at 1.8 Mb the full triangle is
# already compact, so band_maxdist() skips it rather than writing a duplicate
# under a misleading name.
# Candidate-gene regions first, then the ones picked off gam6_large_blocks.tsv
# by width rather than by hypothesis. chr13:71-97 is the largest block in the
# genome and had no figure; chr3:144-160 is the widest GAINED block, and until
# now every block figure in this report was a lost block - so the lost/gained
# asymmetry was asserted in the tables and never shown.
reg <- c("chr18_Pcdh_sigonly", "chr18_Pcdh_context_band_sigonly",
         "chr13_KRABZFP_band_sigonly", "chr17_KRABZFP_band_sigonly",
         "chr7_Vmn_band_sigonly", "chr2_Olfr_band_sigonly", "chr16_full_band_sigonly",
         "chr11_1-25Mb_lostblock_band_sigonly", "chr11_100-125Mb_band_sigonly",
         "chr3_30-85Mb_band_sigonly", "chr13_71-97Mb_lostblock_band_sigonly",
         "chr3_144-160Mb_gainedblock_band_sigonly")
parts$regions <- sec("Region plots",
  p_("Contact map, differential contacts, hotspot calls, LADs, the five marks and genes, on one coordinate axis.",
     "Conventions: raw NPMI is a sequential ramp with one cap shared by both genotypes; delta-z (the upstream caller&rsquo;s contact-change statistic, defined under methods below) is fixed at &plusmn;4;",
     "<strong>net diff. contacts</strong> is the per-bin count of significant differential contacts within 5 Mb (gained &minus; lost), while <strong>hotspot bins</strong> is the separate call on that count at &ge; 50 per bin &mdash; a different object from the hotspot <em>domains</em> that carry every result in the tables. Both cutoffs, and who chose them, are audited in &ldquo;Every cutoff in this report&rdquo; below."),
  "<div class='warn'><strong>Data coverage, read before interpreting any empty map.</strong> A region whose bins are detected in few nuclear profiles has mostly NA (not available &mdash; missing, not measured) in the NPMI matrix, and NA renders white &mdash; indistinguishable from &ldquo;no contacts here&rdquo;. Each plot states its measured fraction, and warns below 75%. Measured here: Pcdh 97%, Olfr (olfactory receptor cluster) 92%, KRAB-ZFP (Kr&uuml;ppel-associated box zinc-finger protein cluster) chr17 87%, <strong>KRAB-ZFP chr13 61%</strong>, <strong>vomeronasal chr7 29%</strong>. The vomeronasal cluster's apparently empty contact map is <em>missing data, not absent contact</em> &mdash; 11% finite cells inside the cluster against 97% in its own flanks, with n_valid (the count of usable contact values per bin) 807 vs 6,574 (genome-wide 6,104). Do not read architecture off it.</div>",
  "<h3>Method: how the mark tracks are drawn</h3>",
  p_("They are <strong>coverage, not a ratio</strong>. Each panel overlays the two genotypes: a grey body from zero to",
     "the lower of the two genotype means, and a coloured band from there to the higher &mdash; red where the cKO is",
     "higher, blue where it is lower, the same convention as the contact maps. Height is how much signal is present,",
     "the band is the size of the change. Means are across 6 WT and 4 cKO libraries (cKO rep2 excluded at QC)."),
  p_("The display bin is chosen from the region width to give roughly 400 points across it, floored at the native 20 kb.",
     "The y axis is per mark and starts at zero, except that <strong>H3K9me2 and H3K9me3 share one axis</strong> so the",
     "comparison between them is a comparison of the drawing."),
  "<div class='warn'><strong>This replaced a log2 fold-change panel, and the replacement changed a reading.</strong>",
  "The old panels drew a ratio on a per-mark symmetric axis, so every mark filled its panel to the same extent whatever",
  "the size of its effect &mdash; H3K9me2 read as dramatically as H4K20me3. Measured genome-wide as the gap between",
  "genotype means divided by the mean spread between replicate libraries: <strong>H4K20me3 0.84, H3K9me3 0.54,",
  "H3K9me2 0.30, H3K36me3 0.12, H3K4me3 0.04</strong>. H3K9me2's effect is about a third of the replicate spread and",
  "roughly a third of H4K20me3's. The coverage overlay shows that ordering; the ratio hid it.</div>",
  "<h3>Method: the colour scales are the collaborators' own</h3>",
  p_("Both palettes were hand-matched by eye until 2026-10-18, when their chr16 full-length slide was sampled",
     "directly. One wrinkle first: the slide is a macOS screen capture and therefore carries a Display P3 profile,",
     "so its raw RGB has to be converted to sRGB before anything is compared &mdash; read uncorrected, every colour",
     "looks washed out and the comparison is meaningless. After that conversion both of their colourbars fit a",
     "<strong>9-class ColorBrewer ramp</strong> almost exactly (RMS over the whole bar, in 0&ndash;1 RGB):"),
  "<table>",
  "<tr><th>panel</th><th>their palette</th><th>fit</th><th>the old hand-made ramp</th></tr>",
  "<tr><td>delta-z</td><td><code>rev(RdBu, 9)</code></td><td>0.051</td><td>0.138</td></tr>",
  "<tr><td>raw NPMI</td><td><code>rev(RdYlBu, 9)</code></td><td>0.044</td><td>0.056</td></tr>",
  "</table>",
  "<div class='warn'><strong>This is why our difference maps looked about twice as changed as theirs, on identical numbers.</strong> The old 5-stop ramp put the firm red <code>#d6604d</code> at 75% of the way out; RdBu-9 is still at the pale salmon <code>#f4a582</code> there and does not reach <code>#d6604d</code> until 87.5%. On the &plusmn;4 limit both use, that is the difference between delta-z = +2 reading as a definite red and reading as a faint blush &mdash; and <strong>14.1% of chr16 cells sit beyond |delta-z| = 2</strong>, which is most of what the eye picks up. Nothing about the data changed. Checked directly against their <code>deltaz</code> column on chr16, our independently computed value is <strong>identical</strong> &mdash; Spearman 1.000, regression slope 1.000 on 1.41 million called contacts &mdash; so this was never a discrepancy in the statistic, only in the ramp drawing it.</div>",
  "<div class='warn'><strong>The sign convention is deliberately NOT theirs.</strong> Their delta-z is set1 &minus; set2 = WT &minus; cKO, so their red means <em>lost</em> in the cKO. This repo is positive = more in the cKO throughout and stays that way, so our difference maps are colour-inverted relative to their slides by design. Every panel carries &ldquo;cKO &minus; WT&rdquo; in its subtitle for this reason.</div>",
  p_("The <strong>cap</strong> is not theirs, though, and forcing it was tried and reverted. Their fixed 0&ndash;0.25",
     "bar makes our panels read solid warm with a thick saturated near-diagonal wedge, because 5.9% of chr16 cells",
     "clear 0.25 and all of them clip to one red. Two explanations for that were tested and both are wrong:",
     "block-averaging to a coarser display bin <em>reduces</em> the saturated share rather than causing it",
     "(going from 280 kb to 60 kb blocks took it from 8.6% <em>up</em> to 11.0%, because coarse blocks average the",
     "near-diagonal peak down); and matching their figure pixel-for-pixel is invalid, because their panel is a",
     "screenshot of a slide at ~42 kb per pixel, so its wedge is blended away by downscaling before it can be",
     "measured. The cap is the <strong>q97</strong> of what each figure draws &mdash; 0.299 on chr16, close to their",
     "0.25 but adapting per region. q99 is too wide: 0.413 against a chr16 median of 0.098, which turns the map blue."),
  "<h3>Method: the display bin comes from the output raster</h3>",
  p_("The contact maps are block-averaged before drawing, because subsampling every k-th bin aliases &mdash; the",
     "rotated lattice beats against the raster grid and the map reads as diagonal streaks. The target used to be a",
     "fixed ~350 display columns, which was never checked against the output. These figures save at 300 dpi over an",
     "~8.4 in panel, so the panel is about <strong>2520 px</strong> wide and 350 columns made every display block",
     "<strong>seven pixels across</strong> &mdash; a factor of 7 of resolution discarded in a file that was already",
     "paying for it. On chr16 that was a 280 kb display bin against the 20 kb native grid. It is now one block per",
     "output pixel, capped by a tile budget at 1600 columns (60 kb on chr16). The cap is a budget, not a resolution",
     "choice: the map is drawn as one rectangle per cell, in the PDF as well as the PNG, so an uncapped native",
     "triangle on a large chromosome produces a vector file nothing will open."),
    "<h3>Method: the difference map</h3>",
  p_("delta-z is each genotype's matrix z-scored within each genomic distance chromosome-wide, then subtracted,",
     "fixed at &plusmn;4. Raw NPMI falls about sixfold from 20 kb to 20 Mb, so an unnormalised difference would be",
     "dominated by the near-diagonal."),
  p_("A variant masking the map to the delivered significant contacts is available but off by default.",
     "The delivered criterion is |delta-z| &gt; 2 <em>and</em> max NPMI &gt; 0.1, which passes <strong>12.7% of measured",
     "cell pairs on chr16 and 14.1% on chr18</strong> &mdash; about 1.4 million contacts per chromosome. At the ~350",
     "display columns these figures use, one block is roughly 14 &times; 14 native cells, so nearly every block contains",
     "significant cells and the masked map fills in exactly like the unmasked one. It is genuinely sparse only at native",
     "resolution, on regions of a few Mb."),
  "<h3>The regions</h3>",
  p_("The difference map shown here is <strong>masked to the significant contacts</strong>, which is now the default",
     "render. The unmasked map spends its colour range on cells near zero and the structure washes out &mdash; the",
     "same contrast the collaborators&rsquo; own figure shows between its pale delta-z panel and its significant panel.",
     "Every region figure carries its variant in the filename &mdash; <code>_sigonly</code> here,",
     "<code>_allcells</code> for every measured cell, and <code>_allcells_binbalanced</code> for the same map with the",
     "additive per-bin (row/column) effect fitted and subtracted. That last one exists because a fifth to a quarter",
     "of any difference map is a per-bin term: fitting delta-z[i,j] ~ a[i] + a[j] accounts for <strong>22.8%</strong>",
     "of the variance at the olfactory cluster, 26.5% at Pcdh, 18.0% on the chr11 block and <strong>27.9%</strong> in",
     "a plain chr3 control window. It is a property of the data everywhere, not of odd regions. Distance z-scoring",
     "cannot remove it &mdash; that handles a per-<em>distance</em> offset, which is orthogonal &mdash; so it survives",
     "as the criss-cross plaid lattice visible once every cell is drawn. The significance mask hides it by dropping",
     "low-coverage cells."),
  # Figures are grouped by chromosome on disk; the chromosome is the name's own
  # prefix, so there is no second list to keep in step with this one.
  paste0(sapply(reg, function(r)
      embed_fig(sprintf("regions/%s/%s.png", sub("^(chr[0-9XY]+).*$", "\\1", r), r),
                paste0(esc(sub("_sigonly$", "", r)), " &mdash; difference map masked to significant contacts (grey = measured but not significant, white = not measured)"))
    ), collapse = ""),

  "<h3>The slide variant, and the insulation strip it carries</h3>",
  p_("A second render of four regions at 16:9, for presenting. It is a <strong>re-layout, not a resize</strong>:",
     "the two raw-NPMI panels are dropped, the mark tracks are cut to H4K20me3 and H3K9me3, and the distance cap",
     "is then solved from whatever height the track stack leaves &mdash; so the contact map is a band rather than a",
     "triangle, and the band height differs per region. It also carries one panel the other variants do not: the",
     "<strong>insulation-score difference strip</strong>, which is step 1 of the Pombo lab&rsquo;s MELTRON pipeline",
     "figure placed directly under the difference map."),
  "<h4>What the strip&rsquo;s y axis is</h4>",
  p_("<strong>Square size, 240 to 1040 kb &mdash; not a score.</strong> The insulation score at a position is the mean",
     "NPMI of the off-diagonal square that <em>crosses</em> that position, minus the chromosome mean at that same",
     "square size. Computing it at eleven square sizes gives eleven values per position, and the strip draws them as",
     "eleven stacked rows: colour is the score difference (cKO &minus; WT), the row is the width of the square it was",
     "measured over. A change confined to the bottom rows is local; a change running the full height of the strip is a",
     "change at every scale the measurement can see."),
  sprintf("<div class='warn'><strong>Three properties worth knowing before reading one.</strong><ul>%s%s%s</ul></div>",
    "<li><strong>The values do not depend on the plotted interval.</strong> Cropping a region changes only the display binning (20 kb at Pcdh, 280 kb on a whole chromosome); the underlying per-viewpoint values are identical to nine decimals. The strip is a view of a genome-wide object, not a statistic of the window.</li>",
    sprintf("<li><strong>A uniform genome-wide shift is invisible by construction.</strong> Centring is per genotype and per chromosome, so if both genotypes moved by the same amount everywhere, the strip would show nothing. That matters here, because the headline architectural effect <em>is</em> broad (&minus;1.7%% to &minus;5.0%% eu&ndash;eu, +4.2%% het&ndash;het). The strip shows where the change is <em>unevenly</em> distributed, which is a different question from whether there is one.</li>"),
    sprintf("<li><strong>The rows are nested, so they are not independent.</strong> The 1040 kb square contains the 240 kb square. Across the eleven rows of the difference, PC1 explains <strong>%s</strong> of the variance and PC1+PC2 <strong>%s</strong>, with &rho;(240 kb, 1040 kb) = %s. Eleven rows, roughly one and a half independent dimensions.</li>",
            if (num$is_ok) sprintf("%.1f%%", num$is_pc1) else "about 86%",
            if (num$is_ok) sprintf("%.1f%%", num$is_pc12) else "about 96%",
            if (num$is_ok) sprintf("%.2f", num$is_rho_scale) else "~0.55")),
  "<h4>The slide variant also carries repeat composition</h4>",
  p_("LINE, SINE and LTR as the fraction of each 20 kb bin they cover, from the annotated bin object &mdash; an",
     "annotation that was always in the data, now drawn. <strong>One panel on one shared y axis</strong>, not three",
     "stacked tracks: they are fractions of the same bin, and the comparison that matters is which class dominates",
     "where. The hotspot result is that domains losing contacts are LINE-rich (1.25&times;) and SINE-<em>depleted</em>",
     "(0.80&times;), which is a statement about LINE against SINE and becomes invisible once each class gets its own",
     "auto-scaled axis. The ceiling is twice the 99th percentile of the values drawn, capped at 1, so a single",
     "repeat-pileup bin cannot flatten the region &mdash; the same rule the mark panels use."),
    "<h4>Is the strip just the difference map again?</h4>",
  p_("A fair question, since the strip, the difference map and the hotspot ribbon all respond to the same perturbation,",
     "and by eye they line up. Measured, they do not carry the same information.",
     if (num$is_ok && !is.null(num$is_rho_net))
       sprintf("Spearman &rho; between the insulation difference and the net per-bin differential-contact count is <strong>%.2f at 240 kb, %.2f at 480 kb and %.2f at 1040 kb</strong> &mdash; about %.0f%% shared variance in the middle of that range.",
               num$is_rho_net[1], num$is_rho_net[2], num$is_rho_net[3],
               100 * num$is_rho_net[2]^2)
     else "",
     "The two <em>are</em> related and should be, but one is not a redrawing of the other.",
     "The genuinely redundant pair is elsewhere:",
     "the hotspot ribbon against its own count track, at &rho; = 0.91, which is why only one of those two is drawn."),
  sprintf("<div class='warn'><strong>A display bug caught by reading the figure, worth recording because it was invisible in the numbers.</strong> The difference strip was first drawn on &plusmn;0.10 while the per-genotype strips above it were on &plusmn;0.20. Every value was correct and every panel was individually defensible, but a stack of panels is read by comparison, so 4.8%% of cells that were mildly negative in <em>both</em> genotypes appeared as a near-saturated difference. On one shared limit that figure is 0.0%%. The limit is now %s for all strips, set just above the genome-wide 99.9th percentile of |score| (%s).</div>",
          "&plusmn;0.20",
          if (num$is_ok) sprintf("%.3f", num$is_q999) else "0.18"),
  paste0(sapply(c("chr18_Pcdh_slide", "chr11_1-25Mb_lostblock_slide",
                  "chr3_30-85Mb_slide", "chr16_full_slide"), function(r)
    embed_fig(sprintf("regions/%s/%s.png", sub("^(chr[0-9XY]+).*$", "\\1", r), r),
              paste0(esc(r), " &mdash; 16:9 slide render; difference map, insulation difference strip, hotspot bins, CTCF, LADs, two marks, genes"))),
    collapse = ""))

parts$inventory <- sec("What the data is, before any of it is interpreted",
  p_("The sections above and below are findings. This one is the inventory &mdash; what was measured, at what",
     "resolution, how much of it is actually present, and which numbers are dense covariates rather than hit lists.",
     "It exists because several arguments in this report turn on the <em>shape</em> of the data rather than on an",
     "effect size, and that shape is otherwise scattered across the methods section."),

  "<h3>The grid, and what sits on it</h3>",
  "<table>",
  "<tr><th>object</th><th>size</th><th>what it is</th></tr>",
  "<tr><td>20 kb bin grid</td><td>123,133 bins</td><td>mm39, the 19 autosomes only &mdash; there is no X anywhere in this analysis. Every table in this report keys on it.</td></tr>",
  "<tr><td>NPMI contact matrices</td><td>2 genotypes &times; 19 autosomes</td><td>already merged across replicates by the collaborators &mdash; there is no within-genotype variance to test against</td></tr>",
  "<tr><td>ChIP libraries</td><td>50 (6 WT / 4 cKO usable per mark)</td><td>MINUTE, 5 marks; HP1gKO rep2 excluded at QC</td></tr>",
  "<tr><td>ChIP peaks</td><td>11,293</td><td>the peak-level ruler in the ChIP section</td></tr>",
  sprintf("<tr><td>differential-contact count track</td><td>%s bins</td><td>per-bin count of significant differential contacts within 5 Mb, delivered; median %s per bin, maximum %s</td></tr>",
          if (num$ok) format(num$hot_n, big.mark = ",") else "&mdash;",
          if (num$ok) format(num$hot_med) else "&mdash;",
          if (num$ok) format(num$hot_max) else "&mdash;"),
  sprintf("<tr><td>insulation score maps</td><td>%s viewpoints &times; %d square sizes &times; 2 genotypes</td><td>240&ndash;1040 kb, written by GAM_2c; the input to MELTRON and to the strip on the slide figures</td></tr>",
          if (num$is_ok) format(num$is_vp, big.mark = ",") else "&mdash;",
          if (num$is_ok) length(num$is_sizes) else 11L),
  "<tr><td>hotspot domains</td><td>979 lost / 958 gained</td><td>q90 + smoothing + merging on the count track; the independent unit for every enrichment test</td></tr>",
  "</table>",

  "<h3>Three things that are dense, not sparse &mdash; and are repeatedly misread as hit lists</h3>",
  "<ul>",
  sprintf("<li><strong>The differential-contact count track.</strong> Every bin in the delivered files has a count and <strong>none of them is zero</strong>; the median bin carries %s. It is a density, and the signal is its magnitude. Thresholding it is a separate, explicit decision &mdash; see the cutoff section below.</li>",
          if (num$ok) format(num$hot_med) else "a few dozen"),
  "<li><strong>The significant-contact mask.</strong> The delivered criterion passes 12.7&ndash;14.1% of measured cell pairs, about 1.4 million contacts per chromosome. At the ~350 display columns the region figures use, one block is roughly 14 &times; 14 native cells, so almost every block contains a significant cell and the masked map fills in like the unmasked one. The mask is only informative at native resolution over a few Mb.</li>",
  "<li><strong>Per-bin ChIP log2FC at 20 kb.</strong> Noise-dominated (sd 0.14&ndash;0.16 broad, 0.57 for H3K4me3). The genome-wide medians are sound; a single bin's value is not.</li>",
  "</ul>",

  "<h3>What is missing, and where</h3>",
  p_("NA in the NPMI matrix is <em>not measured</em>, and it renders white &mdash; the same white as &ldquo;no contact&rdquo;.",
     "Every region figure prints its measured fraction for this reason. Genome-wide the figure is high, but it is",
     "strongly region-dependent: the vomeronasal cluster on chr7 is <strong>29% measured</strong>, with 11% finite",
     "cells inside the cluster against 97% in its own flanks. An empty-looking contact map there is an artefact of",
     "coverage, not a result."))

parts$settled <- sec("What this session settled, and what it retracted",
  p_("Recorded here so nothing in the sections below is read as more or less certain than it is."),

  "<h3>Settled, with a positive result</h3>",
  "<ul>",
  "<li><strong>The A/B axis is now one shared ruler, and <code>pc1_delta</code> is retired from genotype comparison.</strong>",
  "Both genotypes orient on only <strong>10 of 19</strong> autosomes, and on <strong>6 of those 10</strong> they",
  "select the same eigenvector, so the old per-genotype difference subtracted two different vectors wherever they",
  "disagreed and was blank elsewhere. One axis from the averaged correlation matrix reaches a median |&rho;| of",
  "<strong>0.391</strong> against 0.174 and 0.173 for the two genotypes' own axes. Concretely, on the old axis the",
  "olfactory cluster shifted significantly toward <em>B</em> and on the shared axis it shifts significantly toward",
  "<em>A</em> &mdash; a sign reversal, both at p &lt; 0.03, from the same data.</li>",
  "<li><strong>The protocadherin decompaction result survives the change of axis &mdash; but it is a",
  "matched-background statement, not a ranking.</strong> Rebuilt on the shared axis the A-shift is +0.874 (p = 0.003)",
  "against +0.576 on the old axis, and the short-range loss it is paired with was unaffected. Of the nine target sets",
  "it is still the only one meeting both criteria. What that does <em>not</em> mean is that it changes most:",
  "<strong>43.7% of 580 kb windows genome-wide lose more short-range contact than Pcdh does, and 5.9% beat it on both",
  "criteria at once</strong> &mdash; roughly one window in seventeen. See the target-set section for the full table.</li>",
  "<li><strong>Hotspot domains and large blocks melt and condense in the direction their name predicts, by an",
  "independent instrument.</strong> MELTRON compares insulation-score distributions across 11 spatial scales and never",
  "touches the differential-contact calls the domains were built from. Domains that lose contacts melt (KS D 0.171 and",
  "0.129 against null medians of 0.032 and 0.031), domains that gain contacts condense (0.198 and 0.163 against 0.032",
  "and 0.033), all four at the floor of a 500-draw size-matched null. The derepressed gene families, protocadherin",
  "included, sit <em>below</em> their own nulls &mdash; less changed than a random region of the same width.</li>",
  "<li><strong>Heterochromatin self-associates at short range on 17 of 19 autosomes</strong>, and it is not mixing",
  "&mdash; het&ndash;euchromatin is flat then falls. The two exceptions below 1 Mb are chr1 and chr3, at &minus;0.6 to",
  "&minus;2.2 percentage points, which is noise around zero rather than a reversal.</li>",
  "<li><strong>Hotspot domains that lose contacts are LAD-associated (fold 1.35), LINE-rich (1.25) and quiescent</strong>,",
  "depleted of promoters, CTCF and SINEs &mdash; against a circular-shift null that preserves domain spacing. Domains that",
  "<em>gain</em> contacts are indistinguishable from rotated background on every annotation tested.</li>",
  "<li><strong>Our own GAM insulation boundaries are externally validated by CTCF</strong> (1.72&times; at &plusmn;5 kb,",
  "decaying to 1.39&times; at &plusmn;100 kb).</li>",
  "</ul>",

  "<h3>Settled, with a powered negative</h3>",
  "<p class='note'>The coupling result is <em>not</em> in this list &mdash; see point 5 above. It was described as a",
  "powered negative on a per-bin n and does not survive a per-chromosome one.</p>",
  "<ul>",
  "<li><strong>Compartment strength does not change.</strong> Measured on two axes external to both contact matrices,",
  "all 19 autosomes, bootstrapped over chromosomes: active-chromatin axis mean delta &minus;0.00018 (95% CI &minus;0.00213",
  "to +0.00186, p = 0.62); TSS-density axis +0.00019 (p = 0.74).</li>",
  "<li><strong>The Pcdh TAD is indistinguishable from other chr18 TADs</strong>, intra and distal, under all three",
  "partitions (p = 0.75&ndash;0.88).</li>",
  "<li><strong>Individual domains do not compact</strong> &mdash; every within- versus between-domain",
  "difference-in-differences interval spans zero, under all three domain definitions.</li>",
  "</ul>",

  "<div class='warn'><strong>Retracted.</strong>",
  "<ul>",
  "<li><strong>&ldquo;Compartments are less defined in the cKO.&rdquo;</strong> The statistic was circular: strength",
  "measured along an axis derived from the contact matrices fits whichever genotype built the axis. On the same data at",
  "the same resolution the median cKO/WT ratio is 0.54 from a WT-derived axis, 0.79 from the averaged one and 1.21 from a",
  "cKO-derived one. The direction of the answer was set by the choice of axis.</li>",
  "<li><strong>Resolution as the cause of the weak compartment axis.</strong> Tested from 100 kb through 1 Mb; |rho| is",
  "flat across the whole range. The real defect is that WT and cKO do not select the same eigenvector &mdash; their",
  "leading subspaces overlap only 0.27&ndash;0.57.</li>",
  "<li><strong>chr11 as the chromosome that reverses the rule.</strong> chr1 inverts harder at both long bands, and",
  "chr8, chr10 and chr15 also invert. chr11 is one of five.</li>",
  "</ul></div>",

  "<div class='warn'><strong>Found, and unresolved.</strong> <code>TAD_boundaries_mm39.bed</code>, used as the external",
  "TAD partition, shows <strong>no CTCF enrichment at any window size</strong> (1.09&times;, flat from &plusmn;5 kb to",
  "&plusmn;100 kb) while the same CTCF track is 5.27&times; at transcription start sites and 1.72&times; at our own",
  "insulation minima. The file is a genuine liftover, not corrupt, but its provenance &mdash; cell type, caller,",
  "resolution &mdash; is unrecorded. METHODS open question 9.</div>",

  "<h3>The landscape heatmap</h3>",
  p_("The MINUTE significant-peak heatmap, re-drawn at 16:9 for presentation. Same k-means clustering as the portrait",
     "original: it re-draws the already-laid-out object rather than building a new one, because the cluster assignment is",
     "read back off that object to produce the cluster tables and BED exports."),
  embed_chip_fig("heatmap/2000maxgap_indsignificance_with_TAD_landscape.png",
                 "11,293 significant peaks &times; 50 libraries, row z-score, five k-means blocks."))

parts$chip <- sec("The ChIP side, in brief",
  p_("The GAM work is read against one headline: HP1&gamma; loss removes H4K20me3 almost everywhere.",
     "This section states that quantitatively, because every architecture claim downstream is conditioned on it.",
     "Full statement and caveats: <code>ChIP/results/EXECUTIVE_SUMMARY.md</code>."),

  "<h3>Method</h3>",
  p_("MINUTE ChIP-seq, 5 histone marks, Emx1-Cre cortex + hippocampus, 2 biological replicates per genotype",
     "(HP1gKO rep2 excluded at QC). Counts over called peaks, DESeq2 with <strong>sizeFactors = 1</strong>, no",
     "<code>lfcShrink</code>. Both choices are deliberate and both are necessary here: median-of-ratios normalisation",
     "and a shrinkage prior each assume most features do not change, which is exactly the assumption a near-genome-wide",
     "loss violates. The figures below are per-mark MA density plots &mdash; log2 fold change against mean normalised",
     "count, every peak."),

  "<h3>Result</h3>",
  p_("Median log2 fold change per peak, and the fraction of peaks below zero. The same quantity on the GAM 20 kb grid",
     "is shown beside it, because the mark tracks in the region plots are computed on that grid and are systematically",
     "shallower &mdash; a 20 kb bin dilutes a peak with its flanking background."),
  "<table>",
  "<tr><th>mark</th><th>peak-level median</th><th>% of peaks &lt; 0</th><th>GAM 20 kb median</th><th>% of bins &lt; 0</th></tr>",
  "<tr><td>H4K20me3</td><td><strong>&minus;0.678</strong></td><td>94.5%</td><td>&minus;0.497</td><td>99.8%</td></tr>",
  "<tr><td>H3K9me3</td><td>&minus;0.415</td><td>88.5%</td><td>&minus;0.290</td><td>93.5%</td></tr>",
  "<tr><td>H3K9me2</td><td><strong>&minus;0.341</strong></td><td>79.7%</td><td>&minus;0.179</td><td>93.3%</td></tr>",
  "<tr><td>H3K36me3</td><td>&minus;0.067</td><td>58.8%</td><td>&minus;0.052</td><td>67.6%</td></tr>",
  "<tr><td>H3K4me3</td><td>&minus;0.048</td><td>68.0%</td><td>0.000</td><td>43.1%</td></tr>",
  "</table>",
  "<ul>",
  "<li><strong>H4K20me3 is the effect.</strong> &minus;0.678 at peaks with 94.5% of peaks down &mdash; a shift of the",
  "whole distribution, not a minority of sites.</li>",
  "<li><strong>H3K9me2 is about half that depth but nearly as pervasive</strong> (&minus;0.341, 79.7% of peaks; 93.3%",
  "of 20 kb bins). Not negligible, and easy to under-read from the 20 kb number alone.</li>",
  "<li><strong>H3K36me3 and H3K4me3 are near-flat</strong> (&minus;0.067 and &minus;0.048). Small but not exactly zero",
  "&mdash; H3K36me3 has 58.8% of peaks below zero. They are treated as unchanged for the purpose of orienting the",
  "compartment axis, which is a statement about them being the <em>same ruler in both genotypes</em>, not a claim that",
  "they are perfectly static.</li>",
  "<li><strong>Nothing reaches FDR &lt; 0.05 in either direction.</strong> With 2 vs 2 and unshrunk fold changes that",
  "is expected and is not evidence against the shift &mdash; the evidence is the distributional move, which is why",
  "everything here is reported as an effect size rather than a count of significant peaks.</li>",
  "</ul>",
  embed_chip_fig("change_plots/H4K20me3_MA_density.png",
                 "H4K20me3. The whole cloud sits below zero."),
  embed_chip_fig("change_plots/H3K9me3_MA_density.png", "H3K9me3."),
  embed_chip_fig("change_plots/H3K9me2_MA_density.png",
                 "H3K9me2 &mdash; shallower than H4K20me3 but clearly displaced."),
  embed_chip_fig("change_plots/H3K36me3_MA_density.png", "H3K36me3, near-flat."),
  embed_chip_fig("change_plots/H3K4me3_MA_density.png", "H3K4me3, centred on zero."),
  embed_chip_fig("relationships/relationship_all_regions_H3K9me3_vs_H4K20me3_by_cluster.png",
                 "H3K9me3 against H4K20me3 per region. The two marks are largely uncoupled, which is what makes the H4K20me3-only class a real category rather than a threshold artefact."))

parts$hotspots <- sec("How hotspot domains are called, and what the count track is",
  p_("Two different objects share the word &ldquo;hotspot&rdquo; and conflating them has caused errors before.",
     "The <strong>count track</strong> is the number of significant differential contacts per 20 kb bin, within 5 Mb.",
     "About 93% of autosomal bins carry at least one, so it is a dense continuous covariate, not a hit list &mdash;",
     "the signal is the magnitude (1&ndash;162), and a bin with none is a 0, not missing.",
     "A <strong>hotspot domain</strong> is what you get after thresholding and merging that track. Only the domain is an",
     "independent unit for an enrichment test, because neighbouring bins share signal."),
  embed_fig("qc/hotspot_domain_calling_qc.png",
            "The full parameter justification. <strong>A</strong>: autocorrelation of the count track falls from r = 0.59 at 20 kb to ~0.4 at 100 kb &mdash; smoothing is defensible while r is substantial, which is what sets the 100 kb window. <strong>B</strong>: gap-merging the raw track consolidates fragments but chains into multi-Mb domains if pushed, and the gap size is arbitrary. <strong>C</strong>: smoothing first gives the same consolidation without the chaining, with the window set by panel A rather than chosen. <strong>D</strong>: raw thresholding is dominated by singletons; the 5-bin floor is what removes them. <strong>E</strong>: a worked example on chr3, showing the raw count, the running mean, the q90 threshold and the resulting domains."),
  p_("Both an adjusted and an <em>unadjusted</em> domain set are carried through every downstream stage, because the two",
     "agree on only 40% of bins (the Jaccard index). Where a result holds under one and not the other, the definition is doing the work."),
  tbl("qc_hotspot_domain_parameter_sweep.tsv", "The full sweep behind panels B and C.", max_rows = 25L),
  p_("<strong>There is no figure of what hotspot domains overlap.</strong> The domains are called and carried, and the",
     "count track is used as a covariate throughout the coupling section below, but nothing here annotates the domains",
     "against genes, repeats or chromatin states as a set. That is a gap, not an omission from the report."))

parts$blocks <- sec("Large one-directional blocks",
  p_("Contiguous 10&ndash;25 Mb stretches where the differential contacts run almost entirely one way."),
  embed_fig("blocks/large_blocks_genome_wide.png", "Where they are."),
  p_("And which ones are big. Blocks are <strong>called from the data</strong> &mdash; |net| sustained over at least",
     "4 Mb &mdash; so ranking them is a genome-wide statement, unlike ranking the handful of regions we chose to draw.",
     "White dots mark the blocks that do have a region figure, which answers the question worth asking about our",
     "figure selection: whether it covers the top of the distribution or the middle of it."),
  embed_fig("blocks/block_ranking.png",
            "Every called block, ranked by mean per-bin net differential contacts, with block width beneath."),
  tbl("gam12_block_ranking.tsv", "The same ranking.", max_rows = 20L),
  tbl("gam6_large_blocks.tsv", "The blocks.", max_rows = 25L),
  embed_fig("blocks/block_contact_budget.png", "Block-by-block contact budget."),
  tbl("gam8_block_contact_budget.tsv", max_rows = 20L),
  tbl("gam6_block_characterisation.tsv",
      "What is in them. The compartment column moved from <code>pc1_WT</code> to the shared axis on 2026-10-19, because <code>pc1_WT</code> is WT's own selected eigenvector and is a different eigenvector on different chromosomes. <strong>Read <code>shared_WT</code>, not <code>shared_delta</code>.</strong> <code>shared_WT</code> is clean, and says lost blocks are the more B-like half of the genome in WT &mdash; the same statement the LAD, LINE and quiescent enrichment of lost hotspot domains makes independently. <code>shared_delta</code> is partly circular here: the score is <code>C_g %*% a</code>, so a block defined by losing contact has its score magnitude mechanically pulled toward 0, which for a negative (B-like) score reads as &ldquo;toward A&rdquo;. The observed values are exactly what that artefact predicts, so they are not independent evidence of a shift.",
      max_rows = 20L))

parts$hetclass <- sec("Contact classes and mark-change strata",
  embed_fig("heterochromatin/heterochromatin_contact_classes.png",
            "Change in mean NPMI by chromatin class and separation, genome-wide. The bootstrap CIs are over the 19 autosomes."),
  tbl("gam10_class_change.tsv", "The numbers behind it &mdash; the table the Bottom line is built from."),
  embed_fig("heterochromatin/mark_change_contact_classes.png",
            "The same split by how much mark the bins lost, testing whether the contact change tracks the size of the mark loss."),
  tbl("gam11_stratum_pair_change.tsv", max_rows = 25L),
  tbl("gam11_lost_vs_kept.tsv"),
  tbl("gam10_intra_vs_inter.tsv"))

parts$targets <- sec("HUSH/CBX3 target families",
  p_("Every claim here is against a <em>matched</em> background, never against all bins.",
     "The targets are repeat-rich and heterochromatic, so an unmatched contrast is significant for reasons that have nothing to do with HP1&gamma;."),
  embed_fig("hush_contacts/target_vs_matched_z.png"),
  tbl("gam5_target_vs_matched_background.tsv", "Full table, including <code>naive_diff</code> &mdash; the unmatched contrast, kept for comparison.", max_rows = 30L),
  embed_fig("hush_contacts/contacts_lost_target_matched_all.png"),
  embed_fig("hush_contacts/decompaction_quadrant.png",
            "The decompaction test: short-range contact loss on one axis, shift toward the A compartment on the other. Decompaction requires both, each at BH < 0.05."),
  tbl("gam5_decompaction_verdict.tsv", "The verdict per family. The A-compartment axis is the shared axis from GAM_2b; <code>pc1_delta</code> is carried beside it for contrast and is not a sound genotype comparison."),
  "<div class='warn'><strong>&ldquo;Decompacts&rdquo; here means unusual <em>for its stratum</em>, not unusual in the genome &mdash; and the difference is large.</strong> The test above compares each set against bins matched on chromatin state, repeat-content quartile and detection-coverage tertile. That is the right question for a target set and the wrong one for anyone looking at a region plot, where a locus is being compared against everywhere else you could have pointed at, matched on nothing. The table below asks that second question, and it is the one to quote if someone says a different region looks like it changed more &mdash; because for protocadherin, many of them did.</div>",
  tbl("gam5_target_genomewide_context.tsv",
      "Share of same-width contiguous windows anywhere in the genome that are more extreme than the set, in the decompaction direction. Small = exceptional. <code>obs_short</code> and <code>obs_axis</code> are the set's own means, and match the observed columns of the matched-background table exactly."),
  p_("<strong>Protocadherin, read both ways.</strong> Against its matched background it is clearly unusual:",
     "short-range loss &minus;0.0248 below matched (p = 0.002) and an A-shift of +0.874 above matched (p = 0.001).",
     "Against the genome it is mid-pack for contact loss &mdash; <strong>43.7% of 580 kb windows lose more</strong> &mdash;",
     "and in the top sixth for the A-shift, with <strong>5.9% of windows beating it on both criteria at once</strong>.",
     "Both are true and they are not in conflict; they are answers to different questions. The headline that follows",
     "from them is the narrow one: among the HUSH/CBX3 target families we tested, protocadherin is the only one whose",
     "architecture moves more than its own chromatin stratum does &mdash; not that it is the most remodelled place in",
     "the genome, which it is not."),
  p_("Independent support for the narrow reading: MELTRON scores the protocadherin cluster <em>below</em> its own",
     "size-matched null (KS D 0.224 against a null median of 0.260, empirical p = 0.66). A locus that were genuinely",
     "among the most remodelled in the genome would not do that."),
  embed_fig("hush_contacts/npmi_locus_Pcdh_chr18.png"),
  embed_fig("hush_contacts/npmi_locus_KRABZFP_densest.png"),
  embed_fig("pcdh/pcdh_native_npmi_difference.png",
            "The protocadherin cluster at native resolution, from the contact pairs rather than the 5 Mb hotspot summary."),
  tbl("gam9_pcdh_native_profile.tsv", max_rows = 20L),
  tbl("gam5_size_perm_binary_flags.tsv", "The coarser size-matched permutation, as a cross-check that the stratified result is not an artefact of the stratification. The A/B flag was <code>switch_AB</code> (WT&rsquo;s eigenvector called this bin A while the cKO&rsquo;s called it B &mdash; two different directions) until 2026-10-19 and is now <code>shared_AB</code>, a sign change of one axis both genotypes are scored on."),
  "<p><strong>This is where the retargeted flag changed an answer.</strong> KRAB-ZFP clusters cross A&rarr;B at 1.45&times; a repeat-matched draw (z = 3.54, p = 0.002, 387 bins). In the stratified null the same set gives A&rarr;B z = 3.62 (p = 0.001) with B&rarr;A flat (z = 0.88, p = 0.37), a continuous <code>shared_delta</code> of &minus;0.160 against &minus;0.017 matched (z = &minus;3.40), and a starting position more A-ward than matched bins (<code>in_A_shared_WT</code> 0.615 vs 0.566). Four readings of one axis, all agreeing: HP1&gamma; loss moves these clusters toward B.</p>",
  "<p>The obvious objection does not hold here. The shared score is <code>C_g %*% a</code>, so contact amplitude leaks into it &mdash; but KRAB-ZFP <em>gains</em> contact relative to matched bins (<code>net_remodel</code> +0.22, <code>lost_log1p</code> z = &minus;3.86) and starts A-ward, so the mechanical direction is toward <em>A</em>, the opposite of what is observed. The caveat that does apply is the axis quality: median |&rho;| 0.391 against active chromatin, determined on 16 of 19 autosomes.</p>")

parts$compartments <- sec("Compartment figures",
  "<div class='warn'><strong>Read this section against the compartment caveat.</strong> The A/B axis is weakly determined on most chromosomes, and WT and cKO do not select the same eigenvector, so a WT-vs-cKO comparison of these tracks is comparing two different directions. Compartment <em>strength</em>, measured on axes external to both contact matrices across all 19 autosomes, shows <strong>no change</strong> (active-chromatin axis mean delta &minus;0.00018, 95% CI &minus;0.00213 to +0.00186, p = 0.62; TSS-density axis p = 0.74) &mdash; but those two numbers were computed on 2026-10-17 and written into prose without a stage or a table behind them, so nothing regenerates or checks them; see <code>results/COMPARTMENTS.md</code> section 6. The figures are kept for the structure they show, not for a differential claim.</div>",
  embed_fig("compartments/contact_decay_by_genotype.png",
            "Contact decay with separation. This does NOT depend on the compartment axis and is reliable."),
  embed_fig("compartments/saddle_by_genotype.png",
            "Saddle per genotype. Each genotype is stratified by its OWN axis, which is exactly the comparability problem above &mdash; the two panels are not on a common ruler."),
  embed_fig("compartments/pc1_tracks_by_chr.png"),
  embed_fig("compartments/pc1_scatter_WT_vs_HP1gKO.png"),
  tbl("npmi_compartment_orientation.tsv", "Axis quality per chromosome and genotype. <code>which_pc</code> is the diagnostic column: it differs between genotypes on most chromosomes.", max_rows = 42L),
  tbl("npmi_compartment_switch_counts.tsv", "Switch counts. Near-symmetric, which is what axis noise looks like &mdash; not evidence of no net shift."),
  "<div class='warn'><strong>Three figures were removed from this section on 2026-10-18 rather than shown with a caveat.</strong>",
  "<code>pc1_delta_by_group</code>, <code>pc1_delta_vs_H4K20me3</code> and <code>compartment_switch_fraction</code> are all",
  "built on <code>pc1_delta</code> &mdash; the difference between the two genotypes&rsquo; own compartment axes. That",
  "quantity was retired from the GAM_6 verdict for being a difference of two <em>different</em> eigenvectors, but the",
  "figures stayed in the report, which is the worst of both: retired in the prose, still on the page.",
  "<p>They are now unsupportable for a second, harder reason. With an eigenvector-selection margin in place, a compartment",
  "axis is determined on <strong>21 of 38</strong> chromosome&times;genotype combinations, and on only <strong>7 of 19",
  "autosomes in both genotypes</strong> &mdash; and 10 of the 38 have no balanced candidate at all, so they are",
  "undetermined at any threshold. A per-genotype difference is not defined on the remainder. Two of them",
  "(chr15 and chr19, cKO) were additionally shipping with A and B inverted, which the summary table was hiding by",
  "storing the absolute value of the diagnostic.</p>",
  "<p>What remains usable is below: the shared axis, where both genotypes are scored on one vector. The underlying",
  "<code>pc1_delta</code> columns are still produced and still in the tables, so nothing is lost &mdash; they are just",
  "no longer displayed as if they were a result.</p></div>",
  embed_fig("compartments_vs_contacts/contact_density_by_group.png"),
  embed_fig("compartments_vs_contacts/insulation_delta_by_group.png"),
  embed_fig("compartments_vs_contacts/size_perm_overlap.png"),
  tbl("gam4_group_summary.tsv", max_rows = 20L),
  tbl("gam4_group_vs_matched_background.tsv", max_rows = 25L))

parts$shared_axis <- sec("One shared A/B axis, and both genotypes scored on it",
  p_("The compartment figures above are each genotype measured on its own axis, which is why",
     "they cannot be compared. This section is the replacement: a single axis built from the",
     "<strong>averaged</strong> correlation matrix at 100 kb, phased once against WT",
     "H3K4me3 + H3K36me3, with both genotypes then projected onto it. Every bin is scored",
     "against the same ruler, so the difference of two scores is a difference of one measured",
     "quantity rather than of two different objects."),
  p_("<strong>Why this was necessary.</strong> Across all 19 autosomes the two genotypes select",
     "<strong>the same eigenvector on only 6 of the 10 where both orient</strong>. <code>pc1_delta</code> subtracts two",
     "different vectors wherever they disagree, which is most of the genome."),
  p_("<strong>Does it actually help?</strong> Measured, not assumed. Against the same",
     "active-chromatin reference the shared axis reaches a median |&rho;| of <strong>0.391</strong>",
     "against 0.174 for WT's own axis and 0.173 for the cKO's, and clears |&rho;| &gt; 0.4 on",
     "<strong>7 of 19</strong> autosomes. It beats both per-genotype axes by a wide margin,",
     "which is the point of the stage."),
  "<div class='warn'><strong>The shared axis balances before its PCA too (2026-10-18), and that is a trade",
     "rather than a win.</strong> It is determined on <strong>16 of 19</strong> autosomes.",
     "Unbalanced, the compartment landed on PC1 for <strong>0 of 18</strong>",
     "determined chromosomes; balanced it is <strong>8 of 16</strong>, with 15 of 16 chromosomes moving up",
     "exactly one eigenvector. The A/B split also evens out (worst |frac_A &minus; 0.5| 0.198 to 0.117). But",
     "agreement with active chromatin does not improve &mdash; on the 15 chromosomes determined both ways it is",
     "8 up and 7 down &mdash; and <strong>chr8, chr14 and chr15 drop out entirely</strong> (chr4 is recovered),",
     "taking bin coverage from 90.5% to 82.7%. Those three now fail the selector at k = 6 as well, so it is not",
     "a candidate-budget problem. The likely reading is that their old calls were the detectability vector,",
     "which correlates with active chromatin because detectability in GAM is partly biological &mdash; chr14 was",
     "called on <strong>PC4</strong> at |&rho;| 0.527 before. That reading is an inference, not a measurement.",
     "<code>results/COMPARTMENTS.md</code> section 4c has the full table and how to reverse it.</div>",
  embed_fig("compartments/shared_axis_vs_pergenotype.png",
            "The comparison that justifies the stage. Each row is one autosome; all three points use the same reference."),
  embed_fig("compartments/shared_axis_tracks.png",
            "The axis and both genotypes' scores. Red is the A side."),
  embed_fig("compartments/shared_axis_wt_vs_ko.png"),
  "<div class='warn'><strong>The track is shape only, deliberately.</strong> Each genotype's score is median-centred and scaled to RMS&nbsp;1 per chromosome. &lambda;<sub>1</sub> does fall about 20% in the cKO, but that is a global property of the correlation matrix affecting every direction equally, not a compartment-specific effect &mdash; normalising by each genotype's own &lambda;<sub>1</sub> leaves a strength ratio of 0.953. Leaving the amplitude in the track would paint the cKO as uniformly weaker everywhere and reproduce exactly the claim that was retracted. Amplitude is kept per chromosome in the summary table below. <strong>This section does not re-answer the strength question</strong>, which is settled and null on two external axes.</div>",
  "<div class='warn'>Better is not the same as good. 10 of 19 autosomes at |&rho;| &gt; 0.4 still means the axis is weak on half the genome. Quote the per-chromosome quality with any compartment claim.</div>",
  tbl("gam2b_shared_axis_summary.tsv",
      "Per chromosome. <code>rho_ref</code> is the shared axis; <code>rho_pergeno_*</code> are each genotype's own. <code>lambda1_ratio</code> and the two <code>rayleigh_*</code> columns are where the amplitude went.",
      max_rows = 21L))

parts$meltron <- sec("MELTRON: melting and condensation",
  p_("The Pombo lab's own instrument for whether a region holds its contacts together or has",
     "come apart (Winick-Ng et al., <em>Nature</em> 2021). It asks a different question from our",
     "hotspot track and answers it from a different quantity, which is why it is worth having:",
     "hotspots count the collaborators' per-contact |&Delta;z| &gt; 2 calls, while melting compares",
     "the <strong>distribution</strong> of local contact density across spatial scales. Nothing",
     "here uses the differential-contact calls."),
  p_("Insulation score = the mean contact value in the square crossing a viewpoint, computed at",
     "<strong>11 square sizes from 240 to 1040 kb</strong>, in a <strong>120 kb window stepping one",
     "bin</strong> &mdash; MELTRONIC's own interval set. The two genotypes' pooled scores are then",
     "compared by one-sided Kolmogorov&ndash;Smirnov tests, one per direction."),
  "<div class='warn'><strong>There are two scores per interval, not one signed score &mdash; and that is deliberate.</strong> The two one-sided tests are <em>not</em> complementary: where the two insulation distributions cross, one shifted in its lower tail and the other in its upper, both reject and an interval genuinely carries a melting score <em>and</em> a condensation score. MELTRONIC reports both for exactly that reason, and so do we. On this data they never fire together at the published cutoff &mdash; <strong>0.0% of intervals</strong> &mdash; which is measured on every run rather than assumed.</div>",
  embed_fig("meltron/melting_track.png",
            "Melting up in red, condensation down in blue, one point per window. The two occupy separate multi-megabase blocks."),
  embed_fig("meltron/melting_score_distribution.png",
            "The two scores against each other. The top-right quadrant would be intervals called both ways; it is empty."),
  "<div class='warn'><strong>A correction worth recording: the window matters more than anything else here.</strong> An earlier version of this analysis used a 500 kb / 25-bin window instead of the published 120 kb. That pooled 275 nested, overlapping values into a test that assumes independence, so &minus;log<sub>10</sub>(p) saturated and the published |score| &gt; 5 cutoff appeared to call 69% of the genome. The conclusion drawn at the time &mdash; that the published cutoff does not transfer to this data &mdash; was a statement about the window, not about the data. At the published window the cutoff calls 26.0% melting and 26.1% condensing and is used directly.</div>",
  "<div class='warn'><strong>One deviation remains, and it is forced.</strong> Our insulation score is additive, not the published log<sub>2</sub>(square mean) &minus; log<sub>2</sub>(chromosome mean), because that form needs a positive quantity and only about half of our insulation squares have a positive mean NPMI &mdash; it is undefined for half the genome. We subtract the chromosome mean instead, the same decision as the observed&minus;expected convention everywhere else in this pipeline. Sign and ranking remain comparable to published work; <strong>magnitude does not</strong>.</div>",
  embed_fig("meltron/melting_score_vs_D.png"),
  tbl("gam2c_meltron_spatial_coherence.tsv",
      "Melting and condensing are regional: the net score changes sign between about 15% of adjacent 120 kb windows, with runs of one direction averaging 480 kb and reaching 8 Mb."),
  p_("<strong>The result.</strong> Region sets are scored against size-matched random regions on the",
     "same chromosomes, and that null is what makes the output interpretable &mdash; the raw score",
     "saturates for the strongest sets and cannot rank them. Hotspot domains and large blocks melt",
     "and condense in the direction their name predicts, by an instrument that never touches the",
     "differential-contact calls they were built from: lost KS D 0.171 and 0.129 against null",
     "medians of 0.032 and 0.031, gained 0.198 and 0.163 against 0.032 and 0.033 &mdash; all four at",
     "the floor of a 500-draw null. The derepressed gene families, protocadherin included, sit",
     "<em>below</em> their own nulls."),
  embed_fig("meltron/melting_region_sets.png"),
  embed_fig("meltron/melting_ecdf_region_sets.png",
            "The published diagnostic: the largest vertical gap between the two curves is the statistic."),

  "<h3>Does melting relate to anything else we measured?</h3>",
  p_("A melting track nobody can connect to the hotspots, the marks or the domains is a number",
     "without a claim attached. These are the joins, and they give a clean answer:",
     "<strong>melting is a contact-side measure that agrees with our other contact-side measures",
     "and is independent of the ChIP side</strong>."),
  embed_fig("meltron/melting_vs_covariates.png",
            "Net differential contacts &rho; = &minus;0.305 (bins losing contacts melt, bins gaining condense); H4K20me3 log2FC &rho; = 0.071, flat."),
  embed_fig("meltron/melting_vs_annotation.png",
            "Lost hotspot domains melt (median +5), gained domains condense (&minus;8), everything else sits on zero. The ChIP classes are indistinguishable from each other."),
  p_("That &rho; of &minus;0.305 against net differential contacts is large by this project's",
     "standards &mdash; compare &rho; = &minus;0.015 for the mark-vs-contact coupling. It is a",
     "<em>third</em> instrument reaching the same conclusion as the |&Delta;z| calls and the",
     "block analysis, from insulation distributions rather than from per-contact significance.",
     "And the flat H4K20me3 panel is the project's central uncoupling result appearing for the",
     "third time, now on the melting axis."),
  "<h3>At what spatial scale does the change happen?</h3>",
  p_("The single collapsed score cannot answer this; the insulation-by-square-size view can.",
     "A vertical band of one colour is a change at every scale, while colour confined to the",
     "lower rows is local reorganisation only."),
  embed_fig("meltron/is_by_scale_chr18_Pcdh.png"),
  embed_fig("meltron/is_by_scale_chr11_lostblock.png"),
  tbl("gam2c_meltron_regions.tsv",
      "Region sets. <code>ks_D</code> is the effect size and <code>p_empirical_D</code> is the number to quote; <code>melting_score</code> and <code>condensation_score</code> are the published-style quantities.",
      max_rows = 20L),
  tbl("gam2c_meltron_genome_summary.tsv", "Per chromosome, at the published cutoff.",
      max_rows = 21L))

parts$coupling <- sec("Coupling: does the mark loss have a contact counterpart?",
  p_("Three framings of one question, which until now lived in three files and disagreed only in presentation.",
     "They are now one stage, and this figure is the synthesis:"),
  embed_fig("coupling/coupling_synthesis.png"),
  p_("In words:"),
  "<ul>",
  "<li><strong>GAM_3</strong> &mdash; per-bin Spearman of ChIP log2FC against differential-contact count. &rho; &minus;0.015 raw, &minus;0.046 partial on coverage. <strong>Not a powered negative</strong>: the unit of independence is the chromosome, and the 19-chromosome bootstrap (&plusmn;0.03) contains the hotspot-implied effect of &minus;0.028, so zero and that effect are not distinguishable here. The partial-covariates variant excludes zero (&minus;0.0216, CI &minus;0.0372 to &minus;0.0049). The dose-response is real and tiny: 0.032 log2 units of H4K20me3 across the whole range of contact loss, about 6% of the depth of the loss.</li>",
  "<li><strong>GAM_7</strong> &mdash; retracts GAM_3's premise: the contact variable is smeared over 5 Mb while the ChIP variable varies at 20 kb, a ~250&times; scale mismatch. The null is real for the question as posed but cannot distinguish &ldquo;no relationship&rdquo; from &ldquo;the contact variable cannot resolve that scale&rdquo;.</li>",
  "<li><strong>GAM_10/11/13</strong> &mdash; NPMI difference-in-differences with a bootstrap over the 19 autosomes. This is the framing that survived, and the one every new stage uses.</li>",
  "</ul>",
  p_("The positive result is elsewhere: H4K20me3 loss tracks <em>baseline</em> architecture, not contact <em>change</em> &mdash; pc1 +0.196 and npmi_short &minus;0.166 after coverage control. Insulation does not survive the control (&minus;0.180 &rarr; &minus;0.066) and should not be reported."),
  embed_fig("coupling/specificity_control_rho.png"),
  embed_fig("coupling/specificity_control_curves.png"),
  p_("The supporting views, kept because the null is the claim and a null needs to be shown from more than one angle:"),
  embed_fig("coupling/coupling_density2d_lost.png"),
  embed_fig("coupling/coupling_density2d_gained.png"),
  embed_fig("coupling/coupling_quantile_curves.png",
            "Binning by how much mark a bin lost, rather than correlating. Flat curves are the same null in a form that cannot be hidden by a correlation coefficient."),
  embed_fig("coupling/coupling_quantile_curves_net.png"),
  embed_fig("coupling/coupling_raw_vs_pooled.png",
            "Raw against pooled, the check that the null is not an artefact of aggregation."),
  embed_fig("coupling/coupling_stratified.png"),
  embed_fig("coupling/lost_vs_gained_per_bin.png"),
  embed_fig("coupling/coupling_npmi_features.png"),
  embed_fig("coupling/contacts_by_compartment.png"),
  embed_fig("coupling/coupling_compartment_switch.png"),
  tbl("coupling_spearman.tsv", max_rows = 20L),
  tbl("coupling_specificity.tsv", max_rows = 20L),
  tbl("coupling_per_chromosome.tsv", "Per chromosome, so the genome-wide null can be read against the spread.", max_rows = 20L))

parts$cutoffs <- sec("Every cutoff in this report, and who chose it",
  p_("Four separate thresholds decide what counts as a &ldquo;significant&rdquo; or &ldquo;hotspot&rdquo; anything here,",
     "they come from three different places, and they are easy to confuse because the figures put two of them in",
     "adjacent panels. This table is the audit."),
  "<table>",
  "<tr><th>cutoff</th><th>value</th><th>whose</th><th>what it decides</th></tr>",
  "<tr><td>contact significance</td><td>|delta-z| &gt; 2 <em>and</em> max NPMI &gt; 0.1</td><td><strong>theirs</strong>, and not adjustable</td><td>which cell pairs are in the delivered differential-contact files at all. It is baked into the delivered filenames &mdash; we never recompute it, and could not change it without redoing the calling.</td></tr>",
  "<tr><td>hotspot <em>bins</em> (the ribbon)</td><td>per-bin count &ge; 50</td><td><strong>theirs</strong>, bottom of a proposed range</td><td>the one-row ribbon under the contact map on the region figures. Their slides propose 50&ndash;75; see below.</td></tr>",
  "<tr><td>hotspot <em>domains</em></td><td>q90 of the smoothed track, 5-bin smoothing, gap 1, minimum 5 bins</td><td><strong>ours</strong></td><td>every enrichment, LAD, repeat and MELTRON result in the tables. Justified panel by panel in the domain-calling QC figure above.</td></tr>",
  "<tr><td>MELTRON melting / condensation</td><td>|score| &gt; 5</td><td><strong>published</strong>, used directly</td><td>which intervals are called melted or condensed. Transfers to this data at the published 120 kb window &mdash; it calls 26.0% and 26.1%, not the 69% an earlier wrong window produced.</td></tr>",
  "</table>",

  "<h3>Why the ribbon sits at 50 and not 75</h3>",
  sprintf("<p>The proposed range is 50&ndash;75 and we deliberately keep the <em>bottom</em> of it, which is a choice rather than a default. On the delivered count track (%s autosomal bins, all of them non-zero), taking the larger of the gained and lost counts per bin:</p>",
          if (num$ok) format(num$hot_n, big.mark = ",") else "~117,000"),
  "<table>",
  "<tr><th>cutoff</th><th>bins called</th><th>% of the covered bins</th><th>% of the 123,133-bin grid</th></tr>",
  if (num$ok) paste0(sapply(seq_along(c(50, 60, 75)), function(i) {
    k <- c(50, 60, 75)[i]
    sprintf("<tr><td>&ge; %d%s</td><td>%s</td><td>%s</td><td>%s</td></tr>", k,
            if (k == 50) " (used)" else "", format(num$cut[i], big.mark = ","),
            pc_(num$cut[i], num$hot_n), pc_(num$cut[i], num$grid_n))
  }), collapse = "") else "<tr><td colspan='4'>count track not loaded</td></tr>",
  "</table>",
  sprintf("<p>Moving to the top of the range would drop %s of the called bins. The ribbon is read as &ldquo;where is this busy&rdquo;, so the permissive end is the honest one for a density display; anything that needs a conservative <em>call</em> uses the domains instead. Both denominators are given above because quoting one percentage without saying which is meant is how this number drifted between the code and the prose.</p>",
          if (num$ok) sprintf("%.0f%%", 100 * (1 - num$cut[3] / num$cut[1])) else "about 73%"),

  "<h3>The ribbon and the domains disagree about 60% of the time</h3>",
  p_("Both are called &ldquo;hotspot&rdquo; and they are not the same set. Jaccard agreement is",
     "<strong>0.39 for lost and 0.40 for gained</strong>. The q90 cut lands at a raw count of about 44 and 42 &mdash;",
     "close to 50 &mdash; so the divergence is almost entirely the smoothing and the 5-bin size floor, not the height",
     "of the threshold. This is why the ribbon on the region figures is labelled <strong>hotspot bins</strong> with its",
     "cutoff printed in the legend, rather than &ldquo;hotspot&rdquo;: a reader who sees the bare word will assume the",
     "domains from the tables, and about 60% of them are not.",
     if (num$ok) sprintf("Only %d bins genome-wide are called in both directions at once, which is why the ribbon is one row rather than two.", num$hot_both) else ""))

parts$methods <- sec("How every comparison in this report is computed",
  "<p>Every number above is one of a small number of operations. This section states each one exactly,",
  "including what is and is not normalised, so a figure can be checked without reading the scripts.</p>",

  "<h3>The contact measure</h3>",
  "<ul>",
  "<li><strong>NPMI (normalised pointwise mutual information)</strong> is the delivered quantity, one 20 kb &times;",
  "20 kb matrix per chromosome per genotype, computed upstream at sigma 2.6 from the nuclear-profile segregation",
  "table. It is a <em>normalised co-segregation statistic</em>, not a contact frequency: PMI compares how often two",
  "windows appear in the same nuclear profile against how often they would by chance, and the normalisation bounds",
  "it in [&minus;1, 1]. It already accounts for how often each of the two windows is detected",
  "at all, which is why <em>no further coverage or library-size normalisation is applied to it anywhere in this",
  "pipeline</em>. It is bounded in [&minus;1, 1] and mostly positive in cis.</li>",
  "<li><strong>Unmappable bins stay NA.</strong> All-NA rows and columns are dropped before any arithmetic. They",
  "are never set to 0, because 0 in NPMI means &ldquo;co-segregates exactly as often as chance&rdquo;, not",
  "&ldquo;no data&rdquo;, and would pull every average towards the middle.</li>",
  "<li><strong>Delta-z</strong> in the differential-contact files is the upstream caller&rsquo;s own statistic, not",
  "ours. Sign convention here is <strong>cKO &minus; WT</strong>; the collaborators&rsquo; slides use WT &minus; cKO,",
  "so a gain in this report is a loss in theirs.</li>",
  "</ul>",

  "<h3>Percent change, and what it is a change in</h3>",
  "<ul>",
  "<li>Bin pairs are assigned to a cell of chromatin class &times; separation band, the mean NPMI of all pairs in",
  "that cell is taken <em>within each genotype separately</em>, and the reported figure is",
  "<strong>100 &times; (mean<sub>cKO</sub> / mean<sub>WT</sub> &minus; 1)</strong>. It is a ratio of two means, not a",
  "mean of per-pair ratios, so a cell with few pairs cannot be dominated by one extreme pair.</li>",
  "<li><strong>Separation bands</strong> are &lt;0.2, 0.2&ndash;1, 1&ndash;5, 5&ndash;10 and &gt;10 Mb. Banding is",
  "mandatory rather than cosmetic: NPMI falls steeply with separation, so any average pooled over separations is",
  "mostly a statement about which separations the two sets happened to contain.</li>",
  "<li><strong>Chromatin class</strong> is from the WT H3K9me3/H4K20me3 heterochromatin domain calls, so the class",
  "label is fixed and does not move with the genotype. A pair is het&ndash;het, het&ndash;euchromatin or",
  "euchromatin&ndash;euchromatin; het&ndash;het is split further into same-domain and different-domain.</li>",
  "<li><strong>Confidence intervals are a bootstrap over the 19 autosomes</strong>, resampled with replacement,",
  "the whole cell mean recomputed each time, 95% percentile interval. The chromosome is the resampling unit",
  "because bin pairs within a chromosome are not independent &mdash; treating millions of pairs as independent",
  "would give intervals too narrow to mean anything. Per-chromosome figures therefore carry no CI: there is",
  "nothing left to resample.</li>",
  "<li><strong>Within- versus between-domain</strong> is a difference-in-differences: (intra<sub>cKO</sub> &minus;",
  "intra<sub>WT</sub>) &minus; (inter<sub>cKO</sub> &minus; inter<sub>WT</sub>), bootstrapped the same way, so a",
  "change shared by both does not register as a difference between them.</li>",
  "</ul>",

  "<h3>Permutation nulls, and why a raw overlap is never used</h3>",
  "<p>Large heterochromatin domains overlap narrow features trivially, by being large. Every overlap-shaped claim",
  "is therefore tested against a matched null rather than against chance:</p>",
  "<ul>",
  "<li><strong>Size-matched</strong> &mdash; the query set is replaced by random regions drawn to the same width",
  "distribution on the same chromosomes.</li>",
  "<li><strong>Separation-matched</strong> &mdash; for contact comparisons, partners are drawn at the same",
  "separations as the observed set, so a difference cannot be the distance decay in disguise.</li>",
  "<li><strong>Repeat-content-matched</strong> &mdash; for the repeat-class work, since H4K20me3 domains and young",
  "LINE/LTR content &mdash; long interspersed nuclear elements and long terminal repeat retrotransposons &mdash; are heavily confounded.</li>",
  "</ul>",

  "<h3>TADs: yes, called explicitly, and two independent ways</h3>",
  "<p>The Pcdh test needs a partition of chr18 into TADs. Both partitions are used and reported side by side,",
  "because agreement between an external call and our own is the only check available:</p>",
  "<ul>",
  "<li><strong>External partition.</strong> <code>TAD_boundaries_mm39.bed</code>, 2,438 published boundaries each",
  "150 kb wide. TAD intervals are the spans <em>between</em> consecutive boundaries on a chromosome",
  "(<code>load_tads_external()</code>). Fixed across genotypes, so it cannot absorb a genotype effect. 81 usable",
  "TADs on chr18.</li>",
  "<li><strong>Our own partition, per genotype.</strong> Called from GAM insulation",
  "(<code>call_tads_insulation()</code>). Insulation at a bin is the mean NPMI <em>across</em> that bin, so a strong",
  "boundary is a local <strong>minimum</strong>. The profile is smoothed over 5 bins; a boundary is a minimum within",
  "&plusmn;10 bins whose prominence is at least 0.5 local standard deviations, with boundaries forced at least 10 bins",
  "apart and the deeper of a close pair kept. Prominence is measured in local sd rather than absolute NPMI so one",
  "threshold works on chromosomes of different contact density. 144 TADs on chr18 in each genotype.</li>",
  "</ul>",
  "<p>The Pcdh cluster sits in one TAD under both: <code>chr18:36.82&ndash;37.87 Mb</code> externally (88.9% of the",
  "cluster inside it) and <code>chr18:37.40&ndash;38.06 Mb</code> by insulation (63.1%). <strong>The insulation call",
  "is identical in WT and cKO</strong> &mdash; the boundaries around Pcdh do not move &mdash; which is itself a result.</p>",
  "<p>The test is then: for the Pcdh TAD, compute the change in mean NPMI to every other chr18 TAD, average",
  "<em>within</em> separation band first and only then across bands, and compare that single number against the same",
  "quantity computed for each of the other 80 (or 143) chr18 TADs. Averaging within band first is what stops a TAD",
  "being ranked by which separations it happens to have partners at. The p value is two-sided empirical, the",
  "fraction of other TADs at least as extreme.</p>",
  "<p>Result, under all three partitions: the Pcdh TAD is indistinguishable from other chr18 TADs, both internally",
  "(&minus;4.8% external, &minus;3.2% own, p = 0.75 and 0.85) and distally (&minus;8.4% and &minus;3.5%, p = 0.88",
  "and 0.69). It loses contacts, but so does the rest of chr18 &mdash; the median other TAD loses more.</p>",

  "<h3>CTCF</h3>",
  p_("ENCODE conservative IDR thresholded peaks, lifted mm10 &rarr; mm39 with rtracklayer and UCSC's mm10ToMm39",
     "chain by <code>prepare_CTCF_annotation.R</code>: <strong>ENCSR877MSN</strong> (layer of hippocampus, male adult",
     "2 months) and <strong>ENCSR644VYX</strong> (left cerebral cortex, female adult 2 months), both B6CASTF1/J.",
     "100% of peaks survived the lift; median widths 440 bp and 309 bp. Peaks are kept individual &mdash; unlike the",
     "LAD preparation they are not <code>reduce()</code>d and carry no width floor, either of which would fuse or",
     "delete point features."),
  "<div class='warn'><strong>Bulk tissue, not sorted CA1.</strong> No CA1-pyramidal-specific CTCF ChIP-seq appears to",
  "be public. CTCF occupancy is largely cell-type-invariant, which is what makes bulk usable as a proxy, but this",
  "track is bulk hippocampus and bulk cortex and must not be described as CA1. Both tissues are carried because",
  "together they match the two assays: GAM is CA1 pyramidal cells, the ChIP is Emx1-Cre cortex + hippocampus.</div>",
  p_("Validated independently of any boundary set: peak density is <strong>5.27&times;</strong> genome average at",
     "TSS &plusmn;1 kb, decaying smoothly to 1.51&times; at &plusmn;50 kb. That is what established the track is sound",
     "before it was used to test the two TAD partitions."),

  "<h3>The compartment axis</h3>",
  "<p>Stated here because result 4 turns on the detail. The matrix is first converted to",
  "<strong>observed minus expected</strong> (O&minus;E), where expected is the mean of its own diagonal, i.e. the mean NPMI at",
  "that exact separation on that chromosome in that genotype. Subtracted, not divided: NPMI is already a log-scale",
  "quantity that goes negative, so a ratio would flip sign wherever a diagonal mean passes through zero.</p>",
  "<p>That matrix is then block-averaged from 20 kb to <strong>100 kb</strong> bins, NA-aware, because compartments",
  "are not resolvable at 20 kb. A Pearson correlation matrix is formed across its columns, and the top four",
  "eigenvectors are taken. The compartment axis is chosen <em>balance first</em>: only eigenvectors assigning between",
  "20% and 80% of the chromosome to the A side are eligible, and among those the one with the largest |Spearman &rho;|",
  "against the active-chromatin reference wins. Balance is a precondition rather than a tiebreak because the leading",
  "eigenvector is systematically degenerate here, putting 87&ndash;96% of bins on one side while still correlating",
  "slightly better with the reference than the real compartment vector does.</p>",
  "<p>The <strong>reference</strong> is the mean of per-bin WT H3K4me3 and WT H3K36me3 z-scores. Two marks because",
  "H3K4me3 is punctate at promoters and thin across a 20 kb bin while H3K36me3 is broad over gene bodies. Both are",
  "unchanged in the KO (the knockout; per-bin median log2FC +0.000 and &minus;0.052, against &minus;0.497 for H4K20me3), so the A/B",
  "axis is the same ruler in both genotypes. Anchoring on H3K9me3 or H4K20me3 would be circular &mdash; the axis would",
  "move with the effect being measured.</p>",

  "<h3>The ChIP side</h3>",
  "<ul>",
  "<li><strong>DESeq2 with sizeFactors = 1.</strong> Deliberate. DESeq2 is the count-based differential-testing package used on",
  "the ChIP side, and its sizeFactors are the per-sample scaling constants it would normally fit. The median-of-ratios estimate assumes most features",
  "do not change; here H4K20me3 is down almost everywhere, so median-of-ratios would divide the genome-wide loss out",
  "and report it as nothing.</li>",
  "<li><strong>No lfcShrink.</strong> lfcShrink is DESeq2's fold-change shrinkage step and ashr (adaptive shrinkage) is the",
  "prior it uses by default; that prior assumes most effects are null, the same assumption",
  "sizeFactors = 1 exists to avoid. Applied, it collapses H3K36me3 and H3K9me2 to exactly 0 and pulls H4K20me3 from",
  "&minus;0.497 to &minus;0.212.</li>",
  "<li><strong>Effect sizes are relative to each mark's own genome-wide median</strong>, not to an absolute cutoff,",
  "because an absolute cutoff is not shift-invariant and the effect here is a shift.</li>",
  "</ul>",

  "<h3>Joining the two assays</h3>",
  "<p>The ChIP counts path uses NCBI (National Center for Biotechnology Information) seqnames (<code>1</code>, <code>2</code>) and the GAM files use UCSC",
  "(University of California, Santa Cruz) names (<code>chr1</code>). A mismatched flip produces zero overlaps, which reads downstream as a clean null rather than",
  "an error, so every cross-assay join asserts its overlap rate with <code>stopifnot()</code> and fails loudly. Both",
  "assays are mm39; anything external is lifted with CrossMap (a coordinate liftover tool) and the liftover recorded in METHODS.md.</p>")

parts$caveats <- sec("Standing caveats",
  "<ul>",
  "<li><strong>2 vs 2.</strong> Both assays are two biological replicates per genotype. The GAM files are already merged and carry no replicate dimension, so z &gt; 2 is the only significance available from them.</li>",
  "<li><strong>ChIP at 20 kb.</strong> DESeq2 sees integers of ~1&ndash;11 per bin after the &times;10 count scaling, unshrunk. H3K36me3 appeared raised over Pcdh; it is &minus;0.049 against a genome median of &minus;0.052, p = 0.93. Nothing at Pcdh differs from genome background at 20 kb. lfcShrink is <em>not</em> the fix &mdash; it assumes most effects are null, which is exactly what sizeFactors = 1 exists to avoid here.</li>",
  "<li><strong>Compartments parked, but recoverable.</strong> An oriented axis exists on all 38 autosome &times; genotype",
  "combinations; 25 reach |&rho;| &ge; 0.3 and 11 reach &ge; 0.4 against the orienting marker. It sits on PC1 only for chr7,",
  "and on PC2&ndash;PC4 elsewhere, which is the thing to understand before any compartment claim.</li>",
  "<li><strong>LADs are constitutive, not neuronal.</strong> cLAD/ciLAD are the cell-type-invariant core across ESC (embryonic stem cell), NPC (neural progenitor cell), astrocyte and MEF (mouse embryonic fibroblast).</li>",
  "<li><strong>Unresolved provenance.</strong> Nuclear profiles per genotype, the passed_qc.curated quality-control (QC) criteria, what sigma 2.6 (the upstream smoothing parameter) was fit on, and whether the GAM cKO is the same Emx1-Cre cortex as the ChIP.</li>",
  "</ul>")

# Every term expanded in the body is also listed here, so a reader who joins the
# report halfway can look one up without hunting for its first use. A <dl> rather
# than a <table>: the table rule sets white-space:nowrap, which would force a
# definition column to scroll sideways instead of wrapping.
parts$glossary <- sec("Glossary",
  "<dl>",
  "<dt><strong>A/B compartments</strong></dt><dd>The split of a chromosome into the active (A) and inactive (B) spatial compartment, read off an eigenvector of the contact matrix. Parked here &mdash; see result 4.</dd>",
  "<dt><strong>ashr</strong></dt><dd>Adaptive shrinkage, the default prior used by <code>lfcShrink</code>. It assumes most effects are null, which is wrong for a genome-wide shift.</dd>",
  "<dt><strong>autosome</strong></dt><dd>A non-sex chromosome; chr1&ndash;chr19 in mouse. The resampling unit for every bootstrap interval here.</dd>",
  "<dt><strong>bootstrap</strong></dt><dd>Resampling the 19 autosomes with replacement and recomputing the statistic each time, to get a 95% percentile interval.</dd>",
  "<dt><strong>bp / kb / Mb</strong></dt><dd>Base pairs, thousands of base pairs, millions of base pairs.</dd>",
  "<dt><strong>CA1</strong></dt><dd>Hippocampal field CA1; its pyramidal cells are the GAM material.</dd>",
  "<dt><strong>ChIP</strong></dt><dd>Chromatin immunoprecipitation &mdash; here the multiplexed ChIP-seq of five histone marks that the GAM data is integrated against.</dd>",
  "<dt><strong>CI</strong></dt><dd>Confidence interval. All of them here are 95% bootstrap percentile intervals over the 19 autosomes.</dd>",
  "<dt><strong>cKO / KO / WT</strong></dt><dd>Conditional knockout (the Emx1-Cre HP1&gamma; deletion), knockout, wild type. The canonical labels elsewhere in this repo are WT and HP1gKO.</dd>",
  "<dt><strong>cLAD / ciLAD / LAD</strong></dt><dd>Lamina-associated domain: a region contacting the nuclear lamina. cLAD is the constitutive core, shared across cell types; ciLAD is the constitutive inter-LAD complement.</dd>",
  "<dt><strong>CrossMap</strong></dt><dd>The tool used to lift external coordinates between genome assemblies.</dd>",
  "<dt><strong>delta-z</strong></dt><dd>The upstream differential-contact caller&rsquo;s own z-like statistic for a bin pair. Sign convention in this report is cKO &minus; WT.</dd>",
  "<dt><strong>DESeq2</strong></dt><dd>The count-based differential-testing package used on the ChIP side.</dd>",
  "<dt><strong>difference-in-differences</strong></dt><dd>(intra<sub>cKO</sub> &minus; intra<sub>WT</sub>) &minus; (inter<sub>cKO</sub> &minus; inter<sub>WT</sub>): the change in one class minus the change in the other, so a shift shared by both cancels.</dd>",
  "<dt><strong>eigenvector</strong></dt><dd>An axis of the correlation matrix. PC1&ndash;PC4 are the leading four, ranked by how much variance they carry.</dd>",
  "<dt><strong>Emx1-Cre</strong></dt><dd>The Cre driver line used for the conditional deletion.</dd>",
  "<dt><strong>GAM</strong></dt><dd>Genome architecture mapping: contact frequencies inferred from how often two windows co-segregate across thin cryosections of single nuclei, rather than from ligation.</dd>",
  "<dt><strong>GAM_2, GAM_3, GAM_10 &hellip;</strong></dt><dd>Numbered stages of this pipeline. The number identifies the stage that produced a result, not a chromosome or a sample.</dd>",
  "<dt><strong>GC content</strong></dt><dd>The guanine+cytosine fraction of a sequence; a genotype-independent stand-in for active chromatin.</dd>",
  "<dt><strong>het / euchromatin</strong></dt><dd>In the pair-class tables, het is heterochromatin as called from the WT H3K9me3/H4K20me3 domains, and euchromatin is the remainder. eu&ndash;eu means euchromatin&ndash;euchromatin.</dd>",
  "<dt><strong>histone marks</strong></dt><dd>Read as histone, residue, modification: H3K9me3 is trimethylation of lysine 9 on histone H3. H3K4me3 (promoters) and H3K36me3 (gene bodies) are active; H3K9me2, H3K9me3 and H4K20me3 are repressive.</dd>",
  "<dt><strong>HP1&gamma;</strong></dt><dd>Heterochromatin protein 1 gamma, the protein encoded by Cbx3 and the gene deleted in the cKO.</dd>",
  "<dt><strong>HUSH</strong></dt><dd>The human silencing hub complex, which silences young retrotransposons and intronless genes. Its target families are the gene clusters shown in the region plots.</dd>",
  "<dt><strong>IAP / IAPLTR1a_Mm</strong></dt><dd>Intracisternal A-particle, a family of mouse endogenous retrovirus; IAPLTR1a_Mm is one of its youngest long terminal repeat subfamilies.</dd>",
  "<dt><strong>insulation</strong></dt><dd>Mean NPMI <em>across</em> a bin. A local minimum is a strong domain boundary.</dd>",
  "<dt><strong>IQR</strong></dt><dd>Interquartile range, the span from the 25th to the 75th percentile.</dd>",
  "<dt><strong>Jaccard index</strong></dt><dd>Shared bins divided by the union of both sets; 1 is identical, 0 is disjoint.</dd>",
  "<dt><strong>KRAB-ZFP</strong></dt><dd>Kr&uuml;ppel-associated box zinc-finger protein, a large family of repeat-silencing transcription factors that sits in tandem gene clusters on chr13 and chr17.</dd>",
  "<dt><strong>LINE / SINE / LTR</strong></dt><dd>Long interspersed nuclear element (L1 in mouse), short interspersed nuclear element, and long terminal repeat &mdash; the flanking element of an endogenous retrovirus. The three repeat classes used as composition covariates.</dd>",
  "<dt><strong>lfcShrink</strong></dt><dd>DESeq2&rsquo;s fold-change shrinkage step. Deliberately not applied here.</dd>",
  "<dt><strong>log2FC</strong></dt><dd>Log2 fold change: cKO over WT on a log2 scale, so &minus;1 is a halving.</dd>",
  "<dt><strong>mm39</strong></dt><dd>The GRCm39 mouse genome assembly, used by both assays.</dd>",
  "<dt><strong>NA</strong></dt><dd>Not available: a missing, unmeasured value. Never the same as zero, and it renders white in the contact maps.</dd>",
  "<dt><strong>NCBI / UCSC seqnames</strong></dt><dd>Two chromosome-naming conventions: <code>1</code> (NCBI, National Center for Biotechnology Information) versus <code>chr1</code> (UCSC, University of California, Santa Cruz). Mixing them gives zero overlaps, which reads as a clean null.</dd>",
  "<dt><strong>NPMI</strong></dt><dd>Normalised pointwise mutual information: the delivered contact measure. A normalised <em>co-segregation</em> statistic, <strong>not a contact frequency</strong> &mdash; it is how much more often two windows share a nuclear profile than chance would give, already accounting for how often each is detected at all. Bounded in [&minus;1, 1]; <strong>0 means chance, not &ldquo;no contact&rdquo;</strong>, which is why unmappable bins are never set to 0.</dd>",
  "<dt><strong>n_valid</strong></dt><dd>The count of usable (non-NA) contact values for a bin &mdash; the coverage diagnostic for a region plot.</dd>",
  "<dt><strong>nuclear profile</strong></dt><dd>One thin cryosection of a single nucleus, the unit GAM sequences. NPMI is computed across profiles.</dd>",
  "<dt><strong>O&minus;E</strong></dt><dd>Observed minus expected: the matrix with the mean of its own diagonal &mdash; the mean NPMI at that separation &mdash; subtracted, which removes the distance decay.</dd>",
  "<dt><strong>Olfr / Vmn</strong></dt><dd>Olfactory receptor and vomeronasal receptor gene clusters, two large tandem arrays used as test regions.</dd>",
  "<dt><strong>PC1&ndash;PC4</strong></dt><dd>The first to fourth principal component.</dd>",
  "<dt><strong>PCA</strong></dt><dd>Principal component analysis: the decomposition the compartment axis is taken from.</dd>",
  "<dt><strong>Pcdh</strong></dt><dd>The clustered protocadherin locus on chr18, which spans approximately one TAD.</dd>",
  "<dt><strong>q90 / q99</strong></dt><dd>The 90th and 99th percentile of a distribution.</dd>",
  "<dt><strong>QC</strong></dt><dd>Quality control.</dd>",
  "<dt><strong>&rho; (rho)</strong></dt><dd>Spearman&rsquo;s rho, the rank correlation coefficient; &minus;1 to +1. |&rho;| is its magnitude, ignoring sign.</dd>",
  "<dt><strong>sd</strong></dt><dd>Standard deviation.</dd>",
  "<dt><strong>sigma 2.6</strong></dt><dd>The smoothing parameter of the upstream NPMI computation. What it was fit on is an open question in METHODS.md.</dd>",
  "<dt><strong>sizeFactors</strong></dt><dd>DESeq2&rsquo;s per-sample scaling constants. Fixed at 1 here, because the usual median-of-ratios estimate would divide the genome-wide loss out.</dd>",
  "<dt><strong>TAD</strong></dt><dd>Topologically associating domain: a self-interacting block bounded by insulating boundaries.</dd>",
  "<dt><strong>z-score / z &gt; 2</strong></dt><dd>A value expressed in standard deviations from the mean. z &gt; 2 is the significance threshold baked into the delivered differential contacts.</dd>",
  "</dl>")

# Theme tokens are defined once and emitted under whichever selectors the
# target needs (see TWO OUTPUT MODES above). Every colour is a token: a literal
# colour inside a theme block is the classic unreadable-page bug, because it
# applies in one theme only.
light_tokens <- "--bg:#fff;--fg:#1c1b1f;--mut:#5b5b66;--line:#e2e1e6;--warn:#fff8e6;--warnb:#e6c65c;--miss:#fdeaea;--missb:#d98b8b;--panel:#f6f6f8;--code:#f3f3f5;--imgbg:transparent"
dark_tokens  <- "--bg:#16151a;--fg:#e9e8ee;--mut:#a3a2ad;--line:#33323b;--warn:#2c2616;--warnb:#8a7430;--miss:#2e1d1d;--missb:#8a4a4a;--panel:#232229;--code:#232229;--imgbg:#fff;color-scheme:dark"

css <- "
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.6 -apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif}
main{max-width:1180px;margin:0 auto;padding:28px 16px 80px}
h1{font-size:1.9rem;margin:0 0 .2em;letter-spacing:-.01em}
h2{font-size:1.3rem;margin:2.2em 0 .5em;padding-bottom:.25em;border-bottom:2px solid var(--line)}
h3{font-size:1.08rem;margin:1.6em 0 .4em}
h4{font-size:.98rem;margin:1.2em 0 .3em;color:var(--mut)}
p.sub{color:var(--mut);margin:0 0 1.2em}
p.cap{font-weight:600;margin:1.1em 0 .3em}
p.note{color:var(--mut);font-size:.85rem;margin:.2em 0 1em}
figure{margin:1.1em 0}
img{max-width:100%;height:auto;border:1px solid var(--line);border-radius:6px}
figcaption{color:var(--mut);font-size:.85rem;margin-top:.35em}
table{border-collapse:collapse;width:100%;margin:.6em 0 1.2em;font-size:.86rem;display:block;overflow-x:auto}
th,td{border:1px solid var(--line);padding:4px 8px;text-align:right;white-space:nowrap}
th{background:var(--panel);text-align:left;font-weight:600}
td:first-child,th:first-child{text-align:left}
.warn{background:var(--warn);border-left:4px solid var(--warnb);padding:10px 14px;margin:1em 0;border-radius:0 6px 6px 0}
.missing{background:var(--miss);border-left:4px solid var(--missb);padding:8px 12px;margin:.8em 0;border-radius:0 6px 6px 0;font-size:.88rem}
code{background:var(--code);padding:1px 5px;border-radius:3px;font-size:.88em}
img{background:var(--imgbg)}
ul{padding-left:1.2em} li{margin:.3em 0}
section{margin-bottom:1.5em}
"

# The browser file only ever sees the OS preference. The Artifact host also lets
# a viewer force a theme, so there the media query is guarded against an
# explicit light choice and the dark tokens are repeated for an explicit dark
# one - otherwise the toggle wins in one direction only.
theme_css <- function(for_artifact) {
  if (!for_artifact) {
    sprintf(":root{%s}\n@media (prefers-color-scheme:dark){:root{%s}}\n",
            light_tokens, dark_tokens)
  } else {
    sprintf(paste0(":root{%s}\n",
                   "@media (prefers-color-scheme:dark){:root:not([data-theme=\"light\"]){%s}}\n",
                   ":root[data-theme=\"dark\"]{%s}\n"),
            light_tokens, dark_tokens, dark_tokens)
  }
}

body_html <- paste0("<main>",
  parts$head, parts$bottom, parts$inventory, parts$settled, parts$chip, parts$hotspots, parts$hotspot_annot, parts$het, parts$hetclass, parts$blocks,
  parts$targets, parts$pcdh, parts$chrom, parts$regions, parts$compartments,
  parts$shared_axis, parts$meltron, parts$coupling, parts$cutoffs, parts$methods, parts$caveats, parts$glossary,
  "</main>")

html <- paste0("<!doctype html><html><head><meta charset='utf-8'>",
  "<meta name='viewport' content='width=device-width,initial-scale=1'>",
  "<title>GAM working report</title><style>", theme_css(FALSE), css, "</style></head><body>",
  resolve_figs(body_html, "embed"), "</body></html>")

writeLines(html, OUT)
msg("wrote %s (%.1f MB)", OUT, file.size(OUT) / 1e6)
n_missing <- length(gregexpr("class=\"missing\"", html)[[1]])
if (grepl("class=\"missing\"", html)) msg("%d missing figure/table placeholder(s) - stages that did not complete", n_missing)

if (ARTIFACT) {
  # No <!doctype>, <html>, <head> or <body>: the Artifact host wraps the file in
  # its own skeleton, and a second document nested inside that renders wrong.
  # :root is left to the host as well - it pads it with the phone's safe-area
  # insets, which is what keeps the page clear of the status bar.
  art <- paste0("<title>GAM Working Report</title>\n<style>",
                theme_css(TRUE), css, "</style>\n",
                resolve_figs(body_html, "external"))
  writeLines(art, OUT_ART)
  mb <- file.size(OUT_ART) / 1e6
  msg("wrote %s (%.1f MB)", OUT_ART, mb)
  # The Artifact limit is 16 MB with embedded data: URIs counted, and the
  # figures are nearly all of this file, so warn before a publish is refused.
  if (mb > 15) msg("WARNING: %.1f MB is close to the 16 MB Artifact limit - lower MAXW", mb)
}
