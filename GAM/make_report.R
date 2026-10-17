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

embed_fig <- function(rel, caption = NULL) {
  f <- file.path(fig_dir, rel)
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
  tag <- if (ARTIFACT) {
    pub <- file.path("figures", gsub("[/]", "_", rel))
    dest <- file.path(asset_dir, pub)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    file.copy(src, dest, overwrite = TRUE)
    ASSETS$files <- c(ASSETS$files, pub)
    # loading="lazy" matters here: the artifact is read on a phone, and without it
    # every figure in the report is fetched before the first one renders.
    sprintf('<img src="%s" loading="lazy" decoding="async">', pub)
  } else {
    sprintf('<img src="data:image/png;base64,%s">', base64enc::base64encode(src))
  }
  paste0("<figure>", tag,
         if (!is.null(caption)) paste0("<figcaption>", caption, "</figcaption>") else "",
         "</figure>")
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

  "<h3>4. A compartment axis is recovered on every autosome, but not the same eigenvector in both genotypes.</h3>",
  "<p>Per <code>npmi_compartment_orientation.tsv</code>, an oriented compartment axis is recovered on <strong>all 19",
  "autosomes in both genotypes</strong> (38 of 38; only chrX and chrY fail). Its agreement with the orienting",
  "marker &mdash; the mean of per-bin WT H3K4me3 and H3K36me3 z-scores (signal in standard deviations from the genome-wide mean) &mdash; is modest.",
  "Agreement is Spearman&rsquo;s &rho; (rho, the rank correlation coefficient), and |&rho;| is its magnitude, ignoring sign:</p>",
  "<table>",
  "<tr><th>|&rho;| with orienting marker</th><th>combos (of 38)</th></tr>",
  "<tr><td>&ge; 0.3</td><td>25</td></tr><tr><td>&ge; 0.4</td><td>11</td></tr><tr><td>&ge; 0.5</td><td>3</td></tr>",
  "</table>",
  "<p>The more diagnostic fact is <em>which</em> eigenvector carries it. The compartment axis is <strong>PC1 &mdash; the first",
  "principal component &mdash; on only chr7</strong> (both genotypes); elsewhere it is PC2 (the second principal component,",
  "17 combos), PC3 (13) or PC4 (6). GAM_2 (a numbered stage of this pipeline) already subtracts",
  "expected-by-distance before the Pearson correlation and PCA (principal component analysis), so the usual culprit &mdash; the distance-decay gradient",
  "dominating PC1 &mdash; is already handled. Something else is taking the leading variance, and we have not identified",
  "it. Two candidates, both testable and neither tested: per-bin detection frequency, and the large one-directional",
  "GAM_6 blocks.</p>",
  "<p>So the question is <em>answerable</em>, not closed &mdash; but it needs that diagnostic first, because a",
  "compartment-strength comparison built on PC3 of one genotype against PC2 of the other is not measuring one thing.",
  "Each genotype is currently oriented independently, so an apparent A/B switch could be an orientation artefact; orienting",
  "both on a genotype-independent marker (GC content &mdash; the guanine+cytosine fraction &mdash; or gene density, rather than H3K4me3 signal) would remove that",
  "degree of freedom. cLAD/ciLAD &mdash; constitutive lamina-associated domains (LADs) and constitutive inter-LADs &mdash; remain the external",
  "alternative that avoids the axis entirely.</p>",

  "<h3>5. The remodelling does not track the size of the mark loss, bin by bin.</h3>",
  "<p>&rho; = &minus;0.015 raw, &minus;0.046 partialled on coverage, n = 119,952 &mdash; a powered negative. But the contact",
  "variable is smeared over 5 Mb while the ChIP variable varies at 20 kb, a ~250&times; scale mismatch, so this cannot",
  "separate &ldquo;no relationship&rdquo; from &ldquo;cannot resolve it&rdquo;. The architecture does change; it just does not",
  "change in proportion to how much mark a bin lost.</p>",

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
  p_("Two partitions are reported side by side &mdash; published boundaries (fixed, non-circular, different cell type) and insulation minima from our own GAM (self-consistent, free to move between genotypes, but our own call).",
     "A result that holds under only one is a result about the partition."),
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
  embed_fig("chromosomes/contact_budget_by_chromosome.png"),
  tbl("gam15_contact_budget.tsv", "Contact budget by chromosome, chromatin class and separation.", max_rows = 60L),
  paste0(sapply(c("chr3","chr7","chr11","chr16","chr18"), function(ch)
    paste0("<h3>", ch, "</h3>",
           embed_fig(sprintf("chromosomes/%s_overview.png", ch)),
           md_file(file.path(results_dir, "chromosomes", paste0(ch, ".md")),
                   demote = TRUE))), collapse = ""))

# chr18_Pcdh has no _band variant by design: at 1.8 Mb the full triangle is
# already compact, so band_maxdist() skips it rather than writing a duplicate
# under a misleading name.
reg <- c("chr18_Pcdh", "chr18_Pcdh_context_band", "chr13_KRABZFP_band", "chr17_KRABZFP_band",
         "chr7_Vmn_band", "chr2_Olfr_band", "chr16_full_band", "chr11_1-25Mb_lostblock_band",
         "chr3_30-85Mb_band")
parts$regions <- sec("Region plots",
  p_("Contact map, differential contacts, hotspot calls, LADs, the five marks and genes, on one coordinate axis.",
     "Conventions: raw NPMI is a sequential 0&rarr;q99 (the 99th percentile) ramp with one cap shared by both genotypes; delta-z (the upstream caller&rsquo;s contact-change statistic, defined under methods below) is fixed at &plusmn;4;",
     "<strong>net diff. contacts</strong> is the per-bin count of significant differential contacts within 5 Mb (gained &minus; lost), while <strong>hotspot</strong> is the separate call on that count at the 50-per-window cutoff."),
  "<div class='warn'><strong>Data coverage, read before interpreting any empty map.</strong> A region whose bins are detected in few nuclear profiles has mostly NA (not available &mdash; missing, not measured) in the NPMI matrix, and NA renders white &mdash; indistinguishable from &ldquo;no contacts here&rdquo;. Each plot now states its measured fraction, and warns below 75%. Measured here: Pcdh 97%, Olfr (olfactory receptor cluster) 92%, KRAB-ZFP (Kr&uuml;ppel-associated box zinc-finger protein cluster) chr17 87%, <strong>KRAB-ZFP chr13 61%</strong>, <strong>vomeronasal chr7 29%</strong>. The vomeronasal cluster's apparently empty contact map is <em>missing data, not absent contact</em> &mdash; 11% finite cells inside the cluster against 97% in its own flanks, with n_valid (the count of usable contact values per bin) 807 vs 6,574 (genome-wide 6,104). Do not read architecture off it.</div>",
  p_("Mark tracks are aggregated to 100 kb and faded by coverage, because per-bin log2FC at 20 kb is noise-dominated. Fill diverges from each mark's <em>genome median</em>, printed in the axis label &mdash; blue means above that median, which for a negative median is still a loss."),
  paste0(sapply(reg, function(r) embed_fig(sprintf("regions/%s.png", r), esc(r))), collapse = ""))

parts$coupling <- sec("Coupling: does the mark loss have a contact counterpart?",
  p_("Three framings of one question, which until now lived in three files and disagreed only in presentation.",
     "They are now one stage, and this figure is the synthesis:"),
  embed_fig("coupling/coupling_synthesis.png"),
  p_("In words:"),
  "<ul>",
  "<li><strong>GAM_3</strong> &mdash; per-bin Spearman of ChIP log2FC against differential-contact count. Returns a clean null: &rho; &minus;0.015 raw, &minus;0.046 partial on coverage, flat under every variant. Bins losing the most H4K20me3 (median &minus;0.764) and the least (&minus;0.222) have identical contact densities. n = 119,952, so this is a powered negative, not an underpowered one.</li>",
  "<li><strong>GAM_7</strong> &mdash; retracts GAM_3's premise: the contact variable is smeared over 5 Mb while the ChIP variable varies at 20 kb, a ~250&times; scale mismatch. The null is real for the question as posed but cannot distinguish &ldquo;no relationship&rdquo; from &ldquo;the contact variable cannot resolve that scale&rdquo;.</li>",
  "<li><strong>GAM_10/11/13</strong> &mdash; NPMI difference-in-differences with a bootstrap over the 19 autosomes. This is the framing that survived, and the one every new stage uses.</li>",
  "</ul>",
  p_("The positive result is elsewhere: H4K20me3 loss tracks <em>baseline</em> architecture, not contact <em>change</em> &mdash; pc1 +0.196 and npmi_short &minus;0.166 after coverage control. Insulation does not survive the control (&minus;0.180 &rarr; &minus;0.066) and should not be reported."),
  embed_fig("coupling/specificity_control_rho.png"))

parts$methods <- sec("How every comparison in this report is computed",
  "<p>Every number above is one of a small number of operations. This section states each one exactly,",
  "including what is and is not normalised, so a figure can be checked without reading the scripts.</p>",

  "<h3>The contact measure</h3>",
  "<ul>",
  "<li><strong>NPMI (normalised pointwise mutual information)</strong> is the delivered quantity, one 20 kb &times;",
  "20 kb matrix per chromosome per genotype, computed upstream at sigma 2.6 from the nuclear-profile segregation",
  "table. It is a co-segregation frequency already normalised for how often each of the two windows is detected",
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
  "<dt><strong>NPMI</strong></dt><dd>Normalised pointwise mutual information: the delivered contact measure, a co-segregation frequency already normalised for how often each window is detected at all. Bounded in [&minus;1, 1]; 0 means chance.</dd>",
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
  parts$head, parts$bottom, parts$het, parts$pcdh, parts$chrom, parts$regions, parts$coupling, parts$methods, parts$caveats, parts$glossary,
  "</main>")

html <- paste0("<!doctype html><html><head><meta charset='utf-8'>",
  "<meta name='viewport' content='width=device-width,initial-scale=1'>",
  "<title>GAM working report</title><style>", theme_css(FALSE), css, "</style></head><body>",
  body_html, "</body></html>")

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
                theme_css(TRUE), css, "</style>\n", body_html)
  writeLines(art, OUT_ART)
  mb <- file.size(OUT_ART) / 1e6
  msg("wrote %s (%.1f MB)", OUT_ART, mb)
  # The Artifact limit is 16 MB with embedded data: URIs counted, and the
  # figures are nearly all of this file, so warn before a publish is refused.
  if (mb > 15) msg("WARNING: %.1f MB is close to the 16 MB Artifact limit - lower MAXW", mb)
}
