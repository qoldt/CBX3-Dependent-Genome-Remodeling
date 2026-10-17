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
MAXW <- 1100L           # px, downscale wider images before embedding (keeps the
                        # single-file HTML small enough to open on a phone)
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
# is embedded at full size rather than skipped.
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
  b64 <- base64enc::base64encode(src)
  paste0('<figure><img src="data:image/png;base64,', b64, '">',
         if (!is.null(caption)) paste0("<figcaption>", caption, "</figcaption>") else "",
         "</figure>")
}

# --- table rendering ----------------------------------------------------------
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
  "<p class='sub'>HP1&gamma; (Cbx3) cKO vs WT &mdash; GAM CA1 pyramidal cells, ChIP Emx1-Cre cortex+hippocampus. mm39, 20 kb. Generated ",
  Sys.Date(), " by <code>make_report.R</code>.</p>",
  "<div class='warn'><strong>Read first.</strong> GAM is 2 biological replicates per genotype, already merged, so there is no within-genotype variance: ",
  "every number here is an effect size against a matched or bootstrap null, never a per-contact p-value. ",
  "Per-bin ChIP log2FC at 20 kb is noise-dominated (per-bin sd 0.14&ndash;0.16 for the broad marks, 0.57 for H3K4me3); genome-wide medians are sound, local structure is not. ",
  "The A/B compartment axis is weakly determined on most chromosomes and is not used for any claim below.</div>")

parts$bottom <- sec("Bottom line",
  "<p>Everything below this section is the supporting detail. Percentages are change in mean NPMI,",
  "cKO vs WT, with 95% CIs from a bootstrap over the 19 autosomes",
  "(<code>gam10_class_change.tsv</code>, <code>gam13_did.tsv</code>).</p>",

  "<h3>1. H4K20me3 is lost almost everywhere. That part is solid.</h3>",
  "<p>Median log2FC about &minus;0.50 per 20 kb bin, &minus;0.68 over peaks, ~95% of domains down.",
  "A global shift, not a minority-of-sites effect. Everything else here is the architecture question.</p>",

  "<h3>2. Heterochromatin becomes more self-associated at short range. It is NOT mixing.</h3>",
  "<p>This is the central result and the loose phrase &ldquo;heterochromatin gains contacts&rdquo; hides which partners.",
  "Split by partner class:</p>",
  "<table>",
  "<tr><th>pair type</th><th>&lt;0.2 Mb</th><th>0.2&ndash;1 Mb</th><th>5&ndash;10 Mb</th><th>&gt;10 Mb</th></tr>",
  "<tr><td>het&ndash;het, same domain</td><td><strong>+4.3%</strong> (1.8, 6.6)</td><td><strong>+5.9%</strong> (1.3, 18.6)</td><td>&mdash;</td><td>&mdash;</td></tr>",
  "<tr><td>het&ndash;het, different domains</td><td><strong>+4.2%</strong> (3.4, 4.9)</td><td><strong>+6.2%</strong> (3.5, 8.2)</td><td>+5.1% (&minus;5.9, 9.0)</td><td>&minus;0.6% (&minus;10.3, 4.9)</td></tr>",
  "<tr><td>het&ndash;euchromatin</td><td>+2.0% (1.3, 2.8)</td><td>+0.2% (&minus;1.2, 1.6)</td><td><strong>&minus;5.2%</strong> (&minus;10.1, &minus;1.0)</td><td>&minus;6.4% (&minus;13.5, 0.0)</td></tr>",
  "<tr><td>euchromatin&ndash;euchromatin</td><td>+0.9% (0.3, 1.6)</td><td><strong>&minus;1.7%</strong> (&minus;3.0, &minus;0.5)</td><td><strong>&minus;8.1%</strong> (&minus;12.0, &minus;4.2)</td><td><strong>&minus;8.9%</strong> (&minus;13.8, &minus;4.0)</td></tr>",
  "</table>",
  "<p>Read down the columns, not across:</p>",
  "<ul>",
  "<li><strong>The gain is specifically heterochromatin-to-heterochromatin, and only below 1 Mb</strong> &mdash; those are the",
  "only cells whose CIs exclude zero. Beyond 1 Mb the het&ndash;het CIs span zero and nothing can be claimed.</li>",
  "<li><strong>It is not mixing.</strong> Mixing would require het&ndash;euchromatin to rise. It does not: flat at 0.2&ndash;1 Mb",
  "(+0.2%, CI spans zero) and <em>falling</em> at long range (&minus;5.2%, &minus;6.4%). Heterochromatin is keeping more",
  "company with heterochromatin and no more with euchromatin &mdash; segregation sharpening at short range, not dissolving.</li>",
  "<li><strong>At long range everything loses, euchromatin most</strong> (&minus;8 to &minus;9%). The &ldquo;het gains while eu loses,",
  "diverging with distance&rdquo; shorthand conflates two separate things: a short-range het&ndash;het gain, and a",
  "long-range loss that hits euchromatin hardest.</li>",
  "</ul>",

  "<h3>3. Individual domains do not compact. Heterochromatin-to-heterochromatin contact is what rises.</h3>",
  "<p>Within-domain versus between-domain is indistinguishable: every difference-in-differences CI spans zero",
  "(p 0.19&ndash;0.89) under all three domain definitions. So the right statement is &ldquo;heterochromatin contacts",
  "heterochromatin more&rdquo;, not &ldquo;each domain tightens&rdquo;. That null replicates across definitions, which is",
  "what makes it believable rather than a definition artefact.</p>",

  "<h3>4. Whether compartments are better or worse defined: we cannot say.</h3>",
  "<p>The A/B axis (PC1) is usable on only 11 of 38 chromosome &times; genotype combinations, so no compartment-strength",
  "claim is supported either way. The class table above is <em>consistent with</em> sharper short-range segregation, but",
  "that is an inference from contact classes, not a compartment measurement. Answering it properly needs a different",
  "axis &mdash; the cLAD/ciLAD annotation is the obvious external candidate.</p>",

  "<h3>5. The remodelling does not track the size of the mark loss, bin by bin.</h3>",
  "<p>rho = &minus;0.015 raw, &minus;0.046 partialled on coverage, n = 119,952 &mdash; a powered negative. But the contact",
  "variable is smeared over 5 Mb while the ChIP variable varies at 20 kb, a ~250&times; scale mismatch, so this cannot",
  "separate &ldquo;no relationship&rdquo; from &ldquo;cannot resolve it&rdquo;. The architecture does change; it just does not",
  "change in proportion to how much mark a bin lost.</p>",

  "<h3>6. Chromosomes split by direction, and chr11 runs the other way.</h3>",
  "<p>chr3 and chr7 gain; chr11, chr16 and chr18 lose. Genome-wide, het&ndash;het minus eu&ndash;eu is positive.",
  "On chr11 it is <strong>negative</strong> &mdash; &minus;2.3 at 1&ndash;5 Mb and &minus;11.8 at 5&ndash;10 Mb &mdash; i.e. its",
  "heterochromatin loses contacts <em>faster</em> than its euchromatin, which is the inversion.</p>",
  "<div class=\'warn\'><strong>Do not over-read chr11.</strong> It is the least heterochromatic autosome (3.1% of bins",
  "against 7.2% genome-wide), so its het&ndash;het cells rest on ~1,000&ndash;1,200 bin pairs where eu&ndash;eu has millions.",
  "Its dominant feature is a loss across <em>every</em> class (eu&ndash;eu &minus;22% at 5&ndash;10 Mb, &minus;27.8% beyond",
  "10 Mb); the het&ndash;het number is simply more negative still, on thin data. Nothing we have computed attributes this",
  "to the loss of one large LAD-associated domain, and a large heterochromatin domain is unlikely there by construction.</div>",

  "<h3>7. Four things tested and rejected. Recorded so they are not re-derived.</h3>",
  "<ul>",
  "<li><strong>Heterochromatin-rich chromosomes gain contacts.</strong> No: rho = &minus;0.014, p = 0.96 across 19 autosomes.",
  "chr7 and chr11 are simply the extremes of both rankings at once.</li>",
  "<li><strong>The Pcdh TAD changes its distal contacts distinctively.</strong> No: p = 0.75/0.85 intra and 0.88/0.68 distal",
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
  tbl("gam13_definition_agreement.tsv", "Agreement between them (Jaccard on bin membership). For scale, the two hotspot-domain variants agree at 0.40."),
  embed_fig("het_domains/intra_vs_inter_by_definition.png",
            "Within-domain vs between-domain contact change, per definition. Bars are 95% CIs from a bootstrap over the 19 autosomes."),
  embed_fig("het_domains/by_domain_size.png"),
  tbl("gam13_did.tsv", "Difference-in-differences: intra minus inter. A CI excluding zero means the two classes changed by different amounts."),
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
  "<div class='warn'><strong>Retraction.</strong> An earlier agent reported chr11 as having &ldquo;the strongest mark-contact coupling of any chromosome, rho +0.51 to +0.66, rank 1 of 19&rdquo;. That agent produced no output file and the numbers had no source. Directly computed, chr11 <em>is</em> rank 1 of 19 &mdash; but at rho = 0.128, with chr17 at &minus;0.149 (larger, opposite sign) and a genome-wide value of 0.012. The ordering is weak evidence.</div>",
  embed_fig("chromosomes/chromosome_ranking.png"),
  "<div class='warn'><strong>A tempting over-read, tested and rejected.</strong> chr7 is by far the most heterochromatic chromosome (20.2% of bins against 7.2% genome-wide) <em>and</em> gains contacts at every separation, while chr11 is the least (3.1%) and loses most. Across all 19 autosomes that relationship is absent: rho = &minus;0.014, p = 0.96. The two are simply the extremes of both rankings at once.</div>",
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
     "Conventions: raw NPMI is a sequential 0&rarr;q99 ramp with one cap shared by both genotypes; delta-z is fixed at &plusmn;4;",
     "<strong>net diff. contacts</strong> is the per-bin count of significant differential contacts within 5 Mb (gained &minus; lost), while <strong>hotspot</strong> is the separate call on that count at the 50-per-window cutoff."),
  "<div class='warn'><strong>Data coverage, read before interpreting any empty map.</strong> A region whose bins are detected in few nuclear profiles has mostly NA in the NPMI matrix, and NA renders white &mdash; indistinguishable from &ldquo;no contacts here&rdquo;. Each plot now states its measured fraction, and warns below 75%. Measured here: Pcdh 97%, Olfr 92%, KRAB-ZFP chr17 87%, <strong>KRAB-ZFP chr13 61%</strong>, <strong>vomeronasal chr7 29%</strong>. The vomeronasal cluster's apparently empty contact map is <em>missing data, not absent contact</em> &mdash; 11% finite cells inside the cluster against 97% in its own flanks, with n_valid 807 vs 6,574 (genome-wide 6,104). Do not read architecture off it.</div>",
  p_("Mark tracks are aggregated to 100 kb and faded by coverage, because per-bin log2FC at 20 kb is noise-dominated. Fill diverges from each mark's <em>genome median</em>, printed in the axis label &mdash; blue means above that median, which for a negative median is still a loss."),
  paste0(sapply(reg, function(r) embed_fig(sprintf("regions/%s.png", r), esc(r))), collapse = ""))

parts$coupling <- sec("Coupling: does the mark loss have a contact counterpart?",
  p_("Three framings of one question, which until now lived in three files and disagreed only in presentation.",
     "They are now one stage, and this figure is the synthesis:"),
  embed_fig("coupling/coupling_synthesis.png"),
  p_("In words:"),
  "<ul>",
  "<li><strong>GAM_3</strong> &mdash; per-bin Spearman of ChIP log2FC against differential-contact count. Returns a clean null: rho &minus;0.015 raw, &minus;0.046 partial on coverage, flat under every variant. Bins losing the most H4K20me3 (median &minus;0.764) and the least (&minus;0.222) have identical contact densities. n = 119,952, so this is a powered negative, not an underpowered one.</li>",
  "<li><strong>GAM_7</strong> &mdash; retracts GAM_3's premise: the contact variable is smeared over 5 Mb while the ChIP variable varies at 20 kb, a ~250&times; scale mismatch. The null is real for the question as posed but cannot distinguish &ldquo;no relationship&rdquo; from &ldquo;the contact variable cannot resolve that scale&rdquo;.</li>",
  "<li><strong>GAM_10/11/13</strong> &mdash; NPMI difference-in-differences with a bootstrap over the 19 autosomes. This is the framing that survived, and the one every new stage uses.</li>",
  "</ul>",
  p_("The positive result is elsewhere: H4K20me3 loss tracks <em>baseline</em> architecture, not contact <em>change</em> &mdash; pc1 +0.196 and npmi_short &minus;0.166 after coverage control. Insulation does not survive the control (&minus;0.180 &rarr; &minus;0.066) and should not be reported."),
  embed_fig("coupling/specificity_control_rho.png"))

parts$caveats <- sec("Standing caveats",
  "<ul>",
  "<li><strong>2 vs 2.</strong> Both assays are two biological replicates per genotype. The GAM files are already merged and carry no replicate dimension, so z &gt; 2 is the only significance available from them.</li>",
  "<li><strong>ChIP at 20 kb.</strong> DESeq2 sees integers of ~1&ndash;11 per bin after the &times;10 count scaling, unshrunk. H3K36me3 appeared raised over Pcdh; it is &minus;0.049 against a genome median of &minus;0.052, p = 0.93. Nothing at Pcdh differs from genome background at 20 kb. lfcShrink is <em>not</em> the fix &mdash; it assumes most effects are null, which is exactly what sizeFactors = 1 exists to avoid here.</li>",
  "<li><strong>Compartments parked.</strong> A usable axis on 11 of 38 chromosome &times; genotype combinations.</li>",
  "<li><strong>LADs are constitutive, not neuronal.</strong> cLAD/ciLAD are the cell-type-invariant core across ESC/NPC/astrocyte/MEF.</li>",
  "<li><strong>Unresolved provenance.</strong> Nuclear profiles per genotype, the passed_qc.curated criteria, what sigma 2.6 was fit on, and whether the GAM cKO is the same Emx1-Cre cortex as the ChIP.</li>",
  "</ul>")

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
  parts$head, parts$bottom, parts$het, parts$pcdh, parts$chrom, parts$regions, parts$coupling, parts$caveats,
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
