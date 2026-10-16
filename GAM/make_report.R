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
# ================================================================
source("config.R")
suppressPackageStartupMessages({ library(data.table); library(base64enc) })

OUT <- file.path(results_dir, "GAM_REPORT.html")
MAXW <- 1500L           # px, downscale wider images before embedding
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

md_file <- function(path) {
  if (!file.exists(path)) return(sprintf('<div class="missing">missing: <code>%s</code></div>', esc(path)))
  x <- readLines(path, warn = FALSE)
  out <- character(0); intab <- FALSE
  for (ln in x) {
    if (grepl("^\\|", ln)) {
      cells <- trimws(strsplit(sub("^\\|", "", sub("\\|$", "", ln)), "|", fixed = TRUE)[[1]])
      if (all(grepl("^-+$", cells))) next
      tag <- if (!intab) { intab <- TRUE; out <- c(out, "<table>"); "th" } else "td"
      out <- c(out, paste0("<tr>", paste0("<", tag, ">", esc(cells), "</", tag, ">", collapse = ""), "</tr>"))
    } else {
      if (intab) { out <- c(out, "</table>"); intab <- FALSE }
      if (grepl("^### ", ln))      out <- c(out, paste0("<h4>", esc(sub("^### ", "", ln)), "</h4>"))
      else if (grepl("^## ", ln))  out <- c(out, paste0("<h3>", esc(sub("^## ", "", ln)), "</h3>"))
      else if (grepl("^# ", ln))   out <- c(out, paste0("<h2>", esc(sub("^# ", "", ln)), "</h2>"))
      else if (grepl("^- ", ln))   out <- c(out, paste0("<li>", bold(esc(sub("^- ", "", ln))), "</li>"))
      else if (nzchar(trimws(ln))) out <- c(out, paste0("<p>", bold(esc(ln)), "</p>"))
    }
  }
  if (intab) out <- c(out, "</table>")
  paste(out, collapse = "\n")
}
# perl = TRUE for the lazy quantifier; R's default TRE does not support it.
bold <- function(x) gsub("\\*\\*(.+?)\\*\\*", "<strong>\\1</strong>", x, perl = TRUE)

sec <- function(title, ...) paste0("<section><h2>", title, "</h2>", paste0(..., collapse = "\n"), "</section>")
p_   <- function(...) paste0("<p>", paste0(..., collapse = " "), "</p>")

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
  tbl("gam15_chromosome_ranking.tsv", "All 19 autosomes, so any one chromosome can be read against the spread.", max_rows = 20L),
  embed_fig("chromosomes/contact_budget_by_chromosome.png"),
  tbl("gam15_contact_budget.tsv", "Contact budget by chromosome, chromatin class and separation.", max_rows = 60L),
  paste0(sapply(c("chr3","chr7","chr11","chr16","chr18"), function(ch)
    paste0("<h3>", ch, "</h3>",
           embed_fig(sprintf("chromosomes/%s_overview.png", ch)),
           md_file(file.path(results_dir, "chromosomes", paste0(ch, ".md"))))), collapse = ""))

reg <- c("chr18_Pcdh", "chr18_Pcdh_band", "chr13_KRABZFP_band", "chr17_KRABZFP_band",
         "chr7_Vmn_band", "chr2_Olfr_band", "chr16_full_band", "chr11_1-25Mb_lostblock_band",
         "chr3_30-85Mb_band")
parts$regions <- sec("Region plots",
  p_("Contact map, differential contacts, hotspot calls, LADs, the five marks and genes, on one coordinate axis.",
     "Conventions: raw NPMI is a sequential 0&rarr;q99 ramp with one cap shared by both genotypes; delta-z is fixed at &plusmn;4;",
     "<strong>net diff. contacts</strong> is the per-bin count of significant differential contacts within 5 Mb (gained &minus; lost), while <strong>hotspot</strong> is the separate call on that count at the 50-per-window cutoff."),
  p_("Mark tracks are aggregated to 100 kb and faded by coverage, because per-bin log2FC at 20 kb is noise-dominated. Fill diverges from each mark's <em>genome median</em>, printed in the axis label &mdash; blue means above that median, which for a negative median is still a loss."),
  paste0(sapply(reg, function(r) embed_fig(sprintf("regions/%s.png", r), esc(r))), collapse = ""))

parts$coupling <- sec("Coupling: where it stands",
  "<div class='warn'>Three framings of one question currently coexist and disagree in presentation, not in substance. This section states them once; the figures themselves are still the old ones and are the lowest-priority item to redo.</div>",
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

css <- "
:root{--bg:#fff;--fg:#1c1b1f;--mut:#5b5b66;--line:#e2e1e6;--warn:#fff8e6;--warnb:#e6c65c;--miss:#fdeaea;--missb:#d98b8b}
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
th{background:#f6f6f8;text-align:left;font-weight:600}
td:first-child,th:first-child{text-align:left}
.warn{background:var(--warn);border-left:4px solid var(--warnb);padding:10px 14px;margin:1em 0;border-radius:0 6px 6px 0}
.missing{background:var(--miss);border-left:4px solid var(--missb);padding:8px 12px;margin:.8em 0;border-radius:0 6px 6px 0;font-size:.88rem}
code{background:#f3f3f5;padding:1px 5px;border-radius:3px;font-size:.88em}
ul{padding-left:1.2em} li{margin:.3em 0}
section{margin-bottom:1.5em}
@media (prefers-color-scheme:dark){
 :root{--bg:#16151a;--fg:#e9e8ee;--mut:#a3a2ad;--line:#33323b;--warn:#2c2616;--warnb:#8a7430;--miss:#2e1d1d;--missb:#8a4a4a;color-scheme:dark}
 th{background:#232229} code{background:#232229} img{background:#fff}
}
"

html <- paste0("<!doctype html><html><head><meta charset='utf-8'>",
  "<meta name='viewport' content='width=device-width,initial-scale=1'>",
  "<title>GAM working report</title><style>", css, "</style></head><body><main>",
  parts$head, parts$het, parts$pcdh, parts$chrom, parts$regions, parts$coupling, parts$caveats,
  "</main></body></html>")

writeLines(html, OUT)
msg("wrote %s (%.1f MB)", OUT, file.size(OUT) / 1e6)
n_missing <- length(gregexpr("class=\"missing\"", html)[[1]])
if (grepl("class=\"missing\"", html)) msg("%d missing figure/table placeholder(s) - stages that did not complete", n_missing)
