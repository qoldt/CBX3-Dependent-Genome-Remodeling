# CBX3-Dependent Genome Remodeling — repo-level conventions

Multi-assay study of **HP1γ (CBX3)** in neuronal chromatin organization:
**HP1γ conditional knockout (Emx1-Cre) cortex vs wild type**, genome **mm39**
throughout. This file holds the conventions that span assays. Each assay folder
has its own `CLAUDE.md` for its specifics — **read that one before working inside
it**; this file deliberately does not duplicate it.

## Layout

| Path | Assay | Status |
|------|-------|--------|
| `ChIP/` | Multiplexed ChIP-seq (MINUTE), 5 histone marks over called peaks | **Active** — 7 stages, documented in `ChIP/CLAUDE.md` |
| `GAM/` | Genome Architecture Mapping, 20 kb NPMI | **Active** — 9 numbered stages (`GAM_1`–`GAM_7`, plus `GAM_2b`/`GAM_2c`), `run_GAM.R`, `make_report.R`, `verify_GAM.R`, and utilities; documented in `GAM/CLAUDE.md` |
| `shared/` | Cross-assay helpers (`common.R`: seqname flips, `mm39_bins()`, GAM bin-id parsing; `setup.R`) | **Exists** — sourced by both assays' `config.R` |
| `scRNAseq/` | Single-cell RNA-seq | Planned — directory does not exist yet |

The headline result the other assays are being integrated against: HP1γ loss
causes a **broad, near-genome-wide loss of H4K20me3** (median log2FC ≈ −0.68),
deepest over IAP/young-L1 and quiescent chromatin, **largely uncoupled from
H3K9me3** (≈ −0.42), derepressing HUSH/CBX3-target gene clusters. Full statement
with its caveats: `ChIP/results/EXECUTIVE_SUMMARY.md` — read the threshold and QC
caveats there before quoting any number.

## Cross-assay conventions

**Genome build is mm39 everywhere.** External data in mm10 must be lifted with
CrossMap and the liftover recorded in that assay's `METHODS.md`.

**Genotype labels.** Canonical is **`WT`** and **`HP1gKO`**. The GAM files on disk
say **`cKO`** for the same genotype. Translate at the edge — in the assay's
`config.R` and nowhere else — so no analysis code carries two names for one thing.

**Seqname style is the most common silent bug in this repo.** The ChIP
counts/DESeq2 path is **NCBI** (`1`, `2`, …); annotation overlaps and the GAM
files are **UCSC** (`chr1`). A mismatched flip produces zero overlaps, which
reads downstream as a clean null result rather than an error. Any cross-assay
join must assert its overlap rate with `stopifnot()`, not a warning.

**Each assay folder has one `config.R` that is the single source of truth** for
that assay's parameters, paths, thresholds and sample sheet. Stage scripts source
it and must not redefine its values. To change a cutoff, edit `config.R`; to
change samples, edit the sample sheet.

**Stages are numbered, thematic, and independently runnable** after the stage that
writes their input. Each reads its inputs from disk rather than relying on
in-memory state, so any one can be re-run alone. A `run_<ASSAY>.R` at the folder
root executes them in order.

**Figure colours come from the `ltc` palettes defined in config, never raw hex.**
Discrete → `gaby` (≤5 levels) then `hat`; continuous/heatmaps → `heatmap0`, with
the diverging variants for log2FC and z-scores. Greys (`grey65`/`grey75`/`grey80`)
are deliberate for "neutral / background / not-a-category" levels and are *not*
drawn from the palette. Figures save as both `.png` and `.pdf` via `save_fig()` /
`save_base_fig()` into per-analysis subfolders.

**Results trees are git-ignored in full and must stay fully regenerable** from the
committed code plus the externally-hosted raw data. Nothing in `results/` is a
source of truth except the hand-written `.md` files at its top level, which
**are** committed. Currently six, and each is named individually in the
`.gitignore` except-list:

- `ChIP/results/` — `EXECUTIVE_SUMMARY.md`, `METHODS.md`, `COMPARTMENT_CUTOFFS.md`
- `GAM/results/` — `METHODS.md`, `COMPARTMENTS.md`, `BIN_BALANCING.md`

**Add a new one to that list in the same commit that creates it.**
`BIN_BALANCING.md` was silently git-ignored for a commit because the GAM
except-list named only `METHODS.md`, and a committed methods document that is
not actually committed looks identical to one that is.

**A number written in prose is the most common stale thing in this repo**, and
no `stopifnot()` catches it, because each stage is individually correct — the
inconsistency exists only *between* files. The published GAM report once claimed
"38 of 38" in a paragraph citing a table that said 27 of 38.
`GAM/verify_GAM.R` enforces the rule that `results/tables/*.tsv` is the single
source of truth: a statistic is computed *from* the table and every file quoting
it must contain the current value, so the check cannot itself go stale. Run it
before quoting any GAM number (`cd GAM && Rscript verify_GAM.R`); on a failure
the **table is right by definition**. If another assay grows the same problem,
copy the pattern rather than inventing a second one.

**Large raw inputs are never committed** — bigWigs, RepeatMasker annotation, NPMI
matrices, contact files. They are downloaded or staged locally per the assay
README and ignored in `.gitignore`. Small committed inputs (peak BEDs) are the
exception. Before adding any new data folder, check its size and ignore it first.

**Statistical framing carries across assays.** The effect here is genome-wide, so
per-feature testing is underpowered by construction and absolute cutoffs are not
shift-invariant. Prefer **effect sizes expressed relative to each mark's own
genome-wide median**, and **size- or composition-matched permutation nulls** for
anything overlap-shaped (large heterochromatin domains trivially overlap narrow
features). Report enrichments both per-region and size-weighted where they can
diverge.

## Working practice

- **Push directly to `main`** — no feature branches or PRs in this repo.
- Heavy local jobs (bigWig quantification over many regions, parsing the GAM
  matrices) are sustained all-core work on an Apple Silicon Mac — run them under
  the `safe-local-alignment` skill, which handles thermal/swap protection and
  refuses silently-empty output.
- **Check free disk before launching one.** That skill guards thermal pressure
  and swap but *not* disk. On 2026-10-18 the data volume hit 100% (279 MB free
  of 926 GB) while a pipeline was writing, and it surfaced as an unrelated
  `ENOSPC`, not as a pipeline error — R exits 0 after a failed write, so a stage
  can produce a truncated or empty table and still look successful.
  `df -h /System/Volumes/Data`; treat under ~10 GB free as a reason to stop.
- **A stage's own "done" message is not the end of the stage.** The GAM stage
  files are consolidations of an older, finer numbering and print one completion
  message per legacy section, so `GAM_4` announces "GAM_6 complete" and keeps
  going. The reliable end-of-run signals are the thermal guard's
  `FINISHED status=0` and an absent R process
  (`ps -eo comm | grep -cE 'exec/R$'`). Verifying or publishing off a "done"
  line once led to a push made while GAM_4 was still rewriting tables.
- Don't invent data provenance. If the accession, replicate count or QC criterion
  isn't recoverable from the files, write the open question into that assay's
  `METHODS.md` and ask.
