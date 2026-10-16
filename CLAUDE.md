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
| `GAM/` | Genome Architecture Mapping, 20 kb NPMI | **Data only** — raw drop, no pipeline yet |
| `shared/` | Cross-assay helpers (palettes, figure savers, permutation tests) | **Planned** — not yet created |
| `scRNAseq/` | Single-cell RNA-seq | Planned |

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
source of truth except the hand-written `.md` files at its top level
(`EXECUTIVE_SUMMARY.md`, `METHODS.md`, `COMPARTMENT_CUTOFFS.md`), which **are**
committed.

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
- Don't invent data provenance. If the accession, replicate count or QC criterion
  isn't recoverable from the files, write the open question into that assay's
  `METHODS.md` and ask.
