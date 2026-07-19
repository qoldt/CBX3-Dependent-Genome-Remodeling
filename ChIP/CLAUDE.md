# MINUTE — Multiplexed ChIP-seq differential analysis (mm39)

Differential histone-mark ChIP-seq analysis for a HP1γ knockout vs wild-type
comparison, built on scaled bigWigs from the MINUTE multiplexed ChIP-seq pipeline.
Five marks: **H3K4me3, H3K9me2, H3K9me3, H4K20me3, H3K36me3**. Genome **mm39**.

## Two upstream cluster steps (not run by the R runner)

- `Minute Run/` — MINUTE Snakemake pipeline on Slurm → input-scaled bigWigs.
- `Peak Calling/` — MACS3 + `MergePeaks_to_masterPeak.R` → consensus peak BEDs
  (`data/peaks/`).

The R pipeline **starts from** those bigWigs + peaks + repeat/TAD annotation.

## Pipeline (R)

Run from the **`ChIP/` directory** (no `setwd()`; paths are repo-relative).
Each stage reads its inputs from disk, so they can run individually (after
stage 1) or all at once via the runner.

Stages are organised **thematically**. MINUTE_1 writes the annotated `.rds`
handoff; every later stage reads it independently. Stage 6 also reads the
family-exon signal persisted by **stage 5**.

```
config.R + samples.tsv                  # parameters + sample sheet (single sources of truth)
  │
MINUTE_1_Count_and_Annotate.R           # bigWig signal -> counts -> DESeq2 -> annotate
  │  writes: results/rds/annotated_results_..._2000bp_merged.rds  (handoff)
  ├── MINUTE_2_Global_Changes.R         # heatmap + k-means clusters + cluster characterisation;
  │      per-chr change plots + log2FC distributions; mark-vs-mark relationships;
  │      TAD + ChromHMM enrichment. Re-reads bigWigs for the signal matrix.
  │      writes: rds/cluster_analysis_inputs.rds (self-consumed), bed/significant_peaks_*
  ├── MINUTE_3_Differential_loss.R      # co-loss / H4K20me3-only / stable compartments
  ├── MINUTE_4_Repeats.R                # (A) repeat hypergeometric enrichment
  │      (B) direct bw_loci signal change over repeat copies (aggregate + per-copy)
  ├── MINUTE_5_Clustered_Gene_Families.R# HUSH/CBX3 gene clusters; co-loss-vs-background
  │      writes: rds/family_exon_signal.rds (consumed by MINUTE_6), bed/family_*_exons.bed
  └── MINUTE_6_Intersection.R           # KAP1/TRIM28 vs the loss (size-matched null +
         KAP1-bound-vs-unbound genotype effect); skips if the KAP1 BEDs are absent
```

Run everything sequentially:

```r
# from the ChIP/ directory
Rscript run_MINUTE.R        # or interactively: source("run_MINUTE.R")
# store elsewhere?  MINUTE_DATA=/path Rscript run_MINUTE.R
```

Stages 2–5 only need the `.rds` from MINUTE_1; stage 6 needs stage 5's output.

## Where parameters live — edit `config.R` / `samples.tsv` only

`config.R` is sourced by every stage and by `run_MINUTE.R`. It is the **single
source of truth** for parameters; stage scripts must not redefine them. It holds:

- Roots: `MINUTE_DATA` env var → `bigwig_dir` / `annot_dir`; repo-relative
  `peaks_dir`, `results_dir` (+ `counts_dir`/`rds_dir`/`tables_dir`/`fig_dir`/`bed_dir`)
- Loads **`samples.tsv`** into `samples` (sample_id, mark, genotype, replicate,
  bigwig) and derives `marks` + `get_bigwig_files(mark)`
- `regions` — mark → input peak BED
- `loadRepeatBED()` / `load_annotation()` — repeat + TAD annotation
- **`sig_thresholds` + `is_significant()`** — the significance definition
  (`|log2FC| > lfc & pvalue < p`, per mark), shared by every stage that needs a
  "significant" set, so those sets are guaranteed identical.
- Shared plotting constants + helpers: `geno_colors`, `theme_m`, `mark_levels`,
  and `hyper_row()` (one hypergeometric enrichment row) — used across stages.

Genotype/replicate come from `samples.tsv` (matched by sample_id), NOT from
positional `rep()` calls — the DESeq2 design and heatmap annotation cannot drift
from the actual bigWig set. To change a cutoff/peak file, edit `config.R`; to
add/remove/rename a sample, edit `samples.tsv`.

## Data dependencies (NOT in git — download per the README)

bigWigs + annotation are hosted on Google Drive and download into repo-relative
`ChIP/data/bigwig` and `ChIP/data/annotation` (the `MINUTE_DATA` default is now
`data`). Override with the env var to use a copy staged elsewhere.

- `data/bigwig/` — scaled mm39 bigWigs from the MINUTE pipeline (~21 GB)
- `data/annotation/{LINE,SINE,LTR}.mm39.bed` — UCSC RepeatMasker mm39 (~185 MB)
- `data/annotation/TAD_boundaries_mm39.bed` — TAD boundary calls (small)
- `data/peaks/…` — master/consensus peak BEDs referenced in `config.R$regions`
  (small, committed)

## Conventions & gotchas

- **Seqnames:** counts/DESeq2 use **NCBI** style (`1`, `2`, …); annotation
  overlaps switch to **UCSC** (`chr1`) as needed. Exported
  `significant_peaks_<mark>.bed` are **NCBI**. Flip `bed_style <- "UCSC"` in
  MINUTE_2 for `chr`-prefixed output.
- DESeq2 uses `sizeFactors = 1` on purpose — the bigWigs are already scaled.
  Contrast is `HP1gKO` vs `WT` (WT is the reference level).
- Significance uses the **raw** `pvalue`, not BH-adjusted `padj`.
- `annotated_results` (the `.rds`) contains **all** measured peaks with stats +
  annotation; significance is applied downstream, per stage, via `is_significant()`.
- `results/` is git-ignored **except** `results/rds/` (the handoff `.rds` is tracked).

## Not part of the pipeline

- Upstream helpers (`MergePeaks_to_masterPeak.R`) — produce/precompute inputs,
  not run by `run_MINUTE.R`.
