# MINUTE — Multiplexed ChIP-seq differential analysis (mm39)

Differential histone-mark ChIP-seq analysis for a HP1γ knockout vs wild-type
comparison, built on scaled bigWigs from the MINUTE multiplexed ChIP-seq pipeline.
Five marks: **H3K4me3, H3K9me2, H3K9me3, H4K20me3, H3K36me3**. Genome **mm39**.

## Pipeline

Run stages in order. Each reads its inputs from disk, so they can run
individually (after stage 1) or all at once via the runner.

```
config.R                         # single source of truth for ALL parameters
  │
MINUTE_1_Count_and_Annotate.R    # bigWig signal -> counts -> DESeq2 -> annotate
  │  writes: counts/<mark>_bigwig_counts.tsv
  │  writes: annotated_results_H3K9me3_H4K20me3_2000bp_merged.rds   (the handoff)
  ├── MINUTE_2_Hypergeometric_test.R   # enrichment of significant peaks vs annotations
  │      reads:  annotated_results ...rds
  │      writes: odds-ratio / enrichment tables + plots
  └── MINUTE_3_heatmap.R               # heatmap of significant peaks + BED export
         reads:  annotated_results ...rds  (+ re-reads bigWigs for the signal matrix)
         writes: 2000maxgap_indsignificance_with_TAD.pdf
         writes: significant_peaks_with_metadata.bed
         writes: significant_peaks_<mark>.bed   (NCBI seqnames, score = -log10 p)
```

Run everything sequentially:

```r
# from the repo root (the folder containing config.R)
Rscript run_MINUTE.R      # or, interactively:  source("run_MINUTE.R")
```

MINUTE_2 and MINUTE_3 are independent of each other — both only need the
`.rds` produced by MINUTE_1.

## Where parameters live — edit `config.R` only

`config.R` is sourced by every stage and by `run_MINUTE.R`. It is the **single
source of truth**; stage scripts must not redefine these. It holds:

- Paths: `setwd`, `bigwigDir`, `outputDir`, `annotated_rds`
- `chips` — mark → bigWig sample prefix
- `regions` — mark → input peak BED
- `get_bigwig_files()` — builds the 11 sample paths (6 WT + 5 HP1gKO) per mark
- `loadRepeatBED()` / `load_annotation()` — repeat + TAD annotation
- **`sig_thresholds` + `is_significant()`** — the significance definition
  (`|log2FC| > lfc & pvalue < p`, per mark). MINUTE_2 and MINUTE_3 both call
  `is_significant()`, so their "significant" sets are guaranteed identical.

To change a cutoff, a sample, or a peak file, edit `config.R` and rerun.

## Data dependencies (NOT in git — provide locally)

These are too large to commit; place them under the repo root before running:

- `bigwig/` — scaled mm39 bigWigs from the MINUTE pipeline (~21 GB)
- `LINE.mm39.bed`, `SINE.mm39.bed`, `LTR.mm39.bed` — UCSC RepeatMasker mm39
  (~185 MB total)
- `TAD_boundaries_mm39.bed` — small, committed
- `peaks/…` — master/consensus peak BEDs referenced in `config.R$regions`
  (the 4 used files are small and committed; the rest of `peaks/` is not)

## Conventions & gotchas

- **Seqnames:** counts/DESeq2 use **NCBI** style (`1`, `2`, …); annotation
  overlaps switch to **UCSC** (`chr1`) as needed. Exported
  `significant_peaks_<mark>.bed` are **NCBI**. Flip `bed_style <- "UCSC"` in
  MINUTE_3 for `chr`-prefixed output.
- The `*_sig.ucsc.bed` / `all_histone_marks.ucsc.bed` files (chr-prefixed, with
  UCSC `track` headers) are a **separate downstream conversion**, not produced
  by this pipeline.
- DESeq2 uses `sizeFactors = 1` on purpose — the bigWigs are already scaled, so
  no additional normalization is applied.
- `annotated_results` (the `.rds`) contains **all** measured peaks with stats +
  annotation; significance is applied downstream, per stage, via
  `is_significant()`.

## Not part of the pipeline

- `archive_superseded/` — old monolithic `Minute_count_and_annotate_from_scaled_bw*.R`,
  fully replaced by MINUTE_1–3. Kept for reference only.
- `MINUTE_4_deeptools-like.R` — not used.
- Upstream helpers (`MergePeaks_to_masterPeak.R`, `Generate_Countmats.R`,
  `deseq2_minute.R`) — optional; produce/precompute inputs, not run by `run_MINUTE.R`.
