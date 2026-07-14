# MINUTE — Multiplexed ChIP-seq Differential Analysis

Differential histone-mark ChIP-seq analysis comparing **HP1γ knockout vs
wild-type**, built on scaled bigWigs from the
[MINUTE](https://github.com/NBISweden/minute) multiplexed ChIP-seq pipeline.

- **Marks:** H3K4me3, H3K9me2, H3K9me3, H4K20me3, H3K36me3
- **Genome:** mouse mm39
- **Design:** 6 WT vs 5 HP1gKO replicates, tested with DESeq2

The workflow quantifies ChIP signal over called peaks, tests each mark for
genotype-dependent changes, annotates the results (genomic region, gene,
repeats, TAD boundaries), and reports significant regions as heatmaps,
enrichment statistics, and BED tracks.

---

## Where things come from — the two upstream steps

The R analysis **starts from bigWigs + peaks + repeat/TAD annotation**. Those
are produced by two cluster steps that live in this folder but are run
separately on the HPC (Slurm), not by the R runner:

1. **`Minute Run/`** — the [MINUTE](https://github.com/NBISweden/minute)
   Snakemake pipeline. Demultiplexes, aligns to mm39, and emits the
   **input-scaled bigWigs** (`*_<genotype>_rep*.mm39.scaled.bw`). Submit with
   `Minute Run/submit_skeleton.sh`.
2. **`Peak Calling/`** — MACS3 + `MergePeaks_to_masterPeak.R`, producing the
   consensus/master **peak BEDs** in `data/peaks/`.

After those finish on the cluster, stage the outputs locally (see **Data setup**)
and run the R pipeline below.

---

## Overview

```
config.R                         ← all parameters + the sample sheet load (edit this only)
samples.tsv                      ← one row per bigWig: sample_id, mark, genotype, replicate
  │
  ▼
MINUTE_1_Count_and_Annotate.R    Quantify + DESeq2 + annotate
  │                                bigWigs × peaks → counts → DESeq2 → annotated table
  │                                ⇒ results/rds/annotated_results_...rds   (handoff file)
  ├─▶ MINUTE_2_Hypergeometric_test.R   Enrichment of significant peaks vs
  │                                     TAD boundaries / LINE / SINE / LTR (+ subfamilies)
  │
  └─▶ MINUTE_3_heatmap.R               Heatmap of significant peaks + BED export
```

Stage 1 produces one `.rds` that stages 2 and 3 each read. Stages 2 and 3 are
independent — run either, or both, in any order after stage 1.

---

## Directory layout

```
ChIP/
  config.R  run_MINUTE.R  MINUTE_1..3.R      # code
  samples.tsv                                 # sample sheet (single source of truth)
  data/peaks/                                 # consensus/master peak BEDs (committed)
  results/                                    # ALL generated output (git-ignored except rds/)
    counts/  rds/  tables/  figures/  bed/
  Minute Run/  Peak Calling/                  # cluster (Slurm) steps — see above

$MINUTE_DATA/          # large raw inputs, OUTSIDE the repo (default ~/SynologyDrive/MINUTE)
  bigwig/              # *_<genotype>_rep*.mm39.scaled.bw
  annotation/          # LINE.mm39.bed  SINE.mm39.bed  LTR.mm39.bed  TAD_boundaries_mm39.bed
```

**Run the R pipeline from the `ChIP/` directory** — repo-relative paths assume it.
There is no `setwd()`; the only machine-specific path is `$MINUTE_DATA`.

---

## Requirements

R (≥ 4.2) with the following packages:

**CRAN:** `data.table`, `dplyr`, `ggplot2`, `ggrepel`, `circlize`

**Bioconductor:** `DESeq2`, `GenomicRanges`, `rtracklayer`, `ChIPseeker`,
`AnnotationDbi`, `org.Mm.eg.db`, `TxDb.Mmusculus.UCSC.mm39.knownGene`,
`ComplexHeatmap`

**GitHub:** [`wigglescout`](https://github.com/cnluzon/wigglescout) — reads
signal from bigWigs. It needs specific pinned dependencies:

```r
install.packages("remotes")
remotes::install_version("furrr",   version = "0.2.3")
remotes::install_version("future",  version = "1.23.0")
remotes::install_version("globals", version = "0.14.0")
remotes::install_github("cnluzon/wigglescout")
```

---

## Data setup

Small inputs (peak BEDs, `samples.tsv`) are committed. Large raw inputs are
**not** — place them under `$MINUTE_DATA` (default `~/SynologyDrive/MINUTE`,
override with the env var):

| Path | What | Source |
|------|------|--------|
| `$MINUTE_DATA/bigwig/` | scaled mm39 bigWigs (`*_WT_rep*.mm39.scaled.bw`, `*_HP1gKO_rep*…`) | `Minute Run/` pipeline output |
| `$MINUTE_DATA/annotation/LINE.mm39.bed`, `SINE.mm39.bed`, `LTR.mm39.bed` | RepeatMasker annotation | UCSC Table Browser (mm39, rmsk) |
| `$MINUTE_DATA/annotation/TAD_boundaries_mm39.bed` | TAD boundaries | Hi-C boundary calls (mm39) |
| `data/peaks/…` | master/consensus peak BEDs | `Peak Calling/` output (committed) |

The bigWig filename for every sample is listed in **`samples.tsv`** — if your
filenames differ, edit that sheet (nothing else).

---

## Usage

```bash
cd ChIP
Rscript run_MINUTE.R                      # all three stages
# or, if the store is elsewhere:
MINUTE_DATA=/path/to/store Rscript run_MINUTE.R
```

Single stage (stage 1 must have produced the `.rds` first):

```bash
cd ChIP
Rscript MINUTE_1_Count_and_Annotate.R
Rscript MINUTE_3_heatmap.R
```

---

## Configuration

**All parameters live in `config.R`; all samples live in `samples.tsv`.** Between
them they define:

- Roots/paths (`$MINUTE_DATA`, `data/peaks`, `results/…`, the handoff `.rds`)
- `samples.tsv` — mark, genotype, replicate, bigWig filename per sample
- `regions` — mark → input peak BED
- `sig_thresholds` + `is_significant()` — the significance rule

Genotype and replicate are read from the sample sheet by column name, so the
DESeq2 design and the heatmap annotation can never drift from the actual bigWig
set. Significance is `|log2FoldChange| > lfc & pvalue < p`, per mark:

| Mark | log2FC | p-value |
|------|--------|---------|
| H3K4me3 | > 0.5 | < 0.05 |
| H3K9me2 | > 0.5 | < 0.05 |
| H3K9me3 | > 0.5 | < 0.10 |
| H4K20me3 | > 0.5 | < 0.20 |
| H3K36me3 | > 0.5 | < 0.20 |

Both MINUTE_2 and MINUTE_3 call `is_significant()`, so their "significant" sets
are always identical.

---

## How significant peaks are called

Per mark, `MINUTE_1`:

1. **Regions** — the consensus/master peak BED for that mark (from `Peak
   Calling/`). H3K9me3 and H4K20me3 share one **2 kb-merged** peak set
   (`2000bp_merge_H3K9me3_H4K20me3.bed`).
2. **Quantify** — `wigglescout::bw_loci()` reads each of the 11 scaled bigWigs
   (6 WT + 5 HP1gKO) over every peak, giving a peaks × samples matrix. `bw_loci`
   returns the **mean coverage** per interval (i.e. length-normalized), which is
   then `round()`-ed to integer pseudo-counts for DESeq2.
3. **Filter** — drop peaks with `rowSums(counts) ≤ 10`.
4. **Test** — DESeq2 with `sizeFactors = 1` (no renormalization — the bigWigs are
   already input-scaled), design `~condition`, Wald test, contrast **HP1gKO vs
   WT** → a `log2FoldChange` and raw `pvalue` per peak.
5. **Call significance** — `is_significant()`: `|log2FoldChange| > 0.5` **and**
   raw `pvalue <` the per-mark threshold (see table above). The same rule is used
   by MINUTE_2 and MINUTE_3.

Two deliberate but "soft" choices to keep in mind: significance uses the **raw**
`pvalue`, not BH-adjusted `padj`, and the broad heterochromatin marks
(H4K20me3, H3K36me3) get **lenient** `p < 0.20` cutoffs.

### How peak length relates to significance

Length is **not a term in the test**, but it shapes the outcome three ways:

1. **Signal is a mean, so length doesn't inflate the count directly** — a 20 kb
   domain and a 500 bp peak with the same average coverage get the same value.
   Longer intervals do average over more bins, giving a **more stable estimate
   (lower within-group variance)** → DESeq2 estimates lower dispersion → more
   power. So longer peaks are, all else equal, *easier* to call significant — not
   because of bigger numbers, but because of steadier ones.
2. **Peak length is confounded with the mark.** H3K9me3/H4K20me3 use the broad
   2 kb-merged domain set, while H3K4me3 is narrow promoter peaks — so "length"
   and "which mark / which p-threshold" move together. Don't read a
   length↔significance trend as biology without accounting for that.
3. **Length is carried downstream, not into the call.** MINUTE_3 annotates each
   significant peak with `peak_size_kb` (capped at 2 kb for heatmap color), and
   ChIPseeker tends to label very large domains as "Promoter" (a length artifact,
   flagged in a MINUTE_1 comment). Length affects interpretation, not the flag.

> **Caveat:** points above assume `bw_loci` aggregates by **mean**. If the
> quantification were ever switched to summed/total coverage, length *would*
> scale the counts directly and bias long peaks toward significance.

---

## Outputs (all under `results/`)

| File | Produced by | Contents |
|------|-------------|----------|
| `counts/<mark>_bigwig_counts.tsv` | MINUTE_1 | per-peak signal matrix |
| `rds/annotated_results_...rds` | MINUTE_1 | all peaks + DESeq2 stats + annotation (tracked handoff) |
| `tables/enrichment_hypergeometric.tsv` | MINUTE_2 | odds ratios & hypergeometric p-values per annotation |
| `figures/enrichment_dotplot.pdf` | MINUTE_2 | enrichment dot plot |
| `figures/2000maxgap_indsignificance_with_TAD.pdf` | MINUTE_3 | clustered heatmap of significant peaks |
| `bed/significant_peaks_<mark>.bed` | MINUTE_3 | significant regions per mark (NCBI seqnames, score = −log10 p) |
| `bed/significant_peaks_with_metadata.bed` | MINUTE_3 | all significant regions + ChIP/cluster metadata |

---

## Notes

- **Seqnames:** counts/DESeq2 use NCBI style (`1`); exported BEDs are NCBI. Set
  `bed_style <- "UCSC"` in MINUTE_3 for `chr`-prefixed output.
- DESeq2 runs with `sizeFactors = 1` — the bigWigs are already input-scaled, so
  no further normalization is applied.
- Significance uses the **raw** `pvalue` (not BH-adjusted `padj`), with lenient
  per-mark thresholds for the broad heterochromatin marks (see table above).
