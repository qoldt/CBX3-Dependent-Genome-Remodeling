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

## Overview

```
config.R                         ← all parameters live here (edit this only)
  │
  ▼
MINUTE_1_Count_and_Annotate.R    Quantify + DESeq2 + annotate
  │                                bigWigs × peaks → counts → DESeq2 → annotated table
  │                                ⇒ annotated_results_...rds   (handoff file)
  ├─▶ MINUTE_2_Hypergeometric_test.R   Enrichment of significant peaks vs
  │                                     TAD boundaries / LINE / SINE / LTR (+ subfamilies)
  │
  └─▶ MINUTE_3_heatmap.R               Heatmap of significant peaks + BED export
                                        ⇒ significant_peaks_<mark>.bed
                                        ⇒ significant_peaks_with_metadata.bed
```

Stage 1 produces one `.rds` that stages 2 and 3 each read. Stages 2 and 3 are
independent — run either, or both, in any order after stage 1.

---

## Requirements

R (≥ 4.2) with the following packages:

**CRAN:** `data.table`, `dplyr`, `ggplot2`, `ggrepel`, `circlize`

**Bioconductor:** `DESeq2`, `GenomicRanges`, `rtracklayer`, `ChIPseeker`,
`AnnotationDbi`, `org.Mm.eg.db`, `TxDb.Mmusculus.UCSC.mm39.knownGene`,
`ComplexHeatmap`

**GitHub:** [`wigglescout`](https://github.com/cnluzon/wigglescout) — used to
read signal from bigWigs. It needs specific pinned dependencies:

```r
install.packages("remotes")
remotes::install_version("furrr",   version = "0.2.3")
remotes::install_version("future",  version = "1.23.0")
remotes::install_version("globals", version = "0.14.0")
remotes::install_github("cnluzon/wigglescout")
```

---

## Data setup

Large inputs are **not** stored in the repo. Place them under the project root:

| Path | What | Source |
|------|------|--------|
| `bigwig/` | scaled mm39 bigWigs (`*_WT_rep*.mm39.scaled.bw`, `*_HP1gKO_rep*...`) | MINUTE pipeline output |
| `LINE.mm39.bed`, `SINE.mm39.bed`, `LTR.mm39.bed` | RepeatMasker annotation | UCSC Table Browser (mm39, rmsk) |
| `peaks/…` | master/consensus peak BEDs | included for the 4 files used; see `config.R$regions` |
| `TAD_boundaries_mm39.bed` | TAD boundaries | included |

Sample naming is defined by `get_bigwig_files()` in `config.R` — adjust it if
your bigWig filenames differ.

---

## Usage

Run the whole pipeline from the repo root:

```bash
Rscript run_MINUTE.R
```

or interactively:

```r
source("run_MINUTE.R")
```

Run a single stage (stage 1 must have produced the `.rds` first):

```bash
Rscript MINUTE_1_Count_and_Annotate.R
Rscript MINUTE_3_heatmap.R
```

---

## Configuration

**All parameters live in `config.R` — edit it there and nowhere else.** It defines:

- Paths (`bigwigDir`, `outputDir`, the handoff `.rds` name)
- `chips` — mark → bigWig sample prefix
- `regions` — mark → input peak BED
- `sig_thresholds` + `is_significant()` — the significance rule

Significance is `|log2FoldChange| > lfc & pvalue < p`, per mark:

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

## Outputs

| File | Produced by | Contents |
|------|-------------|----------|
| `counts/<mark>_bigwig_counts.tsv` | MINUTE_1 | per-peak signal matrix |
| `annotated_results_...rds` | MINUTE_1 | all peaks + DESeq2 stats + annotation |
| enrichment tables / plots | MINUTE_2 | odds ratios & hypergeometric p-values per annotation |
| `2000maxgap_indsignificance_with_TAD.pdf` | MINUTE_3 | clustered heatmap of significant peaks |
| `significant_peaks_<mark>.bed` | MINUTE_3 | significant regions per mark (NCBI seqnames, score = −log10 p) |
| `significant_peaks_with_metadata.bed` | MINUTE_3 | all significant regions + ChIP/cluster metadata |

---

## Notes

- **Seqnames:** counts/DESeq2 use NCBI style (`1`); exported BEDs are NCBI. Set
  `bed_style <- "UCSC"` in MINUTE_3 for `chr`-prefixed output.
- The `*_sig.ucsc.bed` / `all_histone_marks.ucsc.bed` browser tracks (UCSC
  `track` headers, `chr` prefixes) are a **separate downstream conversion**, not
  produced by this workflow.
- DESeq2 runs with `sizeFactors = 1` — the bigWigs are already scaled, so no
  further normalization is applied.
- `archive_superseded/` holds the old monolithic scripts, fully replaced by
  MINUTE_1–3; `MINUTE_4_deeptools-like.R` is unused.
