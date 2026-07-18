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
  │                                annotation: genomic region, gene, repeats, TAD, ChromHMM state
  │                                ⇒ results/rds/annotated_results_...rds   (handoff file)
  ├─▶ MINUTE_2_Hypergeometric_test.R   Enrichment of significant peaks vs TAD / LINE / SINE /
  │                                     LTR (+ subfamilies); separate ChromHMM-state enrichment
  │
  └─▶ MINUTE_3_heatmap.R               Heatmap + k-means clusters + BED export; per-chromosome
        │                              change plots; size×signal×log2FC relationship plots
        │                              ⇒ results/rds/cluster_analysis_inputs.rds
        └─▶ MINUTE_4_cluster_analysis.R   Characterise clusters: signal/loss profiles,
                                           ChromHMM, repeat/region/gene-family composition
```

Stage 1 produces one `.rds` that stages 2 and 3 each read (independently). Stage 4
reads the cluster inputs persisted by stage 3, so it runs after stage 3.

---

## Directory layout

```
ChIP/
  config.R  run_MINUTE.R  MINUTE_1..3.R      # code
  samples.tsv                                 # sample sheet (single source of truth)
  data/
    peaks/             # consensus/master peak BEDs             (committed)
    bigwig/            # *_<genotype>_rep*.mm39.scaled.bw       (download, git-ignored)
    annotation/        # LINE/SINE/LTR .mm39.bed + TAD BEDs     (download, git-ignored)
  results/                                    # ALL generated output (git-ignored except rds/)
    counts/  rds/  tables/  figures/  bed/
  Minute Run/  Peak Calling/                  # cluster (Slurm) steps — see above
```

**Run the R pipeline from the `ChIP/` directory** — repo-relative paths assume it.
There is no `setwd()`. The large inputs default to `ChIP/data/bigwig` and
`ChIP/data/annotation`; set the `MINUTE_DATA` env var to use a copy staged
elsewhere (e.g. an external mount).

---

## Requirements

R (≥ 4.2) with the following packages:

**CRAN:** `data.table`, `dplyr`, `ggplot2`, `ggrepel`, `circlize`

**Bioconductor:** `DESeq2`, `GenomicRanges`, `rtracklayer`, `ChIPseeker`,
`AnnotationDbi`, `org.Mm.eg.db`, `TxDb.Mmusculus.UCSC.mm39.knownGene`,
`ComplexHeatmap`

**GitHub:** [`wigglescout`](https://github.com/cnluzon/wigglescout) — reads
signal from bigWigs. It is **sensitive to the versions of its `future` stack**;
this is the only combination confirmed working for this pipeline:

| Package | Version | Source |
|---------|---------|--------|
| `furrr` | **0.2.3** | CRAN archive (`install_version`) |
| `future` | **1.23.0** | CRAN archive |
| `globals` | **0.14.0** | CRAN archive |
| `wigglescout` | latest | `cnluzon/wigglescout` (GitHub) |

```r
install.packages("remotes")
remotes::install_version("furrr",   version = "0.2.3")
remotes::install_version("future",  version = "1.23.0")
remotes::install_version("globals", version = "0.14.0")
remotes::install_github("cnluzon/wigglescout")
```

**Restart R after installing** so any newer, already-loaded versions unload.
Verify with `packageVersion("furrr")` etc. before running.

> **Symptom of a version drift:** `bw_loci()` fails inside `MINUTE_1`/`MINUTE_3`
> with `Error: values must be length 1, but FUN(X[[i]]) result is length N`. That
> is `wigglescout`'s internal `future`/`furrr` map returning the wrong shape —
> reinstall the pins above (and restart R) rather than editing the R scripts.
> A stray `future::plan("multisession"/"multicore")` set by another package can
> trigger the same error; `future::plan("sequential")` clears it.

---

## Data setup

Small inputs (peak BEDs, `samples.tsv`) are committed and need no setup. The
large raw inputs (bigWigs ~21 GB, RepeatMasker annotation ~185 MB) are hosted on
Google Drive — **download them into the folders below before running.** They
default to `ChIP/data/…`; to keep them elsewhere, download anywhere and point
`MINUTE_DATA` at that folder instead.

| Download into | What | Google Drive folder |
|---------------|------|---------------------|
| `ChIP/data/bigwig/` | 55 scaled mm39 bigWigs (`*_WT_rep*.mm39.scaled.bw`, `*_HP1gKO_rep*…`, ~21 GB) | **[bigWigs ↗](https://drive.google.com/drive/folders/11embg_Ft3tLefSg3Stj6ZAan8u0YDs5k?usp=sharing)** |
| `ChIP/data/annotation/` | `LINE.mm39.bed`, `SINE.mm39.bed`, `LTR.mm39.bed`, `TAD_boundaries_mm39.bed` | **[annotation ↗](https://drive.google.com/drive/folders/1e_HPy6tVKeAZsVEpnQxG-qSZ_VTCYjTh?usp=sharing)** |
| `ChIP/data/peaks/` | master/consensus peak BEDs | already in the repo (committed) |

### Annotation — `gdown` (4 files)

```bash
cd ChIP
pip install gdown
gdown --folder "https://drive.google.com/drive/folders/1e_HPy6tVKeAZsVEpnQxG-qSZ_VTCYjTh" -O data/annotation
```

### bigWigs — `rclone` (55 files)

The bigWig folder holds **55 files, over `gdown --folder`'s 50-file cap** (it
would silently skip 5), so use [`rclone`](https://rclone.org/drive/), which has no
limit. One-time setup: `rclone config` → add a Google Drive remote named `gdrive`.
Then pull the folder by its ID:

```bash
cd ChIP
rclone copy gdrive: data/bigwig \
  --drive-root-folder-id 11embg_Ft3tLefSg3Stj6ZAan8u0YDs5k -P
```

The files must end up **directly** in `data/bigwig/` and `data/annotation/` (no
extra nested folder). Sanity-check:

```bash
ls ChIP/data/bigwig/*.bw | wc -l       # expect 55
ls ChIP/data/annotation/               # LINE/SINE/LTR + TAD_boundaries .bed
```

The bigWig filename expected for every sample is listed in **`samples.tsv`** — if
your downloaded filenames differ, edit that sheet (nothing else).

> **Note:** the `data/bigwig/` and `data/annotation/` folders are git-ignored, so
> the downloaded files are never committed. Only `data/peaks/` is tracked.

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

### Global changes vs. per-peak significance (H4K20me3, H3K9me3)

Per-peak differential testing assumes a **minority** of peaks change against a
stable background — that assumption fails for the broad heterochromatin marks in
this dataset. In the HP1γ knockout, H4K20me3 is lost almost genome-wide:

| mark | median log2FC | % of domains reduced | median baseMean |
|------|--------------:|---------------------:|----------------:|
| H3K4me3 | +0.03 | 40% | 44.8 |
| H3K36me3 | −0.03 | 54% | 2.9 |
| H3K9me2 | −0.10 | 60% | 3.4 |
| **H3K9me3** | **−0.22** | **80%** | **2.6** |
| **H4K20me3** | **−0.47** | **85%** | **1.6** |

For H4K20me3 the **entire distribution is shifted down** (~0.47 log2 units;
median, 25th *and* 75th percentiles all negative). Two things follow:

1. **Per-peak p-values are the wrong lens.** The finding is "everything is down,"
   not "these peaks differ." Each domain is individually only modestly changed and
   noisy, so each is underpowered — and the counts are tiny (baseMean ~1–3), which
   further guts negative-binomial power. Requiring `p < thr` discards the real,
   genome-wide effect. (This is why the earlier exploratory plots used a
   fold-change cutoff.)
2. **`sizeFactors = 1` is what makes the shift visible.** Default DESeq2
   median-of-ratios normalization assumes most peaks are unchanged and would
   re-center the distribution on 0, *erasing* the global loss. Relying on MINUTE's
   external input-scaling (see [Conventions](#how-significant-peaks-are-called))
   preserves it — so the absolute downshift is real, **provided the bigWigs are on
   a common external scale** (they are: input-scaled).

**Consequence for the change plots** (`*_changes_by_chr_*`): the highlight
(coloured) criterion is **per-mark**, set in `plot_criterion` in `MINUTE_3`:

- **`"effect"`** — colour by `|log2FC| > 0.5` (labelled "effect size", *not*
  "significant") for the globally-shifted marks **H3K9me3, H4K20me3**.
- **`"pvalue"`** — `is_significant()` for the centred marks **H3K4me3, H3K36me3,
  H3K9me2**, where per-peak testing is valid.

The diagnostic `figures/log2FC_distribution_by_mark.png` (section 7 of MINUTE_3)
plots each mark's log2FC density with its median — the evidence behind these
assignments. Two caveats: the global-loss claim is best stated as a
**distributional** result (median shift, % of domains down), not a peak count;
and individual **gene labels** on low-count marks (H4K20me3 baseMean ~1.6) are
unreliable — the distribution is the finding, not any single locus.

> This per-mark criterion governs **only the change plots**. The heatmap, BED
> exports, and MINUTE_2 hypergeometric test still use `is_significant()`
> (p-based) throughout.

### Finding: H3K9me3 / H4K20me3 loss by chromatin type and domain size

In the HP1γ knockout, H4K20me3 is lost broadly across all heterochromatin
clusters (median log2FC ≈ −1.0) whereas H3K9me3 loss is milder and
cluster-specific, and the two marks' per-peak losses are weakly *anti*-correlated
(Spearman ρ ≈ −0.30 in the main heterochromatin clusters) — so they largely do
not fall at the same peaks. Where H3K9me3 *does* co-fall with H4K20me3, the
regions are enriched for constitutive-heterochromatin ChromHMM states (`Het`,
`Quies3`), whereas H4K20me3-only loss dominates generic quiescent (`Quies`,
`Quies2`) and distal-intergenic chromatin. This coupled loss is concentrated in
the **largest domains** — size-weighting the enrichment by domain length
strengthens the `Het`/`Quies3` signal (e.g. `Het` log2 ratio 0.38 → 1.33) —
indicating the biggest constitutive-heterochromatin domains are where H3K9me3
and H4K20me3 decline together.

---

## Outputs (all under `results/`)

| File | Produced by | Contents |
|------|-------------|----------|
| `counts/<mark>_bigwig_counts.tsv` | MINUTE_1 | per-peak signal matrix |
| `rds/annotated_results_...rds` | MINUTE_1 | all peaks + DESeq2 stats + annotation (tracked handoff) |
| `tables/enrichment_hypergeometric.tsv` | MINUTE_2 | odds ratios & hypergeometric p-values per repeat/TAD annotation |
| `figures/enrichment_dotplot.pdf` | MINUTE_2 | repeat/TAD enrichment dot plot |
| `tables/enrichment_chromHMM.tsv` | MINUTE_2 | ChromHMM-state enrichment per mark, both per-region (Wilcoxon) and size-weighted (permutation) |
| `figures/enrichment_chromHMM_dotplot.pdf` | MINUTE_2 | ChromHMM enrichment, per-region (each region equal) |
| `figures/enrichment_chromHMM_sizeweighted_dotplot.pdf` | MINUTE_2 | ChromHMM enrichment, size-weighted by domain length (genomic territory) |
| `tables/enrichment_chromHMM_H3K9me3_loss_vs_unchanged.tsv` + 2 pdfs | MINUTE_2 | chromatin states distinguishing H3K9me3-loss vs -unchanged (within H4K20me3-lost); per-region + size-weighted |
| `figures/2000maxgap_indsignificance_with_TAD.pdf` | MINUTE_3 | clustered heatmap of significant peaks (k-means, seeded) |
| `figures/<mark>_changes_by_chr_coloured_by_{genomic_region,repeat}.png` | MINUTE_3 | per-chromosome domain-size vs log2FC, sized by domain size, coloured by region/repeat |
| `figures/log2FC_distribution_by_mark.png` | MINUTE_3 | per-mark log2FC density + median (global-shift diagnostic) |
| `figures/relationship_total_*.png` | MINUTE_3 | genome-wide size×baseMean×log2FC: binned heatmap, scatter, size×signal interaction |
| `figures/relationship_bycluster_*.png` | MINUTE_3 | the same three relationship views, faceted by k-means cluster |
| `figures/relationship_bycluster_change_H3K9me3_vs_H4K20me3_binned.png` | MINUTE_3 | H3K9me3 vs H4K20me3 DESeq2 log2FC plane, binned & coloured by median domain size, by cluster |
| `figures/relationship_bycluster_mark_change_distribution.png` | MINUTE_3 | H3K9me3 / H4K20me3 DESeq2 log2FC distribution across clusters (uniform vs cluster-specific) |
| `figures/relationship_all_regions_H3K9me3_vs_H4K20me3_by_cluster.png` | MINUTE_3 | all merged-set regions (grey) with significant peaks coloured by cluster, in the ΔH3K9me3 × ΔH4K20me3 plane |
| `rds/cluster_analysis_inputs.rds` | MINUTE_3 | cluster table + per-sample signal matrix (handoff to MINUTE_4) |
| `tables/significant_peaks_clusters.tsv` | MINUTE_3 | significant peaks with cluster + annotation (the metadata BED can't hold it) |
| `bed/significant_peaks_<mark>.bed` | MINUTE_3 | significant regions per mark (NCBI seqnames, score = −log10 p) |
| `bed/significant_peaks_with_metadata.bed` | MINUTE_3 | significant regions, score = log2FC (browser BED; cluster/ChIP are in the TSV) |
| `figures/cluster_signal_profile.png` | MINUTE_4 | mean signal per cluster × mark, WT vs KO (what each cluster is) |
| `figures/cluster_loss_heatmap.png` | MINUTE_4 | log2(KO/WT) signal, cluster × mark |
| `figures/cluster_chromHMM_coverage_{abs,zscore}.png` | MINUTE_4 | mean ChromHMM-state coverage per cluster (absolute + per-state z) |
| `tables/cluster_chromHMM_enrichment.tsv` + `cluster_chromHMM_enrichment{,_sizeweighted}.png` | MINUTE_4 | per-cluster ChromHMM enrichment vs other clusters; per-region + size-weighted |
| `figures/cluster_{repeat,region}_composition.png` | MINUTE_4 | repeat-class / genomic-region composition per cluster |
| `figures/cluster_gene_families.png` | MINUTE_4 | Pcdh / KRAB-Zfp / Vmn / Olfr and pericentromeric fractions per cluster |
| `tables/cluster_summary.tsv` | MINUTE_4 | per-cluster n, median size/log2FC, top state/repeat, family fractions |

---

## Notes

- **Seqnames:** counts/DESeq2 use NCBI style (`1`); exported BEDs are NCBI. Set
  `bed_style <- "UCSC"` in MINUTE_3 for `chr`-prefixed output.
- DESeq2 runs with `sizeFactors = 1` — the bigWigs are already input-scaled, so
  no further normalization is applied.
- Significance uses the **raw** `pvalue` (not BH-adjusted `padj`), with lenient
  per-mark thresholds for the broad heterochromatin marks (see table above).
