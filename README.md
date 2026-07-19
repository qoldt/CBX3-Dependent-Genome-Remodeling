# CBX3-Dependent Genome Remodeling

Analysis code for a study of **HP1γ (CBX3)** in organizing chromatin in neurons,
comparing HP1γ knockout vs wild-type across three complementary assays.

## Structure

| Folder | Assay | Status |
|--------|-------|--------|
| [`ChIP/`](ChIP/) | Multiplexed ChIP-seq (MINUTE pipeline) — differential histone-mark analysis over called peaks | Active |
| `scRNAseq/` | Single-cell RNA-seq | Planned |
| `GAM/` | Genome Architecture Mapping | Planned |

Each folder is self-contained, with its own README describing inputs, usage, and
outputs. Genome build is mouse **mm39** throughout.

## ChIP-seq analysis — methods summary

Differential histone-mark ChIP-seq in **HP1γ (Cbx3) conditional-knockout cortex**
(*Emx1-Cre*) vs wild type. Five marks — **H3K4me3, H3K36me3, H3K9me2, H3K9me3,
H4K20me3** — genome **mm39**, 6 WT vs 5 HP1gKO replicates.

- **Signal.** Input-scaled bigWigs from the **MINUTE** multiplexed ChIP-seq
  pipeline; signal is quantified over consensus (MACS3-merged) peaks.
- **Testing.** DESeq2, contrast **HP1gKO vs WT**, with **`sizeFactors = 1`** —
  deliberate, because the bigWigs are already input-scaled and median-of-ratios
  normalization would erase the genome-wide H4K20me3 loss that is the main
  finding. Significance uses the **raw** p-value with lenient per-mark thresholds;
  because the loss is genome-wide (breaking per-peak power), the broad marks are
  also read via an **effect-size** (|log2FC|) criterion.
- **Annotation.** Each peak carries genomic region, nearest gene, RepeatMasker
  class/family, TAD-boundary overlap, and **per-state ChromHMM coverage
  fractions** (18-state; domains span many states, so fractions beat a single
  dominant label).
- **Key analytical choices.** Enrichments are reported both **per-region** and
  **size-weighted** (by domain length); overlap enrichments use **size-matched
  permutation nulls** (kb–Mb heterochromatin domains would otherwise trivially
  overlap narrow features); repeat and gene-family responses are measured by
  **direct `bw_loci` quantification over the elements/exons themselves**
  (peak-calling-independent). External **KAP1/TRIM28** ChIP peaks (ChIP-Atlas)
  were lifted **mm10 → mm39** with CrossMap.
- **Stages** (`ChIP/run_MINUTE.R`): (1) count & annotate → (2) global changes
  (heatmap, clusters, relationships, TAD/ChromHMM enrichment) → (3) differential
  H3K9me3-vs-H4K20me3 loss → (4) repeats → (5) clustered HUSH/CBX3 gene families
  → (6) intersection with silencing machinery (KAP1).

**Headline result:** HP1γ loss causes a **broad, near-genome-wide loss of
H4K20me3** (deepest over IAP/young-L1 and quiescent chromatin) that is **largely
uncoupled from H3K9me3**, derepressing HUSH/CBX3-target gene clusters
(protocadherins, KRAB-ZFPs, vomeronasal/olfactory receptors).

Full details: [`ChIP/README.md`](ChIP/README.md) (pipeline + outputs),
[`ChIP/results/EXECUTIVE_SUMMARY.md`](ChIP/results/EXECUTIVE_SUMMARY.md)
(findings), [`ChIP/results/METHODS.md`](ChIP/results/METHODS.md) (data sources).

## Data

Large raw inputs (bigWigs, RepeatMasker annotation, etc.) are **not** stored in
this repo — see the per-folder README for what to provide locally and where to
get it.
