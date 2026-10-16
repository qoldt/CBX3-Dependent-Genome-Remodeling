# CBX3-Dependent Genome Remodeling

Analysis code for a study of **HP1γ (CBX3)** in organizing chromatin in neurons,
comparing HP1γ knockout vs wild-type across three complementary assays.

## Structure

| Folder | Assay | Status |
|--------|-------|--------|
| [`ChIP/`](ChIP/) | Multiplexed ChIP-seq (MINUTE pipeline) — differential histone-mark analysis over called peaks | Active |
| `scRNAseq/` | Single-cell RNA-seq | Planned |
| [`GAM/`](GAM/) | Genome Architecture Mapping — 20 kb NPMI contact maps, differential contacts, integrated with the ChIP results | Active |

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

## GAM analysis — methods summary

Genome Architecture Mapping of **HP1γ (Cbx3) conditional-knockout vs wild-type**
nuclei, genome **mm39**, **20 kb** bins. Two biological replicates per genotype,
merged upstream, so the effective n is 2 vs 2 — matching the ChIP side.

- **Inputs.** Dense per-chromosome **NPMI** matrices per genotype (autosomes +
  chrX/chrY), plus delivered **differential contact** pairs and per-bin
  **hotspot** counts (autosomes only; z > 2, NPMI > 0.1, counts within 0-5 Mb).
  Genotype is encoded only by the folder name.
- **Direction.** Contacts stronger in WT are **lost** in the KO; those stronger in
  cKO are **gained**. All downstream columns are named from the KO's point of view.
- **Integration.** ChIP signal and DESeq2 statistics are projected onto the same
  20 kb autosomal grid, with the ChIP compartments (median-centred labels
  primary), repeat content, ChromHMM state and TAD boundaries per bin. PC1
  compartments, insulation and contact decay come from the NPMI matrices, with PC1
  oriented per chromosome to the H3K4me3-high side.
- **Key analytical choices.** The delivered files carry no replicate dimension, so
  no p-value is attached to a contact; GAM-side statistics use size- or
  composition-matched permutation nulls on the bin side. Matrices are processed
  one chromosome at a time.
- **Stages** (`GAM/run_GAM.R`): (1) bins, ChIP signal and annotation -> (2) NPMI
  compartments/insulation/decay -> (3) contact coupling, (4) compartments vs
  contacts, (5) HUSH targets and repeats vs contacts; stages 3-5 read the outputs
  of 1 and 2.

No GAM results have been generated yet. Full details:
[`GAM/README.md`](GAM/README.md) (data and pipeline),
[`GAM/results/METHODS.md`](GAM/results/METHODS.md) (data sources and open
provenance questions).

## Data

Large raw inputs (bigWigs, RepeatMasker annotation, etc.) are **not** stored in
this repo — see the per-folder README for what to provide locally and where to
get it.
