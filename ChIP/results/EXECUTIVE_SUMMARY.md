# Executive summary — CBX3/HP1γ-dependent genome remodeling

Differential histone-mark ChIP-seq in **HP1γ (Cbx3) conditional-knockout cortex**
(*Emx1-Cre*) vs wild type. MINUTE multiplexed ChIP-seq, genome **mm39**, five
marks (H3K4me3, H3K36me3, H3K9me2, H3K9me3, H4K20me3), 6 WT vs 5 HP1gKO
replicates. DESeq2 with `sizeFactors = 1` (bigWigs are input-scaled), contrast
HP1gKO vs WT. Pipeline stages MINUTE_1–5 (see `../README.md`).

## Headline findings

1. **H4K20me3 is lost almost genome-wide — the primary effect.** The entire
   distribution shifts down (median log2FC ≈ **−0.47**, ~85% of domains reduced),
   deepest over broad heterochromatin. This is a global change, not a
   minority-of-peaks change.

2. **H3K9me3 loss is milder and region-specific** (median ≈ −0.22, ~80% down),
   and per peak the two marks' losses are **weakly anti-correlated** (Spearman
   ρ ≈ −0.30 in the main heterochromatin clusters) — the peaks losing the most
   H4K20me3 are *not* the ones losing the most H3K9me3.

3. **Two mechanistically distinct compartments** (MINUTE_5, shared H3K9me3/
   H4K20me3 peak set):
   - **H4K20me3-only loss** (H3K9me3 retained; ~13.7k regions): **large**
     (median 2.3 kb), **distal-intergenic**, **generic-quiescent** chromatin
     (Quies/Quies2 ~82%, little Het), strongly **ERV/IAP-enriched**. H4K20me3
     maintenance here is *uncoupled from H3K9me3*.
   - **Co-loss** (both fall; ~20k regions): **smaller** (1.3 kb), **gene-proximal**,
     **constitutive heterochromatin** (Het/Quies3), including the silenced gene
     families below.

4. **HUSH/CBX3-silenced gene families are derepressed** (MINUTE_4 part 7,
   quantified over each gene's silencing-relevant exon). All four families lose
   H3K9me3 **and** H4K20me3 (paired-Wilcoxon FDR ≤ 1e-11):

   | Family | genes | H3K9me3 log2FC (FDR) | H4K20me3 log2FC (FDR) | H3K4me3 |
   |--------|------:|---------------------:|----------------------:|---------|
   | Olfactory (Or/Olfr) | 1107 | −0.28 (5e-117) | −0.38 (6e-173) | ns |
   | Vomeronasal (Vmn) | 353 | −0.35 (2e-55) | −0.49 (3e-58) | ns |
   | KRAB-ZFP (Zfp) | 428 | −0.18 (1e-34) | −0.32 (1e-49) | ns |
   | Protocadherin (chr18 Pcdha/b/g) | 59 | −0.42 (4e-11) | −0.65 (2e-11) | **+0.36 (up)** |

   Clustered protocadherins additionally **gain H3K4me3** — a full derepression
   signature.

5. **Repeats:** IAP/ERVK elements are strongly enriched in the lost compartments
   (IAPEz-int odds ratio ~21 in the H4K20me3-only group vs stable), consistent
   with HUSH/TRIM28 targets; SINE and ERVL-MaLR are depleted.

6. **Chromatin-state clusters** (MINUTE_3 k-means, 5 groups): an active
   H3K4me3/promoter cluster (which *gains* signal), an H3K36me3/transcription
   cluster, and three heterochromatin flavours that differ in H3K9me2 vs H3K9me3
   content — all of which lose H4K20me3.

## Interpretation

CBX3/HP1γ is required to **maintain H4K20me3 broadly** across quiescent and
ERV-rich chromatin, **largely independently of H3K9me3**. The canonical model
places H4K20me3 (SUV420H1/2) downstream of H3K9me3/HP1; here H4K20me3 falls even
where H3K9me3 is retained, pointing to a more direct HP1γ requirement. Loss of
this maintenance derepresses HUSH-target gene families (KRAB-ZFPs, clustered
protocadherins, vomeronasal/olfactory receptors).

## Statistical framing (read before quoting numbers)

- **Global change breaks per-peak testing.** For H4K20me3/H3K9me3 the effect is
  genome-wide, so per-peak p-values are underpowered; the change plots use an
  **effect-size** criterion (|log2FC|) for these marks and **p-value** for the
  centred marks (H3K4me3/H3K36me3/H3K9me2). `sizeFactors = 1` is deliberate — it
  preserves the global shift that median-of-ratios normalization would erase.
- **Enrichment reported two ways.** ChromHMM/repeat enrichments are given both
  **per-region** (each region equal) and **size-weighted** (by domain length =
  genomic territory); they can diverge, and the divergence is informative
  (e.g. Het co-loss strengthens under size-weighting — driven by large domains).

## QC caveat

**HP1gKO replicate 2** has globally elevated signal across *all five* marks
(within-mark z up to +2.5) — a scaling/depth residual (uniform across marks =
technical, not biological), which `sizeFactors = 1` carries through. Because it
inflates the KO mean, it makes the measured losses **conservative** (likely
underestimated). Recommended: a rep2-excluded sensitivity check. Do **not**
switch to median-of-ratios to "correct" it.

## Where to find things (`results/`)

- **Tables:** `tables/family_exon_summary.tsv`, `tables/diffloss_*` (ChromHMM /
  repeat / region enrichment + per-group gene lists), `tables/enrichment_*`,
  `tables/cluster_*`, `tables/significant_peaks_clusters.tsv`.
- **Figures** (each as `.png` + `.pdf` under `figures/<subanalysis>/`):
  `differential_loss/`, `gene_families/`, `clusters/`, `relationships/`,
  `change_plots/`, `enrichment/`, `heatmap/`.
- **deepTools inputs:** `bed/family_*_exons.bed` (target exons, NCBI seqnames
  matching the bigWigs), `bed/significant_peaks_*.bed`.
- **Handoffs:** `rds/annotated_results_*.rds` (all peaks + stats + annotation +
  ChromHMM), `rds/cluster_analysis_inputs.rds`, `rds/family_exon_signal.rds`.
