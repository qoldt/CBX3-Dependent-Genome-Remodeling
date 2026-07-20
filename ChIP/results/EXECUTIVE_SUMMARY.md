# Executive summary — CBX3/HP1γ-dependent genome remodeling

Differential histone-mark ChIP-seq in **HP1γ (Cbx3) conditional-knockout cortex**
(*Emx1-Cre*) vs wild type. MINUTE multiplexed ChIP-seq, genome **mm39**, five
marks (H3K4me3, H3K36me3, H3K9me2, H3K9me3, H4K20me3), 6 WT vs 5 HP1gKO
libraries — **HP1gKO rep2 is excluded by default** as a technical failure
(discordant with rep3, its own technical replicate), so analyses run 6 vs 4.
Replicates are nested (technical replicates of 2 biological replicates per
genotype), so the effective n is 2 vs 2. DESeq2 with `sizeFactors = 1` (bigWigs are input-scaled), contrast
HP1gKO vs WT. Pipeline stages MINUTE_1–6 (see `../README.md`).

## Headline findings

> **Read compartment sizes as median-relative, not absolute.** MINUTE_3's group
> cutoffs are absolute, so against a genome-wide loss they largely ask "did it
> lose?" — to which nearly everything answers yes — and the answer moves whenever
> the global shift moves. Across sample sets (all libraries vs HP1gKO rep2
> excluded) group sizes change by a mean of **84%** under absolute cutoffs but
> only **11%** when each mark is centred on its own median. Crucially the
> ordering also flips: absolute cutoffs make `co_loss` dominate once rep2 is
> dropped (46,921 vs 8,481), whereas median-centred cutoffs keep
> **`H4K20me3_only` > `co_loss` in both sample sets** (1,798 vs 1,452 with all
> libraries; 1,748 vs 1,171 without rep2). **The uncoupling result below
> therefore holds — but only stated as "loses more than its own genome-wide
> average", not "loses in absolute terms."** Both schemes are reported per run in
> `tables/diffloss_group_definitions.tsv`.


1. **H4K20me3 is lost almost genome-wide — the primary effect.** The entire
   distribution shifts down (median log2FC ≈ **−0.68**, ~95% of domains reduced;
   −0.47 / ~85% with HP1gKO rep2 retained),
   deepest over broad heterochromatin. This is a global change, not a
   minority-of-peaks change.

2. **H3K9me3 loss is milder and region-specific** (median ≈ −0.42, ~89% down;
   −0.22 / ~80% with rep2 retained),
   and per peak the two marks' losses are **weakly anti-correlated** (Spearman
   ρ ≈ −0.30 in the main heterochromatin clusters) — the peaks losing the most
   H4K20me3 are *not* the ones losing the most H3K9me3.

3. **Two mechanistically distinct compartments** (MINUTE_3, shared H3K9me3/
   H4K20me3 peak set):
   - **H4K20me3-only loss** (H3K9me3 retained; ~13.7k regions): **large**
     (median 2.3 kb), **distal-intergenic**, **generic-quiescent** chromatin
     (Quies/Quies2 ~82%, little Het), strongly **ERV/IAP-enriched**. H4K20me3
     maintenance here is *uncoupled from H3K9me3*.
   - **Co-loss** (both fall; ~20k regions): **smaller** (1.3 kb), **gene-proximal**,
     **constitutive heterochromatin** (Het/Quies3), including the silenced gene
     families below.

4. **HUSH/CBX3-silenced gene families are derepressed** (MINUTE_5,
   quantified over each gene's silencing-relevant exon). All four families lose
   H3K9me3 **and** H4K20me3 (paired-Wilcoxon FDR ≤ 1e-11):

   | Family | genes | H3K9me3 log2FC (FDR) | H4K20me3 log2FC (FDR) | H3K4me3 |
   |--------|------:|---------------------:|----------------------:|---------|
   | Olfactory (Or/Olfr) | 1107 | −0.39 (2e-150) | −0.54 (1e-180) | −0.03 (0.03) |
   | Vomeronasal (Vmn) | 353 | −0.47 (5e-58) | −0.59 (2e-59) | −0.04 (0.02) |
   | KRAB-ZFP (Zfp) | 428 | −0.32 (6e-50) | −0.45 (1e-63) | ns |
   | Protocadherin (chr18 Pcdha/b/g) | 59 | −0.59 (3e-11) | −0.81 (2e-11) | **+0.31 (up)** |

   *(HP1gKO rep2 excluded — the previous version reported H3K9me3 losses ~0.1–0.15
   shallower, e.g. KRAB-ZFP −0.18; every family deepened.)*

   Clustered protocadherins additionally **gain H3K4me3** — a full derepression
   signature. Tested against a **matched random-gene background** (MINUTE_5;
   background co-loss rate now 41%), all four families co-lose H3K9me3
   **significantly above that gene background** (Fisher BH ≤ 5e-15): protocadherin
   92% (OR 15), vomeronasal 79% (OR 5.2), KRAB-ZFP 68% (OR 3.0), olfactory 62%
   (OR 2.3).

   **The earlier "KRAB-ZFP is H4K20me3-only" conclusion does not survive rep2
   exclusion** — that rested on KRAB-ZFP H3K9me3 being ~−0.18 and co-losing at
   background (32%, ns), both of which were rep2 artefacts. KRAB-ZFP now loses
   H3K9me3 (−0.32) and co-loses above the gene background. It remains
   H4K20me3-**preferential** but by a graded, reference-dependent margin, not a
   binary: relative to the genome-wide *peak* median (−0.415, repeat/intergenic-
   dominated) KRAB-ZFP loses *less* H3K9me3 than typical (+0.09, p 3e-21), whereas
   relative to other *genes* it loses *more*. So: loses both marks, H4K20me3 more
   than H3K9me3, with H3K9me3 loss below the heterochromatin average but above the
   genic average.

5. **Repeats:** IAP/ERVK elements are strongly enriched in the lost compartments
   (IAPEz-int odds ratio ~21 in the H4K20me3-only group vs stable), consistent
   with HUSH/TRIM28 targets; SINE and ERVL-MaLR are depleted. **Direct signal
   quantification over the repeat copies themselves** (MINUTE_4, `bw_loci`,
   peak-calling-independent) confirms this with magnitude: H4K20me3 is lost over
   every class but **deepest over IAP** (IAPEz-int −0.52, IAPLTR1a_Mm −0.55) and
   **young L1** (L1MdT −0.45, L1MdGf −0.41, L1MdA −0.38 log2FC), shallowest over
   **SINE controls (~−0.18)**; H3K9me3 loss is milder (−0.10 to −0.23, ERV-focused);
   and the active/euchromatic marks (H3K4me3, H3K36me3, H3K9me2) stay flat over all
   classes (|log2FC| < 0.1) — a clean specificity control. This ordering
   (IAP > young L1 > older ERV/LINE > SINE) is the mouse read-out of HUSH/TRIM28
   substrate preference and does **not** depend on peak calling or the significance
   definition (unlike the MINUTE_4 hypergeometric test).

6. **Chromatin-state clusters** (MINUTE_2 k-means, 5 groups): an active
   H3K4me3/promoter cluster (which *gains* signal), an H3K36me3/transcription
   cluster, and three heterochromatin flavours that differ in H3K9me2 vs H3K9me3
   content — all of which lose H4K20me3.

7. **KAP1/HUSH occupancy is coupled to the loss selectively** (MINUTE_6; two mm39
   KAP1/TRIM28 tracks — Neuro-2a neural + ChIP-Atlas all-cell — vs a **size-matched
   permutation null**):
   - The loss *compartments* (co-loss / H4K20me3-only) are **not** KAP1-enriched
     beyond what their large domain size predicts (fold ≈ 1); the naive ~1.9×
     overlap is entirely a domain-length artifact.
   - **KRAB-ZFP and protocadherin exons are strongly KAP1-bound** (fold ~2.5× and
     ~4× vs size-matched exons); vomeronasal/olfactory exons are KAP1-**depleted**
     (fold ≈ 0), so their silencing is KAP1-independent.
   - Integrating occupancy with the genotype effect: **KAP1 binding predicts the
     H4K20me3 loss specifically at KRAB-ZFP** (KAP1-bound KRAB-ZFP exons lose
     Δ ≈ −0.3 log2 *more* H4K20me3 than unbound; Wilcoxon p < 1e-4, both tracks),
     and modestly genome-wide (Δ ≈ −0.13). The protocadherin H4K20me3 loss is
     broad and **not** KAP1-graded.

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

## Pipeline (how these results are produced)

Six thematic stages, run by `ChIP/run_MINUTE.R`; each reads the annotated `.rds`
from stage 1 (details in `../README.md`, exact accessions in `METHODS.md`).

1. **Count & annotate** — input-scaled mm39 bigWigs × consensus peaks → DESeq2
   (`sizeFactors = 1`, HP1gKO vs WT) → per-peak stats + annotation (genomic
   region, gene, repeats, TAD, per-state ChromHMM coverage fractions).
2. **Global changes** — clustered heatmap (k-means), per-chromosome change plots,
   per-mark log2FC distributions, mark-vs-mark relationships, cluster
   characterisation, and TAD + ChromHMM enrichment of the significant set.
3. **Differential loss** — split the shared H3K9me3/H4K20me3 peak set into
   co-loss / H4K20me3-only / stable and characterise each.
4. **Repeats** — hypergeometric repeat-class enrichment **and** direct `bw_loci`
   signal change over the repeat copies themselves (peak-calling-independent).
5. **Clustered gene families** — per-gene silencing-exon signal for the HUSH/CBX3
   families + the co-loss-vs-background test.
6. **Intersection** — KAP1/TRIM28 vs the loss: size-matched permutation overlap
   + KAP1-bound-vs-unbound genotype effect (skips if the KAP1 BEDs are absent).

## Where to find things (`results/`)

- **Tables:** `tables/family_exon_summary.tsv`, `tables/family_coloss_vs_background.tsv`,
  `tables/diffloss_*` (ChromHMM / repeat / region enrichment + per-group gene lists),
  `tables/enrichment_*` (incl. young-L1 subfamilies), `tables/kap1_{intersection,genotype_effect}.tsv`,
  `tables/repeat_signal_*` (direct repeat-class change), `tables/cluster_*`,
  `tables/significant_peaks_clusters.tsv`.
- **Figures** (each as `.png` + `.pdf` under `figures/<subanalysis>/`):
  `differential_loss/` (incl. KAP1), `gene_families/`, `repeats/`, `clusters/`,
  `relationships/`, `change_plots/`, `enrichment/`, `heatmap/`.
- **Methods:** `results/METHODS.md` (KAP1 data source + mm10→mm39 liftover; HUSH-target note).
- **deepTools inputs:** `bed/family_*_exons.bed` (target exons, NCBI seqnames
  matching the bigWigs), `bed/significant_peaks_*.bed`.
- **Handoffs:** `rds/annotated_results_*.rds` (all peaks + stats + annotation +
  ChromHMM), `rds/cluster_analysis_inputs.rds`, `rds/family_exon_signal.rds`.
