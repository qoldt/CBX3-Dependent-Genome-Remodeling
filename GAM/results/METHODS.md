# Methods notes — GAM

## Data sources and processing parameters

Genome Architecture Mapping (GAM) data for **HP1γ (Cbx3) conditional knockout
(cKO) vs wild type (WT)**, aligned to **mm39** and binned at **20 kb**. The
delivered products (about 7.7 GB, git-ignored) are:

- **NPMI matrices**, one dense square matrix per chromosome per genotype (chr1-19,
  chrX, chrY), `HP1g_mm39.segregation_at20000.2026-02-23.passed_qc.curated.sigma.2.6.NPMI.<chr>.squared.txt.gz`.
  Parameters recoverable from the filename: segregation table at 20 kb, dated
  2026-02-23, restricted to `passed_qc.curated` profiles, `sigma.2.6`. Genotype is
  encoded only by the directory (`WT_at20Kb` / `cKO_at20Kb`); the filenames are
  identical between the two.
- **Differential contacts** (`Diff_contacts/`, 38 `.bg2.gz`, autosomes only),
  WT vs cKO contact pairs called significant by the upstream caller with
  **z > 2**, **NPMI > 0.1** (in the selecting genotype) and no distance cap
  (`maxdist = inf`). Observed contact distances span 20 kb to 57.18 Mb. A
  contact is labelled by the genotype in which it is stronger: `selWT`
  (25,558,790 pairs; lost in the KO) and `selcKO` (25,042,198 pairs; gained).
- **Hotspots** (`Hotspots/`, 38 `.bg.gz`, autosomes only): per 20 kb bin, the
  number of significant differential contacts within 0-5 Mb. About 114,793
  (selWT) and 114,651 (selcKO) bins carry a count (range 1-162) out of 123,133
  autosomal bins, so these are dense continuous tracks rather than hit lists.

**Replication.** Two biological replicates per genotype, merged upstream into a
single segregation table per genotype (confirmed by the data owner). The effective
biological n is 2 vs 2, the same as the ChIP side, where 6 WT vs 4 HP1gKO
libraries are nested technical replicates of 2 biological replicates per
genotype. The delivered files have no replicate dimension, so the z > 2 filter is
the only significance carried by the inputs; downstream statistics use
matched-background permutation on the bin side (size-matched nulls via
`size_perm()`), not a p-value on a contact.

**Integration with ChIP.** ChIP signal and DESeq2 statistics (from the MINUTE
pipeline, `../ChIP/`) are projected onto the 20 kb autosomal GAM grid. Seqnames
are flipped between NCBI (ChIP peaks) and UCSC (GAM) with `to_ucsc()`/`to_ncbi()`.
Differential contact direction is expressed from the KO's point of view (`selWT`
= lost, `selcKO` = gained). Chromatin compartments from ChIP use the MINUTE_3
cutoffs, with the median-centred variant as primary. PC1 compartments from NPMI
are oriented per chromosome so that A is the H3K4me3-high side.

No GAM results have been generated at the time of writing.

---

## OPEN QUESTIONS / TODO — provenance not recoverable from the files

These are **not** known and must not be assumed. Each needs an answer from the
data owner before this file can be treated as complete.

1. **How many nuclear profiles (tubes) per genotype** went into the WT and cKO
   segregation tables, before and after QC? Depth per genotype affects NPMI
   noise and could bias a WT-vs-cKO comparison if unequal.
2. **What are the QC criteria behind `passed_qc.curated`?** Which metrics and
   thresholds, and what does "curated" add (manual review?).
3. **What was sigma 2.6 fit on?** Which data or procedure set the value, and was
   it fit per genotype or jointly?
4. **Which GAM protocol variant was used, and what is the reference?** (For
   example standard GAM, multiplex GAM, or another; slice thickness; number of
   nuclear profiles per tube.) Citation to be added.
5. **Is the cKO the same Emx1-Cre Cbx3 conditional-knockout cortex as the ChIP
   data?** Same tissue, age, sex and genotype definition, or a different cohort?
6. **Which tool produced the NPMI matrices and differential contacts**, and what
   does the caller assume (distributions, multiple-testing handling, how `deltaz`
   is defined)? The upstream version and command are not recorded in the files.
7. **Are the two merged biological replicates balanced** in profile count, and
   were they collected in the same batch? Pooling unequal replicates would let
   one replicate dominate a genotype.
8. **The segregation tables themselves — windows x nuclear profiles, per
   genotype.** This is now the top data request, because it is the only thing
   that unlocks GAM's defining capability and we do not have it.

   GAM detects **three-way (and higher) contacts** directly, because a nuclear
   profile either contains a window or does not and no ligation step forces loci
   into pairs. Beagrie et al. 2017 (Nature 543:519) used exactly this to report
   an abundance of three-way contacts genome-wide, enriched between highly
   transcribed regions and super-enhancers, modelled with SLICE. Nothing in our
   delivery supports that analysis:

   - What we hold is **pairwise** — NPMI matrices plus differential contacts and
     hotspots derived from them.
   - Three-way co-segregation **cannot be recovered from pairwise NPMI**, as a
     matter of arithmetic rather than effort: pairwise marginals do not determine
     the probability that all three loci occupy the same slice. Three loci that
     pairwise-contact strongly may never co-occur in one profile.
   - Per GAMtools' own documentation the segregation table is "the crucial and
     most important output of a GAM experiment - all further downstream analysis
     will generally be based on" it; `gamtools matrix` then derives the pairwise
     matrices we were sent (cosegregation / linkage / dprime / npmi). So the
     table exists upstream by construction, and our filenames name it:
     `HP1g_mm39.segregation_at20000.<date>.passed_qc.curated.sigma.2.6.NPMI.chrN.squared.txt.gz`.

   **Ask for:** `segregation_at_20000.table` (and a coarser one, 50 kb or 100 kb,
   for triplet power) for each genotype, unmerged per replicate if possible,
   which would also settle questions 1 and 7. **Also ask** how many NPs passed
   QC per genotype: triplet frequencies scale as the cube of detection
   probability, so power here depends on NP count far more steeply than the
   pairwise analysis does, and a triplet null needs the NP-level data to permute.

   **What could then be asked** (all currently blocked): do heterochromatin
   domains lose multi-way contacts while keeping pairwise ones, and does that
   depend on domain size; are multi-way hubs enhancer-centred here as they are in
   mESC (we have Enh/EnhG/EnhPr/EnhPois from the 18-state cortex ChromHMM,
   though it is P0 cortex, not adult CA1); and does the Pcdh cluster form a
   multi-way hub, which is a good test case since its NPMI coverage is 97% and
   the alpha cluster's HS5-1 enhancer is a known multi-loop organiser.

   **A pairwise proxy exists but does not answer the question.** Triangle
   density in a thresholded NPMI graph ("do a bin's partners also contact each
   other?") is computable from what we have, but it cannot distinguish one cell
   holding three loci together from three cells each holding a pair - which is
   the entire distinction being asked about. If it is ever reported it must be
   labelled a pairwise clustering statistic, never a multi-contact result.

9. **What is `TAD_boundaries_mm39.bed`?** It is used as the "external" TAD
   partition in the Pcdh test, but nothing records its cell type, the caller
   that produced it, its original resolution or its source assembly. The file
   itself looks like a genuine liftover (per-feature width jitter, boundaries
   spaced on a 50 kb grid at a median 950 kb, padded to 150 kb windows), so it
   is not corrupt. But CTCF peaks show **no enrichment at it at any window
   size** (1.09x genome average at +/-5 kb through +/-100 kb), while the same
   CTCF track is 5.27x enriched at TSSs and 1.72x at our own GAM insulation
   minima. Either it is from a cell type far from neurons, or its boundary
   positions are only accurate to the width of the window. Until the provenance
   is known it should be treated as the weaker of the two partitions.

10. **What exactly is floxed in the Cbx3 conditional allele, and where?** The
    Cbx3 locus itself shows a sharp, local perturbation and we cannot currently
    tell whether that is biology or the knockout construct.

    Cbx3 is **chr6:51,447,340-51,460,684** (mm39, 13.3 kb), sharing a
    bidirectional promoter region with Hnrnpa2b1. The 20 kb bin that contains
    that promoter, `chr6:51440000-51460000`, carries **73 lost and 5 gained**
    differential contacts against genome-wide medians of 14 and 13 - the
    **97.8th percentile** for lost contacts, so it is a called hotspot bin at
    the >= 50 cutoff. The bin is TSS-dense (20 TSSs), and even among bins with
    >= 20 TSSs the median is 13.5, so 73 is high for its own class.

    Per-bin ChIP log2FC over that bin: **H3K36me3 -0.438** against a genome
    median of -0.052, with H4K20me3 -0.379 (*less* lost than the -0.497 median),
    H3K9me3 -0.167 (median -0.290) and H3K4me3 +0.036.

    **This is very probably the construct, not chromatin.** Emx1-Cre excises the
    floxed segment in the cells assayed, so sequence is physically missing:
    reads drop out over the deleted interval, ChIP log2FC there goes negative
    whatever the chromatin is doing, and GAM's detection of the altered locus is
    affected too. A collapse of H3K36me3 over a gene body is exactly what a
    deleted or silenced gene looks like, which makes this a **positive control
    that the knockout worked** rather than a finding.

    Everything around it is unremarkable, which supports that reading: the
    +/-1 Mb neighbourhood sits at the **54th percentile** of net contact change,
    no called hotspot domain overlaps it, no large block contains it (it falls
    between the chr6:27.6-46.0 lost block and the chr6:54.5-63.3 gained block),
    and the insulation change over +/-200 kb is at the **31st percentile** of
    |difference| genome-wide. The perturbation is one bin wide.

    To resolve it we need the floxed coordinates (which exons, which loxP
    positions) and ideally the targeting construct, so the affected interval can
    be masked in both assays. Until then the Cbx3 bin should be excluded from
    any per-bin analysis rather than read as a result. Related: question 5.
