# Methods notes

## KAP1 / TRIM28 ChIP-seq data source and preparation

Mouse KAP1 (TRIM28) ChIP-seq peaks were obtained from **ChIP-Atlas**
(https://chip-atlas.org), which provides uniformly reprocessed public ChIP-seq
peaks (MACS2, q < 1×10⁻⁵). No primary adult-mouse-cortex KAP1 ChIP-seq with
public processed peaks was available, so two tracks were used: (i) a
**neural-lineage** track from Neuro-2a neuroblastoma cells (GEO **GSE110032**,
KAP1 ChIP-seq, replicates GSM2986446 / GSM2986447; 15,355 peaks), used as the
primary track for its neural relevance; and (ii) the **all-cell-type aggregate**
of 47 mouse Trim28 experiments (276,439 peaks; ChIP-Atlas `Oth.ALL.05.Trim28`),
used as a well-powered robustness track, since KAP1 occupancy at its
ERV/H3K9me3-heterochromatin targets is largely conserved across cell types.

Peaks were downloaded in **mm10** coordinates
(`https://dbarchive.biosciencedbc.jp/kyushu-u/mm10/assembled/Oth.ALL.05.Trim28.*.bed`)
and lifted to **mm39** with **CrossMap** using the UCSC `mm10ToMm39.over.chain`
chain (`https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz`);
99.7 % (neural) and 99.0 % (all-cell) of peaks mapped. Lifted peaks were
restricted to the primary chromosomes and converted to Ensembl-style seqnames
(1–19, X, Y) to match the input-scaled mm39 bigWigs.

Overlap of the differential-loss compartments (and the HUSH/CBX3 gene-family
exons) with KAP1 peaks was tested against a **size-matched permutation null**
(N = 2000): for each set, regions were resampled from the appropriate universe
(the measured merged-peak set for the loss compartments; all annotated mm39 gene
exons for the gene families) matched to the query's genomic-length distribution
(20 length bins), and the overlap fraction was recomputed to build the null.
Enrichment is reported as observed/expected fold and empirical p (MINUTE_6).

## Note on HUSH targets (Danac et al., Mol Cell 2024)

The HUSH-target data in Danac et al. (Mol Cell 2024; PMID 39013473; GEO
GSE268799) is **human (hg38)** CUT&RUN in K562, and its dominant HUSH targets are
primate-specific young LINE-1 (L1HS/L1PA) with no orthologous mouse locus — so
coordinate liftOver to mm39 is not valid. The mouse-native equivalent is used
instead: enrichment for **young/active mouse L1 families (L1MdA/T/Gf/F)** against
the mm39 RepeatMasker annotation (MINUTE_5).

## Peak / region sets: what is quantified, and where each file came from

Recorded 2026-11-20 after checking every file in `data/peaks/` against the
per-replicate MACS3 calls on the lab Synology (`SynologyDrive/MINUTE/peaks/`).

**What the analysis quantifies** (`config.R` `regions`):

| Mark | File | Provenance |
|---|---|---|
| H3K9me2, H3K9me3, H4K20me3 | `2000bp_merge_H3K9me3_H4K20me3.bed` (126,197) | Reproduced exactly as `reduce(c(H3K9me3, H4K20me3), min.gapwidth = 2000)` of `merged_H3K9me3_masterPeak.bed` and `merged_H4K20me3_masterPeak.bed` (both on the Synology, Feb 2025). The merge step is not coded; the two inputs have **no committed recipe** and do not derive from any call set now on disk. |
| H3K4me3 | `master_peaks_consensus/H3K4me3_consensus_masterPeak.bed` | Byte-identical copy also at `data/peaks/H3K4me3_consensus_masterPeak.bed`. `MergePeaks_to_masterPeak.R` does not list H3K4me3, so no committed code makes it. |
| H3K36me3 | `merged_H3K36me3_masterPeak.bed` | **No committed recipe** (Feb 2025, same origin as the merged H3K9me3/H4K20me3 sets). |

**H3K9me2 moved onto the shared set on 2026-11-20.** MACS3 models H3K9me2
poorly - no per-replicate call set looks right over the signal - and the file
it used before, `merged_H3K9me2_masterPeak.bed` (38,787 regions, median 323 bp,
17.2 Mb), has no recipe and overlaps at most 70% of any call set on disk. The
old file is kept in `data/peaks/` and in the IGV session as
`retired_H3K9me2_merged` for comparison.

**What `Peak Calling/MergePeaks_to_masterPeak.R` actually produced.** The
`master_peaks_consensus/*` sets for H3K9me2, H3K9me3, H4K20me3 and H3K36me3
(>= 2 of the 11 libraries, HP1gKO rep2 included) are reproduced exactly from
the calls in `broad_1000min_5000maxgap/` - **not** `maxgap2000/`, which is the
path hard-coded in the script, and not the `--max-gap 2000`, no-minimum-length
settings in the committed `Snakefile`. So the MACS3 parameters behind the
consensus sets (min length 1000, max gap 5000) are not in the repo either.
These consensus sets are not used by the analysis.

**Open:** where the Feb 2025 `merged_*_masterPeak.bed` files came from (which
calls, which merge). Until that is known, every peak-level H3K9me3, H4K20me3 and
H3K36me3 number rests on a region set that cannot be regenerated from this repo.
The 20 kb bin results (GAM_1) do not depend on any peak file.

**Quantisation caveat (applies to every mark, sharpest for H3K9me2).** MINUTE_1
rounds the mean coverage per region to an integer for DESeq2. H3K9me2's mean
coverage over these regions is ~1.4, so per-region log2FC takes few values
(517 distinct across 77,122 regions; -0.222 alone is 10,592). Per-region
H3K9me2 "gains" are therefore mostly 1<->2 rounding steps; read the
distribution, not individual regions.
