# MINUTE (ChIP) parameters — why each value is what it is

`config.R` holds the values and `utils.R` the code; this file holds the
reasoning that used to sit in `config.R` as comment paragraphs. It was moved
here on 2026-11-23 without changing a value. If you change a parameter, update
its section here in the same commit.

## Input locations

Large raw inputs (bigWigs + repeat/TAD annotation) are not committed. Download
them from the Google Drive links in the README into the repo-relative folders
`ChIP/data/bigwig` and `ChIP/data/annotation` (the defaults; scripts run from
`ChIP/`). To use a copy staged elsewhere (e.g. a Synology/external mount), set
`MINUTE_DATA=/path/to/store Rscript run_MINUTE.R`.

`results_dir` is overridable so a sensitivity analysis can be written to a
separate tree without disturbing the committed run:
`MINUTE_RESULTS=results_norep2 MINUTE_EXCLUDE=... Rscript run_MINUTE.R`.

## Sample exclusion — HP1gKO rep2 by default

**Why.** Replicates are nested, not independent: WT rep1-3 and rep4-6 are
technical replicates of biological replicates 1 and 2; likewise HP1gKO rep2-3
(rep1's library failed) and rep4-6. HP1gKO rep2 is the highest-signal library in
all five marks, and it is badly discordant with rep3 — its own technical twin
(H4K20me3 median coverage 1.92 vs 1.14, while the bio-2 trio sits at
0.99-1.32). Technical replicates should agree; this one does not, so it is
treated as a technical failure rather than biological variation.

**Effect of excluding it:** the global loss deepens (H4K20me3 median log2FC
−0.474 → −0.678) while the control marks stay near zero, and the apparent
"gain" set roughly halves. About 0.04-0.07 of the deepening is a general offset
from dropping the highest KO library, so it is not all correction.

**Override** with `MINUTE_EXCLUDE`: a comma-separated list of sample_id values or
of `<genotype>:<replicate>` pairs that expand across every mark.
`MINUTE_EXCLUDE=none` keeps all samples (reproduces the pre-exclusion run):
`MINUTE_EXCLUDE=none MINUTE_RESULTS=results_allreps Rscript run_MINUTE.R`.

**Caveat that exclusion does not fix:** the DESeq2 design `~condition` still
treats technical replicates as independent, so the effective n is 2 biological
replicates vs 2, not 6 vs 4. Per-peak p-values are anticonservative on that
count — though the p-value distributions (MINUTE_2 section 7) show they are also
conservatively miscalibrated by the global shift, which dominates.

## Regions quantified per mark

H3K9me2 is quantified over the shared H3K9me3/H4K20me3 regions, not its own
peaks (changed 2026-11-20). MACS3 models H3K9me2 poorly — no per-replicate call
set looks right over the signal — and `merged_H3K9me2_masterPeak.bed` (38,787
narrow regions, median 323 bp, 17 Mb) has no committed recipe and does not
derive from any call set on disk. Provenance of every region file:
`results/METHODS.md`, "Peak / region sets".

## Annotation

ChromHMM is the 18-state cortex (P0) segmentation, mm39. It is annotated per
peak in MINUTE_1 (column `chromHMM_state`), enriched in MINUTE_2 and used in
MINUTE_4. KAP1/TRIM28 peaks were prepared by liftOver from ChIP-Atlas mm10 (see
the README methods); MINUTE_6 skips its analysis if they are absent.

## Significance

`sig_thresholds` is the single definition of "significant" — `|log2FoldChange| >
lfc` and raw `pvalue < p`, per mark — applied by `is_significant()` in every
stage that needs a significant set, so those sets are guaranteed identical.
Significance uses the raw p-value, not BH-adjusted padj.

## Gene families

The HUSH/CBX3-target families are quantified directly over their exons, rather
than relying on individual genes passing per-peak significance. The exon rule
picks the exon that carries the H3K9me3/HUSH silencing, per family, so the
biologically relevant exon is quantified per gene rather than diluted across all
exons:

- KRAB-ZFP → `last` (3'-most exon: the ZNF array)
- Protocadherin → `first` (5'-most exon: the large variable exon)
- Vomeronasal / Olfactory → `largest` (the single coding exon)
