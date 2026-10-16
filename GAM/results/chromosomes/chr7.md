# chr7

Generated 2026-10-17 by `GAM_7_Figures.R` (section B). Numbers here are regenerable; do not hand-edit.

## Mark-contact coupling, in context

- rho(H4K20me3 log2FC, net contact change) = **0.017**, rank 10 of 19.
- The spread across all 19 autosomes is -0.149 to 0.128 around a genome-wide 0.0119, so this rank is weak evidence on its own.
- Median H4K20me3 log2FC here -0.462, against -0.497 genome-wide.
- 20.2% of bins are heterochromatic (>=50% coverage by the merged H3K9me3/H4K20me3 set), against 7.2% genome-wide.

## Contact budget

Percent change in mean NPMI, HP1gKO vs WT, by chromatin class and separation.

| class | 0.02-1Mb | 1-5Mb | 5-10Mb | >10Mb |
|---|---|---|---|---|
| het-het | +8.9% | +10.5% | +11.8% | +8.7% |
| het-eu | +3.5% | +6.9% | +7.6% | +13.8% |
| eu-eu | +0.8% | +3.7% | +9.2% | +13.5% |

het-het minus eu-eu, per band (genome-wide this is positive: heterochromatin gains while euchromatin loses):

| 0.02-1Mb | 1-5Mb | 5-10Mb | >10Mb |
|---|---|---|---|
| +8.2 | +6.7 | +2.6 | -4.7 |

Bands where this chromosome runs against the genome-wide direction: >10Mb.

## Large one-directional blocks (GAM_6)

8 block(s) on this chromosome, 51.9 Mb in total.

| start (Mb) | end (Mb) | width (Mb) | direction | mean net |
|---|---|---|---|---|
| 9.58 | 19.44 | 9.86 | gained | -1.29 |
| 34.32 | 40.84 | 6.52 | lost | +1.39 |
| 47.30 | 54.78 | 7.48 | gained | -1.24 |
| 57.64 | 61.66 | 4.02 | lost | +1.36 |
| 90.38 | 95.34 | 4.96 | lost | +1.40 |
| 100.80 | 107.86 | 7.06 | lost | +1.52 |
| 109.46 | 115.94 | 6.48 | gained | -1.70 |
| 128.70 | 134.20 | 5.50 | lost | +1.34 |

## Caveats carried from the pipeline

- GAM is 2 biological replicates per genotype, already merged, so there is no within-genotype variance; every claim here is an effect size against a matched or bootstrap null, never a per-contact p-value.
- Per-bin ChIP log2FC at 20 kb is noise-dominated (per-bin sd 0.14-0.16 for the broad marks). Chromosome-level medians are sound; do not read local structure off a single bin.
- The A/B compartment axis is weakly determined on most chromosomes and is not used here.

