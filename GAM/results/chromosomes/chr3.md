# chr3

Generated 2026-10-17 by `GAM_7_Figures.R` (section B). Numbers here are regenerable; do not hand-edit.

## Mark-contact coupling, in context

- rho(H4K20me3 log2FC, net contact change) = **0.059**, rank 4 of 19.
- The spread across all 19 autosomes is -0.149 to 0.128 around a genome-wide 0.0119, so this rank is weak evidence on its own.
- Median H4K20me3 log2FC here -0.543, against -0.497 genome-wide.
- 6.7% of bins are heterochromatic (>=50% coverage by the merged H3K9me3/H4K20me3 set), against 7.2% genome-wide.

## Contact budget

Percent change in mean NPMI, HP1gKO vs WT, by chromatin class and separation.

| class | 0.02-1Mb | 1-5Mb | 5-10Mb | >10Mb |
|---|---|---|---|---|
| het-het | -0.9% | -0.9% | +9.2% | -0.6% |
| het-eu | -0.2% | +0.9% | +4.5% | -0.5% |
| eu-eu | +0.8% | +1.5% | +1.2% | -5.6% |

het-het minus eu-eu, per band (genome-wide this is positive: heterochromatin gains while euchromatin loses):

| 0.02-1Mb | 1-5Mb | 5-10Mb | >10Mb |
|---|---|---|---|
| -1.7 | -2.4 | +8.0 | +4.9 |

Bands where this chromosome runs against the genome-wide direction: 0.02-1Mb, 1-5Mb.

## Large one-directional blocks (GAM_6)

11 block(s) on this chromosome, 109.5 Mb in total.

| start (Mb) | end (Mb) | width (Mb) | direction | mean net |
|---|---|---|---|---|
| 2.80 | 7.42 | 4.62 | lost | +1.42 |
| 13.74 | 28.52 | 14.78 | lost | +1.85 |
| 28.98 | 36.02 | 7.04 | gained | -1.65 |
| 43.48 | 52.18 | 8.70 | lost | +1.18 |
| 52.90 | 64.10 | 11.20 | gained | -2.07 |
| 66.96 | 75.66 | 8.70 | lost | +1.41 |
| 76.40 | 88.20 | 11.80 | gained | -2.04 |
| 96.22 | 100.24 | 4.02 | lost | +1.64 |
| 108.24 | 120.30 | 12.06 | lost | +2.48 |
| 120.86 | 131.66 | 10.80 | lost | +1.94 |
| 143.96 | 159.75 | 15.79 | gained | -2.10 |

## Caveats carried from the pipeline

- GAM is 2 biological replicates per genotype, already merged, so there is no within-genotype variance; every claim here is an effect size against a matched or bootstrap null, never a per-contact p-value.
- Per-bin ChIP log2FC at 20 kb is noise-dominated (per-bin sd 0.14-0.16 for the broad marks). Chromosome-level medians are sound; do not read local structure off a single bin.
- The A/B compartment axis is weakly determined on most chromosomes and is not used here.

