# chr11

Generated 2026-10-17 by `GAM_7_Figures.R` (section B). Numbers here are regenerable; do not hand-edit.

## Mark-contact coupling, in context

- rho(H4K20me3 log2FC, net contact change) = **0.128**, rank 1 of 19.
- The spread across all 19 autosomes is -0.149 to 0.128 around a genome-wide 0.0119, so this rank is weak evidence on its own.
- Median H4K20me3 log2FC here -0.415, against -0.497 genome-wide.
- 3.1% of bins are heterochromatic (>=50% coverage by the merged H3K9me3/H4K20me3 set), against 7.2% genome-wide.

## Contact budget

Percent change in mean NPMI, HP1gKO vs WT, by chromatin class and separation.

| class | 0.02-1Mb | 1-5Mb | 5-10Mb | >10Mb |
|---|---|---|---|---|
| het-het | +7.4% | -14.7% | -33.8% | -24.7% |
| het-eu | -2.9% | -10.3% | -17.4% | -25.3% |
| eu-eu | -4.8% | -12.4% | -22.0% | -27.8% |

het-het minus eu-eu, per band (genome-wide this is positive: heterochromatin gains while euchromatin loses):

| 0.02-1Mb | 1-5Mb | 5-10Mb | >10Mb |
|---|---|---|---|
| +12.3 | -2.3 | -11.8 | +3.1 |

Bands where this chromosome runs against the genome-wide direction: 1-5Mb, 5-10Mb.

## Large one-directional blocks (GAM_6)

7 block(s) on this chromosome, 65.0 Mb in total.

| start (Mb) | end (Mb) | width (Mb) | direction | mean net |
|---|---|---|---|---|
| 2.72 | 23.30 | 20.58 | lost | +2.39 |
| 42.60 | 50.12 | 7.52 | lost | +1.74 |
| 57.98 | 66.04 | 8.06 | gained | -1.26 |
| 66.18 | 71.46 | 5.28 | gained | -1.43 |
| 83.38 | 97.00 | 13.62 | gained | -1.84 |
| 97.18 | 103.16 | 5.98 | gained | -1.03 |
| 103.48 | 107.48 | 4.00 | gained | -1.20 |

## Caveats carried from the pipeline

- GAM is 2 biological replicates per genotype, already merged, so there is no within-genotype variance; every claim here is an effect size against a matched or bootstrap null, never a per-contact p-value.
- Per-bin ChIP log2FC at 20 kb is noise-dominated (per-bin sd 0.14-0.16 for the broad marks). Chromosome-level medians are sound; do not read local structure off a single bin.
- The A/B compartment axis is weakly determined on most chromosomes and is not used here.

