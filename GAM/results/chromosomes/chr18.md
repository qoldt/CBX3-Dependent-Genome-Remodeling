# chr18

Generated 2026-10-17 by `GAM_7_Figures.R` (section B). Numbers here are regenerable; do not hand-edit.

## Mark-contact coupling, in context

- rho(H4K20me3 log2FC, net contact change) = **-0.015**, rank 14 of 19.
- The spread across all 19 autosomes is -0.149 to 0.128 around a genome-wide 0.0119, so this rank is weak evidence on its own.
- Median H4K20me3 log2FC here -0.507, against -0.497 genome-wide.
- 4.8% of bins are heterochromatic (>=50% coverage by the merged H3K9me3/H4K20me3 set), against 7.2% genome-wide.

## Contact budget

Percent change in mean NPMI, HP1gKO vs WT, by chromatin class and separation.

| class | 0.02-1Mb | 1-5Mb | 5-10Mb | >10Mb |
|---|---|---|---|---|
| het-het | -2.4% | -7.6% | -11.1% | -23.4% |
| het-eu | -2.5% | -8.7% | -13.6% | -21.6% |
| eu-eu | -5.3% | -12.3% | -14.0% | -19.1% |

het-het minus eu-eu, per band (genome-wide this is positive: heterochromatin gains while euchromatin loses):

| 0.02-1Mb | 1-5Mb | 5-10Mb | >10Mb |
|---|---|---|---|
| +2.9 | +4.8 | +2.9 | -4.3 |

Bands where this chromosome runs against the genome-wide direction: >10Mb.

## Large one-directional blocks (GAM_6)

5 block(s) on this chromosome, 38.0 Mb in total.

| start (Mb) | end (Mb) | width (Mb) | direction | mean net |
|---|---|---|---|---|
| 7.90 | 15.68 | 7.78 | lost | +1.25 |
| 21.84 | 32.38 | 10.54 | lost | +1.97 |
| 33.88 | 41.18 | 7.30 | gained | -1.41 |
| 54.86 | 59.52 | 4.66 | gained | -1.00 |
| 82.96 | 90.72 | 7.76 | gained | -2.83 |

## Caveats carried from the pipeline

- GAM is 2 biological replicates per genotype, already merged, so there is no within-genotype variance; every claim here is an effect size against a matched or bootstrap null, never a per-contact p-value.
- Per-bin ChIP log2FC at 20 kb is noise-dominated (per-bin sd 0.14-0.16 for the broad marks). Chromosome-level medians are sound; do not read local structure off a single bin.
- The A/B compartment axis is weakly determined on most chromosomes and is not used here.

