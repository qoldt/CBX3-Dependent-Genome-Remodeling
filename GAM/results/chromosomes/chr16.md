# chr16

Generated 2026-10-17 by `GAM_7_Figures.R` (section B). Numbers here are regenerable; do not hand-edit.

## Mark-contact coupling, in context

- rho(H4K20me3 log2FC, net contact change) = **0.040**, rank 7 of 19.
- The spread across all 19 autosomes is -0.149 to 0.128 around a genome-wide 0.0119, so this rank is weak evidence on its own.
- Median H4K20me3 log2FC here -0.515, against -0.497 genome-wide.
- 4.8% of bins are heterochromatic (>=50% coverage by the merged H3K9me3/H4K20me3 set), against 7.2% genome-wide.

## Contact budget

Percent change in mean NPMI, HP1gKO vs WT, by chromatin class and separation.

| class | 0.02-1Mb | 1-5Mb | 5-10Mb | >10Mb |
|---|---|---|---|---|
| het-het | +0.3% | +1.2% | -8.4% | -13.0% |
| het-eu | -2.9% | -8.6% | -12.9% | -23.6% |
| eu-eu | -3.7% | -13.2% | -24.8% | -28.7% |

het-het minus eu-eu, per band (genome-wide this is positive: heterochromatin gains while euchromatin loses):

| 0.02-1Mb | 1-5Mb | 5-10Mb | >10Mb |
|---|---|---|---|
| +4.0 | +14.3 | +16.4 | +15.7 |

Bands where this chromosome runs against the genome-wide direction: none.

## Large one-directional blocks (GAM_6)

5 block(s) on this chromosome, 48.1 Mb in total.

| start (Mb) | end (Mb) | width (Mb) | direction | mean net |
|---|---|---|---|---|
| 18.30 | 27.76 | 9.46 | gained | -1.97 |
| 31.70 | 40.44 | 8.74 | lost | +2.20 |
| 43.52 | 52.56 | 9.04 | gained | -1.62 |
| 73.14 | 87.08 | 13.94 | lost | +2.61 |
| 91.10 | 98.01 | 6.91 | gained | -2.03 |

## Caveats carried from the pipeline

- GAM is 2 biological replicates per genotype, already merged, so there is no within-genotype variance; every claim here is an effect size against a matched or bootstrap null, never a per-contact p-value.
- Per-bin ChIP log2FC at 20 kb is noise-dominated (per-bin sd 0.14-0.16 for the broad marks). Chromosome-level medians are sound; do not read local structure off a single bin.
- The A/B compartment axis is weakly determined on most chromosomes and is not used here.

