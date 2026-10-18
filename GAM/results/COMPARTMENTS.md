# A/B compartments in this dataset: how they are called, and why it is hard

Hand-written, committed. Written 2026-10-18 because the reasoning here is
spread across four scripts and will not be reconstructable from the code.

**Short version.** The textbook Hi-C recipe does not work out of the box on this
GAM data, because the leading eigenvector of the correlation matrix is a
**per-bin detectability component, not compartments**. The pipeline currently
works around that by hunting down the eigenvector list for one that looks like
an A/B split. That salvage is the source of nearly every compartment problem in
this project. Balancing the matrix before the PCA fixes the root cause on
chromosomes that have a compartment, and shows that several chromosomes simply
do not have one.

---

## 1. What is being computed

A and B are the two spatial compartments of the nucleus - active/gene-rich and
inactive. Loci in the same compartment contact each other more, which appears in
a contact matrix as a checkerboard. The standard recipe:

1. remove the distance effect (observed minus expected, by separation)
2. Pearson correlation matrix of that
3. **PCA - PC1 is the compartment axis**
4. orient PC1 with an active mark so positive = A

Step 4 here uses `active_chromatin_ref()` in `config.R`: the mean of per-bin WT
H3K4me3 and H3K36me3 z-scores. Two marks rather than one because H3K4me3 is
punctate at promoters and thin across a bin, while H3K36me3 is broad over
transcribed gene bodies and better matched to this resolution. Both are
UNCHANGED in the KO (per-bin median log2FC +0.000 and -0.052, against -0.497 for
H4K20me3), which is what makes the axis the same ruler in both genotypes.

Code path: `GAM_2_NPMI_Compartments.R` `reduce_npmi()` -> aggregate to
`gam_compartment_binsize` -> `pc1_corr()` -> `select_compartment_pc()`
(`config.R`) -> `orient_sign()`.

## 2. Where it breaks: step 3, not step 4

**PC1 is the compartment axis on 2 of 24 determined combinations.** It is PC2 on
14, PC3 on 7, PC4 on 1.

The reason is in `config.R`'s own note on `select_compartment_pc()`: the leading
eigenvector is "systematically degenerate (frac_A 0.04-0.13 - almost every bin
on one side)". Measured directly on the unbalanced O-E matrix:

| chromosome | frac_A of PC1, unbalanced |
|---|---|
| chr19 WT | **0.033** |
| chr3 WT | **0.983** |

A vector that puts 97-98% of a chromosome on one side is not an A/B split
whatever it correlates with. So the pipeline selects on **balance first, then
correlation** - it looks for a candidate with `frac_A` in [0.2, 0.8] and takes
the best |rho| among those.

That is a reasonable salvage and it is also the root of the downstream trouble:
we end up using the 2nd-4th component, carrying a few percent of the variance,
where Hi-C would use a dominant 1st. Even the selected axis reaches only a
median |rho| of **0.322** against active chromatin (range 0.034-0.481). A
well-determined Hi-C compartment axis is 0.6-0.8.

**The orienting marks are not the limitation.** The two marks agree with each
other at rho 0.764 at 100 kb, so a perfect axis could reach |rho| ~0.93 against
this reference. We observe 0.322. The ceiling is the axis, not the ruler.

## 3. Why PC1 is degenerate: it is detectability

A window detected in more or fewer nuclear profiles shifts its **entire row and
column** of the matrix. That is a rank-ish per-bin component, and a per-bin
component in a correlation matrix produces exactly the signature observed - a
dominant leading eigenvector with nearly every bin on one side.

The same effect was quantified independently, for a different reason, on the
difference maps: fitting `Dz[i,j] ~ a[i] + a[j]` accounts for **18-28%** of the
variance in every region tested, worst in a plain control window. See
`BIN_BALANCING.md`.

### The experiment

Additive per-bin balancing of the O-E matrix *before* the correlation and PCA,
compared against the current unbalanced path. Same inputs, same aggregation,
WT only:

| | PC1 frac_A | PC1 \|rho\| | PC1 var% | best \|rho\| over PC1-4 |
|---|---|---|---|---|
| chr11 unbalanced | degenerate | 0.458 | 17.6 | 0.518 (PC2) |
| **chr11 balanced** | **0.561** | **0.507** | 11.7 | **0.507 (PC1)** |
| chr19 unbalanced | 0.033 | 0.281 | 16.9 | 0.281 (PC1) |
| chr19 balanced | 0.535 | 0.206 | 10.8 | 0.320 (PC4) |
| chr3 unbalanced | 0.983 | 0.086 | 11.0 | 0.086 (PC1) |
| chr3 balanced | 0.481 | 0.002 | 5.5 | 0.082 (PC2) |

Three things follow, and the third matters as much as the first:

1. **Balancing fixes the degeneracy everywhere.** frac_A moves from 0.033 and
   0.983 to 0.48-0.56. PC1 becomes an actual split rather than a lopsided
   coverage component. PC1's variance share drops (17.6% -> 11.7% on chr11),
   which is the removed per-bin component.
2. **On a chromosome that has a compartment, it works.** chr11 PC1 becomes the
   compartment at |rho| 0.507, slightly better than the 0.518 the eigenvector
   hunt currently recovers from PC2 - and with no hunt needed.
3. **It does not manufacture a compartment where there is none.** chr3 and
   chr19 stay at |rho| <= 0.32 after balancing. Those chromosomes genuinely
   lack a strong compartment axis in this data; the eigenvector hunt was
   dressing that up as a call.

## 4. The downstream damage, and what was done about it

Problems traced to the salvage, all fixed or documented 2026-10-18:

- **The cKO's axis is selected to match WT's.** `GAM_2:421-428` passes WT's own
  PC1 as the cKO's `orient_ref`, and `select_compartment_pc()` ranks candidates
  by correlation with that reference. So the cKO eigenvector is chosen for
  resembling WT's. Consequences: the cKO can never be called where WT fails
  (0 of 19 chromosomes), and the report's "the two genotypes select the same
  eigenvector on only 8 of 19" UNDERSTATES the disagreement - they are not
  independent selections, the cKO is steered toward WT and still disagrees on
  11. **Not yet fixed.**
- **No margin requirement.** argmax|rho| with no margin made the answer depend
  on how many eigenvectors were offered. `gam_compartment_pc_margin` = 0.05
  now requires the winner to beat the runner-up. Result: 24 of 38 combos
  oriented, 10 chromosomes with both genotypes.
- **Two chromosomes shipped inverted.** chr15 and chr19 cKO had A and B swapped
  (rho -0.061 and -0.134 against the independent H3K4me3 check). It was a
  warning nobody acted on, and `GAM_2b:309` stored `abs(rho_h3k4me3)`, so the
  table that exists to catch it was hiding it. Both fixed, plus an assertion at
  the point the table is written.
- **`pc1_delta` is not a usable genotype comparison** and its three figures were
  removed from the report.

## 5. What to do next

**Balance before the PCA.** That is the root-cause fix: it removes the
eigenvector hunt on good chromosomes and stops dressing up weak ones. It should
replace `select_compartment_pc()`'s balance-first heuristic, not sit alongside
it. `normGAM` (PMC6936146) provides ICE and Knight-Ruiz for GAM data if a
multiplicative balancing on the raw matrix is preferred to the additive one on
O-E used here.

**Then fix the cKO reference** to use the same independent ChIP ruler as WT, so
the two genotypes are selected independently and a difference between them means
something.

**Expect fewer chromosomes, not more.** Balancing is likely to confirm that
several chromosomes have no compartment axis in this data. That is the honest
outcome and it is better than the current position, which is a call on every
chromosome at |rho| 0.03-0.48.

Until both are done, use the **shared axis** (`GAM_2b`), where both genotypes
are scored on one vector, and treat per-genotype A/B as a coarse annotation on
the chromosomes listed in `compartments_WT.bed` rather than as a measurement.
