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
  (rho -0.061 and -0.134 against the independent check). It was a warning
  nobody acted on, and `GAM_2b` stored the ABSOLUTE value of that correlation,
  so the table that exists to catch it was hiding it. Both fixed, plus an
  assertion at the point the table is written.
- **`pc1_delta` is not a usable genotype comparison** and its three figures were
  removed from the report.

## 4b. Both root-cause fixes are now implemented (2026-10-18)

**Balancing before the PCA** (`gam2_balance_before_pca`) and **both genotypes
selected on the same external ruler** (the cKO no longer gets WT's own
eigenvector). Measured against the previous state:

| | before | after |
|---|---|---|
| compartment axis is PC1 | 2 of 24 | **16 of 27** |
| oriented | 24 of 38 | 27 of 38 |
| WT median \|rho\| vs the active ruler | 0.340 | 0.358 |
| cKO called on a chromosome where WT fails | **0 of 19** | **3** |

Read these carefully, because two of them say less than they appear to.

- **The eigenvector hunt is largely gone.** PC1 carries the compartment on 16 of
  27 rather than 2 of 24. This is the fix working: the detectability component
  is out of the decomposition, so the leading eigenvector is the compartment.
  `which_pc` is now a DIAGNOSTIC - a chromosome still selecting PC3 is saying
  something about itself.
- **The cKO is genuinely independent.** `WT fails / cKO succeeds` was
  structurally impossible before (0 of 19) because an undetermined WT axis is
  all NA and the cKO had nothing to correlate against. It is now 3.
- **|rho| barely moved**, 0.340 -> 0.358. Balancing did NOT make the axis agree
  better with active chromatin; it made the right eigenvector carry it. The
  ceiling is the data. Compartment structure in this GAM dataset is weak, and
  no amount of better decomposition will change that - which is what the margin
  is for.

**`rho_check`, not `rho_h3k4me3`.** The orientation check is now TSS density
from the annotation, which involves no ChIP at all. It had to change: once both
genotypes are oriented on the active-chromatin ruler, a check against that same
ruler is checking the reference against itself, which is the flaw the WT path
already had. The column was renamed because keeping the old name would have
been a false claim about the source - it briefly fooled a before/after
comparison during this very session.

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

## 4c. The shared axis now balances too — and it is a trade (2026-10-18)

Section 4b fixed the **per-genotype** axes (`GAM_2`). But the BEDs,
`compartment_switch`, `GAM_6`'s decompaction verdict and the figures all run off
the **shared** axis (`GAM_2b`), and that one still formed its correlation matrix
from the raw O-E map. It showed exactly the pathology section 3 describes: the
compartment landed on PC1 for **0 of 18** determined chromosomes (PC2 x8,
PC3 x8, PC4 x2).

`gam2b_balance_before_pca` (params v2) applies the same estimator per genotype,
before the two correlation matrices are averaged.

**The mechanism is confirmed.** 15 of 16 chromosomes moved up *exactly one*
eigenvector (3->2, 2->1, 4->3). One eigenvector was being spent on per-bin
detectability; removing it shifts the whole spectrum by one.

| | unbalanced (v1) | balanced (v2) |
|---|---|---|
| compartment on PC1 | 0 of 18 | **8 of 16** |
| determined | **18 of 19** | 16 of 19 |
| median \|rho\| vs active chromatin | **0.423** | 0.391 |
| median \|frac_A - 0.5\| | 0.061 | **0.043** |
| worst \|frac_A - 0.5\| | 0.198 | **0.117** |
| bin coverage for the BEDs | **90.5%** | 82.7% |

**It does not improve agreement with active chromatin.** On the 15 chromosomes
determined in both versions: 8 up, 7 down, mean 0.383 -> 0.391, median
0.422 -> 0.387. A wash. Same conclusion as 4b — the ceiling is the data. What
balancing buys is which eigenvector carries the compartment, and a more even
A/B split.

**It costs chr8, chr14 and chr15**, which now fail the selector entirely
(chr4 is recovered, so 18 -> 16). Not an eigenvector-budget problem: tested
k = 6, all three still fail with "no balanced eigenvector among the top 6".

**Why the loss is probably honest.** In Hi-C, balancing removes a near-pure
technical bias. In GAM, detection frequency is *partly biological* — open
chromatin is more detectable — so the per-bin term is part artefact and part
compartment. A detectability vector therefore correlates with active chromatin
on its own, and can pass an orientation check while not being the compartment.
chr14 was "determined" before on **PC4** at |rho| 0.527, the second-best score
in the genome, which is what that failure mode looks like. The balanced run
refuses it instead. That reading is an inference, not a measurement, and it is
the part of this change most worth disagreeing with.

**To reverse:** `gam2b_balance_before_pca <- FALSE` in `GAM_2b_Shared_Axis.R`
restores v1 exactly. `GAM2B_K=<n>` re-runs the k sensitivity test.

**One thing to watch.** Shared-axis coverage fell to 82.7%, and `GAM_6` asserts
`ov > 0.80`. That is 2.7 points of headroom. If another chromosome drops out,
`GAM_6` will stop rather than quietly produce a verdict on two-thirds of the
genome — which is the intended behaviour, but it will look like a crash.

## 6. The strength null does not regenerate from committed code (2026-10-19)

Sections 3-5 are about how badly the A/B *axis* behaves. This section is about
the one compartment result that does not depend on an axis at all - and it has a
different problem.

The strength test projects both genotypes onto rulers external to both contact
matrices and finds nothing:

| ruler | mean delta | 95% CI | p |
|---|---|---|---|
| active chromatin (WT H3K4me3 + H3K36me3 z) | -0.00018 | -0.00213 to +0.00186 | 0.62 |
| TSS density | +0.00019 | -0.00057 to +0.00103 | 0.74 |

Both over all 19 autosomes, bootstrap over chromosomes, paired Wilcoxon.
Negative on 10 of 19 for the active-chromatin ruler.

**This is the strongest compartment statement in the analysis and the least
reproducible part of it.** It was computed on 2026-10-17 and recorded in `f5096fa`
as prose in `CLAUDE.md` and `make_report.R`. It was never persisted:

- no file in `results/tables/` contains it;
- no committed script computes it (`grep -rn '0.00018' --include='*.R'` returns
  `make_report.R` only, where it is a literal in a `<div class='warn'>`);
- `verify_GAM.R` therefore cannot check it - its facts are computed *from*
  tables, so a statistic with no table generates no fact and no failure. A
  green 39/39 says nothing about this number.

That breaks the repo rule that everything in `results/` regenerates from
committed code plus the raw data, and it breaks it for exactly the number a
reviewer is most likely to press on.

**What closing it needs.** A stage - `GAM_2d`, or a section of `GAM_2b`, which
already loads both correlation matrices and `active_chromatin_ref()` - that for
each autosome scores both genotypes on each external ruler, computes the
strength contrast, bootstraps over chromosomes, and writes
`gam2d_strength_external_axes.tsv`. Then register it in `verify_GAM.R` as a
fact so the prose is checked like everything else. Until that exists, quote the
two rows above as provisional and say where they came from.
