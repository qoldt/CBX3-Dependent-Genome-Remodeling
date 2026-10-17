# GAM — Genome Architecture Mapping (mm39, 20 kb)

HP1γ cKO vs wild-type GAM, integrated with the MINUTE ChIP-seq results. Read the
root `../CLAUDE.md` first for the cross-assay conventions; this file holds only
what is specific to GAM. Data inventory and formats: `README.md`. Unresolved
provenance: `results/METHODS.md`.

## Pipeline

Run from the **`GAM/` directory** (no `setwd()`):

```
config.R                      parameters, path builders, loaders (single source of truth)
GAM_1_Bins_and_Signal.R       ChIP signal/stats + annotation on the 20 kb grid     (heavy)
  │  writes results/rds/gam_bins_annotated.rds
GAM_2_NPMI_Compartments.R     NPMI -> PC1/compartments/insulation/decay            (heavy)
  │  needs GAM_1 (PC1 orientation); writes results/rds/npmi_bin_features.rds
  ├── GAM_2b_Shared_Axis.R    ONE A/B axis, both genotypes scored on it     (heavy)
  │    needs GAM_1 + NPMI only; writes results/rds/shared_axis_scores.rds
  │    GAM_6 joins its score_delta, so run it before GAM_6
  ├── GAM_2c_Meltron.R        MELTRON melting / condensation                (heavy)
  │    needs NPMI only; writes meltron_is.rds + meltron_scores.rds
  ├── GAM_3_Coupling.R        does contact remodelling coincide with the mark loss?
  ├── GAM_4_Large_Scale.R     10-25 Mb one-directional blocks + block contact budget
  ├── GAM_5_Heterochromatin.R heterochromatin domains, intra vs inter, mark strata  (heavy)
  ├── GAM_6_Targets.R         HUSH families, repeats, protocadherin, the Pcdh TAD
  └── GAM_7_Figures.R         region plots + per-chromosome overviews               (heavy)

utilities (not pipeline stages, run on their own):
  QC_Hotspot_Domains.R        hotspot domain parameter sweep, for the collaborating lab
  prepare_LAD_annotation.R    fetch + liftOver the constitutive LAD sets to mm39
  prepare_CTCF_annotation.R   fetch + liftOver ENCODE CTCF peaks to mm39 (BULK
                              hippocampus + cortex, NOT sorted CA1 - see its header)
  export_igv.R                refresh the minimal IGV session
  make_report.R               assemble results/GAM_REPORT.html from whatever has run
  run_overnight.sh            run the batch under the thermal guard, checking outputs
```

**Stages were consolidated 2026-10-17 (15 scripts -> 7).** Each merged stage keeps
its sections in `local({...})` blocks labelled with the stage they came from, and
**every output file name is unchanged** (`gam6_*`, `gam8_*`, `gam9_*`, `gam10_*`,
`gam11_*`, `gam13_*`, `gam14_*`, `gam15_*`, `coupling_*`), so the prefixes no
longer match the file they live in - grep the prefix, not the filename. The merge
was verified by re-running every stage and diffing all 50 result tables
byte-for-byte against the pre-merge outputs.

`Rscript run_GAM.R` runs them in order; `./run_overnight.sh` runs the analysis
batch under the thermal guard and asserts each stage's outputs are non-empty. It also needs the ChIP pipeline to have
run (`../ChIP/results/rds/annotated_results_*.rds`; `require_chip_rds()` stops
with an actionable message otherwise). GAM_1, GAM_2, GAM_5 and GAM_7 are sustained all-core
jobs: run them under the `safe-local-alignment` skill. GAM_3-7 should `stop()`
with "run GAM_1 first" if an rds is missing; they must not recompute anything
GAM_1/GAM_2 provide, and join on `bin_id`.

## Single source of truth

`config.R` holds every parameter, the genotype bridge, path builders and loaders.
**Stages must not redefine them** (no local bin size, z cutoff, file paths or raw
`fread()` of a data file). To change sigma/z/NPMI, edit `config.R`; the path
builders then resolve to different files. Cross-assay helpers (`save_fig`,
`size_perm`, `hyper_row`, `to_ucsc`/`to_ncbi`, `mm39_bins`, `gam_bin_id`,
`parse_gam_bin`, palettes) live in `../shared/common.R`. If a stage needs a
parameter or helper that does not exist, add it to the right file rather than
inlining it. GAM-only dependency: `R.utils` (for gzipped `fread`).

## Gotchas that silently corrupt results

1. **Genotype is encoded ONLY by the directory name.** The NPMI filenames are
   identical in `WT_at20Kb/` and `cKO_at20Kb/`. Never resolve a matrix by
   basename or copy one out of its folder; always go through
   `npmi_file(genotype, chr)` / `load_npmi(genotype, chr)`.
2. **Direction semantics.** `selWT` = contact stronger in WT = **lost** in the KO;
   `selcKO` = **gained**. Columns are `n_contacts_lost` / `n_contacts_gained`.
   Reversing these inverts the headline; restate the mapping in a comment wherever
   it is used. In the `.bg2` files `set1` is WT and `set2` is cKO.
3. **Hotspots are dense, not sparse.** About 93 % of the 123,133 autosomal bins
   carry at least one differential contact. Do not treat presence as a hit list or
   count "hotspot bins"; the signal is the count magnitude (1-162). Absent bins
   are 0, not NA (`load_hotspots()` already fills them).
4. **Autosomes only for differential products; X/Y exist only in NPMI.**
   `Diff_contacts/` and `Hotspots/` cover chr1-19. NPMI covers chr1-19 + chrX +
   chrY (`gam_npmi_chroms`). Cross-assay integration stays on `gam_autosomes`;
   X/Y are for NPMI-only analyses.
5. **Replication: 2 biological replicates per genotype, already merged.** The
   effective n is 2 vs 2, matching ChIP. The delivered files carry no replicate
   dimension, so the only significance is the baked-in z > 2, and there is no
   within-genotype variance to test against. This forbids: a p-value or FDR on an
   individual contact or bin, any replicate-level test, and treating 25 million
   pairs (or 123k bins) as independent observations. Every statistical claim must
   come from matched-background permutation on the bin side (`size_perm()`,
   composition-matched nulls), reported as effect sizes. State "2 vs 2" in any
   write-up.
6. **One NPMI matrix at a time.** chr1 is about 9.7k x 9.7k, about 760 MB as
   doubles, and both genotypes double that. Load a chromosome, extract the
   per-bin features, `rm()` and `gc()`, then move on. Never hold all chromosomes
   or both genotypes of a large chromosome together. Diff contacts are about 25 M
   pairs per direction: load per chromosome with `select=` for needed columns.
7. **UCSC vs NCBI.** GAM files are **UCSC** (`chr1`); the ChIP peak `.rds` is
   **NCBI** (`1`). Flip with `to_ucsc()` / `to_ncbi()`. A botched flip yields zero
   overlaps, which looks like a clean null: assert overlap rates with
   `stopifnot()`.
8. **Coordinates.** GAM bin ids (`chr19:0-20000`) and file coordinates are
   0-based half-open; GRanges is 1-based. Use `parse_gam_bin()` / `gam_bin_id()`,
   never hand-rolled arithmetic.
9. **Genotype labels.** Canonical is `WT` / `HP1gKO`. Only `config.R` ever sees
   `cKO`. Feature columns are suffixed `_WT` / `_HP1gKO`, never `_cKO`.
10. **Colours** come from `../shared/common.R` (`disc_pal`, `scale_fill_disc`,
    `scale_fill_heat0_div`, `heat_col_fun`, `geno_colors`, `repeat_palette`,
    `theme_m`). No raw hex. Save with `save_fig(plot, name, subdir=)`.

## The ChIP tracks at 20 kb are noise-dominated (2026-10-16)

The genome-wide medians are sound. **Local per-bin structure is not**, and it has
already been misread once.

The bigWigs are input-scaled, so per-bin coverage is 0.1-3, not counts. GAM_1
multiplies by `gam_count_scale` (10) and rounds so DESeq2 has integers, which
means DESeq2 sees roughly **1-11 per 20 kb bin**, and `results()` returns an
**unshrunk** `log2FoldChange`. Measured per-bin sd: 0.14-0.16 for the broad
marks, up to **0.57 for H3K4me3** in its lowest-coverage quartile. 84% of
H3K4me3 bins have `baseMean < 0.5`; 30% of H3K9me3 bins and 17% of H3K36me3.

What this cost: H3K36me3 appeared to be **increased** over the Pcdh cluster. It
is not. Pcdh median is -0.049 against a genome median of -0.052, p = 0.93 versus
the rest of the genome, and the apparent bulge at chr18:37.3-37.8 Mb is -0.061,
p = 0.55. **Nothing at Pcdh differs from genome background at 20 kb, H3K4me3
included** - the peak-level H3K4me3 gain is a promoter-scale effect that a 20 kb
bin dilutes away. Use the ChIP peak-level results for anything focal.

**`lfcShrink` is not the fix.** ashr collapses H3K36me3 and H3K9me2 to exactly
0.000 (sd 0.000) and pulls H4K20me3 from -0.497 to -0.212, because a shrinkage
prior assumes most effects are null - the same assumption as median-of-ratios
normalisation, which `sizeFactors = 1` exists here to avoid. Do not re-test this.

The fix is representational, and GAM_12 implements it: aggregate the DISPLAY to
`gam12_mark_binsize` (100 kb), coverage-weighted; floor the panel y-span
(`gam12_mark_minspan`) so a narrow dynamic range is not stretched to full height;
and set opacity from `baseMean` against the mark's q90 so regions with no signal
visibly fade. Also note the track FILL diverges from each mark's **genome
median**, not from zero - blue means "above that median", which for a negative
median is still a loss.

## A/B compartments: axis weakly determined; NO strength change (2026-10-17)

The compartment call runs and is not artefactual, but it is **not reliable on
most chromosomes**. Do not report an A/B assignment, a net A<->B shift, or
`compartment_switch` without restating the quality numbers below.

**What was wrong and is now fixed.** The first draft ranked candidate
eigenvectors purely by correlation with an active-chromatin reference. That
picks the wrong vector: the LEADING eigenvector is systematically degenerate
(`frac_A` 0.04-0.13 - nearly every bin on one side) yet often correlates
slightly *better* with active marks than the real compartment vector. Balance is
now a precondition over 4 eigenvectors, not a tiebreak. Result: implausible
`frac_A` 13/38 -> 0/38, PC1 chosen 2/38 instead of 19/38.

**What is still wrong.** Usable axis (|rho| > 0.4) only **11 of 38** chromosome x
genotype combinations; median |rho| 0.348. The switch table is near-symmetric
(A->B 21,393 vs B->A 21,078); **do not read it as "no net shift"** - 36% of bins
switching is what axis noise looks like.

**Usable without reservation** (none depend on the PC): insulation,
`npmi_short` / `npmi_long`, `decay_slope`, and the raw NPMI.

### Settled 2026-10-17 - do not re-derive any of this

1. **Resolution is NOT the limitation. The old lead 1 is dead.** Tested 100 kb,
   240 kb, 500 kb and 1 Mb on chr1/7/11/16/19. |rho| is flat across the whole
   range (chr11 0.465 -> 0.425, chr19 0.250 -> 0.306, chr16 0.469 -> 0.548).
   Beagrie 2017 validated GAM compartments at 1 Mb, so this was the obvious fix
   and it does nothing. **The comment in `config.R` attributing the weak axis to
   resolution is wrong** and should not be used as justification.
2. **The real defect is that WT and cKO do not select the same eigenvector.**
   Mean squared cosine between the two genotypes' leading 4-eigenvector subspaces
   is only **0.27-0.57**, the 4th principal angle is near-orthogonal, and the
   selected PC differs between genotypes on 4 of 5 chromosomes tested. The
   eigengap lambda1/lambda2 is 1.2-3.1, so the leading vectors are close to an
   arbitrary rotation of each other. Any per-genotype PC1-vs-PC1 comparison is
   comparing two different objects.
3. **A shared axis is constructible and generalises.** Build it from the averaged
   correlation matrix, fit its direction once against the active-chromatin
   reference, score both genotypes on it. Held-out rho (weights fitted on odd
   coarse bins, evaluated on even) matches in-sample: chr16 0.555 -> 0.565,
   chr11 0.519 -> 0.553, chr1 0.435 -> 0.437. k = 4 and k = 8 agree. This is a
   real improvement over per-genotype selection **as a definition of the axis**.
4. **Never measure compartment STRENGTH along an eigen-derived axis.** It is
   circular: the axis fits whichever genotype's matrix built it. Same data, k=4,
   100 kb, median strength ratio cKO/WT by axis source - **WT axis 0.540
   (28/28 below 1), averaged axis 0.786 (25/28), cKO axis 1.207 (8/28)**. The
   direction of the answer is set by the choice of axis. An earlier
   "compartments are less defined in the cKO" claim came from the averaged axis
   and is **retracted**.
5. **There is no detectable compartment-strength change.** Measured on axes
   EXTERNAL to both contact matrices, so neither genotype can be favoured, all
   19 autosomes at 100 kb, bootstrap over chromosomes:
   - active-chromatin axis (WT H3K4me3 + H3K36me3 z): mean delta **-0.00018**,
     95% CI -0.00213 to +0.00186, negative on 10/19, paired Wilcoxon p = 0.62
   - TSS-density axis (no ChIP at all): mean delta +0.00019, CI -0.00057 to
     +0.00103, negative on 11/19, p = 0.74

   Two independent rulers, both null. Report this as a negative result, not as
   "not yet tested".
6. **lambda1 falls ~20% in the cKO** (ratio 0.58-1.02, median 0.81). That is a
   global property of the correlation matrix affecting every direction equally,
   not a compartment-specific effect, and it is what made the biased test look
   positive. Normalising by each genotype's own lambda1 on the averaged axis
   gives a strength ratio of 0.953 - i.e. nothing.
7. **Detection frequency explains chr7 and chr11, and only those.** A leading
   eigenvector tracks per-bin `n_valid` at |rho| 0.34-0.47 on chr7 and chr11,
   versus 0.01-0.23 on chr1, chr16 and chr19. Those are exactly the two
   chromosomes showing whole-chromosome contact offsets, so treat them as
   coverage-affected rather than biologically opposite.
8. **Our strength statistic is additive and is NOT comparable to published
   values.** The field computes saddle strength as (AA+BB)/(AB+BA) on
   **observed/expected**, which is positive so the ratio is defined. NPMI is
   already log-scale and goes negative, so GAM_2 subtracts the distance
   expectation rather than dividing; on O-E the AB denominator crosses zero and
   the ratio returns nonsense (it was producing negative "strengths"). Use
   **(AA+BB)/2 - AB**. Never quote it against a Hi-C paper's number.

### How the field does this, for when it is picked up again

- **dcHiC** is the standard differential-compartment tool: SVD per sample, PC
  selected and sign-phased against GC content or gene density, **quantile
  normalisation of the compartment scores across all samples** to put them on a
  common scale, then a per-bin Mahalanobis distance with a chi-square p and
  IHW/FDR. Note it assumes the same PC is selected in every sample - an
  assumption **our data violates** (point 2), which is why it cannot simply be
  run as-is here.
- **GENOVA / Imakaev saddle** for strength: O/E matrix, subtract 1, eigen
  decompose, first eigenvector scaled by sqrt(eigenvalue); phase by correlating
  the arm-wise score with an active mark (H3K4me1) and flipping on a negative
  correlation. Strength from the top/bottom 20% of score-ranked quantile groups.
  Our phasing against H3K4me3+H3K36me3 is the same idea.
- **HiCDOC** is the newer alternative that does prediction plus differential
  analysis with replicates.

**Cache gotcha:** `gam2_params$version` must be bumped whenever the compartment
METHOD changes, not only its parameters - a cached per-chromosome reduction
carries the PC it was computed with. Currently version 2.

## The shared A/B axis (GAM_2b) supersedes pc1_delta for WT vs cKO

`pc1_WT` and `pc1_HP1gKO` are each genotype's OWN selected eigenvector, and
settled point 2 above says those are different objects. Now measured on all 19
autosomes rather than the 5 originally tested: **the two genotypes select the
same eigenvector on only 8 of 19**. `pc1_delta` is therefore a difference
between two different vectors and must not be used for a genotype comparison.

`GAM_2b_Shared_Axis.R` implements settled point 3. One axis from the AVERAGED
correlation matrix at 100 kb, phased once against `active_chromatin_ref()`, then
both genotypes scored on it as `s_g = C_g %*% a` (one power-iteration step, so
if `a` were C_g's own eigenvector the score returns that eigenvector back).

**It is measurably better, and the comparison is in the output rather than
asserted.** Against the same active-chromatin reference, autosomes only:

| axis | median rho | reaches rho > 0.4 |
|---|---|---|
| WT's own | 0.315 | 5 of 19, either genotype |
| cKO's own | 0.254 | |
| **shared** | **0.422** | **10 of 19** |

The shared axis beats BOTH per-genotype axes on **16 of 19** autosomes
(exceptions chr3, chr8, chr13). Per chromosome in
`gam2b_shared_axis_summary.tsv`.

**The track is SHAPE ONLY.** Each genotype's score is median-centred and scaled
to RMS 1 per chromosome. Deliberate: lambda1 falls ~20% in the cKO, but that is
global (settled point 6), and leaving it in the track would repaint the
retracted "compartments are weaker in the cKO" claim across the whole genome.
Amplitude is kept in `rayleigh_WT` / `rayleigh_HP1gKO` / `lambda1_ratio` in the
summary, where it reads as what it is. **This stage does not re-answer the
strength question** - that is settled and null on two external axes.

Still weak evidence: 10 of 19 is better, not good. Report the quality numbers
with any compartment claim, as before.

**`GAM_6`'s decompaction verdict now uses `shared_delta`.** The earlier table
called Protocadherin "DECOMPACTS (short-range loss + toward A)" with the A-shift
resting entirely on `pc1_delta`. `pc1_delta` is still computed and still in the
table, labelled "not a sound comparison", so the two can be contrasted - it just
no longer decides a verdict. The short-range-loss half was never affected.

## "Protocadherin decompacts" is a matched-background claim, NOT a ranking

`gam5_decompaction_verdict.tsv` says Protocadherin DECOMPACTS. That is a
statement about its own chromatin stratum - the null matches on ChromHMM state,
repeat-content quartile and detection-coverage tertile - and it is routinely
misread as "Pcdh changes most", including in earlier write-ups here. It does not.

`gam5_target_genomewide_context.tsv` (GAM_6 section 4b) exists to stop that
misreading and should be quoted beside the verdict. For Protocadherin, against
same-width contiguous windows anywhere in the genome:

| | Protocadherin |
|---|---|
| windows losing MORE short-range contact | **43.7%** |
| windows more A-ward | 15.1% |
| windows beating it on BOTH criteria | **5.9%** |

So it is mid-pack for contact loss and in the top sixth for the A-shift, with
roughly one window in seventeen beating it outright. Both readings are true;
they answer different questions. MELTRON agrees with the narrow one - the Pcdh
cluster scores BELOW its own size-matched null (D 0.224 vs 0.260, p = 0.66).

**Gotcha that produced a wrong version of this table once:** `target_sets` holds
ROW INDICES into `d` in its current order (row i == bins[i]). `setorder(d, ...)`
inside an analysis block silently repoints every target set at different bins,
and `d` is used by the sections below. Always scan on a copy. The symptom was a
reported Pcdh mean of -0.019 for a quantity whose value over its actual bins is
+0.802; the two are now cross-checked against the matched-background table's
`observed` column, which is computed independently in the same stage.

## MELTRON (GAM_2c): two scores per interval, at the published window

`GAM_2c_Meltron.R` ports the Pombo lab's MELTRON / MELTRONIC (Winick-Ng et al.
Nature 2021; github.com/pombo-lab/Meltron, github.com/DominikSzabo1/MELTRONIC).
Insulation score = mean contact value in the square `[v-s, v) x (v, v+s]`, over
11 square sizes 240-1040 kb; the two genotypes' pooled IS distributions compared
with one-sided KS tests.

**Two scores per interval, never one signed number.** MELTRONIC emits
`melting_score` AND `condensation_score`, each Bonferroni-corrected separately,
and this follows it. The two one-sided tests are **not complementary**: where the
two distributions cross, both reject and an interval genuinely carries both.
`net_score` is derived for display and is secondary. On our data the two never
fire together at the published cutoff (**0.0%** of intervals), which is printed
every run as a result rather than assumed.

**The window is 120 kb stepping one bin, and getting it wrong breaks the score.**
MELTRONIC's `120kb_genomic_intervals.bed` is 120 kb wide at a 40 kb step - 3 bins
at their resolution. We match the physical width (6 bins at 20 kb). An earlier
version here used 500 kb / 25 bins, which pooled 275 nested overlapping values
into a test assuming independence; `-log10(p_adj)` saturated and the published
`|score| > 5` cutoff appeared to call 69% of the genome. **The conclusion drawn
then - "the published cutoff does not transfer" - was wrong**, and was a
statement about the window, not the data. At the published window the cutoff
calls 26.0% melting and 26.1% condensing and is used directly.

**One deviation remains, and it is forced.** IS is additive here, not
`log2(square mean) - log2(chromosome mean)`: **only ~51% of our insulation
squares have a positive mean NPMI**, so the log form is undefined for half the
genome. We subtract instead - the same decision, for the same reason, as settled
point 8. Sign and ranking stay comparable to published work; the **magnitude does
not**. (20 kb viewpoints rather than 40 kb is the other difference, and is only a
grid choice.)

**Regional, not fine-grained.** `gam2c_meltron_spatial_coherence.tsv`: the net
score changes sign between ~15% of adjacent 120 kb windows, median run 480 kb,
longest 8 Mb. Melting and condensing regions are discrete multi-megabase blocks.

**Quote `p_empirical_D`, not the raw score, for region sets.** Every region set
is scored against size-matched random regions on the same chromosomes. Result:
hotspot domains and large blocks melt and condense in the direction their name
predicts - lost D 0.171 / 0.129 against null medians 0.032 / 0.031, gained 0.198
/ 0.163 against 0.032 / 0.033, all four at the 1/501 floor - by an instrument
that never touches the differential-contact calls they were built from. The
derepressed gene families, protocadherin included, sit BELOW their own nulls.

## Region plots: three variants, and the slide one carries the Pombo strips

`GAM_7_Figures.R` section A renders each region as `<name>` (full triangle) and
`<name>_band` (distance capped so the panel is 0.95 in tall). `GAM7_SLIDE=1`
renders a third, `<name>_slide`, for four regions at exactly 13.33 x 7.5 in.

**Slide mode is a re-layout, not a resize.** A full region figure is ~1.09:1 and
cannot be squashed into 16:9, so slide mode drops the two raw NPMI maps, cuts the
marks to H3K9me3 + H4K20me3, shortens the caption to one line, and then solves
the contact map's distance cap FROM the height left over. The map is therefore
clipped to a band by construction, and how thin that band is depends on what
else is in the stack.

**It also carries the insulation-score map strips** (`ismap = TRUE`, default on
in slide mode and off elsewhere): WT, cKO and their difference, read from
`meltron_is.rds`, drawn directly under the difference map - the Pombo lab's own
display, x = position, y = insulation square size. The strips and the map above
them are the same quantity at two readings, so they should agree; a strip that
contradicts the map above it is a bug, not a finding.

**All three strips share one fixed limit (`gam12_is_lim` 0.20), and the
difference must never get its own.** Drawn at 0.10 first, on the reasoning that
the difference has no regional offset and deserves a tighter range, the panel
lied: at chr18:38.10 Mb, IS_WT -0.161 and IS_cKO -0.127 are both solid blue on
the strip scale while their +0.034 difference came out a third of the way to
saturated red. Stacked panels are read by comparison, and two rulers in one
stack is a display bug however right each number is. Genome-wide, 4.8% of cells
with BOTH genotypes below -0.05 were drawn over two-thirds saturated in the
difference; on the shared limit that is 0.0%, and the difference still reaches
52% of the bar at its own 99th percentile.

**Every row height in these figures is absolute, including the chrome.** That is
deliberate - mixing `"in"` map rows with `"null"` track rows made the layout
non-deterministic - but it means a figure saved shorter than its content does
NOT shrink to fit: it overflows the canvas and the outer rows are cut off
silently. Two consequences:

1. `gam12_h_chrome_slide` must be **measured**, on `patchwork:::patchworkGrob(p)`
   and not on `plot_table(p)`. `plot_table` returns the inner table only and
   omits the figure title, the figure caption and the outer margins - 0.55 in of
   them. Measuring the wrong one is why every slide figure shipped without its
   title and caption until 2026-10-17.
2. The caption is **one line and ggplot never wraps it**. Keep the slide caption
   under ~250 characters or it is cut off mid-word.

Adding the three strips costs 3 x `gam12_h_ismap` out of a fixed total, and the
trade against the map is 3:1 - every 0.1 in of strip costs 0.3 in of map, which
is why the strips are 0.33 in rather than the 0.42 they were first drawn at.

## GAM_5 sections A and B use a domain definition section C disproves

`GAM_5_Heterochromatin.R` section C builds heterochromatin domains properly
(gap-merged, size-swept) and shows that the definition sections A and B use -
single unmerged fragments above a width cut, from a 126,197-fragment peak set of
mean width 3.5 kb - is wrong. A and B were left on the old definition
deliberately, so the change is made and reported on its own rather than silently
restating earlier numbers; the file header says so.

**Consequence for reading the outputs:** `gam10_*` tables come from the old
definition, `gam13_*` from the merged one, and the report presents results
against all three definitions wherever they exist. Do not quote a `gam10_*`
domain-level number as if it were the current definition. This is the largest
outstanding inconsistency in the GAM tree.

## The `.rds` handoffs (schemas)

**`results/rds/gam_bins_annotated.rds`** (GAM_1), a named list:
- `$bins` GRanges (UCSC, autosomes, 123,133, `gam_bins()` order), mcol `bin_id`.
- `$stats` long data.table: `bin_id, mark (factor, mark_levels), baseMean,
  log2FoldChange, pvalue, padj`.
- `$l2fc` matrix [bins x 5 marks]; `$signal` matrix [bins x samples], un-rounded.
- `$annot` data.table keyed `bin_id`: `chr, start, end, n_contacts_lost,
  n_contacts_gained, frac_LINE/SINE/LTR, hmm_<state>, tad_boundary, n_tss,
  compartment_abs, compartment_centred` (levels `stable`, `H4K20me3_only`,
  `co_loss`, NA elsewhere).
- `$meta` parameters used.
- Compartment cutoffs are MINUTE_3's verbatim (`H4_LOST_CUT = -0.5`,
  `K9_LOSS_CUT = -0.3`, `UNCH_CUT = 0.15`) applied per bin. **`compartment_centred`
  (each mark minus its own genome-wide median) is the primary label downstream**;
  the absolute cutoffs are not shift-invariant against a genome-wide loss
  (`../ChIP/results/COMPARTMENT_CUTOFFS.md`).

**`results/rds/npmi_bin_features.rds`** (GAM_2), a named list:
- `$features` data.table keyed `bin_id` (autosomes plus X/Y where available):
  `bin_id, chr, start, end`, then per genotype suffixed `_WT` / `_HP1gKO`:
  `pc1, compartment (A/B), insulation, npmi_short, npmi_long, decay_slope,
  n_valid`; plus `pc1_delta` (HP1gKO minus WT), `insulation_delta`,
  `compartment_switch` (`A_A, A_B, B_A, B_B`).
- `$insulation_windows`, `$saddle` (per genotype), `$meta`.
- PC1 is oriented **per chromosome** so A (positive) is the WT H3K4me3-high side.
  The sign is arbitrary per chromosome, so an unoriented PC1 makes WT vs KO
  meaningless. No silent fallback if GAM_1's rds is absent: `stop()`.

## results/ is git-ignored and regenerable

Every output, the `.rds` files included, is regenerable by `run_GAM.R`. The only
committed files under `results/` are the hand-written `.md` files
(`METHODS.md`). The input data (`NPMI_matrices/`, `Diff_contacts/`, `Hotspots/`,
about 7.7 GB) is git-ignored too. `GAM_RESULTS=<dir>` writes a run elsewhere.

## Do not invent provenance

Number of nuclear profiles, the `passed_qc.curated` criteria, what sigma 2.6 was
fit on, the protocol variant and whether the cKO is the same Emx1-Cre cortex as
the ChIP are unknown. They are open questions in `results/METHODS.md`; ask rather
than assume.
