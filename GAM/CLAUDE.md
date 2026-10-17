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

## MELTRON (GAM_2c): use the KS statistic D, not the published melting score

`GAM_2c_Meltron.R` ports the Pombo lab's MELTRON / MELTRONIC (Winick-Ng et al.
Nature 2021, github.com/pombo-lab/Meltron). Insulation score = mean contact
value in the square `[v-s, v) x (v, v+s]`, over 11 square sizes 240-1040 kb;
the two genotypes' pooled IS distributions compared with a one-sided KS test;
published score = -log10(Bonferroni-corrected p).

**Two deviations, both forced, both in the file header:**

1. **IS is additive here, not `log2(square mean) - log2(chromosome mean)`.** The
   published form needs a positive quantity. **Only ~50% of our insulation
   squares have a positive mean NPMI**, so the log form is undefined for half
   the genome. We subtract instead - the same decision, for the same reason, as
   settled point 8. Sign and ranking stay comparable to published work; the
   **magnitude does not**.
2. 20 kb viewpoints rather than 40 kb, so the output joins the rest of the repo
   without a regrid. Square sizes unchanged.

**The published |score| > 5 cutoff does not transfer and is NOT used.** It calls
about a third of the genome in the published setting; here it calls roughly two
thirds, because each window pools 11 square sizes x 25 viewpoints into a test
that assumes independence - they are nested and overlapping - so
`-log10(p_adj)` saturates. Calls come from the tails of the genome's own
**signed KS D** instead, which is a relative statement ("the most melted 5%")
and is labelled as such. The published figure is printed beside ours so the
discrepancy stays visible.

**The track is regionally coherent over megabases, not at window level.**
Measured every run into `gam2c_meltron_spatial_coherence.tsv`: aggregated to one
point per 500 kb sliding window, the direction changes sign between ~25% of
adjacent windows, the median run of one direction is 3 windows (1.5 Mb) and the
longest is 28 (14 Mb). That is the scale of GAM_4's one-directional blocks - and
those blocks come out melting (lost) and condensing (gained) at the floor of the
null here. So read regions, not single windows.

**Quote `p_empirical_D`, not the raw score.** Every region set is also scored
against size-matched random regions on the same chromosomes, and that null is
what makes the output interpretable at all: **the raw score saturates at the
+/-300 cap for all four significant sets**, so it cannot rank or separate them.
Against the null they separate cleanly. D also shows something the score cannot:
the gene-family sets sit BELOW their own nulls (e.g. Pcdh 0.224 against a null
median of 0.260), i.e. less changed than a random region of the same width, which
a one-sided score reports simply as "not significant".

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
