# GAM — Genome Architecture Mapping (mm39, 20 kb)

HP1γ cKO vs wild-type GAM, integrated with the MINUTE ChIP-seq results. Read the
root `../CLAUDE.md` first for the cross-assay conventions; this file holds only
what is specific to GAM. Data inventory and formats: `README.md`. Unresolved
provenance: `results/METHODS.md`.

## Pipeline

Run from the **`GAM/` directory** (no `setwd()`):

```
config.R                              # parameters, path builders, loaders (single source of truth)
GAM_1_Bins_and_Signal.R               # ChIP signal/stats + annotation on the 20 kb grid  (heavy)
  │  writes results/rds/gam_bins_annotated.rds
GAM_2_NPMI_Compartments.R             # NPMI -> PC1/compartments/insulation/decay           (heavy)
  │  needs GAM_1 (PC1 orientation); writes results/rds/npmi_bin_features.rds
  ├── GAM_3_Contact_Coupling.R
  ├── GAM_4_Compartments_vs_Contacts.R
  └── GAM_5_HUSH_and_Repeats_Contacts.R   # 3-5 read only the two .rds; independent of each other
```

`Rscript run_GAM.R` runs them in order. It also needs the ChIP pipeline to have
run (`../ChIP/results/rds/annotated_results_*.rds`; `require_chip_rds()` stops
with an actionable message otherwise). GAM_1 and GAM_2 are sustained all-core
jobs: run them under the `safe-local-alignment` skill. GAM_3-5 should `stop()`
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

## A/B compartments: how the axis is defined (resolved in code, unverified on data)

GAM_2 does **not** take the leading eigenvector as the compartment axis, and does
**not** median-centre it. Both of those were wrong in the first draft and the
reasoning is worth keeping, because the failure was invisible:

- Taking the top eigenvector unconditionally gave, on chr19, a vector with 97% of
  bins on one side of zero. That is not an A/B split; it is the overall contact
  structure. In Hi-C-like data the compartment signal sits in PC1 **or** PC2
  depending on the chromosome, which is why cooltools-style pipelines select by
  correlation with GC or gene density.
- Median-centring the result forced `frac_A = 0.5` per chromosome per genotype,
  which in turn forced GAM_4's switch table to near-symmetric A→B / B→A counts
  **whatever the biology did** (chr19 draft: 609 vs 612). It masked the first
  problem and made the net-shift question unanswerable by construction.

What it does instead: `select_compartment_pc()` computes two eigenvectors and
picks whichever correlates better (|Spearman|) with an **active-chromatin
reference**, then orients the sign by the same reference. The reference is the
mean of the per-bin z-scores of **WT H3K4me3 and WT H3K36me3**.

Why those two marks:
- They are complementary at this resolution. H3K4me3 is punctate at promoters and
  thin across a 20 kb bin; H3K36me3 is broad over transcribed gene bodies and
  better matched to the bin size. Their average is steadier than either alone.
- They are the marks that **do not change** in the KO (per-bin median log2FC
  +0.000 and −0.052, against −0.497 for H4K20me3). The A/B axis is therefore the
  same ruler in both genotypes. Anchoring on H3K9me3/H4K20me3 would be circular —
  the axis would move with the effect under study.

With the right eigenvector selected, no re-centring is needed and `frac_A` is free
to differ from 0.5, which is what makes a net A↔B shift detectable.

**Still to confirm on the first real run** (the numbers are in
`tables/npmi_compartment_orientation.tsv`): how often `which_pc == 2` — if that is
common, say so in the methods rather than leaving it implicit; that `rho_pc1` vs
`rho_pc2` separate cleanly, since a near-tie means the axis is poorly determined
for that chromosome; and that WT and KO PC correlate strongly genome-wide, because
a wholesale reorganisation is a much bigger claim than this design can support.
Until those are checked, treat `compartment_switch` as diagnostic.

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
