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

## OPEN ISSUE — the A/B compartment call is not yet trustworthy

**Do not report an absolute A/B assignment, or a net A→B compartment shift, from
the current GAM_2 code.** Two coupled problems, both unresolved:

1. **PC1 is median-centred per chromosome per genotype**, so `frac_A ≈ 0.5` by
   construction and a chromosome-wide A/B imbalance cannot be represented. That
   in turn constrains GAM_4's switch table to near-symmetric A→B and B→A counts
   whatever the biology does. GAM_2's chr19 smoke test shows the signature
   exactly: A→B 609 vs B→A 612. The code comment calls this "fine for the
   WT-vs-KO question" — that is the one claim in this stage that is wrong, since
   the WT-vs-KO question in GAM_4 *is* the net shift.
2. **The leading eigenvector may not be the compartment eigenvector.** GAM_2
   computes two (`eigs_sym(k = 2)`) but always takes the first, and the raw PC1
   on chr19 sat ~2 SD above zero with 97% of bins on one side — the signature of
   a coverage/mean component rather than an A/B split. Standard practice is to
   select the eigenvector that best tracks an independent A/B marker, which here
   would be the H3K4me3 reference already loaded for orientation.

Centring was a reasonable response to (2), but it masks (2) rather than fixing
it, and it costs (1).

**What IS usable now:** PC1 as a within-chromosome *ranking*, the WT-vs-KO PC1
*delta* at matched rank, insulation, the distance bands and the decay slope.
**What is not:** `frac_A`, `compartment_<g>`, `compartment_switch`, and any
claim of the form "the KO moves territory from A to B".

**To resolve, on the first real run:** record each eigenvector's correlation with
WT H3K4me3 per chromosome and see whether PC1 or PC2 tracks it; then set the A/B
boundary from that marker rather than from PC1's own median. Until then treat
GAM_4's compartment-switch section as diagnostic output, not a result.

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
