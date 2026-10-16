# GAM — Genome Architecture Mapping, HP1γ cKO vs wild type

Analysis of **Genome Architecture Mapping (GAM)** data comparing **HP1γ (Cbx3)
conditional knockout vs wild type**, at **20 kb** resolution on mouse **mm39**,
and its integration with the MINUTE ChIP-seq results in [`../ChIP/`](../ChIP/).

GAM here supplies three delivered products: dense **NPMI** contact matrices per
genotype, **differential contact** pairs, and per-bin **hotspot** counts of those
differential contacts. The R pipeline projects the ChIP changes onto the same
20 kb grid and asks how 3D contact remodelling relates to the heterochromatin
changes.

> **Status:** data inventoried and the pipeline scaffolded. No GAM results have
> been generated or interpreted yet; this README describes what each stage
> produces, not what it found.

---

## Data inventory

| Folder | Contents | Files | Size | Chromosomes |
|--------|----------|-------|------|-------------|
| `NPMI_matrices/WT_at20Kb/` | dense square NPMI matrices, WT | 21 | 2.7 GB | chr1-19, chrX, chrY |
| `NPMI_matrices/cKO_at20Kb/` | dense square NPMI matrices, cKO | 21 | 2.3 GB | chr1-19, chrX, chrY |
| `Diff_contacts/` | significant differential contact pairs (`.bg2.gz`) | 38 | | autosomes only (19 x 2 directions) |
| `Hotspots/` | per-bin differential contact counts (`.bg.gz`) | 38 | | autosomes only (19 x 2 directions) |

Total about 7.7 GB, **not in git** (see *Data dependency*).

Parameters encoded in the filenames (all set in `config.R`): mm39, bin size
20 kb (`segregation_at20000`), `sigma.2.6`, segregation date `2026-02-23`,
`passed_qc.curated`, prefix `HP1g_mm39`.

### NPMI matrices

`NPMI_matrices/{WT,cKO}_at20Kb/HP1g_mm39.segregation_at20000.2026-02-23.passed_qc.curated.sigma.2.6.NPMI.<chr>.squared.txt.gz`

- Dense **square** text matrices of normalised pointwise mutual information
  between 20 kb bins of one chromosome. chr19 is 3072 x 3072; chr1 is about
  9.7k x 9.7k (about 760 MB as doubles).
- Row/column labels are bin ids of the form `chr19:0-20000` (0-based half-open).
- **The filenames are identical in the two genotype folders — genotype is encoded
  only by the directory name** (`WT_at20Kb` vs `cKO_at20Kb`).
- Autosomes **plus chrX and chrY**.

### Differential contacts (`.bg2.gz`)

`Diff_contacts/significant_WT_vs_cKO_<chr>_<sel>_maxdistinf_z2_npmi0.1.bg2.gz`,
with `<sel>` = `selWT` or `selcKO`. Columns:

```
chrom_x, start_x, end_x, chrom_y, start_y, end_y, bindist,
set1, set2, comparison, set1_npmi, set2_npmi, maxnpmi,
set1_zscore, set2_zscore, deltaz, abs_deltaz
```

(`set1` is WT and `set2` is cKO in these `WT_vs_cKO` files; `deltaz` is the
z-score difference.) Filters baked into the files: **z > 2, NPMI > 0.1,
maxdist = inf**. Observed `bindist` spans 20 kb to 57.18 Mb. Totals:
**selWT = 25,558,790 pairs; selcKO = 25,042,198 pairs.**

### Hotspots (`.bg.gz`)

`Hotspots/significant_WT_vs_cKO_<chr>_<sel>_maxdistinf_z2_npmi0.1_dist0-5000000.bg.gz`.
Columns `chrom, start, end, contact_count`, where `contact_count` is the number of
significant differential contacts within **0-5 Mb** of that bin. About 114,793
bins (selWT) and 114,651 (selcKO), counts 1-162.

**These are dense continuous tracks, not sparse hit lists.** The autosomal 20 kb
grid has 123,133 bins, so about 93 % of bins carry at least one differential
contact. The signal is the count *magnitude*; bins absent from a file have zero.

### Direction semantics

| File tag | Meaning | In the KO |
|----------|---------|-----------|
| `selWT` | contact stronger in WT | **lost** |
| `selcKO` | contact stronger in cKO | **gained** |

Downstream columns are named from the KO's point of view
(`n_contacts_lost`, `n_contacts_gained`). Getting this backwards inverts the
headline result.

---

## Provenance and replication

- **Replication:** two biological replicates per genotype, **already merged** into
  the single segregation table per genotype that is on disk. The effective
  biological n is therefore 2 vs 2 — the same as the ChIP side (6 vs 4 libraries
  that are nested technical replicates of 2 biological replicates per genotype),
  so the assays are replication-matched.
- The delivered files carry **no replicate dimension**. The baked-in z > 2 is the
  only significance available from them, so GAM-side statistical claims come from
  matched-background permutation on the bin side, not from a p-value attached to
  a contact.
- Everything not recoverable from the filenames (nuclear profile counts, QC
  criteria, what sigma was fit on, protocol) is listed as an open question in
  [`results/METHODS.md`](results/METHODS.md).

---

## Pipeline

```
config.R                          parameters, path builders, loaders (edit this only)
  │
GAM_1_Bins_and_Signal.R           ChIP signal/stats + annotation on the 20 kb grid
  │   writes results/rds/gam_bins_annotated.rds
  ▼
GAM_2_NPMI_Compartments.R         NPMI -> PC1 compartments, insulation, decay
  │   writes results/rds/npmi_bin_features.rds
  ├─▶ GAM_3_Contact_Coupling.R
  ├─▶ GAM_4_Compartments_vs_Contacts.R
  └─▶ GAM_5_HUSH_and_Repeats_Contacts.R
```

| Stage | Inputs | Outputs |
|-------|--------|---------|
| **GAM_1** `Bins_and_Signal` | ChIP bigWigs (`MINUTE_DATA`), the ChIP annotated `.rds` (`require_chip_rds()`), Hotspots, repeat/TAD/ChromHMM annotation | `rds/gam_bins_annotated.rds`: the 123,133-bin autosomal grid with per-bin ChIP DESeq2 stats and signal, per-bin differential contact counts (`n_contacts_lost`/`n_contacts_gained`), repeat fractions, ChromHMM coverage, TAD boundary, TSS count, and ChIP compartment labels (`compartment_abs`, `compartment_centred`) |
| **GAM_2** `NPMI_Compartments` | NPMI matrices (both genotypes, autosomes + chrX/chrY), GAM_1 rds | `rds/npmi_bin_features.rds`: per bin and genotype PC1, A/B compartment, insulation, short/long-range NPMI, decay slope, valid-cell count; PC1/insulation deltas and A/B switch class; saddle matrices |
| **GAM_3** `Contact_Coupling` | GAM_1 + GAM_2 rds | how per-bin contact remodelling relates to ChIP change (tables/figures) |
| **GAM_4** `Compartments_vs_Contacts` | GAM_1 + GAM_2 rds | ChIP-defined compartments vs GAM compartments vs contact remodelling |
| **GAM_5** `HUSH_and_Repeats_Contacts` | GAM_1 + GAM_2 rds | HUSH/CBX3 targets and repeat content vs contact remodelling, using size-matched permutation nulls |

GAM_1 and GAM_2 are sustained all-core jobs; GAM_3-5 need only the two `.rds`
files and can be re-run independently. The two `.rds` schemas are summarised in
[`CLAUDE.md`](CLAUDE.md). PC1 is oriented per chromosome so that A (positive) is
the H3K4me3-high side, which is why GAM_2 depends on GAM_1.

---

## Requirements and usage

R (>= 4.2). Packages install themselves on first run via `../shared/setup.R`
(set `MINUTE_AUTO_INSTALL=0` to disable and get an error instead). GAM adds one
dependency of its own, **`R.utils`** (`data.table::fread` needs it to read the
gzipped inputs); `config.R` registers it.

The ChIP pipeline must have been run first:

```bash
cd ChIP && Rscript run_MINUTE.R     # produces the annotated peak .rds + uses the bigWigs
cd ../GAM
Rscript run_GAM.R                   # all stages
```

GAM_1 and GAM_2 are heavy local jobs; run them under the repo's
`safe-local-alignment` guidance (thermal/swap protection on Apple Silicon).
A single stage, once its inputs exist:

```bash
cd GAM
Rscript GAM_3_Contact_Coupling.R
```

### Environment variables

| Variable | Default | Purpose |
|----------|---------|---------|
| `GAM_DATA` | `.` | folder holding `NPMI_matrices/`, `Diff_contacts/`, `Hotspots/` |
| `GAM_RESULTS` | `results` | output tree (use a separate one for a sensitivity run) |
| `CHIP_RESULTS` | `../ChIP/results` | ChIP results tree providing the annotated peak `.rds` |
| `MINUTE_DATA` | `../ChIP/data` | ChIP bigWigs, read by GAM_1 |
| `MINUTE_AUTO_INSTALL` | on | set `0` to disable package auto-install |

---

## Configuration

**All parameters live in `config.R`**: bin size, sigma, segregation date, z and
NPMI cutoffs, hotspot window (5 Mb), the genotype bridge (`WT`/`HP1gKO` to
`WT`/`cKO`), chromosome sets, the `lost`/`gained` mapping, path builders
(`npmi_file`, `diff_file`, `hotspot_file`) and loaders (`load_hotspots`,
`load_diff_contacts`, `load_npmi`). Stages must not redefine them.

---

## Outputs

Everything is written under `results/` (`counts/ rds/ tables/ figures/ bed/`),
which is **git-ignored in full and regenerable** by `run_GAM.R`. The exceptions
are the hand-written top-level `.md` files (`results/METHODS.md`), which are
committed. Figures are saved as both `.png` and `.pdf` into per-analysis
subfolders by `save_fig()`.

---

## Data dependency (NOT in git)

The GAM data drop (about 7.7 GB: `NPMI_matrices/`, `Diff_contacts/`, `Hotspots/`)
is git-ignored. Place it in `GAM/` or point `GAM_DATA` at wherever it is staged.
Sanity check:

```bash
ls GAM/NPMI_matrices/WT_at20Kb | wc -l     # expect 21
ls GAM/NPMI_matrices/cKO_at20Kb | wc -l    # expect 21
ls GAM/Diff_contacts | wc -l               # expect 38
ls GAM/Hotspots | wc -l                    # expect 38
```
