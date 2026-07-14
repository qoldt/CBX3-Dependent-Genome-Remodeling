# CBX3-Dependent Genome Remodeling

Analysis code for a study of **HP1γ (CBX3)** in organizing chromatin in neurons,
comparing HP1γ knockout vs wild-type across three complementary assays.

## Structure

| Folder | Assay | Status |
|--------|-------|--------|
| [`ChIP/`](ChIP/) | Multiplexed ChIP-seq (MINUTE pipeline) — differential histone-mark analysis over called peaks | Active |
| `scRNAseq/` | Single-cell RNA-seq | Planned |
| `GAM/` | Genome Architecture Mapping | Planned |

Each folder is self-contained, with its own README describing inputs, usage, and
outputs. Genome build is mouse **mm39** throughout.

## Data

Large raw inputs (bigWigs, RepeatMasker annotation, etc.) are **not** stored in
this repo — see the per-folder README for what to provide locally and where to
get it.
