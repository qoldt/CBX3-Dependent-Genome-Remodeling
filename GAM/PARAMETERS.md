# GAM parameters — why each value is what it is

`config.R` holds the values and `utils.R` the code; this file holds the
reasoning that used to sit in `config.R` as comment paragraphs. It was moved
here on 2026-11-23 without changing a value. If you change a parameter, update
its section here in the same commit.

## Input file parameters

`gam_binsize`, `gam_sigma`, `gam_seg_date`, `gam_seg_prefix`, `gam_z`,
`gam_npmi`, `gam_maxdist_label` and `gam_hotspot_maxdist` match the file names
on disk. Changing one changes which files the path builders resolve to, so a
re-run at a different sigma / z / NPMI cutoff is a config edit, not a code
edit. `|delta-z| > 2` and NPMI > 0.1 are baked into the delivered files; they
are not ours to choose.

## Genotype labels and direction

Canonical repo-wide labels are WT / HP1gKO (see `../CLAUDE.md`). The GAM files
say "cKO" for HP1gKO. Everything downstream of `config.R` uses the canonical
name; only the path builders ever see "cKO".

The NPMI filenames are **identical** between the two genotype folders — genotype
is encoded only by the directory name (`WT_at20Kb` / `cKO_at20Kb`). A matrix file
copied out of its folder has lost its genotype, so never resolve one by basename.

The upstream caller labels each significant differential contact by the
genotype in which it is SELECTED (i.e. stronger):

    selWT  -> contact stronger in WT   -> LOST   in HP1gKO
    selcKO -> contact stronger in cKO  -> GAINED in HP1gKO

Because the whole study is framed as "what does losing HP1g do?", every
downstream column is named from the KO's point of view: `n_contacts_lost` /
`n_contacts_gained`. This mapping is the single most reversible mistake in the
GAM integration — getting it backwards inverts the headline result.

## Chromosomes

The differential products (`Diff_contacts`, `Hotspots`) are autosome-only. The
NPMI matrices additionally carry chrX and chrY. Keep cross-assay integration on
the autosomes; X/Y are available for NPMI-only analyses.

## Inputs from the ChIP assay

GAM_1 puts ChIP signal on the GAM bin grid, so it needs the ChIP bigWigs, the
ChIP annotation BEDs and the ChIP sample sheet.

**Why not `source("../ChIP/config.R")`:** both configs define `results_dir` /
`counts_dir` / `rds_dir` / `fig_dir` / `tables_dir` by the same names, so
sourcing the ChIP one would silently redirect every GAM output into
`ChIP/results` (and re-attach the whole ChIP library set). The few things GAM
needs are therefore re-derived, and they must track `ChIP/config.R`: same
`MINUTE_DATA` root, same sample sheet, same default exclusion.

The consensus peak BEDs are committed in the ChIP tree (small), unlike the
bigWigs and annotation. GAM_5 (old GAM_10) reads the merged H3K9me3/H4K20me3
set from there to define heterochromatin domains.

**Sample exclusion.** By default HP1gKO rep2 is excluded, exactly as on the ChIP
side — it is the highest-signal library in all five marks and discordant with
rep3, its own technical twin, so it is treated as a technical failure. Keeping
the same default is what makes the per-bin log2FC comparable to the published
peak-level numbers; `MINUTE_EXCLUDE=none` reproduces the all-library run.
Syntax: comma-separated sample_ids and/or `<genotype>:<replicate>` pairs.

## Annotation

**LADs.** Peric-Hupkes 2010 (GSE17051), lifted mm9 → mm39 by
`prepare_LAD_annotation.R`. Constitutive sets: cLAD is LAD in all four profiled
cell types, ciLAD is LAD in none. They are deliberately **not** neuronal — they
are the cell-type-invariant core, which is what makes them an axis independent
of our own marks and of either genotype. Per-cell-type (NPC, astrocyte) calls
exist in the same GEO entry if cell-type specificity ever matters.

**CTCF.** ENCODE conservative IDR thresholded peaks, lifted mm10 → mm39 by
`prepare_CTCF_annotation.R`. The third, orthogonal line of evidence on where a
boundary is: the two TAD partitions are a published boundary set and our own
insulation minima, neither of which is an independent measurement of the
architectural protein itself.

**Bulk tissue, not sorted CA1.** ENCSR877MSN is layer of hippocampus and
ENCSR644VYX is left cerebral cortex, both whole tissue from adult B6CASTF1/J.
No CA1-pyramidal-specific CTCF ChIP-seq is public. CTCF occupancy is largely
cell-type-invariant, which is what makes bulk usable, but label it bulk wherever
it appears — it is not a CA1 track. Both tissues are carried because together
they match the two assays: GAM is CA1 pyramidal cells (hippocampus), the ChIP
is Emx1-Cre cortex + hippocampus.

**TAD partitions.** Two independent partitions, kept side by side on purpose.
Any "this locus talks to distal TADs differently" claim depends on where the
boundaries are, so a result that only holds under one partition is a result
about the partition.

- **External** — published boundaries, fixed across genotypes. Non-circular,
  but a different cell type and blind to anything the KO does.
- **Own** — insulation minima from our own GAM, per genotype
  (`gam_tad_insul`). Self-consistent with the contact data and free to move
  between genotypes, but it is our own call, so it cannot also be evidence
  about boundary strength.

## ChIP quantification on the GAM grid

**`gam_count_scale` (10).** DESeq2 needs integers, but mean coverage over a
20 kb bin is often below 1, so plain rounding would collapse most bins onto 0/1.
Scaling by 10 before rounding preserves the gradation; the factor cancels in
log2FC, and baseMean is divided back out so it stays on the coverage scale.
Same value and reasoning as `ChIP/MINUTE_7_H3K36me3_Genes.R`.

**`gam_cores`.** Leaves one core free; GAM_1 is a sustained all-core job on a
laptop, so lower it (`GAM_CORES=4`) when running alongside anything else.

**Mark-loss group cutoffs (`H4_LOST_CUT`, `K9_LOSS_CUT`, `UNCH_CUT`).** MINUTE_3's,
verbatim. Do not retune them here. They define the stable / H4K20me3_only /
co_loss groups, and the whole point of reusing them is that the GAM bins are
labelled by the same rule as the ChIP peaks. They are absolute and therefore
not shift-invariant against a genome-wide loss, which is why every GAM stage
treats the median-centred variant as primary
(`ChIP/results/COMPARTMENT_CUTOFFS.md`). These are not A/B compartments.

## Heterochromatin domains

`het_domain_bed` is 126,197 fragments of mean width 3.5 kb (444 Mb total). It is
a merged **peak** set, not a domain set: GAM_10 `reduce()`s it and keeps
whatever single fragment happens to exceed a width cut, which is not the same
thing as a domain. Gap-merging first is what turns it into one. Every
domain-level result is reported against the unmerged set too, as the hotspot
domains already are.

`gam_het_merge` was chosen from the sweep in GAM_13. The fragments percolate:
genome coverage runs 1.3% (gap 0) → 4.5% (10 kb) → 17.5% (20 kb) → 64.8%
(50 kb) at min_width 100 kb, because once the gap exceeds the typical spacing
between peaks everything fuses into chromosome-length blocks. 10 kb sits below
that transition and yields 331 domains over 110 Mb, median 160 kb — genuinely
large domains rather than the whole genome.

Worth stating plainly: this BED contains very little large-domain structure. At
gap 0 only 33 Mb sits in pieces ≥ 100 kb, so the 444 Mb it covers is almost
entirely sub-50 kb peaks. GAM_10's "domains ≥ 100 kb" were whichever single
fragments happened to be long, not merged domains.

## Chromatin classes (2026-10-19)

Replaces the old binary is_het / eu_eu split, which called everything that was
not an H3K9me3/H4K20me3 peak "euchromatin". That was wrong twice over: it put
Polycomb (facultative heterochromatin) in the euchromatin bucket, and the bucket
was 85% ChromHMM-quiescent, i.e. unmarked chromatin, not euchromatin. Only about
one bin in nine of the old "euchromatin" was actually active.

    LEVEL 1  euchromatin | heterochromatin        - the headline contrast
    LEVEL 2  TSS, Enhancer, Tx | Polycomb, H3K9me3, Quiescent

Level 2 nests exactly inside level 1, so any level-1 number is the pooled
level-2 numbers and the two can never disagree.

**State → level-2 mapping, and why it is not just the state names.** Mostly the
18-state cortex P0 model's own names are taken at face value. Three groups are
assigned against their name, because our own WT ChIP contradicts it and the
group is large enough for that to matter. Mean WT signal in bins where the state
covers ≥ 20% of the bin (genome-wide mean is ~1.0–1.5 for each):

| state  | bins   | H3K4me3 | H3K36me3 | H3K9me3 | H4K20me3 | assigned |
|--------|-------:|--------:|---------:|--------:|---------:|----------|
| Het    |  1,372 | 0.70 | 0.79 | 4.49 | 3.54 | H3K9me3 |
| Quies3 |  5,331 | 0.39 | 0.68 | 3.03 | 2.45 | H3K9me3 ← not quiescent |
| Quies2 |  9,983 | 0.60 | 0.77 | 2.37 | 1.89 | H3K9me3 ← not quiescent |
| Quies  | 94,902 | 0.97 | 0.81 | 1.05 | 1.15 | Quiescent |
| QuiesG | 22,198 | 3.10 | 1.72 | 0.66 | 0.85 | Tx ← not quiescent |

Quies2 and Quies3 carry H3K9me3 at 2.4–3.0x, against 1.05x for Quies proper —
they are H3K9me3 heterochromatin wearing a quiescent label. QuiesG runs the
other way: H3K4me3 3.10 and H3K9me3 depleted at 0.66, which is gene-associated
chromatin. Leaving these three under "Quiescent" would have put 15,314 bins of
real heterochromatin and 22,198 bins of gene-adjacent chromatin in the wrong
level-1 class.

Four small states are kept on their names although our marks disagree — TssBiv
(H3K9me3 2.16) and EnhPr (H3K9me3 3.54) look heterochromatic, EnhPois and EnhLo
are quiet. All four are under 150 bins at ≥ 20% coverage, so reassigning them
would be fitting the map to a few hundred bins. Recorded so the inconsistency is
a known choice rather than an oversight.

**We do not measure H3K27me3.** The Polycomb class rests entirely on this
external segmentation, whose training marks are not recorded (METHODS open
question). Our own data can only confirm what Polycomb is not: ReprPC/ReprPCWk
bins carry the lowest H3K9me3 (0.37) and H4K20me3 (0.43) of any class, below
active chromatin, so they are certainly not constitutive heterochromatin.

**`gam_het_peak_cut` (0.5).** A bin whose H3K9me3/H4K20me3 peak coverage reaches
this is called H3K9me3 heterochromatin whatever ChromHMM says. This override is
not a hedge — it is the better-evidenced call. 8,266 bins that ChromHMM puts in
plain Quies carry mean WT H3K9me3 of 3.94 against 0.88 for the Quies bins
without a peak, a 4.5x separation measured in this tissue. ChromHMM's own Het
state is only 0.92% of the genome against our 7.2%, so it under-calls badly
here and is folded in rather than relied on.

## Distance bands

One definition, because these were copied into GAM_5 (twice), GAM_7 (twice) and
`make_report.R`, and a band set that disagrees between two tables is the kind of
thing no `stopifnot` catches.

The top band was ">10Mb" until 2026-10-19, open-ended over 10–195 Mb. That
pooled 364 M possible pairs, and nearly half of them were beyond 50 Mb, so
anything happening at 10–25 Mb was swamped by the very long range. Splitting at
50 Mb divides it almost evenly: 189 M pairs in 10–50 Mb against 174 M beyond
50 Mb.

**Read the >50Mb band per chromosome with care.** It is open-ended, so it means a
different separation on every chromosome: chr19 (61 Mb) can only contribute
50–61 Mb and has 0.16 M pairs, while chr1 (195 Mb) contributes 50–195 Mb with
26 M — a 160-fold spread. Pooled genome-wide that is fine; as one column of a
per-chromosome figure it is not comparable across rows, and chr19 is
effectively unpowered there. The same was true of the old ">10Mb", less
visibly.

## Hotspot domains

The hotspot counts are a dense, smooth track, not a sparse hit list, so there
are two legitimate units and they answer different questions:

- **per-bin count** — the continuous covariate. Highest resolution, no threshold
  chosen, right for correlation/coupling (GAM_3). Neighbouring bins are not
  independent (measured autocorrelation r = 0.59 at 20 kb, 0.40 at 100 kb, 0.26
  at 500 kb), so any CI over bins must resample chromosomes, not bins.
- **domain** — the independent unit. Right for enrichment/overlap tests and for
  "how much territory remodels", because counting adjacent bins as independent
  observations pseudo-replicates — the same reason `chromHMM_size_weighted()`
  permutes region labels rather than base pairs.

**Why smooth before thresholding.** Thresholding the raw track fragments it: at
the 90th percentile the median run is 1 bin and only 22% of hotspot bins fall in
domains ≥ 100 kb, because a single bin dipping below the cutoff splits a real
domain in two. Absorbing sub-threshold gaps fixes some of that (gap ≤ 2 → 62%)
but chains through long quiet stretches if pushed (gap ≤ 3 produces a 3.6 Mb
"domain"), and the gap size is arbitrary. A 100 kb (5-bin) running mean first
gives median 4 bins and 80% of hotspot bins in ≥ 100 kb domains without that
pathology. 100 kb is chosen to match the autocorrelation structure: still
r ≈ 0.40, so averaging over it is justified, whereas a 500 kb window (r ≈ 0.26)
would fuse genuinely separate domains.

**The unadjusted variant (`gam_hotspot_params_raw`)** is carried alongside so
that every domain-level result can be reported both ways. It is
threshold-and-nothing-else: no smoothing, strict adjacency, no size floor. It is
the conservative reading — it makes no assumption that a sub-threshold bin
inside a run belongs to the domain — at the cost of fragmenting real domains
into singletons (at q90, median run length 1 bin, 78% of hotspot bins in runs
under 100 kb). Keeping both is the same instinct as reporting enrichments
per-region and size-weighted on the ChIP side: where the two agree the
conclusion does not rest on the smoothing choice, and where they disagree that
disagreement is itself the finding. The truly unadjusted view is the per-bin
count itself, which the pipeline always retains as `n_contacts_lost` /
`n_contacts_gained`; this variant is about domain construction only.

## Large contact-direction blocks

A second, much coarser scale than `gam_hotspot_params`. The hotspot domains are
0.1–2 Mb; the tracks show one-directional territory at 10–20 Mb — e.g. chr11 is
almost entirely contact-losing over 0–20 Mb and almost entirely contact-gaining
over 56–108 Mb. A caller with a 100 kb floor and a q90 cutoff cannot see that:
it fragments a 16 Mb block into dozens of small domains and reports them as
unrelated.

These blocks are called on the signed net direction,
`net = log2((gained + 1) / (lost + 1))`, rather than on either density alone,
because the feature is the imbalance: a block can be busy in both directions and
still be unremarkable. net > 0 means contacts are preferentially gained in the
cKO, net < 0 lost.

This was `log2(lost/gained)` until 2026-10-18 — the collaborators' WT-minus-KO
sign, and the only signed quantity in the repo pointing the opposite way. It
made the genome-wide block figure shade lost blocks red while every contact map
uses red for gained. `mean_net` in `gam6_large_blocks.tsv` changes sign from
that date; `dirn` does not, being a label rather than a number.

Verified not to be a coverage artifact before building on it: the correlation
between net and the WT-minus-KO count of valid NPMI cells is +0.005 at 4 Mb and
+0.025 per bin, i.e. nothing. Blocks are not where data is missing.
