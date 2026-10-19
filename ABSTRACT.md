# Abstract

**HP1γ loss broadly depletes H4K20me3 and increases short-range heterochromatin
self-association without detectable change in A/B compartmentalisation.**

*Draft of 2026-10-19. Every number here is quoted from a table in
`ChIP/results/tables/` or `GAM/results/tables/`, which are the source of truth;
`cd GAM && Rscript verify_GAM.R` checks that this file still agrees with them.
Re-run it after any stage re-run and fix this file, not the table.*

---

We profiled HP1γ (Cbx3) conditional-knockout cortex (*Emx1-Cre*) against wild
type with MINUTE multiplexed ChIP-seq for five histone marks and Genome
Architecture Mapping at 20 kb resolution (mm39; 2 vs 2 biological replicates,
6 vs 4 libraries after excluding HP1gKO rep2 at QC).

**Mark loss is broad and mark-selective.** H4K20me3 falls almost genome-wide
(median log2FC −0.497 on the 20 kb grid, −0.68 at peak level, ~95% of peaks
reduced), ahead of H3K9me3 (−0.290) and H3K9me2 (−0.179), while H3K36me3
(−0.052) and H3K4me3 (0.000) are unchanged. H4K20me3 loss is largely uncoupled
from H3K9me3: regions losing H4K20me3 alone outnumber co-loss regions under both
cutoff schemes.

**Contacts redistribute from long to short range, graded by chromatin class.**
Below 0.2 Mb, heterochromatin–heterochromatin NPMI rises +4.18% (95% CI
3.36–4.92) against +2.04% for heterochromatin–euchromatin and +0.95% for
euchromatin–euchromatin; above 10 Mb the ordering inverts (−0.60%, −6.38%,
−8.91%). The heterochromatin-versus-euchromatin contrast below 1 Mb holds on
**17 of 19 autosomes** in both sub-bands (chr1 and chr3 are the exceptions) and
is invariant to the large per-chromosome NPMI offsets that dominate the raw
values. The monotone gradient across the three classes argues against
heterochromatin–euchromatin mixing, which would compress the classes together
rather than separate them.

**Remodelling is domain-specific.** Domains losing contacts are enriched for
constitutive LADs (1.35×), LINE-rich (1.24×) and quiescent chromatin (1.10×) and
depleted of SINEs (0.80×), TSSs (0.91×) and CTCF peaks (0.88×); all hold against
both a shuffled and a rotated null. Domains gaining contacts are
indistinguishable from background. GAM-derived insulation boundaries are
CTCF-enriched (1.72×), validating the boundary calls independently of the
remodelling result.

**Two powered negatives.** Per-bin mark loss does not predict per-bin contact
change (|ρ| ≤ 0.021 for every mark × direction, n ≈ 120,000), although a ~250×
mismatch between the scale of the mark domains and that of the contact bins
limits what this can exclude. Compartment strength is unchanged on two axes
external to both contact matrices (active-chromatin axis mean Δ −0.00018, 95% CI
−0.00213 to +0.00186, p = 0.62; TSS-density axis p = 0.74) — but see the
provenance note below before relying on those two numbers.

**Conclusion.** The architectural consequence of HP1γ loss is a short-range gain
in heterochromatin self-association over LAD/LINE-rich quiescent domains, not
compartment-scale reorganisation.

---

## What this abstract deliberately does not claim

**No per-genotype A/B comparison.** The A/B axis is weakly determined in this
dataset: against an independent active-chromatin ruler the best eigenvector
reaches median |ρ| 0.39, and the axis is determined on only **16 of 19**
autosomes — chr8, chr14 and chr15 have no usable axis at all. Removing per-bin detectability before the decomposition fixed
*which* eigenvector carries the compartment (PC1 on 0 of 18 → 8 of 16) but did
not improve agreement with active chromatin (0.42 → 0.39; 8 chromosomes up, 7
down), so the ceiling is the data rather than a remaining parameter choice.
Worse, a cKO/WT strength ratio computed along an eigen-derived axis takes the
value 0.540, 0.786 or 1.207 depending only on whether the axis came from WT, the
average, or the cKO — the axis choice sets the answer. An earlier claim that
compartments are less defined in the cKO was retracted for exactly this reason.

**What survives is the negative result**, because the two external-axis tests
project both genotypes onto ChIP-derived and TSS-derived rulers and use no
eigenvector at all. Both are null across all 19 autosomes with a tight CI.
Compartment figures elsewhere in the repo are kept for the structure they show,
not for a differential claim. Full record: `GAM/results/COMPARTMENTS.md`.

**The compartment-strength null is not currently regenerable, and that is a
problem.** It was computed on 2026-10-17 and the result was written into prose
(`GAM/CLAUDE.md`, `GAM/make_report.R`) without being persisted as a stage or a
table. No file in `results/tables/` contains it, no committed script produces
it, and `verify_GAM.R` therefore cannot check it — so the single compartment
claim this abstract does make is the one piece of the analysis that breaks the
repo's own rule that everything in `results/` regenerates from committed code.
The numbers are quoted here as they stand in the commit that settled the
question (`f5096fa`), not re-derived. Until a stage emits them, treat them as
provisional.

**Open provenance questions** are listed in `GAM/results/METHODS.md`; the one
that bears on this abstract is Q9, an external TAD annotation showing no CTCF
enrichment at any window size (1.09×, flat) while the same CTCF track reaches
5.27× at TSSs and 1.72× at our own boundaries. That annotation is therefore not
used for any claim above.
