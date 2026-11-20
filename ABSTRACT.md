# Abstract

**HP1γ loss broadly depletes H4K20me3 and increases short-range heterochromatin
self-association.**

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
Bins are partitioned on two nested levels: euchromatin (TSS, enhancer,
transcribed) against heterochromatin (Polycomb, H3K9me3, quiescent). Below
0.2 Mb, heterochromatin–heterochromatin NPMI rises **+1.60%** while
heterochromatin–euchromatin is +0.23% and euchromatin–euchromatin falls
−0.58%. All three then fall with distance, deepest at 10–50 Mb
(het–het **−9.79%**, het–eu −10.62%, eu–eu −10.93%), and recover partway in
the >50 Mb band (−5.93%, −6.83%, −7.61%) where the confidence intervals all
reach 0 and the band is not comparable across chromosomes — chr19 contributes
only 50–61 Mb against chr1's 50–195 Mb. The classes are strictly ordered
het–het > het–eu > eu–eu in every one of the six distance bands, and the
heterochromatin-versus-euchromatin contrast below 1 Mb now holds on **19 of 19
autosomes** in both sub-bands — offset-invariant, since the large
per-chromosome NPMI offsets cancel in a within-chromosome contrast. The
monotone gradient argues against heterochromatin–euchromatin mixing, which
would compress the classes together rather than separate them.

**The gain is specific to H3K9me3 heterochromatin, not to repression in
general.** At level 2, H3K9me3–H3K9me3 pairs rise **+4.13%** below 0.2 Mb,
peak at +5.91% over 0.2–1 Mb, and are the only class still above zero out to
10–50 Mb (+4.49% at 5–10 Mb, +1.44% at 10–50 Mb, the last with a CI spanning
0). Polycomb–Polycomb (+2.56%) and quiescent–quiescent (+1.43%) rise below
0.2 Mb and then fall away, and TSS–TSS collapses at long range (−24.33% at
5–10 Mb). An earlier
binary partition that called everything without an H3K9me3/H4K20me3 peak
"euchromatin" pooled Polycomb and quiescent chromatin into that class; it
measured the same +4.2% for constitutive heterochromatin but could not show
that the other repressive classes behave differently.

**Remodelling is domain-specific.** Domains losing contacts are enriched for
constitutive LADs (1.35×), LINE-rich (1.24×) and quiescent chromatin (1.10×) and
depleted of SINEs (0.80×), TSSs (0.91×) and CTCF peaks (0.88×); all hold against
both a shuffled and a rotated null. Domains gaining contacts are
indistinguishable from background. GAM-derived insulation boundaries are
CTCF-enriched (1.72×), validating the boundary calls independently of the
remodelling result.

**A powered negative.** Per-bin mark loss does not predict per-bin contact
change (|ρ| ≤ 0.021 for every mark × direction, n ≈ 120,000), although a ~250×
mismatch between the scale of the mark domains and that of the contact bins
limits what this can exclude.

**Conclusion.** HP1γ loss has **two separable architectural consequences**, not
one, and they act on different chromatin. H3K9me3 domains gain self-association
at short range and keep gaining out to 5 Mb (**+4.13%** below 0.2 Mb, **+5.91%**
at 0.2–1 Mb, **+4.59%** at 1–5 Mb). LAD- and LINE-rich **quiescent** chromatin
does the opposite beyond 0.2 Mb (quiescent–quiescent +1.43% below 0.2 Mb, −0.96%
at 0.2–1 Mb with a CI spanning 0, then **−4.34%** at 1–5 Mb), and the mixed
H3K9me3–quiescent pair sits between them (+1.93%, −0.07%, −2.53%). That is the
same asymmetry the domain result shows from the other side: the domains that
*lose* contacts are the LAD-associated (1.35×), LINE-rich (1.24×) quiescent
ones, while those that *gain* are indistinguishable from background. The large
one-directional blocks are where this is visible — chr11:2.7–23.3 Mb (LINE
1.21×, cLAD 1.16×, 74.6% quiescent, 6.6% H3K9me3) and chr16:73.1–87.1 Mb (LINE
1.29×, SINE 0.43×, cLAD 1.78×, 88.8% quiescent, 3.3% H3K9me3) both lose heavily.

An earlier version of this conclusion read "a short-range gain in heterochromatin
self-association over LAD/LINE-rich quiescent domains", which welded the two
findings into one and attributed the gain to the chromatin that is in fact
losing. Quiescent chromatin is 69% of the genome and is what LAD/LINE-rich
domains are made of; only H3K9me3, 9% of bins, gains past 0.2 Mb.

---

## What this abstract deliberately does not claim

**No A/B compartments.** Compartments were called by eigendecomposition of the
NPMI matrices, and that was removed on 2026-11-20. NPMI is already a normalised
quantity, so its eigenvectors do not recover compartments (the leading one was
per-bin detectability), and the collaborating lab confirmed it is the wrong
input. Every compartment result went with it: the per-genotype and shared axes,
the switch counts, the compartment-strength null, and the retracted claim that
compartments are less defined in the cKO. This abstract makes no compartment
claim, positive or negative.

**Open provenance questions** are listed in `GAM/results/METHODS.md`; the one
that bears on this abstract is Q9, an external TAD annotation showing no CTCF
enrichment at any window size (1.09×, flat) while the same CTCF track reaches
5.27× at TSSs and 1.72× at our own boundaries. That annotation is therefore not
used for any claim above.
