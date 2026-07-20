# Why the H3K9me3 / H4K20me3 compartment cutoffs must be read median-relative

**Short version:** MINUTE_3 sorts peaks into `stable` / `H4K20me3_only` /
`co_loss` using fixed log2FC cutoffs. Because *both* marks lose genome-wide,
those cutoffs largely ask *"did it lose?"* — to which most of the genome answers
yes — instead of *"did it lose more than typical?"*. That makes the compartment
sizes depend on where the global shift happens to sit, which in turn depends on
which replicates are included. Centring each mark on its own median fixes it.
**The uncoupling result survives the fix; only its phrasing has to change.**

Figure: `figures/differential_loss/diffloss_cutoff_schemes.png`
Table: `tables/diffloss_group_definitions.tsv` (both schemes, regenerated per run)

---

## 1. What the classification is trying to do

MINUTE_3 asks whether the two heterochromatin marks are lost **together** or
**independently**, by binning each shared peak with fixed thresholds:

| rule | cutoff |
|---|---|
| H4K20me3 "lost" | log2FC < **−0.5** |
| H3K9me3 "lost" | log2FC < **−0.3** |
| either "unchanged" | \|log2FC\| < **0.15** |

- `co_loss` = H4K20me3 lost **and** H3K9me3 lost
- `H4K20me3_only` = H4K20me3 lost **and** H3K9me3 unchanged
- `stable` = both unchanged

## 2. Why that breaks here

Those numbers implicitly assume the average peak doesn't change. In this dataset
it does: **the typical peak has H3K9me3 log2FC = −0.415** (HP1gKO rep2 excluded),
which is already past the −0.3 "lost" line. A completely unremarkable,
bang-average peak is therefore auto-labelled "H3K9me3 lost".

> **63.1% of all peaks pass the "H3K9me3 lost" test.**
> A criterion satisfied by two-thirds of the genome is not identifying anything.

The mirror image does the other damage: "H3K9me3 unchanged" (\|log2FC\| < 0.15)
now requires a peak to sit ~0.4 **above** the genome median — genuinely unusual.
So `H4K20me3_only`, which *requires* H3K9me3 to be unchanged, is starved by
construction while `co_loss` absorbs almost everything.

### A single peak makes it concrete

Take a peak with H4K20me3 = **−1.30**, H3K9me3 = **−0.42**:

| | test | verdict |
|---|---|---|
| absolute | −1.30 < −0.5 ✓ and −0.42 < −0.3 ✓ | **co_loss** — "both marks lost" |
| median-centred | −1.30 − (−0.678) = **−0.62** ✓ lost; −0.42 − (−0.415) = **−0.005** → typical | **H4K20me3_only** |

Its H3K9me3 is *identical to the genome median*. It lost exactly the average
amount — there is no evidence H3K9me3 was specially affected there. The absolute
rule calls it co-loss purely because the genome-wide median already sits past the
cutoff.

This is not a rare edge case:

> **64.7% of the 46,921 `co_loss` peaks have H3K9me3 within ±0.15 of the genome
> median.** Only **17.1%** genuinely lose more H3K9me3 than typical.

The compartment labelled "both marks lost together" is nearly two-thirds peaks
where only one mark did anything distinctive.

## 3. The fix

Subtract each mark's own median before applying the same cutoffs. The question
changes from *"did it lose?"* to *"did it lose **more than typical for this
mark**?"* — which is the comparative question an uncoupling claim actually needs,
and which does not move when the global shift moves.

Selectivity returns to something meaningful: **12.4%** of peaks clear the
median-centred H3K9me3 cutoff, versus **63.1%** for the absolute one.

## 4. The evidence that the centred version is better

A classification that reflects biology should give the same answer regardless of
which replicates happen to be included. Both schemes were run on two sample sets
differing by a single library (all libraries vs HP1gKO rep2 excluded):

| scheme | stable | H4K20me3_only | co_loss | **mean change across sample sets** |
|---|---|---|---|---|
| absolute | 4,105 → 931 | 13,682 → 8,481 | 19,938 → 46,921 | **84%** |
| median-centred | 19,418 → 17,305 | 1,798 → 1,748 | 1,452 → 1,171 | **11%** |

Dropping one library restructures the absolute compartments completely — no
biology changed, only the global shift moved and the fixed cutoffs landed
somewhere different relative to it. The centred groups barely move (~8× more
stable).

## 5. What this means for the result

Under **absolute** cutoffs with rep2 excluded, `co_loss` (46,921) swamps
`H4K20me3_only` (8,481) — which reads as *"the marks lose together"*, the
opposite of the published conclusion.

Under **median-centred** cutoffs, `H4K20me3_only` exceeds `co_loss` in **both**
sample sets (1,798 vs 1,452 with all libraries; 1,748 vs 1,171 without rep2).

> **The uncoupling result holds.** But it must be stated as *"loses more
> H4K20me3 than its own genome-wide average, while losing a typical amount of
> H3K9me3"* — not *"loses H4K20me3 in absolute terms while H3K9me3 is
> unchanged"*. The absolute framing made the claim look as though it had
> collapsed when rep2 was removed, and would have collapsed it for real had
> anyone re-run with a different replicate set.

## 6. Status in the pipeline

The **absolute** scheme is still what MINUTE_3's downstream outputs use, so the
existing tables and figures are unchanged. The median-centred counts are reported
alongside them per run in `tables/diffloss_group_definitions.tsv`, so the
dependence is visible in the outputs rather than only in prose.

Switching the downstream definition would move every table and figure in that
stage — a larger decision, but the median-centred version is the more defensible
one to publish.
