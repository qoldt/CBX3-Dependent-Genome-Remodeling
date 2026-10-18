# Per-bin balancing of the GAM difference map

Hand-written. Covers the `_allcells_binbalanced` region-plot variant: what it
does, why, what it cannot distinguish, and what it does **not** touch.

**It is a display correction. It changes no statistic, no table and no
significance call.** See "Scope" at the bottom before anything else.

---

## 1. The problem

GAM detects whether a 20 kb window is present in a nuclear profile. Windows
differ in how easily they are detected — mappability, repeat content, how the
locus survives cryosectioning — and that detectability is not guaranteed to be
identical between genotypes, which are separate libraries with separate profile
counts.

A per-window effect is not a per-pair effect. If window *i* is detected slightly
less often in the cKO, **every pair involving *i*** shifts: its whole row and its
whole column of the matrix move together.

In the rotated triangle the row and the column of a single bin are two diagonal
lines meeting at that bin, so one biased window draws a **V**. Several hundred
of them produce the criss-cross plaid lattice visible in the `_allcells` maps.

### Why NPMI does not already fix it

It partly does, and this is the obvious objection. NPMI is built on pointwise
mutual information,

    PMI(i,j) = log[ p(i,j) / (p(i) p(j)) ]

so each window's **marginal detection probability** `p(i)` is divided out by
construction, then the result is normalised to [-1, 1]. A first-order per-window
detectability difference is therefore already removed.

What survives:

- the normalisation is exact only under the independence model it is built on;
- **variance** is not normalised — a poorly detected window has a noisier NPMI
  row, and the per-distance z-scoring applied downstream rescales that noise
  differently in each genotype;
- the two genotypes are z-scored **independently**, so any residual per-window
  deviation that differs in magnitude between them reappears in the difference.

So what is measured below is a *residual* per-bin term after NPMI's own
normalisation, not the raw detectability bias.

### Why distance z-scoring does not fix it either

`zbydist()` z-scores within each genomic distance. That removes a function of
`|i - j|`. A per-bin effect is a function of `i` and of `j` separately. The two
are orthogonal apart from a constant, so one cannot remove the other.

---

## 2. The model

Write the displayed difference map as

    D[i,j] = a_i + a_j + e[i,j]

- `a` is one number per bin — everything about bin *i* that is true of *i*
  against **every** partner.
- `e[i,j]` is what depends on the specific pair: the pairwise signal plus noise.

`D` is symmetric and the diagonal is excluded throughout.

### Identifiability

`a` is determined only up to an additive constant: replacing `a_i` with
`a_i + c` adds `2c` to every cell. The model therefore needs a constraint, and
we use

    mean(a) = 0

which fixes the constant and leaves the overall level of the map alone.

### Estimator

With the constraint above, the row means of the partial residual recover `a`
directly:

    mean_j( D[i,j] - a_j )  =  a_i + mean(a) + mean_j(e[i,j])
                            =  a_i                              (+ noise)

so the estimate is a fixed-point iteration, damped by averaging with the
previous iterate for stability, NA-aware so unmeasured cells contribute nothing:

```r
A <- D; diag(A) <- NA
aa <- rep(0, nrow(A))
for (it in 1:60)
  aa <- (aa + rowMeans(sweep(A, 2, aa, "-"), na.rm = TRUE)) / 2
aa[!is.finite(aa)] <- 0
D <- D - outer(aa, aa, "+")
```

60 iterations is well past convergence at these sizes; the damping factor of 1/2
is there because the undamped iteration oscillates.

### How much of the map it is

Variance explained, `1 - var(D - fitted) / var(D)`, over the pairs actually
drawn:

| region | width | per-bin term |
|---|---|---|
| chr3:30-40 Mb (plain control window) | 10 Mb | **27.9%** |
| chr18 Pcdh | 1.8 Mb | 26.5% |
| chr2 Olfr | 7.6 Mb | 22.8% |
| chr11:1-25 Mb lost block | 24 Mb | 18.0% |

The control window is the **highest**. This is a property of the data
everywhere, not of unusual loci.

---

## 3. Relation to published practice

Per-bin bias correction is established for GAM. `normGAM`
(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6936146/) provides five
normalisations for GAM matrices, two of which are matrix balancing: **ICE**
(iterative correction and eigenvector decomposition) and **Knight-Ruiz
2-norm**, alongside normalised linkage disequilibrium, vanilla coverage and
sequential component normalisation. In Hi-C the same idea is ICE
(Imakaev et al. 2012) and KR (as used in Juicer).

Two differences from that prior art, stated so this is not mistaken for an
off-the-shelf method:

1. **Those balance a single raw matrix multiplicatively**, `M[i,j] / (b_i b_j)`,
   forcing equal row sums. This balances a **difference of two
   distance-z-scored matrices**, which is a different object.
2. **Additive, not multiplicative.** A multiplicative correction on counts is an
   additive one in log space; `D` here is already a difference of z-scores, so
   the additive form is the matching analogue. Dividing z-scores by a per-bin
   factor would not be meaningful.

---

## 4. What it cannot distinguish

**A genuine bin-level biological effect and a technical detectability
difference are the same thing to this model.** Both are "bin *i* is off by a
constant against every partner".

A window that truly loses contact with everything — which is what a LAD bin
detaching from the lamina would look like — is removed by this correction along
with the artefact. That is why it is a **separate variant and not applied to the
other maps**.

Read the pair:

- a feature present in `_allcells` **and** `_allcells_binbalanced` is pairwise
  structure;
- a feature that disappears was either a detection artefact **or** a real
  whole-bin effect, and this map cannot tell you which. Use the per-bin tracks
  (differential-contact counts, mark coverage, LAD/repeat annotation) to decide.

---

## 5. Scope

Balancing is applied in **two** places, for two different reasons. When this
document was first written it was only the first, and the list below said no
statistic used it. That stopped being true on 2026-10-18 and the claim is
corrected here.

**(a) Figures — `gam_region_plot()`**, to the display matrix `D` only, after
aggregation and before the tile data is built. Purely cosmetic: it removes the
criss-cross plaid lattice so the region maps are readable.

**(b) Compartments — before the PCA, in `GAM_2_NPMI_Compartments.R`
(`gam2_balance_before_pca`, params v8) and `GAM_2b_Shared_Axis.R`
(`gam2b_balance_before_pca`, params v2).** This one is not cosmetic and it does
change results. The same per-bin term that draws the plaid also loads onto the
leading eigenvector of the O-E correlation matrix, which pushes the compartment
down the spectrum and makes the A/B axis a mixture of compartment and
detectability. Removing it first takes the compartment back to PC1 on much of
the genome: in `GAM_2`, PC1 went from 2 of 24 to 16 of 27; in `GAM_2b`, from
0 of 18 to 8 of 16, with 15 of 16 chromosomes moving up by exactly one
eigenvector — the signature of one eigenvector having been spent on
detectability.

**Read the trade honestly.** Balancing does not make the axis agree better with
active chromatin. On the chromosomes determined both with and without it, median
|rho| moves 0.423 -> 0.387 and the sign test is 8 up / 7 down: a wash. What it
buys is *which* eigenvector carries the compartment and a more even A/B split
(median |frac_A - 0.5| 0.061 -> 0.043, worst case 0.198 -> 0.117). The ceiling
is the data, as `COMPARTMENTS.md` section 4b says. There is also a cost in
coverage, recorded there.

**A specific caveat for GAM.** In Hi-C, balancing targets a near-pure technical
bias. In GAM, detection frequency is *itself partly biological* — open chromatin
is more detectable — so `a[i] + a[j]` is part artefact and part compartment, and
subtracting it removes some real signal along with the bias. `GAM/CLAUDE.md`
settled point 7 records a leading eigenvector tracking per-bin `n_valid` at
|rho| 0.34-0.47 on chr7 and chr11. This is the most likely reason balancing
helps the eigenvector ordering while leaving |rho| flat.

**Still not touched:**

- the delivered `Diff_contacts` files and their significance criterion
  (`|delta-z| > 2` and `max NPMI > 0.1`), which are the collaborators' and are
  used exactly as delivered;
- the per-bin differential-contact counts and the hotspot bins and domains
  derived from them;
- every statistic in `results/tables/` **except** the compartment outputs named
  in (b) above and anything downstream of them (`npmi_compartment_*`,
  `gam2b_shared_axis_summary`, the compartment BEDs, `compartment_switch` and
  its consumers in `GAM_3` and `GAM_6`);
- the `_sigonly` maps, and the slide figures, which use `_sigonly`.

Only the `_allcells_binbalanced` variant sets `balance = TRUE`, and that variant
is unmasked by construction. Balancing a significance-masked map would be
meaningless — it would fit row effects on a sparse, selected subset — and the
function guards against the combination.
