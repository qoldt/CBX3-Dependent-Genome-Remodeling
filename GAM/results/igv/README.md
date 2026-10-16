# Minimal IGV view - GAM/ChIP integration (mm39)

Open `gam_chip_minimal.xml` in IGV (File > Open Session), or drag the
individual files onto an mm39 genome. One track per data type:

| Track | What it is |
|---|---|
| H4K20me3 log2FC | the ChIP effect, per 20 kb bin, KO vs WT |
| GAM contacts lost | differential contacts stronger in WT (lost in the KO) |
| Hotspot domains | called remodelled domains, adjusted variant |
| GAM PC1 (WT) | compartment axis; gaps = no usable axis on that chromosome |
| ChIP compartments | stable / H4K20me3_only / co_loss |
| GAM large blocks | 10-25 Mb one-directional territory (red lost, blue gained) |

Default locus is the whole of chr11, which shows the dominant feature most
clearly: contacts are lost across essentially all of 0-25 Mb and gained
across 55-110 Mb. That large-block structure covers 41.6% of the genome.

The PC1 track is NOT reliable genome-wide - only 11 of 38 chromosome x
genotype axes reach |rho| > 0.4. See `../../CLAUDE.md` before using it.

Read together: H4K20me3 falls almost everywhere, while contact direction is
organised into very large blocks - the two are uncoupled at bin level
(GAM_3: rho -0.015). Zoom out to see the blocks; zoom in and they dissolve.
Regenerate with `cd GAM && Rscript export_igv.R`.
