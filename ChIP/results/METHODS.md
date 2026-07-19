# Methods notes

## KAP1 / TRIM28 ChIP-seq data source and preparation

Mouse KAP1 (TRIM28) ChIP-seq peaks were obtained from **ChIP-Atlas**
(https://chip-atlas.org), which provides uniformly reprocessed public ChIP-seq
peaks (MACS2, q < 1×10⁻⁵). No primary adult-mouse-cortex KAP1 ChIP-seq with
public processed peaks was available, so two tracks were used: (i) a
**neural-lineage** track from Neuro-2a neuroblastoma cells (GEO **GSE110032**,
KAP1 ChIP-seq, replicates GSM2986446 / GSM2986447; 15,355 peaks), used as the
primary track for its neural relevance; and (ii) the **all-cell-type aggregate**
of 47 mouse Trim28 experiments (276,439 peaks; ChIP-Atlas `Oth.ALL.05.Trim28`),
used as a well-powered robustness track, since KAP1 occupancy at its
ERV/H3K9me3-heterochromatin targets is largely conserved across cell types.

Peaks were downloaded in **mm10** coordinates
(`https://dbarchive.biosciencedbc.jp/kyushu-u/mm10/assembled/Oth.ALL.05.Trim28.*.bed`)
and lifted to **mm39** with **CrossMap** using the UCSC `mm10ToMm39.over.chain`
chain (`https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz`);
99.7 % (neural) and 99.0 % (all-cell) of peaks mapped. Lifted peaks were
restricted to the primary chromosomes and converted to Ensembl-style seqnames
(1–19, X, Y) to match the input-scaled mm39 bigWigs.

Overlap of the differential-loss compartments (and the HUSH/CBX3 gene-family
exons) with KAP1 peaks was tested against a **size-matched permutation null**
(N = 2000): for each set, regions were resampled from the appropriate universe
(the measured merged-peak set for the loss compartments; all annotated mm39 gene
exons for the gene families) matched to the query's genomic-length distribution
(20 length bins), and the overlap fraction was recomputed to build the null.
Enrichment is reported as observed/expected fold and empirical p (MINUTE_6).

## Note on HUSH targets (Danac et al., Mol Cell 2024)

The HUSH-target data in Danac et al. (Mol Cell 2024; PMID 39013473; GEO
GSE268799) is **human (hg38)** CUT&RUN in K562, and its dominant HUSH targets are
primate-specific young LINE-1 (L1HS/L1PA) with no orthologous mouse locus — so
coordinate liftOver to mm39 is not valid. The mouse-native equivalent is used
instead: enrichment for **young/active mouse L1 families (L1MdA/T/Gf/F)** against
the mm39 RepeatMasker annotation (MINUTE_5).
