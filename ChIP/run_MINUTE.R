# ================================================================
# MINUTE pipeline runner — executes all stages in sequence.
# Run from the ChIP/ directory (the folder containing config.R):
#   cd ChIP && Rscript run_MINUTE.R
# or, in an interactive session with the ChIP/ working directory:
#   source("run_MINUTE.R")
# bigWigs/annotation default to ChIP/data/bigwig and ChIP/data/annotation
# (download per the README). To use a copy staged elsewhere, set MINUTE_DATA:
#   MINUTE_DATA=/path/to/store Rscript run_MINUTE.R
#
# Stages are organised thematically; all read the annotated .rds from stage 1.
# Stage 6 additionally reads the family-exon signal persisted by stage 5.
# ================================================================
source("config.R")

message("=== MINUTE_1: quantify + DESeq2 + annotate ===")
source("MINUTE_1_Count_and_Annotate.R")

message("=== MINUTE_2: global changes (heatmap, clusters, relationships, TAD/ChromHMM) ===")
source("MINUTE_2_Global_Changes.R")

message("=== MINUTE_3: differential H3K9me3 vs H4K20me3 loss ===")
source("MINUTE_3_Differential_loss.R")

message("=== MINUTE_4: repeats (enrichment + direct signal change) ===")
source("MINUTE_4_Repeats.R")

message("=== MINUTE_5: clustered HUSH/CBX3 gene families ===")
source("MINUTE_5_Clustered_Gene_Families.R")

message("=== MINUTE_6: intersection with silencing machinery (KAP1/TRIM28) ===")
source("MINUTE_6_Intersection.R")

message("=== MINUTE pipeline complete ===")
