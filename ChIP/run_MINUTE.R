# ================================================================
# MINUTE pipeline runner — executes all stages in sequence.
# Run from the ChIP/ directory (the folder containing config.R):
#   cd ChIP && Rscript run_MINUTE.R
# or, in an interactive session with the ChIP/ working directory:
#   source("run_MINUTE.R")
# bigWigs/annotation default to ChIP/data/bigwig and ChIP/data/annotation
# (download per the README). To use a copy staged elsewhere, set MINUTE_DATA:
#   MINUTE_DATA=/path/to/store Rscript run_MINUTE.R
# ================================================================
source("config.R")

message("=== MINUTE_1: quantify + DESeq2 + annotate ===")
source("MINUTE_1_Count_and_Annotate.R")

message("=== MINUTE_2: hypergeometric enrichment ===")
source("MINUTE_2_Hypergeometric_test.R")

message("=== MINUTE_3: heatmap + clusters + change/relationship plots ===")
source("MINUTE_3_heatmap.R")

message("=== MINUTE_4: k-means cluster characterisation ===")
source("MINUTE_4_cluster_analysis.R")

message("=== MINUTE pipeline complete ===")
