# ================================================================
# MINUTE pipeline runner — executes all stages in sequence.
# Run from the project root (the folder containing config.R):
#   Rscript run_MINUTE.R
# or, in an interactive session:
#   source("run_MINUTE.R")
# ================================================================
source("config.R")

message("=== MINUTE_1: quantify + DESeq2 + annotate ===")
source("MINUTE_1_Count_and_Annotate.R")

message("=== MINUTE_2: hypergeometric enrichment ===")
source("MINUTE_2_Hypergeometric_test.R")

message("=== MINUTE_3: heatmap + significant-region BED export ===")
source("MINUTE_3_heatmap.R")

message("=== MINUTE pipeline complete ===")
