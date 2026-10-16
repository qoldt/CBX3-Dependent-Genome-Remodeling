# ================================================================
# GAM pipeline runner - executes all stages in sequence.
# Run from the GAM/ directory (the folder containing config.R):
#   cd GAM && Rscript run_GAM.R
# or, in an interactive session with the GAM/ working directory:
#   source("run_GAM.R")
#
# Environment variables (all optional):
#   GAM_DATA=/path          GAM inputs (NPMI_matrices/, Diff_contacts/, Hotspots/);
#                           default "." (the data ships inside GAM/, git-ignored)
#   GAM_RESULTS=dir         write the whole output tree elsewhere; default "results"
#   CHIP_RESULTS=dir        ChIP results tree holding the annotated peak .rds;
#                           default ../ChIP/results
#   MINUTE_DATA=/path       ChIP bigWigs, read by GAM_1 via bw_loci(); default ../ChIP/data
#   MINUTE_AUTO_INSTALL=0   do not auto-install missing packages (error instead)
#
# Prerequisite: the ChIP pipeline must have run first (cd ../ChIP && Rscript
# run_MINUTE.R) - GAM_1 reads its annotated_results_*.rds and the bigWigs.
#
# Dependency chain (NOT a flat list):
#   GAM_1  ->  GAM_2  ->  { GAM_3, GAM_4, GAM_5 }
#   GAM_1 writes results/rds/gam_bins_annotated.rds
#   GAM_2 writes results/rds/npmi_bin_features.rds (needs GAM_1: PC1 is oriented
#         using the WT H3K4me3 signal from GAM_1's rds)
#   GAM_3/4/5 read those two .rds and are independent of each other, so any of
#   them can be re-run alone once GAM_1 and GAM_2 have produced their output.
#
# Cost: GAM_1 (bw_loci over every bigWig on the 20 kb grid) and GAM_2 (parses
# the dense NPMI matrices, one chromosome at a time) are sustained all-core
# local jobs. Run them under the repo's `safe-local-alignment` guidance
# (thermal/swap protection on this Apple Silicon Mac). GAM_3-5 need only the two
# .rds files from GAM_1/GAM_2 and are comparatively light.
# ================================================================
source("config.R")

message("=== GAM_1: bin grid + ChIP signal/stats + annotation (heavy) ===")
source("GAM_1_Bins_and_Signal.R")

message("=== GAM_2: NPMI matrices -> PC1 compartments, insulation, decay (heavy) ===")
source("GAM_2_NPMI_Compartments.R")

message("=== GAM_3: differential-contact coupling to ChIP changes (needs GAM_1, GAM_2) ===")
source("GAM_3_Contact_Coupling.R")

message("=== GAM_4: GAM compartments vs ChIP compartments vs contacts (needs GAM_1, GAM_2) ===")
source("GAM_4_Compartments_vs_Contacts.R")

message("=== GAM_5: HUSH targets and repeats vs contact remodelling (needs GAM_1, GAM_2) ===")
source("GAM_5_HUSH_and_Repeats_Contacts.R")

message("=== GAM pipeline complete ===")
