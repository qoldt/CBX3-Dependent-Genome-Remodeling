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
#   GAM_1  ->  GAM_2  ->  { GAM_2b, GAM_2c, GAM_3, GAM_4, GAM_5, GAM_6, GAM_7 }
#   GAM_1 writes results/rds/gam_bins_annotated.rds
#   GAM_2 writes results/rds/npmi_bin_features.rds (needs GAM_1: PC1 is oriented
#         using the WT H3K4me3 signal from GAM_1's rds)
#   GAM_2b writes results/rds/shared_axis_scores.rds and needs only GAM_1 plus
#         the NPMI matrices; GAM_6 joins its score_delta, so it runs before GAM_6.
#   GAM_2c writes results/rds/meltron_is.rds + meltron_scores.rds from the NPMI
#         matrices alone, and reads GAM_1's hotspot domain tables if present.
#   GAM_3-7 read those .rds and are independent of each other, so any of
#   them can be re-run alone once GAM_1 and GAM_2 have produced their output.
#   (GAM_6 section C and GAM_7 section B also read tables GAM_4 writes.)
#
# Cost: GAM_1 (bw_loci over every bigWig on the 20 kb grid) and GAM_2 (parses
# the dense NPMI matrices, one chromosome at a time) are sustained all-core
# local jobs. Run them under the repo's `safe-local-alignment` guidance
# (thermal/swap protection on this Apple Silicon Mac). GAM_3-7 need the two .rds
# files from GAM_1/GAM_2; GAM_5 and GAM_7 also re-read the NPMI matrices, so they
# are not light either - run the whole thing under the guard, or use
# run_overnight.sh, which does that and checks each stage's outputs.
# ================================================================
source("config.R")

message("=== GAM_1: bin grid + ChIP signal/stats + annotation (heavy) ===")
source("GAM_1_Bins_and_Signal.R")

message("=== GAM_2: NPMI matrices -> PC1 compartments, insulation, decay (heavy) ===")
source("GAM_2_NPMI_Compartments.R")

message("=== GAM_2b: one shared A/B axis, both genotypes scored on it ===")
source("GAM_2b_Shared_Axis.R")

message("=== GAM_2c: MELTRON - chromatin melting and condensation ===")
source("GAM_2c_Meltron.R")

message("=== GAM_4: large one-directional blocks and the block contact budget ===")
source("GAM_4_Large_Scale.R")

message("=== GAM_5: heterochromatin domains and their contacts ===")
source("GAM_5_Heterochromatin.R")

# GAM_3 AFTER GAM_5, not before. Its synthesis panel reads
# gam10_class_change.tsv, which GAM_5 writes - in the old order that read always
# picked up the PREVIOUS run's table, so a full pipeline run produced a figure
# one generation out of date and nothing said so.
message("=== GAM_3: does contact remodelling coincide with the mark loss? ===")
source("GAM_3_Coupling.R")

message("=== GAM_6: the derepressed targets (families, repeats, protocadherin) ===")
source("GAM_6_Targets.R")

# GAM_8 after GAM_4 (it needs gam6_large_blocks.tsv) and after GAM_1. It is
# placed here rather than next to GAM_4 only so the heavy NPMI passes do not
# run back to back.
message("=== GAM_8: chromatin classes of the contacts inside the large blocks ===")
source("GAM_8_Block_Classes.R")

message("=== GAM_7: region plots and per-chromosome overviews ===")
source("GAM_7_Figures.R")

message("=== GAM pipeline complete ===")

# Consistency check. Runs in a SEPARATE process on purpose: it must read the
# tables as they now sit on disk, not as this session happens to hold them in
# memory, and a non-zero exit here must not kill a pipeline that has already
# produced its outputs. It is advisory at this point - the failure it reports is
# always "a file quotes a number the tables no longer say", which is for a human
# to fix, not for the pipeline to retry.
message("=== verify: are the outputs and the prose that cites them consistent? ===")
vrc <- system2("Rscript", "verify_GAM.R")
if (!identical(vrc, 0L))
  message("\n*** verify_GAM.R reported inconsistencies - see above. ***\n",
          "The tables are right by definition; fix whatever quotes them.\n")

message("Report:  Rscript make_report.R      ->  results/GAM_REPORT.html")
message("Verify:  Rscript verify_GAM.R       ->  exit 1 if anything disagrees")
