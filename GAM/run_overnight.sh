#!/usr/bin/env bash
# ==============================================================================
# Unattended overnight driver for the GAM analysis batch.
#
#     cd GAM && ./run_overnight.sh [logfile]
#
# WHY A DRIVER AND NOT AGENTS. Three of four subagents in the session that
# produced this batch hung after emitting one sentence and returned nothing, and
# unsourced numbers from one of them reached a summary and had to be retracted.
# A script leaves artifacts on disk whether or not anything is watching, and is
# re-runnable. Everything here is therefore a committed R stage.
#
# WHY SEQUENTIAL. Several stages parse the NPMI matrices (7.7 GB, chr1 alone is
# ~760 MB per genotype as doubles) and this machine is thermally limited - see
# the safe-local-alignment skill, which exists because a parallel local job once
# forced a thermal emergency shutdown. Everything runs under thermal_guard.sh,
# one stage at a time.
#
# FAILURE POLICY. A stage that fails does NOT abort the night: the remaining
# stages still run, and the run ends with a PASS/FAIL table. Each stage declares
# the outputs it must produce, and they are checked for existence AND non-zero
# size - "exit 0" is not evidence a stage did anything, which is the same
# mistake the skill's verify_star_logs.sh exists to catch.
# ==============================================================================
set -uo pipefail
cd "$(dirname "$0")"

LOG="${1:-$PWD/overnight.log}"
GUARD="$HOME/.claude/skills/safe-local-alignment/scripts/thermal_guard.sh"
[ -x "$GUARD" ] || { echo "thermal guard not found at $GUARD" >&2; exit 1; }

RESULTS="${GAM_RESULTS:-results}"
declare -a NAMES STATUS

say() { printf '%s %s\n' "$(date '+%H:%M:%S')" "$*" | tee -a "$LOG"; }

# run <label> <script> <required output> [<required output> ...]
run() {
  local label="$1" script="$2"; shift 2
  local t0 rc=0 missing=""
  say "=== START $label ($script)"
  t0=$(date +%s)
  "$GUARD" -l "$LOG.$label" -p "R" -- Rscript "$script" >>"$LOG" 2>&1 || rc=$?
  for out in "$@"; do
    if [ ! -s "$out" ]; then missing="$missing $out"; fi
  done
  local dt=$(( $(date +%s) - t0 ))
  if [ $rc -ne 0 ]; then
    say "=== FAIL  $label (exit $rc, ${dt}s)"; NAMES+=("$label"); STATUS+=("FAIL exit=$rc")
  elif [ -n "$missing" ]; then
    say "=== EMPTY $label (${dt}s) missing:$missing"; NAMES+=("$label"); STATUS+=("EMPTY$missing")
  else
    say "=== PASS  $label (${dt}s)"; NAMES+=("$label"); STATUS+=("PASS ${dt}s")
  fi
}

say "overnight batch starting; log=$LOG"

# Priority 1 - heterochromatin domains and the Pcdh TAD question
run het_domains GAM_13_Het_Domains_Expanded.R \
  "$RESULTS/tables/gam13_merge_sweep.tsv" \
  "$RESULTS/tables/gam13_domain_definitions.tsv" \
  "$RESULTS/tables/gam13_intra_vs_inter.tsv"

run pcdh_tads GAM_14_Pcdh_TAD_Interactions.R \
  "$RESULTS/tables/gam14_pcdh_tad_summary.tsv" \
  "$RESULTS/tables/gam14_pcdh_by_separation.tsv"

# Priority 2 - per-chromosome characterisation
run chromosomes GAM_15_Chromosome_Reports.R \
  "$RESULTS/tables/gam15_chromosome_ranking.tsv" \
  "$RESULTS/tables/gam15_contact_budget.tsv" \
  "$RESULTS/chromosomes/chr11.md"

# Priority 3 - region plots, now including the other clustered families
run region_plots GAM_12_Region_Plots.R \
  "$RESULTS/figures/regions/chr18_Pcdh.pdf" \
  "$RESULTS/figures/regions/chr7_Vmn.pdf" \
  "$RESULTS/figures/regions/chr2_Olfr.pdf"

# The report goes last so it picks up everything that succeeded.
run report make_report.R "$RESULTS/GAM_REPORT.html"

say ""
say "================ overnight summary ================"
for i in "${!NAMES[@]}"; do printf '%-14s %s\n' "${NAMES[$i]}" "${STATUS[$i]}" | tee -a "$LOG"; done
say "==================================================="
fails=0
for s in "${STATUS[@]}"; do case "$s" in PASS*) ;; *) fails=$((fails+1));; esac; done
say "$fails stage(s) did not pass"
exit 0
