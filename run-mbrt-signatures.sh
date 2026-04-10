#!/bin/bash
# run-mbrt-signatures.sh — Multi-agent autonomous MBRT signature analysis
#
# VM PREREQUISITES (manual, one-time — see Step 2 setup):
#   1. Git: git config --global credential.helper store && clone repos
#   2. Claude: install, `claude login`, overlay config from jeszyman/claude
#   3. Mounts: /mnt/gcs (gcsfuse, read-only), Box data pre-staged to local
#   4. Conda: conda env create -f ~/repos/spatial-rads/config/spatial-rads-conda-env.yaml -n spatial-rads
#   5. Data staged to $SCRATCH (see Data Staging section in plan)
#   6. ~/.claude/settings.local.json with {"hooks": {}}
#
# LAUNCH:
#   cd ~/repos/spatial-rads && nohup bash run-mbrt-signatures.sh > orchestrator.log 2>&1 &
#
set -euo pipefail

# Source nvm so `claude` (nvm-installed) is on PATH under non-interactive bash.
# Without this, nohup/login shells miss nvm's shims and the claude CLI appears missing.
export NVM_DIR="$HOME/.nvm"
[ -s "$NVM_DIR/nvm.sh" ] && . "$NVM_DIR/nvm.sh"

REPO="$HOME/repos/spatial-rads"
RESULTS="$REPO/results/signatures"
PLAN="$REPO/plan-mbrt-signatures.md"
SCRATCH="/mnt/data/spatial-rads"
TIMESTAMP=$(date +%Y%m%d-%H%M)
BRANCH="autonomous/mbrt-signatures-$TIMESTAMP"

mkdir -p "$RESULTS/data" "$RESULTS/plots" "$RESULTS/tables"

# ============================================================
# Preflight checks
# ============================================================
echo "[$(date)] Preflight checks..."

[[ -f "$PLAN" ]] || { echo "FATAL: $PLAN not found"; exit 1; }
[[ -d "$SCRATCH" ]] || { echo "FATAL: $SCRATCH not found — stage data first"; exit 1; }
command -v claude >/dev/null || { echo "FATAL: claude not installed"; exit 1; }
command -v conda >/dev/null || { echo "FATAL: conda not found"; exit 1; }

# Verify staged data
[[ -f "$SCRATCH/seurat_clustered.rds" ]] || { echo "FATAL: Mutter_01 Seurat object not staged"; exit 1; }
M02_COUNT=$(find "$SCRATCH" -name "seuratObject_*_Mutter_02_*.RDS" | wc -l)
[[ "$M02_COUNT" -ge 1 ]] || { echo "FATAL: No Mutter_02 RDS files found in $SCRATCH"; exit 1; }
echo "  Found $M02_COUNT Mutter_02 RDS files"

# Verify Mutter_01 analysis outputs
[[ -f "$REPO/dev/peak_valley_analysis/data/pv_degs_bulk.tsv" ]] || { echo "FATAL: pv_degs_bulk.tsv missing"; exit 1; }
[[ -f "$REPO/dev/peak_valley_analysis/data/mbrt4h_peak_valley.tsv" ]] || { echo "FATAL: mbrt4h_peak_valley.tsv missing"; exit 1; }
[[ -f "$REPO/dev/peak_valley_analysis/data/stripe_model.rds" ]] || { echo "FATAL: stripe_model.rds missing"; exit 1; }

# Conda env
conda run -n spatial-rads R --version >/dev/null 2>&1 || { echo "FATAL: spatial-rads conda env not working"; exit 1; }

# Disk space (need >= 15 GB free)
AVAIL_GB=$(df --output=avail -BG "$SCRATCH" | tail -1 | tr -d ' G')
[[ "$AVAIL_GB" -ge 15 ]] || { echo "FATAL: Only ${AVAIL_GB}GB free on scratch, need 15GB"; exit 1; }

echo "[$(date)] Preflight OK — ${AVAIL_GB}GB free, $M02_COUNT Mutter_02 files staged"

# ============================================================
# Feature branch
# ============================================================
cd "$REPO"
git checkout -b "$BRANCH"
echo "[$(date)] Working on branch: $BRANCH" | tee -a "$RESULTS/orchestrator.log"

# ============================================================
# Watchdog (background bash — not Claude)
# ============================================================
(
  while true; do
    {
      echo "=== WATCHDOG $(date) ==="
      df -h "$SCRATCH" 2>/dev/null || echo "Scratch disk not mounted"
      echo "Scratch usage: $(du -sh "$SCRATCH" 2>/dev/null | cut -f1)"
      echo "Results files: $(find "$RESULTS" -type f | wc -l)"
      pgrep -fa "claude" 2>/dev/null | head -3 || echo "No claude process"
      # Check for R processes (pipeline running)
      pgrep -fa "Rscript" 2>/dev/null | head -3 || echo "No R process"
      echo ""
    } >> "$RESULTS/watchdog.log"
    sleep 300
  done
) &
WATCHDOG_PID=$!
trap "kill $WATCHDOG_PID 2>/dev/null" EXIT

# ============================================================
# Validation gate functions
# ============================================================
validate_executor_output() {
  local ok=0
  [[ -f "$RESULTS/SUMMARY.md" ]] || { echo "GATE FAIL: SUMMARY.md missing"; ok=1; }
  [[ -f "$RESULTS/data/data_inventory.tsv" ]] || { echo "GATE FAIL: data_inventory.tsv missing"; ok=1; }
  local n_plots=$(find "$RESULTS/plots" -name "*.png" 2>/dev/null | wc -l)
  [[ "$n_plots" -gt 5 ]] || { echo "GATE FAIL: only $n_plots plots (need >5)"; ok=1; }
  local n_data=$(find "$RESULTS/data" -name "*.tsv" 2>/dev/null | wc -l)
  [[ "$n_data" -gt 3 ]] || { echo "GATE FAIL: only $n_data data files (need >3)"; ok=1; }
  # Check SUMMARY.md has substance
  local wc=$(wc -w < "$RESULTS/SUMMARY.md" 2>/dev/null || echo 0)
  [[ "$wc" -gt 200 ]] || { echo "GATE FAIL: SUMMARY.md only $wc words"; ok=1; }
  return $ok
}

validate_reviewer_output() {
  [[ -f "$RESULTS/REVIEW.md" ]] || { echo "GATE FAIL: REVIEW.md missing"; return 1; }
  # Check for numbered findings
  grep -q "R[0-9]" "$RESULTS/REVIEW.md" || { echo "GATE FAIL: no numbered findings in REVIEW.md"; return 1; }
  return 0
}

# ============================================================
# Shared system prompt rules
# ============================================================
SYSTEM_RULES="MANDATORY RULES:
- Never ask questions. Debug and fix problems yourself.
- /mnt/gcs/ is READ-ONLY. NEVER write to /mnt/gcs/. NEVER use gsutil to write.
- /mnt/rclone/box/ is READ-ONLY. NEVER write to Box.
- Staged data is on $SCRATCH (local disk). Use these copies, not FUSE mounts.
- Pipeline intermediates go to $SCRATCH.
- Deliverables (SUMMARY.md, figures/, tables/, data/) go to $RESULTS.
- Use conda run -n spatial-rads Rscript <script> for all R execution.
- Do NOT modify any files in dev/peak_valley_analysis/.
- Check for existing R processes before launching new ones.
- Set random seeds (seed=42) for all stochastic steps.
- Commit to the current branch when done. Do not push."

# ============================================================
# Phase 1: EXECUTOR (8h max)
# ============================================================
echo "[$(date)] === PHASE 1: EXECUTOR ===" | tee -a "$RESULTS/orchestrator.log"

timeout 8h claude -p "$(cat "$PLAN")

You are the EXECUTOR agent. Follow Executor Tasks 1-7 in the plan above, in order. All setup (data staging, conda env) is already complete — do not redo it.

Key paths:
- Mutter_01 Seurat object: $SCRATCH/seurat_clustered.rds
- Mutter_01 analysis data: $REPO/dev/peak_valley_analysis/data/
- Mutter_02 RDS files: $SCRATCH/mutter02/
- Mutter_02 sample layout: $SCRATCH/mutter02/sample_layout.tsv (slide → sample_id → treatment; per Jenn Fazzari 2026-04-09)
- Gene list XLSX: $REPO/config/CosMx-Mouse-Universal-Cell-Characterization-Gene-List-(1).XLSX
- Write deliverables to: $RESULTS/

Write R scripts to $RESULTS/ and execute them. Save output data to $RESULTS/data/ and plots to $RESULTS/plots/.

FOV-TO-SAMPLE INFERENCE: Mutter_02 has NO FOV-level map. Each slide holds 3 samples (0 Gy / MBRT 2d / SBRT 20 Gy 2d). Recover the 3-sample grouping from Seurat metadata (FOV spatial coordinates or fov field), document the method in SUMMARY.md, and sanity-check that FOV counts per inferred sample are roughly balanced.

INTERPRETIVE BOUNDARY: Mutter_01 4h MBRT is the ONLY dataset with validated peak/valley ground truth (H2AX IHC). Mutter_02 has NO H2AX — stripe/peak-valley results there are exploratory pattern-detection only, never validated biology. At all non-4h timepoints (Mutter_01 later timepoints AND all of Mutter_02), you SCORE cells with UCell (primary) and AddModuleScore (secondary, seed=42) — do NOT classify them as peak/valley. Frame results as 'signature persistence' or 'consistent with persistent spatial patterning' — never claim cells were in peaks.

CIRCULARITY: Exclude ALL 18 DDR genes from peak/valley signatures (listed in plan under Circularity Rules). These were used to fit the stripe model.

STRIPE MODEL (Mutter_02 2d, exploratory): The 15-deg tilt is tissue mounting angle, NOT beam angle. Conserve 1.02mm spacing (collimator-fixed), but RE-FIT tilt and offset per Mutter_02 slide using Cdkn1a/p21 only (H2AX unavailable). Expect weaker signal than Mutter_01 4h (p21 decays by 2d). Include a random-null baseline (100 permutations of angle/offset) to judge whether any detected stripes exceed chance. A negative result is informative — report it plainly.

$SYSTEM_RULES

Your final output MUST include results/signatures/SUMMARY.md answering the scientific intention." \
  --dangerously-skip-permissions \
  -n "spatial-sigs-executor-$TIMESTAMP" \
  < /dev/null 2>&1 | tee "$RESULTS/executor.log"

EXEC_EXIT=$?
echo "[$(date)] Executor exit code: $EXEC_EXIT" | tee -a "$RESULTS/orchestrator.log"

# Executor validation gate
if validate_executor_output; then
  echo "[$(date)] Executor gate: PASSED" | tee -a "$RESULTS/orchestrator.log"
else
  echo "[$(date)] Executor gate: FAILED — skipping reviewer and reviser" | tee -a "$RESULTS/orchestrator.log"
  kill $WATCHDOG_PID 2>/dev/null || true
  echo "[$(date)] === ABORTED (executor gate) ===" | tee -a "$RESULTS/orchestrator.log"
  exit 1
fi

# ============================================================
# Phase 2: REVIEWER (2h max)
# ============================================================
echo "[$(date)] === PHASE 2: REVIEWER ===" | tee -a "$RESULTS/orchestrator.log"

timeout 2h claude -p "You are an ADVERSARIAL REVIEWER of an MBRT spatial transcriptomics analysis.

An executor agent processed Mutter_02 CosMx data (4 slides, ~1000-gene panel, 2-day timepoint ONLY, 3 samples per slide: 0 Gy / MBRT 2d / SBRT 20 Gy 2d, NO H2AX ground truth), integrated it with existing Mutter_01 analysis, and tested whether MBRT peak/valley transcriptomic signatures (derived from spatially validated Mutter_01 4h data) persist at 2d and replicate in independent biological samples.

Read these files:
- ~/repos/spatial-rads/plan-mbrt-signatures.md (the plan — see Reviewer Checklist section)
- ~/repos/spatial-rads/results/signatures/SUMMARY.md (findings)
- ~/repos/spatial-rads/results/signatures/plots/ (all figures)
- ~/repos/spatial-rads/results/signatures/data/ (all data tables)
- ~/repos/spatial-rads/results/signatures/executor.log (scan for errors, skipped steps)
- ~/repos/spatial-rads/dev/peak_valley_analysis/data/pv_degs_bulk.tsv (original signature source)

Find problems. Be ruthless and specific.

Key concerns to investigate:
- Circularity: are zone-classification genes in the scoring signature? All 18 DDR genes excluded?
- Cell-type-stratified signatures: used as primary (not bulk)? Could apparent persistence be a cell-composition artifact?
- UCell primary: is UCell the primary scorer? Is AddModuleScore reported as a secondary comparison?
- Batch effects: did Mutter_01 vs Mutter_02 integration distort biology?
- Pseudobulk independence: are FOVs from the same slide treated as independent biological replicates? (They should NOT be — they are technical, not biological, replicates.)
- FOV-to-sample inference: was the 3-samples-per-slide recovery method documented and sanity-checked?
- H2AX-absent stripe caveat: are Mutter_02 stripe results framed as exploratory/pattern-detection only, never validated peak-vs-valley biology?
- p21-only stripe fitting at 2d: was a random-null baseline used? Was the weaker-signal expectation acknowledged?
- Interpretive boundary: was all non-Mutter_01-4h data framed as 'scoring' not 'classification'?
- Multiple testing: FDR correction on all gene-level tests?
- n=1 caveat: did any Mutter_01 comparisons use formal p-values?

Write numbered findings (R1, R2, ...) to ~/repos/spatial-rads/results/signatures/REVIEW.md. Commit.

$SYSTEM_RULES" \
  --dangerously-skip-permissions \
  -n "spatial-sigs-reviewer-$TIMESTAMP" \
  < /dev/null 2>&1 | tee "$RESULTS/reviewer.log"

REV_EXIT=$?
echo "[$(date)] Reviewer exit code: $REV_EXIT" | tee -a "$RESULTS/orchestrator.log"

# Reviewer validation gate
if validate_reviewer_output; then
  echo "[$(date)] Reviewer gate: PASSED" | tee -a "$RESULTS/orchestrator.log"
else
  echo "[$(date)] Reviewer gate: FAILED — skipping reviser" | tee -a "$RESULTS/orchestrator.log"
  kill $WATCHDOG_PID 2>/dev/null || true
  echo "[$(date)] === ABORTED (reviewer gate) ===" | tee -a "$RESULTS/orchestrator.log"
  exit 1
fi

# ============================================================
# Phase 3: REVISER (3h max)
# ============================================================
echo "[$(date)] === PHASE 3: REVISER ===" | tee -a "$RESULTS/orchestrator.log"

timeout 3h claude -p "You are the REVISER agent.

Two agents already ran on an MBRT spatial transcriptomics signature analysis:
1. EXECUTOR: processed Mutter_02, scored signatures, wrote SUMMARY.md
2. REVIEWER: adversarial critique with numbered findings in REVIEW.md

Read:
- ~/repos/spatial-rads/plan-mbrt-signatures.md (Reviser Protocol section)
- ~/repos/spatial-rads/results/signatures/SUMMARY.md
- ~/repos/spatial-rads/results/signatures/REVIEW.md

For each R# in REVIEW.md:
- Valid + fixable: fix it. Re-run the affected R script, update figures/tables/SUMMARY.md.
- Valid + unfixable: add to Limitations in SUMMARY.md with impact assessment.
- Wrong: rebut under 'Reviewer Response' section with evidence.

After making fixes, verify that previously-passing validation still passes:
- SUMMARY.md still >500 words and answers the scientific intention
- All referenced figures still exist
- No new R errors introduced

Add a 'Revision Notes' section to SUMMARY.md documenting every R# and what was done.
Commit when done.

$SYSTEM_RULES" \
  --dangerously-skip-permissions \
  -n "spatial-sigs-reviser-$TIMESTAMP" \
  < /dev/null 2>&1 | tee "$RESULTS/reviser.log"

REVIS_EXIT=$?
echo "[$(date)] Reviser exit code: $REVIS_EXIT" | tee -a "$RESULTS/orchestrator.log"

# ============================================================
# Done
# ============================================================
kill $WATCHDOG_PID 2>/dev/null || true

echo "[$(date)] === COMPLETE ===" | tee -a "$RESULTS/orchestrator.log"
echo "Exits: executor=$EXEC_EXIT reviewer=$REV_EXIT reviser=$REVIS_EXIT" | tee -a "$RESULTS/orchestrator.log"
echo "Branch: $BRANCH" | tee -a "$RESULTS/orchestrator.log"
echo "Results: $RESULTS/"
echo ""
echo "Next: review with review-mbrt-signatures.md"
