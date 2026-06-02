# Retro: MBRT Peak/Valley Signature Persistence and Replication

Pilot retrospective for the virtual-scientist **Step 4b** workflow. Companion
to `plan-mbrt-signatures.md`, `review-mbrt-signatures.md`, and
`run-mbrt-signatures.sh`. Filled out during the Step 4 review session; the
**Action items** section is the payload — everything above it is evidence.

Run: `autonomous/mbrt-signatures-20260410-1540`
Launcher start: 2026-04-10 15:40:09 UTC
Launcher end: 2026-04-10 23:40:09 UTC (phase 1 SIGKILL by `timeout 8h`)
Terminal state: **FAILED** — phase 1 hit wall-clock timeout. Phases 2 and 3 never ran.
VM: jeff-frag-test (48 cores, 613 GiB RAM)

---

## Wall clock

- Expected (from plan): not formally budgeted. Plan used the skill default `timeout 8h` for phase 1 without justification. No wall-clock estimate existed in the plan's Executor Tasks.
- Actual per phase:
  - Phase 1 (executor): **8h 00m** — SIGKILL'd by `timeout 8h` at 23:40:09 UTC exactly (start 15:40:09 UTC). Executor did not voluntarily exit.
  - Phase 2 (reviewer): **never ran** — launcher's `set -e` + timeout exit code 124 aborted the script immediately.
  - Phase 3 (reviser): **never ran** — same reason.
- Largest time sink: `01_process_mutter02.R`.
  - Attempt #1: 52 min, crashed at FindNeighbors on `future.globals.maxSize` (see Failures and recoveries).
  - Attempt #2: **5h 48m and counting when SIGKILL hit.** Ran from 17:44 → 23:40 UTC. Never reached `mutter02_qc_summary.tsv` (line 232 of the script) or `saveRDS(mutter02_processed.rds)` (line 239). All compute lost — no checkpoint.
  - Script 03 (signature definition on Mutter_01): completed successfully at 16:01 UTC (well under 20 min). These are the only real deliverables this run produced.

## Resource utilization

- Peak CPU cores used vs available:
  - Attempt #1: **1.0 core effective / 48 available** (2%). Single-threaded throughout.
  - Attempt #2: peaked at ~1.4 cores effective (CPU:wall ratio 1.37×), declined to **0.98× over the final 38 min** as the script entered single-threaded phases (UMAP finalize, TransferData coordinator). Over 20 R worker threads were alive but most were blocked on coordination, not computing. **Never exceeded ~3% of the box's CPU capacity.**
- Peak RSS vs available: 86 GB / 613 GB (**14%**). Never close to memory limits.
- Swap: **zero pages used** throughout. Not a memory-pressure failure.
- Scratch high-water: **0 bytes** in `/mnt/data/spatial-rads/signatures/`. The merged-object pipeline never wrote a single intermediate checkpoint — everything was held in the R session's RAM and died with it.
- Under-utilized windows: **the entire run.** At no point was the VM resource-constrained. We paid for a 48-core / 613 GiB box and used ~3% of the CPU and 14% of the RAM, then ran out of wall-clock time. This is the single most important data point in the retro: **the failure mode was scheduling, not capacity.**

## Agent decisions worth revisiting

(Not wrong, just suboptimal. Add entries live during the run.)

- **2026-04-10 16:27–17:43 UTC — `future.globals.maxSize` crash at FindNeighbors after 52 min of compute.**
  `01_process_mutter02.R` ran at 100% CPU single-core (RSS 20–22 GB) on a
  48-core, 613 GiB box. **Reached FindNeighbors on merged ~2.8M-cell object
  after ~52 min of per-slide QC/Normalize/ScaleData/PCA work**, then died
  with:
  ```
  Computing nearest neighbor graph
  Error in getGlobalsAndPackages(expr, envir = envir, globals = globals) :
    The total size of the 7 globals exported for future expression
    ('FUN()') is 1.17 GiB. This exceeds the maximum allowed size of ...
  ```
  The script had `options(future.globals.maxSize = 8 * 1024^3)` at the top,
  but the crash suggests a scoping / plan-override issue: `library(UCell)`
  pulls in `future` and may flip the plan; Seurat's internal FindNeighbors
  can launch its own future expression that doesn't honor the option set in
  the top-level session. **52 min of wall clock lost to a retryable error.**

  **Executor recovered autonomously**: caught the error in the tool-result
  stream, edited the script (9042 → 9096 bytes, mtime 17:43), relaunched as
  PID 119574 at 17:32 UTC. Live at check time, 11 min in, still in the
  per-slide loop.

- **2026-04-10 16:27+ UTC — Single-threaded Seurat on ~2.8M cells (background tax, attempt #1).**
  Independent of the crash: the agent did not set
  `future::plan("multicore", workers=N)` before the heavy Seurat steps.
  Attempt #1 ran single-threaded end-to-end on a 48-core box. The executor
  recognized this for attempt #2 and added multicore, but the fix came with
  its own problem (next bullet).

- **2026-04-10 17:44–23:40 UTC — Attempt #2 multicore collapse.**
  For the retry, the executor set `future::plan("multicore")` and launched
  20 R worker threads. CPU utilization started at 1.37× effective (~140%),
  held at 1.22× through the per-slide loop, then **collapsed to 0.98× over
  the final 38 min of wall clock**. The multicore advantage evaporates when
  Seurat hits phases that don't parallelize (`RunUMAP`, TransferData's k-NN
  projection step, etc.). Over-allocation of workers doesn't help and adds
  coordination overhead that can *slow* the single-threaded sections. The
  real lesson: **merging all four Mutter_02 slides into a 2.8M-cell
  monolith was the wrong architecture**, not "use multicore." Per-slide
  `TransferData` from Mutter_01 as a reference would have kept every
  operation on ~700K cells per slide — tractable with or without multicore.

- **2026-04-10 15:40–23:40 UTC — 8h phase 1 timeout was miscalibrated for
  the work.** The plan used the skill's default `timeout 8h` with no
  justification. A 2.8M-cell merged Seurat pipeline with
  `FindNeighbors`/`RunUMAP`/`TransferData` needs more than 8h even when
  parallelized well, and the plan knew the object size in advance. The
  timeout was set blindly. Correct approach: the plan should compute an
  expected phase 1 wall clock from the data size and the slowest known
  operations, then size the timeout at ~2× that estimate with a floor of 6h
  and a ceiling that matches the day's budget.

- **2026-04-10 15:40–23:40 UTC — No intermediate checkpoints in the R
  script.** `01_process_mutter02.R` ran as a single monolithic Rscript
  with one `saveRDS` at the very end (line 239). When the timeout hit, all
  5h 48m of attempt #2 compute — per-slide QC, normalize, scale, PCA,
  merge, and a deep portion of the merged-object NN/UMAP/TransferData work
  — evaporated. A checkpoint after per-slide processing, and another after
  merge, would have let a rerun resume from the midpoint instead of
  starting over. This is a generic lesson for any multi-hour executor
  script: **checkpoint at every stable intermediate state**, not just at
  the end.

## Failures and recoveries

- **2026-04-10 17:43 UTC — Executor caught `future.globals.maxSize` crash and
  retried autonomously.** `01_process_mutter02.R` died at FindNeighbors after
  52 min of compute with a `getGlobalsAndPackages` error. Executor diagnosed
  from tool-result stderr, edited the script (9042 → 9096 bytes), relaunched
  in background. No nanny intervention needed. This is the executor prompt
  rule #2 ("root cause before retry") working as designed — cost: 52 min
  wall clock.
- **2026-04-10 17:53:29–17:59:05 UTC — CLI retry storm (NOT a single
  transient — CORRECTED post-run).** At the time I reported to the user
  that there was "one transient `UND_ERR_SOCKET` around 18:00" and that CLI
  internal retry handled it. **That was wrong.** Post-mortem jsonl
  inspection shows **10 consecutive `subtype:"api_error"` events** between
  17:53:29 and 17:59:05 UTC, with `retryAttempt` climbing from 7 through 10
  out of `maxRetries: 10`. The CLI came within a single attempt of giving
  up entirely. My Monitor filter only surfaced the assistant-text tail of
  the storm (the one event at 18:00), not the earlier `api_error` system
  messages. **This is the most important nanny lesson of the run**: the
  jq filter must match `subtype == "api_error"` and surface *all* of them,
  not just assistant-text errors. A single error at retry 10/10 is not a
  "one transient blip" — it is the last gasp of a six-minute failure burst
  I failed to see.
- **2026-04-10 23:40:09 UTC — SIGKILL by `timeout 8h`. Terminal failure.**
  Phase 1 hit wall clock, claude PID 88592 died, R worker 119574 died
  with it, launcher 199660 exited on `set -e`. `orchestrator.log` ends
  after the phase 1 start banner because no further line was ever written
  — the script was killed mid-line in the phase 1 invocation. Phases 2 and
  3 never ran. 5h 48m of attempt #2 compute lost with zero intermediate
  output. This is the run's actual terminal state.
- Retry loops triggered in `run_claude_phase`: 0 (the launcher-level retry
  is for CLI exits; CLI-internal retries are handled by the CLI itself and
  don't bubble up to the launcher).
- Interventions needed by nanny: 0 actual; 1 attempted. The user asked for
  "kill claude externally, preserve R" (option 2) at ~23:44 UTC, but by
  then the SIGKILL had already fired at 23:40. My timing arithmetic was
  off by ~5 minutes and the window to act closed during deliberation. See
  Nanny process lessons.

## Plan-template gaps

Things `plan-mbrt-signatures.md` should have specified but didn't.

- **No phase 1 wall-clock estimate.** Plan inherited `timeout 8h` from the
  skill default without sanity-checking it against the known data size
  (2.8M merged cells, 4 Seurat objects to process, FindNeighbors + UMAP +
  TransferData all on the merged object). A plan that knows its input
  dimensions must size its timeout. Candidate rule: every plan's Executor
  Tasks section must include a line `Expected phase 1 wall clock: N hours
  (justification: ...)` and the launcher's `timeout` must be ≥ 2× that
  estimate.
- **No "merge vs per-slide" architectural decision.** Plan did not specify
  whether to process Mutter_02 as 4 independent slides or as a merged
  object. Executor chose merge, which was the expensive path and also the
  unnecessary one (per-slide `TransferData` from Mutter_01 reference is
  scientifically equivalent for signature projection and computationally
  tractable). Plans for scRNA / spatial work on multiple samples must
  explicitly state the unit-of-processing and why.
- **No checkpoint requirement.** Plan did not require the executor to
  write intermediate RDS checkpoints in long-running scripts. A checkpoint
  after per-slide processing and another after merge would have made
  attempt #2's 5h 48m recoverable instead of total loss. Candidate
  executor prompt rule (see below).
- **No separation of "cheap/fast sanity" and "expensive/full" passes.** The
  Mutter_02 processing script went directly to full-object processing with
  no smoke test on one slide first. The virtual-scientist skill already
  mandates a single-sample smoke test gate in Step 2 (setup), but that's
  for the data staging phase — not for the analysis script itself. A
  downstream smoke test (run the script on one slide, confirm
  saveRDS-worthy output, *then* run on all slides) would have caught the
  merged-object blowup before it consumed phase 1.

## Skill-body and module gaps

Things `virtual-scientist` skill and its org-mode module
(`private.org` id `e243c532-65ab-4434-a844-ab386869a47a`) should teach but
don't.

- **Step 4b itself** — this retro is the pilot; the skill does not currently
  have a workflow retrospective step. Step 4 is scientific-review only, and
  workflow lessons evaporate between runs. (Meta: the very act of writing this
  file is the first Step 4b action item.)
- **Parallelization preflight** — no guidance in Step 1 plan template or
  Executor prompt rules about matching worker counts to available cores. On a
  48-core / 613 GiB box, a single-threaded Seurat job burns wall clock
  silently. Candidate rule: "For any R / python worker that could be
  parallelized and will run on ≥1M rows/cells, set the parallel backend in the
  plan's Executor Tasks explicitly. For Seurat, that means
  `future::plan('multicore', workers=N)` + `options(future.globals.maxSize=...)`
  before FindNeighbors/TransferData/IntegrateLayers." **Caveat learned in
  attempt #2**: multicore is necessary but not sufficient. Seurat has large
  single-threaded sections (UMAP finalize, TransferData coordinator) where
  more workers don't help and coordination overhead can slow things down.
  The rule should also say: **prefer per-sample work to merged-object work
  whenever the science allows**, because per-sample work scales linearly
  and merged-object work hits single-threaded ceilings.
- **Phase 1 timeout sizing.** Skill currently offers `timeout 8h` as the
  canonical phase 1 wrapper without any sizing guidance. Candidate rule:
  "Phase 1 timeout = max(expected wall clock × 2, 6h, data-size-dependent
  floor). For Seurat/scanpy work on merged objects > 1M cells: start at
  12h. For > 2M cells: 16h or don't merge."
- **Retry-burst detection in the nanny jq filter.** The filter at
  `~/jsonl-filter.jq` only surfaces assistant-text and tool events. It
  silently drops `type:"system"` entries with `subtype:"api_error"`. A
  six-minute retry storm is invisible. Candidate filter addition: a branch
  for `type == "system" and subtype == "api_error"` that emits a timestamped
  retry-attempt line. Candidate nanny rule: "Three `api_error` events
  within 10 minutes is an automatic alarm, not a transient."
- **Checkpoint discipline in executor prompt rules.** New candidate rule #7:
  "Any R or python script expected to run more than 1 hour must save
  intermediate RDS / pickle checkpoints at every stable state: after data
  load, after per-sample processing, after merge, after each heavy
  downstream transformation. If a phase timeout SIGKILLs the script, a
  rerun must be able to resume, not restart."
- **Launcher behavior on phase 1 timeout.** Currently the launcher's
  `set -e` causes immediate exit when `timeout` returns 124, so phase 2
  and phase 3 never run even if phase 1 produced usable partial outputs.
  Candidate rule: the launcher should treat "timeout but deliverables
  exist" differently from "timeout and empty results directory." If
  `results/signatures/data/` is non-empty, run the validation gate anyway
  — it may still pass and allow the reviewer phase to salvage what's
  there. (In this run, the data/ directory had Mutter_01 signatures from
  16:01 UTC that the reviewer could have critiqued.)

## Action items (the payload)

Each action item is specific: which file, which section, what text.

**Skill / module updates (Step 4b pilot deliverables):**

- [ ] **Update virtual-scientist skill** (`~/.claude/skills/virtual-scientist/SKILL.md`):
  add **Step 4b — Workflow retrospective** section with the `retro-<name>.md`
  template. Amend the Step 1 artifacts list to include `retro-<name>.md`.
- [ ] **Update virtual-scientist module** (`private.org` id
  `e243c532-65ab-4434-a844-ab386869a47a`): mirror the Step 4b addition so the
  conceptual home and the operational skill stay in sync.

**Executor prompt rules (add to the rule block the launcher injects):**

- [ ] **Rule: parallelization preflight.** For R / python work on ≥1M
  rows/cells on a multi-core box, set the parallel backend explicitly at
  the top of the script, AND prefer per-sample work to merged-object work
  when the science allows (per-sample parallelizes cleanly; merged-object
  hits single-threaded ceilings).
- [ ] **Rule: intermediate checkpoints (new rule #7).** Any script
  expected to run > 1 hour must save intermediate `saveRDS` / `pickle`
  checkpoints at every stable state (after data load, after per-sample
  processing, after merge, after each heavy transformation). A timeout
  SIGKILL must leave a resumable state on disk.
- [ ] **Rule: smoke-test before full pass.** Before running a heavy
  analysis script on all samples, run it on one sample end-to-end and
  verify a non-degenerate checkpoint. This applies inside the executor
  phase, not just to the Step 2 data-staging smoke test.

**Plan template (virtual-scientist Step 1 plan artifacts):**

- [ ] **Require an explicit phase 1 wall-clock estimate** in every plan's
  Executor Tasks section, with justification, used to size the launcher's
  `timeout` wrapper.
- [ ] **Require a "merge vs per-sample" decision statement** in every
  plan that processes multiple scRNA / spatial samples.
- [ ] **Require a checkpoint plan** (named intermediate files, disk
  budget for them) in every plan with an expected phase 1 > 2h.

**Launcher script conventions (skill's `run-<name>.sh` guidance):**

- [ ] **Phase 1 timeout sizing rule.** Default from `timeout 8h` to
  `timeout 12h` for merged-object scRNA / spatial work > 1M cells.
  Canonical example in the skill body.
- [ ] **Timeout-tolerant validation gate.** When phase 1 times out,
  inspect `results/*/data/` — if non-empty, run the validation gate
  anyway. Allow partial-success handoff to phase 2. Code sketch in skill.

**Nanny tooling (jq filter + Monitor discipline):**

- [ ] **Update `~/jsonl-filter.jq`** to emit `type:"system" +
  subtype:"api_error"` events with retry attempt numbers. File lives at
  `~/.claude/skills/virtual-scientist/jsonl-filter.jq` (staged to VM at
  `~/jsonl-filter.jq` during Step 2).
- [ ] **Nanny rule: three `api_error` events within 10 min = alarm.** Add
  to skill's nanny section.
- [ ] **Nanny rule: elapsed-vs-timeout tracker.** On every status check,
  compute `claude elapsed time` vs the phase timeout and surface "time
  remaining" to the user. If < 30 min, warn; if < 15 min, intervene
  first and ask second.

**Project memory:**

- [ ] Add `feedback_vm_runs.md` capturing the combined rules:
  Seurat-multicore + per-sample preference + checkpoint discipline +
  timeout sizing + elapsed tracker. Reference this retro as the source
  incident.

**Rerun-specific action items:**

- [ ] **Rewrite `plan-mbrt-signatures.md`** for a rerun:
  - Per-slide `TransferData` from Mutter_01 reference (no merged object).
  - Checkpoint after each slide's processing (`mutter02_slide_{1..4}_processed.rds`).
  - Explicit phase 1 budget: 12h.
  - Use `future::plan("multicore", workers=8)` — not 20. Over-allocation
    hurt attempt #2.
  - Smoke test: process slide 1 fully, verify `saveRDS` output, then loop.
- [ ] **Update `run-mbrt-signatures.sh`**: `timeout 12h` on phase 1.
- [ ] **Keep the Mutter_01 signatures from 16:01 UTC** — those are valid
  deliverables (Task 3 of the original plan) and should be cherry-picked
  to the rerun's feature branch rather than recomputed.

## Nanny process lessons (my-side)

- **Attach `Monitor` to the jsonl at the first status check, not the fifth.**
  I polled the VM via one-shot `gcloud compute ssh` calls for ~2 hours before
  switching to `Monitor` with a `tail -F` + jq filter on the current session
  jsonl. Every poll was an IAP SSH round trip; the Monitor-based stream is
  push-based and free. The virtual-scientist skill already documents the jq
  filter and the jsonl location, but does not explicitly tell the nanny to
  set up Monitor as the **first** action of a status check. Candidate rule:
  "On the first nanny check of a run, stage `~/jsonl-filter.jq` on the VM
  and launch `Monitor` with `tail -n 0 -F <current jsonl> | jq -f
  ~/jsonl-filter.jq`. Subsequent checks become targeted one-shots (proc
  state, disk, file listings) and the turn stream arrives free."
- **The jq filter was too narrow.** I used a filter that only emitted
  assistant-text and tool events. It silently skipped `type:"system" +
  subtype:"api_error"` messages. That's why I saw "one transient socket
  error at 18:00" instead of the actual six-minute, 10-attempt retry
  storm between 17:53 and 17:59. **Every status update I gave the user
  about the API error was wrong because the filter lied by omission.**
  Action item already listed above: update the canonical filter to emit
  `api_error` events too.
- **Elapsed-vs-timeout was not on my status dashboard.** I tracked
  "claude elapsed time" as a curiosity but never compared it against the
  8h phase 1 timeout until hour 5+. A simple calculation — "claude has
  been up 7h 51m, timeout is 8h, you have ~9 minutes before SIGKILL" —
  was absent from every status update I gave the user until the very
  last one. Candidate dashboard line for every nanny status:
  `phase 1: NhNm / 8h (Nm remaining)`.
- **Timing arithmetic failure under time pressure.** At 23:31 UTC I told
  the user "~7 minutes left" before SIGKILL, then spent those minutes
  investigating the launcher's retry logic and process tree. The SIGKILL
  fired at 23:40:09 UTC while I was mid-investigation; the user's "option
  2" call at 23:44 couldn't be acted on because there was nothing left
  to kill. **Rule: when the timeout clock is under 15 minutes, intervene
  first and investigate second.** A premature kill of the R worker would
  at worst have been a clean cancel; waiting to understand the tree
  meant losing the decision window entirely.
- **I conflated "jsonl frozen" with "bash tool blocked on long Rscript."**
  That interpretation was charitable and also correct for much of the
  run, but it masked the fact that I had no independent signal of R
  progress other than CPU burn rate and open libraries. Next time,
  attach `strace -p` or at minimum `py-spy`-equivalent to see what the R
  process is actually doing, or insist on intermediate file writes from
  the script itself.
- **I under-estimated how slow merged-object Seurat is.** At the first
  5h status check I still believed `mutter02_processed.rds` was "due in
  60-90 minutes." It wasn't — the script never produced it. A rougher,
  more cynical estimate ("merged-object FindNeighbors/UMAP on 2.8M cells
  with partial multicore can easily take 8+ hours") would have triggered
  an intervention recommendation 3+ hours earlier.

## Cross-project harvest notes

_For later use when ≥2 projects have retros. Glob
`~/repos/*/retro-*.md`, look for action items that repeat across projects,
and promote those to skill/module updates._
