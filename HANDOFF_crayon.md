# ClrAttention:  handoff to the Linux workstation (crayon)

*First written 2026-10-13 by the cloud session that ran the project on the Mac.  Updated 2026-10-13 22:05 PT by the session that moved it to crayon, for a Claude Code session running locally on crayon.  The master copy is the claude.ai project doc `claude/clr_handoff_workstation.md`.*

## State as of 2026-10-13 22:05 PT

- **Code.**  The local clone is at `~/Projects/ClrAttention`, on branch `accuracy-review`.  Commit 25318d1 (the job-runner changes for two machines) has been applied as 6ac82d1.  A further commit, cbaa2bc, adds `jobs_pc/` and `runner_pc.out` to `.gitignore`.  Both are pushed to GitHub `rimorob/ClrAttention`, whose `accuracy-review` head is cbaa2bc.  In the clone, the push credentials come from `gh` (the repo-local `credential.helper` is `!gh auth git-credential`), and the repo-local git identity is Claude <noreply@anthropic.com>, matching earlier commits.
- **Decision record.**  `CITATION_LOG.md` (D1–D29) holds the decisions.  `docs/COMPUTE_CONSTRAINTS.md` lists every choice made to fit the Mac.  `docs/whitepapers/*.md` are the findings logs.
- **Mac results.**  The Mac finished every queued run on 2026-09-24.  Its outputs are in `results/` in the local clone (copied from iCloud):  `chips_bootstrap` (200 cluster and 100 replicate draws), `beeline_v2` (graded ChIP and filtered STRING truths), `perturbation_B100v2` (variance computed without the perturbation's own experiments) and `tfnode_attention`.  The white papers have not yet been updated with `beeline_v2`, `perturbation_B100v2` or the extended chips bootstrap.
- **Machine.**  crayon has two Xeon E5-2690 v4 CPUs (28 physical cores, 56 threads), 125 GB of RAM, R 4.6.1, and `gh` logged in as rimorob.  `analysis/parallel.R` sizes the workers automatically.  With the 3.3 GB chips budget, that comes to about 26 workers without hyper-threading and about 31 with it.

## How the jobs run (updated 2026-10-13 22:30 PT)

Everything runs from the local clone `~/Projects/ClrAttention`.  Nothing reads from or writes to iCloud while jobs run, to avoid sync races and network errors.  `data/` and `results/` were copied once from `~/iCloud/Prototypes/ClrAttention/`, and that folder is now only a stale archive.  The old iCloud runner and the `210_bootstrap_pc` job were stopped.

1. **The runner.**  Run it as a capped user service, so it survives editor sessions and logouts (lingering is enabled) and a runaway job cannot trigger the global OOM killer:  `systemd-run --user --unit=clr-jobs --working-directory=$HOME/Projects/ClrAttention -p MemoryMax=95G -p OOMPolicy=continue -p CPUWeight=10 -p Nice=19 bash tools/job_runner.sh jobs_pc`.  OOMPolicy=continue keeps the runner alive when one job hits the cap; CPUWeight and Nice keep the desktop responsive when jobs use every hardware thread.  No MemoryHigh: its reclaim pressure tripped systemd-oomd on 2026-10-16, which killed both lanes and the desktop session.  A second, PIDC-only lane runs Julia directly (NetworkInference.jl on Julia 1.6 segfaults from R and with threads): `--unit=clr-jobs2 -p MemoryMax=10G ... jobs_pc2`; its edges are scored later in the main lane with `baselines.R --methods pidc`.  Keep the caps' sum at or under about 105 GB, size `--workers` from measured peaks, and run no heavy interactive tests while jobs run.  Check it with `systemctl --user status clr-jobs`, or watch both lanes with `watch -n 30 tools/progress.sh` (progress bars; expected minutes per job in `<lane>/estimates.tsv` or a `# expect:` line, log-based steps via `# progress: <total> <regex>`).  It also starts the control channel (`jobs_pc/control/`).  Do not hand-start jobs alongside it, and never edit a script in place while a job runs it.
2. **The `jobs_pc` queue.**  The jobs run in name order.  Each one sources `jobs_pc/common.sh`, which changes into the repo root (relative to its own location), sets `JOBS_ROOT=jobs_pc`, and loads `CLR_HT` from `jobs_pc/ht.env`.  Results stay in `results/`.
   - `010_setup.sh` runs `tools/setup_linux.sh`.
   - `020_bench_ht.sh` runs `tools/bench_ht.sh`.  It writes `CLR_HT=1` to `jobs_pc/ht.env` only if hyper-threading gives at least 5% more draws per hour, and `CLR_HT=0` otherwise.
   - `030_chips_bootstrap_1000.sh` runs `chips_bootstrap.R --Bcluster 1000 --Brep 200`, which resumes from the Mac's finished draws in `results/chips_bootstrap`.
   - `040_perturbation_B500.sh` runs `perturbation.R --B 500 --out results/perturbation_B500v2`.
   - `050_avg_null_B1000.sh` runs `run_m3d_regulondb.R --set avg --B 1000 --threads (nproc - 4) --out results/avg_B1000`.

**Fixes made on crayon (2026-10-13).**  `~/.bashrc` exported `R_LIBS` pointing at the R 4.4 library, whose Rcpp does not load under R 4.6.  The line is commented out (backup `~/.bashrc.bak-2026-10-13`), and `jobs_pc/common.sh` also unsets `R_LIBS`.  `analysis/parallel.R` now counts physical cores from `/proc/cpuinfo` (commit 9e6c87c), because `detectCores(logical = FALSE)` returned 56 here instead of 28.  The first benchmark run had used the wrong count and was discarded.

## Checks

1. `jobs_pc/runner.log` exists, and `ps -eo pid,args | grep -E "job_runner|control_runner"` shows the `jobs_pc` runner and its control channel.
2. Watch `jobs_pc/done/010_setup.log`.  If R packages fail to build, the system libraries in the header of `tools/setup_linux.sh` may be missing, and installing them needs sudo.
3. Do not queue the five jobs above again.

## Queue after the first five jobs (added 2026-10-14)

All five original jobs finished.  The remaining COMPUTE_CONSTRAINTS analyses were queued in priority order (`jobs_pc/queue/`, one at a time):  060 lab splits and cap sweep (`analysis/lab_splits.R`), 070 nested chips bootstrap with TF-node metrics (`analysis/nested_bootstrap.R`, 1,000 draws), 080 lab-matched perturbation null (`perturbation.R --null lab`), 085/087 BEELINE with 5,000-gene and all-gene CLR context (`beeline.R --context`), 090 BEELINE nested cell bootstrap in the 5,000-gene context (`analysis/beeline_boot.R`), 095 SSEM-Lasso (`perturbation.R --ssem 1 --reuse 1`), 100 FDR operators in the chips bootstrap (`chips_bootstrap.R --fdr_B 100`, 200 draws).  The user's guidance:  compute is effectively free and a week is fine; size runs for what they need to settle.  The draft paper and journal proposal are in `docs/whitepapers/paper_draft.md`.

## What is left

- **Analyses that need new code** (from `docs/COMPUTE_CONSTRAINTS.md`):
  - re-selecting the depth inside each bootstrap draw (a nested bootstrap);
  - 20–50 lab splits, and cell splits for BEELINE, instead of one;
  - full gene context for BEELINE, with the 5,000 most variable genes as the intermediate step;
  - a perturbation null matched by lab or condition;
  - bootstraps for BEELINE and the TF-node benchmark;
  - a softmax-cap sweep (caps of 100 and 200, and no cap);
  - the FDR-thresholded operators inside the chips bootstrap.
- **White papers.**  Update them with `beeline_v2`, `perturbation_B100v2` and the extended chips bootstrap, then with the full-scale runs as they finish.
- **Results through git.**  `results/` is git-ignored.  Commit result summaries, not the large `.rds` fits, and push from crayon.
- **Validation.**  Have a separate QC agent that did not produce the analysis check each new result before it goes into a white paper.

## Model choice (agreed with the user, 2026-10-16)

Opus at high effort for analysis design, result interpretation, independent QC reviews and paper writing; Opus at medium as the default; Sonnet for routine work (queue babysitting, status tables, plotting from a settled spec, boilerplate scripts, mechanical ports of designed code). Apply the same split when spawning subagents.

## Scheduler (2026-10-20): the standing rule is ~70-80% machine utilization

`tools/scheduler.sh` (service `clr-sched`, queue `jobs_pc/sched/`) replaces the one-at-a-time runner. Each job declares `# mem: <GB>`, `# cpus: <N>` and optionally `# after: job1,job2`; jobs run in parallel within 90 GB / 50 threads with at least 15 GB left free, all inside `clrjobs.slice` (`~/.config/systemd/user/clrjobs.slice`: hard cap 95 GB, no swap, which protects the desktop), each in its own scope `clrjob-<name>` capped at min(1.3 x declared, 95 GB) with OOMPolicy=kill. Restart the scheduler with `-p Restart=on-failure`; running jobs survive a restart and are finalized from `running/<job>.rc`. Stop a job with `systemctl --user stop clrjob-<name>.scope` (stops its R workers too). The END line in `jobs_pc/runner.log` records the measured peak; use it to recalibrate declarations. Watch with `watch -n 30 tools/progress.sh` or the auto-refreshing `docs/status/progress.html`.
