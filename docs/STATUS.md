# Project status

*Thursday 22 October 2026, 21:45, before a second shutdown (GPUs back in place and visible, both x16, both report NUMA node 0 where node 1 was expected). Living file: update at each checkpoint. Everything is committed and pushed; the jobs resume from their caches (GRNBoost2 8 runs, K562 keff splits 3-4 saved, BEELINE 17%, RPE1 t\* splits in `results/clr_logpt_test/rpe1/splits/`). The CLR(log P^t) checklist below says where to continue; the next step is the EcoCyc truth build. GPU work is set up (`~/venvs/scfm`, checkpoints hashed in the scfm preregistration change log) but the lamindb setup, leakage audit and smoke test are not started: run them only when the CPU jobs are light.*

## Restarting after the shutdown

1. The scheduler (`clr-sched`, a boot-time user service; lingering is on) starts by itself, finds the interrupted jobs in `jobs_pc/running/` without exit codes, and puts them back in `jobs_pc/sched/` (`REQUEUE` lines in `jobs_pc/runner.log`). All three resume from their caches.
2. Check: `systemctl --user status clr-sched`, `tools/progress.sh`, `grep -E "REQUEUE|START" jobs_pc/runner.log | tail`.
3. GPU check after the move: `nvidia-smi --query-gpu=index,name,pci.bus_id --format=csv` and, per card, `cat /sys/bus/pci/devices/<bus>/current_link_width` (expect 16) and `numa_node` (expect 1 for the moved card).
4. Start the CORUM order-of-operations expert (planned for noon; see "Next").

## Jobs

| Job | Data | What it tests | State at shutdown | After restart |
|---|---|---|---|---|
| 0893b keff selections, K562 | Human K562 Perturb-seq pseudobulk | Preregistered keff test, phase 1: held-out depth per width variant (no truths) | 2/10 splits done (cached) | resumes; about 1 day at 2 workers |
| 0894 GRNBoost2 diagnostic | hESC_1000, mHSC-E_500 (BEELINE), E. coli M3D | Whether arboreto's early stopping (broken by scikit-learn ≥ 1.3) explains GRNBoost2's last place | 8 runs done (kept) | resumes; about 15 h left |
| 110c BEELINE bootstrap | Human scRNA-seq (hESC, hHep) + 5 mouse | Bootstrap over cells for CLR-NS vs CLR and \|r\| and attention, 5,000-gene context | 1,577/3,080 tasks (51%) | resumes at 5 workers; raise to 12 when the keff jobs end |

## Results in hand (all reviewed; details in `docs/paper_plan.md`)

- **Headline (E. coli M3D):** attention beats CLR by +0.022 AUPR [+0.019, +0.026], in 1,000/1,000 nested bootstrap draws; held-out depth t = 9.
- **Hub control:** the lead survives on strong-evidence truths.
- **Per regulon:** ties CLR and GENIE3.
- **Mechanism:** the walk trades pair-level resolution for module recovery. It loses on tight complexes (K562 CORUM 1.8× lower than CLR).
- **Hard FDR operators:** soft attention beats the best one in 200/200 paired cluster-bootstrap draws.
- **BEELINE single cell:** CLR is comparable to GENIE3 and PIDC. CLR-NS has the highest mean but was chosen in-sample, so it is not confirmatory. On human data the result is split.
- **keff selections:** M3D and EcoMAC are committed (`docs/keff_selections/`). Held-out keff selection runs to the widest width (40), as the review predicted.

## Documents

| Document | Path | State |
|---|---|---|
| White paper and supplement | `docs/whitepapers/paper.md`, `supplement.md` | updated 22 Oct, fact-checked |
| Method table | `docs/whitepapers/sota_table.md` | reviewed; Bonferroni-corrected ties |
| Abstracts (150-word Nature Machine Intelligence; Nature Communications) and bullet outline | `docs/whitepapers/abstract_draft.md`, `abstract_points.md` | drafted; register of Hayete, Gardner & Collins 2007 |
| Graphical abstract (code-built) and diffusion prompt | `docs/whitepapers/figures/graphical_abstract.png`, `graphical_abstract_diffusion_prompt.md` | redesigned after fresh-eyes reviews; panel 2 being rebuilt as the CLR criss-cross (PI request) |
| Checkpoint for Jim Collins | `docs/checkpoint_collins.md`; private page https://claude.ai/artifact/5sbGiidecw7ZqkjbBz4eBT | fact-checked; share from the page's menu |
| Contributions and AI-use disclosure | `docs/contributions.md` | draft; PI to confirm marked items |
| Preregistrations | `docs/preregistration_keff.md`, `docs/preregistration_scfm.md` | committed before any scoring |

## Open items as of 22 Oct, 16:45

- **CORUM fix found (exploratory):** CLR applied again to log Pᵗ rescues complexes (K562 CORUM 7.95× random vs CLR 6.39×) without hurting TF regulons. PI decisions recorded in `docs/plans/corum_analysis.md`. **Next step:** preregistered and committed (`docs/preregistration_clr_logpt.md`: confirmatory = E. coli M3D + EcoMAC, EcoCyc complexes superiority and RegulonDB SC non-inferiority vs attention at t = 9; RPE1 merged CORUM key secondary; K562 in-sample only). `analysis/clr_logpt_test.R` written (gates, a3r operator, K562 reproduction mode) and the K562 reproduction check PASSED (22 Oct: AUPR 0.105435083 vs 0.1054351; `results/clr_logpt_test/k562_repro/repro_check.csv`). Remaining: the label-free t\* phase (selections committed to `docs/clr_logpt_selections/`) and the EcoCyc truth build (commit each), code review, then score E. coli (truth modes are gate-checked stubs so far).
- **K562 keff selections:** a one-shot service (`clr-k562-boost`) switches the job to 4 workers as soon as splits 3 and 4 are saved; nothing is lost.
- **GRNBoost2 diagnostic:** early sign that no early stopping fixes GRNBoost2 (mHSC-E AUPRC ratio 2.75); wait for the full run, then a reviewer.
- **Thermal watch** (`clr-thermal`) runs until midnight; CPU 2 sits at 82-86 °C under load (no warnings).
- **Model use:** subagents default to Sonnet (`worker` agent; user setting); Opus high only for design, statistical interpretation and claims.

## CLR(log Pᵗ) confirmatory checklist (resume here after any shutdown)

Design is fixed in `docs/preregistration_clr_logpt.md` (PI chose option A, 22 Oct: E. coli is the only confirmatory set; K562 in-sample; RPE1 key secondary). Tick each step only when committed. Jobs run through `jobs_pc/sched` and resume from caches; check `tools/progress.sh` first.

- [x] Preregistration committed (90b0e3d)
- [x] `analysis/clr_logpt_test.R` with gates; K562 reproduction check passed (e891f99, 5881335)
- [x] Label-free t\* phase, M3D (a3r t\* 6) and EcoMAC (a3r t\* 5), committed; attention t\* 9 / 8 match K10; none below 4. RPE1 still running (job 0902, then commit `rpe1_*` files): `--mode selections`; commit to `docs/clr_logpt_selections/` (script committed; jobs 0900 m3d and 0901 ecomac done, 0902 rpe1 in `jobs_pc/sched`, results in `results/clr_logpt_test/<dataset>/`)
  - How to resume: if a job died, requeue `jobs_pc/sched/0902_clrlogpt_sel_rpe1.sh` (the scheduler restarts interrupted jobs; finished splits in `results/clr_logpt_test/<dataset>/splits/` are reused).
- [x] E. coli complex truth built and its SHA-256 committed (source: EBI Complex Portal 83333, the preregistered fallback, because EcoCyc proteins.dat needs a BioCyc login; no fallback to operon-included needed; M3D 1,150 and EcoMAC 1,148 positive pairs, 54 and 53 complexes; log `docs/clr_logpt_ecoli_truth_log.md`)
- [x] Scoring implemented in `analysis/clr_logpt_test.R` (c64d2a5: `--mode truth|boot|verdict`, jackknife, cluster bootstrap, Holm, outcome table; smoke-tested on synthetic data only, no truth read). Staged jobs (not queued) in `jobs_pc/staged/` (untracked: jobs_pc is gitignored): 0910 m3d, 0911 ecomac, 0912 boot, 0913 rpe1, 0914 k562, 0915 verdict
- [x] Code review done (`docs/clr_logpt_review_report.md`, 26 Oct; reviewer saw no E. coli result). Script ready; 5 fixes committed (bdcc05b human_truths.R not sourced, 825b801 Holm for the small-loss note + three-valued adopt, 93e5bc8 verdict inputs, 1f45cfc cache robustness); synthetic smoke test passed
- [x] Before queuing the staged jobs, all done 26 Oct: (1) RPE1 t\* job 0902 finished, `rpe1_*` selections committed (a3r t\*=12, attention t\*=9); (2) change-log entry in the preregistration (complex source, release, SHA-256, counts, jackknife seed stream, K562 t\* note); (3) `--mode repro` rerun at the new preregistration hash passed (AUPR 0.105435); (4) `after: 0914` added to 0915; headers of 0910-0915 checked against measured peaks
- [ ] Truth phase, once: E. coli, then RPE1, then K562 (in-sample). QUEUED 26 Oct in `jobs_pc/sched`: 0910 truth M3D (12 GB), 0911 truth EcoMAC (12 GB, after 0910), 0912 M3D bootstrap (14 GB, after 0910), 0913 truth RPE1 (32 GB, after 0911), 0914 truth K562 (16 GB, after 0913), 0915 verdict (after 0912, 0913, 0914). **How to resume:** the truth phase refuses reruns of a finished dataset; if a job dies, requeue the same script unchanged (caches resume from `results/clr_logpt_test/truth/<dataset>/cache`, the bootstrap per draw); never edit the preregistration or selections in between.
- [ ] Write-up: verdict per the preregistration's outcome table

If a job died mid-step, rerun that step; do not skip a gate. No truth may be read before the three boxes above it are ticked.


## CLR(log P^t) confirmatory verdict (26 Oct 13:41; script `analysis/clr_logpt_test.R`, outputs `results/clr_logpt_test/verdict/`, committed copies in `docs/clr_logpt_verdict/`)

**NOT ADOPTED: failure to confirm; attention stays the main method** (preregistered rule: P1 and P2 must each reach one-sided p < 0.025). Complex superiority vs attention at t = 9: M3D +0.0257 AUPR (jackknife SE 0.0182, p = 0.079), EcoMAC +0.0285 (SE 0.0174, p = 0.051). Both point the predicted way and exceed the +0.005 minimum, but neither passes. Regulon non-inferiority (P3, P4) passed (M3D +0.0043, EcoMAC +0.0033; Holm p < 1e-5). The M3D cluster bootstrap interval for P1 is [+0.013, +0.043] (rule 3 satisfied), which does not override the jackknife result: the two resamplings measure different variance (see `docs/plans/clr_logpt_interpretation.md`), and the preregistration made the bootstrap an extra condition, not an alternative. Reproduction check passed. S1-S5 are descriptive (Holm family only applies if adopted). Interpretation and write-up wording: to be reviewed (Opus) before any claim; K562 and RPE1 are in-sample or secondary and cannot rescue the verdict. K562 keff selections job 0893c finished and are committed (`docs/keff_selections/k562ess_selections.csv`); keff phase 2 on K562 is queued as 0896.

- **Interpretation memo (27 Oct): `docs/plans/clr_logpt_interpretation.md`.** It covers what may and may not be claimed, why the jackknife and the bootstrap disagree (truth-sampling vs experiment-sampling variance; the bootstrap is a necessary condition, not an override), the power gap (complex SE 0.018, about 9× the regulon SE, which was the only arm in the power analysis; minimum detectable effect about 0.05), exact paper wording, how S1–S5 may be mentioned (intervals only, no p-values) and what a follow-up would need (a new preregistration with a larger, untouched complex truth).

## Next

1. **CORUM order-of-operations expert (PI, 22 Oct), now including the asymmetry assessment (26 Oct).** Members of a complex share one latent variable, so all their pairs are strongly coupled. The softmax and the walk may spread attention across the clique and then into the wider translation program. The question is whether collapsing strongly coupled cliques into latent variables before attending, or attending on residuals after removing a clique's shared factor, fixes the loss. A max-effort expert analyses the CORUM results and interviews the PI (questions relayed in the session), starting from `docs/plans/corum_analysis.md` and `docs/plans/asymmetry_expert.md`.
   - **Asymmetry (PI, 26 Oct; columns as regulators or latents, rows as observed genes): partial pass, on synthetic data only** (`docs/plans/asymmetry_expert.md`; scripts in `analysis/synthetic_asymmetry/`; no real data or truth read).
     - It is not a complex fix. Mass conservation keeps the 1/(n−1) dilution under any row, column or Sinkhorn normalization. Column-CLR on log Pᵗ ties a3r on complexes (×1.01 / ×0.97) and costs regulons (×0.94 / ×0.86). Sinkhorn costs regulons ×0.74–0.85. Residual or units-first attention deletes regulon programs (×0.30). a3r (already preregistered) remains the complex fix.
     - Passing part: latent-column CLR (C3). Each column becomes the attention-weighted latent of its attenders (leave-one-out), and genes are scored against it with CLR's row and column z-scores. It beats a3r on direct TF–target edges by ×1.12–1.25 in 12/12 synthetic seeds; on complexes it scores ×0.68.
     - Direction is not identifiable from co-expression: there are no direction claims.
     - **Follow-up (memo section 6):** only after the CLR(log Pᵗ) verdict, preregister C3 as an edge-level readout on E. coli RegulonDB direct edges (M3D, EcoMAC; vs CLR and attention, plus a3r if adopted). Merge it with the planned directed DREAM5-style evaluation (item 6), so that the edge truth is touched once. t_a is chosen label-free; keff 10 is fixed. No C3 computation on real data before that.
2. **Figures:**
   - finish the revisions from the fresh-eyes reviews (Shapley, MI vs ρ, lexA scatters; one colour key and finding titles);
   - rebuild graphical-abstract panel 2 as the CLR criss-cross;
   - build the paper figures that have legends but no images (Figs 1, 2, 4–7).
3. **keff test, phase 2:** K562 truth scoring finished (job 0896 retry, peak 47.9 GB). M3D (0897), EcoMAC (0898, after 0897) and `--dataset combine` (0899, after both) are queued, 3 workers each, mem 30 (phase 1 peaks were 18-20 GB at 5 workers); each is behind the byte-identity gate, so if one dies, requeue the same script (never edit selections). Original plan: run `--dataset m3d` and `--dataset ecomac` truth phases with `--selections_ref docs/keff_selections/<dataset>_selections.csv`, then `--dataset combine`.
4. **Learned-attention comparator (scPRINT, scGPT):** preregistered. Install a separate environment (`~/venvs/scfm`, PyTorch for CUDA 12.1 to support Pascal) after the GPU move; then run the leakage audit, the smoke test and the runs.
   - **Install done 22 Oct (installs only, nothing run; PI-approved, nice 19, 1 job at a time):** `~/venvs/scfm` (Python 3.10.22, uv). Exact versions in `docs/scfm_env_freeze.txt`. torch 2.3.1+cu121 (sees GTX 1080 and 1060, capability 6.1; `get_arch_list` has sm_60 but no sm_61, so Pascal runs on the sm_60 binary, which only the smoke test can confirm), scprint 1.0.13, scgpt 0.2.4 (installed `--no-deps`), torchtext 0.18.0, datasets 2.21.0, galore-torch 1.0, lightning 2.6.6.
   - **Resolved conflicts:** torch 2.5.1 was tried first but torchtext 0.18 (needed by scGPT) has an ABI mismatch with it, so torch was pinned to 2.3.1. scGPT's own pins (pytorch-lightning 1.9.5, scvi-tools<1, orbax, cell-gears) would downgrade lightning for scPRINT, so scGPT is installed without them; `pip check` lists those four as unmet. `import scgpt` works; scGPT's perturbation/scVI parts are not usable (not needed for the preregistered attention and embedding readouts). owlready2 needed `--no-build-isolation` (setuptools/packaging pin).
   - **Blocker for scPRINT:** `import scprint` fails until a lamindb instance is configured (django settings; needs `lamin init`, which writes a database and may fetch ontologies). Left undone; do it at the smoke-test step. Also unverified: scPRINT on torch 2.3.1 (scPRINT does not pin torch).
   - **Checkpoints:** scPRINT `medium-v1.5.ckpt` (222 MB) downloaded to `~/scfm_models/scprint/`, HF repo `jkobject/scPRINT` commit 8e1ad1a1bc40c1fc2f4f1700fe27923eea9a2bed, SHA-256 96ff7c5d9e1ad5076e818bde3c38f498c95cd354cc42029d90691b5ba2964ff8 (log in the preregistration change log before scoring). Not downloaded: `large-v1.ckpt` 1.2 GB and `small-v1.ckpt` 111 MB (secondary, can fetch on request); scGPT `scGPT_human` (source is a Google Drive link, not pinned in the preregistration; fetch after the PI confirms the source).
5. **GRNBoost2 (rerun running since 27 Oct 19:58):** PI approved the rerun everywhere with legacy early stopping (rule committed before any rerun, top of `docs/plans/grnboost2_diagnostic.md`).
   - **Diagnostic read (27 Oct):** `docs/plans/grnboost2_diagnostic.md`. Early stopping explains most of the gap: M3D SC 1.18 → 1.36, BEELINE +0.12/+0.27 log2 vs CLR; wrapper and seed are not at fault.
   - **Done:** `tools/grnboost2.py` defaults to `--early_stop legacy`; equivalence gate vs stock arboreto 0.1.6 + scikit-learn 1.2.2 passed (mDC_500, 3 seeds; `docs/grnboost2_equivalence/`). Old GRNBoost2 edges, scores and result tables moved or copied (never deleted) to `results/grnboost2_superseded/`; nothing reuses them. GENIE3 reproduces exactly from its seed (smoke test, mDC_500).
   - **Queue (in `jobs_pc/sched`, strictly sequential via `# after`, 16 threads each):** 0920 M3D production seed (`baselines.R`, gate on non-GRNBoost2 rows) → 0921 M3D seeds 1-3 (edges, must equal job 0894's legacy_es) → 0922 BEELINE GENIE3 cache → 0923 BEELINE 14 cases (gate) → 0924 BEELINE seeds 1-3 → 0925 K562 truth-gene run → 0926 EcoMAC → 0927 `macro_aupr.R` into `results/macro_aupr_v3` → 0928 `sota_table.R` into `results/sota_table_gb2_legacy.md`. Estimated 25-40 h wall in total.
   - **Resume:** if a job dies, requeue the same script from `jobs_pc/done/` unchanged (caches resume: GRNBoost2 edges per case or seed, GENIE3 and PIDC scores). Do not edit the paper, Table 2, `docs/whitepapers/sota_table.md` or the Collins checkpoint until 0928 finishes and the gates passed. A failed gate (rc 1) stops the chain on purpose: read its output in the job log.
   - **Open:** the K562 all-genes secondary run (`baselines_k562ess_genie3`, Holm family of 40, 8.5k genes) is not rerun (its macro tables hit the memory cap before; would need about 12-20 h more). Its GRNBoost2 is superseded and `macro_aupr_v3` excludes that family; PI to decide whether to rerun or drop it.
6. **Other planned analyses:**
   - log-scale depth vs module size;
   - directed DREAM5-style evaluation with ARACNE;
   - a method-neutral hub stratifier;
   - the EcoMAC core result, after checking array overlap with M3D;
   - diffCLR (paper 2).

## Machine rules (see `HANDOFF_crayon.md`)

- **Target utilization:** keep the machine at about 70–80% (standing rule).
- **Memory protection:** jobs run in `clrjobs.slice` (hard cap 95 GB). The scheduler's budget is 90 GB and 56 threads, with 15 GB always left free.
- **Auto-retry:** a job killed for memory is retried once at 1.5× its memory.
- **Status requests:** answered with the four-part table (running / queued / finished / planned, with data and a description of what each job tests).
