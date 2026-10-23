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
- [ ] EcoCyc complex truth build; commit its SHA-256
- [ ] Implement the scoring (jackknife, bootstrap, Holm per the preregistration); reviewer who has seen no E. coli result checks script against the file
- [ ] Truth phase, once: E. coli, then RPE1, then K562 (in-sample)
- [ ] Write-up: verdict per the preregistration's outcome table

If a job died mid-step, rerun that step; do not skip a gate. No truth may be read before the three boxes above it are ticked.

## Next

1. **CORUM order-of-operations expert (PI, 22 Oct).** Members of a complex share one latent variable, so all their pairs are strongly coupled. The softmax and the walk may spread attention across the clique and then into the wider translation program. The question is whether collapsing strongly coupled cliques into latent variables before attending, or attending on residuals after removing a clique's shared factor, fixes the loss. A max-effort expert analyses the CORUM results and interviews the PI (questions relayed in the session).
2. **Figures:**
   - finish the revisions from the fresh-eyes reviews (Shapley, MI vs ρ, lexA scatters; one colour key and finding titles);
   - rebuild graphical-abstract panel 2 as the CLR criss-cross;
   - build the paper figures that have legends but no images (Figs 1, 2, 4–7).
3. **keff test, phase 2:** commit the K562 selections, then score truths (the scoring step refuses to run unless the selections match the committed copies).
4. **Learned-attention comparator (scPRINT, scGPT):** preregistered. Install a separate environment (`~/venvs/scfm`, PyTorch for CUDA 12.1 to support Pascal) after the GPU move; then run the leakage audit, the smoke test and the runs.
   - **Install done 22 Oct (installs only, nothing run; PI-approved, nice 19, 1 job at a time):** `~/venvs/scfm` (Python 3.10.22, uv). Exact versions in `docs/scfm_env_freeze.txt`. torch 2.3.1+cu121 (sees GTX 1080 and 1060, capability 6.1; `get_arch_list` has sm_60 but no sm_61, so Pascal runs on the sm_60 binary, which only the smoke test can confirm), scprint 1.0.13, scgpt 0.2.4 (installed `--no-deps`), torchtext 0.18.0, datasets 2.21.0, galore-torch 1.0, lightning 2.6.6.
   - **Resolved conflicts:** torch 2.5.1 was tried first but torchtext 0.18 (needed by scGPT) has an ABI mismatch with it, so torch was pinned to 2.3.1. scGPT's own pins (pytorch-lightning 1.9.5, scvi-tools<1, orbax, cell-gears) would downgrade lightning for scPRINT, so scGPT is installed without them; `pip check` lists those four as unmet. `import scgpt` works; scGPT's perturbation/scVI parts are not usable (not needed for the preregistered attention and embedding readouts). owlready2 needed `--no-build-isolation` (setuptools/packaging pin).
   - **Blocker for scPRINT:** `import scprint` fails until a lamindb instance is configured (django settings; needs `lamin init`, which writes a database and may fetch ontologies). Left undone; do it at the smoke-test step. Also unverified: scPRINT on torch 2.3.1 (scPRINT does not pin torch).
   - **Checkpoints:** scPRINT `medium-v1.5.ckpt` (222 MB) downloaded to `~/scfm_models/scprint/`, HF repo `jkobject/scPRINT` commit 8e1ad1a1bc40c1fc2f4f1700fe27923eea9a2bed, SHA-256 96ff7c5d9e1ad5076e818bde3c38f498c95cd354cc42029d90691b5ba2964ff8 (log in the preregistration change log before scoring). Not downloaded: `large-v1.ckpt` 1.2 GB and `small-v1.ckpt` 111 MB (secondary, can fetch on request); scGPT `scGPT_human` (source is a Google Drive link, not pinned in the preregistration; fetch after the PI confirms the source).
5. **GRNBoost2:** read the diagnostic, then decide how to rerun GRNBoost2 fairly (legacy early stopping or none) everywhere it appears.
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
