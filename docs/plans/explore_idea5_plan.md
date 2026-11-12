# Exploratory plan: self-supervised learned attention on the M3D compendium (ideation memo, section 2 idea 5)

*12 November 2026. Written and committed BEFORE any training or scoring. Not a preregistration. Exploratory: M3D is the project's tuning set and its RegulonDB truths are consumed, so nothing here is confirmatory, no p-values are computed, and nothing goes into the paper or the pre-read. After this commit the file is not edited; every deviation is logged in `docs/explore_idea5_results.md`.*

**Question.** The paper's thesis is "estimated, not learned" attention. The scFM arms used foreign pretraining, so they could not isolate the kernel. Here one attention layer over genes is trained self-supervised on the same compendium (M3D) and read out exactly like the CLR attention. If learning the kernel from the same arrays beats the CLR kernel on RegulonDB co-membership, the thesis is wrong at this scale; a tie supports it.

## 0. Choices that interpret the brief (stated so they can be overruled)

1. **Matrix: 4,297 genes x 466 replicate-averaged experiments**, not the 907-chip matrix. This is the matrix of `analysis/baselines.R`, `lab_splits.R`, `robust_z.R`, the truths, and the equivalence gate (`avg_E_coli_v4_Build_6_exps466probes4297.tab`, constant rows dropped; `docs/COMPUTE_CONSTRAINTS.md` D keeps the 466 as primary to avoid pseudo-replication). The brief's "907 arrays" is the chip-level file of the same genes.
2. **Gain unit = absolute pooled AUPR** (as `docs/plans/ecomac_nonoverlap_rescore.md`; the pass/fail thresholds of +/-0.005 are on this scale). Enrichment (AUPR / base rate; SC base rate 0.0943, all 0.1205) is reported alongside.
3. **Hub-stratified** = within-stratum percentile in 10 strata of H_i + H_j (`rz_hub_stratify` of `analysis/robust_z/robust_z_lib.R`). Two stratifiers, both fixed for every method, not model-specific: H = pi, the stationary distribution of the full-data CLR attention walk (the stratifier of the sota table; **primary**), and H = row mean of raw MI (method-neutral; sensitivity). Both are computed once from the full-data CLR.
4. **"Zero-step" read-out** = t = 1, the attention operator before any further walk, S = (P + P')/2. "Step 0" always means the training checkpoint before any optimiser step. **Pairwise read-out** = the score matrix C below, before the softmax and the walk (for CLR it is the CLR score).
5. **Matched comparators.** Primary comparator = full-data CLR attention at t = 9 (the project's attention, V0 of `robust_z`). Secondary (half-trained models): the CLR attention built from the same training half (the model's own step 0), same code path.

## 1. Data, splits, preprocessing

- Genes x experiments matrix X as in `robust_z.R` (M3D branch). 39 labs from the `experimenter:` field (unknown-lab experiments are their own lab).
- **One balanced lab split** by `lab_split()` of `lab_splits.R` (shuffle labs, assign each to the half with fewer experiments), `set.seed(20261112)`. Halves: 254 and 212 experiments. Both directions are run: H1 train / H2 validation and H2 train / H1 validation. Every seed uses this same split (the seed varies initialisation and batch order only).
- Per-gene standardisation uses the **training** samples' mean and SD (n - 1); validation samples are standardised with the training constants. For the full-data refit (section 4) the mean and SD come from all 466 experiments.
- Truths are never read by any training or stopping code.

## 2. Model

One attention layer over genes, G = 4,297 keys and queries (every gene is both). For a minibatch of samples with standardised expression `Xb` (G x B):

- **Logits** `L_ij = a_i * Z_ij / tau0 + D_ij` for `j != i`, with `D = Q K' / sqrt(d)`, then row-centred over `j != i` (`D_ij <- D_ij - mean_{j' != i} D_ij'`). `Q, K` are G x d learned embeddings; **d = 32 (primary) and d = 128**. `a_i = exp(gamma_i)`, `gamma` learned, one scalar per query gene (a query-specific temperature on the CLR kernel). `Z` is the package's CLR score matrix (Stouffer, B-spline MI with `bins = "hg"`, `transform = "none"`, truncated at 0, symmetrised, zero diagonal) computed on the **same samples the model is trained on**; `tau0` is the median over rows of the per-row softmax temperature that the readout's bisection assigns to the row of Z (about 10 effective partners over the top 50 positive scores). `Z/tau0` is a fixed buffer, not learned.
- **Two initialisations.** *CLR-initialised*: `gamma = 0`, `Q ~ N(0, 1)`, `K = 0` (so `D = 0` and `L = Z/tau0` exactly at step 0; gradients reach `K` immediately). *Scratch*: no `Z` term (`L = D`), `Q, K ~ N(0, 1)`, no `gamma`.
- **Attention** `A = softmax_j (L_ij)` over all `j != i` (diagonal masked to `-inf`; a standard full softmax). **Values = expression, signed:** `v_ij = s_ij x_j` with `s_ij = sign(Pearson r_ij)` on the training samples (+1 at r = 0), fixed, not learned (the package's `diffuse(values = "signed")`; MI and the readout are sign-blind, so an unsigned value would handicap any gene whose best partner is anti-correlated).
- **Prediction** `xhat_i(s) = sum_{j != i} A_ij s_ij x_j(s)`.

**Masking and loss (precise).** Gene i's own value in a sample is masked from its own prediction by the diagonal exclusion, in every sample, for every gene simultaneously (leave-one-gene-out; no other masking, no dropout). The prediction of gene i's profile across the samples of the batch uses only the other genes' values in the same samples. **Loss = mean squared error over all genes and all samples of the batch**, `mean_{i,s} (x_is - xhat_is)^2` on standardised expression. Held-out score: the same MSE on the validation samples (other lab half), reported also as R^2 = 1 - MSE / mean(x_val^2) (the predict-the-training-mean baseline).

**Readout of learned logits (identical to the project's).** The exported quantities are `D` (row-centred, float32) and `gamma`; the readout reassembles `L_ij = a_i Z_ij/tau0 + D_ij` in float64 (so that step 0 is exact).

1. Score matrix `C = pmax(tau0 * (L + L') / 2, 0)`, diagonal 0 (the CLR convention: symmetrise, truncate at 0; for the CLR-initialised model at step 0 `C = Z` exactly).
2. Per-gene softmax by the package's `clr:::.softmax_rows(C, keff = 10, cap = 50)` (top 50 positive scores, per-row temperature by bisection to exp(entropy) = 10), then the row normalisation and self-loop of `ClrAttention$.finish_operator` (`rz_operator`).
3. Lazy walk `P = 0.5 I + 0.5 A`; symmetric read-out `S_t = (P^t + P^t')/2` (`rz_lazy_P`, `rz_walk`).
4. Read-outs scored: **pairwise** (`C`), **t = 1**, **t = 9**, **t\*** (below).

**Label-free depth rule t\*** (`lab_splits.R`): for a model trained on half h, the criterion at each integer t = 1..32 is the AUPR of `S_t` (built from the model's operator) against the top 1% of the upper-triangular CLR scores computed on the **other** half (target = `Z` of the held-out half, 99th-percentile cut). t\* per config and seed = round of the geometric mean of the two directions' argmax t. The CLR comparators use t = 9 (project depth; the same rule gives 9 on M3D in `robust_z`). The refit models (section 4) use the t\* of their config and seed.

## 3. Training protocol (fixed here)

- **Optimiser** Adam (betas 0.9, 0.999, eps 1e-8, no weight decay), fp32, batch 64 samples (shuffled each epoch with a seeded generator; last batch partial), linear warm-up over the first 10 steps then constant learning rate. **Epoch cap 600; early stopping:** validation MSE after every epoch (and at step 0); stop after 30 consecutive epochs without a relative improvement of 1e-4 over the best; the checkpoint of the best validation MSE is exported (ties: earliest). If the best epoch is 0 the exported model is the init (reported as such).
- **Learning rate** selected per config (init x d) from {1e-3, 3e-3, 1e-2} by the lowest best-validation MSE averaged over the two directions, on a **tuning seed 0** that is not used afterwards; ties go to the smaller rate. Validation MSE only; no truth.
- **Seeds 1-5** (final), x 2 directions x 4 configs (CLR/scratch x d 32/128) = 40 half-trained models.
- Hardware: GPU device 1 (GTX 1060) only, torch of `~/venvs/scfm` (2.3.1+cu121); no new environment. Bit-identity across reruns is not promised (reported: seed-0 rerun difference).

## 4. Models that are scored

| Group | Models | Read-outs scored |
|---|---|---|
| CLR references | full-data CLR (V0 of `robust_z`); CLR of each training half (the step 0 of the half models) | pairwise, t = 1, t = 9 (half-CLR also its own t\*) |
| **Refits (primary)** | per config and seed: trained on **all 466 experiments** for `E` epochs, `E` = round(mean of the two directions' best epochs); CLR-initialised models start from the full-data CLR (exactly the paper's attention at step 0); 20 models | pairwise, t = 1, t = 9, t\* |
| Half-trained (secondary) | the 40 early-stopped models of section 3 (each in its own direction) | t = 9, t\* |
| Untrained control | scratch layer at step 0 (random init, same readout), 2 d x 5 seeds | pairwise, t = 1, t = 9 |

## 5. Equivalence gates (fail = stop and report, before any truth is scored)

G1. For the full data and for each training half: CLR-initialised step 0 (exported `D` = 0, `gamma` = 0) through the learned-logit readout reproduces the package (`ClrAttention$calibrate()`, `$build_operator(softmax_keff = 10, softmax_cap = 50)`, and the P^9 walk mass of `analysis/baselines.R`): max absolute difference of C, of the operator and of S_9 below 1e-8.
G2. The scoring code on the full-data CLR attention at t = 9 reproduces `results/robust_z/m3d/pooled_metrics.csv` (V0: pairwise and attention t = 9; raw, pi and MI-row-mean stratified; SC and all) to 1e-10 in AUPR, and the 100 jackknife values of `results/robust_z/m3d`'s stream.
G3. The R port of the per-row temperature reproduces `clr:::.softmax_rows` weights (max difference below 1e-10).

## 6. Truths and metrics

- M3D RegulonDB truths as in `robust_z.R`: **SC** (strong evidence C, S; base rate 0.0943) and **all** evidence (0.1205), regulator-agnostic co-membership pairs, same-operon pairs removed (`comembership_pairs`).
- **Pooled AUPR** (`ap_fast`) raw, hub-stratified by pi (primary) and by MI row mean. **Paired delete-half gene jackknife, 100 replicates, seed 20261016, same halves as `robust_z.R` for every method** (SC drawn first, then all); the interval of a contrast is the difference plus/minus 1.96 SD of the paired replicate differences (the `robust_z` convention; an uncalibrated, gene-resampling-only interval). Contrasts: each model minus full-data CLR attention at t = 9; half-trained models also minus their half-CLR t = 9.
- Seed-to-seed spread: mean, SD and range of the primary contrast over the 5 seeds, reported next to the SD of the jackknife replicates.

## 7. Pass/fail, stated before compute

Primary cell: **CLR-initialised, d = 32, refit, read out at its t\*; hub-stratified (pi) strong-evidence (SC) pooled AUPR; contrast = learned minus full-data CLR attention at t = 9; per seed with its own jackknife interval.**

- **Per seed:** *win* if the contrast is >= +0.005 and the lower end of its 95% jackknife interval is above 0; *loss* if the contrast is below -0.005; *tie* if it lies within +/-0.005 (|contrast| <= 0.005); a contrast >= +0.005 with an interval that includes 0 is *unresolved gain* (not a win).
- **Overall:** **learned wins** if at least 3 of the 5 seeds are wins; **learned loses** if at least 3 of 5 are losses; **tie** if at least 3 of 5 are ties; otherwise **mixed**, with the counts reported.
- The same rule is reported, descriptively and not as the verdict, for d = 128, for the scratch models, for the read-out at t = 9, for the MI-row-mean stratifier, for the all-evidence truth, and for the half-trained models against their half-CLR comparator. If the pi and MI-row-mean verdicts differ, the report says so; the verdict stays with pi.
- The expectation in the memo (learned ~ CLR attention on regulons) is a prediction, not a threshold.

## 8. Smoke test (before the real run; no truth)

On 600 randomly chosen genes (seeded), end to end: prep (MI, CLR, tau0), a short training run (a few epochs, both inits, both d), export, the learned-logit readout, gate G1 and G3 on the subset, the depth curve. No truth file is opened, no AUPR is computed. Timing and GPU memory are read from it to set the declarations below.

## 9. Jobs, resources, disk

- Scheduler (`tools/scheduler.sh`, queue `jobs_pc/sched`, declarations as the scFM jobs): prep (R, `# cpus: 8`, `# mem` about 8), training (python, **GPU device 1 only**, `# cpus: 2`, `# mem` about 8, `# after` prep), readout + depth curves + scoring (R, `# cpus: 8`, `# mem` about 14, `# after` training). Each job runs in its own capped scope at 1.3 x the declared memory; the declarations are set from the smoke test and from `jobs_pc/peaks.tsv` after the run; peak memory and runtime go in the results document. Other workers run CPU jobs concurrently; nice 19 throughout.
- Disk: exported logits are G x G float32 (74 MB) per model, about 5 GB transient under `results/idea5/logits/`, **deleted after scoring**; kept: small CSVs, logs and the document (a row for the transient space is added to `docs/disk_budget.md`). No checkpoint is committed.

## 10. Post-hoc check on EcoMAC non-overlap arrays (labelled post hoc; only if cheap)

The trained models are gene-specific (embeddings of the 4,297 M3D probes), so after training their operator does not depend on any array and scoring it against RegulonDB on EcoMAC would repeat the M3D scoring on fewer genes and add no information. The well-defined check is **label-free held-out reconstruction**: on the EcoMAC arrays without an M3D counterpart (subset B of `docs/plans/ecomac_nonoverlap_rescore.md`; shared genes only; genes standardised on those arrays), the MSE / R^2 of the M3D-trained refit models against the M3D full-data CLR attention at step 0. It is done only if the subset-B matrix can be built by calling the existing pipeline unchanged (no edits to `analysis/baselines.R`) in well under an hour of one core; otherwise it is skipped and the results document says so. No truth is scored on EcoMAC.

## 11. Outputs and reporting

Scripts under `analysis/idea5/`; tables under `results/idea5/` (small CSVs force-added); `docs/explore_idea5_results.md`: method as run, gates, training-curve summary (best epochs, held-out R^2 at step 0 and at the best epoch, fraction of rows that end up isolated), jobs and peak memory, all tables, the literal pass/fail, deviations, and a reading limited to the tables. `docs/STATUS.md`, other workers' files, DREAM5, the paper and the pre-read are not touched.
