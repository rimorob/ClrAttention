# Spectral-embedding exploratory pass: code, jobs, resources

Implements `docs/plans/spectral_embedding_pass_plan.md` (read it first: definitions, ledger, descriptive success criterion, multiplicity). Built and smoke-tested on 21 Nov 2026. **Nothing was submitted to the scheduler, no real truth was scored, Jurkat was not touched, no download or deletion was made, nothing was committed.** Job files are staged here, not in `jobs_pc/`.

## Scripts

| Script | Role |
|---|---|
| `sp_lib.R` | Shared code. Sources (does not redefine) the set-up of `analysis/beeline.R`, `score_truth`, `ap_fast`, `analysis/scgpt_ext/ext_lib.R` (`ext_edge`, `ext_complex_truth`, `ext_score_complex`, `ext_plant_complexes`), `gl_make_P`, `ad_Psym`, `ad_eigs`, `sk_c2z`. New: unit loader (copy of the expression construction of `score_ext.R`, gated by the stored `tcaclr2_edge_t2.csv`), spectral embedding, the seven read-outs, the two ring families, the persistent-homology validator, the 35-variant registry (pgmm_k50 not run). Jurkat cannot be opened through `sp_dataset` (only the development datasets are allowed unless `sp_score.R --heldout YES`). |
| `sp_score.R` | Evaluation driver (same units, scGPT gene set, metrics as `score_ext.R`). `--truth real|synthetic` (synthetic reads no truth file), `--datasets`, `--units`, `--only` (read-outs), `--krules`, `--variants` (explicit list; required with `--heldout YES`), `--expr`. Output `results/spectral_pass/<dataset>/`: `metrics_spectral.csv` (same columns as `metrics_ext.csv`, `method` = variant name), `variant_info.csv` (k, components, model, BIC, `g1_wins`, `at_max`, seconds), `eig_info.csv` (edge, t2, eigenvalues above the edge, the first 60 eigenvalues), `ring_planes.csv` (per plane: both ring families, BIC differences, fitted radii, H1 persistence and its permuted-data null), `gates.csv`, `unit_timing.csv`, `complex_truth/`. Never writes into `results/scgpt_ext`. |
| `sp_controls.R` | `--mode planes` (plane-level tests of the two ring families), `positive` (planted clusters, subspace clusters, ring; synthetic embedding and real embedding), `negative` (per-gene permuted expression), `both`. No truth file. |
| `sp_summary.R` | Descriptive tables: per variant, dataset, truth, the mean-over-units AUPR/base of the variant, of tcaCLR2, scGPT token cosine, scGPT attention (and aCLR, tcaCLR, CLR-NS, scGPT hidden for context) read from `results/scgpt_ext/<dataset>/metrics_ext.csv`, and the ratios; plus the mean log2 ratio over units and the three primary truths. No test, no threshold. Jurkat only with `--datasets Jurkat`. |
| `jobs/` | 7900-7904 (one per development dataset), 7905 (plane tests + positive control), 7906 (negative controls, five datasets), 7910 (summary), **7990 (Jurkat, HELD, sealed)**. |

## Model list (35 listed variants = 5 embeddings x 7 read-outs, 34 scored: pgmm is not run for k > 30; exact definitions in plan sections 3-4)
Embedding: `gl_make_P(ns_clr(X), 10, 50)` -> `ad_Psym` -> top eigenpairs, trivial one dropped, `Y[, j] = lambda_j^(t2/2) v_j`, t2 and edge as tcaCLR2. k: edge rule (`kE`: eigenvalues above the permuted-data edge among the first 60, cap 50) and 5, 10, 20, 50. Read-outs: `cos`; `cos_clr2` (`sk_c2z` on the cosine); `mclust` (BIC, G = 2..30, 14 models); `hddc` (HDclassif, K = 2..20, 14 models, BIC for K and each component's dimension); `pgmm` (G = 2..4, q = 1..3, 8 families, only k <= 30); `ring` (family A: polar transform, bulk + thin rings at the origin, adjacent planes (j, j+1), BIC against mclust, accept at +10); `rmix` (family B: mixture of 1-2 radial Gaussians with fitted centre, proper planar density, BIC against mclust, accept at +10). Validator, not a read-out: H1 persistence (TDAstats) per plane against the permuted-data null.

## Inputs
Expression: as `score_ext.R` (`results/scfm_*/` prep directories and `expr/` files; K562 and hESC/hHep from their h5ad / BEELINE inputs); the scGPT gene list of each unit (`~/scfm_models/step4_*/scgpt_A/<case>.genes.txt`, gene list only); stored `results/scgpt_ext/<dataset>/tcaclr2_edge_t2.csv` (gate); truths only with `--truth real` (`data/human_truth/{dorothea,collectri}.tsv`, Complex Portal via `analysis/scgpt_ext/build_complex_truth.R`, KD bulk file); `results/scgpt_ext/<dataset>/metrics_ext.csv` (comparison columns, summary only).

## Declared resources
No forked workers (a forked child of a process that ran clr can hang; nothing here forks), so the process tree is the main process. Measured peak RSS (`/usr/bin/time -v`, whole run): 0.81 GB (negative control hESC_500), 0.87 GB (plane tests), 0.90-0.96 GB (hESC_500 score), 1.04 GB (negative control HepG2), 1.09 GB (positive control). Declared 4 GB (3.7x the largest), 1 cpu (single-threaded: BLAS and OpenMP set to 1), 7910 2 GB.

| Job | mem (GB) | cpus | after |
|---|---|---|---|
| 7900 hESC, 7901 hHep, 7902 K562, 7903 HepG2, 7904 RPE1 (score, real truths) | 4 | 1 | none |
| 7905 plane tests + positive control | 4 | 1 | none |
| 7906 negative controls (five datasets, unit 1) | 4 | 1 | none |
| 7910 summary | 2 | 1 | 7900,...,7904 (comma-separated) |
| 7990 Jurkat (HELD, frozen variants only, refuses to run until `FROZEN` is set) | 4 | 1 | 7910 |

## Smoke tests (synthetic truths only; outputs in `results/spectral_pass/smoke/`)
Run with up to 6 concurrent single-threaded processes on a shared machine (load average above 6), so times are pessimistic.

| Run | Wall time | Peak RSS |
|---|---|---|
| plane tests (9 planted 2-D cases, both families, validator) | 98 s | 0.87 GB |
| score, hESC_500, all 35 variants, synthetic truths | 43.5 min (gate: edge 0.837826 = stored, t2 17) | 0.96 GB |
| score, hESC_500, 4 variants (`--variants`), synthetic truths | 65 s | 0.90 GB |
| positive control (synthetic + real HepG2 unit 1 embedding) | 40.5 min | 1.09 GB |
| negative control, HepG2 unit 1 / hESC_500 | 23.7 min / 15.9 min | 1.04 / 0.81 GB |
| score, HepG2 unit 1, all 35 variants, synthetic truths (before the pgmm k <= 30 rule) | 98.8 min (gate: edge difference 1e-16); 34 variants scored, `pgmm_k50` ran 43 min and returned no model | 1.15 GB |

Where the time goes (hESC_500, 841 genes, k_edge 15; seconds, loaded machine): CLR-NS, edge and eigenpairs about 30; ring planes 1-49 (both families, mclust per plane, validator) 309; permuted-data validator null about 150; mclust 42-228 per k; hddc 19-66 per k; pgmm 81 (k = 5), 158 (10), 198 (15), 259 (20), 956 (k = 50, hESC_500). HepG2 unit 1 (1,253 genes, k_edge 22): ring planes 444, mclust 82-354, hddc 23-82, pgmm 201 (k = 5), 295 (10), 673 (22), 597 (20), **2,588 and no model at k = 50**. After the smoke run pgmm was limited to k <= 30 (plan section 3), which removes the k = 50 call; pgmm is then still about 40% of the per-unit time.

## Estimated full-size time
Measured per unit with the k = 50 pgmm call: 43.5 min (hESC_500, 841 genes) and 98.8 min (HepG2, 1,253 genes), on a machine with load average above 6. With pgmm limited to k <= 30 (about 16 min and 43 min less): about 28 min and 55 min per unit, perhaps less on an idle machine. Units: hESC 2 (the 1,000-gene case is larger), hHep 2, K562 5, HepG2 5, RPE1 5. Estimated wall time per job: 7900 about 1.5 h, 7901 about 1.5 h, 7902 to 7904 about 4.5 h each (the five jobs are independent and can run concurrently, 1 cpu each). Controls: 7905 about 45 min, 7906 about 2 h. Summary seconds. Jurkat (frozen subset only): scales with the number of frozen variants (the ring families need all 49 planes at about 8 min per unit plus about 3 min for the validator null; pgmm_kE and pgmm_k20 about 10 min each). To shorten: `--only` / `--variants` drop read-outs, but any dropped variant must be recorded in the plan's ledger before scoring (the 35-variant family is the pre-specified one); pgmm is the obvious candidate (the plan allows "HDDC and/or pgmm").

## Controls found in the smoke tests (descriptive; for the lead to read, no interpretation made here)
See the report of the build and `results/spectral_pass/smoke/controls/`. Items to know before reading real scores: (1) at the edge rule the per-gene permuted-data control has 0 eigenvalues above the edge (HepG2 unit 1 and hESC_500), so the edge-rule variants are not scored there; (2) at forced k = 5, 10, 20 the mixtures still return components on permuted data (mclust G = 2-17, hddc K = 3-19, pgmm at its upper limit G = 4) with `g1_wins` FALSE, and ring family A accepts 0-5 planes of 19; family B accepts 0-1; (3) planted ring: found by both families in the synthetic embedding (plane 7, validator above the null), not accepted by either family when planted in the real HepG2 embedding (validator above the null, H1 0.94 vs null maximum 0.41); mclust, hddc and pgmm recover planted clusters and the ring pairs in the real embedding at kE and k = 10.

## Things that did not work as first written, or are odd
- `pgmm::pgmmEM` is far slower than the other mixtures (above); G is at its upper limit (4) in every smoke run (`at_max` TRUE).
- mclust and hddc return many components on the real embeddings (G 6-22 at k <= 22, K 9-19); mclust hit its limit G = 30 at k = 50 on HepG2 (`at_max` TRUE); the BIC prefers many non-Gaussian-tail components.
- The literal radial density N(r; mu, s) / (2 pi r) of the PI is not integrable at r = 0, so family B uses the normalised form f(x) = exp(-(r - mu)^2 / 2 s^2) / Z (plan section 3); equal to it to first order for thin rings.
- Family B's likelihood can fail to produce any valid ring component (`B_dbic = -Inf` in `ring_planes.csv`), including in the plane of a planted ring in the real embedding.
- A running R script reads its file incrementally: editing `sp_score.R` while a smoke run was going broke the trailing `writeLines`/`say` lines of that run (the metric files had already been written inside the dataset loop). Do not edit the scripts while jobs run.
- Smoke metric values come from synthetic truths (planted from CLR-NS and the spectral cosine) and mean nothing; `smoke/summary_codecheck/` only checks that the summary script runs (variant values synthetic, comparison columns real).
